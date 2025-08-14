from __future__ import annotations
from typing import TYPE_CHECKING
import panel as pn
import holoviews as hv
import param
import pandas as pd
import numpy as np
from holoviews.streams import Selection1D
import networkx as nx

if TYPE_CHECKING:
    from clusterf.app import ClusterFApp


class SuperClusterViewer(param.Parameterized):
    """
    A modular viewer for displaying cluster graphs using HoloViews.
    Handles visualization of compound clusters with interactive features.
    """
    
    # Parameters
    app: "ClusterFApp" = param.Parameter(default=None, doc="Reference to main ClusterF application")
    visible = param.Boolean(default=True, doc="Whether the cluster viewer is visible")
    
    # Plot data
    cluster_data = param.DataFrame(default=None, allow_None=True, doc="Current cluster data for plotting")
    color_dict = param.Dict(default={}, doc="Color mapping for cluster categories")
    selected_nodes = param.List(default=[], doc="Currently selected cluster nodes")
    
    # Plot configuration
    plot_width = param.Integer(default=600, bounds=(400, 1600), doc="Width of the cluster plot")
    plot_height = param.Integer(default=600, bounds=(300, 1200), doc="Height of the cluster plot")
    point_size = param.Integer(default=25, bounds=(5, 50), doc="Size of cluster nodes")
    
    def __init__(self, **params):
        super().__init__(**params)
        
        # Graph components
        self.G = None  # NetworkX graph
        self.pos = None  # Node positions
        self.cluster_color_map = {}  # Mapping of nodes to colors
        
        # HoloViews elements
        self.graph = None
        self.selection = None
        
        # Main plot pane
        self.plot = pn.pane.HoloViews(
            object=None, 
            sizing_mode="stretch_both", 
            margin=0,
            name="Cluster Graph"
        )
        
    def initialize_graph(self, graph):
        """
        Initialize the graph visualization with a NetworkX graph.
        
        Args:
            graph: NetworkX graph object
        """
        self.G = graph
        if not self.G or len(self.G.nodes()) == 0:
            self.plot.object = hv.Text(0.5, 0.5, "No graph data available").opts(
                width=self.plot_width, height=self.plot_height, xaxis=None, yaxis=None
            )
            return
            
        # Generate positions for nodes using networkx
        self.pos = nx.spring_layout(self.G)
        
        # Update colors based on current color dictionary
        self._update_node_colors()
        
        # Create the visualization using hv.Graph.from_networkx
        self._create_graph_visualization()
        
        # Initialize with all nodes selected
        if self.selection:
            self.selection.update(index=list(range(len(self.G.nodes()))))
        
    def _update_node_colors(self):
        """Update node colors based on current color dictionary and library data."""
        if not self.app or not hasattr(self.app, 'library') or not self.app.library:
            # Default colors if no library available
            self.cluster_color_map = {node: "#999999" for node in self.G.nodes()}
            return
            
        library = self.app.library
        
        # Only get color dictionary from color picker if we don't already have one
        if not self.color_dict and hasattr(self.app, 'color_picker') and self.app.color_picker.color_dict:
            self.color_dict = self.app.color_picker.color_dict.copy()
        
        if not hasattr(library, 'df') or "Cluster" not in library.df.columns:
            self.cluster_color_map = {node: "#999999" for node in self.G.nodes()}
            return
            
        try:
            # Create a mapping of clusters to their category counts
            categories = library.df[library.df["Category"] != "Miss"]
            cluster_category_counts = (
                categories.groupby("Cluster")["Category"]
                .value_counts()
                .unstack(fill_value=0)
            )
            
            # Apply the color determination function to each cluster
            self.cluster_color_map = {
                cluster: self._determine_cluster_color(counts)
                for cluster, counts in cluster_category_counts.iterrows()
            }
            
        except Exception as e:
            print(f"Error updating node colors: {e}")
            self.cluster_color_map = {node: "#999999" for node in self.G.nodes()}
    
    def _determine_cluster_color(self, category_counts):
        """Determine cluster color based on the category with the most compounds."""
        try:
            # Find the category with the maximum count
            if category_counts.sum() == 0:
                return "#999999"

            max_category = category_counts.idxmax()
            return self.color_dict.get(max_category, "#999999")

        except (AttributeError, TypeError, ValueError):
            return "#999999" 
    
    def _create_graph_visualization(self):
        """Create the HoloViews graph visualization using hv.Graph.from_networkx."""
        if not self.G:
            return
            
        # Add color attributes to the graph nodes
        for node in self.G.nodes():
            self.G.nodes[node]['color'] = self.cluster_color_map.get(node, "#999999")
        
        # Create the graph using HoloViews
        self.graph = hv.Graph.from_networkx(self.G, self.pos).opts(
            # Node styling
            node_size=self.point_size,
            node_color='color',
            node_line_color='black',
            node_line_width=0,
            
            # Edge styling
            edge_line_width=0.5,
            edge_line_color='gray',
            edge_line_alpha=0.6,
            
            # General plot styling
            width=self.plot_width,
            height=self.plot_height,
            min_width=300,
            min_height=300,
            max_width=600,
            responsive=True,
            xaxis=None,
            yaxis=None,
            
            # Tools
            tools=["tap", "box_select", "lasso_select"],
            active_tools=["tap"]
        )
        
        # Add labels as a separate overlay
        node_positions = np.array([self.pos[n] for n in self.G.nodes()])
        node_labels = list(self.G.nodes())
        self.labels = hv.Labels(
            {("x", "y"): node_positions, "text": node_labels},
            ["x", "y"],
            "text"
        ).opts(
            text_color="black", 
            text_font_size="10pt"
        )
        
        # Combine graph and labels
        self.graph = self.graph * self.labels
        
        # Set up selection stream
        self.selection = Selection1D(source=self.graph)
        self.selection.param.watch(self._on_selection_change, "index")
        
        # Update the plot
        self._update_plot()
    
    def _on_selection_change(self, event):
        """Handle changes in node selection."""
        indices = self.selection.index if self.selection else []
        
        if indices:
            node_list = list(self.G.nodes())
            self.selected_nodes = [node_list[i] for i in indices]
            # Notify app of selection change
            if hasattr(self.app, '_on_cluster_selection_change'):
                self.app._on_cluster_selection_change(self.selected_nodes)
        else:
            self.selected_nodes = []
            # Clear selection in app
            if hasattr(self.app, '_on_cluster_selection_change'):
                self.app._on_cluster_selection_change([])
        
        # Update the plot visualization
        self._update_plot()
    
    def _update_plot(self):
        """Update the graph plot with current selection state."""
        if not self.graph:
            return
            
        if len(self.selected_nodes) > 0:
            # For selection styling, we'll recreate the graph with updated styling
            # Add selection attributes to the graph nodes
            for node in self.G.nodes():
                if node in self.selected_nodes:
                    self.G.nodes[node]['selected'] = True
                    self.G.nodes[node]['size'] = self.point_size * 2
                    self.G.nodes[node]['line_width'] = 2
                else:
                    self.G.nodes[node]['selected'] = False
                    self.G.nodes[node]['size'] = self.point_size
                    self.G.nodes[node]['line_width'] = 0
            
            # Update edge attributes based on selection
            for edge in self.G.edges():
                if edge[0] in self.selected_nodes or edge[1] in self.selected_nodes:
                    self.G.edges[edge]['connected'] = True
                    self.G.edges[edge]['line_width'] = 2
                    self.G.edges[edge]['line_color'] = 'black'
                    self.G.edges[edge]['alpha'] = 1.0
                else:
                    self.G.edges[edge]['connected'] = False
                    self.G.edges[edge]['line_width'] = 1
                    self.G.edges[edge]['line_color'] = 'gray'
                    self.G.edges[edge]['alpha'] = 0.1
            
            # Create the styled graph
            styled_graph = hv.Graph.from_networkx(self.G, self.pos).opts(
                # Node styling
                node_size=hv.dim('size'),
                node_color='color',
                node_line_color='black',
                node_line_width=hv.dim('line_width'),
                
                # Edge styling
                edge_line_width=hv.dim('line_width'),
                edge_line_color=hv.dim('line_color'),
                edge_line_alpha=hv.dim('alpha'),
                
                # General plot styling
                width=self.plot_width,
                height=self.plot_height,
                min_width=300,
                min_height=300,
                max_width=600,
                responsive=True,
                xaxis=None,
                yaxis=None,
                
                # Tools
                tools=["tap", "box_select", "lasso_select"],
                active_tools=["tap"]
            )
            
            # Add labels overlay
            styled_graph = styled_graph * self.labels
            
            self.plot.object = styled_graph
        else:
            # No selection - show default view
            self.plot.object = self.graph
    
    def update_colors(self, color_dict):
        """
        Update the color mapping and refresh the visualization.
        
        Args:
            color_dict: Dictionary mapping categories to colors
        """
        self.color_dict = color_dict.copy()
        
        if self.G:
            self._update_node_colors()
            self._create_graph_visualization()
    
    def set_selection(self, node_indices):
        """
        Set the selected nodes by index.
        
        Args:
            node_indices: List of node indices to select
        """
        if self.selection and node_indices:
            self.selection.update(index=node_indices)
        elif self.selection:
            self.selection.update(index=[])
    
    def clear_selection(self):
        """Clear all node selections."""
        self.set_selection([])
    
    def get_view(self):
        """Return the Panel component for embedding in layouts."""
        return self.plot
    
    @property
    def view(self):
        """Property accessor for the Panel component."""
        return self.get_view()