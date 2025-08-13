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
        self.node_positions = None  # Array of node positions
        self.node_labels = []  # List of node labels
        self.edges_data = []  # Edge data for visualization
        self.edges_array = None  # Array of edge coordinates
        self.cluster_color_map = {}  # Mapping of nodes to colors
        
        # HoloViews elements
        self.points = None
        self.labels = None
        self.initial_edges = None
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
        
        # Convert node positions and edges to lists of coordinates for HoloViews
        self.node_positions = np.array([self.pos[n] for n in self.G.nodes()])
        self.node_labels = list(self.G.nodes())
        self.edges_data = [
            (self.pos[edge[0]], self.pos[edge[1]]) for edge in self.G.edges()
        ]
        
        # Update colors based on current color dictionary
        self._update_node_colors()
        
        # Extract x and y coordinates for edges
        self.edges_array = np.array(
            [(x1, y1, x2, y2) for ((x1, y1), (x2, y2)) in self.edges_data]
        )
        
        # Create the visualization elements
        self._create_visualization_elements()
        
        # Initialize with all nodes selected
        self.selection.update(index=list(range(len(self.node_labels))))
        
    def _update_node_colors(self):
        """Update node colors based on current color dictionary and library data."""
        if not self.app or not hasattr(self.app, 'library') or not self.app.library:
            # Default colors if no library available
            self.cluster_color_map = {node: "#999999" for node in self.node_labels}
            return
            
        library = self.app.library
        
        # Only get color dictionary from color picker if we don't already have one
        if not self.color_dict and hasattr(self.app, 'color_picker') and self.app.color_picker.color_dict:
            self.color_dict = self.app.color_picker.color_dict.copy()
        
        if not hasattr(library, 'df') or "Cluster" not in library.df.columns:
            self.cluster_color_map = {node: "#999999" for node in self.node_labels}
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
            self.cluster_color_map = {node: "#999999" for node in self.node_labels}
    
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
    
    def _create_visualization_elements(self):
        """Create the HoloViews visualization elements."""
        # Get colors for nodes
        cmap = [self.cluster_color_map.get(node, "#999999") for node in self.node_labels]
        
        # Create data frame for points
        data = pd.DataFrame({
            "x": self.node_positions[:, 0],
            "y": self.node_positions[:, 1],
            "color": cmap,
        })
        
        # Create HoloViews Points for the nodes
        self.points = hv.Points(data, ["x", "y"]).opts(
            size=self.point_size,
            width=self.plot_width,
            height=self.plot_height,
            responsive=True,
            xaxis=None,
            yaxis=None,
            color="color",
            tools=["tap", "box_select", "lasso_select"],
            active_tools=["tap"],
        )
        
        # Create labels for nodes
        self.labels = hv.Labels(
            {("x", "y"): self.node_positions, "text": self.node_labels},
            ["x", "y"],
            "text",
        ).opts(text_color="black", text_font_size="10pt")
        
        # Create HoloViews Segments for the edges
        self.initial_edges = hv.Segments(self.edges_array).opts(
            line_width=0.5, color="gray"
        )
        
        # Set up selection stream
        self.selection = Selection1D(source=self.points)
        self.selection.param.watch(self._on_selection_change, "index")
        
        # Update the plot
        self._update_plot()
    
    def _on_selection_change(self, event):
        """Handle changes in node selection."""
        indices = self.selection.index if self.selection else []
        
        if indices:
            self.selected_nodes = [self.node_labels[i] for i in indices]
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
        if not self.points or not self.initial_edges or not self.labels:
            return
            
        if len(self.selected_nodes) > 0:
            # Get selected and non-selected nodes
            n_nodes = len(self.node_positions)
            selected_indices = self.selection.index if self.selection else []
            non_selected_indices = [
                i for i in range(n_nodes) if i not in selected_indices
            ]
            
            # Create separate point objects for selected and non-selected nodes
            if selected_indices:
                selected_data = self.points.data.loc[selected_indices]
                selected_points = hv.Points(selected_data).opts(
                    size=self.point_size * 2, 
                    color="color", 
                    line_color="black",
                    line_width=2,
                    responsive=True,
                    xaxis=None,
                    yaxis=None,
                    tools=["tap", "box_select"],
                    active_tools=["tap"],
                )
            else:
                selected_points = hv.Points([]).opts(size=0)
            
            if non_selected_indices:
                non_selected_data = self.points.data.loc[non_selected_indices]
                non_selected_points = hv.Points(non_selected_data).opts(
                    size=self.point_size, 
                    color="color",
                    alpha=0.6
                )
            else:
                non_selected_points = hv.Points([]).opts(size=0)
            
            # Get connected and non-connected edges
            connected_edges = []
            non_connected_edges = []
            
            for edge in self.G.edges():
                if any(node in edge for node in self.selected_nodes):
                    x1, y1 = self.pos[edge[0]]
                    x2, y2 = self.pos[edge[1]]
                    connected_edges.append((x1, y1, x2, y2))
                else:
                    x1, y1 = self.pos[edge[0]]
                    x2, y2 = self.pos[edge[1]]
                    non_connected_edges.append((x1, y1, x2, y2))
            
            # Create edge segments
            if connected_edges:
                connected_segments = hv.Segments(connected_edges).opts(
                    line_width=2, color="black"
                )
            else:
                connected_segments = hv.Segments([]).opts(line_width=0)
                
            if non_connected_edges:
                faded_segments = hv.Segments(non_connected_edges).opts(
                    line_width=1, color="gray", alpha=0.1
                )
            else:
                faded_segments = hv.Segments([]).opts(line_width=0)
            
            # Combine all elements
            self.plot.object = (
                connected_segments * faded_segments * 
                non_selected_points * selected_points * self.labels
            )
        else:
            # No selection - show default view
            self.plot.object = self.initial_edges * self.points * self.labels
    
    def update_colors(self, color_dict):
        """
        Update the color mapping and refresh the visualization.
        
        Args:
            color_dict: Dictionary mapping categories to colors
        """
        self.color_dict = color_dict.copy()
        
        if self.G and self.node_labels:
            self._update_node_colors()
            self._create_visualization_elements()
    
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