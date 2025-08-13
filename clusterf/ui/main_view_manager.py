from __future__ import annotations
from typing import TYPE_CHECKING
import panel as pn
import param
from clusterf.ui.cluster_viewer import SuperClusterViewer
from clusterf.ui.responsive_layout_manager import ResponsiveLayoutManager
from clusterf.ui.compound_grid import CompoundGrid
from clusterf.ui.category_histogram import CategoryHistogram
from clusterf.ui.compound_data_chart import CompoundDataChart
from clusterf.ui.compound_table import CompoundTable

if TYPE_CHECKING:
    from clusterf.app import ClusterFApp


class MainViewManager(param.Parameterized):
    """
    Manages the main content area of the ClusterF application.
    Coordinates between different visualization components like cluster viewer,
    compound tables, histograms, etc.
    """
    
    app: "ClusterFApp" = param.Parameter(default=None, doc="Reference to main ClusterF application")
    visible = param.Boolean(default=True, doc="Whether the main view is visible")
    clusters_built = param.Boolean(default=False, doc="Whether super clusters have been built")
    
    def __init__(self, app: "ClusterFApp", **params):
        super().__init__(**params)
        self.app = app
        
        # Initialize visualization components
        self.cluster_viewer = SuperClusterViewer(app=app)
        self.compound_grid = CompoundGrid(app=app)
        self.category_histogram = CategoryHistogram(app=app)
        self.compound_data_chart = CompoundDataChart(app=app)
        self.compound_table = CompoundTable(app=app)
        
        # Initialize responsive layout manager
        self.responsive_layout = ResponsiveLayoutManager(app=app)
        
        # Main content layout - will be populated after clustering
        self.main_content = pn.Column(
            pn.pane.Markdown(
                "## Welcome to ClusterF\n\n"
                "1. Select a library and dataset from the sidebar\n"
                "2. Click 'Build Super Clusters' to generate visualizations\n"
                "3. Explore your compound clusters interactively",
                margin=(20, 20)
            ),
            sizing_mode="stretch_both",
            margin=0
        )
        
        # Watch for clustering completion
        if hasattr(app, 'sc_builder'):
            # Monitor when clusters are built
            app.sc_builder.param.watch(self._on_clusters_built, 'clusters_built')
    
    def _on_clusters_built(self, event):
        """Handle when super clusters are built and update the main view."""
        if event.new and self.app.library:
            self.clusters_built = True
            # Don't automatically set up visualizations - let the super cluster selector handle this
            self._setup_main_layout()
        else:
            self.clusters_built = False
            self._reset_main_view()
    
    def _setup_main_layout(self):
        """Set up the main layout structure after clustering is complete."""
        if not self.app.library or not hasattr(self.app.library, 'super_clusters'):
            return
        
        # Update responsive layout with all components
        components = {
            'cluster_graph': self.cluster_viewer.view,
            'compound_grid': self.compound_grid.view,
            'category_histogram': self.category_histogram.view,
            'compound_data_chart': self.compound_data_chart.view,
            'compound_table': self.compound_table.view
        }
        self.responsive_layout.update_components(components)
        
        # Create the main layout with responsive layout manager
        self.main_content.objects = [
            self.responsive_layout.get_layout_controls(),
            self.responsive_layout.get_current_layout()
        ]
    
    def _reset_main_view(self):
        """Reset the main view to the welcome message."""
        self.main_content.objects = [
            pn.pane.Markdown(
                "## Welcome to ClusterF\n\n"
                "1. Select a library and dataset from the sidebar\n"
                "2. Click 'Build Super Clusters' to generate visualizations\n"
                "3. Explore your compound clusters interactively",
                margin=(20, 20)
            )
        ]
    
    def _on_cluster_selection_change(self, selected_nodes):
        """Handle cluster selection changes from the cluster viewer."""
        # This method will be called by the cluster viewer when nodes are selected
        # Here we can update other views based on the selection
        print(f"Cluster nodes selected: {selected_nodes}")
        
        # Update compound table with selected clusters
        if selected_nodes:
            # Get compounds from selected clusters
            compounds = []
            for cluster in selected_nodes:
                cluster_compounds = self.app.library.df[
                    self.app.library.df["Cluster"] == cluster
                ]["Compound"].tolist()
                compounds.extend(cluster_compounds)
            
            # Update components with selected compounds
            self.compound_table.update_table(compounds=compounds)
            self.compound_grid.update_compounds(compounds)
            self.compound_data_chart.update_chart(compounds)
        else:
            # Clear selections
            self.compound_table.update_table()
            self.compound_grid.update_compounds([])
            self.compound_data_chart.update_chart([])
    
    def _on_compound_selection_change(self, selected_compounds):
        """Handle compound selection changes from the compound table."""
        # Update compound grid and data chart with selected compounds
        self.compound_grid.update_compounds(selected_compounds)
        self.compound_data_chart.update_chart(selected_compounds)
    
    def update_colors(self, color_dict):
        """Update colors across all visualizations."""
        if self.clusters_built and hasattr(self, 'cluster_viewer'):
            self.cluster_viewer.update_colors(color_dict)
            self.compound_grid.update_colors(color_dict)
            self.category_histogram.update_colors(color_dict)
            self.compound_data_chart.update_colors(color_dict)
    
    def update_cluster_view(self, super_cluster_number):
        """Update the cluster view with a specific super cluster."""
        if not self.clusters_built or not self.app.library or not hasattr(self.app.library, 'super_clusters'):
            return
        
        try:
            # Get the super cluster data (1-indexed to 0-indexed)
            super_cluster_idx = super_cluster_number - 1
            if super_cluster_idx >= len(self.app.library.super_clusters):
                return
            
            super_cluster_data = self.app.library.super_clusters[super_cluster_idx]
            cluster_nodes = super_cluster_data[2]  # [super_cluster_number, compound_count, cluster_nodes]
            
            if not cluster_nodes:
                return
            
            # Get the first member cluster from the super cluster
            member_cluster = cluster_nodes[0]
            
            # Build subgraph for this cluster
            self.app.library.build_subgraph(member_cluster)
            
            # Initialize cluster viewer with the subgraph
            if hasattr(self.app.library, 'sub_graph'):
                self.cluster_viewer.initialize_graph(
                    self.app.library.sub_graph
                )
                
                # Update colors from color picker
                if hasattr(self.app, 'color_picker') and self.app.color_picker.color_dict:
                    self.cluster_viewer.update_colors(
                        self.app.color_picker.color_dict
                    )
                
                # Update other components with the new super cluster
                self.category_histogram.update_histogram(super_cluster_number)
                self.compound_table.update_table(super_cluster=super_cluster_number)
                
                print(f"Main view updated to super cluster {super_cluster_number} with {len(cluster_nodes)} member clusters")
            
        except Exception as e:
            print(f"Error updating main view for super cluster {super_cluster_number}: {e}")
    
    @property
    def view(self):
        """Return the main content Panel component."""
        return self.main_content
    
    def get_view(self):
        """Return the main content Panel component."""
        return self.main_content
