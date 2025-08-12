from __future__ import annotations
from typing import TYPE_CHECKING
import panel as pn
import param
from clusterf.ui.cluster_viewer import SuperClusterViewer

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
            self._setup_main_visualizations()
        else:
            self.clusters_built = False
            self._reset_main_view()
    
    def _setup_main_visualizations(self):
        """Set up the main visualization layout after clustering is complete."""
        if not self.app.library or not hasattr(self.app.library, 'super_clusters'):
            return
            
        # Get the first super cluster and build subgraph for it
        if not self.app.library.super_clusters:
            return
            
        # Get the first cluster node from the first super cluster (following original clusterf.py pattern)
        first_super_cluster = self.app.library.super_clusters[0]
        cluster_nodes = first_super_cluster[2]  # [super_cluster_number, compound_count, cluster_nodes]
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
        
        # Create the main layout with cluster viewer
        self.main_content.objects = [
            pn.Card(
                self.cluster_viewer.view,
                title="Cluster Network Graph",
                collapsed=False,
                margin=(5, 5)
            ),
            # Placeholder for additional visualizations
            pn.Card(
                pn.pane.Markdown(
                    "Additional visualizations (compound table, histograms, etc.) "
                    "will be added here."
                ),
                title="Additional Views",
                collapsed=True,
                margin=(5, 5)
            )
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
        
        # Future: Update compound table, histograms, etc. based on selection
        # if selected_nodes:
        #     self._update_compound_table(selected_nodes)
        #     self._update_category_histogram(selected_nodes)
    
    def update_colors(self, color_dict):
        """Update colors across all visualizations."""
        if self.clusters_built and hasattr(self, 'cluster_viewer'):
            self.cluster_viewer.update_colors(color_dict)
    
    @property
    def view(self):
        """Return the main content Panel component."""
        return self.main_content
    
    def get_view(self):
        """Return the main content Panel component."""
        return self.main_content
