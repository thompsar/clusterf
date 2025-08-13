from __future__ import annotations
from typing import TYPE_CHECKING
import panel as pn
import param

if TYPE_CHECKING:
    from clusterf.app import ClusterFApp


class SuperClusterSelector(param.Parameterized):
    """
    A discrete slider component for selecting super clusters with real-time and throttled update modes.
    """
    
    app: "ClusterFApp" = param.Parameter(default=None, doc="Reference to main ClusterF application")
    current_super_cluster = param.Selector(objects=[1], default=1, doc="Currently selected super cluster")
    total_super_clusters = param.Integer(default=1, doc="Total number of super clusters available")
    update_mode = param.Selector(objects=["Real-time", "Throttled"], default="Throttled", doc="Update mode for the slider")
    
    def __init__(self, app: "ClusterFApp", **params):
        super().__init__(**params)
        self.app = app
        
        # Create the discrete slider widget
        self.slider_widget = pn.widgets.DiscreteSlider.from_param(
            self.param.current_super_cluster,
            name="Super Cluster",
            disabled=True,
            width=200
        )
        
        # Create update mode selector
        self.update_mode_widget = pn.widgets.RadioButtonGroup.from_param(
            self.param.update_mode,
            name="Update Mode",
            width=200
        )
        
        # Create the controls panel
        self.controls = pn.Card(
            self.slider_widget,
            self.update_mode_widget,
            title="Super Cluster Selector",
            width=220,
            collapsed=False,
            margin=(5, 5)
        )
        
        # Set up watchers for different update modes
        self.slider_widget.param.watch(self._on_slider_change, "value")
        self.slider_widget.param.watch(self._on_slider_throttled, "value_throttled")
        
        # Watch for when clusters are built to enable the slider
        if hasattr(app, 'sc_builder'):
            app.sc_builder.param.watch(self._on_clusters_built, 'clusters_built')
    
    def _on_clusters_built(self, event):
        """Handle when super clusters are built and update the slider."""
        if event.new and self.app.library and hasattr(self.app.library, 'super_clusters'):
            total_clusters = len(self.app.library.super_clusters)
            self.total_super_clusters = total_clusters
            
            # Update the selector options to include all super cluster numbers
            cluster_options = list(range(1, total_clusters + 1))
            self.param.current_super_cluster.objects = cluster_options
            
            self.slider_widget.disabled = False
            self.current_super_cluster = 1  # Start with first super cluster
            
            # Trigger the first super cluster view update
            self._update_cluster_view(1)
        else:
            self.slider_widget.disabled = True
            self.total_super_clusters = 1
            self.param.current_super_cluster.objects = [1]
    
    def _on_slider_change(self, event):
        """Handle real-time slider changes."""
        if self.update_mode == "Real-time":
            self._update_cluster_view(event.new)
    
    def _on_slider_throttled(self, event):
        """Handle throttled slider changes (when user stops moving the slider)."""
        if self.update_mode == "Throttled":
            self._update_cluster_view(event.new)
    
    def _update_cluster_view(self, super_cluster_number):
        """Update the cluster view with the selected super cluster."""
        if hasattr(self.app, 'main_view'):
            self.app.main_view.update_cluster_view(super_cluster_number)
    
    def set_super_cluster(self, super_cluster_number):
        """Programmatically set the super cluster selection."""
        if super_cluster_number in self.param.current_super_cluster.objects:
            self.current_super_cluster = super_cluster_number
