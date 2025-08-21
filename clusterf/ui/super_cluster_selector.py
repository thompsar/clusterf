from __future__ import annotations
from typing import TYPE_CHECKING
import panel as pn
import param
import holoviews as hv

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
        self._counts = []  # list of (super_cluster_number, total_compounds)
        self._base_scatter = None
        self._current_overlay = None
        
        # Create the discrete slider widget
        self.slider_widget = pn.widgets.DiscreteSlider.from_param(
            self.param.current_super_cluster,
            name="Super Cluster",
            disabled=True,
            width=200
        )
        
        # Small scatter chart above the slider
        self.chart_pane = pn.pane.HoloViews(
            object=None,
            height=120,
            sizing_mode="stretch_width",
            margin=(0, 0, 5, 0),
        )

        # Create update mode selector
        self.update_mode_widget = pn.widgets.RadioButtonGroup.from_param(
            self.param.update_mode,
            name="Update Mode",
            width=200
        )
        
        # Create the controls panel
        self.controls = pn.Card(
            pn.Column(
                self.chart_pane,
                self.slider_widget,
                self.update_mode_widget,
                sizing_mode="stretch_width",
                margin=0,
            ),
            title="Super Cluster Selector",
            width=220,
            collapsed=False,
            visible=False,  # Initially hidden until clusters are built
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
            
            # Show the super cluster selector card
            self.controls.visible = True
            
            # Build counts data and chart
            try:
                sc_list = self.app.library.super_clusters
                # Expected structure: [super_cluster_number, compound_count, cluster_nodes]
                self._counts = [
                    (int(sc[0]) if len(sc) > 0 else idx + 1, int(sc[1]) if len(sc) > 1 else 0)
                    for idx, sc in enumerate(sc_list)
                ]
            except Exception:
                # Fallback: infer 1..N with zero counts if structure differs
                self._counts = [(i + 1, 0) for i in range(total_clusters)]
            self._build_counts_chart()
            
            # Trigger the first super cluster view update
            self._update_cluster_view(1)
        else:
            self.slider_widget.disabled = True
            self.total_super_clusters = 1
            self.param.current_super_cluster.objects = [1]
            # Hide the super cluster selector card
            self.controls.visible = False
            self._counts = []
            if hasattr(self, "chart_pane"):
                self.chart_pane.object = None
    
    def _on_slider_change(self, event):
        """Handle real-time slider changes."""
        # Always update the highlight dot in real time
        try:
            self._update_highlight(event.new)
        except Exception:
            pass
        # Update main view only in real-time mode
        if self.update_mode == "Real-time":
            self._update_cluster_view(event.new)
    
    def _on_slider_throttled(self, event):
        """Handle throttled slider changes (when user stops moving the slider)."""
        # Ensure highlight matches the final position as well
        try:
            self._update_highlight(event.new)
        except Exception:
            pass
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
            try:
                self._update_highlight(super_cluster_number)
            except Exception:
                pass

    def reset(self):
        """Reset the selector to the initial disabled/hidden state."""
        self.total_super_clusters = 1
        self.param.current_super_cluster.objects = [1]
        self.current_super_cluster = 1
        if hasattr(self, "slider_widget"):
            self.slider_widget.disabled = True
        if hasattr(self, "controls"):
            self.controls.visible = False
        if hasattr(self, "chart_pane"):
            self.chart_pane.object = None
        self._counts = []

    # ---- Internal helpers for chart ----
    def _build_counts_chart(self):
        """Build the base scatter plot of total compounds per super cluster and add highlight."""
        if not self._counts:
            self.chart_pane.object = None
            return
        try:
            base = hv.Scatter(self._counts, kdims=["SuperCluster"], vdims=["Count"]).opts(
                color="#1f77b4",  # blue
                size=5,
                alpha=0.8,
                tools=["hover"],
                active_tools=[],
                responsive=True,
                min_height=100,
                xlabel="Super Cluster",
                ylabel='# Cmpds',
                show_grid=False,
                xlim=(-10, max(k for k, _ in self._counts)+10),
            )
            self._base_scatter = base
            # Compose with initial highlight
            self._update_highlight(self.current_super_cluster)
        except Exception as e:
            print(f"Error building super cluster counts chart: {e}")
            self.chart_pane.object = None

    def _update_highlight(self, sc_number: int):
        """Update the red highlight dot to the given super cluster number."""
        if not self._counts or not self._base_scatter:
            return
        try:
            # Find y value for the given x
            y = 0
            for x, cnt in self._counts:
                if int(x) == int(sc_number):
                    y = cnt
                    break
            highlight = hv.Scatter([(int(sc_number), y)], kdims=["SuperCluster"], vdims=["Count"]).opts(
                color="red",
                size=8,
                alpha=1.0,
            )
            overlay = self._base_scatter * highlight
            self._current_overlay = overlay
            self.chart_pane.object = overlay
        except Exception as e:
            print(f"Error updating highlight point: {e}")
