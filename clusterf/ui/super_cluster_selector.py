from __future__ import annotations
from typing import TYPE_CHECKING
import panel as pn
import holoviews as hv
from holoviews import streams
import param

if TYPE_CHECKING:
    from clusterf.app import ClusterFApp


class SuperClusterSelector(param.Parameterized):
    """
    A discrete slider component for selecting super clusters with real-time and throttled update modes.
    """

    app: "ClusterFApp" = param.Parameter(
        default=None, doc="Reference to main ClusterF application"
    )
    current_super_cluster = param.Selector(
        objects=[1], default=1, doc="Currently selected super cluster"
    )
    total_super_clusters = param.Integer(
        default=1, doc="Total number of super clusters available"
    )

    def __init__(self, app: "ClusterFApp", **params):
        super().__init__(**params)
        self.app = app
        self._counts = []
        self._base_scatter = None

        # Create the discrete slider widget
        self.slider_widget = pn.widgets.DiscreteSlider.from_param(
            self.param.current_super_cluster,
            name="Super Cluster",
            disabled=True,
            width=200,
        )

        self.chart_pane = pn.pane.HoloViews(
            object=None, height=200, sizing_mode="stretch_width", margin=0
        )

        # Create the controls panel
        self.controls = pn.Card(
            self.chart_pane,
            self.slider_widget,
            title="Super Cluster Selector",
            width=220,
            collapsed=False,
            visible=False,
        )

        # Set up watchers for different update modes
        self.slider_widget.param.watch(self._on_slider_change, "value")
        self.slider_widget.param.watch(
            self._on_slider_change_throttled, "value_throttled"
        )

        # Watch for when clusters are built to enable the slider
        if hasattr(app, "sc_builder"):
            app.sc_builder.param.watch(self._on_clusters_built, "clusters_built")

    def _on_clusters_built(self, event):
        """Handle when super clusters are built and update the slider."""
        if (
            event.new
            and self.app.library
            and hasattr(self.app.library, "super_clusters")
        ):
            total_clusters = len(self.app.library.super_clusters)
            self.total_super_clusters = total_clusters

            # Update the selector options to include all super cluster numbers
            cluster_options = list(range(1, total_clusters + 1))
            self.param.current_super_cluster.objects = cluster_options

            self.slider_widget.disabled = False
            self.current_super_cluster = 1  # Start with first super cluster

            # Show the super cluster selector card
            self.controls.visible = True

            self._counts = [
                (super_clust_id, counts)
                for super_clust_id, counts, __ in self.app.library.super_clusters
            ]
            # Fast lookup for highlight
            self._count_by_id = {k: v for k, v in self._counts}

            # Build chart once
            self._build_counts_chart()

            # Initialize highlight at current slider
            self._update_highlight(self.current_super_cluster)

            # Trigger the first super cluster view update
            self._update_cluster_view(self.current_super_cluster)
        else:
            self.slider_widget.disabled = True
            self.total_super_clusters = 1
            self.param.current_super_cluster.objects = [1]
            self.controls.visible = False
            self._counts = []
            if hasattr(self, "chart_pane"):
                self.chart_pane.object = None

    def _on_slider_change(self, event):
        """Handle real-time slider changes."""
        self._update_highlight(event.new)

    def _on_slider_change_throttled(self, event):
        """Handle throttled slider changes (when user stops moving the slider)."""
        # self._update_highlight(event.new)
        self._update_cluster_view(event.new)

    def _update_cluster_view(self, super_cluster_number):
        """Update the cluster view with the selected super cluster."""
        if hasattr(self.app, "main_view"):
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

    def _build_counts_chart(self):
        """Build the base scatter plot of total compounds per super cluster and add highlight."""
        if not self._counts:
            self.chart_pane.object = None
            return
        try:
            # Small hook to ensure logo is removed even if a toolbar sneaks in
            def _hide_toolbar_logo(plot, element):
                plot.state.toolbar_location = None
                tb = getattr(plot.state, "toolbar", None)
                if tb is not None:
                    tb.logo = None

            base = hv.Scatter(
                self._counts, kdims=["SuperCluster"], vdims=["Count"]
            ).opts(
                color="#1f77b4",  # blue
                size=5,
                alpha=0.8,
                tools=["hover"],
                active_tools=[],
                responsive=True,
                height=200,
                min_height=100,
                xlabel="Super Cluster",
                ylabel="# Cmpds",
                show_grid=False,
                xlim=(-10, max(k for k, _ in self._counts) + 10),
            )
            self._base_scatter = base

            # Streamed highlight
            self._highlight_stream = streams.Pipe(data=[])

            def _highlight(data):
                # data is a list of (x, y) pairs
                return hv.Scatter(data, kdims=["SuperCluster"], vdims=["Count"]).opts(
                    color="red", size=10
                )

            highlight_dm = hv.DynamicMap(_highlight, streams=[self._highlight_stream])

            # Compose once
            overlay = (base * highlight_dm).opts(
                toolbar=None,
                hooks=[_hide_toolbar_logo],
                shared_axes=False,
                axiswise=True,
                framewise=True,
            )
            self.chart_pane.object = overlay

        except Exception as e:
            print(f"Error building super cluster counts chart: {e}")
            self.chart_pane.object = None

    def _update_highlight(self, sc_number: int):
        """Update only the highlight via stream to avoid re-render (preserves slider focus)."""
        try:
            if not getattr(self, "_highlight_stream", None):
                return
            sc_id = int(sc_number)
            # Prefer dict lookup; fall back to list index if needed
            if hasattr(self, "_count_by_id") and sc_id in self._count_by_id:
                y = self._count_by_id[sc_id]
            else:
                # Fallback for older state
                idx = sc_id - 1
                if not self._counts or idx < 0 or idx >= len(self._counts):
                    self._highlight_stream.send([])
                    return
                y = self._counts[idx][1]

            self._highlight_stream.send([(sc_id, y)])
        except Exception as e:
            print(f"Error updating highlight: {e}")
