from __future__ import annotations
from typing import TYPE_CHECKING
import os
import panel as pn
import param
from clusterf.ui.cluster_builder import SuperClusterBuilder
from clusterf.ui.color_picker import CategoryColorPicker
from clusterf.ui.main_view_manager import MainViewManager
from clusterf.ui.super_cluster_selector import SuperClusterSelector
from clusterf.ui.dataset_stats import DatasetStatsCard

if TYPE_CHECKING:
    from clusterf.core.chemistry import ChemLibrary


class ClusterFApp(param.Parameterized):
    library: ChemLibrary | None = None
    LIBRARIES_DIR = os.path.join("data", "libraries")
    DATASETS_DIR = os.path.join("data", "datasets")

    # Add a parameter to track dataset changes
    dataset_loaded = param.Boolean(
        default=False, doc="Whether a dataset has been loaded"
    )

    def __init__(self, **params):
        super().__init__(**params)
        pn.extension("tabulator")

        self.sc_builder = SuperClusterBuilder(app=self)
        self.color_picker = CategoryColorPicker(app=self)
        self.main_view = MainViewManager(app=self)
        self.super_cluster_selector = SuperClusterSelector(app=self)
        self.dataset_stats = DatasetStatsCard(app=self)

        self.sidebar = pn.Column(
            self.sc_builder.controls,
            self.color_picker.controls,
            self.super_cluster_selector.controls,
            self.dataset_stats.view,
            width=250,
            margin=(5, 5),
            styles={"padding": "10px", "border-right": "1px solid #ddd"},
        )

        # Watch for color changes and propagate to main view
        self.color_picker.param.watch(self._on_color_change, "color_dict")
    # Also watch dataset_loaded to propagate resets if needed (handled via _on_dataset_loaded)

    def _on_color_change(self, event):
        """Handle color changes from the color picker."""
        if event.new and self.main_view:
            self.main_view.update_colors(event.new)

    def _on_cluster_selection_change(self, selected_nodes):
        """Handle cluster selection changes from visualizations."""
        # This method can be called by visualization components
        # to notify the app of selection changes
        if hasattr(self.main_view, "_on_cluster_selection_change"):
            self.main_view._on_cluster_selection_change(selected_nodes)

    def _on_compound_selection_change(self, selected_compounds):
        """Handle compound selection changes from the compound table."""
        # This method can be called by the compound table
        # to notify the app of compound selection changes
        if hasattr(self.main_view, "_on_compound_selection_change"):
            self.main_view._on_compound_selection_change(selected_compounds)

    def _on_miss_compounds_toggle(self, show_miss_compounds):
        """Handle miss compounds toggle changes from the compound table."""
        # This method can be called by the compound table
        # to notify the app of miss compounds toggle changes
        if hasattr(self.main_view, "_on_miss_compounds_toggle"):
            self.main_view._on_miss_compounds_toggle(show_miss_compounds)

    def _on_dataset_loaded(self):
        """Handle dataset loading completion."""
        # Reset app UI to launch state when a new dataset is loaded
        self.reset_app_state()
        # Toggle the dataset_loaded parameter to trigger watchers (color picker, etc.)
        self.dataset_loaded = not self.dataset_loaded

        # Notify main view of dataset change (if main view exists)
        if hasattr(self, "main_view") and hasattr(self.main_view, "_on_dataset_changed"):
            self.main_view._on_dataset_changed()

    def reset_app_state(self):
        """Reset UI components and clustering state to initial launch state."""
        # Reset clustering state on the library if present
        if hasattr(self, "library") and self.library and hasattr(self.library, "reset_clustering_state"):
            self.library.reset_clustering_state()
        # Reset main view and children
        if hasattr(self, "main_view") and hasattr(self.main_view, "reset"):
            self.main_view.reset()
        # Reset super cluster selector
        if hasattr(self, "super_cluster_selector") and hasattr(self.super_cluster_selector, "reset"):
            self.super_cluster_selector.reset()
        # Reset color picker UI; it will rebuild after dataset_loaded toggles
        if hasattr(self, "color_picker") and hasattr(self.color_picker, "reset"):
            self.color_picker.reset()

    def _on_color_picker_rebuilt(self):
        """Handle color picker rebuild completion."""
        # Now that the color picker has finished rebuilding, update all components
        if hasattr(self, "main_view"):
            # Update all components with the new colors
            if hasattr(self.main_view, "category_histogram"):
                self.main_view.category_histogram.refresh_histogram()
            if hasattr(self.main_view, "compound_grid"):
                self.main_view.compound_grid.update_colors(self.color_picker.color_dict)
            if hasattr(self.main_view, "compound_data_chart"):
                self.main_view.compound_data_chart.update_colors(
                    self.color_picker.color_dict
                )
            if hasattr(self.main_view, "cluster_viewer"):
                self.main_view.cluster_viewer.update_colors(
                    self.color_picker.color_dict
                )
            if hasattr(self.main_view, "compound_table"):
                self.main_view.compound_table.update_colors(
                    self.color_picker.color_dict
                )
            # Refresh stats in case category names/colors changed
            if hasattr(self, "dataset_stats"):
                self.dataset_stats.refresh()

    def serve(self):
        return pn.Row(
            self.sidebar,
            pn.Spacer(width=10),  # Add spacing between sidebar and main view
            self.main_view.view,
            sizing_mode="stretch_both",
            margin=0,
        )


if __name__ == "__main__":
    pn.serve(ClusterFApp().serve(), show=True)
