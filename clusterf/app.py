from __future__ import annotations
from typing import TYPE_CHECKING
import os
import panel as pn
import param
from clusterf.ui.cluster_builder import SuperClusterBuilder
from clusterf.ui.color_picker import CategoryColorPicker
from clusterf.ui.main_view_manager import MainViewManager

if TYPE_CHECKING:
    from clusterf.core.chemistry import ChemLibrary


class ClusterFApp(param.Parameterized):
    library: ChemLibrary | None = None
    LIBRARIES_DIR = os.path.join("data", "libraries")
    DATASETS_DIR = os.path.join("data", "datasets")

    def __init__(self, **params):
        super().__init__(**params)
        pn.extension("tabulator")

        self.sc_builder = SuperClusterBuilder(app=self)
        self.color_picker = CategoryColorPicker(app=self)
        self.main_view = MainViewManager(app=self)

        self.sidebar = pn.Column(
            self.sc_builder.controls, 
            self.color_picker.controls, 
            width=200,
            margin=(5, 5)
        )
        
        # Watch for color changes and propagate to main view
        self.color_picker.param.watch(self._on_color_change, 'color_dict')

    def _on_color_change(self, event):
        """Handle color changes from the color picker."""
        if event.new and self.main_view:
            self.main_view.update_colors(event.new)
    
    def _on_cluster_selection_change(self, selected_nodes):
        """Handle cluster selection changes from visualizations."""
        # This method can be called by visualization components
        # to notify the app of selection changes
        if hasattr(self.main_view, '_on_cluster_selection_change'):
            self.main_view._on_cluster_selection_change(selected_nodes)

    def serve(self):
        return pn.Row(
            self.sidebar,
            self.main_view.view,
            sizing_mode="stretch_both"
        )


if __name__ == "__main__":
    pn.serve(ClusterFApp().serve(), show=True)
