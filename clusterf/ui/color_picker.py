from __future__ import annotations
from typing import TYPE_CHECKING, Dict, List

import panel as pn
import param

if TYPE_CHECKING:
    from clusterf.app import ClusterFApp


def _category_sort_key(category: str) -> tuple:
    """Sorts categories in a stable, meaningful order similar to clusterf.py."""
    c = str(category).lower()
    if "selective" in c:
        return (0, c)
    if "universal" in c:
        return (1, c)
    if "increaser" in c or "decreaser" in c:
        return (2, c)
    if "bidirectional" in c:
        return (3, c)
    if "interfering" in c:
        return (4, c)
    return (5, c)


class CategoryColorPicker(param.Parameterized):
    """
    Builds a panel of color pickers for dataset categories once a dataset is loaded.

    Usage (patterned after clustering.py/app.py):
        app = ClusterFApp()
        color_picker = CategoryColorPicker(app=app)
        sidebar = pn.Column(color_picker.controls, ...)

    The widget panel becomes visible once `app.library.subset_df` exists and has a
    "Category" column. It rebuilds automatically when the selected dataset changes.
    """

    app: "ClusterFApp"

    _color_widgets: Dict[str, pn.widgets.ColorPicker] = param.Dict(
        default={}, precedence=-1
    )
    color_dict: Dict[str, str] = param.Dict(default={}, precedence=-1)

    # Wong 2011 palette (same spirit as clusterf.py default_colors)
    _default_palette: List[str] = [
        "#009E73",  # Bluish Green
        "#98FB98",  # Pale Green
        "#56B4E9",  # Sky Blue
        "#F0E442",  # Yellow
        "#E69F00",  # Orange
        "#CC79A7",  # Reddish Purple
        "#D55E00",  # Vermillion
        "#FF6347",  # Tomato
        "#FFB6C1",  # Light Pink
        "#DDA0DD",  # Plum
        "#0072B2",  # Blue
        "#F0E68C",  # Khaki
    ]

    def __init__(self, app: "ClusterFApp", **params) -> None:
        super().__init__(**params)
        self.app = app

        # Initial placeholder control; will be replaced when data is available
        self.controls = pn.Card(
            pn.pane.Markdown("*No categories loaded. Load a dataset first.*"),
            title="Category Colors",
            collapsed=True,
            visible=False,
            width=220
        )

        # If the app has a dataset already loaded, build immediately
        self._maybe_build()

        # Rebuild when dataset selection changes (pattern: observe app dataset_loaded param)
        # This ensures the color picker rebuilds when a dataset is actually loaded
        if hasattr(self.app, "param") and hasattr(self.app.param, "dataset_loaded"):
            self.app.param.watch(self._on_dataset_changed, "dataset_loaded")

    def get_colors(self) -> Dict[str, str]:
        """Return the current category→color mapping."""
        return dict(self.color_dict)

    def _on_dataset_changed(self, _):
        self._maybe_build()

    def _maybe_build(self) -> None:
        library = getattr(self.app, "library", None)
        subset_df = getattr(library, "subset_df", None)
        if subset_df is None or "Category" not in getattr(subset_df, "columns", []):
            # Hide controls if dataset isn’t ready
            self.controls.visible = False
            self._color_widgets.clear()
            self.color_dict = {}
            return

        # Extract categories present in the subset (excluding NaN)
        try:
            categories = sorted(
                [c for c in subset_df["Category"].dropna().unique().tolist()],
                key=_category_sort_key,
            )
        except Exception:
            categories = []
        if not categories:
            self.controls.visible = False
            self._color_widgets.clear()
            self.color_dict = {}
            return

        # Build widgets (reuse existing colors where possible)
        widgets: List[pn.widgets.ColorPicker] = []
        new_color_dict: Dict[str, str] = {}
        for i, category in enumerate(categories):
            prev = self.color_dict.get(category)
            default = (
                prev if prev else self._default_palette[i % len(self._default_palette)]
            )
            widget = pn.widgets.ColorPicker(
                name=str(category), value=default, width=120, height=50, margin=(5, 10)
            )
            widget.param.watch(self._on_color_change, "value")
            self._color_widgets[str(category)] = widget
            new_color_dict[str(category)] = default
            widgets.append(widget)
        self.color_dict = new_color_dict

        # Layout: 2 per row for readability
        rows: List[pn.Row] = []
        for i in range(0, len(widgets), 2):
            rows.append(pn.Row(*widgets[i : i + 2]))
        body = pn.Column(*rows)

        # Update card
        self.controls.objects = [body]
        self.controls.visible = True
        
        # Trigger a color change event to notify other components of the new color dictionary
        # This ensures components like the category histogram get updated colors immediately
        self.param.trigger('color_dict')
        
        # Also notify the app that the color picker has finished rebuilding
        if hasattr(self.app, '_on_color_picker_rebuilt'):
            self.app._on_color_picker_rebuilt()

    def _on_color_change(self, _):
        # Sync dict from widget values
        self.color_dict = {name: w.value for name, w in self._color_widgets.items()}
