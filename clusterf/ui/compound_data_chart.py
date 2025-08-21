from __future__ import annotations
from typing import TYPE_CHECKING
import panel as pn
import param
import holoviews as hv

if TYPE_CHECKING:
    from clusterf.app import ClusterFApp


# Bokeh hook: push Bars below other glyphs so overlays like ErrorBars render on top
# see https://github.com/holoviz/holoviews/issues/1968
def _bars_to_underlay(plot, element):
    """Set the Bars glyph renderer to 'underlay' level in Bokeh to ensure
    other overlays (e.g., ErrorBars) are drawn above the bars.
    """
    try:
        renderer = plot.handles.get("glyph_renderer")
        if renderer is not None:
            # Valid Bokeh levels: 'image', 'underlay', 'glyph', 'annotation', 'overlay'
            # renderer.level = "underlay"
            renderer.level = "underlay"
    except Exception:
        # Best-effort; ignore if backend/structure differs
        pass


def _push_grid_to_bottom(plot, element):
    """Force the plot's grid renderers to the lowest 'image' level so the grid
    is always the most underlaid item.
    """
    try:
        fig = getattr(plot, "state", None)
        if fig is None:
            return
        # xgrid and ygrid are lists of Grid renderers
        for grid in list(getattr(fig, "xgrid", [])) + list(getattr(fig, "ygrid", [])):
            try:
                grid.level = "image"
            except Exception:
                continue
    except Exception:
        # Ignore if backend changes or attributes differ
        pass


class CompoundDataChart(param.Parameterized):
    """
    A component for displaying compound data visualizations.
    Shows compound data charts and analysis.
    """

    app: "ClusterFApp" = param.Parameter(
        default=None, doc="Reference to main ClusterF application"
    )
    color_dict = param.Dict(default={}, doc="Color mapping for categories")
    selected_compounds = param.List(default=[], doc="Currently selected compounds")
    show_miss_compounds = param.Boolean(
        default=True, doc="Whether to show 'Miss' compounds"
    )

    def __init__(self, app: "ClusterFApp", **params):
        super().__init__(**params)
        self.app = app

        # Create a responsive container for charts
        self.chart_pane = pn.Column(
            sizing_mode="stretch_both",
            margin=0,
        )

        # Create the main view
        self.view = pn.Card(
            self.chart_pane,
            title="Compound Data Analysis",
            collapsed=False,
            margin=(5, 5),
            sizing_mode="stretch_both",
        )

    def update_chart(self, compound_ids: list = None, metric: str = "Delta Lifetime Z"):
        """Update the chart with new compound data."""
        if compound_ids is not None:
            self.selected_compounds = compound_ids

        # Get current colors from the app's color picker if available
        if hasattr(self.app, "color_picker") and self.app.color_picker.color_dict:
            self.color_dict = self.app.color_picker.color_dict.copy()

        if not self.selected_compounds or not self.app.library:
            self.chart_pane.objects = [
                pn.pane.HoloViews(
                    hv.Text(0.5, 0.5, "No compounds selected").opts(
                        width=600, height=400, xaxis=None, yaxis=None
                    ),
                    sizing_mode="stretch_width",
                )
            ]
            return

        try:
            # Get compound data from dataset
            compound_df = self.app.library.dataset_df[
                self.app.library.dataset_df["Compound"].isin(self.selected_compounds)
            ]

            if compound_df.empty:
                self.chart_pane.objects = [
                    pn.pane.HoloViews(
                        hv.Text(
                            0.5, 0.5, "No data available for selected compounds"
                        ).opts(width=600, height=400, xaxis=None, yaxis=None),
                        sizing_mode="stretch_width",
                    )
                ]
                return

            # Filter out "Miss" compounds if show_miss_compounds is False
            if not self.show_miss_compounds:
                # Get category information from the main library df
                category_info = self.app.library.df[
                    self.app.library.df["Compound"].isin(self.selected_compounds)
                ][["Compound", "Category"]]

                # Filter out Miss compounds
                non_miss_compounds = category_info[category_info["Category"] != "Miss"][
                    "Compound"
                ].tolist()

                # Update compound_df to only include non-Miss compounds
                compound_df = compound_df[
                    compound_df["Compound"].isin(non_miss_compounds)
                ]

                if compound_df.empty:
                    self.chart_pane.objects = [
                        pn.pane.HoloViews(
                            hv.Text(0.5, 0.5, "No non-Miss compounds available").opts(
                                width=600, height=400, xaxis=None, yaxis=None
                            ),
                            sizing_mode="stretch_width",
                        )
                    ]
                    return

            # Create statistics (include count to determine if error bars are needed)
            stats_df = (
                compound_df.groupby(["Compound", "Construct"])[metric]
                .agg(["mean", "std", "count"])
                .reset_index()
                .rename(columns={"mean": "Mean", "std": "Std", "count": "N"})
            )

            # Add category information
            stats_df = stats_df.merge(
                compound_df[["Compound", "Category"]].drop_duplicates(),
                on="Compound",
                how="left",
            )

            # Add colors using color_dict mapping (matching clusterf.py)
            stats_df["Color"] = stats_df["Category"].map(self.color_dict)
            stats_df.loc[stats_df["Color"].isna(), "Color"] = (
                "#999999"  # Default color for missing categories
            )

            # Create title based on number of compounds (matching clusterf.py)
            if len(self.selected_compounds) == 1:
                title = f"Delta Lifetime Z for Compound {self.selected_compounds[0]}"
            else:
                title = f"Delta Lifetime Z for {len(self.selected_compounds)} Compounds"

            # Calculate symmetrical y-axis limits based on data including std
            # Use max(|mean| + std) to ensure error bars fit
            max_abs_value = 0.0
            if not stats_df.empty:
                max_abs_value = float(
                    (stats_df["Mean"].abs() + stats_df["Std"].fillna(0)).max()
                )
            # Add some padding and ensure we can see ±4 lines
            y_limit = max(max_abs_value * 1.1, 4.5)

            # Build one bar chart per Construct (Compound on x-axis), fill available width
            chart_panes = []
            constructs = list(stats_df["Construct"].drop_duplicates())
            for construct in constructs:
                construct_df = stats_df[stats_df["Construct"] == construct]

                bars = hv.Bars(
                    construct_df,
                    kdims=["Compound"],
                    vdims=["Mean", "Category", "Color"],
                ).opts(
                    color="Color",
                    responsive=True,
                    width=600,
                    min_height=250,
                    title=f"{title} – Construct: {construct}",
                    ylabel=metric,
                    xlabel="Compound",
                    ylim=(-y_limit, y_limit),
                    xrotation=45,
                    tools=["hover"],
                    show_grid=True,
                    show_legend=True,
                    legend_position="right",
                    hooks=[_bars_to_underlay, _push_grid_to_bottom],
                )

                # Add symmetric error bars (Mean ± Std) only when N>1 and Std is present
                error_df = construct_df[["Compound", "Mean", "Std", "N"]].copy()
                error_df = error_df[error_df["N"] > 1]
                error_df = error_df[error_df["Std"].notna() & (error_df["Std"] > 0)]

                # Add horizontal reference lines at ±4
                hline_pos4 = hv.HLine(4).opts(
                    color="red", line_dash="dashed", line_width=2, alpha=0.7
                )
                hline_neg4 = hv.HLine(-4).opts(
                    color="red", line_dash="dashed", line_width=2, alpha=0.7
                )
                overlay = bars * hline_pos4 * hline_neg4
                if not error_df.empty:
                    errorbars = hv.ErrorBars(
                        error_df,
                        kdims=["Compound"],
                        vdims=["Mean", "Std"],
                    ).opts(
                        color="black",
                        line_width=2,
                        line_alpha=0.9,
                    )
                    # Raw replicate points for compounds with N>1 (same condition as error bars)
                    raw_df = self.app.library.dataset_df[
                        (self.app.library.dataset_df["Compound"].isin(self.selected_compounds))
                        & (self.app.library.dataset_df["Construct"] == construct)
                    ][["Compound", metric]].copy()
                    # Limit raw points to compounds that have error bars (N>1)
                    raw_df = raw_df[raw_df["Compound"].isin(error_df["Compound"])]

                    points = hv.Scatter(
                        raw_df,
                        kdims=["Compound"],
                        vdims=[metric],
                    ).opts(
                        color="black",
                        size=6,
                        alpha=0.6,
                    )

                    overlay = overlay * errorbars * points

                chart_panes.append(
                    pn.pane.HoloViews(
                        overlay,
                        sizing_mode="stretch_width",
                        min_height=300,
                        margin=0,
                    )
                )

            # Populate the column with one responsive pane per construct
            self.chart_pane.objects = chart_panes

        except Exception as e:
            print(f"Error creating compound data chart: {e}")
            self.chart_pane.objects = [
                pn.pane.HoloViews(
                    hv.Text(0.5, 0.5, f"Error: {str(e)}").opts(
                        width=600, height=400, xaxis=None, yaxis=None
                    ),
                    sizing_mode="stretch_width",
                )
            ]

    def update_colors(self, color_dict: dict):
        """Update the color scheme and refresh the chart."""
        self.color_dict = color_dict.copy()
        # Force refresh the chart with new colors
        self.update_chart()

    def refresh_chart(self):
        """Force refresh the chart with current data and colors."""
        self.update_chart()

    def get_show_miss_compounds(self) -> bool:
        """Get the current state of show_miss_compounds."""
        return self.show_miss_compounds

    def set_show_miss_compounds(self, show_miss_compounds: bool):
        """Set the show_miss_compounds state and update the chart if needed."""
        if self.show_miss_compounds != show_miss_compounds:
            self.show_miss_compounds = show_miss_compounds
            # Update the chart with current compounds to apply the new filter
            if self.selected_compounds:
                self.update_chart()

    def reset(self):
        """Reset the chart to the initial state."""
        self.color_dict = {}
        self.selected_compounds = []
        self.show_miss_compounds = True
        if hasattr(self, "chart_pane"):
            self.chart_pane.objects = []
