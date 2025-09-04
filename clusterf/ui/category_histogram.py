from __future__ import annotations
from typing import TYPE_CHECKING
import panel as pn
import param
import holoviews as hv
import pandas as pd

if TYPE_CHECKING:
    from clusterf.app import ClusterFApp


class CategoryHistogram(param.Parameterized):
    """
    A component for displaying category distribution histograms.
    Shows the distribution of compound categories in the current super cluster.
    """

    app: "ClusterFApp" = param.Parameter(
        default=None, doc="Reference to main ClusterF application"
    )
    color_dict = param.Dict(default={}, doc="Color mapping for categories")
    current_super_cluster = param.Integer(default=1, doc="Current super cluster number")

    def __init__(self, app: "ClusterFApp", **params):
        super().__init__(**params)
        self.app = app

        # Create the histogram pane
        self.histogram_pane = pn.pane.HoloViews(
            object=None, sizing_mode="stretch_both", margin=0
        )

        # Create the main view
        self.view = pn.Card(
            self.histogram_pane,
            title="Category Distribution",
            collapsed=False,
            margin=(5, 5),
            sizing_mode="stretch_both",
        )

    def update_histogram(self, super_cluster_number: int = None):
        """Update the histogram for a specific super cluster."""
        if super_cluster_number is not None:
            self.current_super_cluster = super_cluster_number

        # Get current colors from the app's color picker if available
        if hasattr(self.app, "color_picker") and self.app.color_picker.color_dict:
            self.color_dict = self.app.color_picker.color_dict.copy()

        if not self.app.library or not hasattr(self.app.library, "subset_df"):
            self.histogram_pane.object = hv.Text(
                0.5, 0.5, "No subset data available for histogram."
            ).opts(width=600, height=250, xaxis=None, yaxis=None)
            return

        try:
            # Get compounds in current super cluster
            cluster_compounds = self._get_current_super_cluster_compounds()

            if not cluster_compounds:
                self.histogram_pane.object = hv.Text(
                    0.5, 0.5, "No compounds in super cluster"
                ).opts(width=600, height=250, xaxis=None, yaxis=None)
                return

            # Get categories for these compounds
            subset_data = self.app.library.subset_df[
                self.app.library.subset_df["Compound"].isin(cluster_compounds)
            ]

            if subset_data.empty:
                self.histogram_pane.object = hv.Text(
                    0.5, 0.5, "No category data for current super cluster."
                ).opts(width=600, height=250, xaxis=None, yaxis=None)
                return

            # Get all available categories from color picker to ensure completeness
            all_categories = list(self.color_dict.keys()) if self.color_dict else []

            # Sort categories using the same logic as clusterf.py
            def category_sort_key(category):
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

            all_categories = sorted(all_categories, key=category_sort_key)

            # Count categories in current super cluster (excluding "Miss")
            category_counts = subset_data[subset_data["Category"] != "Miss"][
                "Category"
            ].value_counts()

            # Create complete DataFrame including all categories (even with 0 counts)
            hist_data = pd.DataFrame(
                {
                    "Category": all_categories,
                    "Count": [category_counts.get(cat, 0) for cat in all_categories],
                    "%Total": [
                        f"{int(100 * category_counts.get(cat, 0) / len(cluster_compounds))}%"
                        for cat in all_categories
                    ],
                }
            )

            # Map colors from color_dict
            colors = []
            for category in hist_data["Category"]:
                colors.append(self.color_dict.get(category, "#999999"))

            hist_data["Color"] = colors

            # Create HoloViews bar chart
            bars = hv.Bars(hist_data, ["Category"], ["Count", "%Total", "Color"]).opts(
                color="Color",
                min_height=300,
                responsive=True,
                title=f"Category Distribution - Super Cluster {self.current_super_cluster}",
                ylabel="Count",
                xlabel="",
                xrotation=45,
                tools=["hover"],
                active_tools=[],
                shared_axes=False,
                axiswise=True,
                framewise=True,
                ylim=(
                    0,
                    1.2 * hist_data["Count"].max()
                    if hist_data["Count"].max() > 0
                    else 1,
                ),
            )

            self.histogram_pane.object = bars

        except Exception as e:
            print(f"Error creating category histogram: {e}")
            self.histogram_pane.object = hv.Text(0.5, 0.5, f"Error: {str(e)}").opts(
                width=600, height=250, xaxis=None, yaxis=None
            )

    def update_colors(self, color_dict: dict):
        """Update the color scheme and refresh the histogram."""
        self.color_dict = color_dict.copy()
        # Force refresh the histogram with new colors
        self.update_histogram()

    def refresh_histogram(self):
        """Force refresh the histogram with current data and colors."""
        # Get updated colors from the app's color picker if available
        if hasattr(self.app, "color_picker") and self.app.color_picker.color_dict:
            self.color_dict = self.app.color_picker.color_dict.copy()
        else:
            # Clear old color dict if no new colors available
            self.color_dict = {}

        # If no categories available, show appropriate message
        if not self.color_dict:
            self.histogram_pane.object = hv.Text(
                0.5, 0.5, "Load a dataset to see category distribution."
            ).opts(width=600, height=250, xaxis=None, yaxis=None)
            return

        self.update_histogram()

    def reset(self):
        """Reset the histogram to the initial state."""
        self.color_dict = {}
        self.current_super_cluster = 1
        if hasattr(self, "histogram_pane"):
            self.histogram_pane.object = None

    def _get_current_super_cluster_compounds(self):
        """Get list of compounds in the current super cluster."""
        if not self.app.library or not hasattr(self.app.library, "df"):
            return []

        # Check if SuperCluster column exists
        if "SuperCluster" in self.app.library.df.columns:
            return self.app.library.df[
                self.app.library.df["SuperCluster"] == self.current_super_cluster
            ]["Compound"].tolist()
        else:
            # Fallback to old method if SuperCluster column doesn't exist
            super_cluster_idx = self.current_super_cluster - 1
            if hasattr(self.app.library, "super_clusters") and super_cluster_idx < len(
                self.app.library.super_clusters
            ):
                cluster_list = self.app.library.super_clusters[super_cluster_idx][2]
                return self.app.library.df[
                    self.app.library.df["Cluster"].isin(cluster_list)
                ]["Compound"].tolist()
            else:
                return []
