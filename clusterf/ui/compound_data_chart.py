from __future__ import annotations
from typing import TYPE_CHECKING
import panel as pn
import param
import holoviews as hv
import pandas as pd

if TYPE_CHECKING:
    from clusterf.app import ClusterFApp


class CompoundDataChart(param.Parameterized):
    """
    A component for displaying compound data visualizations.
    Shows compound data charts and analysis.
    """
    
    app: "ClusterFApp" = param.Parameter(default=None, doc="Reference to main ClusterF application")
    color_dict = param.Dict(default={}, doc="Color mapping for categories")
    selected_compounds = param.List(default=[], doc="Currently selected compounds")
    
    def __init__(self, app: "ClusterFApp", **params):
        super().__init__(**params)
        self.app = app
        
        # Create the chart pane
        self.chart_pane = pn.pane.HoloViews(
            object=None,
            sizing_mode="stretch_both",
            margin=0
        )
        
        # Create the main view
        self.view = pn.Card(
            self.chart_pane,
            title="Compound Data Analysis",
            collapsed=False,
            margin=(5, 5)
        )
    
    def update_chart(self, compound_ids: list = None, metric: str = "Delta Lifetime Z"):
        """Update the chart with new compound data."""
        if compound_ids is not None:
            self.selected_compounds = compound_ids
        
        # Get current colors from the app's color picker if available
        if hasattr(self.app, 'color_picker') and self.app.color_picker.color_dict:
            self.color_dict = self.app.color_picker.color_dict.copy()
        
        if not self.selected_compounds or not self.app.library:
            self.chart_pane.object = hv.Text(0.5, 0.5, "No compounds selected").opts(
                width=600, height=400, xaxis=None, yaxis=None
            )
            return
        
        try:
            # Get compound data from dataset
            compound_df = self.app.library.dataset_df[
                self.app.library.dataset_df["Compound"].isin(self.selected_compounds)
            ]
            
            if compound_df.empty:
                self.chart_pane.object = hv.Text(0.5, 0.5, "No data available for selected compounds").opts(
                    width=600, height=400, xaxis=None, yaxis=None
                )
                return
            
            # Create statistics
            stats_df = (
                compound_df.groupby(["Compound", "Construct"])[metric]
                .agg(["mean", "std"])
                .reset_index()
                .rename(columns={"mean": "Mean", "std": "Std"})
            )
            
            # Add category information
            stats_df = stats_df.merge(
                compound_df[["Compound", "Category"]].drop_duplicates(),
                on="Compound",
                how="left"
            )
            
            # Add colors using color_dict mapping (matching clusterf.py)
            stats_df["Color"] = stats_df["Category"].map(self.color_dict)
            stats_df.loc[stats_df["Color"].isna(), "Color"] = "#999999"  # Default color for missing categories
            
            # Create title based on number of compounds (matching clusterf.py)
            if len(self.selected_compounds) == 1:
                title = f"Delta Lifetime Z for Compound {self.selected_compounds[0]}"
            else:
                title = f"Delta Lifetime Z for {len(self.selected_compounds)} Compounds"
            
            # Calculate symmetrical y-axis limits based on data (matching clusterf.py)
            max_abs_value = max(abs(stats_df["Mean"].min()), abs(stats_df["Mean"].max()))
            # Add some padding and ensure we can see ±4 lines
            y_limit = max(max_abs_value * 1.1, 4.5)
            
            # Create bar chart (matching clusterf.py dimensions and options)
            bars = hv.Bars(
                stats_df,
                kdims=["Construct", "Compound"],
                vdims=["Mean", "Category", "Color"]
            ).opts(
                color="Color",
                width=600,
                height=400,
                title=title,
                ylabel=metric,
                xlabel="Construct",
                ylim=(-y_limit, y_limit),
                xrotation=45,
                tools=["hover"],
                show_grid=True,
                show_legend=True,
                legend_position="right"
            )
            
            # Add horizontal reference lines at ±4 (matching clusterf.py)
            hline_pos4 = hv.HLine(4).opts(
                color="red", line_dash="dashed", line_width=2, alpha=0.7
            )
            hline_neg4 = hv.HLine(-4).opts(
                color="red", line_dash="dashed", line_width=2, alpha=0.7
            )
            
            # Combine bars with reference lines
            chart = bars * hline_pos4 * hline_neg4
            self.chart_pane.object = chart
            
        except Exception as e:
            print(f"Error creating compound data chart: {e}")
            self.chart_pane.object = hv.Text(0.5, 0.5, f"Error: {str(e)}").opts(
                width=600, height=400, xaxis=None, yaxis=None
            )
    
    def update_colors(self, color_dict: dict):
        """Update the color scheme and refresh the chart."""
        self.color_dict = color_dict.copy()
        # Force refresh the chart with new colors
        self.update_chart()
        
    def refresh_chart(self):
        """Force refresh the chart with current data and colors."""
        self.update_chart()
