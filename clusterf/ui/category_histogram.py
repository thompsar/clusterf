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
    
    app: "ClusterFApp" = param.Parameter(default=None, doc="Reference to main ClusterF application")
    color_dict = param.Dict(default={}, doc="Color mapping for categories")
    current_super_cluster = param.Integer(default=1, doc="Current super cluster number")
    
    def __init__(self, app: "ClusterFApp", **params):
        super().__init__(**params)
        self.app = app
        
        # Create the histogram pane
        self.histogram_pane = pn.pane.HoloViews(
            object=None,
            sizing_mode="stretch_both",
            margin=0
        )
        
        # Create the main view
        self.view = pn.Card(
            self.histogram_pane,
            title="Category Distribution",
            collapsed=False,
            margin=(5, 5)
        )
    
    def update_histogram(self, super_cluster_number: int = None):
        """Update the histogram for a specific super cluster."""
        if super_cluster_number is not None:
            self.current_super_cluster = super_cluster_number
        
        if not self.app.library or not hasattr(self.app.library, 'super_clusters'):
            self.histogram_pane.object = hv.Text(0.5, 0.5, "No data available").opts(
                width=400, height=300, xaxis=None, yaxis=None
            )
            return
        
        try:
            # Get compounds in current super cluster
            cluster_compounds = self._get_current_super_cluster_compounds()
            
            if not cluster_compounds:
                self.histogram_pane.object = hv.Text(0.5, 0.5, "No compounds in super cluster").opts(
                    width=400, height=300, xaxis=None, yaxis=None
                )
                return
            
            # Get categories for these compounds
            subset_data = self.app.library.subset_df[
                self.app.library.subset_df["Compound"].isin(cluster_compounds)
            ]
            
            if subset_data.empty:
                self.histogram_pane.object = hv.Text(0.5, 0.5, "No category data available").opts(
                    width=400, height=300, xaxis=None, yaxis=None
                )
                return
            
            # Count categories (excluding "Miss")
            category_counts = subset_data[subset_data["Category"] != "Miss"]["Category"].value_counts()
            
            if category_counts.empty:
                self.histogram_pane.object = hv.Text(0.5, 0.5, "No valid categories found").opts(
                    width=400, height=300, xaxis=None, yaxis=None
                )
                return
            
            # Create histogram data
            hist_data = pd.DataFrame({
                "Category": category_counts.index,
                "Count": category_counts.values,
                "Percentage": [f"{int(100 * count / len(cluster_compounds))}%" 
                             for count in category_counts.values]
            })
            
            # Add colors
            hist_data["Color"] = [self.color_dict.get(cat, "#999999") for cat in hist_data["Category"]]
            
            # Create HoloViews bar chart
            bars = hv.Bars(hist_data, ["Category"], ["Count", "Percentage", "Color"]).opts(
                color="Color",
                min_width=300,
                min_height=200,
                responsive=True,
                title=f"Category Distribution - Super Cluster {self.current_super_cluster}",
                ylabel="Count",
                xlabel="",
                xrotation=45,
                tools=["hover"],
                ylim=(0, 1.2 * hist_data["Count"].max())
            )
            
            self.histogram_pane.object = bars
            
        except Exception as e:
            print(f"Error creating category histogram: {e}")
            self.histogram_pane.object = hv.Text(0.5, 0.5, f"Error: {str(e)}").opts(
                width=400, height=300, xaxis=None, yaxis=None
            )
    
    def update_colors(self, color_dict: dict):
        """Update the color scheme and refresh the histogram."""
        self.color_dict = color_dict.copy()
        self.update_histogram()
    
    def _get_current_super_cluster_compounds(self):
        """Get list of compounds in the current super cluster."""
        if not self.app.library or not hasattr(self.app.library, 'df'):
            return []
        
        # Check if SuperCluster column exists
        if "SuperCluster" in self.app.library.df.columns:
            return self.app.library.df[
                self.app.library.df["SuperCluster"] == self.current_super_cluster
            ]["Compound"].tolist()
        else:
            # Fallback to old method if SuperCluster column doesn't exist
            super_cluster_idx = self.current_super_cluster - 1
            if hasattr(self.app.library, 'super_clusters') and super_cluster_idx < len(self.app.library.super_clusters):
                cluster_list = self.app.library.super_clusters[super_cluster_idx][2]
                return self.app.library.df[
                    self.app.library.df["Cluster"].isin(cluster_list)
                ]["Compound"].tolist()
            else:
                return []
