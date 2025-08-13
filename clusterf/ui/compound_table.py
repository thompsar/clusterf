from __future__ import annotations
from typing import TYPE_CHECKING
import panel as pn
import param
import pandas as pd

if TYPE_CHECKING:
    from clusterf.app import ClusterFApp


class CompoundTable(param.Parameterized):
    """
    A component for displaying compound data in a tabular format.
    Provides filtering, selection, and data management capabilities.
    """
    
    app: "ClusterFApp" = param.Parameter(default=None, doc="Reference to main ClusterF application")
    show_miss_compounds = param.Boolean(default=True, doc="Whether to show 'Miss' compounds")
    selected_compounds = param.List(default=[], doc="Currently selected compounds")
    current_super_cluster = param.Integer(default=None, doc="Current super cluster being displayed")
    
    def __init__(self, app: "ClusterFApp", **params):
        super().__init__(**params)
        self.app = app
        
        # Create the table widget
        self.table_widget = pn.widgets.Tabulator(
            pd.DataFrame(),
            sizing_mode="stretch_both",
            selectable="checkbox",
            disabled=True,
            show_index=False,
            pagination=None,
            margin=0,
            styles={"padding": "0px"}
        )
        
        # Create table controls
        self.controls = self._create_table_controls()
        
        # Create the main view
        self.view = pn.Card(
            pn.Column(
                self.controls,
                self.table_widget,
                sizing_mode="stretch_both",
                margin=0
            ),
            title="Compound Data Table",
            collapsed=False,
            margin=(5, 5)
        )
        
        # Set up event handlers
        self.table_widget.param.watch(self._on_table_selection, "selection")
    
    def _create_table_controls(self):
        """Create table control buttons."""
        # Toggle miss compounds button
        self.toggle_miss_button = pn.widgets.Button(
            name="Hide Misses",
            button_type="light",
            width=120,
            height=28
        )
        self.toggle_miss_button.on_click(self._toggle_miss_compounds)
        
        # Retest selected button
        self.retest_button = pn.widgets.Button(
            name="Retest Selected",
            button_type="light",
            width=120,
            height=28
        )
        self.retest_button.on_click(self._retest_selected)
        
        return pn.Row(
            self.toggle_miss_button,
            self.retest_button,
            sizing_mode="fixed",
            width=250,
            margin=(0, 2, 0, 0)
        )
    
    def update_table(self, compounds: list = None, super_cluster: int = None):
        """Update the table with new compound data."""
        if not self.app.library:
            return
        
        # Update current super cluster if provided
        if super_cluster is not None:
            self.current_super_cluster = super_cluster
        
        try:
            if compounds is not None:
                # Filter by specific compounds
                table_df = self.app.library.df[
                    self.app.library.df["Compound"].isin(compounds)
                ].copy()
            elif self.current_super_cluster is not None:
                # Filter by current super cluster
                if "SuperCluster" in self.app.library.df.columns:
                    table_df = self.app.library.df[
                        self.app.library.df["SuperCluster"] == self.current_super_cluster
                    ].copy()
                else:
                    # Fallback to old method
                    super_cluster_idx = self.current_super_cluster - 1
                    if hasattr(self.app.library, 'super_clusters') and super_cluster_idx < len(self.app.library.super_clusters):
                        cluster_list = self.app.library.super_clusters[super_cluster_idx][2]
                        table_df = self.app.library.df[
                            self.app.library.df["Cluster"].isin(cluster_list)
                        ].copy()
                    else:
                        table_df = pd.DataFrame()
            else:
                # Show all compounds
                table_df = self.app.library.df.copy()
            
            # Filter out "Miss" compounds if toggle is disabled
            if not self.show_miss_compounds:
                table_df = table_df[table_df["Category"] != "Miss"]
            
            # Select only the desired columns
            desired_columns = ["Compound", "clogP", "Cluster", "Category", "Retest"]
            available_columns = [col for col in desired_columns if col in table_df.columns]
            
            # If some columns are missing, add them with default values
            for col in desired_columns:
                if col not in table_df.columns:
                    if col == "Retest":
                        table_df[col] = False  # Default retest status
                    else:
                        table_df[col] = ""  # Default empty string
            
            # Select only the desired columns in the correct order
            table_df = table_df[desired_columns]
            
            # Reset index
            table_df = table_df.reset_index(drop=True)
            
            # Update the table
            self.table_widget.value = table_df
            self.table_widget.disabled = False
            
        except Exception as e:
            print(f"Error updating compound table: {e}")
            self.table_widget.value = pd.DataFrame()
    
    def _toggle_miss_compounds(self, event):
        """Toggle the display of 'Miss' compounds."""
        self.show_miss_compounds = not self.show_miss_compounds
        
        # Update button label
        if self.show_miss_compounds:
            self.toggle_miss_button.name = "Hide Misses"
        else:
            self.toggle_miss_button.name = "Show Misses"
        
        # Refresh table with current super cluster context
        self.update_table()
    
    def _retest_selected(self, event):
        """Toggle retest status for selected compounds."""
        if not self.table_widget.selection:
            return
        
        try:
            # Get selected row indices
            selected_indices = self.table_widget.selection
            selected_compounds = self.table_widget.value.loc[
                selected_indices, "Compound"
            ].values
            
            # Toggle retest status in library data
            compound_mask = self.app.library.df["Compound"].isin(selected_compounds)
            self.app.library.df.loc[compound_mask, "Retest"] = ~self.app.library.df.loc[
                compound_mask, "Retest"
            ]
            
            # Update dataset dataframe if available
            if hasattr(self.app.library, 'dataset_df'):
                dataset_mask = self.app.library.dataset_df["Compound"].isin(selected_compounds)
                self.app.library.dataset_df.loc[
                    dataset_mask, "Retest"
                ] = ~self.app.library.dataset_df.loc[dataset_mask, "Retest"]
            
            # Update subset dataframe if available
            if hasattr(self.app.library, 'subset_df'):
                subset_mask = self.app.library.subset_df["Compound"].isin(selected_compounds)
                self.app.library.subset_df.loc[subset_mask, "Retest"] = ~self.app.library.subset_df.loc[
                    subset_mask, "Retest"
                ]
            
            # Refresh table
            self.update_table()
            
            # Restore selection
            self.table_widget.selection = selected_indices
            
            print(f"Toggled Retest status for {len(selected_compounds)} compounds")
            
        except Exception as e:
            print(f"Error toggling retest status: {e}")
    
    def _on_table_selection(self, event):
        """Handle table selection changes."""
        if not event.new:
            self.selected_compounds = []
            return
        
        # Get selected compounds
        selected_indices = event.new
        self.selected_compounds = self.table_widget.value.loc[
            selected_indices, "Compound"
        ].values.tolist()
        
        # Notify app of selection change
        if hasattr(self.app, '_on_compound_selection_change'):
            self.app._on_compound_selection_change(self.selected_compounds)
    
    def get_selected_compounds(self):
        """Get currently selected compounds."""
        return self.selected_compounds.copy()
    
    def clear_super_cluster_context(self):
        """Clear the current super cluster context."""
        self.current_super_cluster = None
