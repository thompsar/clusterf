import pandas as pd
import panel as pn
import param
from typing import Dict, List, Optional, Callable


class CompoundTableManager(param.Parameterized):
    """
    A class to manage compound table functionality including display, styling,
    filtering, and interactions. This class encapsulates all table-related
    operations that were previously part of the main ClusterF class.
    """
    
    # Parameters for controlling table state
    show_miss_compounds = param.Boolean(default=True)
    table_visible = param.Boolean(default=False)
    visible_columns = param.List(default=["Compound"])
    
    # Action parameters
    toggle_miss_button = param.Action(
        lambda x: x.param.trigger("toggle_miss_button"), label="Hide Misses"
    )
    retest_selected_button = param.Action(
        lambda x: x.param.trigger("retest_selected_button"), label="Retest Selected"
    )
    save_button = param.Action(
        lambda x: x.param.trigger("save_button"), label="Save Spreadsheet"
    )
    
    def __init__(self, library=None, color_dict=None, **params):
        super().__init__(**params)
        
        self.library = library
        self.color_dict = color_dict or {}
        self.subcategory_columns = {}
        
        # Initialize the compound table
        self.compound_table = pn.widgets.Tabulator(
            pd.DataFrame(),  # Start with empty dataframe
            sizing_mode="stretch_both",
            selectable="checkbox",  # Enable multiple row selection
            disabled=True,
            visible=False,
            show_index=False,
            pagination=None,
            margin=0,
            styles={"padding": "0px"},
        )
        self.compound_table.editable = True
        
        # Enable in-cell toggling for Retest and show symbols instead of True/False
        self.compound_table.editors = {"Retest": {"type": "tickCross"}}
        self.compound_table.formatters = {"Retest": {"type": "tickCross"}}
        
        # Set up table event handlers
        self.compound_table.on_click(self.on_table_click)
        self.compound_table.param.watch(self.on_table_selection, "selection")
        
        # Cache the previous table dataframe to detect cell edits
        self._table_cache = None
        # Watch for any edits in the table value
        self.compound_table.param.watch(self.on_table_edited, "value")
        
        # Callbacks for external interactions
        self.on_compound_selected_callback = None
        self.on_compounds_selected_callback = None
        self.on_retest_changed_callback = None
        
        # Watch parameter changes
        self.param.watch(self.toggle_miss_compounds, "toggle_miss_button")
        self.param.watch(self.retest_selected_compounds, "retest_selected_button")
        self.param.watch(self.save_spreadsheet, "save_button")
    
    def set_library(self, library):
        """Set the library instance and update the table."""
        self.library = library
        if hasattr(library, 'df'):
            # Reset table with library data
            self.update_visible_columns()
    
    def set_color_dict(self, color_dict: Dict[str, str]):
        """Set the color dictionary for styling."""
        self.color_dict = color_dict
        self.apply_styling()
    
    def set_subcategory_columns(self, subcategory_columns: Dict[str, Dict]):
        """Set subcategory columns data."""
        self.subcategory_columns = subcategory_columns
        self.update_visible_columns()
    
    def update_visible_columns(self, fine_threshold: Optional[float] = None):
        """Update the visible columns list based on current state."""
        if not self.library or not hasattr(self.library, 'df'):
            return
            
        base_columns = ["Compound"]
        if fine_threshold is not None:
            base_columns.append(str(fine_threshold))
        
        # Add subcategory columns if they exist
        if self.subcategory_columns:
            additional_columns = list(self.subcategory_columns.keys()) + ["Category", "Retest"]
            self.visible_columns = base_columns + additional_columns
        else:
            self.visible_columns = base_columns
    
    def build_table(self, clusters_or_super_cluster, fine_threshold: float) -> pd.DataFrame:
        """
        Build table from either a super cluster number or list of cluster nodes.

        Args:
            clusters_or_super_cluster: Either:
                - int/float: Super cluster number (1-indexed) to show all compounds in that super cluster
                - list: List of cluster node IDs to show compounds from specific selected clusters
            fine_threshold: The fine threshold value for clustering

        Returns:
            pd.DataFrame: Filtered dataframe with compounds from specified clusters
        """
        if not self.library or not hasattr(self.library, 'df'):
            return pd.DataFrame(columns=self.visible_columns)
            
        table_df = self.library.df.copy()

        # Add subcategory columns if they exist
        if self.subcategory_columns:
            for col_name, col_data in self.subcategory_columns.items():
                # Map compound values to the column
                table_df[col_name] = table_df["Compound"].map(col_data)

        # Check if input is a single integer (super cluster number)
        if isinstance(clusters_or_super_cluster, (int, float)):
            # Filter by super cluster number using the SuperCluster column (if available)
            if "SuperCluster" in table_df.columns:
                table_df = table_df[
                    table_df["SuperCluster"] == clusters_or_super_cluster
                ]
            else:
                # Fallback to old method if SuperCluster column doesn't exist
                super_cluster_idx = int(clusters_or_super_cluster) - 1
                if hasattr(self.library, "super_clusters") and super_cluster_idx < len(
                    self.library.super_clusters
                ):
                    cluster_list = self.library.super_clusters[super_cluster_idx][2]
                    table_df = table_df[
                        table_df[str(fine_threshold)].isin(cluster_list)
                    ]
                else:
                    return pd.DataFrame(
                        columns=self.visible_columns
                    )  # Return empty table
        else:
            # Input is a list of cluster nodes - filter by those specific clusters
            if not isinstance(clusters_or_super_cluster, list):
                clusters_or_super_cluster = [clusters_or_super_cluster]
            table_df = table_df[
                table_df[str(fine_threshold)].isin(clusters_or_super_cluster)
            ]

        # Filter out "Miss" compounds if toggle is disabled
        if not self.show_miss_compounds:
            table_df = table_df[table_df["Category"] != "Miss"]

        # Ensure all visible columns exist in table_df
        available_columns = [
            col for col in self.visible_columns if col in table_df.columns
        ]
        table_df = table_df[available_columns].reset_index()

        return table_df
    
    def update_table(self, clusters_or_super_cluster, fine_threshold: float, 
                    preserve_selection: bool = False):
        """Update the table with new data."""
        # Store current parameters for potential refresh
        self._current_clusters = clusters_or_super_cluster
        self._current_fine_threshold = fine_threshold
        
        # Preserve current selection if requested
        selected_compounds = []
        if preserve_selection and self.compound_table.selection:
            selected_compounds = self.compound_table.value.loc[
                self.compound_table.selection, "Compound"
            ].values
        
        # Build new table
        new_table = self.build_table(clusters_or_super_cluster, fine_threshold)
        self.compound_table.value = new_table
        
        # Restore selection if preserved
        if preserve_selection and len(selected_compounds) > 0:
            self.compound_table.selection = new_table[
                new_table["Compound"].isin(selected_compounds)
            ].index.tolist()
        
        # Apply styling
        self.apply_styling()
        
        return new_table
    
    def apply_styling(self):
        """Apply styles to the compound table based on cluster colors."""
        if not hasattr(self, "compound_table") or not self.color_dict:
            return

        # Define a function to style each row based on category
        def style_row(row):
            if hasattr(row, 'Category'):
                color = self.color_dict.get(row.Category, "#FFFFFF")  # Default to white
                return [f"background-color: {color};"] * len(row)
            return [""] * len(row)

        try:
            self.compound_table.style.apply(style_row, axis=1)
        except Exception:
            # Silently handle styling errors
            pass
    
    def toggle_miss_compounds(self, event=None):
        """Update table display when Miss compounds toggle button is pressed."""
        # Toggle the boolean value
        self.show_miss_compounds = not self.show_miss_compounds

        # Update button label based on current state
        if self.show_miss_compounds:
            self.param.toggle_miss_button.label = "Hide Misses"
        else:
            self.param.toggle_miss_button.label = "Show Misses"

        # Trigger callback if set
        if hasattr(self, '_on_miss_toggle_callback') and self._on_miss_toggle_callback:
            self._on_miss_toggle_callback()
    
    def retest_selected_compounds(self, event=None):
        """Toggle Retest status for all selected compounds in the table."""
        if not self.compound_table.selection or not self.library:
            return

        # Get selected row indices
        selected_indices = self.compound_table.selection
        selected_compounds = self.compound_table.value.loc[
            selected_indices, "Compound"
        ].values

        # Update library dataframes
        self._update_retest_status(selected_compounds)
        
        # Refresh the table to show updated values
        self._refresh_table_after_retest_change(selected_indices)
        
        # Trigger callback if set
        if self.on_retest_changed_callback:
            self.on_retest_changed_callback(selected_compounds)

        print(f"Toggled Retest status for {len(selected_compounds)} compounds")
    
    def _update_retest_status(self, compound_ids: List[str]):
        """Update Retest status in all relevant dataframes."""
        if not self.library:
            return
            
        # Update main library df
        compound_mask = self.library.df["Compound"].isin(compound_ids)
        self.library.df.loc[compound_mask, "Retest"] = ~self.library.df.loc[
            compound_mask, "Retest"
        ]

        # Update the dataset dataframe if it exists
        if hasattr(self.library, 'dataset_df'):
            dataset_mask = self.library.dataset_df["Compound"].isin(compound_ids)
            self.library.dataset_df.loc[
                dataset_mask, "Retest"
            ] = ~self.library.dataset_df.loc[dataset_mask, "Retest"]

        # Update the subset dataframe if it exists
        if hasattr(self.library, 'subset_df'):
            subset_mask = self.library.subset_df["Compound"].isin(compound_ids)
            self.library.subset_df.loc[subset_mask, "Retest"] = ~self.library.subset_df.loc[
                subset_mask, "Retest"
            ]
    
    def save_spreadsheet(self, event=None):
        """Save the current table as an Excel file."""
        if self.compound_table.value is not None and not self.compound_table.value.empty:
            self.compound_table.value.to_excel("compound_list.xlsx", index=False)
            print("Compound table saved as compound_list.xlsx")
    
    def on_table_click(self, event):
        """Handle table cell clicks."""
        # Determine column name from event
        try:
            col = event.column
            if isinstance(col, int):
                if 0 <= col < len(self.compound_table.value.columns):
                    col_name = self.compound_table.value.columns[col]
                else:
                    col_name = None
            else:
                col_name = col
        except Exception:
            col_name = None

        # If clicking the Retest cell, toggle its value in-place and propagate
        if col_name == "Retest":
            self._handle_retest_cell_click(event)
            return  # Don't treat this click as a row selection

        # Otherwise: normal behavior — select the row and trigger callback
        if hasattr(event, 'row') and event.row is not None:
            compound = self.compound_table.value.loc[event.row, "Compound"]
            self.compound_table.selection = [event.row]
            
            # Trigger callback for single compound selection
            if self.on_compound_selected_callback:
                self.on_compound_selected_callback(compound, event.row)
    
    def _handle_retest_cell_click(self, event):
        """Handle clicks specifically on Retest cells."""
        try:
            row = event.row
            df = self.compound_table.value
            if (
                row is None
                or "Retest" not in df.columns
                or "Compound" not in df.columns
            ):
                return

            compound = df.at[row, "Compound"]
            curr_val = (
                bool(df.at[row, "Retest"])
                if pd.notna(df.at[row, "Retest"])
                else False
            )
            new_val = not curr_val

            # Update the visible table cell immediately
            df.at[row, "Retest"] = new_val
            
            # Propagate change to backing DataFrames first
            self._update_retest_status([compound])
            
            # Now rebuild table from source to ensure consistency
            if hasattr(self, '_current_clusters') and hasattr(self, '_current_fine_threshold'):
                new_table = self.build_table(self._current_clusters, self._current_fine_threshold)
                self.compound_table.value = new_table
                # Find and restore selection
                try:
                    new_row = new_table[new_table["Compound"] == compound].index[0]
                    self.compound_table.selection = [new_row]
                except (IndexError, KeyError):
                    pass
            else:
                # Fallback to simple update
                self.compound_table.value = df
            
            # Re-apply table styling
            self.apply_styling()
            
            # Trigger callback if set
            if self.on_retest_changed_callback:
                self.on_retest_changed_callback([compound])

        except Exception as e:
            print(f"Error handling Retest cell click: {e}")
    
    def on_table_selection(self, event):
        """Handle multiple row selection in the compound table."""
        if not event.new:
            # No selection - trigger callback with empty list
            if self.on_compounds_selected_callback:
                self.on_compounds_selected_callback([])
            return
            
        # Get selected row indices
        selected_indices = event.new  # always a list
        selected_compounds = self.compound_table.value.loc[
            selected_indices, "Compound"
        ].values

        # Trigger callback for multiple compound selection
        if self.on_compounds_selected_callback:
            self.on_compounds_selected_callback(selected_compounds)
    
    def on_table_edited(self, event):
        """Handle table edits, particularly for Retest column changes."""
        try:
            import pandas as _pd  # local alias
        except Exception:
            _pd = pd

        new_df = event.new
        old_df = (
            event.old if event.old is not None else getattr(self, "_table_cache", None)
        )
        # Always update cache for subsequent diffs
        if hasattr(new_df, "copy"):
            self._table_cache = new_df.copy()
        else:
            self._table_cache = new_df

        # Validate frames/columns
        if new_df is None or not isinstance(new_df, _pd.DataFrame):
            return
        if "Compound" not in new_df.columns or "Retest" not in new_df.columns:
            return
        if old_df is None or not isinstance(old_df, _pd.DataFrame):
            # First render; nothing to compare yet
            return

        # Build a diff of Retest values keyed by Compound
        try:
            left = old_df[["Compound", "Retest"]].rename(
                columns={"Retest": "Retest_old"}
            )
            right = new_df[["Compound", "Retest"]].rename(
                columns={"Retest": "Retest_new"}
            )
            merged = left.merge(right, on="Compound", how="outer")
        except Exception:
            return

        changed = merged[merged["Retest_old"] != merged["Retest_new"]]
        if changed.empty:
            return

        # Apply changes to all backing DataFrames
        changed_compounds = []
        for _, row in changed.iterrows():
            comp = row["Compound"]
            changed_compounds.append(comp)

        if changed_compounds:
            self._update_retest_status(changed_compounds)
            
            # Re-apply table styling
            self.apply_styling()
            
            # Trigger callback if set
            if self.on_retest_changed_callback:
                self.on_retest_changed_callback(changed_compounds)
    
    def set_callbacks(self, 
                     on_compound_selected: Optional[Callable] = None,
                     on_compounds_selected: Optional[Callable] = None,
                     on_retest_changed: Optional[Callable] = None,
                     on_miss_toggle: Optional[Callable] = None):
        """Set callback functions for table interactions."""
        self.on_compound_selected_callback = on_compound_selected
        self.on_compounds_selected_callback = on_compounds_selected
        self.on_retest_changed_callback = on_retest_changed
        self._on_miss_toggle_callback = on_miss_toggle
    
    def show(self):
        """Show the table."""
        self.compound_table.visible = True
        self.table_visible = True
    
    def hide(self):
        """Hide the table."""
        self.compound_table.visible = False
        self.table_visible = False
    
    def clear_selection(self):
        """Clear the current table selection."""
        self.compound_table.selection = []
    
    def get_selected_compounds(self) -> List[str]:
        """Get the currently selected compound IDs."""
        if not self.compound_table.selection:
            return []
        return self.compound_table.value.loc[
            self.compound_table.selection, "Compound"
        ].values.tolist()
    
    def get_table_controls_panel(self):
        """Create a panel with table control buttons."""
        if not self.table_visible:
            return pn.Spacer(width=0, height=0)

        return pn.Column(
            pn.Param(
                self.param,
                parameters=["toggle_miss_button"],
                widgets={
                    "toggle_miss_button": {
                        "type": pn.widgets.Button,
                        "width": 120,
                        "height": 28,
                    }
                },
                margin=(0, 0, 2, 0),
                show_name=False,
                sizing_mode="fixed",
                width=122,
            ),
            pn.Param(
                self.param,
                parameters=["retest_selected_button"],
                widgets={
                    "retest_selected_button": {
                        "type": pn.widgets.Button,
                        "width": 120,
                        "height": 28,
                    }
                },
                margin=(0, 0, 0, 0),
                show_name=False,
                sizing_mode="fixed",
                width=122,
            ),
            sizing_mode="fixed",
            width=122,
            margin=(0, 2, 0, 0),
        )
    
    def _refresh_table_after_retest_change(self, selected_indices):
        """Refresh the table display after retest changes to show updated symbols."""
        try:
            # Get the current data source and rebuild the table to reflect backing dataframe changes
            if hasattr(self, '_current_clusters') and hasattr(self, '_current_fine_threshold'):
                # Rebuild the table from source data
                new_table = self.build_table(self._current_clusters, self._current_fine_threshold)
                self.compound_table.value = new_table
            else:
                # Fallback: just force a refresh of the current data
                current_data = self.compound_table.value
                self.compound_table.value = current_data.copy()
            
            # Restore the selection
            self.compound_table.selection = selected_indices
            
            # Re-apply styling
            self.apply_styling()
        except Exception as e:
            print(f"Error refreshing table after retest change: {e}")
    
    @property
    def widget(self):
        """Get the underlying Tabulator widget."""
        return self.compound_table
