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

    app: "ClusterFApp" = param.Parameter(
        default=None, doc="Reference to main ClusterF application"
    )
    show_miss_compounds = param.Boolean(
        default=False, doc="Whether to show 'Miss' compounds"
    )
    selected_compounds = param.List(default=[], doc="Currently selected compounds")
    current_super_cluster = param.Integer(
        default=None, doc="Current super cluster being displayed"
    )

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
            styles={"padding": "0px"},
        )

        self.table_widget.formatters = {"Retest": {"type": "tickCross"}}

        # Create table controls
        self.controls = self._create_table_controls()

        # Create the main view
        self.view = pn.Card(
            pn.Column(
                self.controls, self.table_widget, sizing_mode="stretch_both", margin=0
            ),
            title="Compound Data Table",
            collapsed=False,
            margin=(5, 5),
        )

        # Set up event handlers
        self.table_widget.param.watch(self._on_table_selection, "selection")

        # Watch for table value changes (when cells are edited)
        self.table_widget.param.watch(self._on_table_value_changed, "value")

    def _create_table_controls(self):
        """Create table control buttons."""
        # Toggle miss compounds button
        # Label will be set based on current state below
        self.toggle_miss_button = pn.widgets.Button(
            name="Show Misses", button_type="light", width=120, height=28
        )
        self.toggle_miss_button.on_click(self._toggle_miss_compounds)

        # Retest selected button
        self.retest_button = pn.widgets.Button(
            name="Retest Selected", button_type="light", width=120, height=28
        )
        self.retest_button.on_click(self._retest_selected)

        # Save Retest button
        self.save_retest_button = pn.widgets.Button(
            name="Save Retest", button_type="primary", width=120, height=28
        )
        self.save_retest_button.on_click(self._save_retest_changes)

        # Ensure button label reflects current state without recreating the widget
        self.toggle_miss_button.name = (
            "Hide Misses" if self.show_miss_compounds else "Show Misses"
        )
        return pn.Row(
            self.retest_button,
            self.save_retest_button,
            self.toggle_miss_button,
            sizing_mode="fixed",
            width=380,
            margin=(0, 2, 0, 0),
        )

    def _save_retest_changes(self, event=None):
        """Sync Retest changes and save only the primary dataset."""
        try:
            lib = getattr(self, "app", None)
            lib = getattr(lib, "library", None) if lib else None
            if lib is None:
                return
            # Ensure subset_df/dataset_df exist
            if not hasattr(lib, "subset_df") or not hasattr(lib, "dataset_df"):
                return

            # First, update working dataset for UI consistency
            if hasattr(lib, "sync_retest_to_dataset"):
                lib.sync_retest_to_dataset()

            # Then, sync and save only the primary dataset
            if hasattr(lib, "sync_retest_to_primary_dataset"):
                lib.sync_retest_to_primary_dataset()
            if hasattr(lib, "save_primary_dataset"):
                lib.save_primary_dataset()

        except Exception as e:
            print(f"Error saving retest changes: {e}")

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
                        self.app.library.df["SuperCluster"]
                        == self.current_super_cluster
                    ].copy()
                else:
                    # Fallback to old method
                    super_cluster_idx = self.current_super_cluster - 1
                    if hasattr(
                        self.app.library, "super_clusters"
                    ) and super_cluster_idx < len(self.app.library.super_clusters):
                        cluster_list = self.app.library.super_clusters[
                            super_cluster_idx
                        ][2]
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
            # TODO: add a column picker
            # available_columns = [col for col in desired_columns if col in table_df.columns]

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
            # create editors for the table from the column names, set all to none unless the column name is Retest
            self.table_widget.editors = {
                col: {"type": "tickCross"} if col == "Retest" else None
                for col in table_df.columns
            }
            # Update the table
            self.table_widget.value = table_df
            self.table_widget.disabled = False
            self._apply_category_styling()

        except Exception as e:
            print(f"Error updating compound table: {e}")
            self.table_widget.value = pd.DataFrame()

    def _apply_category_styling(self):
        """Apply styling to the Category column based on color dictionary."""
        if not hasattr(self.app, "color_picker") or not hasattr(
            self.app.color_picker, "color_dict"
        ):
            return

        color_dict = self.app.color_picker.color_dict
        if not color_dict:
            return

        def style_row(row):
            # Create a list of empty styles for all columns
            styles = [""] * len(row)

            # Only style the Category column
            if "Category" in row.index and pd.notna(row["Category"]):
                category = row["Category"]
                color = color_dict.get(category, "#FFFFFF")  # Default to white
                category_index = row.index.get_loc("Category")
                styles[category_index] = f"background-color: {color};"

            return styles

        self.table_widget.style.apply(style_row, axis=1)

    def update_colors(self, color_dict: dict):
        """Update colors and reapply styling to the table."""
        # The color_dict is already available through self.app.color_picker.color_dict
        # so we just need to reapply the styling
        self._apply_category_styling()
        self.update_table()

    def _toggle_miss_compounds(self, event):
        """Toggle the display of 'Miss' compounds."""
        self.show_miss_compounds = not self.show_miss_compounds
        self.clear_selection()
        # Update button label
        if self.show_miss_compounds:
            self.toggle_miss_button.name = "Hide Misses"
        else:
            self.toggle_miss_button.name = "Show Misses"

        # Refresh table with current super cluster context
        self.update_table()

        # Notify app of miss compounds toggle change
        if hasattr(self.app, "_on_miss_compounds_toggle"):
            self.app._on_miss_compounds_toggle(self.show_miss_compounds)

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
            if hasattr(self.app.library, "dataset_df"):
                dataset_mask = self.app.library.dataset_df["Compound"].isin(
                    selected_compounds
                )
                self.app.library.dataset_df.loc[
                    dataset_mask, "Retest"
                ] = ~self.app.library.dataset_df.loc[dataset_mask, "Retest"]

            # Update subset dataframe if available
            if hasattr(self.app.library, "subset_df"):
                subset_mask = self.app.library.subset_df["Compound"].isin(
                    selected_compounds
                )
                self.app.library.subset_df.loc[
                    subset_mask, "Retest"
                ] = ~self.app.library.subset_df.loc[subset_mask, "Retest"]

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

        # Get selected compounds from current table view
        selected_indices = event.new
        self.selected_compounds = self.table_widget.value.loc[
            selected_indices, "Compound"
        ].values.tolist()

        # Notify app of selection change
        if hasattr(self.app, "_on_compound_selection_change"):
            self.app._on_compound_selection_change(self.selected_compounds)

    def _get_full_context_compounds(self):
        """Get the full list of compounds from the current context (before any filtering)."""
        if not self.app.library:
            return []

        # Check if we have a current super cluster context
        if self.current_super_cluster is not None:
            return self.app.library.df[
                self.app.library.df["SuperCluster"] == self.current_super_cluster
            ]["Compound"].tolist()

        # If no specific context, return all compounds
        return self.app.library.df["Compound"].tolist()

    def get_selected_compounds(self):
        """Get currently selected compounds."""
        return self.selected_compounds.copy()

    def clear_super_cluster_context(self):
        """Clear the current super cluster context."""
        self.current_super_cluster = None

    def clear_selection(self):
        """Clear the table selection and reset selected compounds."""
        self.selected_compounds = []
        if hasattr(self, "table_widget"):
            self.table_widget.selection = []

    def _on_table_value_changed(self, event):
        """Handle table value changes when cells are edited in-place."""
        if event.new is None or event.old is None:
            return

        # Check if Retest column values have changed
        if "Retest" not in event.new.columns or "Retest" not in event.old.columns:
            return

        # Find compounds where Retest values have changed
        old_df = event.old
        new_df = event.new

        # Create a comparison DataFrame
        comparison = pd.merge(
            old_df[["Compound", "Retest"]].rename(columns={"Retest": "Retest_old"}),
            new_df[["Compound", "Retest"]].rename(columns={"Retest": "Retest_new"}),
            on="Compound",
            how="outer",
        )

        # Find rows where Retest values changed
        changed = comparison[comparison["Retest_old"] != comparison["Retest_new"]]

        if changed.empty:
            return

        # Propagate changes to backing DataFrames
        for _, row in changed.iterrows():
            compound = row["Compound"]
            new_val = bool(row["Retest_new"]) if pd.notna(row["Retest_new"]) else False

            # Update main library df
            if hasattr(self.app.library, "df") and isinstance(
                self.app.library.df, pd.DataFrame
            ):
                mask = self.app.library.df["Compound"] == compound
                if mask.any():
                    self.app.library.df.loc[mask, "Retest"] = new_val

            # Update dataset_df if present
            if hasattr(self.app.library, "dataset_df") and isinstance(
                self.app.library.dataset_df, pd.DataFrame
            ):
                dmask = self.app.library.dataset_df["Compound"] == compound
                if dmask.any():
                    self.app.library.dataset_df.loc[dmask, "Retest"] = new_val

            # Update subset_df if present
            if hasattr(self.app.library, "subset_df") and isinstance(
                self.app.library.subset_df, pd.DataFrame
            ):
                smask = self.app.library.subset_df["Compound"] == compound
                if smask.any():
                    self.app.library.subset_df.loc[smask, "Retest"] = new_val

        # Keep dtypes consistent
        try:
            self.app.library.df["Retest"] = self.app.library.df["Retest"].astype(
                "boolean"
            )
            if hasattr(self.app.library, "dataset_df") and isinstance(
                self.app.library.dataset_df, pd.DataFrame
            ):
                self.app.library.dataset_df["Retest"] = self.app.library.dataset_df[
                    "Retest"
                ].astype("boolean")
            if hasattr(self.app.library, "subset_df") and isinstance(
                self.app.library.subset_df, pd.DataFrame
            ):
                self.app.library.subset_df["Retest"] = self.app.library.subset_df[
                    "Retest"
                ].astype("boolean")
        except Exception:
            pass

    def reset(self):
        """Reset the table to the initial state."""
        self.show_miss_compounds = False
        self.selected_compounds = []
        self.current_super_cluster = None
        if hasattr(self, "table_widget"):
            self.table_widget.value = pd.DataFrame()
            self.table_widget.disabled = True
            self.table_widget.selection = []
        # Reset toggle button label to reflect hidden misses
        if hasattr(self, "toggle_miss_button"):
            self.toggle_miss_button.name = "Show Misses"
        # Ensure downstream components pick up the default state
        try:
            self.update_table()
        except Exception:
            pass
