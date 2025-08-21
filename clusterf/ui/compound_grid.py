from __future__ import annotations
from typing import TYPE_CHECKING
import panel as pn
import param

if TYPE_CHECKING:
    from clusterf.app import ClusterFApp


class CompoundGrid(param.Parameterized):
    """
    A component for displaying compound structures in a grid layout.
    Handles compound visualization with carousel functionality.
    """
    
    app: "ClusterFApp" = param.Parameter(default=None, doc="Reference to main ClusterF application")
    compounds = param.List(default=[], doc="List of compound IDs to display")
    color_dict = param.Dict(default={}, doc="Color mapping for compound categories")
    show_miss_compounds = param.Boolean(default=False, doc="Whether to show 'Miss' compounds")
    
    def __init__(self, app: "ClusterFApp", **params):
        super().__init__(**params)
        self.app = app
        
        # Create the carousel component (reusing existing implementation)
        from clusterf.ui.carrousel import Carrousel
        self.carousel = Carrousel()
        
        # Create the main view
        self.view = pn.Card(
            self.carousel.view(),
            title="Compound Structures",
            collapsed=False,
            margin=(5, 5)
        )
    
    def update_compounds(self, compound_ids: list):
        """Update the compound grid with new compounds."""
        if not compound_ids:
            self.carousel.svgs = []
            self.common_substructure = None
            return
        
        try:
            
            # Filter out "Miss" compounds if show_miss_compounds is False
            if not self.show_miss_compounds:
                # Get compound data to check categories
                compound_data = self.app.library.get_compounds(compound_ids)
                # Keep only non-Miss compounds
                non_miss_mask = compound_data["Category"] != "Miss"
                filtered_compound_ids = compound_data.loc[non_miss_mask, "Compound"].tolist()
                
                if not filtered_compound_ids:
                    self.carousel.svgs = []
                    self.compounds = []
                    return
            else:
                filtered_compound_ids = compound_ids
            
            # Generate compound grid images
            # Determine grid parameters based on number of compounds
            if len(filtered_compound_ids) > 1:
                mols_per_row = 4
                max_rows = 3
                orient = True
            else:
                mols_per_row = 4
                max_rows = 1
                orient = False

            # Generate the grid with current colors
            try:
                # Get compound data for categories and reuse if already fetched
                if not self.show_miss_compounds and 'compound_data' in locals():
                    # Reuse the compound data we already fetched for filtering
                    pass
                else:
                    # Get the actual categories from the compounds
                    compound_data = self.app.library.get_compounds(filtered_compound_ids)
                
                actual_categories = set(compound_data["Category"].unique())
                
                # Get the current color dictionary from the app's color picker
                app_color_dict = {}
                if hasattr(self.app, 'color_picker') and hasattr(self.app.color_picker, 'color_dict'):
                    app_color_dict = self.app.color_picker.color_dict.copy()
                
                # Create a complete color dictionary
                complete_color_dict = {}
                for category in actual_categories:
                    if category in app_color_dict:
                        complete_color_dict[category] = app_color_dict[category]
                    elif category in self.color_dict:
                        complete_color_dict[category] = self.color_dict[category]
                    else:
                        complete_color_dict[category] = "#999999"  # Default color
                
                # Ensure "Miss" category is handled
                if "Miss" in actual_categories:
                    complete_color_dict["Miss"] = "none"
                
                svgs, self.common_substructure = self.app.library.draw_compound_grid(
                    filtered_compound_ids,
                    mols_per_row=mols_per_row,
                    max_rows=max_rows,
                    color_dict=complete_color_dict,
                    orient=orient,
                    legend=False,
                )
                self.carousel.svgs = svgs
                self.compounds = filtered_compound_ids
            except Exception as grid_error:
                print(f"Error updating compound grid: {grid_error}")
                self.carousel.svgs = []
                self.compounds = []
            
        except Exception as e:
            print(f"Error updating compound grid: {e}")
            self.carousel.svgs = []
    
    def update_colors(self, color_dict: dict):
        """Update the color scheme for compound visualization."""
        self.color_dict = color_dict.copy()
        # Re-render compounds with new colors if needed
        if self.compounds:
            self.update_compounds(self.compounds)
    
    def get_show_miss_compounds(self) -> bool:
        """Get the current state of show_miss_compounds."""
        return self.show_miss_compounds
    
    def set_show_miss_compounds(self, show_miss_compounds: bool):
        """Set the show_miss_compounds state and update the grid if needed."""
        if self.show_miss_compounds != show_miss_compounds:
            self.show_miss_compounds = show_miss_compounds
            # Update the grid with current compounds to apply the new filter
            if self.compounds:
                self.update_compounds(self.compounds)

    def reset(self):
        """Reset the grid to the initial state."""
        self.compounds = []
        self.color_dict = {}
        self.show_miss_compounds = False
        if hasattr(self, "carousel"):
            self.carousel.svgs = []
