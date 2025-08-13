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
        
        if not self.app.library:
            return
        
        try:
            # Filter out non-numeric compound IDs (category names, etc.)
            filtered_compound_ids = []
            for item in compound_ids:
                try:
                    # Try to convert to int to ensure it's a valid compound ID
                    compound_id = int(item)
                    filtered_compound_ids.append(compound_id)
                except (ValueError, TypeError):
                    continue
            
            if not filtered_compound_ids:
                self.carousel.svgs = []
                self.compounds = []
                return
            
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
                # Ensure color_dict has all necessary categories
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
