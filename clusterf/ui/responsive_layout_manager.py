from __future__ import annotations
from typing import TYPE_CHECKING
import panel as pn
import param

if TYPE_CHECKING:
    from clusterf.app import ClusterFApp


class ResponsiveLayoutManager(param.Parameterized):
    """
    Manages responsive layout that adapts to screen size and user preferences.
    Provides both tab-based (compact) and grid-based (expanded) layouts.
    """
    
    # Layout modes
    layout_mode = param.Selector(
        objects=["Auto", "Compact", "Expanded"], 
        default="Auto",
        doc="Layout mode: Auto adapts to screen size, Compact uses tabs, Expanded uses grid"
    )
    
    # Screen size breakpoints (in pixels)
    COMPACT_BREAKPOINT = param.Integer(default=1200, doc="Screen width below which to use compact layout")
    EXPANDED_BREAKPOINT = param.Integer(default=1600, doc="Screen width above which to use expanded layout")
    
    # Layout state
    current_screen_width = param.Integer(default=1200, doc="Current screen width for auto mode")
    
    def __init__(self, app: "ClusterFApp", **params):
        super().__init__(**params)
        self.app = app
        
        # Layout containers
        self.tab_layout = None
        self.grid_layout = None
        self.current_layout = None
        
        # Initialize layouts
        self._create_tab_layout()
        self._create_grid_layout()
        
        # Set initial layout
        self._update_layout()
    
    def _create_tab_layout(self):
        """Create tab-based layout for compact screens."""
        # Create placeholder content for now - will be populated when components are available
        self.tab_layout = pn.Tabs(
            ("Network Graph", pn.pane.Markdown("Network graph will appear here")),
            ("Compound Grid", pn.pane.Markdown("Compound grid will appear here")),
            ("Category Analysis", pn.pane.Markdown("Category histogram will appear here")),
            ("Compound Data", pn.pane.Markdown("Compound data chart will appear here")),
            ("Compound Table", pn.pane.Markdown("Compound table will appear here")),
            sizing_mode="stretch_both"
        )
    
    def _create_grid_layout(self):
        """Create grid layout for expanded screens."""
        self.grid_layout = pn.GridSpec(
            nrows=3, ncols=3,
            sizing_mode="stretch_both",
            margin=0,
            name="Expanded Layout",
            styles={"gap": "2px", "padding": "2px", "box-sizing": "border-box"}
        )
        
        # Placeholder content for now - using non-overlapping regions
        self.grid_layout[0:2, 0:2] = pn.pane.Markdown("Network Graph (2x2)")
        self.grid_layout[0:2, 2] = pn.pane.Markdown("Compound Grid (2x1)")
        self.grid_layout[2, 0] = pn.pane.Markdown("Category Histogram")
        self.grid_layout[2, 1] = pn.pane.Markdown("Compound Data Chart")
        self.grid_layout[2, 2] = pn.pane.Markdown("Compound Table")
    
    def update_components(self, components: dict):
        """
        Update the layout with actual components.
        
        Args:
            components: Dictionary with keys: 'cluster_graph', 'compound_grid', 
                       'category_histogram', 'compound_table', 'compound_data_chart'
        """
        # Update tab layout
        if 'cluster_graph' in components:
            self.tab_layout[0] = ("Network Graph", components['cluster_graph'])
        if 'compound_grid' in components:
            self.tab_layout[1] = ("Compound Grid", components['compound_grid'])
        if 'category_histogram' in components:
            self.tab_layout[2] = ("Category Analysis", components['category_histogram'])
        if 'compound_data_chart' in components:
            self.tab_layout[3] = ("Compound Data", components['compound_data_chart'])
        if 'compound_table' in components:
            self.tab_layout[4] = ("Compound Table", components['compound_table'])
        
        # Recreate grid layout to avoid overlapping issues
        self._create_grid_layout()
        
        # Update grid layout with new components
        if 'cluster_graph' in components:
            self.grid_layout[0:2, 0:2] = components['cluster_graph']
        if 'compound_grid' in components:
            self.grid_layout[0:2, 2] = components['compound_grid']
        if 'category_histogram' in components:
            self.grid_layout[2, 0] = components['category_histogram']
        if 'compound_data_chart' in components:
            self.grid_layout[2, 1] = components['compound_data_chart']
        if 'compound_table' in components:
            self.grid_layout[2, 2] = components['compound_table']
        
        # Update current layout
        self._update_layout()
    
    def _update_layout(self):
        """Update the current layout based on mode and screen size."""
        if self.layout_mode == "Auto":
            if self.current_screen_width < self.COMPACT_BREAKPOINT:
                self.current_layout = self.tab_layout
            else:
                self.current_layout = self.grid_layout
        elif self.layout_mode == "Compact":
            self.current_layout = self.tab_layout
        else:  # Expanded
            self.current_layout = self.grid_layout
    
    def set_screen_width(self, width: int):
        """Update screen width for auto mode calculations."""
        self.current_screen_width = width
        if self.layout_mode == "Auto":
            self._update_layout()
    
    def get_current_layout(self):
        """Get the current active layout."""
        return self.current_layout
    
    def get_layout_controls(self):
        """Get layout mode selector controls."""
        return pn.Row(
            pn.widgets.RadioButtonGroup.from_param(
                self.param.layout_mode,
                name="Layout Mode",
                button_type="light",
                width=200
            ),
            pn.Spacer(),
            margin=(5, 5)
        )
