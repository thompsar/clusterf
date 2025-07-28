import re
import os
import panel as pn
import numpy as np
import pandas as pd
import holoviews as hv
from holoviews.streams import Selection1D
import param
import networkx as nx
from src.chemistry import ChemLibrary
from src.carrousel import Carrousel

pn.extension("tabulator", sizing_mode="stretch_width")

"""
TODO LIST:
- [x] Orient selected single compounds with those in the cluster
- [ ] Add more detail to the compound table (LogP, MW, etc.), definitely add Category!
    - [ ] Add histogram of logP and MW to the graph plot?
- [ ] Add plots accoriding to values in the table (e.g. MW vs LogP)
    - [ ] Plot values should be highlighted to to reveal compounds of interest. Interesting patterns may appear.
- [x] Convert compound ID text box into "compound search" box, have that adjust to appropraite cluster/graph
    - [x] Currently only seeks super cluster. Add selection of exact cluster within super cluster.
- [ ] Add a "reset" button to reset the graph to the original view
- [x] Results of reclustering should be saved to the CSV of the subset file
- [ ] Cluster hover should provide basic stats about cluster (number of compounds per category, total compounds, etc.)
- [x] Subset select should have a drop down menu that shows the subset categories, with color pickers for colorizing the clusters
- [x] Selection of a cluster(s) draws a grid of compounds, all aligned
    - [x] Breaks for larger clusters. Fix this! is there a way to make a carousel of images?
        - See https://discourse.holoviz.org/t/is-there-any-widget-equivalent-to-a-carousel/3431/3
        -[x] BUG: Cluster 155 is getting a color swatch in an empty grid box at the very end. Fix this!!!
        -[x] BUG: Switching clusters needs to reset the selected index of the carrosel. Fix this!
        - [x] TODO: draw_compounds in chemistry.py doesnt return a list for a single compound, but the carrosel expects a list.
    - [x] Issue: chemistry.py:310: RuntimeWarning: More than 20 figures have been opened. Figures created through the pyplot interface 
          (`matplotlib.pyplot.figure`) are retained until explicitly closed and may consume too much memory. (To control this warning, 
          see the rcParam `figure.max_open_warning`). Consider using `matplotlib.pyplot.close()`.
        - Dealt with by switching over to SVG, which is a much better solution for several reasons.
- [x] Colorize table based on category
- [ ] Orienting the compounds in cluster 1772 creates some atom collisions and weird orientations. Fix this.
- [x] Restore the ability to sort/filter the table
        - solving the table selection bug that kept popping up seemed to fix this. fix was setting pagination=None for some reason?
- [ ] Clique highlighing? Common structure?
    - [ ] Would be really cool to have a philogenetic tree of the variations in chemical strucutre from the core structure
- [ ] Size nodes according, roughly to the number of compounds in the cluster? Or perhaps different shape or outline for singletons?
- [ ] Pull rendering of cluster stats chart out of chemistry.py and into here (or make own class)
    - [ ] Label axes accordingly
    - [ ] Convert to a bar chart?
    - [ ] Make clickable to select super cluster
    - [ ] Colorize based on category
- [x] BUG: deselecting a node does not reset the compound grid or table
- [ ] BUG: ? Where is cluster 5030? It is not a singleton but it is not in the cluster charts.
    - 5030 gets tossed out during building of graph due to the fact that none of the compounds included match with other clusters under
        the conditions ive mostly tested (0.2, 0.4). Modifying fingerprinting scheme may help, but hasn't worked yet.
    - [ ] FIX: Re-introduce pruned clusters to graph so that Compound Search works (fails currently). Do this by checking to see if they contain any compounds of interest!
-[ ] Add a way to triage compounds in the table (radio boxes that persist through selections, for example).

NEW:
-[ ] Switching subsets causes key error with Category, likely need to reset the color widgets and histogram.
-[ ] Rendering compound data causes table to scroll back to the top. Seems like a persistent bug over the years in tabulator, may have to live with it for now.
-[x] Filter out Misses from the table
-[ ] Allow for selection of multiple compounds from the table. Draw those particular compounds in the grid. Plot their data in the compound data chart.
    - [x] Implement multi-select functionality in the compound table
    - [x] Update the compound grid to display selected compounds
    - [x] Update the compound data chart to display data for selected compounds
    - [x] Remove rendering of single compounds in the control panel, instead render the selected compounds in the compound grid?
    - [ ] Plot scatter and error bars for selected compounds in the compound data chart

"""


# def category_sort_key(category):
#     category_str = str(category).lower()
#     if (
#         "hit" in category_str
#         and "(nr)" not in category_str
#         and "questionable" not in category_str
#     ):
#         return (0, category)  # Highest priority
#     elif "(nr)" in category_str or "questionable" in category_str:
#         return (1, category)  # Second priority
#     elif "interfering" in category_str:
#         return (2, category)  # Third priority
#     else:
#         return (3, category)  # Lowest priority (includes "Miss" and others)

def category_sort_key(category):
    category_str = str(category).lower()
    if (
        "universal" in category_str
    ):
        return (0, category)  # Highest priority
    elif "selective" in category_str:
        return (1, category)  # Second priority
    elif "bidirectional" in category_str:
        return (2, category)  # Third priority
    elif "interfering" in category_str:
        return (3, category)  # Third priority
    else:
        return (4, category)  # Lowest priority (includes "Miss" and others)


class ClusterF(param.Parameterized):
    libraries = ["cns", "cns_kinase"]  # 'diverset': ChemLibrary('diverset')}

    # Widgets
    lib_select = param.Selector(objects=libraries, default=libraries[1])
    fine_threshold = param.Selector(objects=[0.2], default=0.2)
    coarse_threshold = param.Number(0.4)
    dataset_select = param.FileSelector(path=os.path.join("compound_subsets", "*.csv*"))
    compound_input = param.String()
    selected_compound = param.String()
    selected_nodes = param.List()  # Track the selected node
    recluster_button = param.Action(
        lambda x: x.param.trigger("recluster_button"), label="Recluster Library"
    )
    cluster_slider = param.Selector()
    search_button = param.Action(
        lambda x: x.param.trigger("search_button"), label="Search"
    )
    save_button = param.Action(
        lambda x: x.param.trigger("save_button"), label="Save Spreadsheet"
    )
    save_dataset_button = param.Action(
        lambda x: x.param.trigger("save_dataset_button"), label="Save Dataset"
    )
    show_miss_compounds = param.Boolean(default=True)
    toggle_miss_button = param.Action(
        lambda x: x.param.trigger("toggle_miss_button"), label="Hide Miss Compounds"
    )
    retest_selected_button = param.Action(
        lambda x: x.param.trigger("retest_selected_button"), label="Retest Selected"
    )

    # Dynamic color parameters will be added programmatically
    categories = param.List(default=[])
    color_widgets_visible = param.Boolean(default=False)

    def __init__(self, **params):
        super().__init__(**params)

        # Color blind friendly palette - based on Wong 2011 palette
        self.default_colors = [
            "#009E73",  # Bluish Green
            "#98FB98",  # Pale Green
            "#56B4E9",  # Sky Blue
            "#F0E442",  # Yellow
            "#E69F00",  # Orange
            "#CC79A7",  # Reddish Purple
            "#D55E00",  # Vermillion
            "#FF6347",  # Tomato
            "#FFB6C1",  # Light Pink
            "#DDA0DD",  # Plum
            "#0072B2",  # Blue
            "#F0E68C",  # Khaki
        ]

        self.library = ChemLibrary(self.lib_select)
        
        # Update fine threshold options based on loaded library
        self.update_fine_threshold_options()
        
        self.color_widgets = {}
        self.color_collapse = None
        self.visible_columns = ["Compound", str(self.fine_threshold)]
        # SuperCluster will be added to visible columns after it's created
        self.compound_table = pn.widgets.Tabulator(
            self.library.df[self.visible_columns],
            sizing_mode="stretch_both",
            selectable="checkbox",  # Enable multiple row selection
            disabled=True,
            visible=False,
            show_index=False,
            pagination=None,
            margin=0,
            styles={"padding": "0px"},
        )

        # Slider widget takes advantage of param for real time updating, but uses param.watch
        # with throttling to update only when slider stops for more complicated operations
        # e.g. buliding a graph plot
        self.slider_widget = pn.widgets.DiscreteSlider.from_param(
            self.param.cluster_slider, name="Super Cluster", disabled=True
        )
        self.slider_widget.param.watch(self.update_throttled, "value_throttled")
        self.compound_table.on_click(self.on_table_click)
        self.compound_table.param.watch(self.on_table_selection, "selection")

        # Network of sub-clusters in super cluster
        self.cluster_graph = pn.pane.HoloViews(
            object=None, sizing_mode="stretch_both", margin=0
        )
        # Scatter plot of super cluster sizes (in terms of sub-clusters)
        # TODO: make it number of compounds instead?
        self.cluster_chart = pn.pane.HoloViews(
            object=None, sizing_mode="stretch_width", margin=0
        )
        self.compound_grid = Carrousel()
        self.compound_image = pn.pane.SVG(
            object=None, name="Compound Image", sizing_mode="stretch_both"
        )
        self.category_histogram = pn.pane.HoloViews(
            object=None, sizing_mode="stretch_both", margin=0
        )
        self.histogram_selection = None  # Will be initialized when histogram is created
        self.compound_data_chart = pn.pane.HoloViews(
            object=None, sizing_mode="stretch_both", margin=0
        )
        self.common_substructure = None

    # Update graph view and table after slider release
    def update_throttled(self, event):
        """This function updates only when the slider stops moving (throttled)."""
        # Get the first cluster node from the super cluster
        # TODO: get rid of use of super_clusters attribute and
        # reference dataframe instead?
        cluster_nodes = self.library.super_clusters[event.new - 1][2]

        member_cluster = cluster_nodes[0]
        self.common_substructure = None
        self.library.build_subgraph(member_cluster)
        self.initialize_graph_plot(self.library.sub_graph)
        self.graph_view()
        self.compound_table.value = self.build_table(
            event.new
        )  # Use super cluster number directly
        self.selected_compound = ""
        self.compound_image.object = None
        # Update the category histogram for the new super cluster
        if hasattr(self, "color_dict"):
            self.category_histogram.object = self.create_category_histogram()
        # self.compound_grid.svgs=[]
        # clear selection in table
        self.compound_table.selection = []

    @param.depends("cluster_slider", watch=True)
    def update_realtime(self):
        """This function updates in real-time as the slider moves."""
        # Find the compound count for the current super cluster
        super_cluster_data = self.library.super_clusters[self.cluster_slider - 1]
        compound_count = super_cluster_data[1]
        chart = self.draw_cluster_chart()
        self.cluster_chart.object = chart * hv.Scatter(
            [(self.cluster_slider, compound_count)]
        ).opts(color="red", size=12)

    @param.depends("fine_threshold", watch=True)
    def update_fine_threshold(self):
        self.visible_columns = ["Compound", "SMILES", str(self.fine_threshold)]
        # Remove SuperCluster from visible columns since it will be recreated
        if "SuperCluster" in self.visible_columns:
            self.visible_columns.remove("SuperCluster")
        self.compound_table.visible = False
        self.param.trigger("toggle_miss_button")  # Refresh button panel visibility
        self.compound_image.object = None
        self.compound_grid.svgs = []
        self.cluster_graph.object = None
        self.cluster_chart.object = None
        self.slider_widget.disabled = True

    @param.depends("toggle_miss_button", watch=True)
    def toggle_miss_compounds(self):
        """Update table display when Miss compounds toggle button is pressed."""
        # Toggle the boolean value
        self.show_miss_compounds = not self.show_miss_compounds

        # Update button label based on current state
        if self.show_miss_compounds:
            self.param.toggle_miss_button.label = "Hide Misses"
        else:
            self.param.toggle_miss_button.label = "Show Misses"

        if self.compound_table.visible:
            # Preserve the current selection
            selected_compounds = self.compound_table.value.loc[
                self.compound_table.selection, "Compound"
            ].values

            if self.selected_nodes:
                # If specific nodes are selected, rebuild table for those nodes
                new_table = self.build_table(self.selected_nodes)
                self.compound_table.value = new_table
            else:
                # Otherwise rebuild table for current super cluster
                new_table = self.build_table(self.slider_widget.value)
                self.compound_table.value = new_table

            # Reapply the selection to the table
            self.compound_table.selection = new_table[
                new_table["Compound"].isin(selected_compounds)
            ].index.tolist()

            # Update the compound grid to reflect the filtered compounds
            if (
                not self.compound_table.value.empty
                and "Compound" in self.compound_table.value.columns
            ):
                compounds = self.compound_table.value["Compound"].values
                self.update_compound_grid(compounds)
            else:
                self.update_compound_grid([])

            # Re-apply table styling
            self.style_compound_table()

    @param.depends("retest_selected_button", watch=True)
    def retest_selected_compounds(self):
        """Toggle Retest status for all selected compounds in the table."""
        if not self.compound_table.selection:
            return

        # Get selected row indices
        selected_indices = self.compound_table.selection
        selected_compounds = self.compound_table.value.loc[
            selected_indices, "Compound"
        ].values

        compound_mask = self.library.df["Compound"].isin(selected_compounds)
        self.library.df.loc[compound_mask, "Retest"] = ~self.library.df.loc[
            compound_mask, "Retest"
        ]

        # update the dataset dataframe
        dataset_mask = self.library.dataset_df["Compound"].isin(selected_compounds)
        self.library.dataset_df.loc[
            dataset_mask, "Retest"
        ] = ~self.library.dataset_df.loc[dataset_mask, "Retest"]

        # update the subset dataframe
        subset_mask = self.library.subset_df["Compound"].isin(selected_compounds)
        self.library.subset_df.loc[subset_mask, "Retest"] = ~self.library.subset_df.loc[
            subset_mask, "Retest"
        ]

        # Also update dataset_df with current clustering and Retest information
        if hasattr(self.library, 'dataset_df') and hasattr(self.library, 'super_clusters'):
            self.library.update_dataset_with_clustering(self.fine_threshold, self.coarse_threshold)

        # Refresh the table to show updated values
        if self.selected_nodes:
            new_table = self.build_table(self.selected_nodes)
        else:
            new_table = self.build_table(self.slider_widget.value)

        self.compound_table.value = new_table
        self.style_compound_table()

        # Restore the selection
        self.compound_table.selection = selected_indices

        print(f"Toggled Retest status for {len(selected_compounds)} compounds")

    @param.depends("lib_select", watch=True)
    def update_compound_df(self):
        # BUG: Changing library results in key error due to missing columns in self.visible_columns
        self.library = ChemLibrary(self.lib_select)
        
        # Update fine threshold options based on new library
        self.update_fine_threshold_options()
        
        self.compound_table.value = self.library.df[self.visible_columns]
        self.compound_input = ""
        self.selected_compound = ""
        self.compound_table.visible = False
        self.param.trigger("toggle_miss_button")  # Refresh button panel visibility
        self.compound_image.object = None

    @param.depends("dataset_select", watch=True)
    def load_dataset_df(self):
        """Load dataset and create subset_df from unique compounds with categories."""

        if not self.dataset_select:
            return

        # Load the dataset
        self.library.load_dataset(self.dataset_select)
        self.visible_columns = self.visible_columns + ["Category", "Retest"]

        # Create color pickers for categories found in the subset
        if (
            hasattr(self.library, "subset_df")
            and "Category" in self.library.subset_df.columns
        ):
            unique_categories = self.library.subset_df["Category"].dropna().unique()

            self.categories = sorted(unique_categories)
            self.create_color_widgets(self.categories)
            self.color_widgets_visible = True

        # Clear any existing clustering results from UI
        self.compound_table.visible = False
        self.param.trigger("toggle_miss_button")  # Refresh button panel visibility
        self.compound_image.object = None
        self.compound_grid.svgs = []
        self.cluster_graph.object = None
        self.cluster_chart.object = None
        self.slider_widget.disabled = True
        self.compound_data_chart.object = None

    @param.depends("recluster_button", watch=True)
    def recluster_library(self):
        if not hasattr(self.library, "subset_df"):
            # Necessary if dataset_select has only a single dataset or if the first
            # dataset in the list is selected and wont trigger load_dataset_df
            # probably should fix this...someday...
            self.load_dataset_df()
            

        self.slider_widget.disabled = True
        # ensure coarse threshold is greater than fine threshold
        if self.coarse_threshold <= self.fine_threshold:
            self.coarse_threshold = np.round(self.fine_threshold + 0.1, 2)
        # check to see if coarse_threshold is already in library.subdf
        if str(self.coarse_threshold) not in self.library.subset_df.columns:
            self.library.cluster_subset(self.coarse_threshold)
            # previously self.library.subset_df.to_csv(self.subset_select, index=False)
            # save clustering results back to dataset
            # implement some version of below
            # self.library.save_dataset(self.dataset_select)
        self.library.build_graph(self.fine_threshold, self.coarse_threshold)

        # Update dataset_df with clustering information and current Retest values
        self.library.update_dataset_with_clustering(self.fine_threshold, self.coarse_threshold)

        self.param.cluster_slider.objects = list(
            range(1, len(self.library.super_clusters) + 1)
        )
        self.slider_widget.value = 1
        self.slider_widget.disabled = False
        self.update_cluster_colors()
        # draw the super cluster scatter plot
        chart = self.draw_cluster_chart()

        # highlight the first cluster
        first_super_cluster = self.library.super_clusters[0]
        self.cluster_chart.object = chart * hv.Scatter([first_super_cluster]).opts(
            color="red", size=12
        )

        # Get the first cluster node from the first super cluster
        cluster_nodes = self.library.super_clusters[0][2]
        member_cluster = cluster_nodes[0]
        self.library.build_subgraph(member_cluster)
        self.initialize_graph_plot(self.library.sub_graph)
        self.graph_view()

        # Add SuperCluster to visible columns if not already there
        if "SuperCluster" not in self.visible_columns:
            self.visible_columns.insert(-2, "SuperCluster")
        # self.visible_columns = self.visible_columns + ["Category", "Retest"]
        new_table = self.build_table(1)  # First super cluster
        self.compound_table.value = new_table
        self.style_compound_table()
        self.compound_table.visible = True
        self.param.trigger("toggle_miss_button")  # Refresh button panel visibility

        compounds = new_table["Compound"].values
        self.update_compound_grid(compounds)

        # Initialize the category histogram for the first super cluster
        if hasattr(self, "color_dict"):
            self.category_histogram.object = self.create_category_histogram()

    def get_available_fine_thresholds(self):
        """Get available fine threshold values from the library dataframe."""
        if not hasattr(self.library, 'df') or self.library.df is None:
            return [0.2]  # Default threshold if no data loaded
        
        available_thresholds = []
        for column in self.library.df.columns:
            try:
                # Check if column name can be converted to float (clustering threshold)
                threshold_value = float(column)
                # Only include reasonable threshold values (between 0.05 and 0.8)
                if 0.05 <= threshold_value <= 0.8:
                    available_thresholds.append(threshold_value)
            except ValueError:
                continue
        
        # Sort thresholds and return, or default if none found
        if available_thresholds:
            return sorted(available_thresholds)
        else:
            return [0.2]  # Default threshold if no clustering columns found

    def update_fine_threshold_options(self):
        """Update the fine threshold selector options based on available data."""
        available_thresholds = self.get_available_fine_thresholds()
        
        # Update the parameter objects
        self.param.fine_threshold.objects = available_thresholds
        
        # Set default to first available threshold if current value is not available
        if self.fine_threshold not in available_thresholds:
            self.fine_threshold = available_thresholds[0]

    def build_table(self, clusters_or_super_cluster):
        """
        Build table from either a super cluster number or list of cluster nodes.

        Args:
            clusters_or_super_cluster: Either:
                - int/float: Super cluster number (1-indexed) to show all compounds in that super cluster
                - list: List of cluster node IDs to show compounds from specific selected clusters

        Returns:
            pd.DataFrame: Filtered dataframe with compounds from specified clusters
        """
        table_df = self.library.df

        # Check if input is a single integer (super cluster number)
        if isinstance(clusters_or_super_cluster, (int, float)):
            # Filter by super cluster number using the SuperCluster column
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
                        table_df[str(self.fine_threshold)].isin(cluster_list)
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
                table_df[str(self.fine_threshold)].isin(clusters_or_super_cluster)
            ]

        # Filter out "Miss" compounds if toggle is disabled
        if not self.show_miss_compounds:
            table_df = table_df[table_df["Category"] != "Miss"]

        table_df = table_df[self.visible_columns].reset_index()

        return table_df

    @param.depends("search_button", watch=True)
    def search_compound(self):
        compounds = re.findall(r"\d+", self.compound_input)
        if compounds:
            # NOTE: cheesy way to force search for single compounds
            # not great for usability FIX.
            compound = [int(comp) for comp in compounds][0]
            try:
                # Get cluster information for the compound
                compound_data = self.library.get_compounds(compound)
                if compound_data.empty:
                    self.compound_input = "Compound not found in library, check ID"
                    return

                cluster = compound_data[str(self.fine_threshold)].values[0]

                # Check if we have SuperCluster column for more efficient lookup
                if "SuperCluster" in self.library.df.columns:
                    super_cluster_data = self.library.df[
                        self.library.df["Compound"] == str(compound)
                    ]["SuperCluster"]
                    if not super_cluster_data.empty:
                        super_cluster = int(super_cluster_data.iloc[0])
                    else:
                        super_cluster = None
                else:
                    # Fallback to old method
                    super_cluster = None
                    for super_cluster_info in self.library.super_clusters:
                        super_cluster_number = super_cluster_info[0]
                        cluster_nodes = super_cluster_info[2]
                        if cluster in cluster_nodes:
                            super_cluster = super_cluster_number
                            break

                if super_cluster:
                    self.slider_widget.value = super_cluster
                    # TODO: Below seems like a janky fix...but it works.
                    event = MockEvent(super_cluster)
                    self.update_throttled(event)
                    # find the index of the cluster containing the compound, update the selection
                    index = list(np.argwhere(self.node_labels == cluster)[0])
                    self.selection.update(index=index)
                    # find the index of the compound in the table, update the selection
                    index = self.compound_table.value[
                        self.compound_table.value["Compound"] == str(compound)
                    ].index[0]
                    self.compound_table.selection = [
                        int(index)
                    ]  # highlights the row in table, but doesnt click
                    event = MockEvent(row=int(index))
                    # simulate click on table row
                    self.on_table_click(event)
                else:
                    # TODO: Everything below here is a temporary fix to handle compounds that are in the hit list
                    # but do not form super clusters (e.g. singletons, orphans, etc.)
                    # note that the display of the graph and slider will not change in this case
                    # up above  there was no if statement
                    # NOTE: --autoreload breaks displaying on svg file for some reason

                    # update table
                    new_table = self.build_table(
                        [cluster]
                    )  # Pass as list for single cluster
                    self.compound_table.value = new_table
                    self.selected_nodes = []

                    # draw grid for compounds in this cluster
                    compounds = new_table["Compound"].values
                    self.update_compound_grid(compounds)

            except IndexError:
                # TODO: Should return default looking text, not input style text
                self.compound_input = "Compound not found in library, check ID"

    @param.depends("save_button", watch=True)
    def save_spreadsheet(self):
        # NOTE: see link (1) at bottom of this file if you want to implment saving
        # the compound image to the excel file
        self.compound_table.value.to_excel("compound_list.xlsx", index=False)

    @param.depends("save_dataset_button", watch=True)
    def save_dataset(self):
        """Save the current dataset with all modifications back to the original file."""
        if not hasattr(self.library, "dataset_df") or self.dataset_select is None:
            print("No dataset loaded to save")
            return
        
        try:
            # Ensure dataset is up to date with latest clustering and Retest info
            if hasattr(self.library, 'super_clusters'):
                self.library.update_dataset_with_clustering(self.fine_threshold, self.coarse_threshold)
            
            # Save the dataset_df back to the original file
            self.library.dataset_df.to_csv(self.dataset_select, index=False)

            
        except Exception as e:
            print(f"Error saving dataset: {e}")

    def on_table_click(self, event):
        compound = self.compound_table.value.loc[event.row, "Compound"]
        image = self.library.draw_compound(
            compound,
            common_substructure=self.common_substructure,
            legend=False,
        )
        self.selected_compound = f"### Compound ID: {compound}"
        self.compound_image.object = image
        # add the row to the selection
        self.compound_table.selection = [event.row]  # Highlight the clicked row

    def on_table_selection(self, event):
        """Handle multiple row selection in the compound table."""
        if not event.new:
            # No selection - clear everything
            self.compound_image.object = None
            self.compound_data_chart.object = None
            return
        # Get selected row indices
        selected_indices = event.new  # always a list
        selected_compounds = self.compound_table.value.loc[
            selected_indices, "Compound"
        ].values

        # Update compound grid with selected compounds
        self.update_compound_grid(selected_compounds)
        # Update compound data chart with selected compounds
        # self.compound_data_chart.object = self.create_compound_data_scatter(
        #     selected_compounds
        # )
        self.compound_data_chart.object = self.create_compound_data_chart(
            selected_compounds
        )

        # For single selection, also update the compound image
        if len(selected_compounds) == 1:
            image = self.library.draw_compound(
                selected_compounds[0],
                common_substructure=self.common_substructure,
                legend=False,
            )
            self.selected_compound = f"### Compound ID: {selected_compounds[0]}"
            self.compound_image.object = image
        else:
            # Multiple selections - clear single compound display
            self.selected_compound = f"### {len(selected_compounds)} compounds selected"
            self.compound_image.object = None

    def initialize_graph_plot(self, G):
        self.G = G

        # Generate positions for nodes using networkx
        self.pos = nx.spring_layout(self.G)  # Positions of nodes

        # Convert node positions and edges to lists of coordinates for HoloViews
        self.node_positions = np.array([self.pos[n] for n in self.G.nodes()])
        self.node_labels = list(self.G.nodes())
        self.edges_data = [
            (self.pos[edge[0]], self.pos[edge[1]]) for edge in self.G.edges()
        ]

        # use node_labels and self.library.df['Node Color'] to create a color map
        cmap = [self.cluster_color_map[node] for node in self.node_labels]

        # Extract x and y coordinates for edges
        self.edges_array = np.array(
            [(x1, y1, x2, y2) for ((x1, y1), (x2, y2)) in self.edges_data]
        )
        data = pd.DataFrame(
            {
                "x": self.node_positions[:, 0],
                "y": self.node_positions[:, 1],
                "color": cmap,
            }
        )

        # Create HoloViews Points for the nodes
        self.points = hv.Points(data, ["x", "y"]).opts(
            size=15,
            min_width=300,
            min_height=300,
            responsive=True,
            xaxis=None,
            yaxis=None,
            color="color",
            tools=["tap", "box_select", "lasso_select"],
            active_tools=["tap"],
        )
        # create cmap by mapping self.library.subset_df to with color_dict

        self.labels = hv.Labels(
            {("x", "y"): self.node_positions, "text": self.node_labels},
            ["x", "y"],
            "text",
        ).opts(text_color="black", text_font_size="10pt")
        # Create HoloViews Segments for the edges
        self.initial_edges = hv.Segments(self.edges_array).opts(
            line_width=0.5, color="gray"
        )

        self.selection = Selection1D(source=self.points)

        # Watch for node selection changes
        self.selection.param.watch(self.update_selection, "index")

        # Select all nodes by default
        self.selection.update(index=list(range(len(self.node_labels))))

    def update_selection(self, event):
        """Update the graph plot, compound plots and table when a node is selected."""
        # clear compound core structure from previous selection
        self.common_substructure = None
        indices = self.selection.index
        if indices:
            self.selected_nodes = [self.node_labels[i] for i in indices]

            # update table
            new_table = self.build_table(self.selected_nodes)
            self.compound_table.value = new_table
            self.style_compound_table()
            # draw grid
            compounds = new_table["Compound"].values
            self.update_compound_grid(compounds)
            # clear selected compound
            self.selected_compound = ""
            self.compound_image.object = None
            self.compound_data_chart.object = None

        else:
            # node has been deselected, clear everything and reset table
            # Note: This can also go below in graph_view? Perhaps here is better
            # because it is more explicit for the action being taken.
            self.update_compound_grid([])  # Clear the grid
            self.compound_image.object = None
            self.compound_data_chart.object = None
            # Reset table to show current super cluster
            self.compound_table.value = self.build_table(self.slider_widget.value)
            self.style_compound_table()
            self.selected_nodes = []

    def on_histogram_click(self, event):
        """Handle clicks on histogram bars to select all subclusters with compounds of that category."""
        if not event.new or not hasattr(self, "library"):
            return

        try:
            # Get the index of the clicked bar

            clicked_index = event.new[0] if event.new else None

            if clicked_index is None:
                return

            # Get all available categories (sorted same way as histogram)
            all_categories = (
                list(self.color_widgets.keys()) if self.color_widgets else []
            )
            all_categories = sorted(all_categories, key=category_sort_key)

            if clicked_index >= len(all_categories):
                return

            # Get the clicked category
            clicked_category = all_categories[clicked_index]

            # Find all compounds in current super cluster
            cluster_compounds = self.get_current_super_cluster_compounds()

            # Get compounds of the clicked category
            category_compounds = self.library.subset_df[
                (self.library.subset_df["Compound"].isin(cluster_compounds))
                & (self.library.subset_df["Category"] == clicked_category)
            ]["Compound"].tolist()

            if not category_compounds:
                return

            # Find which subclusters contain these compounds
            subclusters_with_category = set()
            for compound in category_compounds:
                compound_data = self.library.df[self.library.df["Compound"] == compound]
                if not compound_data.empty:
                    subcluster = compound_data[str(self.fine_threshold)].iloc[0]
                    subclusters_with_category.add(subcluster)

            # Find the node indices for these subclusters
            node_indices = []
            for i, node_label in enumerate(self.node_labels):
                if node_label in subclusters_with_category:
                    node_indices.append(i)

            # Update the graph selection
            if node_indices:
                self.selection.update(index=node_indices)

        except Exception as e:
            print(f"Error in histogram click handler: {e}")

    @param.depends("selected_nodes", watch=True)
    def graph_view(self):
        """Return the graph plot with updated segments and edges."""
        if len(self.selected_nodes) > 0:
            # get selected nodes
            n_nodes = len(self.node_positions)
            selected_nodes = self.points.data.loc[self.selection.index]
            non_selected_index = [
                i for i in range(n_nodes) if i not in self.selection.index
            ]
            non_selected_nodes = self.points.data.loc[non_selected_index]

            selected_points = hv.Points(selected_nodes).opts(
                size=30, color="color", line_color="k", tools=["tap", "box_select"]
            )
            non_selected_points = hv.Points(non_selected_nodes).opts(
                size=15, color="color"
            )

            # Get all connected edges for the selected node
            connected_edges = []
            non_connected_edges = []
            for edge in self.G.edges():
                if any([node in edge for node in self.selected_nodes]):
                    x1, y1 = self.pos[edge[0]]
                    x2, y2 = self.pos[edge[1]]
                    connected_edges.append((x1, y1, x2, y2))
                else:
                    x1, y1 = self.pos[edge[0]]
                    x2, y2 = self.pos[edge[1]]
                    non_connected_edges.append((x1, y1, x2, y2))

            connected_segments = hv.Segments(connected_edges).opts(
                line_width=2, color="black"
            )
            # Faded segments (e.g., red with reduced opacity)
            faded_segments = hv.Segments(non_connected_edges).opts(
                line_width=1, color="gray", alpha=0.1
            )

            self.cluster_graph.object = (
                connected_segments
                * faded_segments
                * self.points
                * selected_points
                * non_selected_points
                * self.labels
            )
        else:
            # Reset to default view (all edges normal)
            self.cluster_graph.object = self.initial_edges * self.points * self.labels

    def create_category_histogram(self):
        """Create a histogram showing category distribution for the current super cluster."""
        if not hasattr(self, "library") or not hasattr(self.library, "subset_df"):
            return hv.Text(0.5, 0.5, "No subset data available for histogram.").opts(
                width=600, height=250, xaxis=None, yaxis=None
            )

        try:
            # Get compounds in current super cluster
            cluster_compounds = self.get_current_super_cluster_compounds()

            # Get categories for these compounds
            subset_data = self.library.subset_df[
                self.library.subset_df["Compound"].isin(cluster_compounds)
            ]

            if subset_data.empty:
                return hv.Text(
                    0.5, 0.5, "No category data for current super cluster."
                ).opts(width=600, height=250, xaxis=None, yaxis=None)

            # Get all available categories from color widgets to ensure completeness
            all_categories = (
                list(self.color_widgets.keys()) if self.color_widgets else []
            )

            all_categories = sorted(all_categories, key=category_sort_key)

            # Count categories in current super cluster
            category_counts = subset_data[subset_data["Category"] != "Miss"][
                "Category"
            ].value_counts()

            # Create complete DataFrame including all categories (even with 0 counts)
            hist_data = pd.DataFrame(
                {
                    "Category": all_categories,
                    "Count": [category_counts.get(cat, 0) for cat in all_categories],
                    "%Total": [
                        str(
                            int(
                                np.round(
                                    100
                                    * category_counts.get(cat, 0)
                                    / len(cluster_compounds)
                                )
                            )
                        )
                        + "%"
                        for cat in all_categories
                    ],
                }
            )

            # Map colors from color widgets
            colors = []
            for category in hist_data["Category"]:
                colors.append(self.color_dict.get(category, "#999999"))

            hist_data["Color"] = colors
            # Create HoloViews bar chart
            bars = hv.Bars(hist_data, ["Category"], ["Count", "%Total", "Color"]).opts(
                color="Color",
                min_width=300,
                min_height=300,
                responsive=True,
                title=f"Category Distribution - Super Cluster {self.slider_widget.value}",
                ylabel="Count",
                xlabel="",
                xrotation=45,
                tools=["hover", "tap"],  # Add tap tool for interactivity
                active_tools=["tap"],
                ylim=(
                    0,
                    1.2 * hist_data["Count"].max(),
                ),  # Set y-axis limits to auto scale
            )

            # Set up selection stream for histogram bar clicks
            # Recreate the selection stream each time to attach it to the new bars object
            self.histogram_selection = Selection1D(source=bars)
            self.histogram_selection.param.watch(self.on_histogram_click, "index")

            return bars

        except Exception as e:
            return hv.Text(0.5, 0.5, f"Error creating histogram: {str(e)}").opts(
                width=600, height=250, xaxis=None, yaxis=None
            )

    def update_category_histogram(self):
        """Update the category histogram with current data."""
        if hasattr(self, "category_histogram"):
            self.category_histogram.object = self.create_category_histogram()

    def create_color_widgets(self, categories):
        """Create color picker widgets for each category."""
        self.color_widgets = {}

        categories = sorted(categories, key=category_sort_key)

        for i, category in enumerate(categories):
            # Assign default color from palette
            default_color = self.default_colors[i % len(self.default_colors)]

            # Create color picker widget
            color_widget = pn.widgets.ColorPicker(
                name=category, value=default_color, width=120, height=50, margin=(5, 10)
            )

            # Watch for color changes
            color_widget.param.watch(self.on_color_change, "value")
            self.color_widgets[category] = color_widget

    def on_color_change(self, event):
        """Handle color picker changes and update visualization."""
        if hasattr(self, "library") and hasattr(self.library, "subset_df"):
            self.update_cluster_colors()
            if hasattr(self, "cluster_color_map"):
                self.refresh_visualizations()
            # Update the histogram with new colors
            if hasattr(self, "category_histogram"):
                self.category_histogram.object = self.create_category_histogram()
            # Update the compound grid with new colors - get compounds directly
            if hasattr(self, "selected_nodes") and self.selected_nodes:
                current_table = self.compound_table.value
                if not current_table.empty and "Compound" in current_table.columns:
                    compounds = current_table["Compound"].values
                    self.update_compound_grid(compounds)
            # Update the compound data chart with new colors
            if hasattr(self, "compound_data_chart") and self.compound_data_chart.object is not None:
                # Get currently selected compounds from the table
                if hasattr(self, "compound_table") and self.compound_table.selection:
                    selected_indices = self.compound_table.selection
                    selected_compounds = self.compound_table.value.loc[selected_indices, "Compound"].values
                    self.compound_data_chart.object = self.create_compound_data_chart(selected_compounds)
            # style the compound table with the new colors
            if hasattr(self, "compound_table"):
                self.compound_table.value = (
                    self.compound_table.value
                )  # Trigger re-render
                self.style_compound_table()

    def update_cluster_colors(self):
        """Update cluster colors based on current color picker values."""
        if not hasattr(self.library, "subset_df"):
            return

        # Create color dictionary from current widget values
        self.color_dict = {
            category: widget.value for category, widget in self.color_widgets.items()
        }

        fine_threshold = str(self.fine_threshold)

        # Create a mapping of clusters to their category counts
        categories = self.library.df[self.library.df["Category"] != "Miss"]
        cluster_category_counts = (
            categories.groupby(fine_threshold)["Category"]
            .value_counts()
            .unstack(fill_value=0)
        )

        self.cluster_category_counts = cluster_category_counts

        # Apply the color determination function to each cluster
        self.cluster_color_map = {
            cluster: self.determine_cluster_color(counts)
            for cluster, counts in cluster_category_counts.iterrows()
        }

    def style_compound_table(self):
        """Apply styles to the compound table based on cluster colors."""
        if not hasattr(self, "compound_table") or not hasattr(
            self, "cluster_color_map"
        ):
            return

        # Define a function to style each row based on category
        def style_row(row):
            color = self.color_dict.get(row.Category, "#FFFFFF")  # Default to white
            return [f"background-color: {color};"] * len(row)

        self.compound_table.style.apply(style_row, axis=1)

    def determine_cluster_color(self, category_counts):
        """Determine cluster color based on the category with the most compounds."""
        try:
            # Find the category with the maximum count
            if category_counts.sum() == 0:
                return "#999999"  # Default color if no compounds

            max_category = category_counts.idxmax()
            return self.color_dict.get(max_category, "#999999")

        except (AttributeError, TypeError, ValueError):
            return "#999999"  # Default color if no categories found

    def refresh_visualizations(self):
        """Refresh graph and other visualizations with new colors."""
        if hasattr(self, "G") and hasattr(self, "cluster_color_map"):
            # Update graph colors
            cmap = [
                self.cluster_color_map.get(node, "#999999") for node in self.node_labels
            ]

            # Update points data
            data = pd.DataFrame(
                {
                    "x": self.node_positions[:, 0],
                    "y": self.node_positions[:, 1],
                    "color": cmap,
                }
            )

            # Recreate points with new colors
            self.points = hv.Points(data, ["x", "y"]).opts(
                size=15,
                min_width=300,
                min_height=300,
                responsive=True,
                xaxis=None,
                yaxis=None,
                color="color",
                tools=["tap", "box_select", "lasso_select"],
                active_tools=["tap"],
            )

            # Update selection source
            self.selection.source = self.points

            # Refresh the graph view
            self.graph_view()

    def update_compound_grid(self, compounds):
        """Update the compound grid display with given compounds."""
        if len(compounds) == 0:
            self.compound_grid.svgs = []
            self.common_substructure = None
            return

        # Determine grid parameters based on number of compounds
        if len(compounds) > 1:
            mols_per_row = 4
            max_rows = 4
            orient = True
        else:
            mols_per_row = 4
            max_rows = 1
            orient = False

        # Generate the grid with current colors
        grid_image, self.common_substructure = self.library.draw_compound_grid(
            compounds,
            mols_per_row=mols_per_row,
            max_rows=max_rows,
            color_dict=self.color_dict,
            orient=orient,
            legend=False,
        )

        # Update the compound grid display
        self.compound_grid.svgs = grid_image

    def get_color_widgets_panel(self):
        """Create a collapsible panel containing all color widgets."""
        if not self.color_widgets:
            return pn.pane.Markdown("*No categories loaded. Load a subset file first.*")

        # Create rows of color widgets (2 per row for better layout)
        widget_rows = []
        widgets_list = list(self.color_widgets.values())

        for i in range(0, len(widgets_list), 2):
            row_widgets = widgets_list[i : i + 2]
            widget_rows.append(pn.Row(*row_widgets))

        return pn.Column(*widget_rows, margin=(10, 5))

    def get_current_super_cluster_compounds(self):
        """Get list of compounds in the current super cluster."""
        current_super_cluster_num = self.slider_widget.value

        return self.library.df[
            self.library.df["SuperCluster"] == current_super_cluster_num
        ]["Compound"].tolist()

    def create_compound_data_chart(self, compound_ids, metric="Delta Lifetime Z"):
        """Create a bar chart showing Delta Lifetime values for one or more compounds across constructs."""

        compound_df = self.library.dataset_df[
            self.library.dataset_df["Compound"].isin(compound_ids)
        ]

        if compound_df.empty:
            return hv.Text(0.5, 0.5, "No data available for selected compounds").opts(
                width=600, height=400
            )

        stats_df = (
            compound_df.groupby(["Compound", "Construct"])[metric]
            .agg(["mean", "std"])
            .reset_index()
            .rename(columns={"mean": "Mean", "std": "Std"})
        )

        stats_df = stats_df.merge(
            compound_df[["Compound", "Category"]].drop_duplicates(),
            on="Compound",
            how="left",
        )
        # make a Color column in compound_df by using self.color_dict to map
        stats_df["Color"] = stats_df["Category"].map(self.color_dict)
        stats_df.loc[stats_df["Color"] == "none", "Color"] = (
            "#999999"  # Default color for missing categories
        )

        # Create title based on number of compounds
        if len(compound_ids) == 1:
            title = f"Delta Lifetime Z for Compound {compound_ids[0]}"
        else:
            title = f"Delta Lifetime Z for {len(compound_ids)} Compounds"

        # Calculate symmetrical y-axis limits based on data
        max_abs_value = max(abs(stats_df["Mean"].min()), abs(stats_df["Mean"].max()))
        # Add some padding and ensure we can see ±4 lines
        y_limit = max(max_abs_value * 1.1, 4.5)

        # Create bar chart
        bars = hv.Bars(
            stats_df,
            kdims=["Construct", "Compound"],
            vdims=["Mean", "Category", "Color"],
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
            legend_position="right",
        )

        # Add horizontal reference lines at ±4
        hline_pos4 = hv.HLine(4).opts(
            color="red", line_dash="dashed", line_width=2, alpha=0.7
        )
        hline_neg4 = hv.HLine(-4).opts(
            color="red", line_dash="dashed", line_width=2, alpha=0.7
        )

        # Combine bars with reference lines
        chart = bars * hline_pos4 * hline_neg4
        return chart

    def create_compound_data_scatter(self, compound_ids, metric="Delta Lifetime Z"):
        """Create a scatter plot showing Delta Lifetime values for one or more compounds across constructs."""

        compound_df = self.library.dataset_df[
            self.library.dataset_df["Compound"].isin(compound_ids)
        ]

        if compound_df.empty:
            return hv.Text(0.5, 0.5, "No data available for selected compounds").opts(
                width=600, height=400
            )

        # Color palette for constructs
        color_palette = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]

        # Create DataFrames
        stats_df = (
            compound_df.groupby(["Compound", "Construct"])[metric]
            .agg(["mean", "std"])
            .reset_index()
            .rename(columns={"mean": "Mean", "std": "Std"})
        )

        # Create title based on number of compounds
        if len(compound_ids) == 1:
            title = f"Delta Lifetime Z for Compound {compound_ids[0]}"
        else:
            title = f"Delta Lifetime Z for {len(compound_ids)} Compounds"

        # Create bar chart with error bars
        scatter = hv.Scatter(
            stats_df, kdims=["Compound"], vdims=["Mean", "Construct"]
        ).opts(
            color="Construct",
            cmap=color_palette,
            size=8,
            width=600,
            height=400,
            title=title,
            ylabel=metric,
            xlabel="Construct",
            ylim=(-7, 7),
            xrotation=45,
            tools=["hover"],
            show_grid=False,
            show_legend=False,
        )

        return scatter

    def draw_cluster_chart(self):
        """
        Draws hv scatter plot of supercluster sizes with hover and tap tools
        """
        data = [
            (super_clust_id, compound_count)
            for super_clust_id, compound_count, _ in self.library.super_clusters
        ]
        # Use the new super_clusters format which already contains [super_cluster_number, compound_count]
        plot = hv.Scatter(data).opts(
            tools=["hover", "tap"],
            width=300,
            height=200,
            size=10,
            xlabel="Super Cluster",
            ylabel="# Compounds",
        )

        return plot

    @param.depends("toggle_miss_button")
    def get_table_controls_panel(self):
        """Create a panel with table control buttons, only visible when table is visible."""
        if not self.compound_table.visible:
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


class MockEvent:
    def __init__(self, new=None, row=None):
        self.new = new
        self.row = row


clusterF = ClusterF()

sidebar = pn.Column(
    pn.Param(
        clusterF.param,
        parameters=[
            "lib_select",
            "dataset_select",
            "fine_threshold",
            "coarse_threshold",
            "recluster_button",
        ],
        widgets={
            "lib_select": pn.widgets.RadioButtonGroup,
            "dataset_select": {"name": "Select Dataset", "value": ""},
        },
        show_name=False,
    ),
    # Color picker controls in a collapsible widget
    pn.Card(
        clusterF.get_color_widgets_panel,
        title="Category Colors",
        collapsed=True,
        visible=pn.bind(lambda x: x, clusterF.param.color_widgets_visible),
        margin=(5, 2),
    ),
    clusterF.cluster_chart,
    clusterF.slider_widget,
    pn.Param(
        clusterF.param,
        parameters=["compound_input"],
        widgets={
            "compound_input": {
                "type": pn.widgets.TextAreaInput,
                "placeholder": "Search for Super Cluster containing a specific compound ID",
                "name": "Compound ID Search",
            },
        },
        show_name=False,
    ),
    pn.Param(
        clusterF.param,
        parameters=["search_button", "save_button", "save_dataset_button"],
        default_layout=pn.Row,
        margin=(2, 2),
        show_name=False,
    ),
    pn.pane.Markdown(clusterF.param.selected_compound, margin=(0, 5)),
    clusterF.compound_image,
    margin=0,
    styles={"padding": "0px"},
)

main = pn.GridSpec(
    nrows=3,
    ncols=2,
    sizing_mode="stretch_both",
    margin=0,
    name="Main Layout",
    styles={"gap": "2px", "padding": "2px", "box-sizing": "border-box"},
)

# Layout:
# Row 0, Col 0: Network graph
main[0, 0] = clusterF.cluster_graph
# Row 0, Col 1: Compound grid (carousel)
main[0:2, 1] = clusterF.compound_grid.view()
# Row 1, Col 0: Category histogram
main[1, 0] = clusterF.category_histogram
# Row 2, Col 1: Compound table (Tabulator)
main[2, 0] = pn.Row(
    clusterF.get_table_controls_panel,
    clusterF.compound_table,
    sizing_mode="stretch_width",
    margin=0,
)
# Row 0, Col 2: Compound lifetime chart
main[2, 1] = clusterF.compound_data_chart


pn.template.FastListTemplate(
    site="ClusterF",
    title="A ChemBridge Compound Viewer",
    sidebar=sidebar,
    main=main,
).servable()


# References:
# (1) https://stackoverflow.com/questions/39563065/how-to-insert-an-image-in-an-excel-sheet-using-openpyxl
