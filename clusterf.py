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
- [ ] Add plots accoriding to values in the table (e.g. MW vs LogP)
    - [ ] Plot values should be highlighted to to reveal compounds of interest. Interesting patterns may appear.
- [x] Convert compound ID text box into "compound search" box, have that adjust to appropraite cluster/graph
    - [x] Currently only seeks super cluster. Add selection of exact cluster within super cluster.
- [ ] Add a "reset" button to reset the graph to the original view
- [x] Results of reclustering should be saved to the CSV of the subset file
- [ ] Cluster hover should provide basic stats about cluster (number of compounds per category, total compounds, etc.)
- [ ] Subset select should have a drop down menu that shows the subset categories, with color pickers for colorizing the clusters
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
- [ ] Colorize table based on category
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

"""

def category_sort_key(category):
                category_str = str(category).lower()
                if 'hit' in category_str:
                    return (0, category)  # Highest priority
                elif 'questionable' in category_str:
                    return (1, category)  # Second priority
                elif 'interfering' in category_str:
                    return (2, category)  # Third priority
                else:
                    return (3, category)  # Lowest priority (includes "Miss" and others)
    

class ClusterF(param.Parameterized):
    libraries = ["cns", "cns_kinase"]  # 'diverset': ChemLibrary('diverset')}

    # Widgets
    lib_select = param.Selector(objects=libraries, default=libraries[1])
    fine_threshold = param.Number(0.2)
    coarse_threshold = param.Number(0.4)
    subset_select = param.FileSelector(path=os.path.join("compound_subsets", "*.csv*"))
    compound_input = param.String()
    selected_compound = param.String()
    selected_nodes = param.List()  # Track the selected node
    selected_cluster = param.List()
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
        self.color_widgets = {}
        self.color_collapse = None
        self.visible_columns = ["Compound", "SMILES", str(self.fine_threshold)]
        self.compound_table = pn.widgets.Tabulator(
            self.library.df[self.visible_columns],
            height=400,
            width=800,
            selectable=1,
            disabled=True,
            visible=False,
            show_index=False,
            pagination=None,
        )

        # Slider widget takes advantage of param for real time updating, but uses param.watch
        # with throttling to update only when slider stops for more complicated operations
        # e.g. buliding a graph plot
        self.slider_widget = pn.widgets.DiscreteSlider.from_param(
            self.param.cluster_slider, name="Super Cluster", disabled=True
        )
        self.slider_widget.param.watch(self.update_throttled, "value_throttled")
        self.compound_table.on_click(self.on_click)
        self.cluster_chart = pn.pane.HoloViews(object=None)
        self.cluster_graph = pn.pane.HoloViews(object=None)
        self.compound_grid = Carrousel()
        self.compound_image = pn.pane.SVG(object=None, width=300, name="Compound Image")
        self.category_histogram = pn.pane.HoloViews(
            object=None, sizing_mode="stretch_width"
        )
        self.common_substructure = None

    # Update graph view and table after slider release
    def update_throttled(self, event):
        """This function updates only when the slider stops moving (throttled)."""

        cluster_set = self.library.super_clusters[event.new - 1]
        member_cluster = cluster_set[0]  # TODO: fix this
        self.common_substructure = None
        self.library.build_subgraph(member_cluster)
        self.initialize_graph_plot(self.library.sub_graph)
        self.graph_view()
        self.compound_table.value = self.build_table(cluster_set)
        self.selected_compound = ""
        self.compound_image.object = None
        # Update the category histogram for the new super cluster
        if hasattr(self, "color_dict"):
            self.category_histogram.object = self.create_category_histogram()
        # self.compound_grid.svgs=[]

    @param.depends("cluster_slider", watch=True)
    def update_realtime(self):
        """This function updates in real-time as the slider moves."""
        # add a red point to the cluster chart plot corresponding to the value of self.cluster_slider
        cluster_sizes = self.library.cluster_sizes
        chart = self.library.cluster_chart
        self.cluster_chart.object = chart * hv.Scatter(
            cluster_sizes[self.cluster_slider - 1]
        ).opts(color="red", size=12)

    @param.depends("fine_threshold", watch=True)
    def update_fine_threshold(self):
        self.visible_columns = ["Compound", "SMILES", str(self.fine_threshold)]
        self.compound_table.visible = False
        self.compound_image.object = None
        self.compound_grid.svgs = []
        self.cluster_graph.object = None
        self.cluster_chart.object = None
        self.slider_widget.disabled = True

    @param.depends("lib_select", watch=True)
    def update_compound_df(self):
        self.library = ChemLibrary(self.lib_select)
        self.compound_table.value = self.library.df[self.visible_columns]
        self.compound_input = ""
        self.selected_compound = ""
        self.compound_table.visible = False
        self.compound_image.object = None

    @param.depends("subset_select", watch=True)
    def load_subset_df(self):
        # BUG: this only works if there are multiple csv files in the directory
        # fixing temporarily with if hasattr below in recluster_library
        self.library.load_subset(self.subset_select)

        # Create color pickers for categories found in the subset
        if (
            hasattr(self.library, "subset_df")
            and "Category" in self.library.subset_df.columns
        ):
            unique_categories = self.library.subset_df["Category"].dropna().unique()

            self.categories = sorted(unique_categories)
            self.create_color_widgets(self.categories)
            self.color_widgets_visible = True

    @param.depends("recluster_button", watch=True)
    def recluster_library(self):
        if not hasattr(self.library, "subset_df"):
            self.library.load_subset(self.subset_select)

        self.slider_widget.disabled = True
        # ensure coarse threshold is greater than fine threshold
        if self.coarse_threshold <= self.fine_threshold:
            self.coarse_threshold = np.round(self.fine_threshold + 0.1, 2)
        # check to see if coarse_threshold is already in library.subdf
        if str(self.coarse_threshold) not in self.library.subset_df.columns:
            self.library.cluster_subset(self.coarse_threshold)
            # save clustering results
            self.library.subset_df.to_csv(self.subset_select, index=False)
        self.library.build_graph(self.fine_threshold, self.coarse_threshold)
        self.colorize_clusters()

        self.param.cluster_slider.objects = list(
            range(1, len(self.library.super_clusters) + 1)
        )
        self.slider_widget.value = 1
        self.slider_widget.disabled = False
        # draw the histogram
        chart = self.library.draw_cluster_chart()
        # highlight the first cluster
        cluster_sizes = self.library.cluster_sizes
        self.cluster_chart.object = chart * hv.Scatter(cluster_sizes[0]).opts(
            color="red", size=12
        )

        member_cluster = self.library.super_clusters[0][0]  # TODO: fix this
        self.library.build_subgraph(member_cluster)
        self.initialize_graph_plot(self.library.sub_graph)
        self.graph_view()

        self.compound_table.value = self.build_table(self.library.super_clusters[0])
        self.compound_table.visible = True

        # Initialize the category histogram for the first super cluster
        if hasattr(self, "color_dict"):
            self.category_histogram.object = self.create_category_histogram()

    def build_table(self, clusters):
        if not isinstance(clusters, list):
            clusters = [clusters]
        table_df = self.library.df
        table_df = table_df[table_df[str(self.fine_threshold)].isin(list(clusters))]
        table_df = table_df[self.visible_columns].reset_index()
        return table_df

    # below currently does nothing!!!
    @param.depends("selected_cluster", watch=True)
    def update_table(self):
        nodes = self.library.nodes.data.loc[self.selected_cluster, "index"].values
        new_df = self.library.df
        new_df = new_df[new_df[self.fine_threshold].isin(nodes)][
            self.visible_columns
        ].reset_index()
        self.compound_table.value = new_df
        self.compound_table.visible = True

    @param.depends("search_button", watch=True)
    def search_compound(self):
        compounds = re.findall(r"\d+", self.compound_input)
        if compounds:
            # NOTE: cheesy way to force search for single compounds
            # not great for usability FIX.
            compound = [int(comp) for comp in compounds][0]
            try:
                cluster = self.library.get_compounds(compound)[
                    str(self.fine_threshold)
                ].values[0]
                # super_cluster = [
                #     idx
                #     for idx, super_cluster in enumerate(self.library.super_clusters)
                #     if cluster in super_cluster
                # ][0] + 1
                super_cluster = [
                    idx
                    for idx, super_cluster in enumerate(self.library.super_clusters)
                    if cluster in super_cluster
                ]
                if super_cluster:
                    super_cluster = super_cluster[0] + 1
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
                    self.on_click(event)
                else:
                    # TODO: Everything below here is a temporary fix to handle compounds that are in the hit list
                    # but do not form super clusters (e.g. singletons, orphans, etc.)
                    # note that the display of the graph and slider will not change in this case
                    # up above  there was no if statement
                    # NOTE: --autoreload breaks displaying on svg file for some reason

                    # update table
                    new_table = self.build_table(cluster)
                    self.compound_table.value = new_table
                    self.selected_nodes = []

                    # code below was taken from update_selection to allow drawing of grid
                    compounds = new_table["Compound"].values

                    if len(compounds) > 1:
                        mols_per_row = 4
                        max_rows = 3
                        orient = True
                    else:
                        mols_per_row = 4
                        max_rows = 1
                        orient = False

                    grid_image, self.common_substructure = (
                        self.library.draw_compound_grid(
                            compounds,
                            mols_per_row=mols_per_row,
                            max_rows=max_rows,
                            color_dict=self.color_dict,
                            orient=orient,
                            legend=False,  # note this actually does nothing right now
                        )
                    )

                    self.compound_grid.svgs = grid_image

            except IndexError:
                # TODO: Should return default looking text, not input style text
                self.compound_input = "Compound not found in library, check ID"

    @param.depends("save_button", watch=True)
    def save_spreadsheet(self):
        # NOTE: see link (1) at bottom of this file if you want to implment saving
        # the compound image to the excel file
        self.compound_table.value.to_excel("compound_list.xlsx", index=False)

    def on_click(self, event):
        compound = self.compound_table.value.loc[event.row, "Compound"]
        image = self.library.draw_compound(
            compound,
            common_substructure=self.common_substructure,
            legend=False,
        )
        self.selected_compound = "### Compound ID: " + compound
        self.compound_image.object = image

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
            width=800,
            height=700,
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
        ).opts(text_color="black")
        # Create HoloViews Segments for the edges
        self.initial_edges = hv.Segments(self.edges_array).opts(
            line_width=1, color="gray"
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
            # draw grid
            compounds = new_table["Compound"].values

            if len(compounds) > 1:
                mols_per_row = 4
                max_rows = 3
                orient = True
            else:
                mols_per_row = 4
                max_rows = 1
                orient = False

            grid_image, self.common_substructure = self.library.draw_compound_grid(
                compounds,
                mols_per_row=mols_per_row,
                max_rows=max_rows,
                color_dict=self.color_dict,
                orient=orient,
                legend=False,  # note this actually does nothing right now
            )

            self.compound_grid.svgs = grid_image
            # clear selected compound
            self.selected_compound = ""
            self.compound_image.object = None

        else:
            # node has been deselected, clear everything and reset table
            # Note: This can also go below in graph_view? Perhaps here is better
            # because it is more explicit for the action being taken.
            # self.compound_grid.object = None
            self.compound_grid.svgs = []
            self.compound_image.object = None
            super_cluster = self.library.super_clusters[self.slider_widget.value - 1]
            self.compound_table.value = self.build_table(super_cluster)
            self.selected_nodes = []

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
                size=50, color="color", line_color="k", tools=["tap", "box_select"]
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

    def colorize_clusters(self):
        self.update_cluster_colors()
        # TODO: update grid colors
        #         
    def create_category_histogram(self):
        """Create a histogram showing category distribution for the current super cluster."""
        if not hasattr(self, "library") or not hasattr(self.library, "subset_df"):
            return hv.Text(0.5, 0.5, "No subset data available for histogram.").opts(
                width=600, height=250, xaxis=None, yaxis=None
            )

        try:
            # Get current super cluster
            current_super_cluster = self.library.super_clusters[
                self.slider_widget.value - 1
            ]

            # Get compounds in current super cluster
            cluster_compounds = self.library.df[
                self.library.df[str(self.fine_threshold)].isin(current_super_cluster)
            ]["Compound"].tolist()

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
            category_counts = subset_data["Category"].value_counts()
            # print(category_counts)

            # Create complete DataFrame including all categories (even with 0 counts)
            hist_data = pd.DataFrame(
                {
                    "Category": all_categories,
                    "Count": [category_counts.get(cat, 0) for cat in all_categories],
                    "%Total Comps": [
                        str(int(np.round(
                            100 * category_counts.get(cat, 0) / len(cluster_compounds)
                        )))+'%'
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
            bars = hv.Bars(
                hist_data, ["Category"], ["Count", "%Total Comps", "Color"]
            ).opts(
                color="Color",
                width=800,
                height=400,
                title=f"Category Distribution - Super Cluster {self.slider_widget.value}",
                xlabel="Category",
                ylabel="Count",
                xrotation=45,
                tools=["hover"],
                toolbar="above",
                ylim=(
                    0,
                    1.2 * hist_data[hist_data["Category"] != "Miss"]["Count"].max(),
                ),  # Set y-axis limits to auto scale
            )

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

    def update_cluster_colors(self):
        """Update cluster colors based on current color picker values."""
        if not hasattr(self.library, "subset_df"):
            return

        # Create color dictionary from current widget values
        self.color_dict = {
            category: widget.value for category, widget in self.color_widgets.items()
        }

        # Rebuild cluster color mapping
        category_df = self.library.df.merge(
            self.library.subset_df[["Compound", "Category"]], on="Compound", how="outer"
        )
        fine_threshold = str(self.fine_threshold)

        # Create a mapping of clusters to their category counts
        cluster_category_counts = category_df.groupby(fine_threshold)["Category"].value_counts().unstack(fill_value=0)
        self.cluster_category_counts = cluster_category_counts

        # Apply the color determination function to each cluster
        self.cluster_color_map = {
            cluster: self.determine_cluster_color(counts)
            for cluster, counts in cluster_category_counts.iterrows()
        }

        # Update histogram if it exists
        if hasattr(self, "category_histogram") and hasattr(self, "slider_widget"):
            self.category_histogram.object = self.create_category_histogram()

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
                width=800,
                height=700,
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
            "subset_select",
            "fine_threshold",
            "coarse_threshold",
            "recluster_button",
        ],
        widgets={
            "lib_select": pn.widgets.RadioButtonGroup,
            "subset_select": {"name": "Select Compound Subset", "value": ""},
        },
        show_name=False,
    ),
    clusterF.cluster_chart,
    clusterF.slider_widget,
    # Color picker controls in a collapsible widget
    pn.Card(
        clusterF.get_color_widgets_panel,
        title="Category Colors",
        collapsed=True,
        visible=pn.bind(lambda x: x, clusterF.param.color_widgets_visible),
        margin=(10, 5),
    ),
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
        parameters=["search_button", "save_button"],
        default_layout=pn.Row,
        margin=(-2, 5),
        show_name=False,
    ),
    pn.pane.Markdown(clusterF.param.selected_compound, margin=(0, 10)),
    clusterF.compound_image,
)
main = pn.Column(
    pn.Row(
        clusterF.cluster_graph,
        pn.Column(clusterF.compound_grid.view(), height=700, scroll=True),
        height=700,
    ),
    pn.Row(
        clusterF.category_histogram,
        clusterF.compound_table,
        height=400,
    ),
    width=1600,
)

pn.template.FastListTemplate(
    site="ClusterF",
    title="A ChemBridge Compound Viewer",
    sidebar=sidebar,
    main=main,
).servable()


# References:
# (1) https://stackoverflow.com/questions/39563065/how-to-insert-an-image-in-an-excel-sheet-using-openpyxl
