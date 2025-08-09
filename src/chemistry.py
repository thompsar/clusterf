import warnings
import pandas as pd
import numpy as np
import networkx as nx
import holoviews as hv
from holoviews import opts
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, AllChem, Draw
from rdkit.Chem import rdFMCS

# Below needed for some reason to get the draw functions to work
# even though it is not called. If not imported, draw_compound
# will throw  AttributeError: 'str' object has no attribute 'data'
# when img.data is called.
# Its not sufficient to change img.data to img.
from rdkit.Chem.Draw import IPythonConsole # noqa: F401


class ChemLibrary:
    """
    Class for loading info about ChemBridge compound library
    """

    def __init__(self, library_name):
        self.library_name = library_name.lower()
        self.load_library()

    def load_library(self):
        """
        Loads chemical library file
        """
        path = "compound_libraries/"
        self.library_file = path + self.library_name + ".csv"
        try:
            self.df = pd.read_csv(self.library_file)
            self.df = self.standardize_df(self.df)

        except FileNotFoundError:
            raise FileNotFoundError(
                "Library not found, supported libraries are cns and diverset"
            )

    def load_dataset(self, path):
        """
        Loads dataset containing potentially multiple copies of compounds across different constructs.
        Creates both dataset_df (full dataset) and subset_df (unique compounds with categories).
        Expects dataset to have Category column pre-populated, including "Miss".
        Supports both CSV and Parquet file formats.
        """
        if path.lower().endswith(".parquet"):
            self.dataset_df = pd.read_parquet(path)
        else:
            self.dataset_df = pd.read_csv(path, dtype={"Sub Categories": str})
        # drop super cluster column if it exists
        # since dataset will have na for some values and dont want to deal with those for now.
        if "SuperCluster" in self.dataset_df.columns:
            self.dataset_df = self.dataset_df.drop(columns=["SuperCluster"])

        self.dataset_df = self.standardize_df(self.dataset_df)
        if "Retest" not in self.dataset_df.columns:
            self.dataset_df["Retest"] = False

        # Create subset_df with unique compounds and their categories
        self.create_subset_df()

        # Merge categories into main library df
        self.df = self.df.merge(
            self.subset_df[["Compound", "Category", "Retest"]],
            on="Compound",
            how="left",
        )
        self.df["Category"] = self.df["Category"].fillna("Miss")
        # Note below, converting to 'boolean' and not 'bool' was necessary to
        # avoid FutureWarning: Downcasting object dtype arrays on .fillna, .ffill, .bfill is deprecated and
        # will change in a future version. Call result.infer_objects(copy=False) instead.
        # To opt-in to the future behavior, set `pd.set_option('future.no_silent_downcasting', True)`
        # except below line also throws FutureWarning:
        # self.df["Retest"] = self.df["Retest"].fillna(False).astype(bool).infer_objects(copy=False)
        self.df["Retest"] = self.df["Retest"].fillna(False)
        self.df["Retest"] = self.df["Retest"].astype("boolean").fillna(False)

        # save for later debugging
        # print('retest column value counts:')
        # print(self.df['Retest'].value_counts())

    def create_subset_df(self):
        """
        Creates subset_df from dataset_df with unique compounds and consolidated categories.
        If clustering columns exist in dataset_df, they will be preserved in subset_df.
        """
        if not hasattr(self, "dataset_df"):
            raise ValueError("dataset_df not loaded. Call load_dataset first.")

        # Get clustering columns that might already exist in dataset_df
        clustering_columns = []
        for column in self.dataset_df.columns:
            try:
                # Check if column name can be converted to float (clustering threshold)
                float(column)
                clustering_columns.append(column)
            except ValueError:
                continue
        columns = ["Compound"] + clustering_columns + ["Category", "Retest"]
        # Create subset_df with unique compounds
        self.subset_df = (
            self.dataset_df[self.dataset_df.Category != "Miss"][columns]
            .drop_duplicates(subset="Compound")
            .reset_index(drop=True)
        )
        self.subset_df["Retest"] = (
            self.subset_df["Retest"].astype("boolean").fillna(False)
        )

        self.subset_df = self.standardize_df(self.subset_df)

    def standardize_df(self, df):
        """
        Standardizes the columns of a dataframe
        """
        # Convert Compound column to string in both
        df.Compound = df.Compound.astype(str)
        # Need to fix this, since we wont necessarily have all the columns in the df
        for column in np.arange(0.05, 0.8, 0.05):
            col = str(np.round(column, 2))
            if col in df.columns:
                # Handle NaN values before converting to int
                df[col] = df[col].fillna(-1).astype(int)
                df[col] = df[col].replace(-1, np.nan)
        return df

    def cluster_subset(self, coarse_thresh, radius=2, fpSize=2048):
        # check to see if subset_df already has SMILES strings
        if "SMILES" not in self.subset_df.columns:
            subset_smiles = self.df[
                self.df.Compound.isin(self.subset_df.Compound.values)
            ][["Compound", "SMILES"]]
            # merge subset smiles with subset_df
            self.subset_df = subset_smiles.merge(self.subset_df, on="Compound")
        mols = [MolFromSmiles(smiles) for smiles in self.subset_df.SMILES.values]
        # see https://greglandrum.github.io/rdkit-blog/posts/2023-01-18-fingerprint-generator-tutorial.html
        # for more info on finger printing
        fpgen = AllChem.GetRDKitFPGenerator()
        # use morgan fingerprints instead? Corresponds to ECFP4
        # fpgen = AllChem.GetMorganGenerator(radius=2,fpSize=2048)
        #ECFP6 is radius=3
        #fpgen = AllChem.GetMorganGenerator(radius=radius,fpSize=fpSize)
        # fpgen = AllChem.GetMorganGenerator(radius=2,fpSize=2048,atomInvariantsGenerator=AllChem.GetMorganFeatureAtomInvGen())

        # fpgen = AllChem.GetTopologicalTorsionGenerator()
        # see also https://www.researchgate.net/post/Two_similar_compounds_with_a_low_smaller_than_085_Tanimoto_coefficient2
        # for more info on how fingerprinting might not capture two near identical looking molecules.
        # different type of fingerprint can be used. see singleton clusters 3950,6353 from CNS
        fps = [fpgen.GetFingerprint(mol) for mol in mols]
        
    
        clusters = ClusterFps(fps, coarse_thresh)

        for idx, cluster in enumerate(clusters, start=1):
            self.subset_df.loc[cluster, str(np.round(coarse_thresh, 2))] = int(idx)
            # all_smiles.Cluster = all_smiles.Cluster.astype(int)
        self.subset_df = self.standardize_df(self.subset_df)

    def build_graph(self, fine_thresh, coarse_thresh):
        """
        Builds graph of subclusters
        """
        # convert fine_thresh and coarse_thresh to string if not already
        fine_thresh = str(fine_thresh)
        coarse_thresh = str(coarse_thresh)

        # Precompute cluster to compound and compound to cluster mappings
        fine_clusters = self.df.groupby(fine_thresh).Compound.apply(list).to_dict()
        fine_compounds = self.df.set_index("Compound")[fine_thresh].to_dict()
        coarse_clusters = (
            self.subset_df.groupby(coarse_thresh).Compound.apply(list).to_dict()
        )
        coarse_compounds = self.subset_df.set_index("Compound")[coarse_thresh].to_dict()

        # Initialize edge array
        dims = len(fine_clusters)
        edge_array = np.zeros((dims + 1, dims + 1))

        # Build edge array
        for cluster, comps in fine_clusters.items():
            coarse_mapping = set(
                [coarse_compounds[comp] for comp in comps if comp in coarse_compounds]
            )
            if len(coarse_mapping) > 0:
                compound_mapping = set(
                    [
                        comp
                        for coarse_clust in coarse_mapping
                        for comp in coarse_clusters[coarse_clust]
                    ]
                )
                related_clusters = list(
                    set([fine_compounds[comp] for comp in compound_mapping])
                )
                edge_array[cluster, related_clusters] = 1
        # find rows that only have 1 element
        self.singletons = np.where(edge_array.sum(axis=1) == 1)[0]
        # set diagonal to 0
        np.fill_diagonal(edge_array, 0)
        # add back singleton nodes to edge_array
        edge_array[self.singletons, self.singletons] = 1
        self.graph = nx.from_numpy_array(edge_array)
        # drop 0 node from graph
        self.graph.remove_node(0)

        # Get superclusters
        # subgraph are sets, convert to list
        super_cluster_groups = [
            list(subgraph)
            for subgraph in nx.connected_components(self.graph)
            if len(subgraph) > 1
        ]

        # add singletons to super_cluster_groups
        super_cluster_groups.extend([[node] for node in self.singletons])

        super_cluster_counts = []
        for subgraph in super_cluster_groups:
            mask = self.df[fine_thresh].isin(subgraph)
            compound_count = len(self.df[mask]["Compound"])
            super_cluster_counts.append([compound_count, subgraph])

        # sort super_cluster_counts by compound count
        super_cluster_counts.sort(key=lambda x: x[0], reverse=True)

        # Map each cluster to its super cluster number (1-indexed) and count compounds
        self.super_clusters = [
            [super_clust_id, compound_count, sub_clusters]
            for super_clust_id, (compound_count, sub_clusters) in enumerate(
                super_cluster_counts, start=1
            )
        ]

        # Initialize super cluster column with NaN
        self.df["SuperCluster"] = np.nan
        # Map each cluster to its super cluster number (1-indexed)
        for super_clust_id, compound_count, sub_clusters in self.super_clusters:
            self.df.loc[self.df[fine_thresh].isin(sub_clusters), "SuperCluster"] = (
                super_clust_id
            )
        print(
            f"Total compounds in super clusters: {np.sum([comp for _, comp, __ in self.super_clusters])}"
        )

    def extract_sub_categories(self):
        """
        Extract sub-categories from the dataset_df into individual columns with symbols.
        This processes the 'Sub Categories' column if it exists.
        """
        if (
            not hasattr(self, "dataset_df")
            or "Sub Categories" not in self.dataset_df.columns
        ):
            return {}

        # Get unique compounds and their sub-categories (convert to string for hashing)
        temp_df = self.dataset_df[["Compound", "Sub Categories"]].copy()
        temp_df["Sub Categories"] = temp_df["Sub Categories"].astype(str)
        unique_compounds = temp_df.drop_duplicates()

        # Parse sub-categories and get all unique keys
        all_keys = set()
        parsed_subcats = {}

        for _, row in unique_compounds.iterrows():
            compound = row["Compound"]
            subcats_raw = row["Sub Categories"]

            if pd.isna(subcats_raw) or subcats_raw == "nan":
                parsed_subcats[compound] = {}
                continue

            try:
                # Parse the string representation of the dictionary
                subcats_dict = eval(str(subcats_raw))
                if isinstance(subcats_dict, dict):
                    parsed_subcats[compound] = subcats_dict
                    all_keys.update(subcats_dict.keys())
                else:
                    parsed_subcats[compound] = {}
            except (SyntaxError, NameError, TypeError, ValueError):
                parsed_subcats[compound] = {}

        # Convert values to symbols
        def value_to_symbol(value):
            if value == "+":
                return "↑"
            elif value == "-":
                return "↓"
            elif value == "Interfering":
                return "!"
            else:
                return "—"  # em dash for None or other values

        # Create columns for each sub-category key
        subcategory_columns = {}
        for key in sorted(all_keys):
            column_data = {}
            for compound in self.df["Compound"]:
                compound_subcats = parsed_subcats.get(compound, {})
                column_data[compound] = value_to_symbol(compound_subcats.get(key))
            subcategory_columns[key] = column_data

        return subcategory_columns

    def update_dataset_with_clustering(self, fine_thresh, coarse_thresh):
        """
        Update dataset_df with clustering information and current Retest values.
        This should be called after clustering is complete.
        """
        if not hasattr(self, "dataset_df"):
            return
        # Get clustering and category information from the main library df
        cluster_info = self.df[["Compound", str(fine_thresh), "Retest"]].copy()

        # Remove existing clustering columns from dataset_df if they exist
        columns_to_drop = []
        if str(fine_thresh) in self.dataset_df.columns:
            columns_to_drop.append(str(fine_thresh))
        if str(coarse_thresh) in self.dataset_df.columns:
            columns_to_drop.append(str(coarse_thresh))
        if "Retest" in self.dataset_df.columns:
            columns_to_drop.append("Retest")

        if columns_to_drop:
            self.dataset_df = self.dataset_df.drop(columns=columns_to_drop)

        # Merge clustering information into dataset_df
        self.dataset_df = self.dataset_df.merge(
            cluster_info,
            on="Compound",
            how="left",
        )

        # # Add coarse clustering information from subset_df if it exists
        # if hasattr(self, 'subset_df') and str(coarse_thresh) in self.subset_df.columns:
        #     coarse_info = self.subset_df[['Compound', str(coarse_thresh)]].copy()
        #     self.dataset_df = self.dataset_df.merge(
        #         coarse_info,
        #         on='Compound',
        #         how='left'
        # )

    def compute_saturation_metrics(self, fine_thresh):
        """
        Compute and attach saturation metrics to subset_df:
        - SCS (Super Cluster Saturation): for each super cluster, fraction of non-miss compounds

        Metric is merged back onto subset_df per Compound so it can be saved/exported.
        """
        if not hasattr(self, "df") or not hasattr(self, "subset_df"):
            return

        fine_col = str(fine_thresh)
        # Ensure required columns exist
        required_cols = {"Compound", fine_col, "Category"}
        if not required_cols.issubset(set(self.df.columns)):
            return

        # If SuperCluster hasn't been computed yet, skip gracefully
        has_super = "SuperCluster" in self.df.columns

        df_main = self.df[["Compound", fine_col, "Category"]].copy()
        if has_super:
            df_main["SuperCluster"] = self.df["SuperCluster"]

        # Helper to compute saturation: non-miss / total
        def _saturation(group: pd.DataFrame) -> float:
            total = len(group)
            if total == 0:
                return np.nan
            non_miss = (group["Category"] != "Miss").sum()
            return float(non_miss) / float(total)


        # Super cluster saturation (SCS)
        if has_super:
            scs_series = df_main.dropna(subset=["SuperCluster"]).groupby("SuperCluster").apply(_saturation)
            scs_map = scs_series.to_dict()
        else:
            scs_map = {}

        # Map metrics back per compound using df_main (compound -> cluster ids)
        df_metrics = df_main[["Compound", fine_col]].copy()

        if has_super:
            df_metrics["SCS"] = df_main["SuperCluster"].map(scs_map)
        else:
            df_metrics["SCS"] = np.nan

        df_metrics = df_metrics.drop(columns=[fine_col])

        # Merge into subset_df per compound
        self.subset_df = self.subset_df.merge(df_metrics, on="Compound", how="left")


    def build_subgraph(self, member_cluster):
        """
        Builds subgraph of related clusters
        """

        for sub_graph in nx.connected_components(self.graph):
            if member_cluster in sub_graph:
                self.sub_graph = self.graph.subgraph(sub_graph).copy()
                break
        self.node_pos = nx.spring_layout(self.sub_graph, iterations=100)

    def draw_subgraph(self, member_cluster):
        """
        Draws subgraph of related clusters
        """
        cb_colors = [
            (230, 159, 0),
            (86, 180, 233),
            (0, 158, 115),
            (240, 228, 66),
            (0, 114, 178),
            (213, 94, 0),
            (204, 121, 167),
        ]
        cb_colors = ["#%02x%02x%02x" % color for color in cb_colors]

        self.build_subgraph(member_cluster)

        kwargs = dict(width=800, height=800, xaxis=None, yaxis=None)
        opts.defaults(opts.Nodes(**kwargs), opts.Graph(**kwargs))

        plot = hv.Graph.from_networkx(
            self.sub_graph,
            positions=self.node_pos,
        ).opts(
            tools=["tap", "lasso_select"],
            node_color="type",
            edge_alpha=1,
            node_alpha=0.9,
            node_size=50,
        )
        self.nodes = plot.nodes
        labels = hv.Labels(plot.nodes, ["x", "y"], "index")
        self.graph_plot = plot * labels.opts(text_font_size="12pt", text_color="black")
        return self.graph_plot

    def get_compounds(self, compound_ids):
        """Gets compound(s) info from library dataframe

        Parameters
        ----------
        compound_id : list of str, int or single str, int
            Chembridge compound ID(s)

        Returns
        -------
        if list:
        pd.DataFrame: chem info including compound name, SMILES string, other details and molfile
        if str or int:
        pd.Series: chem info including compound name, SMILES string, other details and molfile
        """
        # convert compound_id to string
        if not isinstance(compound_ids, (list, np.ndarray)):
            compound_ids = [compound_ids]
        compound_ids = [str(compound_id) for compound_id in compound_ids]

        subset = self.df[self.df["Compound"].isin(compound_ids)]

        # check for missing compounds
        missing = set(compound_ids) - set(subset["Compound"])
        if len(missing) > 0:
            msg = "Compound(s) " + ", ".join(missing) + " not found in library"
            warnings.warn(msg)
        return subset

    def draw_compound_grid(
        self,
        compound_ids,
        mols_per_row=6,
        max_rows=3,
        common_substructure=None,
        color_dict=None,
        orient=False,
        transparent=False,
        legend=True,
    ):
        """Draws Chembridge compound(s)

        Parameters
        ----------
        compound_ids : list of str, int, or single str, int
            Chembridge compound IDs
        mols_per_row : int, optional
            number of molecules drawn per row, by default 6
        max_rows : int, optional
            maximum number of rows, by default 3
        transparent : bool, optional
            draw grid with transparent background, by default False
        legend : bool, optional
            include compound id as legend, by default True

        Returns
        -------
        image : SVG of compound or compounds in grid
        """

        img_size = 300  # Smaller size for better performance in web apps
        chem_info = self.get_compounds(compound_ids)
        # redefine compound_ids below to handle input of bad ID
        compound_ids = chem_info["Compound"].values
        all_categories = chem_info["Category"].values
        all_mols = [
            Chem.MolFromSmiles(smiles_str) for smiles_str in chem_info["SMILES"]
        ]

        if orient:
            all_mols, highlight_atoms, common_substructure = orient_mols(
                all_mols, return_pattern=True
            )

            highlight_colors = [
                {idx: (86 / 255, 180 / 255, 233 / 255, 0.5) for idx in highlight}
                for highlight in highlight_atoms
            ]
        else:
            highlight_atoms = None
            highlight_colors = None

        ngrids = int(np.ceil(len(all_mols) / (mols_per_row * max_rows)))
        mols_per_grid = mols_per_row * max_rows

        def raggedify(var, mols_per_grid):
            if var is None:
                return [None] * ngrids
            else:
                return [
                    var[i : i + mols_per_grid]
                    for i in range(0, len(var), mols_per_grid)
                ]

        vars = [
            all_mols,
            compound_ids,
            highlight_atoms,
            highlight_colors,
            all_categories,
        ]

        vars = map(lambda x: raggedify(x, mols_per_grid), vars)

        imgs = []
        for mols, ids, atoms, colors, categories in zip(*vars):
            # Render the molecules in SVG format
            # legend below used to be [str(id) for id in ids],
            img = Draw.MolsToGridImage(
                mols,
                molsPerRow=mols_per_row,
                subImgSize=(img_size, img_size),
                highlightAtomLists=atoms,
                highlightAtomColors=colors,
                legends=[
                    str(id) + "\n\n" + str(cat) if cat != "Miss" else str(id)
                    for id, cat in zip(ids, categories)
                ],
                useSVG=True,  # Use SVG rendering for efficiency
            )

            # add gridlines (make into function later)
            delimiter = "</rect>\n"
            head, tail = img.data.split(delimiter)
            # head, tail = img.split(delimiter)
            head = head + delimiter
            std_bg_style = "opacity:1.0;fill:#FFFFFF;stroke:none"

            # note stroke-width is 4 since gridlines are 2 wide
            # but double up with neighboring cells
            if transparent:
                new_bg_style = "opacity:0.0;fill:none;stroke:black;stroke-width:4"
                head = head.replace(std_bg_style, new_bg_style)
            else:
                new_bg_style = "opacity:1.0;fill:#FFFFFF;stroke:black;stroke-width:4"
                head = head.replace(std_bg_style, new_bg_style)

            if color_dict is None:
                color_dict = {category: "none" for category in set(categories)}
            # override color_dict miss to none
            color_dict["Miss"] = "none"
            grid_fill = [
                (
                    f'<rect width="{img_size}" '
                    f'height="{img_size}" '
                    f'x="{x}" y="{y}" '
                    'style="fill:none;'
                    "fill-opacity:0.3;"
                    "stroke:black;"
                    "stroke-width:2;"
                    'stroke-opacity:1"/>'
                )
                for y in range(0, img_size * (len(mols) // mols_per_row + 1), img_size)
                for x in range(0, img_size * mols_per_row, img_size)
            ]
            # below requires the grid_fill be draw row wise as is done above with y,x in for loop
            # i dont really like this, but its less clunky than trying to cram it in above
            # find a better way!
            grid_fill[: len(categories)] = [
                line.replace("fill:none", f"fill:{color_dict[category]}")
                for category, line in zip(categories, grid_fill[: len(categories)])
            ]

            grid_fill = "\n".join(grid_fill)
            img.data = head + grid_fill + tail
            # img = head + grid_fill + tail
            imgs.append(img)
        if orient:
            return imgs, common_substructure
        else:
            return imgs, None

    def draw_compound(
        self,
        compound_id,
        common_substructure=None,
        orient=False,
        transparent=False,
        legend=True,
    ):
        """Draws Chembridge compound

        Parameters
        ----------
        compound_ids : list of str, int, or single str, int
            Chembridge compound IDs
        mols_per_row : int, optional
            number of molecules drawn per row, by default 6
        max_rows : int, optional
            maximum number of rows, by default 3
        transparent : bool, optional
            draw grid with transparent background, by default False
        legend : bool, optional
            include compound id as legend, by default True

        Returns
        -------
        image : SVG of compound or compounds in grid
        """

        img_size = 300  # Smaller size for better performance in web apps
        chem_info = self.get_compounds(compound_id)
        mols = [Chem.MolFromSmiles(smiles_str) for smiles_str in chem_info["SMILES"]]

        if common_substructure is not None:
            mols, highlight_atoms = orient_mols(
                mols, common_substructure, return_pattern=False
            )
        # There is a reason I'm using grid image below but I don't remember why
        # perhaps due to useSVG?
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=1,
            subImgSize=(img_size, img_size),
            useSVG=True,  # Use SVG rendering for efficiency
        )

        if orient:
            return img, common_substructure
        else:
            return img


def ClusterFps(fps, cutoff=0.2):
    # adapted from https://github.com/tsudalab/ChemGE/blob/master/results/clustering.py
    # https://www.tsudalab.org
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina

    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs


def orient_mols(mols, common_substructure=None, return_pattern=False):
    # Find the Maximum Common Substructure (MCS)
    # https://greglandrum.github.io/rdkit-blog/posts/2023-10-27-mcswhatsnew.html
    params = rdFMCS.MCSParameters()
    params.BondCompareParameters.RingMatchesRingOnly = True
    params.BondCompareParameters.CompleteRingsOnly = True
    # params.StoreAll = False

    if common_substructure is None:
        mcs_result = rdFMCS.FindMCS(mols, params)
        common_substructure = Chem.MolFromSmarts(mcs_result.smartsString)

    # Ensure 2D coordinates for the common substructure
    AllChem.Compute2DCoords(common_substructure)
    # Align molecules to the common substructure
    for mol in mols:
        AllChem.Compute2DCoords(mol)  # Ensure the molecule has 2D coordinates
        AllChem.GenerateDepictionMatching2DStructure(mol, common_substructure)

    # Highlight atoms involved in the MCS
    highlight_atoms = [mol.GetSubstructMatch(common_substructure) for mol in mols]

    if return_pattern:
        return mols, highlight_atoms, common_substructure
    else:
        return mols, highlight_atoms

    # References:
    # (1) https://projects.volkamerlab.org/teachopencadd/talktorials/T005_compound_clustering.html
