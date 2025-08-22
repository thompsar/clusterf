import os
import warnings
import pandas as pd
import numpy as np
import networkx as nx
from scipy.sparse import coo_matrix
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
from rdkit.Chem.Draw import IPythonConsole  # noqa: F401


class ChemLibrary:
    """
    Class for loading info about chemical compound library
    Supports parquet file format for library data
    Clustering information is stored in separate files
    """

    def __init__(
        self,
        library_name,
        fine_threshold=0.2,
        method="RDKit",
        libraries_dir="data/libraries",
    ):
        self.library_name = library_name
        self.fine_threshold = fine_threshold
        self.method = method
        self.libraries_dir = libraries_dir
        self.load_library()
        self.load_clustering_info()

    def load_library(self):
        """
        Loads chemical library file (Parquet format)
        Separates platemap from chemical data
        """
        library_path = os.path.join(self.libraries_dir, f"{self.library_name}.parquet")
        try:
            # Load the full library file
            full_library = pd.read_parquet(library_path)

            # Convert Compound column to string
            full_library.Compound = full_library.Compound.astype(str)

            # Separate platemap (all compounds including DMSO) from chemical data (compounds only)
            self.platemap = full_library[["Plate", "Col", "Row", "Compound"]].copy()

            # Create chemical data by filtering out DMSO and other non-compound entries
            # Keep only rows where Compound is not DMSO and has valid SMILES
            self.df = (
                full_library[(full_library.Compound != "DMSO")]
                .copy()
                .reset_index(drop=True)
            )

        except FileNotFoundError:
            raise FileNotFoundError(f"Library file not found: {library_path}")

    def load_clustering_info(self):
        """
        Loads clustering information from separate file
        Expected format: {library_name}_clusters.parquet with columns: Compound, Cluster, Method, Threshold
        """
        # Determine clustering file path based on library name
        clustering_path = os.path.join(
            self.libraries_dir, f"{self.library_name}_clusters.parquet"
        )

        try:
            self.clustering_df = pd.read_parquet(clustering_path)

            # Filter for the current method and threshold
            self.current_clustering = self.clustering_df[
                (self.clustering_df.Method == self.method)
                & (self.clustering_df.Threshold == self.fine_threshold)
            ]

            # Add clustering info to chemical data (self.df)
            # Remove existing Cluster column if it exists
            if "Cluster" in self.df.columns:
                self.df = self.df.drop(columns=["Cluster"])

            # Merge clustering info into chemical data
            self.df = self.df.merge(
                self.current_clustering[["Compound", "Cluster"]],
                on="Compound",
                how="left",
            )
            # ensure that the cluster column is an integer
            self.df["Cluster"] = self.df["Cluster"].astype(int)

        except FileNotFoundError:
            warnings.warn(f"Clustering file not found: {clustering_path}")
            self.clustering_df = None
            self.current_clustering = None
            # Add empty Cluster column to chemical data
            self.df["Cluster"] = np.nan

    def set_clustering_parameters(self, fine_threshold=None, method=None):
        """
        Update clustering parameters and reload clustering info
        """
        if fine_threshold is not None:
            self.fine_threshold = fine_threshold
        if method is not None:
            self.method = method

        # Reload clustering info with new parameters
        self.load_clustering_info()

    def get_available_clustering_parameters(self):
        """
        Get available clustering methods and thresholds from the clustering file

        Returns
        -------
        dict: Dictionary with 'methods' and 'thresholds' keys containing lists of available values
        """
        if self.clustering_df is None:
            return {"methods": [], "thresholds": []}

        methods = sorted(self.clustering_df.Method.unique().tolist())
        thresholds = sorted(self.clustering_df.Threshold.unique().tolist())

        return {"methods": methods, "thresholds": thresholds}

    def save_clustering_info(self, clustering_data, method, threshold):
        """
        Save clustering information to the clustering parquet file

        Parameters
        ----------
        clustering_data : pd.DataFrame
            DataFrame with columns: Compound, Cluster
        method : str
            Clustering method used
        threshold : float
            Clustering threshold used
        """
        # Determine clustering file path based on library name
        clustering_path = os.path.join(
            self.libraries_dir, f"{self.library_name}_clusters.parquet"
        )

        # Prepare new clustering data
        new_clustering = clustering_data.copy()
        new_clustering["Method"] = method
        new_clustering["Threshold"] = threshold

        # Convert Compound column to string
        new_clustering.Compound = new_clustering.Compound.astype(str)

        try:
            # Load existing clustering data if file exists
            if os.path.exists(clustering_path):
                existing_clustering = pd.read_parquet(clustering_path)
                existing_clustering.Compound = existing_clustering.Compound.astype(str)

                # Remove existing entries for this method and threshold
                existing_clustering = existing_clustering[
                    ~(
                        (existing_clustering.Method == method)
                        & (existing_clustering.Threshold == threshold)
                    )
                ]

                # Combine existing and new clustering data
                combined_clustering = pd.concat(
                    [existing_clustering, new_clustering], ignore_index=True
                )
            else:
                combined_clustering = new_clustering

            # Save to parquet file
            combined_clustering.to_parquet(clustering_path, index=False)

            # Reload clustering info
            self.load_clustering_info()

        except Exception as e:
            raise Exception(f"Failed to save clustering info: {str(e)}")

    def load_dataset(self, path):
        """
        Loads dataset containing potentially multiple copies of compounds across different constructs.
        Creates both dataset_df (full dataset) and subset_df (unique compounds with categories).
        Expects dataset to have Category column pre-populated, including "Miss".
        Supports both CSV and Parquet file formats.
        """
        # Backwards-compatible: treat as primary dataset load
        self.load_primary_dataset(path)

    def _read_dataset_file(self, path):
        """Read a dataset file and remove stale clustering columns."""
        if path.lower().endswith(".parquet"):
            df = pd.read_parquet(path)
        else:
            df = pd.read_csv(path, dtype={"Compound": str, "Sub Categories": str})
        if "SuperCluster" in df.columns:
            df = df.drop(columns=["SuperCluster"])
        return df

    def _standardize_dataset_df(self, df):
        if "Compound" in df.columns:
            df["Compound"] = df["Compound"].astype(str)
        if "Retest" not in df.columns:
            df["Retest"] = False
        try:
            df["Retest"] = df["Retest"].astype("boolean").fillna(False)
        except Exception:
            pass
        return df

    def _reset_df_category_retest_columns(self):
        if "Category" in self.df.columns:
            self.df = self.df.drop(columns=["Category"])
        if "Retest" in self.df.columns:
            self.df = self.df.drop(columns=["Retest"])

    def _merge_subset_into_df(self):
        subset_columns = ["Compound", "Category", "Retest"]
        existing_columns = [col for col in subset_columns if col in self.df.columns]
        new_columns = [col for col in subset_columns if col not in self.df.columns]
        if new_columns:
            self.df = self.df.merge(
                self.subset_df[["Compound"] + new_columns],
                on="Compound",
                how="left",
            )
        else:
            for col in existing_columns:
                if col in self.subset_df.columns:
                    try:
                        self.df[col] = self.df["Compound"].map(
                            self.subset_df.set_index("Compound")[col]
                        )
                    except KeyError:
                        continue
        if self.df.columns.duplicated().any():
            self.df = self.df.loc[:, ~self.df.columns.duplicated()]
        self.df["Category"] = self.df["Category"].fillna("Miss")
        try:
            self.df["Retest"] = self.df["Retest"].fillna(False)
            self.df["Retest"] = self.df["Retest"].astype("boolean").fillna(False)
        except Exception:
            pass

    def load_primary_dataset(self, path):
        """Load primary dataset into primary_dataset_df and initialize dataset_df copy."""
        self.primary_dataset_path = path
        self.dataset_path = path  # legacy alias
        primary_df = self._standardize_dataset_df(self._read_dataset_file(path))
        self.primary_dataset_df = primary_df
        self.dataset_df = self.primary_dataset_df.copy()
        self._reset_df_category_retest_columns()
        self.create_subset_df()
        self._merge_subset_into_df()

    def load_secondary_datasets(self, paths):
        """Load multiple secondary datasets and merge into working dataset_df."""
        if not hasattr(self, "primary_dataset_df"):
            return
        merged = [self.primary_dataset_df.copy()]
        for p in paths or []:
            try:
                df = self._standardize_dataset_df(self._read_dataset_file(p))
                merged.append(df)
            except Exception as e:
                print(f"Failed to load secondary dataset {p}: {e}")
        self.dataset_df = pd.concat(merged, ignore_index=True, sort=False)
        self._reset_df_category_retest_columns()
        self.create_subset_df()
        self._merge_subset_into_df()

    def sync_retest_to_dataset(self):
        """Update dataset_df['Retest'] from subset_df['Retest'] aligned on Compound.

        Efficient alignment via map on the Compound index. Keeps boolean dtype stable.
        """
        if not hasattr(self, "dataset_df") or not hasattr(self, "subset_df"):
            return
        try:
            # Build map Compound -> Retest from subset_df
            retest_map = (
                self.subset_df.set_index("Compound")["Retest"].astype("boolean").to_dict()
            )
            # Update dataset_df where Compound exists in map
            mask = self.dataset_df["Compound"].isin(retest_map.keys())
            if mask.any():
                self.dataset_df.loc[mask, "Retest"] = (
                    self.dataset_df.loc[mask, "Compound"].map(retest_map)
                )
                # Normalize dtype
                self.dataset_df["Retest"] = self.dataset_df["Retest"].astype("boolean")
        except Exception as e:
            print(f"Error syncing Retest to dataset: {e}")

    def sync_retest_to_primary_dataset(self):
        """Update primary_dataset_df['Retest'] only, from subset_df mapping."""
        if not hasattr(self, "primary_dataset_df") or not hasattr(self, "subset_df"):
            return
        try:
            retest_map = (
                self.subset_df.set_index("Compound")["Retest"].astype("boolean").to_dict()
            )
            mask = self.primary_dataset_df["Compound"].isin(retest_map.keys())
            if mask.any():
                self.primary_dataset_df.loc[mask, "Retest"] = (
                    self.primary_dataset_df.loc[mask, "Compound"].map(retest_map)
                )
                self.primary_dataset_df["Retest"] = self.primary_dataset_df["Retest"].astype(
                    "boolean"
                )
        except Exception as e:
            print(f"Error syncing Retest to primary dataset: {e}")

    def save_dataset(self):
        """Persist current dataset_df back to its source path (parquet or csv)."""
        if not hasattr(self, "dataset_df"):
            return
        path = getattr(self, "dataset_path", None)
        if not path:
            print("No dataset_path set; cannot save dataset.")
            return
        try:
            # Save using original format by extension
            if path.lower().endswith(".parquet"):
                self.dataset_df.to_parquet(path, index=False)
            else:
                self.dataset_df.to_csv(path, index=False)
            print(f"Saved dataset to {path}")
        except Exception as e:
            print(f"Failed to save dataset to {path}: {e}")

    def save_primary_dataset(self):
        """Persist primary_dataset_df back to its original path only."""
        if not hasattr(self, "primary_dataset_df"):
            return
        path = getattr(self, "primary_dataset_path", None)
        if not path:
            print("No primary_dataset_path set; cannot save primary dataset.")
            return
        try:
            if path.lower().endswith(".parquet"):
                self.primary_dataset_df.to_parquet(path, index=False)
            else:
                self.primary_dataset_df.to_csv(path, index=False)
            print(f"Saved primary dataset to {path}")
        except Exception as e:
            print(f"Failed to save primary dataset to {path}: {e}")

    def create_subset_df(self):
        """
        Creates subset_df from dataset_df with unique compounds and consolidated categories.
        Gets clustering information from the chemical data (self.df).
        """
        if not hasattr(self, "dataset_df"):
            raise ValueError("dataset_df not loaded. Call load_dataset first.")

        # Get clustering columns that might already exist in dataset_df
        clustering_columns = []
        for column in self.dataset_df.columns:
            if column in ["Cluster", "SuperCluster"]:
                clustering_columns.append(column)

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

        # Get clustering information from chemical data (self.df)
        if "Cluster" in self.df.columns:
            # Remove existing Cluster column if it exists
            if "Cluster" in self.subset_df.columns:
                self.subset_df = self.subset_df.drop(columns=["Cluster"])

            # Merge clustering info from chemical data into subset_df
            # cluster_info = self.df[["Compound", "Cluster"]].copy()
            # copy not needed
            cluster_info = self.df[["Compound", "Cluster"]]
            self.subset_df = self.subset_df.merge(
                cluster_info, on="Compound", how="left"
            )
        else:
            # Add empty Cluster column if no clustering data available
            self.subset_df["Cluster"] = np.nan

    def cluster_subset_df(self, coarse_thresh, radius=2, fpSize=2048):
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
        # ECFP6 is radius=3
        # fpgen = AllChem.GetMorganGenerator(radius=radius,fpSize=fpSize)
        # fpgen = AllChem.GetMorganGenerator(radius=2,fpSize=2048,atomInvariantsGenerator=AllChem.GetMorganFeatureAtomInvGen())

        # fpgen = AllChem.GetTopologicalTorsionGenerator()
        # see also https://www.researchgate.net/post/Two_similar_compounds_with_a_low_smaller_than_085_Tanimoto_coefficient2
        # for more info on how fingerprinting might not capture two near identical looking molecules.
        # different type of fingerprint can be used. see singleton clusters 3950,6353 from CNS
        fps = [fpgen.GetFingerprint(mol) for mol in mols]

        clusters = ClusterFps(fps, coarse_thresh)

        for idx, cluster in enumerate(clusters, start=1):
            self.subset_df.loc[cluster, "CoarseCluster"] = int(idx)

    def build_graph(self):
        """
        Build a graph connecting fine clusters that co-occur in at least one
        coarse cluster (defined on the smaller subset of compounds).

        Speed notes:
        - Single merge to align fine/coarse on the subset of compounds.
        - Sparse incidence B (fine x coarse), then A = B @ B.T (boolean).
        - No Python loops over clusters or compounds.
        """

        # 1) Fine clusters for all compounds (df is larger; each compound appears once)
        # BUG: potential future bug by removing the notna filtering
        # df_fine = self.df.loc[self.df["Cluster"].notna(), ["Compound", "Cluster"]]
        df_fine = self.df[["Compound", "Cluster"]]

        # 2) Coarse clusters for subset only (each compound appears once)
        df_coarse = self.subset_df[["Compound", "CoarseCluster"]]

        # 3) Restrict to compounds present in the coarse subset
        #    (this also guarantees we only relate via the coarse partition on the subset)
        df_fc = df_fine.merge(df_coarse, on="Compound", how="inner")
        if df_fc.empty:
            # No overlap => empty graph
            self.graph = nx.Graph()
            return

        # 4) Factorize to get compact integer ids (fast & memory friendly)
        fine_codes, fine_uniques = pd.factorize(df_fc["Cluster"], sort=True)
        coarse_codes, coarse_uniques = pd.factorize(df_fc["CoarseCluster"], sort=True)

        n_fine = len(fine_uniques)
        n_coarse = len(coarse_uniques)

        # 5) Build a boolean incidence matrix B (fine x coarse), 1 if any compound maps (F,C)
        # Using duplicates is fine; coo_matrix will sum them; we then binarize.
        data = np.ones(len(df_fc), dtype=np.uint8)
        B = coo_matrix(
            (data, (fine_codes, coarse_codes)), shape=(n_fine, n_coarse)
        ).tocsr()
        B.data[:] = 1  # binarize in case of multiples

        # 6) Fine–fine adjacency via sparse boolean multiplication
        # A[i,j] > 0 if fine cluster i shares any coarse cluster with j
        A = B @ B.T
        A.setdiag(0)  # remove self-loops
        A.eliminate_zeros()

        # 7) Build the graph (nodes are 0..n_fine-1 for now)
        # Note networkx includes singletons by default
        G = nx.from_scipy_sparse_array(A)  # undirected graph

        # 8) Relabel nodes to original fine cluster labels for external consistency
        mapping = {i: int(fine_uniques[i]) for i in range(n_fine)}
        G = nx.relabel_nodes(G, mapping)

        # Store results
        self.graph = G

        super_cluster_counts = []
        for subgraph in nx.connected_components(G):
            mask = self.df["Cluster"].isin(subgraph)
            compound_count = len(self.df[mask]["Compound"])
            super_cluster_counts.append([compound_count, list(subgraph)])

        # sort super_cluster_counts by compound count
        super_cluster_counts.sort(key=lambda x: x[0], reverse=True)

        # Map each cluster to its super cluster number (1-indexed) and count compounds
        self.super_clusters = [
            [super_clust_id, compound_count, sub_clusters]
            for super_clust_id, (compound_count, sub_clusters) in enumerate(
                super_cluster_counts, start=1
            )
        ]

        # Initialize super cluster column with NaN in chemical data
        self.df["SuperCluster"] = np.nan
        # Map each cluster to its super cluster number (1-indexed)
        # Note: would it be beneficial to store compound count for later saturation calculations?
        for super_clust_id, (compound_count, sub_clusters) in enumerate(
            super_cluster_counts, start=1
        ):
            self.df.loc[self.df["Cluster"].isin(sub_clusters), "SuperCluster"] = (
                super_clust_id
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

    def update_dataset_with_clustering(self):
        """
        Update dataset_df with clustering information and current Retest values.
        This should be called after clustering is complete.
        """
        if not hasattr(self, "dataset_df"):
            return

        # Get clustering information from chemical data (self.df) and Retest from subset_df
        cluster_info = self.df[["Compound", "Cluster", "SuperCluster"]].copy()
        retest_info = self.subset_df[["Compound", "Retest"]].copy()

        # Merge retest info into cluster info
        cluster_info = cluster_info.merge(retest_info, on="Compound", how="left")

        # Remove existing clustering columns from dataset_df if they exist
        columns_to_drop = []
        if "Cluster" in self.dataset_df.columns:
            columns_to_drop.append("Cluster")
        if "SuperCluster" in self.dataset_df.columns:
            columns_to_drop.append("SuperCluster")
        if "CoarseCluster" in self.dataset_df.columns:
            columns_to_drop.append("CoarseCluster")
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

    def reset_clustering_state(self):
        """Clear any attributes derived from building super clusters/graphs.

        This is used when loading a new dataset or switching libraries to ensure
        no stale state persists.
        """
        for attr in [
            "graph",
            "sub_graph",
            "node_pos",
            "super_clusters",
        ]:
            if hasattr(self, attr):
                try:
                    delattr(self, attr)
                except Exception:
                    setattr(self, attr, None)
        # Also clear SuperCluster column from df if present (it will be recomputed later)
        try:
            if hasattr(self, "df") and "SuperCluster" in self.df.columns:
                self.df = self.df.drop(columns=["SuperCluster"])
        except Exception:
            pass

    def compute_saturation_metrics(self):
        """
        Compute and attach saturation metrics to subset_df:
        - SCS (Super Cluster Saturation): for each super cluster, fraction of non-miss compounds

        Metric is merged back onto subset_df per Compound so it can be saved/exported.
        """
        if not hasattr(self, "subset_df") or not hasattr(self, "df"):
            return

        # Ensure required columns exist in chemical data
        required_cols = {"Compound", "Cluster"}
        if not required_cols.issubset(set(self.df.columns)):
            return

        # If SuperCluster hasn't been computed yet, skip gracefully
        has_super = "SuperCluster" in self.df.columns

        # Use chemical data for saturation calculations (includes all compounds)
        df_main = self.df[["Compound", "Cluster"]].copy()
        if has_super:
            df_main["SuperCluster"] = self.df["SuperCluster"]

        # Add category information from subset_df
        category_info = self.subset_df[["Compound", "Category"]].copy()
        df_main = df_main.merge(category_info, on="Compound", how="left")
        # Fill missing categories as "Miss" (compounds not in dataset)
        df_main["Category"] = df_main["Category"].fillna("Miss")

        # Helper to compute saturation: non-miss / total
        def _saturation(group: pd.DataFrame) -> float:
            total = len(group)
            if total == 0:
                return np.nan
            non_miss = (group["Category"] != "Miss").sum()
            return float(non_miss) / float(total)

        # Super cluster saturation (SCS)
        if has_super:
            scs_series = (
                df_main.dropna(subset=["SuperCluster"])
                .groupby("SuperCluster")
                .apply(_saturation)
            )
            scs_map = scs_series.to_dict()
        else:
            scs_map = {}

        # Map metrics back per compound using df_main (compound -> cluster ids)
        df_metrics = df_main[["Compound", "Cluster"]].copy()

        if has_super:
            df_metrics["SCS"] = df_main["SuperCluster"].map(scs_map)
        else:
            df_metrics["SCS"] = np.nan

        df_metrics = df_metrics.drop(columns=["Cluster"])

        # Merge into subset_df per compound
        self.subset_df = self.subset_df.merge(df_metrics, on="Compound", how="left")

    def build_subgraph(self, member_cluster):
        """
        Builds subgraph of related clusters
        """
        nodes = nx.node_connected_component(self.graph, member_cluster)
        self.sub_graph = self.graph.subgraph(nodes)
        self.node_pos = nx.spring_layout(self.sub_graph, iterations=50)

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

        subset = self.df[self.df["Compound"].isin(compound_ids)].copy()

        # Respect input order: sort subset according to the ordered compound_ids list
        try:
            order = pd.Categorical(
                subset["Compound"], categories=compound_ids, ordered=True
            )
            subset = (
                subset.assign(__order=order)
                .sort_values("__order")
                .drop(columns=["__order"])
            )
        except Exception:
            # If anything goes wrong, fall back to original subset order
            pass

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
        # redefine compound_ids below to handle input of bad ID, preserving order from get_compounds
        compound_ids = chem_info["Compound"].values.tolist()
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
                    (str(cid) + "\n\n" + str(cat)) if cat != "Miss" else str(cid)
                    for cid, cat in zip(ids, categories)
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
