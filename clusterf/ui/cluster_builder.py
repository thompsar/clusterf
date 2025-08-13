from __future__ import annotations
from typing import TYPE_CHECKING
import os
import panel as pn
import param

from clusterf.core.chemistry import ChemLibrary

if TYPE_CHECKING:
    from clusterf.app import ClusterFApp


class SuperClusterBuilder(param.Parameterized):
    app: ClusterFApp
    library_select = param.Selector(default=None, doc="Select a compound library")
    dataset_select = param.Selector(default=None, doc="Select a dataset")
    method = param.Selector(objects=["RDKit"], default="RDKit", doc="Clustering method")
    fine_threshold = param.Selector(objects=[0.2], default=0.2)
    coarse_threshold = param.Selector(
        objects=[0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55], default=0.4
    )
    cluster_button = param.Action(lambda self: self.param.trigger("cluster_button"))
    clusters_built = param.Boolean(default=False, doc="Whether super clusters have been built")

    def __init__(self, app: ClusterFApp, **params) -> None:
        super().__init__(**params)
        self.app = app

        # Get available libraries (strip .parquet extension)
        compound_libraries = []
        if os.path.exists(app.LIBRARIES_DIR):
            for f in os.listdir(app.LIBRARIES_DIR):
                if f.endswith(".parquet") and not f.endswith("_clusters.parquet"):
                    library_name = f.replace(".parquet", "")
                    compound_libraries.append(library_name)
        
        # --- Create widgets first---
        self.library_select_widget = pn.widgets.Select.from_param(
            self.param.library_select, name="Library", width=200
        )
        self.dataset_select_widget = pn.widgets.Select.from_param(
            self.param.dataset_select, name="Dataset", width=200
        )
        
        # Create side-by-side method and threshold selectors
        self.method_widget = pn.widgets.Select.from_param(
            self.param.method, name="Method", width=95
        )
        self.fine_threshold_widget = pn.widgets.Select.from_param(
            self.param.fine_threshold, name="Threshold", width=95
        )
        
        # Combine method and threshold widgets side by side
        self.clustering_params_row = pn.Row(
            self.method_widget,
            self.fine_threshold_widget,
            width=200
        )
        
        self.coarse_threshold_widget = pn.widgets.Select.from_param(
            self.param.coarse_threshold, name="Coarse threshold", width=200
        )
        self.cluster_button_widget = pn.widgets.Button.from_param(
            self.param.cluster_button,
            name="Build Super Clusters",
            button_type="primary",
            width=200,
        )

        # Create controls first
        self.controls = pn.Card(
            self.library_select_widget,
            self.dataset_select_widget,
            self.clustering_params_row,
            self.coarse_threshold_widget,
            self.cluster_button_widget,
            title="Cluster Library",
            width=220,
            collapsed=False,
        )

        # Now set the library and dataset selections (this will trigger the watchers)
        compound_libraries = sorted(compound_libraries)
        self.param.library_select.objects = compound_libraries
        if compound_libraries:
            self.library_select = compound_libraries[0]

        datasets = sorted(
            f for f in os.listdir(app.DATASETS_DIR) if f.endswith((".csv", ".parquet"))
        )
        self.param.dataset_select.objects = datasets
        if datasets:
            self.dataset_select = datasets[0]

    @param.depends("library_select", watch=True)
    def _load_library(self):
        if self.library_select:
            self.app.library = ChemLibrary(
                self.library_select,
                fine_threshold=self.fine_threshold,
                method=self.method,
                libraries_dir=self.app.LIBRARIES_DIR
            )
            
            # Update clustering parameter options based on available clustering data
            self._update_clustering_options()
            
            # Reset clusters when library changes
            self.clusters_built = False
            
            # Expand the cluster builder card when library changes
            self.controls.collapsed = False
            # print(self.app.library.df.head())

    def _update_clustering_options(self):
        """Update clustering parameter options based on available clustering data"""
        if hasattr(self.app, "library") and self.app.library:
            available_params = self.app.library.get_available_clustering_parameters()
            
            # Update method options
            if available_params["methods"]:
                self.param.method.objects = available_params["methods"]
                # Set default to first available method if current is not available
                if self.method not in available_params["methods"]:
                    self.method = available_params["methods"][0]
            
            # Update threshold options
            if available_params["thresholds"]:
                self.param.fine_threshold.objects = available_params["thresholds"]
                # Set default to first available threshold if current is not available
                if self.fine_threshold not in available_params["thresholds"]:
                    self.fine_threshold = available_params["thresholds"][0]

    @param.depends("dataset_select", watch=True)
    def _load_dataset(self):
        if hasattr(self.app, "library") and self.app.library and self.dataset_select:
            self.app.library.load_dataset(
                os.path.join(self.app.DATASETS_DIR, self.dataset_select)
            )
            # Reset clusters when dataset changes
            self.clusters_built = False
            
            # Expand the cluster builder card when dataset changes
            self.controls.collapsed = False
            # print(self.app.library.dataset_df.head())
    
    @param.depends("method", "fine_threshold", "coarse_threshold", watch=True)
    def _reset_clusters(self):
        """Reset clusters when clustering parameters change."""
        if hasattr(self.app, "library") and self.app.library:
            # Update clustering parameters and reload clustering info
            self.app.library.set_clustering_parameters(
                fine_threshold=self.fine_threshold,
                method=self.method
            )
        self.clusters_built = False
        
        # Expand the cluster builder card when parameters change
        self.controls.collapsed = False

    @param.depends("cluster_button", watch=True)
    def _build_super_clusters(self):
        if getattr(self.app, "library", None):
            try:
                self.app.library.cluster_subset_df(self.coarse_threshold)
                self.app.library.build_graph(self.coarse_threshold)
                
                # Collapse the cluster builder card after successful clustering
                self.controls.collapsed = True
                # Signal that clusters have been built successfully
                self.clusters_built = True
                print(f"Super clusters built successfully! "
                      f"Found {len(getattr(self.app.library, 'super_clusters', []))} super clusters.")
                
                
                
            except Exception as e:
                print(f"Error building super clusters: {e}")
                self.clusters_built = False
        else:
            print("No library loaded. Please select a library first.")
            self.clusters_built = False
