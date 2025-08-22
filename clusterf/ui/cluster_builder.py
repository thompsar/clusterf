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
    # Primary/Secondary dataset widgets
    primary_dataset_select = param.Selector(default=None, doc="Primary dataset")
    secondary_datasets_select = param.ListSelector(default=[], doc="Secondary datasets")
    method = param.Selector(objects=["RDKit"], default="RDKit", doc="Clustering method")
    fine_threshold = param.Selector(objects=[0.2], default=0.2)
    coarse_threshold = param.Selector(
        objects=[0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55], default=0.4
    )
    cluster_button = param.Action(lambda self: self.param.trigger("cluster_button"))
    clusters_built = param.Boolean(
        default=False, doc="Whether super clusters have been built"
    )

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
        self.primary_dataset_widget = pn.widgets.Select.from_param(
            self.param.primary_dataset_select, name="Primary Dataset", width=200
        )
        self.secondary_datasets_widget = pn.widgets.MultiSelect.from_param(
            self.param.secondary_datasets_select, name="Secondary Datasets", width=200, size=6
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
            self.method_widget, self.fine_threshold_widget, width=200
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
            self.primary_dataset_widget,
            self.secondary_datasets_widget,
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
        self.param.primary_dataset_select.objects = datasets
        if datasets:
            self.primary_dataset_select = datasets[0]
            # Secondary initially excludes primary
            self._update_secondary_options()

    @param.depends("library_select", watch=True)
    def _load_library(self):
        if self.library_select:
            # Reset the app state before loading a new library
            if hasattr(self.app, "reset_app_state"):
                self.app.reset_app_state()
            self.app.library = ChemLibrary(
                self.library_select,
                fine_threshold=self.fine_threshold,
                method=self.method,
                libraries_dir=self.app.LIBRARIES_DIR,
            )

            # Update clustering parameter options based on available clustering data
            self._update_clustering_options()

            # Reset clusters when library changes
            self.clusters_built = False

            # Expand the cluster builder card when library changes
            self.controls.collapsed = False
            
            # Trigger dataset change/reset sequence for UI widgets that depend on categories
            if hasattr(self.app, "_on_dataset_loaded"):
                self.app._on_dataset_loaded()
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

    def _update_secondary_options(self):
        """Refresh secondary options to exclude the currently selected primary."""
        all_datasets = list(self.param.primary_dataset_select.objects or [])
        primary = self.primary_dataset_select
        secondary_options = [d for d in all_datasets if d != primary]
        self.param.secondary_datasets_select.objects = secondary_options
        # Drop any selected secondaries that now conflict with primary
        self.secondary_datasets_select = [d for d in self.secondary_datasets_select if d in secondary_options]

    @param.depends("primary_dataset_select", watch=True)
    def _on_primary_dataset_changed(self):
        # Adjust secondary options when primary changes
        self._update_secondary_options()
        if hasattr(self.app, "library") and self.app.library and self.primary_dataset_select:
            # Reset state before loading new primary
            if hasattr(self.app, "reset_app_state"):
                self.app.reset_app_state()
            self.app.library.load_primary_dataset(
                os.path.join(self.app.DATASETS_DIR, self.primary_dataset_select)
            )
            self.clusters_built = False
            self.controls.collapsed = False
            if hasattr(self.app, "_on_dataset_loaded"):
                self.app._on_dataset_loaded()

    @param.depends("secondary_datasets_select", watch=True)
    def _on_secondary_datasets_changed(self):
        if hasattr(self.app, "library") and self.app.library:
            # Load secondaries into working dataset
            paths = [os.path.join(self.app.DATASETS_DIR, f) for f in (self.secondary_datasets_select or [])]
            self.app.library.load_secondary_datasets(paths)
            self.clusters_built = False
            self.controls.collapsed = False
            if hasattr(self.app, "_on_dataset_loaded"):
                self.app._on_dataset_loaded()

    @param.depends("method", "fine_threshold", "coarse_threshold", watch=True)
    def _reset_clusters(self):
        """Reset clusters when clustering parameters change."""
        if hasattr(self.app, "library") and self.app.library:
            # Update clustering parameters and reload clustering info
            self.app.library.set_clustering_parameters(
                fine_threshold=self.fine_threshold, method=self.method
            )
        self.clusters_built = False
        # Clear any computed clustering state on library
        if hasattr(self.app, "library") and hasattr(self.app.library, "reset_clustering_state"):
            self.app.library.reset_clustering_state()

        # Expand the cluster builder card when parameters change
        self.controls.collapsed = False

    @param.depends("cluster_button", watch=True)
    def _build_super_clusters(self):
        if getattr(self.app, "library", None):
            try:
                # Clear previous clustering state before building anew
                if hasattr(self.app.library, "reset_clustering_state"):
                    self.app.library.reset_clustering_state()
                self.app.library.cluster_subset_df(self.coarse_threshold)
                self.app.library.build_graph()

                # Collapse the cluster builder card after successful clustering
                self.controls.collapsed = True
                # Signal that clusters have been built successfully
                self.clusters_built = True

            except Exception as e:
                print(f"Error building super clusters: {e}")
                self.clusters_built = False
        else:
            print("No library loaded. Please select a library first.")
            self.clusters_built = False
