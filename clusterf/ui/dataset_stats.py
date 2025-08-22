from __future__ import annotations
from typing import TYPE_CHECKING, Dict, Any

import pandas as pd
import panel as pn
import param

if TYPE_CHECKING:
    from clusterf.app import ClusterFApp


class DatasetStatsCard(param.Parameterized):
    app: "ClusterFApp"

    def __init__(self, app: "ClusterFApp", **params) -> None:
        super().__init__(**params)
        self.app = app

        self._primary_pane = pn.pane.Markdown("", sizing_mode="stretch_width")
        self._full_pane = pn.pane.Markdown("", sizing_mode="stretch_width")

        # Build card
        self.view = pn.Card(
            pn.pane.Markdown("**Primary Dataset**", margin=(0, 0, 2, 0)),
            self._primary_pane,
            pn.layout.Divider(),
            pn.pane.Markdown("**Full Dataset (with secondaries)**", margin=(4, 0, 2, 0)),
            self._full_pane,
            title="Dataset Stats",
            width=220,
            collapsed=False,
        )

        # Refresh initially if data exists
        self.refresh()

        # Watch for dataset loaded toggles
        if hasattr(self.app, "param") and hasattr(self.app.param, "dataset_loaded"):
            self.app.param.watch(self._on_dataset_loaded, "dataset_loaded")

    # Public API to refresh from other components
    def refresh(self) -> None:
        primary_stats = self._compute_primary_stats()
        full_stats = self._compute_full_stats()

        self._primary_pane.object = self._render_stats_md(primary_stats)
        self._full_pane.object = self._render_stats_md(full_stats)

    # Watcher callback
    def _on_dataset_loaded(self, _):
        self.refresh()

    # ---- Internals ----
    def _stats_from_df(self, df: pd.DataFrame, dedupe_by_compound: bool = True) -> Dict[str, Any]:
        if df is None or df.empty or "Category" not in df.columns:
            return {"total_hits": 0, "interfering": 0, "retest": 0, "by_category": {}}

        value_counts = df[df.Category!="Miss"].drop_duplicates(subset="Compound").Category.value_counts()
        

        # Interfering total
        try:
            interfering_total = value_counts['Interfering']
        except KeyError:
            interfering_total = 0


        value_counts.drop(index=['Interfering'], inplace=True)
        non_interfering_hit_counts = value_counts

        total_hits = int(value_counts.sum())

        # Retest total (unique compounds with any True)
        
        retest_total = int(df.groupby("Compound")["Retest"].any().sum())
        
        #     retest_total = int(bool(tmp["Retest"].any()))

        return {
            "total_hits": total_hits,
            "interfering": interfering_total,
            "retest": retest_total,
            "by_category": {k: int(v) for k, v in non_interfering_hit_counts.to_dict().items()},
        }

    def _compute_primary_stats(self) -> Dict[str, Any]:
        lib = getattr(self.app, "library", None)
        if lib is None or not hasattr(lib, "primary_dataset_df"):
            return {"total_hits": 0, "interfering": 0, "retest": 0, "by_category": {}}
        return self._stats_from_df(lib.primary_dataset_df, dedupe_by_compound=True)

    def _compute_full_stats(self) -> Dict[str, Any]:
        lib = getattr(self.app, "library", None)
        if lib is None:
            return {"total_hits": 0, "interfering": 0, "retest": 0, "by_category": {}}

        # Prefer subset_df (already unique per compound); else use working dataset and dedupe
        if hasattr(lib, "subset_df") and isinstance(lib.subset_df, pd.DataFrame) and not lib.subset_df.empty:
            return self._stats_from_df(lib.subset_df, dedupe_by_compound=False)
        elif hasattr(lib, "dataset_df"):
            return self._stats_from_df(lib.dataset_df, dedupe_by_compound=True)
        else:
            return {"total_hits": 0, "interfering": 0, "retest": 0, "by_category": {}}

    def _render_stats_md(self, stats: Dict[str, Any]) -> str:
        if not stats:
            return "_No data_"
        lines = [
            f"- **Hits (total)**: {stats.get('total_hits', 0)}",
            f"- **Interfering**: {stats.get('interfering', 0)}",
            f"- **Retest selected**: {stats.get('retest', 0)}",
        ]
        by_cat = stats.get("by_category", {}) or {}
        if by_cat:
            lines.append("- **By category**:")
            for cat, cnt in by_cat.items():
                lines.append(f"  - {cat}: {cnt}")
        return "\n".join(lines)


