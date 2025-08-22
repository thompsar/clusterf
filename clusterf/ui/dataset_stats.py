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
    def _aggregate_categories(self, df: pd.DataFrame) -> pd.DataFrame:
        if df is None or df.empty or "Compound" not in df.columns or "Category" not in df.columns:
            return pd.DataFrame(columns=["Compound", "Category", "Retest"]).astype({"Retest": "boolean"})

        # Normalize types
        temp = df.copy()
        if "Compound" in temp.columns:
            temp["Compound"] = temp["Compound"].astype(str)
        if "Retest" not in temp.columns:
            temp["Retest"] = False
        try:
            temp["Retest"] = temp["Retest"].astype("boolean").fillna(False)
        except Exception:
            pass

        # Exclude Miss rows for aggregation
        non_miss = temp[temp["Category"] != "Miss"].copy()
        if non_miss.empty:
            out = (
                temp[["Compound"]]
                .drop_duplicates()
                .assign(Category="Miss", Retest=False)
            )
            return out.astype({"Retest": "boolean"})

        def _aggregate_category_label(group: pd.DataFrame) -> str:
            cats = group["Category"].astype(str)
            # Interfering overrides everything
            if cats.str.contains("Interfering", case=False, na=False).any():
                return "Interfering"

            plus_mask = cats == "Hit(+)"
            minus_mask = cats == "Hit(-)"

            # Build construct-based labels if Construct is present
            if "Construct" in group.columns:
                constructs_plus = (
                    group.loc[plus_mask, "Construct"].astype(str).unique().tolist()
                )
                constructs_minus = (
                    group.loc[minus_mask, "Construct"].astype(str).unique().tolist()
                )
                constructs_plus = sorted([c for c in constructs_plus if pd.notna(c)])
                constructs_minus = sorted([c for c in constructs_minus if pd.notna(c)])
                if constructs_plus and constructs_minus:
                    constructs = sorted(set(constructs_plus) | set(constructs_minus))
                    return f"{'/'.join(constructs)} Hit(+/-)"
                if constructs_plus:
                    return f"{'/'.join(constructs_plus)} Hit(+)"
                if constructs_minus:
                    return f"{'/'.join(constructs_minus)} Hit(-)"

            # Fallback to first non-miss category
            non_miss_cats = cats[cats != "Miss"].dropna().astype(str)
            return non_miss_cats.iloc[0] if not non_miss_cats.empty else "Miss"

        agg_cat = non_miss.groupby("Compound").apply(_aggregate_category_label)
        out = agg_cat.reset_index().rename(columns={0: "Category"})

        # Retest aggregated per compound (any True)
        try:
            retest_map = (
                temp.groupby("Compound")["Retest"].apply(lambda s: bool(pd.Series(s).astype(bool).any())).to_dict()
            )
        except Exception:
            retest_map = {}
        out["Retest"] = out["Compound"].map(retest_map).fillna(False)
        try:
            out["Retest"] = out["Retest"].astype("boolean")
        except Exception:
            pass
        return out

    def _compute_primary_stats(self) -> Dict[str, Any]:
        lib = getattr(self.app, "library", None)
        if lib is None or not hasattr(lib, "primary_dataset_df"):
            return {"total_hits": 0, "interfering": 0, "retest": 0, "by_category": {}}

        agg = self._aggregate_categories(lib.primary_dataset_df)
        return self._stats_from_agg(agg)

    def _compute_full_stats(self) -> Dict[str, Any]:
        lib = getattr(self.app, "library", None)
        if lib is None:
            return {"total_hits": 0, "interfering": 0, "retest": 0, "by_category": {}}

        # Prefer subset_df if present, else aggregate from dataset_df
        if hasattr(lib, "subset_df") and isinstance(lib.subset_df, pd.DataFrame) and not lib.subset_df.empty:
            agg = lib.subset_df[["Compound", "Category", "Retest"]].copy()
            if "Retest" not in agg.columns:
                agg["Retest"] = False
            try:
                agg["Retest"] = agg["Retest"].astype("boolean").fillna(False)
            except Exception:
                pass
        elif hasattr(lib, "dataset_df"):
            agg = self._aggregate_categories(lib.dataset_df)
        else:
            agg = pd.DataFrame(columns=["Compound", "Category", "Retest"])  # empty

        return self._stats_from_agg(agg)

    def _stats_from_agg(self, agg: pd.DataFrame) -> Dict[str, Any]:
        if agg is None or agg.empty:
            return {"total_hits": 0, "interfering": 0, "retest": 0, "by_category": {}}

        cats = agg["Category"].astype(str)
        # Total hits: any category indicating a hit (exclude Miss and Interfering)
        is_hit = cats.str.contains("Hit", na=False)
        is_interfering = cats.str.contains("Interfering", case=False, na=False)
        total_hits = int(((is_hit) & (~is_interfering)).sum())

        # Interfering total
        total_interfering = int(is_interfering.sum())

        # Retest total
        total_retest = int(pd.Series(agg.get("Retest", False)).astype(bool).sum())

        # Breakdown by hit category (exclude Interfering and Miss)
        mask = (~is_interfering) & (cats != "Miss")
        by_cat = (
            cats[mask]
            .value_counts()
            .sort_index()
            .to_dict()
        )

        return {
            "total_hits": total_hits,
            "interfering": total_interfering,
            "retest": total_retest,
            "by_category": by_cat,
        }

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


