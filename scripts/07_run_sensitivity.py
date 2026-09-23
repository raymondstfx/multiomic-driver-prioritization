"""Run Phase 7 robustness analyses without changing the primary ranking."""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _common import run_stage
from scipy import sparse
from scipy.stats import spearmanr

from multiomic_driver.evaluation.sensitivity import (
    combine_sensitivity_scores,
    common_guide_sets,
    component_depth_diagnostics,
    effect_dimension_contributions,
    guide_composition,
    latent_centroid_effects,
)
from multiomic_driver.utils.artifacts import figure_dir, table_dir
from multiomic_driver.utils.io import require_path

LOGGER = logging.getLogger("07_run_sensitivity")


def _raw_atac_depth(path: Path, n_obs: int) -> tuple[np.ndarray, np.ndarray]:
    source = ad.read_h5ad(path, backed="r")
    try:
        if source.n_obs != n_obs:
            raise ValueError(
                "Raw and processed ATAC objects have different cell counts"
            )
        total = np.zeros(n_obs, dtype=float)
        detected = np.zeros(n_obs, dtype=float)
        for start in range(0, n_obs, 2000):
            stop = min(start + 2000, n_obs)
            block = source.X[start:stop]
            total[start:stop] = np.asarray(block.sum(axis=1)).ravel()
            detected[start:stop] = (
                block.getnnz(axis=1)
                if sparse.issparse(block)
                else np.count_nonzero(block, axis=1)
            )
        return total, detected
    finally:
        source.file.close()


def _scores(summary: pd.DataFrame) -> pd.Series:
    return summary.set_index("candidate")["score"]


def _rank_correlation(primary: pd.DataFrame, alternative: pd.DataFrame) -> float:
    merged = primary[["candidate", "multiomic_rank"]].merge(
        alternative[["candidate", "multiomic_rank"]], on="candidate"
    )
    return float(
        spearmanr(merged["multiomic_rank_x"], merged["multiomic_rank_y"]).statistic
    )


def main(config: dict) -> None:
    processed = Path(config["paths"]["processed"])
    interim = Path(config["paths"]["interim"])
    phase04 = table_dir(config, "phase04")
    phase05 = table_dir(config, "phase05")
    phase035 = table_dir(config, "phase035")
    tables = table_dir(config, "phase07")
    figures = figure_dir(config, "phase07")
    settings = config["sensitivity"]
    primary = pd.read_csv(require_path(phase05 / "candidate_ranking.csv"))
    eligibility = pd.read_csv(
        require_path(phase035 / "candidate_phase4_eligibility.csv")
    )
    candidates = primary["candidate"].astype(str).tolist()

    rna = ad.read_h5ad(require_path(processed / "rna_processed.h5ad"), backed="r")
    atac = ad.read_h5ad(require_path(processed / "atac_processed.h5ad"), backed="r")
    try:
        if not rna.obs_names.equals(atac.obs_names):
            raise ValueError("RNA and ATAC cells are not aligned")
        metadata = rna.obs.copy()
        rna_matrix = np.asarray(rna.obsm["X_pca"])
        atac_matrix = np.asarray(atac.obsm["X_lsi"])

        atac_total, atac_detected = _raw_atac_depth(
            require_path(interim / "atac_merged.h5ad"), atac.n_obs
        )
        diagnostics = pd.concat(
            [
                component_depth_diagnostics(
                    rna_matrix,
                    metadata["total_counts"].to_numpy(),
                    metadata["n_genes_by_counts"].to_numpy(),
                    prefix="PC",
                    explained_variance_ratio=np.asarray(
                        rna.uns["pca"]["variance_ratio"]
                    ),
                ).assign(modality="RNA"),
                component_depth_diagnostics(
                    atac_matrix,
                    atac_total,
                    atac_detected,
                    prefix="LSI",
                    singular_values=np.asarray(atac.uns["lsi"]["singular_values"]),
                    explained_variance_ratio=np.asarray(
                        atac.uns["lsi"]["explained_variance_ratio"]
                    ),
                ).assign(modality="ATAC"),
            ],
            ignore_index=True,
        )
        diagnostics.to_csv(tables / "latent_component_diagnostics.csv", index=False)

        contributions = pd.concat(
            [
                effect_dimension_contributions(
                    pd.read_csv(phase04 / "rna_effect_summary.csv"),
                    prefix="PC",
                    modality="RNA",
                ),
                effect_dimension_contributions(
                    pd.read_csv(phase04 / "atac_effect_summary.csv"),
                    prefix="LSI",
                    modality="ATAC",
                ),
            ],
            ignore_index=True,
        )
        contributions.to_csv(tables / "effect_dimension_contributions.csv", index=False)

        minimum = int(config["phase4"]["min_cells_per_group"])
        methods = ["mean", "winsorised_mean", "trimmed_mean", "median"]
        method_rankings = []
        method_effects = []
        method_summary = []
        for method in methods:
            _, rna_summary = latent_centroid_effects(
                rna_matrix,
                metadata,
                candidates,
                min_cells_per_group=minimum,
                method=method,
                winsor_fraction=float(settings["winsor_fraction"]),
                trim_fraction=float(settings["trim_fraction"]),
            )
            _, atac_summary = latent_centroid_effects(
                atac_matrix,
                metadata,
                candidates,
                min_cells_per_group=minimum,
                method=method,
                winsor_fraction=float(settings["winsor_fraction"]),
                trim_fraction=float(settings["trim_fraction"]),
            )
            ranking = combine_sensitivity_scores(
                _scores(rna_summary), _scores(atac_summary)
            )
            if method == "mean":
                for column, observed in (
                    ("rna_score", _scores(rna_summary)),
                    ("atac_score", _scores(atac_summary)),
                ):
                    expected = primary.set_index("candidate")[column].reindex(
                        observed.index
                    )
                    if not np.allclose(
                        observed.to_numpy(),
                        expected.to_numpy(),
                        rtol=1e-5,
                        atol=1e-8,
                    ):
                        raise ValueError(
                            f"Sensitivity mean recomputation does not match {column}"
                        )
            ranking.insert(0, "analysis", method)
            ranking.insert(1, "is_primary_estimator", method == "mean")
            ranking["primary_rank_spearman"] = _rank_correlation(primary, ranking)
            method_rankings.append(ranking)
            method_effects.extend(
                [
                    rna_summary.assign(modality="RNA"),
                    atac_summary.assign(modality="ATAC"),
                ]
            )
            method_summary.append(
                {
                    "analysis": method,
                    "is_primary_estimator": method == "mean",
                    "candidate_count": len(ranking),
                    "rank_spearman_vs_primary": _rank_correlation(primary, ranking),
                    "hic2_rank": int(
                        ranking.loc[
                            ranking["candidate"].eq("HIC2"), "multiomic_rank"
                        ].iloc[0]
                    ),
                }
            )
        pd.concat(method_rankings, ignore_index=True).to_csv(
            tables / "robust_centroid_rankings.csv", index=False
        )
        pd.concat(method_effects, ignore_index=True).to_csv(
            tables / "robust_centroid_effects.csv", index=False
        )

        primary_effects = {
            "RNA": pd.read_csv(phase04 / "rna_effect_summary.csv"),
            "ATAC": pd.read_csv(phase04 / "atac_effect_summary.csv"),
        }
        loo_rankings = []
        for modality, frame in primary_effects.items():
            prefix = "PC" if modality == "RNA" else "LSI"
            dimensions = [column for column in frame if column.startswith(prefix)]
            other_modality = "ATAC" if modality == "RNA" else "RNA"
            other_scores = primary_effects[other_modality].set_index("target_gene")[
                "mean_effect_norm"
            ]
            for removed in dimensions:
                retained = [column for column in dimensions if column != removed]
                scores = pd.Series(
                    np.linalg.norm(frame[retained].to_numpy(), axis=1),
                    index=frame["target_gene"],
                )
                rna_scores = scores if modality == "RNA" else other_scores
                atac_scores = scores if modality == "ATAC" else other_scores
                ranking = combine_sensitivity_scores(rna_scores, atac_scores)
                ranking.insert(0, "removed_dimension", removed)
                ranking.insert(0, "modality", modality)
                ranking["rank_spearman_vs_primary"] = _rank_correlation(
                    primary, ranking
                )
                loo_rankings.append(ranking)
        pd.concat(loo_rankings, ignore_index=True).to_csv(
            tables / "leave_one_dimension_out_rankings.csv", index=False
        )

        ntc = (
            metadata["guide_assignment_status"].astype(str).eq("singlet")
            & metadata["is_non_targeting"].astype(str).str.lower().eq("true")
        ).to_numpy()
        scaled = {}
        for modality, matrix, frame, prefix in (
            ("RNA", rna_matrix, primary_effects["RNA"], "PC"),
            ("ATAC", atac_matrix, primary_effects["ATAC"], "LSI"),
        ):
            sd = matrix[ntc].std(axis=0, ddof=0)
            sd[sd == 0] = 1.0
            dimensions = [column for column in frame if column.startswith(prefix)]
            scaled[modality] = pd.Series(
                np.linalg.norm(frame[dimensions].to_numpy() / sd, axis=1),
                index=frame["target_gene"],
            )
        ntc_ranking = combine_sensitivity_scores(scaled["RNA"], scaled["ATAC"])
        ntc_ranking["rank_spearman_vs_primary"] = _rank_correlation(
            primary, ntc_ranking
        )
        ntc_ranking.to_csv(tables / "ntc_sd_scaled_ranking.csv", index=False)

        composition = guide_composition(metadata)
        composition.to_csv(tables / "guide_composition.csv", index=False)
        composition.loc[
            composition["target_gene"].astype(str).str.startswith("NTC")
        ].to_csv(tables / "ntc_guide_composition.csv", index=False)
        common = common_guide_sets(composition, min_cells_per_group=minimum)
        common_candidates = [
            candidate for candidate in candidates if common.get(candidate)
        ]
        if len(common_candidates) >= 2:
            _, common_rna = latent_centroid_effects(
                rna_matrix,
                metadata,
                common_candidates,
                min_cells_per_group=minimum,
                allowed_guides=common,
            )
            _, common_atac = latent_centroid_effects(
                atac_matrix,
                metadata,
                common_candidates,
                min_cells_per_group=minimum,
                allowed_guides=common,
            )
            common_ranking = combine_sensitivity_scores(
                _scores(common_rna), _scores(common_atac)
            )
            common_ranking.to_csv(tables / "common_guide_ranking.csv", index=False)

        ntc_targets = composition.loc[
            composition["target_gene"].astype(str).str.startswith("NTC")
        ]
        ntc_group_counts = (
            ntc_targets.groupby(
                ["target_gene", "condition", "replicate"], observed=True
            )["cell_count"]
            .sum()
            .unstack(["condition", "replicate"], fill_value=0)
        )
        expected_ntc_groups = pd.MultiIndex.from_product(
            [["DMSO", "Dasatinib"], [1, 2]]
        )
        ntc_group_counts = ntc_group_counts.reindex(
            columns=expected_ntc_groups, fill_value=0
        )
        balanced_ntcs = set(
            ntc_group_counts.index[
                ntc_group_counts.min(axis=1).ge(minimum)
            ].astype(str)
        )
        if balanced_ntcs:
            _, balanced_rna = latent_centroid_effects(
                rna_matrix,
                metadata,
                candidates,
                min_cells_per_group=minimum,
                allowed_ntc_targets=balanced_ntcs,
            )
            _, balanced_atac = latent_centroid_effects(
                atac_matrix,
                metadata,
                candidates,
                min_cells_per_group=minimum,
                allowed_ntc_targets=balanced_ntcs,
            )
            balanced_ranking = combine_sensitivity_scores(
                _scores(balanced_rna), _scores(balanced_atac)
            )
            balanced_ranking.insert(
                0, "included_ntc_targets", ";".join(sorted(balanced_ntcs))
            )
            balanced_ranking["rank_spearman_vs_primary"] = _rank_correlation(
                primary, balanced_ranking
            )
            balanced_ranking.to_csv(
                tables / "balanced_ntc_ranking.csv", index=False
            )

        guide_flags = (
            composition.groupby("target_gene", observed=True)
            .agg(
                exact_guide_classes=("guide", "nunique"),
                all_classes_present_in_four_groups=(
                    "present_in_all_four_groups",
                    "all",
                ),
                minimum_class_group_count=("cell_count", "min"),
                maximum_class_fraction=("within_target_group_fraction", "max"),
            )
            .reset_index()
        )
        guide_flags["composition_sensitive_flag"] = (
            ~guide_flags["all_classes_present_in_four_groups"]
            | guide_flags["maximum_class_fraction"].gt(0.8)
        )
        guide_flags.to_csv(tables / "guide_composition_flags.csv", index=False)

        replicate_rankings = []
        for replicate in (1, 2):
            modality_scores = {}
            for modality in ("rna", "atac"):
                effects = pd.read_csv(phase04 / f"{modality}_replicate_effects.csv")
                dims = [
                    c
                    for c in effects
                    if c.startswith("PC" if modality == "rna" else "LSI")
                ]
                selected = effects.loc[effects["replicate"].eq(replicate)]
                modality_scores[modality] = pd.Series(
                    np.linalg.norm(selected[dims].to_numpy(), axis=1),
                    index=selected["target_gene"],
                )
            ranking = combine_sensitivity_scores(
                modality_scores["rna"], modality_scores["atac"]
            )
            ranking.insert(0, "replicate", replicate)
            ranking["rank_spearman_vs_primary"] = _rank_correlation(primary, ranking)
            replicate_rankings.append(ranking)
        pd.concat(replicate_rankings, ignore_index=True).to_csv(
            tables / "replicate_specific_rankings.csv", index=False
        )

        threshold_rankings = []
        count_columns = [
            "dmso_r1_cells",
            "dmso_r2_cells",
            "dasatinib_r1_cells",
            "dasatinib_r2_cells",
        ]
        for threshold in settings["cell_thresholds"]:
            eligible = (
                eligibility.loc[
                    eligibility[count_columns].notna().all(axis=1)
                    & eligibility[count_columns].min(axis=1).ge(int(threshold)),
                    "target_gene",
                ]
                .astype(str)
                .tolist()
            )
            if len(eligible) < 2:
                continue
            _, rs = latent_centroid_effects(
                rna_matrix, metadata, eligible, min_cells_per_group=int(threshold)
            )
            _, ats = latent_centroid_effects(
                atac_matrix, metadata, eligible, min_cells_per_group=int(threshold)
            )
            ranking = combine_sensitivity_scores(_scores(rs), _scores(ats))
            ranking.insert(0, "min_cells_per_group", int(threshold))
            threshold_rankings.append(ranking)
        pd.concat(threshold_rankings, ignore_index=True).to_csv(
            tables / "cell_threshold_rankings.csv", index=False
        )

        weight_rankings = []
        rna_primary = primary.set_index("candidate")["rna_score"]
        atac_primary = primary.set_index("candidate")["atac_score"]
        for weight in settings["modality_rna_weights"]:
            ranking = combine_sensitivity_scores(
                rna_primary, atac_primary, rna_weight=float(weight)
            )
            weight_rankings.append(ranking)
        pd.concat(weight_rankings, ignore_index=True).to_csv(
            tables / "modality_weight_rankings.csv", index=False
        )

        summary = pd.DataFrame(method_summary)
        summary.to_csv(tables / "sensitivity_summary.csv", index=False)
        pivot = pd.concat(method_rankings).pivot(
            index="candidate", columns="analysis", values="multiomic_rank"
        )
        ax = pivot.plot(marker="o", figsize=(10, 6))
        ax.invert_yaxis()
        ax.set_ylabel("Multi-omic rank (lower is better)")
        ax.set_title("Rank sensitivity to latent-centroid aggregation")
        plt.tight_layout()
        plt.savefig(figures / "robust_centroid_rank_stability.png", dpi=180)
        plt.close()
    finally:
        rna.file.close()
        atac.file.close()
    LOGGER.info("Phase 7 sensitivity outputs written to %s", tables)


if __name__ == "__main__":
    run_stage("07_run_sensitivity", main)
