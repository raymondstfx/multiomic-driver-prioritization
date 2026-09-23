"""Phase 7 robustness analyses that never overwrite the primary ranking."""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, trim_mean

from multiomic_driver.utils.design import DEFAULT_DESIGN

CONDITIONS = DEFAULT_DESIGN.conditions
REPLICATES = DEFAULT_DESIGN.replicates


def _boolean(values: pd.Series) -> np.ndarray:
    if pd.api.types.is_bool_dtype(values.dtype):
        return values.fillna(False).to_numpy(dtype=bool)
    return values.astype(str).str.lower().eq("true").to_numpy()


def _centroid(
    block: np.ndarray,
    method: str,
    *,
    winsor_fraction: float,
    trim_fraction: float,
) -> np.ndarray:
    if method == "mean":
        return block.mean(axis=0)
    if method == "winsorised_mean":
        lower = np.quantile(block, winsor_fraction, axis=0)
        upper = np.quantile(block, 1.0 - winsor_fraction, axis=0)
        return np.clip(block, lower, upper).mean(axis=0)
    if method == "trimmed_mean":
        return np.asarray(trim_mean(block, trim_fraction, axis=0))
    if method == "median":
        return np.median(block, axis=0)
    raise ValueError(f"Unsupported latent-centroid method: {method}")


def latent_centroid_effects(
    matrix: np.ndarray,
    metadata: pd.DataFrame,
    candidates: Iterable[str],
    *,
    min_cells_per_group: int,
    method: str = "mean",
    winsor_fraction: float = 0.01,
    trim_fraction: float = 0.05,
    allowed_guides: Mapping[str, set[str]] | None = None,
    allowed_ntc_targets: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute robust replicate-specific effects and candidate summaries."""
    if matrix.shape[0] != len(metadata):
        raise ValueError("Matrix rows must match metadata rows")
    singlet = metadata["guide_assignment_status"].astype(str).eq("singlet").to_numpy()
    ntc = _boolean(metadata["is_non_targeting"])
    target = metadata["target_gene"].astype(str).to_numpy()
    guide = metadata["guide"].astype(str).to_numpy()
    condition = metadata["condition"].astype(str).to_numpy()
    replicate = metadata["replicate"].astype(int).to_numpy()
    replicate_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for candidate in candidates:
        effects = []
        for rep in REPLICATES:
            vectors: dict[tuple[str, str], np.ndarray] = {}
            counts: dict[tuple[str, str], int] = {}
            for cond in CONDITIONS:
                base = singlet & (condition == cond) & (replicate == rep)
                candidate_mask = base & (target == candidate) & ~ntc
                if allowed_guides is not None:
                    candidate_mask &= np.isin(
                        guide, sorted(allowed_guides.get(candidate, set()))
                    )
                ntc_mask = base & ntc
                if allowed_ntc_targets is not None:
                    ntc_mask &= np.isin(target, sorted(allowed_ntc_targets))
                for population, mask in (
                    ("candidate", candidate_mask),
                    ("NTC", ntc_mask),
                ):
                    count = int(mask.sum())
                    if count < min_cells_per_group:
                        raise ValueError(
                            f"{candidate}:{cond}:R{rep}:{population} has {count} "
                            f"cells; requires {min_cells_per_group}"
                        )
                    vectors[(cond, population)] = _centroid(
                        matrix[mask],
                        method,
                        winsor_fraction=winsor_fraction,
                        trim_fraction=trim_fraction,
                    )
                    counts[(cond, population)] = count
            effect = (
                vectors[("Dasatinib", "candidate")]
                - vectors[("Dasatinib", "NTC")]
                - vectors[("DMSO", "candidate")]
                + vectors[("DMSO", "NTC")]
            )
            effects.append(effect)
            row: dict[str, object] = {
                "candidate": candidate,
                "replicate": rep,
                "aggregation_method": method,
                "effect_norm": float(np.linalg.norm(effect)),
                "dmso_candidate_cells": counts[("DMSO", "candidate")],
                "dmso_ntc_cells": counts[("DMSO", "NTC")],
                "dasatinib_candidate_cells": counts[("Dasatinib", "candidate")],
                "dasatinib_ntc_cells": counts[("Dasatinib", "NTC")],
            }
            row.update({f"dim_{i + 1}": value for i, value in enumerate(effect)})
            replicate_rows.append(row)
        stacked = np.vstack(effects)
        mean_effect = stacked.mean(axis=0)
        denominator = float(np.linalg.norm(stacked[0]) * np.linalg.norm(stacked[1]))
        cosine = (
            float(np.dot(stacked[0], stacked[1]) / denominator)
            if denominator > 0
            else float("nan")
        )
        row = {
            "candidate": candidate,
            "aggregation_method": method,
            "score": float(np.linalg.norm(mean_effect)),
            "replicate_cosine": cosine,
        }
        row.update({f"dim_{i + 1}": value for i, value in enumerate(mean_effect)})
        summary_rows.append(row)
    return pd.DataFrame(replicate_rows), pd.DataFrame(summary_rows)


def combine_sensitivity_scores(
    rna_scores: pd.Series,
    atac_scores: pd.Series,
    *,
    rna_weight: float = 0.5,
) -> pd.DataFrame:
    """Combine alternative modality scores using the fixed z-score definition."""
    if not 0.0 <= rna_weight <= 1.0:
        raise ValueError("rna_weight must be between zero and one")
    joined = pd.concat(
        [rna_scores.rename("rna_score"), atac_scores.rename("atac_score")], axis=1
    ).dropna()
    for modality in ("rna", "atac"):
        score = joined[f"{modality}_score"]
        sd = score.std(ddof=0)
        if np.isclose(sd, 0.0):
            raise ValueError(f"{modality} sensitivity scores have zero variance")
        joined[f"{modality}_z"] = (score - score.mean()) / sd
        joined[f"{modality}_rank"] = score.rank(method="min", ascending=False).astype(
            int
        )
    joined["rna_weight"] = rna_weight
    joined["atac_weight"] = 1.0 - rna_weight
    joined["multiomic_score"] = (
        rna_weight * joined["rna_z"] + (1.0 - rna_weight) * joined["atac_z"]
    )
    joined["multiomic_rank"] = (
        joined["multiomic_score"].rank(method="min", ascending=False).astype(int)
    )
    joined = joined.sort_values(
        ["multiomic_rank", "multiomic_score"], ascending=[True, False]
    )
    joined.insert(0, "candidate", joined.index)
    return joined.reset_index(drop=True)


def component_depth_diagnostics(
    embedding: np.ndarray,
    total_counts: np.ndarray,
    detected_features: np.ndarray,
    *,
    prefix: str,
    singular_values: np.ndarray | None = None,
    explained_variance_ratio: np.ndarray | None = None,
) -> pd.DataFrame:
    """Summarise component distributions and correlation with library depth."""
    rows = []
    quantiles = (0.001, 0.01, 0.5, 0.99, 0.999)
    for index in range(embedding.shape[1]):
        values = embedding[:, index]
        row = {
            "component": f"{prefix}{index + 1}",
            "pearson_with_total_counts": float(np.corrcoef(values, total_counts)[0, 1]),
            "spearman_with_total_counts": float(
                spearmanr(values, total_counts).statistic
            ),
            "pearson_with_detected_features": float(
                np.corrcoef(values, detected_features)[0, 1]
            ),
            "spearman_with_detected_features": float(
                spearmanr(values, detected_features).statistic
            ),
            "minimum": float(values.min()),
            "q001": float(np.quantile(values, quantiles[0])),
            "q01": float(np.quantile(values, quantiles[1])),
            "median": float(np.quantile(values, quantiles[2])),
            "q99": float(np.quantile(values, quantiles[3])),
            "q999": float(np.quantile(values, quantiles[4])),
            "maximum": float(values.max()),
        }
        if singular_values is not None:
            row["singular_value"] = float(singular_values[index])
        if explained_variance_ratio is not None:
            row["explained_variance_ratio"] = float(explained_variance_ratio[index])
        rows.append(row)
    return pd.DataFrame(rows)


def effect_dimension_contributions(
    effect_summary: pd.DataFrame, *, prefix: str, modality: str
) -> pd.DataFrame:
    """Return each latent dimension's contribution to squared effect magnitude."""
    dimensions = [column for column in effect_summary if column.startswith(prefix)]
    rows = []
    for row in effect_summary.itertuples(index=False):
        values = np.asarray(
            [getattr(row, column) for column in dimensions], dtype=float
        )
        squared = values**2
        total = squared.sum()
        for dimension, value, contribution in zip(
            dimensions, values, squared, strict=True
        ):
            rows.append(
                {
                    "candidate": row.target_gene,
                    "modality": modality,
                    "dimension": dimension,
                    "effect_value": value,
                    "squared_effect": contribution,
                    "fraction_of_squared_norm": contribution / total
                    if total
                    else np.nan,
                }
            )
    return pd.DataFrame(rows)


def guide_composition(metadata: pd.DataFrame) -> pd.DataFrame:
    """Count exact guide-call composition across experimental groups."""
    singlets = metadata.loc[
        metadata["guide_assignment_status"].astype(str).eq("singlet")
        & metadata["target_gene"].notna()
    ].copy()
    counts = (
        singlets.groupby(
            ["target_gene", "guide", "condition", "replicate"], observed=True
        )
        .size()
        .rename("cell_count")
        .reset_index()
    )
    totals = counts.groupby(["target_gene", "condition", "replicate"], observed=True)[
        "cell_count"
    ].transform("sum")
    counts["within_target_group_fraction"] = counts["cell_count"] / totals
    observed = counts.groupby(["target_gene", "guide"], observed=True).size()
    counts = counts.merge(
        observed.rename("observed_experimental_groups"),
        on=["target_gene", "guide"],
    )
    counts["present_in_all_four_groups"] = counts["observed_experimental_groups"].eq(4)
    return counts.sort_values(
        ["target_gene", "guide", "condition", "replicate"]
    ).reset_index(drop=True)


def common_guide_sets(
    composition: pd.DataFrame, *, min_cells_per_group: int
) -> dict[str, set[str]]:
    """Find exact guide-call classes meeting the minimum in all four groups."""
    groups = composition.pivot_table(
        index=["target_gene", "guide"],
        columns=["condition", "replicate"],
        values="cell_count",
        fill_value=0,
        observed=True,
    )
    expected = pd.MultiIndex.from_product([CONDITIONS, REPLICATES])
    groups = groups.reindex(columns=expected, fill_value=0)
    eligible = groups.min(axis=1).ge(min_cells_per_group)
    result: dict[str, set[str]] = {}
    for target, guide in groups.index[eligible]:
        result.setdefault(str(target), set()).add(str(guide))
    return result
