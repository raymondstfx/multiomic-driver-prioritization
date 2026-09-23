"""Build committed, human-facing summaries and a machine-readable manifest."""

from __future__ import annotations

import logging
import platform
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from _common import run_stage

from multiomic_driver.utils.artifacts import (
    PHASE_DIRECTORIES,
    figure_dir,
    file_sha256,
    results_root,
    summary_dir,
    table_dir,
)
from multiomic_driver.utils.io import require_path

LOGGER = logging.getLogger("08_finalize_results")

DESCRIPTIONS = {
    "candidate_ranking.csv": "Primary fixed Phase 5 candidate ranking.",
    "final_candidate_evaluation.csv": "Phase 6 annotations of the fixed ranking.",
    "zfpm2_validation_summary.csv": "Metric-specific downstream ZFPM2 validation.",
    "sensitivity_summary.csv": "Robust-centroid rank stability summary.",
    "latent_component_diagnostics.csv": "Latent components versus depth and distribution diagnostics.",
}


def _copy(source: Path, destination: Path) -> None:
    require_path(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def _git_commit(root: Path) -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else "unavailable"


def main(config: dict) -> None:
    root = results_root(config)
    summary = summary_dir(config)
    phase05 = table_dir(config, "phase05")
    phase06 = table_dir(config, "phase06")
    phase07 = table_dir(config, "phase07")
    phase035 = table_dir(config, "phase035")
    version = str(config["sensitivity"]["primary_analysis_version"])

    ranking_source = phase05 / "candidate_ranking.csv"
    evaluation_source = phase06 / "final_candidate_evaluation.csv"
    zfpm2_source = phase06 / "zfpm2_validation_summary.csv"
    sensitivity_source = phase07 / "sensitivity_summary.csv"
    _copy(evaluation_source, summary / "final_candidate_evaluation.csv")
    _copy(sensitivity_source, summary / "final_sensitivity_summary.csv")

    for filename in ("candidate_ranking.png", "modality_rank_comparison.png"):
        _copy(
            figure_dir(config, "phase05") / filename,
            summary / "figures" / filename,
        )
    for filename in ("top5_candidate_evidence.png", "zfpm2_hic2_vs_ntc.png"):
        _copy(
            figure_dir(config, "phase06") / filename,
            summary / "figures" / filename,
        )
    for filename in ("robust_centroid_rank_stability.png",):
        _copy(figure_dir(config, "phase07") / filename, summary / "figures" / filename)

    ranking = pd.read_csv(ranking_source)
    comparison = pd.read_csv(phase06 / "modality_rank_comparison.csv")
    reliability = pd.read_csv(phase06 / "replicate_reliability_summary.csv")
    eligibility = pd.read_csv(phase035 / "candidate_phase4_eligibility.csv")
    literature = pd.read_csv(phase06 / "top_candidate_literature_evidence.csv")
    canonical = (
        ranking.merge(
            comparison[
                [
                    "candidate",
                    "cross_modal_profile",
                    "rna_to_multiomic_shift",
                    "atac_to_multiomic_shift",
                ]
            ],
            on="candidate",
            validate="one_to_one",
        )
        .merge(
            reliability[
                [
                    "candidate",
                    "rna_consistency_category",
                    "atac_consistency_category",
                ]
            ],
            on="candidate",
            validate="one_to_one",
        )
        .merge(
            eligibility[["target_gene", "min_group_cells"]],
            left_on="candidate",
            right_on="target_gene",
            validate="one_to_one",
        )
        .drop(columns="target_gene")
    )
    canonical["experimental_validation_status"] = canonical["candidate"].map(
        lambda value: (
            "supported benchmark" if value == "HIC2" else "not evaluated"
        )
    )
    evidence = literature.set_index("candidate")["evidence_category"]
    canonical["literature_evidence_category"] = canonical["candidate"].map(
        evidence
    ).fillna("not included in focused review")
    canonical["primary_analysis_version"] = version
    canonical.to_csv(summary / "final_candidate_ranking.csv", index=False)

    validation_rows = []
    hic2_row = ranking.loc[ranking["candidate"].eq("HIC2")].iloc[0]
    for modality in ("rna", "atac", "multiomic"):
        validation_rows.append(
            {
                "validation_id": f"hic2_{modality}_rank",
                "target": "HIC2",
                "metric": f"{modality}_rank",
                "replicate": "combined",
                "value": hic2_row[f"{modality}_rank"],
                "direction": "lower_is_better",
                "interpretation": "post-ranking supported benchmark",
                "analysis_role": "primary_validation",
            }
        )
    for modality in ("rna", "atac"):
        value = hic2_row[f"{modality}_replicate_cosine"]
        validation_rows.append(
            {
                "validation_id": f"hic2_{modality}_replicate_cosine",
                "target": "HIC2",
                "metric": f"{modality}_replicate_cosine",
                "replicate": "R1_vs_R2",
                "value": value,
                "direction": "higher_is_more_consistent",
                "interpretation": "diagnostic; not included in primary score",
                "analysis_role": "diagnostic",
            }
        )
    zfpm2_replicates = pd.read_csv(phase06 / "zfpm2_replicate_did.csv")
    for row in zfpm2_replicates.itertuples(index=False):
        validation_rows.append(
            {
                "validation_id": f"zfpm2_{row.metric}_R{row.replicate}",
                "target": "ZFPM2 downstream of HIC2",
                "metric": row.metric,
                "replicate": f"R{row.replicate}",
                "value": row.replicate_specific_did,
                "direction": "positive" if row.replicate_specific_did > 0 else "negative",
                "interpretation": "within-dataset downstream consistency; not direct regulation",
                "analysis_role": "secondary_validation",
            }
        )
    pd.DataFrame(validation_rows).to_csv(
        summary / "final_validation_summary.csv", index=False
    )

    hic2 = ranking.loc[ranking["candidate"].eq("HIC2")].iloc[0]
    sensitivity = pd.read_csv(sensitivity_source)
    zfpm2 = pd.read_csv(zfpm2_source)
    robust = sensitivity.set_index("analysis")
    zfpm2_lines = []
    for row in zfpm2.itertuples(index=False):
        zfpm2_lines.append(
            f"- `{row.metric}`: mean replicate DiD = {row.mean_replicate_did:.4g}; "
            f"direction agreement = {bool(row.replicate_direction_agreement)}."
        )
    summary_text = f"""# Pipeline summary

Generated from analysis version `{version}` at {datetime.now(timezone.utc).isoformat()}.
Git commit: `{_git_commit(Path(config["_project_root"]))}`. Analysis-config SHA-256:
`{file_sha256(Path(config["_project_root"]) / "configs" / "analysis.yaml")}`.
Data-config SHA-256: `{file_sha256(Path(config["_project_root"]) / "configs" / "data.yaml")}`.
Runtime: Python `{platform.python_version()}` on `{platform.platform()}`. Expected
guide-call archive SHA-256: `{config["guide_calls"]["sha256"]}`.

Dataset: GSE288996, 53,898 paired K562 cells; DMSO and Dasatinib; biological
replicates R1 and R2. Twelve perturbations pass the configured 20-cell minimum in
all four groups; missing groups remain NA and CNOT2 remains below threshold.

## Primary result

The primary Phase 5 ranking is unchanged by Phase 7. The top five candidates are
{", ".join(ranking.nsmallest(5, "display_order")["candidate"])}. HIC2 is RNA rank
{int(hic2.rna_rank)}, ATAC rank {int(hic2.atac_rank)}, and combined rank
{int(hic2.multiomic_rank)}.

For candidate $i$, replicate $r$, and modality $m$, the primary interaction is

$$
\\mathbf E_{{i,r}}^{{(m)}}=(\\bar{{\\mathbf x}}_{{i,r,DASA}}^{{(m)}}-
\\bar{{\\mathbf x}}_{{NTC,r,DASA}}^{{(m)}})-(\\bar{{\\mathbf x}}_{{i,r,DMSO}}^{{(m)}}-
\\bar{{\\mathbf x}}_{{NTC,r,DMSO}}^{{(m)}}).
$$

The score is the Euclidean norm of the replicate-mean effect and is relative to
the current eligible candidate cohort. Its magnitude does not encode a biological
resistance direction and is not causal proof.

## Robustness

The primary estimator is the mean latent-space centroid. Alternative estimators are
reported as sensitivity analyses, not replacements. Their Spearman correlations with
the primary combined rank range from {robust.rank_spearman_vs_primary.min():.3f} to
{robust.rank_spearman_vs_primary.max():.3f}. HIC2 ranks from
{int(robust.hic2_rank.min())} to {int(robust.hic2_rank.max())} across these estimators.
The detailed tables also cover leave-one-dimension-out, NTC-SD scaling, guide
composition, replicate-specific ranks, cell thresholds, and RNA/ATAC weights.

## ZFPM2 downstream validation

ZFPM2 is absent from the perturbation library and is not ranked. HIC2-versus-NTC
ZFPM2 responses are descriptive and metric-dependent:

{chr(10).join(zfpm2_lines)}

These results do not establish direct HIC2-to-ZFPM2 regulation or causal mechanism.

## Navigation

- `final_candidate_ranking.csv`: canonical primary result.
- `final_candidate_evaluation.csv`: Phase 6 interpretation and diagnostics.
- `final_validation_summary.csv`: multi-metric ZFPM2 validation.
- `final_sensitivity_summary.csv`: compact robustness summary.
- `../manifest.csv`: provenance, dimensions, hashes, and roles for all artifacts.
"""
    (summary / "pipeline_summary.md").write_text(summary_text, encoding="utf-8")

    generated = datetime.now(timezone.utc).isoformat()
    repo_root = Path(config["_project_root"])
    rows = []
    phase_lookup = {value: key for key, value in PHASE_DIRECTORIES.items()}
    for path in sorted(root.rglob("*")):
        if (
            not path.is_file()
            or path.name == "manifest.csv"
            or path.name == ".gitkeep"
            or "logs" in path.parts
            or "legacy_flat" in path.parts
        ):
            continue
        relative = path.relative_to(root).as_posix()
        phase = "summary"
        for part in path.parts:
            if part in phase_lookup:
                phase = phase_lookup[part]
                break
        suffix = path.suffix.lower()
        artifact_type = (
            "table"
            if suffix == ".csv"
            else "figure"
            if suffix in {".png", ".pdf", ".svg"}
            else "document"
        )
        n_rows = n_columns = None
        if suffix == ".csv":
            frame = pd.read_csv(path)
            n_rows, n_columns = frame.shape
        is_primary = relative in {
            "summary/final_candidate_ranking.csv",
            "summary/final_candidate_evaluation.csv",
            "summary/final_validation_summary.csv",
            "summary/pipeline_summary.md",
        }
        rows.append(
            {
                "phase": phase,
                "artifact_path": relative,
                "artifact_type": artifact_type,
                "role": "primary" if is_primary else (
                    "sensitivity" if phase == "phase07" else "diagnostic"
                ),
                "primary": is_primary,
                "primary_analysis_version": version,
                "generated_by": "scripts/08_finalize_results.py"
                if phase == "summary"
                else f"phase_{phase}",
                "description": DESCRIPTIONS.get(path.name, path.stem.replace("_", " ")),
                "row_count": n_rows,
                "column_count": n_columns,
                "sha256": file_sha256(path),
                "created_at": generated,
                "git_commit": _git_commit(repo_root),
            }
        )
    pd.DataFrame(rows).to_csv(root / "manifest.csv", index=False)
    LOGGER.info(
        "Final summaries and %d manifest entries written under %s", len(rows), root
    )


if __name__ == "__main__":
    run_stage("08_finalize_results", main)
