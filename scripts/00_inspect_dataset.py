"""Inspect real GSE288996 inputs without performing preprocessing."""

from pathlib import Path

from _common import run_stage

from multiomic_driver.data.align import (
    barcode_overlap_summary,
    classify_guide_calls,
)
from multiomic_driver.data.discovery import discover_experiments
from multiomic_driver.data.loader import (
    load_guide_call_matrix,
    load_guide_count_summary,
    read_barcodes,
    read_features,
    validate_matrix_dimensions,
)
from multiomic_driver.utils.io import require_path


def main(config: dict) -> None:
    raw = require_path(config["paths"]["raw"], "raw GSE288996 directory")
    archive = require_path(config["paths"]["guide_calls"], "official guide calls")
    experiments = discover_experiments(raw)
    files = sorted(path for path in Path(raw).rglob("*") if path.is_file())
    total_bytes = sum(path.stat().st_size for path in files)

    logs = Path(config["paths"]["results"]) / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    inventory = logs / "dataset_files.txt"
    inventory.write_text(
        "\n".join(str(path.relative_to(raw)) for path in files) + "\n",
        encoding="utf-8",
    )

    lines = [
        "GSE288996 Dataset Inspection",
        "============================",
        f"Total dataset size: {total_bytes} bytes ({total_bytes / 1024**3:.3f} GiB)",
        f"Total number of files: {len(files)}",
        "",
        "K562 samples discovered:",
        "RNA: 4",
        "ATAC: 4",
        "guideRNA: 4",
    ]
    for experiment in experiments:
        rna_barcodes = read_barcodes(experiment.rna_barcodes)
        atac_barcodes = read_barcodes(experiment.atac_barcodes)
        features = read_features(experiment.features)
        rna_shape = validate_matrix_dimensions(
            experiment.rna_matrix, rna_barcodes, features
        )
        atac_shape = validate_matrix_dimensions(
            experiment.atac_matrix, atac_barcodes, features
        )
        guide_summary = load_guide_count_summary(experiment.guide_summary)
        calls = load_guide_call_matrix(archive, experiment.group_key)
        assignments = classify_guide_calls(calls)
        assigned = assignments.index[
            assignments["guide_assignment_status"].eq("singlet")
        ]
        overlap = barcode_overlap_summary(rna_barcodes, atac_barcodes, assigned)
        feature_counts = features["feature_type"].value_counts().to_dict()
        status_counts = assignments["guide_assignment_status"].value_counts().to_dict()
        called_examples = assignments.loc[
            assignments["guide_assignment_status"].eq("singlet"),
            ["guide", "target_gene"],
        ].head(3)

        lines.extend(
            [
                "",
                "-" * 60,
                f"{experiment.condition} - replicate {experiment.replicate}",
                "-" * 60,
                "RNA / shared ARC matrix",
                f"matrix: {experiment.rna_matrix.name}",
                f"shape (features x cells): {rna_shape[:2]}",
                f"nonzero entries: {rna_shape[2]}",
                f"barcodes: {len(rna_barcodes)}",
                f"features: {len(features)} {feature_counts}",
                "ATAC",
                f"matrix: {experiment.atac_matrix.name}",
                f"shape (features x cells): {atac_shape[:2]}",
                f"barcodes: {len(atac_barcodes)}",
                "feature/peak annotations: shared RNA feature table",
                "guideRNA",
                f"GEO aggregate file: {experiment.guide_summary.name}",
                f"aggregate guide rows: {len(guide_summary)}",
                f"official call matrix: {calls.shape[0]} cells x {calls.shape[1]} guides",
                f"assignment status: {status_counts}",
                "Barcode examples",
                f"RNA: {rna_barcodes[:3].tolist()}",
                f"ATAC: {atac_barcodes[:3].tolist()}",
                f"guideRNA: {calls.index[:3].tolist()}",
                "Guide assignment examples",
                called_examples.to_string(),
                "Overlap",
                f"RNA cells: {overlap['rna_cells']}",
                f"ATAC cells: {overlap['atac_cells']}",
                f"Guide-assigned cells: {overlap['guide_cells']}",
                f"RNA intersection ATAC: {overlap['rna_atac']}",
                f"RNA intersection guideRNA: {overlap['rna_guide']}",
                f"ATAC intersection guideRNA: {overlap['atac_guide']}",
                (
                    "RNA intersection ATAC intersection guideRNA: "
                    f"{overlap['rna_atac_guide']}"
                ),
            ]
        )
    lines.extend(
        [
            "",
            "Barcode normalization: surrounding whitespace only; '-1' is preserved.",
            "Guide rule: one call, or two calls to the same target, is a singlet.",
            "Non-targeting controls: official targets whose identifiers start with NTC.",
            "GEO guideRNA TXT files are aggregate counts, not cell assignments.",
        ]
    )
    report = "\n".join(lines) + "\n"
    (logs / "dataset_inspection.txt").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    run_stage("00_inspect_dataset", main)
