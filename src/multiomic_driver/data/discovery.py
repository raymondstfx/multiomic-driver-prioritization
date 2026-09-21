"""Discover and validate the flat GEO supplementary-file layout."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from multiomic_driver.data.metadata import parse_sample_name


@dataclass(frozen=True)
class ExperimentFiles:
    """Files for one matched condition-replicate multiome experiment."""

    condition: str
    replicate: int
    rna_sample_id: str
    atac_sample_id: str
    guide_sample_id: str
    rna_matrix: Path
    rna_barcodes: Path
    features: Path
    atac_matrix: Path
    atac_barcodes: Path
    fragments: Path | None
    guide_summary: Path

    @property
    def group_key(self) -> str:
        """Return the author repository's condition-replicate key."""
        condition = "DASA" if self.condition == "Dasatinib" else "DMSO"
        return f"{condition}{self.replicate}"


def _single(paths: list[Path], description: str) -> Path:
    if len(paths) != 1:
        rendered = ", ".join(path.name for path in paths) or "none"
        raise ValueError(f"Expected one {description}, found {len(paths)}: {rendered}")
    return paths[0]


def discover_experiments(raw_dir: str | Path) -> list[ExperimentFiles]:
    """Discover the four K562 experiments and their 12 GEO samples."""
    root = Path(raw_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"Raw dataset directory does not exist: {root}")

    records: dict[tuple[str, int, str], dict[str, object]] = {}
    for path in root.iterdir():
        if not path.is_file() or "K562" not in path.name:
            continue
        try:
            metadata = parse_sample_name(path)
        except ValueError:
            continue
        key = (metadata["condition"], metadata["replicate"], metadata["modality"])
        record = records.setdefault(key, {"metadata": metadata, "paths": []})
        record["paths"].append(path)  # type: ignore[union-attr]

    expected = {
        (condition, replicate, modality)
        for condition in ("Dasatinib", "DMSO")
        for replicate in (1, 2)
        for modality in ("RNA", "ATAC", "guideRNA")
    }
    missing = expected - set(records)
    if missing:
        raise ValueError(f"Missing expected K562 samples: {sorted(missing)}")

    experiments: list[ExperimentFiles] = []
    for condition in ("Dasatinib", "DMSO"):
        for replicate in (1, 2):
            rna = records[(condition, replicate, "RNA")]
            atac = records[(condition, replicate, "ATAC")]
            guide = records[(condition, replicate, "guideRNA")]
            rna_paths = rna["paths"]
            atac_paths = atac["paths"]
            guide_paths = guide["paths"]
            assert isinstance(rna_paths, list)
            assert isinstance(atac_paths, list)
            assert isinstance(guide_paths, list)
            experiments.append(
                ExperimentFiles(
                    condition=condition,
                    replicate=replicate,
                    rna_sample_id=rna["metadata"]["sample_id"],  # type: ignore[index]
                    atac_sample_id=atac["metadata"]["sample_id"],  # type: ignore[index]
                    guide_sample_id=guide["metadata"]["sample_id"],  # type: ignore[index]
                    rna_matrix=_single(
                        [path for path in rna_paths if "matrix.mtx" in path.name],
                        "RNA matrix",
                    ),
                    rna_barcodes=_single(
                        [path for path in rna_paths if "barcodes.tsv" in path.name],
                        "RNA barcode file",
                    ),
                    features=_single(
                        [path for path in rna_paths if "features.tsv" in path.name],
                        "multiome feature file",
                    ),
                    atac_matrix=_single(
                        [path for path in atac_paths if "matrix.mtx" in path.name],
                        "ATAC matrix",
                    ),
                    atac_barcodes=_single(
                        [path for path in atac_paths if "barcodes.tsv" in path.name],
                        "ATAC barcode file",
                    ),
                    fragments=next(
                        (path for path in atac_paths if "fragments.tsv" in path.name),
                        None,
                    ),
                    guide_summary=_single(guide_paths, "guide count summary"),
                )
            )
    return experiments

