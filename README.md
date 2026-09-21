# Background-Aware Multi-omic Driver Prioritisation

## Overview

This repository is a reproducible proof-of-concept for prioritising candidate
molecular drivers from single-cell perturbation multi-omics. It asks whether
treatment-specific effects that agree across RNA and chromatin accessibility
and reproduce across biological replicates provide stronger driver evidence
than a single modality alone.

The method uses MPL-inspired background-aware reasoning, but it is **not** an
implementation of Marginal Path Likelihood, does not estimate selection
coefficients, and does not claim causal gene discovery.

## Dataset

The first dataset is GEO accession **GSE288996**: K562 cells under Dasatinib
and DMSO conditions with RNA, ATAC, and CRISPR guide information. This project
uses processed matrices; raw FASTQ and fragment-level processing are out of
scope. Large data files are intentionally excluded from Git.

Place downloaded GEO files under `data/raw/GSE288996/`. The GEO guideRNA TXT
files are aggregate guide counts rather than cell-level assignments. Download
the authors' official
[`guide_caller_calls.zip`](https://github.com/ucsf-lgr/catatac_public/raw/refs/heads/master/guide_caller_calls.zip)
to the same directory; its expected SHA-256 is recorded in `configs/data.yaml`.

## Method

1. Align RNA, ATAC, and guide assignments by cell barcode.
2. Apply simple modality-specific preprocessing.
3. Aggregate by perturbation, treatment, and biological replicate.
4. Estimate background-corrected treatment effects with Difference-in-Differences.
5. Calculate RNA and ATAC effect magnitudes.
6. Combine standardised scores into an interpretable multi-omic ranking.
7. Evaluate HIC2 and ZFPM2 post hoc; they are never ranking inputs.

For perturbation `i`, the core contrast is:

```text
E_i = (Dasatinib_i - Dasatinib_NT) - (DMSO_i - DMSO_NT)
```

## Repository structure

```text
configs/       Versioned data and analysis settings
data/          Ignored raw, interim, and processed data
docs/          Dataset and methodology notes
scripts/       Command-line pipeline entry points
src/           Reusable Python package
tests/         Unit tests for core logic
results/       Ignored generated tables, figures, and logs
```

## Installation

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
python -m pip install -e ".[dev]"
```

## Usage

Run each stage from the repository root:

```bash
python scripts/00_inspect_dataset.py --config configs/default.yaml
python scripts/01_prepare_data.py --config configs/default.yaml
python scripts/02_preprocess_rna.py --config configs/default.yaml
python scripts/03_preprocess_atac.py --config configs/default.yaml
python scripts/04_compute_effects.py --config configs/default.yaml
python scripts/05_rank_candidates.py --config configs/default.yaml
python scripts/06_evaluate.py --config configs/default.yaml
python scripts/07_generate_figures.py --config configs/default.yaml
```

The inspection stage records the file inventory, real matrix dimensions,
feature types, guide-call structure, barcode examples, and per-sample overlap
in `results/logs/dataset_inspection.txt`. The preparation stage creates:

```text
data/interim/cell_metadata.csv
data/interim/rna_merged.h5ad
data/interim/atac_merged.h5ad
```

RNA and ATAC are split from the shared Cell Ranger ARC matrix without dense
conversion. Barcode suffixes such as `-1` are preserved because they match
across the real files. Run tests with `pytest`.

## Outputs

Generated tables are written to `results/tables/`, figures to
`results/figures/`, and logs to `results/logs/`. These outputs are excluded
from Git apart from directory placeholders.

## Limitations

Guide calls follow the authors' published singlet rule: retain one called guide,
or two guides only when both target the same vector. Cells with other positive
combinations are labelled as multiplets and are not assigned a target. The
score remains evidence prioritisation, not causal proof, and conclusions based
on two externally supported candidates must remain modest.

## Future work

Potential extensions include replicate-stability weighting, gene-to-peak
linking, pathway-level influence analogues, and more advanced representations
after the basic alignment and background correction have been validated.
