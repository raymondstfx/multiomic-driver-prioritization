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
7. Evaluate HIC2 post hoc as a rankable library perturbation.
8. Treat ZFPM2 only as downstream/external biological validation; it is not a
   pooled perturbation and cannot receive a perturbation rank.

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
python scripts/03b_phase4_readiness.py --config configs/default.yaml
python scripts/04_compute_effects.py --config configs/default.yaml
python scripts/05_rank_candidates.py --config configs/default.yaml
python scripts/06_evaluate.py --config configs/default.yaml
python scripts/06b_zfpm2_external_validation.py --config configs/default.yaml
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

Phase 3 writes `rna_processed.h5ad` and `atac_processed.h5ad` under
`data/processed/`. RNA retains raw counts in `layers["counts"]`, stores
log-normalised expression in `X`, and stores 30 PCs in `obsm["X_pca"]`. Because
the four Cell Ranger matrices use independently called peak coordinates, ATAC
first collapses overlapping called peaks into shared consensus intervals (this
does not call new peaks or use fragments). It then stores sparse TF-IDF in `X`
and 30 LSI dimensions in `obsm["X_lsi"]`; LSI1 is retained and documented. QC
tables and condition/replicate diagnostic plots are written to `results/`.

The Phase 3.5 readiness stage characterises mitochondrial-count distributions,
requires each primary targeting perturbation to have at least 20 singlet cells
in every condition × replicate group, verifies NTC/HIC2 coverage, and writes an
explicit eligible-candidate list. Missing perturbation groups remain missing;
they are never converted into zero effects.

Phase 4 keeps biological replicates separate during estimation. For each
candidate and modality it computes an NTC-corrected Difference-in-Differences
vector independently for R1 and R2. Only afterwards are the two vectors
averaged. Replicate cosine similarity is reported as a separate consistency
metric; it is not folded into an effect or ranking score at this stage.

Phase 5 ranks candidates by the magnitude of each replicate-averaged effect
vector. RNA and ATAC magnitudes are standardized independently using population
SD (`ddof=0`), then combined with equal weight. Replicate consistency remains a
separate diagnostic. HIC2 is inspected only after the label-free ranking has
been finalized; ZFPM2 is not a ranked perturbation.

Phase 6 leaves all Phase 5 scores and ranks unchanged. It compares modality
ranks, assigns sign-based cross-modal profiles, and reports replicate cosine
categories as separate annotations. The fixed top five are YEATS4, GPBP1L1,
ZBED6, HIC2, and KMT2B. HIC2 ranks second by RNA, fifth by ATAC, and fourth in
the combined ranking; integration therefore does not outperform RNA alone for
this supported candidate.

The targeted external-validation analysis reads log-normalised ZFPM2 gene
expression from the RNA object and compares HIC2 and NTC singlets within every
condition and replicate. Replicate-specific ZFPM2 Difference-in-Differences
effects are calculated before their mean is reported. This descriptive result
is downstream biological consistency, not a ZFPM2 perturbation score, direct
HIC2-to-ZFPM2 regulation, or causal proof.

## Outputs

Generated tables are written to `results/tables/`, figures to
`results/figures/`, and logs to `results/logs/`. These outputs are excluded
from Git apart from directory placeholders.

## Limitations

Guide calls follow the authors' published singlet rule: retain one called guide,
or two guides only when both target the same vector. Cells with other positive
combinations are labelled as multiplets and are not assigned a target. The
score remains evidence prioritisation, not causal proof. HIC2 is used only for
post-ranking evaluation. ZFPM2 is absent from the guide library and is reserved
for separate downstream biological validation.

This proof-of-concept uses one K562 dataset, one Dasatinib/DMSO comparison, two
biological replicates, and 12 primary perturbations passing a configurable
20-cell-per-group threshold. Scores are latent-space effect magnitudes with
simple equal modality weights; there are no learned weights, full regulatory
network reconstruction, MPL implementation, or selection-coefficient
estimates. A multi-omic rank does not necessarily improve on the best single
modality for each known candidate.

## Future work

Potential extensions include threshold-sensitivity analysis, replicate-
stability weighting, alternative modality weights, feature-level RNA-ATAC
linking, gene/peak mapping, pathway analysis, GRN inference, added drugs, cell
lines and perturbation screens, MultiVI-like joint models, and MPL-inspired
selection modelling. These are not implemented in Phase 6.
