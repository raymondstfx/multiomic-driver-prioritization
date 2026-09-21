# Phase 3 preprocessing QC

Phase 3 was executed on the 53,898 aligned cells produced by Phase 2.

## RNA

- Cells: 53,898 before and after QC (minimum 200 detected genes).
- Genes: 36,601 input; 28,632 detected in at least three cells.
- Median counts per cell: 7,055.
- Median detected genes per cell: 2,863.
- Mitochondrial annotations: 13 genes detected; median mitochondrial fraction
  16.52%. This metric is reported but was not used for filtering because the
  execution plan did not specify a supported threshold.
- Normalisation: target sum 10,000 followed by log1p.
- Raw filtered counts: preserved in `layers["counts"]`.
- Representation: 3,000 HVGs and 30 PCs.

## ATAC

The 792,914 exact input peak coordinates were almost entirely sample-specific:
792,017 occurred in one condition × replicate group, 897 occurred in two or
three groups, and none occurred in all four. Direct LSI on that exact-coordinate
union separated the four samples and was rejected as a technical artefact.

Overlapping existing peak calls were therefore collapsed into 278,457 shared,
non-overlapping consensus intervals. This is coordinate harmonisation, not new
peak calling; no fragments were used. All consensus intervals occurred in at
least three cells. Of these intervals, 133,983 were detected in all four sample
groups. Sparse TF-IDF used
`TF_ij * log(1 + N / (1 + DF_j))`, followed by 30-component truncated SVD.
LSI1 was retained. The corrected diagnostic plots no longer show mutually
exclusive sample axes, although several high-LSI2 outliers remain and should
not be interpreted as biological mechanisms without further QC.

## Cell and perturbation retention

RNA and ATAC each retained all 53,898 cells, with all cell IDs in the same order
and all 53,898 shared between modalities. Counts by group were:

| Condition | Replicate | Cells |
|---|---:|---:|
| DMSO | 1 | 12,621 |
| DMSO | 2 | 13,547 |
| Dasatinib | 1 | 11,920 |
| Dasatinib | 2 | 15,810 |

HIC2 and pooled non-targeting-control singlets were retained identically in
both modalities:

| Condition | Replicate | HIC2 | NTC |
|---|---:|---:|---:|
| DMSO | 1 | 183 | 757 |
| DMSO | 2 | 230 | 994 |
| Dasatinib | 1 | 325 | 1,316 |
| Dasatinib | 2 | 320 | 1,476 |

ZFPM2 remains an RNA feature but is not a pooled guide target. It is excluded
from perturbation ranking and reserved for downstream/external validation.

## Verification

Both processed matrices and their 30-dimensional latent representations contain
only finite values. Source metadata values and cell ordering were preserved.
The synthetic test suite passes (16 tests), and Ruff reports no violations.
