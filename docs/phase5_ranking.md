# Phase 5 multi-omic candidate ranking

Phase 5 ranks exactly the 12 primary candidates approved in Phase 3.5. It uses
the verified Phase 4 `mean_effect_norm` values for RNA PCA and ATAC LSI; raw
modality magnitudes are not compared or added directly.

Each modality is standardized independently with population SD (`ddof=0`). The
primary multi-omic score is `(rna_z + atac_z) / 2`, without learned weights,
known-candidate labels, or replicate-cosine weighting. Negative scores mean
below-average evidence within this 12-candidate set, not a negative biological
effect.

## Results

| Multi-omic rank | Candidate | RNA rank | ATAC rank | Multi-omic score |
|---:|---|---:|---:|---:|
| 1 | YEATS4 | 1 | 3 | 1.778 |
| 2 | GPBP1L1 | 3 | 2 | 1.065 |
| 3 | ZBED6 | 10 | 1 | 0.611 |
| 4 | HIC2 | 2 | 5 | 0.573 |
| 5 | KMT2B | 4 | 4 | -0.066 |
| 6 | CHD2 | 5 | 6 | -0.319 |
| 7 | ADNP | 7 | 8 | -0.377 |
| 8 | PIAS1 | 9 | 7 | -0.407 |
| 9 | TSC22D4 | 6 | 9 | -0.414 |
| 10 | SIN3A | 8 | 10 | -0.579 |
| 11 | SLTM | 11 | 11 | -0.853 |
| 12 | MNT | 12 | 12 | -1.011 |

RNA and ATAC prioritization have a descriptive Spearman correlation of 0.601.
With only 12 candidates this is a compact concordance summary, not a formal
biological significance claim.

## Post-hoc HIC2 evaluation

HIC2 was inspected only after the ranking was complete. It ranks second in RNA,
fifth in ATAC, and fourth in the equal-weight multi-omic ranking. Its RNA and
ATAC replicate cosine similarities are 0.906 and 0.574 respectively. These
values were not used to tune the ranking method.

ZFPM2 is absent from the ranking because it is not a pooled perturbation. It
remains reserved for downstream/external validation.

## Replicate diagnostics

Raw replicate cosine values remain separate from ranking scores. In particular,
MNT has a slightly negative RNA cosine, while ADNP and SIN3A have negative ATAC
cosines. These disagreements remain visible and are neither absolutized nor
silently excluded.

Generated tables and figures under `results/` are ignored by Git. The primary
table is `candidate_ranking.csv`; the additional diagnostics include modality
ranks, replicate norms/cosines, the descriptive rank correlation, and five
ranking figures.
