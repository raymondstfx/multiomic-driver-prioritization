# Phase 7 sensitivity analysis

Phase 7 audits the fixed `primary_v1` ranking; it does not choose a replacement
estimator after inspecting HIC2. Recomputed mean-centroid scores match the Phase 5
primary scores within floating-point tolerance.

## Main findings

- Robust-centroid ranks have Spearman correlations of 0.755–1.000 with the primary
  combined ranking. HIC2 moves from primary rank 4 to rank 3 under 1% winsorisation
  and rank 2 under 5% trimming or median centroids.
- ATAC LSI2, rather than LSI1, is the main instability. Removing LSI1 leaves the
  ranking unchanged; removing LSI2 lowers rank correlation to 0.741 and moves HIC2
  to rank 2. LSI2 accounts for 97.1%, 96.3%, and 85.3% of squared ATAC effect
  magnitude for ZBED6, GPBP1L1, and YEATS4, respectively.
- LSI2 has a heavy right tail (99th percentile 0.00489; maximum 0.204), while its
  raw-depth Spearman correlation is only -0.083. This is a stronger concern than
  the modest depth association of LSI1 because a few candidates' effects are
  dominated by the LSI2 direction.
- Replicate-specific combined ranks correlate 0.951 (R1) and 0.874 (R2) with the
  primary ranking. HIC2 is rank 4 in R1 and rank 2 in R2.
- HIC2 remains rank 4 for thresholds 10, 15, 20, 25, and 30 cells per group. The
  primary threshold remains 20; CNOT2 is not added to the primary analysis.
- HIC2 ranges from rank 5 (ATAC-only) to rank 2 (RNA-only or 75% RNA weight), making
  the dependence on modality weighting explicit.
- Count-based NTC-SD scaling places HIC2 at rank 2 (rank correlation 0.720 with the
  primary ranking). This is a diagnostic, not a retrospectively selected score.

Guide and NTC composition tables preserve exact guide-call classes by condition and
replicate. The common-guide sensitivity retains only guide classes meeting the
configured minimum in every group. These outputs expose, rather than silently average
over, imbalances such as guide classes restricted to one treatment/replicate group.

All detailed outputs are under `results/tables/phase07_sensitivity/`; the compact
summary is versioned in `results/summary/final_sensitivity_summary.csv`.
