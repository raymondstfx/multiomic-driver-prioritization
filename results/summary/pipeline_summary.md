# Pipeline summary

Generated from analysis version `primary_v1` at 2026-09-23T12:01:59.343056+00:00.
Git commit: `1951acfadf6347a7b51e004c84b27ee80bfa937d`. Analysis-config SHA-256:
`229b37e1f6895c0c1a446d124a7f3aed2572742d5428f2ddd1e8d02968fc6b57`.
Data-config SHA-256: `891e26e5ad9eee55660fb8c9a8bb2460858efa6adb847f8d304e27e0a48d3a52`.
Runtime: Python `3.13.0` on `Windows-11-10.0.26200-SP0`. Expected
guide-call archive SHA-256: `2e5992e57d9ec0c1575ab5316e50ab311a0724a846f4f65035c4033e8d1fd1bf`.

Dataset: GSE288996, 53,898 paired K562 cells; DMSO and Dasatinib; biological
replicates R1 and R2. Twelve perturbations pass the configured 20-cell minimum in
all four groups; missing groups remain NA and CNOT2 remains below threshold.

## Primary result

The primary Phase 5 ranking is unchanged by Phase 7. The top five candidates are
YEATS4, GPBP1L1, ZBED6, HIC2, KMT2B. HIC2 is RNA rank
2, ATAC rank 5, and combined rank
4.

For candidate $i$, replicate $r$, and modality $m$, the primary interaction is

$$
\mathbf E_{i,r}^{(m)}=(\bar{\mathbf x}_{i,r,DASA}^{(m)}-
\bar{\mathbf x}_{NTC,r,DASA}^{(m)})-(\bar{\mathbf x}_{i,r,DMSO}^{(m)}-
\bar{\mathbf x}_{NTC,r,DMSO}^{(m)}).
$$

The score is the Euclidean norm of the replicate-mean effect and is relative to
the current eligible candidate cohort. Its magnitude does not encode a biological
resistance direction and is not causal proof.

## Robustness

The primary estimator is the mean latent-space centroid. Alternative estimators are
reported as sensitivity analyses, not replacements. Their Spearman correlations with
the primary combined rank range from 0.755 to
1.000. HIC2 ranks from
2 to 4 across these estimators.
The detailed tables also cover leave-one-dimension-out, NTC-SD scaling, guide
composition, replicate-specific ranks, cell thresholds, and RNA/ATAC weights.

## ZFPM2 downstream validation

ZFPM2 is absent from the perturbation library and is not ranked. HIC2-versus-NTC
ZFPM2 responses are descriptive and metric-dependent:

- `mean_expression`: mean replicate DiD = -0.2692; direction agreement = True.
- `median_expression`: mean replicate DiD = -0.103; direction agreement = True.
- `fraction_expressing`: mean replicate DiD = 0.1091; direction agreement = True.
- `mean_positive_expression`: mean replicate DiD = -0.3422; direction agreement = True.
- `pseudobulk_log1p_cpm`: mean replicate DiD = -0.1012; direction agreement = True.

These results do not establish direct HIC2-to-ZFPM2 regulation or causal mechanism.

## Navigation

- `final_candidate_ranking.csv`: canonical primary result.
- `final_candidate_evaluation.csv`: Phase 6 interpretation and diagnostics.
- `final_validation_summary.csv`: multi-metric ZFPM2 validation.
- `final_sensitivity_summary.csv`: compact robustness summary.
- `../manifest.csv`: provenance, dimensions, hashes, and roles for all artifacts.
