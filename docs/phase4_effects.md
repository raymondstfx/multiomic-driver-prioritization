# Phase 4 replicate-aware effects

Phase 4 uses the 12 primary candidates approved by the Phase 3.5 gate and keeps
the configured minimum at 20 cells per condition × replicate group. CNOT2 is
not included merely to enlarge the candidate set. The five candidates with
genuinely missing groups also remain excluded from the primary analysis.

For RNA PCA and ATAC LSI independently, cell-level latent coordinates are
averaged within candidate/NTC × condition × replicate groups. Difference-in-
Differences is then computed separately for each replicate:

```text
E_i,R1 = (DASA_i,R1 - DASA_NTC,R1) - (DMSO_i,R1 - DMSO_NTC,R1)
E_i,R2 = (DASA_i,R2 - DASA_NTC,R2) - (DMSO_i,R2 - DMSO_NTC,R2)
```

The two vectors are never pooled before these contrasts. Each modality writes
24 replicate rows (12 candidates × 2 replicates) and 12 summary rows. Summary
fields contain the R1 and R2 vector norms, the norm of their elementwise mean,
and cosine similarity between R1 and R2 as a separate replicate-consistency
metric.

## Missing versus insufficient groups

Before calculating effects, the eligibility table is validated. A missing
group must have an NA count, NA minimum-group count, and a `missing_*` reason.
A present group below 20 retains its observed non-zero count and a
`below_min_*(count<20)` reason. Turning an NA count into zero fails validation.
Runtime checks also raise different exception types for absent and undersized
groups. No absent group is represented by a zero vector or zero effect.

## Outputs

Generated, Git-ignored tables are:

```text
results/tables/rna_replicate_effects.csv
results/tables/rna_effect_summary.csv
results/tables/atac_replicate_effects.csv
results/tables/atac_effect_summary.csv
results/tables/phase4_effect_manifest.csv
```

The manifest records 12 candidates, 24 replicate-specific rows, 30 dimensions
per modality, and `replicates_merged_before_effect_estimation=False`. Candidate
ranking and cross-modal score integration remain separate later stages.
