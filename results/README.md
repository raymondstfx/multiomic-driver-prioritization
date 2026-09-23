# Results guide

Start with [`summary/pipeline_summary.md`](summary/pipeline_summary.md). The canonical
ranking is [`summary/final_candidate_ranking.csv`](summary/final_candidate_ranking.csv),
and [`manifest.csv`](manifest.csv) describes every generated artifact and records its
hash and provenance.

Generated tables and figures are grouped by analysis phase:

- `phase03_qc`: RNA/ATAC preprocessing QC.
- `phase035_readiness`: candidate eligibility and explicit exclusion reasons.
- `phase04_effects`: replicate-specific and replicate-averaged DiD effects.
- `phase05_ranking`: the fixed `primary_v1` ranking.
- `phase06_validation`: ranking evaluation, literature context, and ZFPM2 validation.
- `phase07_sensitivity`: robustness analyses that never overwrite the primary result.

Only the compact `summary/` and manifest are versioned. Full generated phase outputs
remain reproducible and are excluded from Git.
