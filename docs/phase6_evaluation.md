# Phase 6 evaluation and biological validation

## Scope

Phase 6 evaluates the fixed Phase 5 ranking without changing its scores or
ranks. It separates three questions: computational agreement between
modalities, post-hoc recovery of HIC2, and downstream ZFPM2 expression. HIC2
was not used to tune the ranking, and ZFPM2 remains outside the perturbation
ranking because it is absent from the guide library.

## Ranking evaluation

All input checks passed: 12 unique primary candidates were present, HIC2 was
present, ZFPM2 was absent, numeric fields were finite, and each rank column was
a complete permutation of 1 through 12. The fixed top five are YEATS4,
GPBP1L1, ZBED6, HIC2, and KMT2B.

The descriptive Spearman correlations were 0.601 for RNA versus ATAC ranks,
0.755 for RNA versus multi-omic ranks, and 0.958 for ATAC versus multi-omic
ranks. With only 12 candidates these values describe this screen and should not
be generalised. YEATS4 and GPBP1L1 have positive standardized evidence in both
modalities; ZBED6 is ATAC-dominant; HIC2 and KMT2B are RNA-dominant under the
predefined sign rule.

HIC2 is RNA rank 2, ATAC rank 5, and multi-omic rank 4, with RNA and ATAC
replicate cosines of 0.906 and 0.574. Thus HIC2 is strongly prioritised by
RNA-based treatment-specific effects and receives additional ATAC evidence,
but multi-omic integration does not improve its position relative to RNA alone.
The framework nevertheless places this experimentally supported pooled
perturbation among the top candidates.

Replicate consistency remains separate from the ranking. Both modalities meet
the configured moderate/high category for YEATS4, GPBP1L1, HIC2, KMT2B, and
CHD2. ZBED6 has weak positive RNA agreement but moderate/high ATAC agreement.
ADNP and SIN3A have directionally inconsistent ATAC effects, and MNT has a
directionally inconsistent RNA effect.

## Literature annotation

The focused evidence table uses source-linked qualitative categories only.
The authoritative CAT-ATAC BioProject record directly connects HIC2 with a
Dasatinib-resistance regulatory network in the source K562 experiment. YEATS4,
ZBED6, and KMT2B have primary mechanistic evidence for transcriptional or
chromatin roles, but the cited experiments do not directly validate Dasatinib
response in K562/CML. No clear direct primary-literature support for GPBP1L1 in
Dasatinib response, K562, or CML was found; that absence remains explicit.

## ZFPM2 downstream validation

ZFPM2 was found once in processed RNA as Ensembl feature `ENSG00000169946`.
The [authoritative source-dataset record](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1220572)
reports independent ZFPM2 loss-of-function evidence for Dasatinib resistance;
the present analysis asks only whether the pooled HIC2 perturbation has a
replicate-consistent downstream ZFPM2 expression signal.
The analysis included HIC2 and NTC singlets in all four condition-by-replicate
groups: HIC2 group sizes were 183, 230, 325, and 320; NTC group sizes were 757,
994, 1,316, and 1,476 for DMSO R1, DMSO R2, Dasatinib R1, and Dasatinib R2,
respectively.

HIC2-minus-NTC mean log-expression contrasts were -0.573 and -0.558 in DMSO
R1/R2, and -0.833 and -0.836 in Dasatinib R1/R2. The treatment-specific DiD was
-0.260 in R1 and -0.278 in R2, with a mean of -0.269 and matching directions.
This is reproducible descriptive evidence that the HIC2 perturbation is
associated with a more negative ZFPM2 contrast under Dasatinib than DMSO. It is
not proof of direct HIC2-to-ZFPM2 regulation, does not make individual cells
independent biological replicates, and is not a ZFPM2 perturbation score.

## Outputs

Phase 6 generates rank comparisons, correlations, reliability annotations,
HIC2 and top-five summaries, literature evidence, final candidate evaluation,
ZFPM2 group summaries, HIC2-minus-NTC contrasts, replicate-specific ZFPM2 DiD,
and two evaluation figures under `results/`. Generated results remain ignored
by Git; the code, tests, configuration, and documentation are versioned.

## Limitations and readiness

The study uses one K562 dataset, one drug/control comparison, two biological
replicates, 12 primary eligible perturbations, a 20-cell/group eligibility
threshold, latent-space effect magnitude, and equal modality weights. It does
not reconstruct a full GRN, infer selection coefficients, implement MPL, or
prove causal drivers. Multi-omic integration need not outperform the best
single modality, as HIC2 demonstrates.

Potential future work includes threshold sensitivity, stability-aware or
alternative modality weighting, feature-level RNA-ATAC linking, regulatory and
pathway analysis, new drugs/cell lines/screens, joint latent models, and
MPL-inspired selection modelling. None is part of Phase 6.
