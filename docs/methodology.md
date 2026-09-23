# Methodology

The first implementation ranks candidate perturbations from treatment-specific,
background-corrected effects in RNA PCA and ATAC LSI space. Latent centroids are computed
by perturbation, condition, and biological replicate. For each perturbation, a
Difference-in-Differences contrast subtracts its DMSO effect from its Dasatinib
effect, each relative to non-targeting controls.

RNA and ATAC effect vectors are scored by Euclidean magnitude. Their scores are
standardised separately and averaged. This intentionally simple integration is
auditable and does not use known biology as a label or score input. HIC2 is a
genuine library perturbation used only for post-ranking evaluation. ZFPM2 is
not in the pooled library and is reserved for downstream/external biological
validation rather than perturbation scoring.

## Phase 3 preprocessing

RNA cells are conservatively filtered at 200 detected genes and genes at three
cells. Counts are preserved in `layers["counts"]`; `X` is library-size
normalised to 10,000 counts per cell and log1p transformed. Three thousand HVGs
are selected with Scanpy's Seurat-compatible method and 30 PCs are calculated
from the HVGs while all retained genes remain in the object.

The source ATAC matrices contain independently called, nearly sample-specific
peak coordinates. Before feature filtering, overlapping existing calls are
collapsed into non-overlapping consensus intervals so every sample uses the
same genomic feature space. This is coordinate harmonisation, not new peak
calling, and does not use fragment files. Consensus peaks are retained when
detected in at least three cells. Sparse TF-IDF is defined as
`TF_ij * log(1 + N / (1 + DF_j))`, where TF is the within-cell count fraction,
N is the number of cells, and DF is the peak's non-zero cell count. Truncated
SVD produces 30 LSI dimensions. LSI1 is retained; any later decision to exclude
it must be justified during downstream modelling.

## Phase 3.5 eligibility gate

The primary Phase 4 set contains targeting singlets only. A perturbation must
occur in all four condition × replicate groups and contain at least 20 cells in
each group. The threshold and all-group requirement are configured under
`phase4` in `configs/analysis.yaml`. Missing groups are represented as missing
and mean that the four-group effect is not estimable; they are not zero-valued
effects. NTC cells form the shared reference and are never ranked candidates.

RNA mitochondrial fractions are reported overall and by experimental group.
No mt% filter is applied because the observed distribution is continuous and
the condition-associated shift should not be erased without stronger evidence.

## Phase 4 replicate-aware effects

Biological replicates are never pooled before effect estimation. For candidate
`i` and replicate `r`, the latent-space vector is:

$$
E_{i,r}=(\bar{x}_{\mathrm{DASA},i,r}-\bar{x}_{\mathrm{DASA},NTC,r})
       -(\bar{x}_{\mathrm{DMSO},i,r}-\bar{x}_{\mathrm{DMSO},NTC,r}).
$$

This is calculated independently in the 30-dimensional RNA PCA space and the
30-dimensional ATAC LSI space. The reported mean effect is `(E_i,R1 + E_i,R2) /
2`. Replicate consistency is the cosine similarity between the two replicate
effect vectors. If either vector has zero magnitude, consistency is missing
rather than being imputed as zero.

The effect stage accepts only candidates marked eligible by Phase 3.5. It
validates that truly absent groups remain NA with `missing_*` reasons and that
present but undersized groups retain `below_min_*(count<threshold)` reasons.
Dedicated exceptions distinguish a missing group from an insufficient group.
No multi-omic ranking is calculated during this stage.

## Phase 5 multi-omic ranking

The primary modality scores are the existing Phase 4 `mean_effect_norm`
values—not the average of R1 and R2 norms. RNA and ATAC scores are independently
standardized over the complete 12-candidate primary population using z-scores
with `ddof=0`. The multi-omic score is the unweighted mean:

$$
S_i=\frac{z(S_i^{RNA})+z(S_i^{ATAC})}{2}.
$$

Ranks are descending and ties use deterministic minimum rank. Zero-variance,
non-finite, duplicate-candidate, and cross-modality candidate-set errors fail
explicitly. Replicate cosine similarity is retained in the output but is not
used as a score multiplier or exclusion rule. Known biology and validation
labels are not ranking inputs.

This method provides treatment-specific perturbation evidence and cross-modal
support. It is not MPL inference and is not, by itself, causal proof.

## Phase 6 evaluation and external validation

Phase 6 consumes the fixed `candidate_ranking.csv`; it does not recompute or
modify Phase 5 scores. The input gate requires 12 unique primary candidates,
finite primary scores, valid (possibly tied) integer ranks, HIC2 presence, and ZFPM2
absence. Rank shifts are `single_modality_rank - multiomic_rank`, so a positive
value means improvement under integration. Spearman correlations are reported
descriptively because the candidate set is small.

Cross-modal profiles are derived only from the signs of the two modality
z-scores. Replicate cosines are classified using configurable thresholds
(`>=0.5` moderate/high, `0` to `<0.5` weak positive, `<0` directional
inconsistency) and remain annotations rather than score multipliers. HIC2 is
looked up only after the ranking is fixed. Qualitative literature evidence is
source-linked and never converted to a numerical score.

ZFPM2 is present in processed RNA as `ENSG00000169946` but absent from the guide
library. Its log-normalised expression is summarised for HIC2 and pooled NTC
singlet cells separately within DMSO/Dasatinib and R1/R2. For replicate `r`, the
targeted downstream contrast is:

$$
\mathrm{ZFPM2\ DiD}_r=(y_{\mathrm{DASA},HIC2,r}-y_{\mathrm{DASA},NTC,r})
 -(y_{\mathrm{DMSO},HIC2,r}-y_{\mathrm{DMSO},NTC,r}).
$$

Only after the two replicate effects are estimated are their mean and direction
agreement reported. Individual cells are not treated as independent biological
replicates for inferential claims. This analysis tests downstream consistency;
it does not assign ZFPM2 a perturbation rank or establish direct regulation.

The final implemented sequence is: verified RNA/ATAC/guide alignment; RNA and
ATAC preprocessing; eligibility filtering; replicate-specific latent centroids;
replicate-specific background-corrected DiD; replicate-mean effect vectors;
modality-specific effect norms; within-modality z-scoring; equal-weight
integration; and post-hoc experimental, literature, and downstream validation.

## Phase 7 robustness and reporting

The fixed `primary_v1` ranking remains the main analysis. Phase 7 separately
tests winsorised, trimmed, and median latent centroids; leave-one-dimension-out
scores; NTC-SD scaling; exact-guide composition; replicate-specific rankings;
cell-count thresholds; and RNA/ATAC weights. These analyses diagnose dependence
on heavy-tailed latent dimensions and design imbalance; they do not select a
more favourable ranking after observing HIC2. Phase 8 creates a compact summary
and a hash-based artifact manifest.
