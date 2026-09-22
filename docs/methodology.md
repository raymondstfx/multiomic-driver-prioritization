# Methodology

The first implementation ranks candidate perturbations from treatment-specific,
background-corrected effects in RNA PCA and ATAC LSI space. Cells are aggregated
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

```text
E_i,r = (Dasatinib_i,r - Dasatinib_NTC,r)
      - (DMSO_i,r       - DMSO_NTC,r)
```

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

```text
multiomic_score_i = (rna_z_i + atac_z_i) / 2
```

Ranks are descending and ties use deterministic minimum rank. Zero-variance,
non-finite, duplicate-candidate, and cross-modality candidate-set errors fail
explicitly. Replicate cosine similarity is retained in the output but is not
used as a score multiplier or exclusion rule. Known biology and validation
labels are not ranking inputs.

This method provides treatment-specific perturbation evidence and cross-modal
support. It is not MPL inference and is not, by itself, causal proof.
