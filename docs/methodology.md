# Methodology

The first implementation ranks candidate perturbations from treatment-specific,
background-corrected effects in RNA PCA and ATAC LSI space. Cells are aggregated
by perturbation, condition, and biological replicate. For each perturbation, a
Difference-in-Differences contrast subtracts its DMSO effect from its Dasatinib
effect, each relative to non-targeting controls.

RNA and ATAC effect vectors are scored by Euclidean magnitude. Their scores are
standardised separately and averaged. This intentionally simple integration is
auditable and does not use HIC2 or ZFPM2 as labels. Those genes are used only for
post-hoc evaluation.

This method provides treatment-specific perturbation evidence and cross-modal
support. It is not MPL inference and is not, by itself, causal proof.

