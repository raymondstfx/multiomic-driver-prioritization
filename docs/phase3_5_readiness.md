# Phase 3.5 QC and Phase 4 readiness

Phase 3.5 was executed on the paired processed RNA and ATAC metadata for 53,898
cells. No Phase 3 filtering threshold was changed and no Phase 4 effects were
calculated.

## Mitochondrial QC

The overall RNA mitochondrial fraction had median 16.52%, 90th percentile
26.54%, 95th percentile 29.79%, and 99th percentile 36.40%. The range was
0–68.08%.

Group medians were:

| Condition | Replicate | Cells | Median mt% | 95th percentile mt% |
|---|---:|---:|---:|---:|
| DMSO | 1 | 12,621 | 19.00 | 32.18 |
| DMSO | 2 | 13,547 | 18.13 | 30.73 |
| Dasatinib | 1 | 11,920 | 14.55 | 25.72 |
| Dasatinib | 2 | 15,810 | 15.27 | 28.08 |

DMSO cells show a systematic upward shift relative to Dasatinib, while the two
replicates within each condition are broadly consistent. The overall histogram
has a continuous right tail rather than a clearly separated failure population.
Applying a fixed 10%, 15%, or 20% cutoff would remove condition-associated
signal unevenly. No additional mt% filter is justified for the initial
proof-of-concept; mitochondrial content remains a documented QC characteristic.

## Candidate eligibility

The primary rule requires a targeting singlet perturbation, all four condition
× replicate groups, and at least 20 cells in every group. Counts absent from a
group are stored as missing values, not zeros.

Twelve candidates pass:

```text
ADNP, CHD2, GPBP1L1, HIC2, KMT2B, MNT,
PIAS1, SIN3A, SLTM, TSC22D4, YEATS4, ZBED6
```

Six remain documented but are excluded from the primary set:

| Candidate | Reason |
|---|---|
| BRD2 | Missing Dasatinib R1 |
| CNOT2 | DMSO R1 has 18 cells, below 20 |
| PLAGL2 | Missing DMSO R1 and Dasatinib R2 |
| PQBP1 | Missing DMSO R2 |
| ZNF330 | Missing Dasatinib R1 |
| ZNF669 | Missing Dasatinib R1 |

## Controls and validation genes

Pooled NTC singlet counts are 757, 994, 1,316, and 1,476 across DMSO R1, DMSO
R2, Dasatinib R1, and Dasatinib R2 respectively. Every reference group exceeds
the configured minimum. HIC2 has 183, 230, 325, and 320 cells and is explicitly
eligible. ZFPM2 is absent from the candidate table and remains downstream or
external biological validation only.

## Gate result

The Phase 4 gate passes: NTC coverage is sufficient, 12 targeting
perturbations are eligible, HIC2 is eligible, and mitochondrial QC does not show
an unresolved failure requiring Phase 3 to be rerun. The project is ready for
replicate-aware latent-centroid and background-corrected effect estimation.
