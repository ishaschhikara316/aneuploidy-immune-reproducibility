# Calibrated test: cohort swap vs re-sampling the same cohort

Both arms use matched group size **n = 330** on PCAWG RNA (pipeline constant), 200 draws each, so sampling noise is matched by construction and no histology restriction is needed.

| arm | mean \|Δβ\| across 60 associations | 1 − r |
|---|--:|--:|
| NULL — two disjoint random halves of TCGA | 0.0654 [0.0432, 0.1050] | 0.465 |
| OBSERVED — TCGA vs non-TCGA | **0.0768** [0.0549, 0.1049] | **0.719** |

- **Cohort swaps are 1.17x more disruptive** than re-sampling the same cohort (95% CI 1.10–1.22).
- permutation-style p (null >= observed): **p = 0.2985** (mean \|Δβ\|), **p = 0.0896** (1−r).

**Not significant:** a cohort swap is not measurably more disruptive than re-sampling the same cohort at this sample size.

