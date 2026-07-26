# Aneuploidy–immune associations are robust to RNA pipeline but fragile across cohorts

*Status: analysis complete, manuscript in preparation. Public data only, CPU-only.*

## Motivation

The literature on genomic instability and tumour immunity contains directly contradictory
reports: whole-genome doubling (WGD) is said to drive immune evasion by silencing antigen
presentation, yet other work reports that tetraploidy *increases* immunogenicity. Almost all of
these findings come from TCGA. We asked a simple question with a clean design: **when an
aneuploidy–immune association fails to replicate, is it the sequencing/processing pipeline, or the
cohort?**

## Design — the two-arm decomposition

Exposures (WGD, ploidy, chromothripsis high/any, percent-genome-altered, number of chromothripsis
events) and tumour purity come from **one shared donor table** (Cortés-Ciriano et al., *Nat Genet*
2020, ShatterSeek/PCAWG), so exposure definitions are **identical in every comparison**. Only one
thing changes at a time:

| setting | donors | RNA pipeline | n |
|---|---|---|--:|
| **S1** | TCGA | TCGA-Xena (RSEM) | 716 |
| **S2** | TCGA — *the same 716 donors* | PCAWG (tophat/STAR FPKM-UQ) | 716 |
| **S3** | non-TCGA (ICGC) | PCAWG | 423 |

- **S1 → S2 isolates the PIPELINE effect** (donors held constant — a genuinely paired comparison)
- **S2 → S3 isolates the COHORT effect** (pipeline held constant)

6 exposures × 10 immune outcomes = **60 associations**, each estimated in all three settings
(model: `outcome ~ exposure + purity + proliferation`; panel scores z-standardised within
histology and within RNA source; FDR-corrected within setting).

## Headline result

| comparison | associations | sign agreement | replication rate | β correlation | median effect ratio |
|---|--:|--:|--:|--:|--:|
| **PIPELINE** (same donors, Xena → PCAWG) | 29 | **100%** | **100%** | **0.97** | 1.08 |
| **COHORT** (same pipeline, TCGA → non-TCGA) | 27 | 100% | **19%** | 0.43 | **0.50** |

Significant associations (FDR<0.05): **S1 = 29/60, S2 = 27/60, S3 = 0/60**.

**Changing the entire RNA pipeline — different aligner, quantifier and normalisation — barely
perturbs these associations (β correlation 0.97, every association replicates). Changing the
cohort halves the effect sizes and eliminates statistical significance entirely.**

## Is the non-TCGA null just lower power? No

- **Down-sampling:** TCGA (S2) subsampled to n=423 (200 draws) still yields a median of **18
  significant associations [2–29]**; the actual non-TCGA cohort at that same n yields **0**.
- **Effect attenuation:** median |β_nonTCGA| / |β_TCGA| = **0.58 [95% CI 0.35–0.74]** — the CI
  excludes 1, so effects are genuinely smaller, not merely less significant.

## Is it tumour-type composition? No

The non-TCGA cohort is 38% lymphoid and has a different histology mix, so we matched on histology:

| setting | donors | significant | median \|β\| |
|---|--:|--:|--:|
| TCGA (all histologies) | 717 | 29/60 | **0.082** |
| TCGA (shared histologies only) | 139 | 1/60 | **0.084** |
| non-TCGA (all) | 423 | 0/60 | 0.056 |
| non-TCGA (shared histologies only) | 260 | 0/60 | **0.045** |
| non-TCGA (solid only) | 262 | 0/60 | **0.031** |

Restricting TCGA to the shared histologies **leaves its effect sizes unchanged** (0.082 → 0.084);
only significance falls, and that is pure sample-size loss (717 → 139). In the non-TCGA cohort the
effect sizes are consistently smaller (0.031–0.056) in every stratum. So the difference is in the
**effect magnitude**, not the tumour types.

## Interpretation

Aneuploidy–immune associations reported in TCGA are **not** artefacts of expression processing —
they survive a complete pipeline swap on the same donors. But their magnitude is cohort-dependent,
shrinking roughly by half in independent ICGC cohorts. This provides a concrete, quantitative
explanation for why the field's reports conflict, and a practical recommendation:

> **Assess replication by effect size, not by significance.** Directions are 100% consistent
> everywhere; it is the magnitude that moves, and significance mostly tracks sample size.

## Honest limitations

- Per-association β differences between S2 and S3 are **not individually significant (0/60)** —
  the attenuation is a consistent *aggregate* shift rather than dramatic per-association reversals.
- Only **4 histologies are shared** between the cohorts, and per-type n is small.
- WGD prevalence differs (TCGA 39% vs non-TCGA 23%), which contributes to both power and possibly
  the effect estimate.
- Cross-sectional bulk expression; the immune panels are marker-score proxies, not cell counts.
- We cannot fully distinguish "different populations/biology" from residual differences in
  specimen handling and library preparation between consortia — only that it is not the
  *computational* pipeline.

## Files

`results/associations.tsv` (all 180 fits) · `results/reproducibility_metrics.md` ·
`results/power_check.md` · `results/histology_matched.md` · `results/beta_difference.tsv` ·
`figures/fig1_pipeline_vs_cohort.png` (key figure) · `fig2_replication_metrics.png` ·
`fig3_effect_sizes.png`
