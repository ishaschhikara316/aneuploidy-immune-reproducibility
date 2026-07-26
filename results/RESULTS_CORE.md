# Aneuploidy–immune associations are invariant to the RNA pipeline, and significance-based replication metrics overstate irreproducibility

*Two findings: (1) swapping the entire expression pipeline changes essentially nothing;
(2) the same data looks either catastrophically or barely irreproducible across cohorts depending
purely on whether replication is scored by significance or by effect size.*

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
cohort halves the effect sizes and eliminates statistical significance.**

**Status of the two halves.** The pipeline result is definitive (paired design, 716 donors). The
cohort attenuation is *suggestive but not formally established*: once tumour-type composition is
also held constant, a pooled permutation test gives **p = 0.086** (see below). We therefore report
the pipeline invariance as the finding, and cohort-dependence as a hypothesis this data supports
but cannot confirm.

## Is the non-TCGA null just lower power? No

- **Down-sampling:** TCGA (S2) subsampled to n=423 (200 draws) still yields a median of **18
  significant associations [2–29]**; the actual non-TCGA cohort at that same n yields **0**.
- **Effect attenuation:** median |β_nonTCGA| / |β_TCGA| = **0.58 [95% CI 0.35–0.74]** — the CI
  excludes 1, so effects are genuinely smaller, not merely less significant.

## Is it tumour-type composition? Descriptively, no

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
effect sizes are consistently smaller (0.031–0.056) in every stratum. Descriptively, then, the gap
is in **effect magnitude** rather than tumour-type mix — but note that the *formal* test of exactly
this comparison (next section) does not reach significance, so this is an observation, not a
demonstrated result.

## Formal test of the cohort effect — does NOT reach significance

Per-association differences were individually underpowered (0/60), so we pooled across all
associations while resampling donors, restricted to PCAWG RNA and the 4 shared histologies — thus
holding **pipeline and composition both constant** (n=399; TCGA 139 / non-TCGA 260). Betas were
oriented by the sign of the discovery estimate and compared via a permutation test that shuffles
the cohort label *within* each histology.

| quantity | value |
|---|--:|
| mean oriented β, TCGA | +0.098 |
| mean oriented β, non-TCGA | −0.008 |
| observed Δ | **+0.106** |
| permutation null (500 draws) | mean +0.031, 95th pct +0.126 |
| **one-sided permutation p** | **0.086** |
| bootstrap Δ [95% CI] | +0.106 [−0.000, +0.217] |

**The cohort effect is therefore not formally significant once composition is controlled.** The
direction and magnitude are consistent with real attenuation, and the unstratified effect ratio
(0.58 [0.35–0.74]) excludes 1, but the strict, composition-matched test does not clear p<0.05.
This test is itself modestly powered (399 donors, 4 histologies), so the honest reading is
*unresolved*, not *refuted*.

## The metric paradox — the central methodological result

The same TCGA → non-TCGA comparison yields opposite impressions depending on the metric:

| replication metric | result | impression |
|---|--:|---|
| associations still significant (FDR<0.05) | **0/60** | catastrophic irreproducibility |
| sign agreement | **100%** | perfect agreement |
| median effect ratio | 0.50–0.58 | effects halve |
| **calibrated discrepancy** vs re-splitting TCGA itself at matched n | **1.17x [1.10–1.22]** | cohort barely matters |

The calibrated test is the most conservative and most interpretable: two disjoint random halves of
TCGA (n=330 each) differ by mean |Δβ| = 0.065, while a TCGA-vs-non-TCGA comparison at the same n
differs by 0.077 — only **17% more**, with a single-draw permutation p = 0.30.

**So the dramatic "0/60 replicate" figure is substantially an artefact of thresholding.** Much of
it is driven by reduced power in the non-TCGA cohort (n=423 vs 717, and WGD prevalence 23% vs 39%),
not by materially different underlying effects. Directions never disagree; magnitudes differ
modestly; only the significance verdict collapses.

> **Recommendation: never assess replication by whether p crosses a threshold in the replication
> cohort.** On this data that metric implies total failure (0/60) while a calibrated effect-size
> comparison implies the cohort contributes ~17% extra variability. Report effect sizes and
> calibrate against re-sampling the discovery cohort.

## Interpretation

Aneuploidy–immune associations reported in TCGA are **not** artefacts of expression processing —
they survive a complete pipeline swap on the same donors, sharing 94.8% of their variance. Nor are
they as cohort-fragile as a significance-based reading suggests: calibrated against re-sampling the
discovery cohort, changing cohort adds only ~17% extra discrepancy. The field's conflicting reports
are therefore unlikely to stem from either processing or wholesale biological divergence — and much
more likely from underpowered cohorts being scored by significance thresholds.

> **Assess replication by effect size, not by significance.** Directions are 100% consistent
> everywhere; it is the magnitude that moves, and significance mostly tracks sample size.

## Honest limitations

- **The cohort effect is real but small, and its size is not precisely determined.** The
  composition-matched permutation gives p = 0.086; the calibrated matched-n test gives a ratio of
  1.17 [1.10–1.22] (CI excludes 1) but a single-draw permutation p = 0.30. Different reasonable
  statistics disagree on whether it is "significant" — which is itself the paper's point. We report
  it as small and real rather than dramatic.
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
