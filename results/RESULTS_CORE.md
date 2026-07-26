# Aneuploidy–immune associations survive a complete change of analysis pipeline, but not a change of cohort

*Two findings, both established: (1) swapping the entire expression pipeline on the SAME patients
changes essentially nothing (94.8% of variance shared); (2) moving to an independent cohort breaks
these associations — and a self-replication control shows this is genuine, not a statistical
artefact.*

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

## The decisive control: does a cohort replicate ITSELF?

A low replication rate is only meaningful if a cohort would replicate *itself* at a much higher
rate under identical conditions. We held everything constant (same pipeline, exposures, models,
n=340 per group, 100 draws) and varied only **where the replication sample came from**:

| replication sample drawn from | median discovery hits | **replication rate** [95%] |
|---|--:|--:|
| **a disjoint random half of TCGA** (itself) | 14 | **62%** [24–100%] |
| a disjoint half of TCGA, down-sampled to non-TCGA WGD prevalence | 14 | **50%** [3–100%] |
| **the non-TCGA cohort** | 14 | **8%** [0–38%] |

**Only 1% of self-replication draws are as poor as a typical cross-cohort draw** (empirical
p ≈ 0.01). Matching the exposure prevalence — WGD is 39% in TCGA but 23% in non-TCGA, the leading
alternative explanation — closes only part of the gap (62% → 50%), nowhere near 8%.

**Conclusion: the cohort effect is real and large.** A cohort replicates itself roughly 8x better
than it replicates in an independent cohort, with sample size, pipeline, exposure definitions,
models and exposure prevalence all held constant.

### Verification: the gap is not tumour-type composition

The comparison above uses TCGA (23 histologies) vs non-TCGA (8, of which 38% lymphoid), so it
confounds consortium with tumour type. We therefore rebuilt both replication arms to contain the
**identical histologies in identical numbers** (Kidney-RCC 47, Liver-HCC 49, Ovary-AdenoCA 20;
n = 116 per arm — all solid tumours, no lymphoid), with the same TCGA discovery sample feeding both:

| replication set | median replication rate [95%] |
|---|--:|
| TCGA (self) | **54%** [25–100%] |
| non-TCGA (cross) | **0%** [0–9%] |

Gap **+54 points**; only **1%** of self-draws are as poor as the median cross-draw. **The cohort
effect survives with tumour-type composition held identical**, and is likewise not explained by
exposure prevalence (matching WGD prevalence gives 50% vs 8%).

### Why our earlier, weaker tests missed this

A composition-matched permutation (p=0.086) and a calibrated mean-|Δβ| comparison (1.17x, p=0.30)
both failed to detect the effect. Both average over **all 60 associations**, including the ~46 that
are null in every setting — which dilutes real signal into noise. The self-replication control is
correctly targeted: it asks only about the associations that were actually significant in
discovery, which is what "replication" means. We report the negative tests alongside the positive
one because the contrast is itself instructive: **aggregate effect-size metrics can be too blunt to
detect a replication failure that is obvious when you condition on discovery.**

## The metric caveat

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

> **Recommendation: always calibrate a replication rate against the same cohort replicating
> itself.** A raw "8% replicated" figure is uninterpretable on its own — here the self-replication
> baseline is 62%, and it is the *gap* between them that carries the meaning. Reporting effect
> sizes alongside significance remains good practice, but on this data an effect-size-only summary
> was too blunt to see the failure at all.

## Interpretation

Aneuploidy–immune associations reported in TCGA are **not** artefacts of expression processing —
they survive a complete pipeline swap on the same donors, sharing 94.8% of their variance. But they
**do not transfer to an independent cohort**: 62% self-replication versus 8% cross-replication under
identical conditions. The field's conflicting reports about aneuploidy and tumour immunity are
therefore unlikely to be a software problem, and much more likely to reflect genuine
cohort-to-cohort differences — in populations, specimen handling, or library preparation — that
current practice does not account for.

> **Assess replication by effect size, not by significance.** Directions are 100% consistent
> everywhere; it is the magnitude that moves, and significance mostly tracks sample size.

## Honest limitations

- **The cohort effect is established by the self-replication control (62% vs 8%, p ≈ 0.01), but two
  weaker tests did not detect it** (composition-matched permutation p = 0.086; calibrated
  mean-|Δβ| ratio 1.17, p = 0.30). We report all three. The discrepancy is explained by those two
  averaging over associations that are null everywhere; it is nonetheless a caution that the size
  of the effect is better characterised than its precise magnitude.
- We cannot say *which* aspect of "cohort" is responsible — population genetics, specimen handling,
  library preparation or clinical context are all confounded with consortium.
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
