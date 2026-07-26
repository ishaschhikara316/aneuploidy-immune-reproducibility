# aneuploidy-immune-reproducibility

**Aneuploidy–immune associations are robust to the RNA pipeline but fragile across cohorts.**

A cross-cohort reproducibility audit of genomic-instability → tumour-immune associations, using a
two-arm design that separates *technical* from *cohort* causes of replication failure.

## The design

Exposures (WGD, ploidy, chromothripsis, PGA…) and purity come from **one shared donor table**
(Cortés-Ciriano et al. *Nat Genet* 2020), so they are defined identically everywhere. Only one
thing changes at a time:

| | donors | RNA pipeline | n |
|---|---|---|--:|
| **S1** | TCGA | TCGA-Xena | 716 |
| **S2** | TCGA — *same donors* | PCAWG | 716 |
| **S3** | non-TCGA (ICGC) | PCAWG | 423 |

`S1→S2` = **pipeline effect** (paired, same donors) · `S2→S3` = **cohort effect** (same pipeline)

## Headline

6 exposures × 10 immune outcomes = 60 associations, each fit in all three settings:

| comparison | sign agreement | replication rate | β correlation | median effect ratio |
|---|--:|--:|--:|--:|
| **pipeline** | 100% | **100%** | **0.97** | 1.08 |
| **cohort** | 100% | **19%** | 0.43 | **0.50** |

Significant (FDR<0.05): S1 = 29/60, S2 = 27/60, **S3 = 0/60**.

Swapping the entire RNA pipeline barely moves these associations. Swapping the cohort halves them.
Ruled out: **power** (down-sampled TCGA at n=423 still gives ~18 hits vs 0) and **tumour-type
composition** (TCGA effect sizes are unchanged when restricted to shared histologies).

**Practical recommendation: assess replication by effect size, not significance** — directions are
100% consistent; magnitude is what moves.

See `results/RESULTS_CORE.md` for full tables, the power check, and honest limitations.

## Reproduce

```bash
python src/01_build_matrix.py       # extract panels from both RNA sources -> analysis matrix
python src/02_audit.py             # 60 associations x 3 settings + reproducibility metrics
python src/03_power_check.py       # down-sampling + effect-attenuation tests
python src/04_histology_matched.py # is it tumour-type composition?
python src/05_figures.py
```

## Data (public, CPU-only)

Chromothripsis/ploidy/purity: Cortés-Ciriano 2020 supp. · TCGA RNA: PanCanAtlas (Xena) ·
PCAWG RNA: open ICGC-25K bucket (`tophat_star_fpkm_uq.v2_aliquot_gl.tsv.gz`) — no access
application required. `data/` holds symlinks to the already-downloaded copies.

*Manuscript in preparation.*
