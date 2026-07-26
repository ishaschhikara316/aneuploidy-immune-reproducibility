# Cross-cohort reproducibility of aneuploidy-immune associations

Systematic audit of **6 genomic-instability exposures x 10 immune outcomes = 60 associations**, each estimated in three settings.

| setting | donors | RNA pipeline | n | significant (FDR<0.05) |
|---|---|---|--:|--:|
| S1 | TCGA | TCGA-Xena | 716 | 29/60 |
| S2 | TCGA (same as S1) | PCAWG | 716 | 27/60 |
| S3 | non-TCGA | PCAWG | 423 | 0/60 |

## Replication of the associations that are significant in S1 (FDR<0.05)

| comparison | n | sign agreement | replication rate | beta correlation | median effect ratio |
|---|--:|--:|--:|--:|--:|
| PIPELINE  (same donors: Xena -> PCAWG RNA) | 29 | 100% | 100% | 0.97 | 1.08 |
| COHORT    (same pipeline: TCGA -> non-TCGA) | 27 | 100% | 19% | 0.43 | 0.50 |
| BOTH      (Xena/TCGA -> PCAWG/non-TCGA) | 29 | 97% | 17% | 0.46 | 0.58 |

*replication rate = same sign AND p<0.05 in the replication setting.*
*median effect ratio = |beta_replication| / |beta_discovery|; <1 means effects shrink.*
