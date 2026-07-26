# Does TCGA replicate itself?

Everything held constant (same pipeline, exposures, models, n=340 per group, 100 draws); only the SOURCE of the replication sample changes.

Replication rule (field standard): of associations at FDR<0.05 in discovery, the fraction with p<0.05 and the same sign in the replication set.

| arm | draws | median discovery hits | **median replication rate** [95%] |
|---|--:|--:|--:|
| SELF  (TCGA -> disjoint TCGA) | 95 | 14 | **62%** [24%–100%] |
| CROSS (TCGA -> non-TCGA) | 95 | 14 | **8%** [0%–38%] |
| SELF, prevalence-matched | 95 | 14 | **50%** [3%–100%] |

- TCGA replicating **itself**: **62%**
- TCGA replicating in **non-TCGA**: **8%**
- TCGA replicating itself at non-TCGA **exposure prevalence**: **50%**
- fraction of SELF draws that look at least as 'failed' as the typical CROSS draw: **1%**

Self-replication is markedly higher than cross-cohort replication, so the cohort difference is not purely a metric artefact.

