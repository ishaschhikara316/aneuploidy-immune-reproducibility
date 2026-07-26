# How robust is pipeline invariance?

Paired arm only: the same 716 TCGA donors scored through TCGA-Xena vs PCAWG RNA.

## 1. Does it hold for WEAK effects too?

| effect-size tertile | n | median \|β\| | β correlation | sign agreement |
|---|--:|--:|--:|--:|
| weak | 20 | 0.033 | **0.746** | 95% |
| medium | 20 | 0.076 | **0.534** | 100% |
| strong | 20 | 0.154 | **0.953** | 100% |

Sign agreement stays 95-100% in every stratum, so the headline is not carried by a few large
effects. **Caveat:** within-tertile correlations are mechanically attenuated by range restriction
(each tertile spans a narrow band of effect sizes), so they are not directly comparable to the
overall r = 0.97 — sign agreement is the fairer within-stratum metric.

## 2. Does it hold within individual tumour types?

| histology | donors | associations | β correlation | sign agreement |
|---|--:|--:|--:|--:|
| ColoRect-AdenoCA | 51 | 60 | 0.966 | 93% |
| Head-SCC | 42 | 60 | 0.955 | 92% |
| Breast-AdenoCA | 83 | 60 | 0.936 | 90% |
| Skin-Melanoma | 36 | 60 | 0.921 | 80% |
| Lung-SCC | 47 | 60 | 0.921 | 87% |
| Liver-HCC | 50 | 60 | 0.890 | 88% |
| Kidney-RCC | 62 | 60 | 0.777 | 85% |
| Lung-AdenoCA | 36 | 60 | 0.768 | 82% |
| Uterus-AdenoCA | 40 | 60 | 0.759 | 73% |
| Kidney-ChRCC | 39 | 60 | 0.506 | 58% |

Median within-histology β correlation: **0.906** (10 tumour types with >=30 donors).

## 3. Practical agreement interval (Bland-Altman)

- mean difference (PCAWG − Xena): **-0.0037** (no systematic bias)
- 95% limits of agreement: **[-0.035, +0.028]** in SD units

Any single association estimated in one pipeline can be expected to fall within roughly ±0.03 SD units of its value in the other.

## 4. Variance decomposition

- fraction of effect-estimate variance SHARED between pipelines: **94.8%**
- pipeline-specific / residual: **5.2%**

