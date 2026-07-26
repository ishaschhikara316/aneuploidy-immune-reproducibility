# Is the cohort effect just tumour-type composition?

Shared histologies between TCGA and non-TCGA (PCAWG RNA): **['Kidney-RCC', 'Liver-HCC', 'Lymph-BNHL', 'Ovary-AdenoCA']**

## Associations significant (FDR<0.05) per setting

| setting | donors | associations | significant | median \|beta\| |
|---|--:|--:|--:|--:|
| TCGA (all histologies) | 717 | 60 | 29 | 0.082 |
| TCGA (shared histologies only) | 139 | 60 | 1 | 0.084 |
| nonTCGA (all) | 423 | 60 | 0 | 0.056 |
| nonTCGA (shared histologies only) | 260 | 60 | 0 | 0.045 |
| nonTCGA (solid only) | 262 | 60 | 0 | 0.031 |

**Reading:** if the non-TCGA null were driven by its different histology mix, restricting to shared histologies (or removing lymphoid tumours) should recover signal.

## Canonical association (WGD -> cytolytic) per shared histology

| histology | cohort | n | WGD+ | beta | p |
|---|---|--:|--:|--:|--:|
| Kidney-RCC | TCGA | 62 | 9 | -0.112 | 0.754 |
| Kidney-RCC | nonTCGA | 47 | 7 | +0.076 | 0.837 |
| Liver-HCC | TCGA | 50 | 20 | -0.423 | 0.0787 |
| Liver-HCC | nonTCGA | 48 | 11 | -0.232 | 0.458 |
| Lymph-BNHL | nonTCGA | 96 | 9 | -0.130 | 0.644 |
| Ovary-AdenoCA | nonTCGA | 68 | 43 | -0.467 | 0.0584 |
