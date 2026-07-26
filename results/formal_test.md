# Formal test of the cohort effect (pooled, composition-controlled)

Per-association tests were individually underpowered (0/60 significant). This pools information across all associations while resampling **donors**, so the correlation between associations is respected.

Restricted to histologies shared by both cohorts (**Kidney-RCC, Liver-HCC, Lymph-BNHL, Ovary-AdenoCA**) and to PCAWG RNA only, so **tumour-type composition and pipeline are both held constant**. n = 399 (TCGA 139 / non-TCGA 260).

Betas oriented by the sign of the TCGA-Xena discovery estimate; Delta = mean(oriented beta, TCGA) − mean(oriented beta, non-TCGA).

## Primary — stratified permutation test

- observed mean oriented beta: **TCGA +0.0983**, **non-TCGA -0.0081**
- **observed Delta = +0.1064**
- permutation null (500 draws, cohort label shuffled *within* histology): mean +0.0308, 95th percentile +0.1257, max +0.2352
- **one-sided permutation p = 0.0858**

**The cohort effect is NOT formally significant** once composition is held constant — the aggregate attenuation cannot be distinguished from chance by this test.

## Secondary — donor bootstrap

- Delta = +0.1064 [95% CI -0.0000, +0.2165]
- effect ratio (non-TCGA / TCGA) = -0.08 [95% CI -1.46, 0.98]

