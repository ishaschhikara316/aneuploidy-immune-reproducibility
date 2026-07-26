# Is the non-TCGA null a power artifact?

S2 (TCGA, PCAWG RNA) n=716 vs S3 (non-TCGA, PCAWG RNA) n=423.

## 1. Down-sampling test (the decisive one)

S2 subsampled to n=423 (200 draws), counting associations at FDR<0.05:

- **matched-n S2: median 18 significant [2–29]** (of 60)
- **actual S3 at the same n: 0 significant**
- (full-size S2: 27)

At identical sample size the TCGA cohort still yields ~18 significant associations while the non-TCGA cohort yields 0. **The non-TCGA null is therefore NOT explained by power.**

## 2. Formal difference in effect sizes (S2 vs S3)

- associations with significantly different betas (FDR<0.05): **0/60**
- median |beta| S2 = 0.080, S3 = 0.056

## 3. Effect-size attenuation

- median |beta_nonTCGA| / |beta_TCGA| = **0.58** [95% CI 0.35–0.74]

