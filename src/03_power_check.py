"""
03 — Is the non-TCGA null a POWER artifact or genuine effect attenuation?

S3 (non-TCGA) has n=423 vs n=716 in S1/S2, and its effects are ~half the size. Both reduce
significance, so the audit's headline is only trustworthy if the two are separated.

Three tests:
  1. DOWN-SAMPLING: repeatedly subsample S2 (TCGA/PCAWG) to n=|S3| and count how many
     associations stay significant. If matched-n S2 still yields many hits, S3's null is NOT power.
  2. FORMAL beta-difference test between S2 and S3 (z on the difference of coefficients).
  3. Bootstrap CI on the median effect-size ratio.

Outputs results/power_check.md
"""
import warnings, numpy as np, pandas as pd, yaml, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from pathlib import Path
from scipy.stats import norm
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((ROOT / "config" / "panels.yaml").read_text())
EXPOSURES, OUTCOMES = CFG["exposures"], [p for p in CFG["panels"] if p != "proliferation"]
rng = np.random.default_rng(20260726)

d = pd.read_parquet(ROOT / "data" / "analysis_matrix.parquet")
xena = set(d[d.rna_source == "TCGA_Xena"].donor_unique_id)
pt = set(d[(d.rna_source == "PCAWG") & (d.cohort == "TCGA")].donor_unique_id)
PAIRED = xena & pt
S2 = d[(d.rna_source == "PCAWG") & d.donor_unique_id.isin(PAIRED)]
S3 = d[(d.rna_source == "PCAWG") & (d.cohort == "nonTCGA")]
N3 = len(S3)


def expo(frame, spec):
    x = pd.to_numeric(frame[spec["col"]], errors="coerce")
    if spec["kind"] == "binary_gt": return (x > spec["cut"]).astype(float)
    if spec["kind"] == "binary":    return x.astype(float)
    return (x - x.mean()) / x.std(ddof=0)


def fit_all(frame):
    out = {}
    for en, es in EXPOSURES.items():
        f = frame.copy(); f["X"] = expo(f, es)
        for oc in OUTCOMES:
            sub = f.dropna(subset=["X", oc, "purity", "proliferation"])
            if len(sub) < 40 or sub["X"].nunique() < 2: continue
            m = smf.ols(f"{oc} ~ X + purity + proliferation", data=sub).fit()
            out[(en, oc)] = (m.params["X"], m.bse["X"], m.pvalues["X"])
    return out


full2, full3 = fit_all(S2), fit_all(S3)
keys = sorted(set(full2) & set(full3))

# ---- 1. down-sampling S2 to n = |S3| ----
counts = []
for _ in range(200):
    samp = S2.iloc[rng.choice(len(S2), N3, replace=False)]
    r = fit_all(samp)
    ps = [r[k][2] for k in keys if k in r]
    if ps:
        counts.append(int((multipletests(ps, method="fdr_bh")[1] < 0.05).sum()))
counts = np.array(counts)

p3 = [full3[k][2] for k in keys]
n_sig3 = int((multipletests(p3, method="fdr_bh")[1] < 0.05).sum())
p2 = [full2[k][2] for k in keys]
n_sig2 = int((multipletests(p2, method="fdr_bh")[1] < 0.05).sum())

# ---- 2. formal beta-difference test S2 vs S3 ----
zrows = []
for k in keys:
    b2, s2, _ = full2[k]; b3, s3, _ = full3[k]
    se = np.sqrt(s2 ** 2 + s3 ** 2)
    z = (b2 - b3) / se
    zrows.append(dict(exposure=k[0], outcome=k[1], beta_S2=b2, beta_S3=b3, z=z,
                      p_diff=2 * (1 - norm.cdf(abs(z)))))
zdf = pd.DataFrame(zrows)
zdf["q_diff"] = multipletests(zdf["p_diff"], method="fdr_bh")[1]

# ---- 3. bootstrap CI on median effect ratio (S1-significant subset) ----
assoc = pd.read_csv(ROOT / "results" / "associations.tsv", sep="\t")
w = assoc.pivot_table(index=["exposure", "outcome"], columns="setting", values=["beta", "q"]).dropna()
sig = w[w[("q", "S1_TCGA_Xena")] < 0.05]
ratios = (sig[("beta", "S3_nonTCGA_PCAWG")].abs() / sig[("beta", "S1_TCGA_Xena")].abs()).values
boot = [np.median(rng.choice(ratios, len(ratios), replace=True)) for _ in range(4000)]

lines = ["# Is the non-TCGA null a power artifact?", "",
         f"S2 (TCGA, PCAWG RNA) n={len(S2)} vs S3 (non-TCGA, PCAWG RNA) n={N3}.", "",
         "## 1. Down-sampling test (the decisive one)", "",
         f"S2 subsampled to n={N3} (200 draws), counting associations at FDR<0.05:", "",
         f"- **matched-n S2: median {int(np.median(counts))} significant "
         f"[{int(np.percentile(counts,2.5))}–{int(np.percentile(counts,97.5))}]** (of {len(keys)})",
         f"- **actual S3 at the same n: {n_sig3} significant**",
         f"- (full-size S2: {n_sig2})", "",
         f"At identical sample size the TCGA cohort still yields ~{int(np.median(counts))} significant "
         f"associations while the non-TCGA cohort yields {n_sig3}. "
         + ("**The non-TCGA null is therefore NOT explained by power.**"
            if n_sig3 < np.percentile(counts, 2.5) else
            "The difference is within the down-sampling range, so power cannot be excluded."), "",
         "## 2. Formal difference in effect sizes (S2 vs S3)", "",
         f"- associations with significantly different betas (FDR<0.05): "
         f"**{int((zdf.q_diff<0.05).sum())}/{len(zdf)}**",
         f"- median |beta| S2 = {zdf.beta_S2.abs().median():.3f}, S3 = {zdf.beta_S3.abs().median():.3f}", "",
         "## 3. Effect-size attenuation", "",
         f"- median |beta_nonTCGA| / |beta_TCGA| = **{np.median(ratios):.2f}** "
         f"[95% CI {np.percentile(boot,2.5):.2f}–{np.percentile(boot,97.5):.2f}]", ""]
(ROOT / "results" / "power_check.md").write_text("\n".join(lines) + "\n")
zdf.to_csv(ROOT / "results" / "beta_difference.tsv", sep="\t", index=False)
print("\n".join(lines))
