"""
02 — The cross-cohort reproducibility audit.

Every genomic-instability exposure x immune outcome association is estimated in three settings:

  S1  TCGA donors  | TCGA-Xena RNA   <- the "classic TCGA discovery" setting
  S2  TCGA donors  | PCAWG RNA       <- SAME donors, SAME exposures, different RNA pipeline
  S3  non-TCGA     | PCAWG RNA       <- SAME pipeline, different donors

  S1 vs S2  isolates the PIPELINE effect   (donors held constant)
  S2 vs S3  isolates the COHORT effect     (pipeline held constant)

For the paired pipeline arm the analysis is restricted to donors present in BOTH S1 and S2.
Model: outcome ~ exposure + purity + proliferation (scores already within-histology z-scored).

Outputs results/associations.tsv and results/reproducibility_metrics.md
"""
import warnings, numpy as np, pandas as pd, yaml, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy import stats
from pathlib import Path
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((ROOT / "config" / "panels.yaml").read_text())
EXPOSURES = CFG["exposures"]
OUTCOMES = [p for p in CFG["panels"] if p != "proliferation"]

d = pd.read_parquet(ROOT / "data" / "analysis_matrix.parquet")

# paired donor set for the pipeline arm
xena_ids = set(d[(d.rna_source == "TCGA_Xena")].donor_unique_id)
pcawg_tcga_ids = set(d[(d.rna_source == "PCAWG") & (d.cohort == "TCGA")].donor_unique_id)
PAIRED = xena_ids & pcawg_tcga_ids
print(f"paired TCGA donors (in both RNA sources): {len(PAIRED)}")

SETTINGS = {
    "S1_TCGA_Xena":   d[(d.rna_source == "TCGA_Xena") & d.donor_unique_id.isin(PAIRED)],
    "S2_TCGA_PCAWG":  d[(d.rna_source == "PCAWG") & d.donor_unique_id.isin(PAIRED)],
    "S3_nonTCGA_PCAWG": d[(d.rna_source == "PCAWG") & (d.cohort == "nonTCGA")],
}
for k, v in SETTINGS.items():
    print(f"  {k}: n={len(v)}")


def exposure_series(frame, spec):
    col, kind = spec["col"], spec["kind"]
    x = pd.to_numeric(frame[col], errors="coerce")
    if kind == "binary_gt":
        return (x > spec["cut"]).astype(float)
    if kind == "binary":
        return x.astype(float)
    return (x - x.mean()) / x.std(ddof=0)          # standardise continuous exposures


rows = []
for sname, sdf in SETTINGS.items():
    for ename, espec in EXPOSURES.items():
        s = sdf.copy()
        s["X"] = exposure_series(s, espec)
        for oc in OUTCOMES:
            sub = s.dropna(subset=["X", oc, "purity", "proliferation"])
            if len(sub) < 40 or sub["X"].nunique() < 2:
                continue
            m = smf.ols(f"{oc} ~ X + purity + proliferation", data=sub).fit()
            ci = m.conf_int().loc["X"]
            rows.append(dict(setting=sname, exposure=ename, outcome=oc, n=len(sub),
                             beta=m.params["X"], lo=ci[0], hi=ci[1], p=m.pvalues["X"]))

res = pd.DataFrame(rows)
res["q"] = np.nan
for s in res.setting.unique():                      # FDR within setting
    m = res.setting == s
    res.loc[m, "q"] = multipletests(res.loc[m, "p"], method="fdr_bh")[1]
res.to_csv(ROOT / "results" / "associations.tsv", sep="\t", index=False)
print(f"\nfitted {len(res)} associations "
      f"({res.exposure.nunique()} exposures x {res.outcome.nunique()} outcomes x {res.setting.nunique()} settings)")

# ---------------- reproducibility metrics ----------------
w = res.pivot_table(index=["exposure", "outcome"], columns="setting",
                    values=["beta", "p", "q"]).dropna()
S1, S2, S3 = "S1_TCGA_Xena", "S2_TCGA_PCAWG", "S3_nonTCGA_PCAWG"


def compare(a, b, label, disc_sig_only=True):
    """replication metrics going from setting a (discovery) to setting b (replication)"""
    sub = w if not disc_sig_only else w[w[("q", a)] < 0.05]
    if len(sub) == 0:
        return None
    ba, bb = sub[("beta", a)], sub[("beta", b)]
    sign = (np.sign(ba) == np.sign(bb)).mean()
    rep = ((np.sign(ba) == np.sign(bb)) & (sub[("p", b)] < 0.05)).mean()
    r = stats.pearsonr(ba, bb)[0] if len(sub) > 2 else np.nan
    shrink = (bb.abs() / ba.abs()).median()
    return dict(comparison=label, n_assoc=len(sub), sign_agreement=sign,
                replication_rate=rep, beta_correlation=r, median_effect_ratio=shrink)


met = [compare(S1, S2, "PIPELINE  (same donors: Xena -> PCAWG RNA)"),
       compare(S2, S3, "COHORT    (same pipeline: TCGA -> non-TCGA)"),
       compare(S1, S3, "BOTH      (Xena/TCGA -> PCAWG/non-TCGA)")]
met = [m for m in met if m]
mdf = pd.DataFrame(met)

# how many associations are significant in each setting?
sig_counts = {s: int((res[res.setting == s]["q"] < 0.05).sum()) for s in SETTINGS}
tot = len(w)

lines = ["# Cross-cohort reproducibility of aneuploidy-immune associations", "",
         f"Systematic audit of **{res.exposure.nunique()} genomic-instability exposures x "
         f"{res.outcome.nunique()} immune outcomes = {tot} associations**, each estimated in three settings.",
         "", "| setting | donors | RNA pipeline | n | significant (FDR<0.05) |", "|---|---|---|--:|--:|",
         f"| S1 | TCGA | TCGA-Xena | {len(SETTINGS[S1])} | {sig_counts[S1]}/{tot} |",
         f"| S2 | TCGA (same as S1) | PCAWG | {len(SETTINGS[S2])} | {sig_counts[S2]}/{tot} |",
         f"| S3 | non-TCGA | PCAWG | {len(SETTINGS[S3])} | {sig_counts[S3]}/{tot} |", "",
         "## Replication of the associations that are significant in S1 (FDR<0.05)", "",
         "| comparison | n | sign agreement | replication rate | beta correlation | median effect ratio |",
         "|---|--:|--:|--:|--:|--:|"]
for m in met:
    lines.append(f"| {m['comparison']} | {m['n_assoc']} | {m['sign_agreement']:.0%} | "
                 f"{m['replication_rate']:.0%} | {m['beta_correlation']:.2f} | {m['median_effect_ratio']:.2f} |")
lines += ["", "*replication rate = same sign AND p<0.05 in the replication setting.*",
          "*median effect ratio = |beta_replication| / |beta_discovery|; <1 means effects shrink.*"]
(ROOT / "results" / "reproducibility_metrics.md").write_text("\n".join(lines) + "\n")

print("\n" + "=" * 74)
print(f"significant (FDR<0.05): S1={sig_counts[S1]}/{tot}  S2={sig_counts[S2]}/{tot}  S3={sig_counts[S3]}/{tot}")
print("=" * 74)
print(mdf.to_string(index=False))
print("\nwrote results/associations.tsv + results/reproducibility_metrics.md")
