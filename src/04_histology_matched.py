"""
04 — Is the "cohort effect" really a TUMOUR-TYPE COMPOSITION effect?

The non-TCGA cohort is 38% lymphoid (Lymph-BNHL/CLL), a lineage absent from TCGA here, and its
histology mix differs entirely. If the replication failure is just composition, then comparing
TCGA and non-TCGA donors *within the same histologies* should remove it.

Comparisons (all on PCAWG RNA, so pipeline is held constant):
  - S2 vs S3 restricted to SHARED histologies
  - S3 solid-only
  - per-shared-histology TCGA vs non-TCGA

Outputs results/histology_matched.md
"""
import warnings, numpy as np, pandas as pd, yaml, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from pathlib import Path
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((ROOT / "config" / "panels.yaml").read_text())
EXPOSURES, OUTCOMES = CFG["exposures"], [p for p in CFG["panels"] if p != "proliferation"]

d = pd.read_parquet(ROOT / "data" / "analysis_matrix.parquet")
p = d[d.rna_source == "PCAWG"]
S2, S3 = p[p.cohort == "TCGA"], p[p.cohort == "nonTCGA"]

h2, h3 = set(S2.histo), set(S3.histo)
shared = sorted(h2 & h3)
print(f"TCGA histologies: {len(h2)} | non-TCGA: {len(h3)} | SHARED: {len(shared)} -> {shared}")
print(f"\nnon-TCGA composition: {S3.histo.value_counts().to_dict()}")


def expo(f, s):
    x = pd.to_numeric(f[s["col"]], errors="coerce")
    if s["kind"] == "binary_gt": return (x > s["cut"]).astype(float)
    if s["kind"] == "binary":    return x.astype(float)
    return (x - x.mean()) / x.std(ddof=0)


def audit(frame, label):
    res = []
    for en, es in EXPOSURES.items():
        f = frame.copy(); f["X"] = expo(f, es)
        for oc in OUTCOMES:
            sub = f.dropna(subset=["X", oc, "purity", "proliferation"])
            if len(sub) < 40 or sub["X"].nunique() < 2: continue
            m = smf.ols(f"{oc} ~ X + purity + proliferation", data=sub).fit()
            res.append(dict(setting=label, exposure=en, outcome=oc, n=len(sub),
                            beta=m.params["X"], p=m.pvalues["X"]))
    r = pd.DataFrame(res)
    if len(r): r["q"] = multipletests(r["p"], method="fdr_bh")[1]
    return r


sets = {
    "TCGA (all histologies)": S2,
    "TCGA (shared histologies only)": S2[S2.histo.isin(shared)],
    "nonTCGA (all)": S3,
    "nonTCGA (shared histologies only)": S3[S3.histo.isin(shared)],
    "nonTCGA (solid only)": S3[~S3.histo.str.startswith("Lymph")],
}
allr, summary = [], []
for lab, frame in sets.items():
    r = audit(frame, lab)
    allr.append(r)
    if len(r):
        summary.append(dict(setting=lab, n_donors=len(frame), n_assoc=len(r),
                            n_sig=int((r["q"] < 0.05).sum()),
                            median_abs_beta=r["beta"].abs().median()))
sm = pd.DataFrame(summary)
print("\n" + sm.to_string(index=False))

# per shared histology, TCGA vs non-TCGA (WGD -> cytolytic as the canonical association)
rows = []
for h in shared:
    for lab, frame in [("TCGA", S2), ("nonTCGA", S3)]:
        sub = frame[frame.histo == h].copy()
        sub["X"] = (pd.to_numeric(sub["ploidy"], errors="coerce") > 2.5).astype(float)
        sub = sub.dropna(subset=["X", "cytolytic", "purity", "proliferation"])
        if len(sub) >= 30 and sub["X"].nunique() > 1:
            m = smf.ols("cytolytic ~ X + purity + proliferation", data=sub).fit()
            rows.append(dict(histology=h, cohort=lab, n=len(sub), wgd_pos=int(sub["X"].sum()),
                             beta=m.params["X"], p=m.pvalues["X"]))
ph = pd.DataFrame(rows)

lines = ["# Is the cohort effect just tumour-type composition?", "",
         f"Shared histologies between TCGA and non-TCGA (PCAWG RNA): **{shared}**", "",
         "## Associations significant (FDR<0.05) per setting", "",
         "| setting | donors | associations | significant | median \\|beta\\| |", "|---|--:|--:|--:|--:|"]
for _, r in sm.iterrows():
    lines.append(f"| {r.setting} | {r.n_donors} | {r.n_assoc} | {r.n_sig} | {r.median_abs_beta:.3f} |")
lines += ["", "**Reading:** if the non-TCGA null were driven by its different histology mix, "
          "restricting to shared histologies (or removing lymphoid tumours) should recover signal.", ""]
if len(ph):
    lines += ["## Canonical association (WGD -> cytolytic) per shared histology", "",
              "| histology | cohort | n | WGD+ | beta | p |", "|---|---|--:|--:|--:|--:|"]
    for _, r in ph.iterrows():
        lines.append(f"| {r.histology} | {r.cohort} | {r.n} | {r.wgd_pos} | {r.beta:+.3f} | {r.p:.3g} |")
(ROOT / "results" / "histology_matched.md").write_text("\n".join(lines) + "\n")
pd.concat(allr).to_csv(ROOT / "results" / "histology_matched.tsv", sep="\t", index=False)
if len(ph):
    print("\nWGD -> cytolytic per shared histology:")
    print(ph.to_string(index=False))
print("\nwrote results/histology_matched.md")
