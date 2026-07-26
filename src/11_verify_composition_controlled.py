"""
11 — VERIFICATION: is the self-vs-cross replication gap really about COHORT, or about the two
cohorts containing different TUMOUR TYPES?

src/10 compared TCGA (23 histologies) with non-TCGA (8 histologies, 38% lymphoid), so its
62% vs 8% gap confounds "different consortium" with "different tumour types". This removes that
confound in the strictest available way:

  Both replication sets are built from the SAME histologies with the SAME number of donors per
  histology (per-histology 1:1 matching). Discovery is identical for both arms. The ONLY difference
  between the two arms is which consortium the replication donors came from.

If the gap survives, the cohort effect is real. If it collapses, the effect was composition.

Outputs results/verify_composition_controlled.md
"""
import warnings, numpy as np, pandas as pd, yaml, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from pathlib import Path
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((ROOT / "config" / "panels.yaml").read_text())
EXPOSURES, OUTCOMES = CFG["exposures"], [p for p in CFG["panels"] if p != "proliferation"]
rng = np.random.default_rng(20260726)
N_DRAW, N_DISC = 100, 340

d = pd.read_parquet(ROOT / "data" / "analysis_matrix.parquet")
p = d[d.rna_source == "PCAWG"].copy()
TC, NT = p[p.cohort == "TCGA"], p[p.cohort == "nonTCGA"]
shared = sorted(set(TC.histo) & set(NT.histo))

# per-histology 1:1 matched pool sizes
cnt = {h: min((TC.histo == h).sum(), (NT.histo == h).sum()) for h in shared}
cnt = {h: c for h, c in cnt.items() if c >= 15}
n_rep = sum(cnt.values())
print(f"shared histologies usable: {cnt}")
print(f"matched replication-set size per arm: {n_rep}")

TCs, NTs = TC[TC.histo.isin(cnt)], NT[NT.histo.isin(cnt)]
TC_rest = TC[~TC.index.isin(TCs.index)]          # discovery pool: TCGA donors NOT in the matched pool
print(f"discovery pool (TCGA, disjoint from replication pool): {len(TC_rest)}")


def expo(f, s):
    x = pd.to_numeric(f[s["col"]], errors="coerce")
    if s["kind"] == "binary_gt": return (x > s["cut"]).astype(float)
    if s["kind"] == "binary":    return x.astype(float)
    sd = x.std(ddof=0); return (x - x.mean()) / sd if sd else x * 0


def fits(frame):
    out = {}
    for en, es in EXPOSURES.items():
        f = frame.copy(); f["X"] = expo(f, es)
        for oc in OUTCOMES:
            sub = f.dropna(subset=["X", oc, "purity", "proliferation"])
            if len(sub) < 25 or sub["X"].nunique() < 2: continue
            try:
                m = smf.ols(f"{oc} ~ X + purity + proliferation", data=sub).fit()
                out[(en, oc)] = (m.params["X"], m.pvalues["X"])
            except Exception: pass
    return out


def rep_rate(bd, rep):
    br = fits(rep)
    ks = sorted(set(bd) & set(br))
    if len(ks) < 10: return np.nan
    q = multipletests([bd[k][1] for k in ks], method="fdr_bh")[1]
    sig = [k for k, qq in zip(ks, q) if qq < 0.05]
    if not sig: return np.nan
    return sum(1 for k in sig if br[k][1] < 0.05 and np.sign(br[k][0]) == np.sign(bd[k][0])) / len(sig)


def matched_sample(frame):
    """draw cnt[h] donors from each shared histology -> identical composition in both arms"""
    parts = []
    for h, c in cnt.items():
        sub = frame[frame.histo == h]
        parts.append(sub.iloc[rng.choice(len(sub), c, replace=False)])
    return pd.concat(parts)


self_r, cross_r = [], []
for _ in range(N_DRAW):
    disc = TC_rest.iloc[rng.choice(len(TC_rest), min(N_DISC, len(TC_rest)), replace=False)]
    bd = fits(disc)
    a = rep_rate(bd, matched_sample(TCs))
    b = rep_rate(bd, matched_sample(NTs))
    if np.isfinite(a): self_r.append(a)
    if np.isfinite(b): cross_r.append(b)

self_r, cross_r = np.array(self_r), np.array(cross_r)
gap = np.median(self_r) - np.median(cross_r)
p_overlap = (self_r <= np.median(cross_r)).mean() if len(self_r) else np.nan

lines = ["# Verification: cohort effect with tumour-type composition held identical", "",
         f"Both replication arms are drawn from the **same histologies with the same number of "
         f"donors each** ({cnt}), total n = {n_rep} per arm. Discovery is the same TCGA sample in "
         f"both arms (n = {min(N_DISC, len(TC_rest))}, disjoint from the replication pool). "
         f"The only difference between arms is the consortium.", "",
         "| replication set | draws | median replication rate [95%] |", "|---|--:|--:|",
         f"| TCGA (self) | {len(self_r)} | **{np.median(self_r):.0%}** "
         f"[{np.percentile(self_r,2.5):.0%}–{np.percentile(self_r,97.5):.0%}] |" if len(self_r) else "| TCGA (self) | 0 | n/a |",
         f"| non-TCGA (cross) | {len(cross_r)} | **{np.median(cross_r):.0%}** "
         f"[{np.percentile(cross_r,2.5):.0%}–{np.percentile(cross_r,97.5):.0%}] |" if len(cross_r) else "| non-TCGA (cross) | 0 | n/a |",
         "", f"- gap (self − cross): **{gap:+.0%}**",
         f"- fraction of self-draws as poor as the median cross-draw: **{p_overlap:.0%}**", "",
         ("**The gap survives with composition held identical — the cohort effect is real and is "
          "not explained by differing tumour types.**" if gap > 0.15 and p_overlap < 0.10 else
          "**The gap does NOT survive composition matching at this sample size** — with identical "
          "tumour-type composition the two arms are not clearly distinguishable, so 'cohort' and "
          "'tumour-type composition' cannot be separated in this data."), ""]
(ROOT / "results" / "verify_composition_controlled.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
