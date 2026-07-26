"""
06 — A single, properly-powered FORMAL test of the cohort effect.

Weakness this fixes: per-association beta differences (S2 vs S3) were individually
non-significant (0/60) because each test uses only its own association. Pooling across the 60
associations gives one well-powered test — but the associations are NOT independent (they share
donors), so the pooling must resample/permute DONORS, not associations.

Primary test — STRATIFIED PERMUTATION (also controls tumour-type composition):
  Working only on PCAWG RNA (pipeline held constant) and only within histologies shared by both
  cohorts, permute the TCGA / non-TCGA label *within each histology*, recompute all associations,
  and recompute the pooled statistic. This tests the cohort effect with composition held fixed.

Statistic: betas are oriented by the sign of the TCGA-Xena (S1) discovery estimate so the
"expected direction" is positive, then
      Delta = mean(oriented beta | TCGA) - mean(oriented beta | non-TCGA)
Delta > 0 means effects are systematically larger in TCGA.

Secondary: donor-level bootstrap CI on Delta and on the effect ratio.

Outputs results/formal_test.md
"""
import warnings, numpy as np, pandas as pd, yaml, statsmodels.formula.api as smf
from pathlib import Path
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((ROOT / "config" / "panels.yaml").read_text())
EXPOSURES, OUTCOMES = CFG["exposures"], [p for p in CFG["panels"] if p != "proliferation"]
rng = np.random.default_rng(20260726)
N_PERM, N_BOOT = 500, 1000

d = pd.read_parquet(ROOT / "data" / "analysis_matrix.parquet")
p = d[d.rna_source == "PCAWG"].copy()                     # pipeline held constant throughout
shared = sorted(set(p[p.cohort == "TCGA"].histo) & set(p[p.cohort == "nonTCGA"].histo))
ps = p[p.histo.isin(shared)].copy()
print(f"shared histologies: {shared}")
print(f"stratified cohort: n={len(ps)} "
      f"(TCGA {int((ps.cohort=='TCGA').sum())} / non-TCGA {int((ps.cohort=='nonTCGA').sum())})")

# orientation from the S1 discovery setting
a = pd.read_csv(ROOT / "results" / "associations.tsv", sep="\t")
s1 = a[a.setting == "S1_TCGA_Xena"].set_index(["exposure", "outcome"])["beta"]
ORIENT = {k: (1.0 if v >= 0 else -1.0) for k, v in s1.items()}


def expo(f, s):
    x = pd.to_numeric(f[s["col"]], errors="coerce")
    if s["kind"] == "binary_gt": return (x > s["cut"]).astype(float)
    if s["kind"] == "binary":    return x.astype(float)
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if sd else x * 0


def oriented_mean(frame):
    """mean oriented beta across all estimable associations in `frame`"""
    vals = []
    for en, es in EXPOSURES.items():
        f = frame.copy(); f["X"] = expo(f, es)
        for oc in OUTCOMES:
            sub = f.dropna(subset=["X", oc, "purity", "proliferation"])
            if len(sub) < 30 or sub["X"].nunique() < 2:
                continue
            try:
                b = smf.ols(f"{oc} ~ X + purity + proliferation", data=sub).fit().params["X"]
            except Exception:
                continue
            vals.append(b * ORIENT.get((en, oc), 1.0))
    return np.mean(vals) if vals else np.nan


def delta(frame, labels):
    t = oriented_mean(frame[labels == "TCGA"])
    n = oriented_mean(frame[labels == "nonTCGA"])
    return t - n, t, n


obs_delta, obs_t, obs_n = delta(ps, ps["cohort"].values)
print(f"\nOBSERVED (shared histologies): mean oriented beta "
      f"TCGA={obs_t:+.4f}  nonTCGA={obs_n:+.4f}  Delta={obs_delta:+.4f}")

# ---------- primary: stratified permutation within histology ----------
null = []
hist = ps["histo"].values
for i in range(N_PERM):
    lab = ps["cohort"].values.copy()
    for h in shared:                                       # permute WITHIN each histology
        m = hist == h
        lab[m] = rng.permutation(lab[m])
    dv, _, _ = delta(ps, lab)
    if np.isfinite(dv):
        null.append(dv)
null = np.array(null)
p_perm = (1 + (null >= obs_delta).sum()) / (1 + len(null))
print(f"permutation null: mean={null.mean():+.4f}  95th pct={np.percentile(null,95):+.4f}  "
      f"max={null.max():+.4f}")
print(f"one-sided permutation p = {p_perm:.4f}  ({len(null)} permutations)")

# ---------- secondary: donor bootstrap CI on Delta and ratio ----------
bt = []
tc = ps[ps.cohort == "TCGA"]; nt = ps[ps.cohort == "nonTCGA"]
for _ in range(N_BOOT):
    a1 = tc.iloc[rng.choice(len(tc), len(tc), replace=True)]
    a2 = nt.iloc[rng.choice(len(nt), len(nt), replace=True)]
    t, n = oriented_mean(a1), oriented_mean(a2)
    if np.isfinite(t) and np.isfinite(n):
        bt.append((t - n, n / t if t else np.nan))
bt = np.array(bt)
d_lo, d_hi = np.nanpercentile(bt[:, 0], [2.5, 97.5])
r_lo, r_hi = np.nanpercentile(bt[:, 1], [2.5, 97.5])
ratio = obs_n / obs_t if obs_t else np.nan

lines = ["# Formal test of the cohort effect (pooled, composition-controlled)", "",
         "Per-association tests were individually underpowered (0/60 significant). This pools "
         "information across all associations while resampling **donors**, so the correlation "
         "between associations is respected.", "",
         f"Restricted to histologies shared by both cohorts (**{', '.join(shared)}**) and to PCAWG "
         f"RNA only, so **tumour-type composition and pipeline are both held constant**. "
         f"n = {len(ps)} (TCGA {int((ps.cohort=='TCGA').sum())} / non-TCGA "
         f"{int((ps.cohort=='nonTCGA').sum())}).", "",
         "Betas oriented by the sign of the TCGA-Xena discovery estimate; "
         "Delta = mean(oriented beta, TCGA) − mean(oriented beta, non-TCGA).", "",
         "## Primary — stratified permutation test", "",
         f"- observed mean oriented beta: **TCGA {obs_t:+.4f}**, **non-TCGA {obs_n:+.4f}**",
         f"- **observed Delta = {obs_delta:+.4f}**",
         f"- permutation null ({len(null)} draws, cohort label shuffled *within* histology): "
         f"mean {null.mean():+.4f}, 95th percentile {np.percentile(null,95):+.4f}, max {null.max():+.4f}",
         f"- **one-sided permutation p = {p_perm:.4f}**", "",
         ("**The cohort effect is formally significant** with composition and pipeline held constant."
          if p_perm < 0.05 else
          "**The cohort effect is NOT formally significant** once composition is held constant — "
          "the aggregate attenuation cannot be distinguished from chance by this test."), "",
         "## Secondary — donor bootstrap", "",
         f"- Delta = {obs_delta:+.4f} [95% CI {d_lo:+.4f}, {d_hi:+.4f}]",
         f"- effect ratio (non-TCGA / TCGA) = {ratio:.2f} [95% CI {r_lo:.2f}, {r_hi:.2f}]", ""]
(ROOT / "results" / "formal_test.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines[-6:]))
print("\nwrote results/formal_test.md")
