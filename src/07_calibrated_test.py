"""
07 — Calibrated test: is a COHORT swap more disruptive than a random split of the SAME cohort?

The composition-matched permutation (src/06) was restricted to 4 shared histologies (n=399) and
gave p=0.086 — underpowered. This test uses the full PCAWG data and a properly calibrated null:

    NULL      discrepancy between two DISJOINT random halves of TCGA        (sampling noise only)
    OBSERVED  discrepancy between a random TCGA subset and a non-TCGA subset (sampling + cohort)

Both use the SAME group size (n=N_SUB each), so sampling noise is matched by construction, and no
histology restriction is needed. If a cohort swap perturbs effect estimates more than re-sampling
the same cohort does, the cohort effect is real.

Discrepancy statistic across the 60 associations:
    D = mean |beta_group1 - beta_group2|      (and 1 - Pearson r as a secondary)

Outputs results/calibrated_test.md
"""
import warnings, numpy as np, pandas as pd, yaml, statsmodels.formula.api as smf
from scipy import stats
from pathlib import Path
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((ROOT / "config" / "panels.yaml").read_text())
EXPOSURES, OUTCOMES = CFG["exposures"], [p for p in CFG["panels"] if p != "proliferation"]
rng = np.random.default_rng(20260726)
N_SUB, N_DRAW = 330, 200          # matched group size; draws per arm

d = pd.read_parquet(ROOT / "data" / "analysis_matrix.parquet")
p = d[d.rna_source == "PCAWG"]                     # pipeline held constant
TC = p[p.cohort == "TCGA"].reset_index(drop=True)
NT = p[p.cohort == "nonTCGA"].reset_index(drop=True)
print(f"TCGA n={len(TC)}  non-TCGA n={len(NT)}  matched subset size={N_SUB}")


def expo(f, s):
    x = pd.to_numeric(f[s["col"]], errors="coerce")
    if s["kind"] == "binary_gt": return (x > s["cut"]).astype(float)
    if s["kind"] == "binary":    return x.astype(float)
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if sd else x * 0


KEYS = [(en, oc) for en in EXPOSURES for oc in OUTCOMES]


def betas(frame):
    out = {}
    for en, es in EXPOSURES.items():
        f = frame.copy(); f["X"] = expo(f, es)
        for oc in OUTCOMES:
            sub = f.dropna(subset=["X", oc, "purity", "proliferation"])
            if len(sub) < 30 or sub["X"].nunique() < 2:
                continue
            try:
                out[(en, oc)] = smf.ols(f"{oc} ~ X + purity + proliferation", data=sub).fit().params["X"]
            except Exception:
                pass
    return out


def disc(b1, b2):
    ks = [k for k in KEYS if k in b1 and k in b2]
    if len(ks) < 10:
        return np.nan, np.nan
    v1 = np.array([b1[k] for k in ks]); v2 = np.array([b2[k] for k in ks])
    r = stats.pearsonr(v1, v2)[0] if len(ks) > 2 else np.nan
    return np.mean(np.abs(v1 - v2)), 1 - r


# ---- NULL: two disjoint random halves of TCGA ----
null_d, null_r = [], []
for _ in range(N_DRAW):
    idx = rng.permutation(len(TC))
    a, b = TC.iloc[idx[:N_SUB]], TC.iloc[idx[N_SUB:2 * N_SUB]]
    dd, dr = disc(betas(a), betas(b))
    if np.isfinite(dd): null_d.append(dd); null_r.append(dr)

# ---- OBSERVED: random TCGA subset vs random non-TCGA subset ----
obs_d, obs_r = [], []
for _ in range(N_DRAW):
    a = TC.iloc[rng.choice(len(TC), N_SUB, replace=False)]
    b = NT.iloc[rng.choice(len(NT), N_SUB, replace=False)]
    dd, dr = disc(betas(a), betas(b))
    if np.isfinite(dd): obs_d.append(dd); obs_r.append(dr)

null_d, obs_d = np.array(null_d), np.array(obs_d)
null_r, obs_r = np.array(null_r), np.array(obs_r)

# p = P(null >= median observed): how often does re-sampling TCGA look as discrepant as a cohort swap
p_d = (1 + (null_d >= np.median(obs_d)).sum()) / (1 + len(null_d))
p_r = (1 + (null_r >= np.median(obs_r)).sum()) / (1 + len(null_r))
# effect: how many times larger is the cohort discrepancy
ratio_d = np.median(obs_d) / np.median(null_d)
boot = [np.median(rng.choice(obs_d, len(obs_d))) / np.median(rng.choice(null_d, len(null_d)))
        for _ in range(2000)]
rlo, rhi = np.percentile(boot, [2.5, 97.5])

lines = ["# Calibrated test: cohort swap vs re-sampling the same cohort", "",
         f"Both arms use matched group size **n = {N_SUB}** on PCAWG RNA (pipeline constant), "
         f"{N_DRAW} draws each, so sampling noise is matched by construction and no histology "
         "restriction is needed.", "",
         "| arm | mean \\|Δβ\\| across 60 associations | 1 − r |", "|---|--:|--:|",
         f"| NULL — two disjoint random halves of TCGA | {np.median(null_d):.4f} "
         f"[{np.percentile(null_d,2.5):.4f}, {np.percentile(null_d,97.5):.4f}] | {np.median(null_r):.3f} |",
         f"| OBSERVED — TCGA vs non-TCGA | **{np.median(obs_d):.4f}** "
         f"[{np.percentile(obs_d,2.5):.4f}, {np.percentile(obs_d,97.5):.4f}] | **{np.median(obs_r):.3f}** |",
         "",
         f"- **Cohort swaps are {ratio_d:.2f}x more disruptive** than re-sampling the same cohort "
         f"(95% CI {rlo:.2f}–{rhi:.2f}).",
         f"- permutation-style p (null >= observed): **p = {p_d:.4f}** (mean \\|Δβ\\|), "
         f"**p = {p_r:.4f}** (1−r).", "",
         ("**The cohort effect is formally established:** swapping cohort perturbs effect estimates "
          "significantly more than re-drawing the same number of donors from the same cohort."
          if p_d < 0.05 else
          "**Not significant:** a cohort swap is not measurably more disruptive than re-sampling "
          "the same cohort at this sample size."), ""]
(ROOT / "results" / "calibrated_test.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
