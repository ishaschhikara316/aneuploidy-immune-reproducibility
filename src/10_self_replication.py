"""
10 — THE DECISIVE TEST: does TCGA replicate ITSELF?

The alarming headline of any replication study is the rate at which discovery findings stay
significant in the replication cohort (here: 19%, and 0/60 at full size). That number is only
meaningful if a cohort would replicate ITSELF at a much higher rate under the same conditions.

So we hold everything constant (same RNA pipeline, same exposures, same models, same sample sizes)
and vary only WHERE the replication sample comes from:

  ARM 1  SELF   : TCGA discovery  ->  a DISJOINT random TCGA replication set
  ARM 2  CROSS  : TCGA discovery  ->  a non-TCGA replication set
  ARM 3  SELF, prevalence-matched: as ARM 1, but the TCGA replication set is down-sampled to the
         non-TCGA whole-genome-doubling prevalence (23%), isolating the effect of how common the
         exposure is from the effect of which cohort it came from.

If ARM 1 ~ ARM 2, the "cohort effect" is an illusion of the metric plus sampling.
Replication rule (the field standard): of associations significant at FDR<0.05 in discovery,
the fraction with p<0.05 AND the same sign in the replication set.

Outputs results/self_replication.md
"""
import warnings, numpy as np, pandas as pd, yaml, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from pathlib import Path
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((ROOT / "config" / "panels.yaml").read_text())
EXPOSURES, OUTCOMES = CFG["exposures"], [p for p in CFG["panels"] if p != "proliferation"]
rng = np.random.default_rng(20260726)
N_SUB, N_DRAW = 340, 100

d = pd.read_parquet(ROOT / "data" / "analysis_matrix.parquet")
p = d[d.rna_source == "PCAWG"].copy()
p["WGDbin"] = (pd.to_numeric(p["ploidy"], errors="coerce") > 2.5).astype(float)
TC = p[p.cohort == "TCGA"].reset_index(drop=True)
NT = p[p.cohort == "nonTCGA"].reset_index(drop=True)
prev_nt = NT["WGDbin"].mean()
print(f"TCGA n={len(TC)} (WGD+ {TC.WGDbin.mean():.1%})   non-TCGA n={len(NT)} (WGD+ {prev_nt:.1%})")


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
            if len(sub) < 30 or sub["X"].nunique() < 2: continue
            try:
                m = smf.ols(f"{oc} ~ X + purity + proliferation", data=sub).fit()
                out[(en, oc)] = (m.params["X"], m.pvalues["X"])
            except Exception: pass
    return out


def rep_rate(disc, rep):
    bd, br = fits(disc), fits(rep)
    ks = sorted(set(bd) & set(br))
    if len(ks) < 10: return np.nan, 0
    q = multipletests([bd[k][1] for k in ks], method="fdr_bh")[1]
    sig = [k for k, qq in zip(ks, q) if qq < 0.05]
    if not sig: return np.nan, 0
    ok = sum(1 for k in sig if br[k][1] < 0.05 and np.sign(br[k][0]) == np.sign(bd[k][0]))
    return ok / len(sig), len(sig)


def match_prev(frame, target, n):
    """subsample `frame` to <=n rows with WGD+ prevalence ~= target.
    Shrinks n to whatever the available positives/negatives can support."""
    pos, neg = frame[frame.WGDbin == 1], frame[frame.WGDbin == 0]
    if target <= 0 or target >= 1: return None
    n = int(min(n, len(pos) / target, len(neg) / (1 - target)))
    if n < 120: return None
    npos = int(round(n * target)); nneg = n - npos
    if len(pos) < npos or len(neg) < nneg: return None
    return pd.concat([pos.iloc[rng.choice(len(pos), npos, replace=False)],
                      neg.iloc[rng.choice(len(neg), nneg, replace=False)]])


arms = {"SELF  (TCGA -> disjoint TCGA)": [],
        "CROSS (TCGA -> non-TCGA)": [],
        "SELF, prevalence-matched": []}
nsig = {k: [] for k in arms}

for _ in range(N_DRAW):
    idx = rng.permutation(len(TC))
    disc = TC.iloc[idx[:N_SUB]]
    rep_self = TC.iloc[idx[N_SUB:2 * N_SUB]]
    r, k = rep_rate(disc, rep_self)
    if np.isfinite(r): arms["SELF  (TCGA -> disjoint TCGA)"].append(r); nsig["SELF  (TCGA -> disjoint TCGA)"].append(k)

    rep_cross = NT.iloc[rng.choice(len(NT), N_SUB, replace=False)]
    r, k = rep_rate(disc, rep_cross)
    if np.isfinite(r): arms["CROSS (TCGA -> non-TCGA)"].append(r); nsig["CROSS (TCGA -> non-TCGA)"].append(k)

    rp = match_prev(TC.iloc[idx[N_SUB:]], prev_nt, N_SUB)
    if rp is not None:
        r, k = rep_rate(disc, rp)
        if np.isfinite(r): arms["SELF, prevalence-matched"].append(r); nsig["SELF, prevalence-matched"].append(k)

rows = []
for k, v in arms.items():
    v = np.array(v)
    if len(v) == 0:
        rows.append(dict(arm=k, draws=0, median_rep_rate=np.nan, lo=np.nan, hi=np.nan,
                         median_n_disc_sig=0))
        continue
    rows.append(dict(arm=k, draws=len(v), median_rep_rate=np.median(v),
                     lo=np.percentile(v, 2.5), hi=np.percentile(v, 97.5),
                     median_n_disc_sig=int(np.median(nsig[k])) if nsig[k] else 0))
df = pd.DataFrame(rows)

self_r = df[df.arm.str.startswith("SELF  ")]["median_rep_rate"].iloc[0]
cross_r = df[df.arm.str.startswith("CROSS")]["median_rep_rate"].iloc[0]
sm = df[df.arm.str.startswith("SELF, prev")]["median_rep_rate"]
selfm_r = sm.iloc[0] if len(sm) else np.nan
# how often does a SELF draw look as bad as the median CROSS draw?
sa = np.array(arms["SELF  (TCGA -> disjoint TCGA)"])
p_overlap = (sa <= cross_r).mean()

lines = ["# Does TCGA replicate itself?", "",
         f"Everything held constant (same pipeline, exposures, models, n={N_SUB} per group, "
         f"{N_DRAW} draws); only the SOURCE of the replication sample changes.", "",
         "Replication rule (field standard): of associations at FDR<0.05 in discovery, the fraction "
         "with p<0.05 and the same sign in the replication set.", "",
         "| arm | draws | median discovery hits | **median replication rate** [95%] |",
         "|---|--:|--:|--:|"]
for _, r in df.iterrows():
    lines.append(f"| {r.arm} | {r.draws} | {r.median_n_disc_sig} | " +
                 ("n/a |" if not np.isfinite(r.median_rep_rate) else
                  f"**{r.median_rep_rate:.0%}** [{r.lo:.0%}–{r.hi:.0%}] |"))
lines += ["",
          f"- TCGA replicating **itself**: **{self_r:.0%}**",
          f"- TCGA replicating in **non-TCGA**: **{cross_r:.0%}**",
          f"- TCGA replicating itself at non-TCGA **exposure prevalence**: " + ("n/a" if not np.isfinite(selfm_r) else f"**{selfm_r:.0%}**"),
          f"- fraction of SELF draws that look at least as 'failed' as the typical CROSS draw: "
          f"**{p_overlap:.0%}**", "",
          ("**A cohort replicating ITSELF fails the standard replication test at a similar rate to "
           "replicating in a different cohort.** The alarming replication rate is therefore mostly a "
           "property of the metric and sample size, not of cohort differences."
           if self_r - cross_r < 0.20 else
           "Self-replication is markedly higher than cross-cohort replication, so the cohort "
           "difference is not purely a metric artefact."), ""]
(ROOT / "results" / "self_replication.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
