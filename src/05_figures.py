"""05 — figures for the cross-cohort reproducibility audit."""
import warnings, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)
BLUE, RED, GREY = "#4C6EF5", "#E03131", "#adb5bd"

a = pd.read_csv(ROOT / "results" / "associations.tsv", sep="\t")
w = a.pivot_table(index=["exposure", "outcome"], columns="setting",
                  values=["beta", "p", "q"]).dropna()
S1, S2, S3 = "S1_TCGA_Xena", "S2_TCGA_PCAWG", "S3_nonTCGA_PCAWG"
sig1 = w[("q", S1)] < 0.05

# ---- Fig 1: the decomposition — pipeline vs cohort ----
fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4), sharex=True, sharey=True)
for ax, (xk, yk, title, note) in zip(axes, [
        (S1, S2, "PIPELINE effect\nsame donors, different RNA pipeline", "Xena → PCAWG"),
        (S2, S3, "COHORT effect\nsame pipeline, different donors", "TCGA → non-TCGA")]):
    x, y = w[("beta", xk)], w[("beta", yk)]
    ax.scatter(x[~sig1], y[~sig1], s=22, color=GREY, alpha=.55, label="ns in TCGA")
    ax.scatter(x[sig1], y[sig1], s=30, color=BLUE if xk == S1 else RED,
               alpha=.85, label="significant in TCGA")
    lim = max(abs(x).max(), abs(y).max()) * 1.12
    ax.plot([-lim, lim], [-lim, lim], color="#333", lw=1, ls="--", zorder=0)
    ax.axhline(0, color=GREY, lw=.8, zorder=0); ax.axvline(0, color=GREY, lw=.8, zorder=0)
    r = stats.pearsonr(x, y)[0]
    ratio = (y[sig1].abs() / x[sig1].abs()).median()
    ax.set_title(f"{title}\nr = {r:.2f}   median effect ratio = {ratio:.2f}", fontsize=10)
    ax.set_xlabel(f"beta ({note.split(' → ')[0]})"); ax.set_ylabel(f"beta ({note.split(' → ')[1]})")
    ax.legend(fontsize=7, loc="upper left")
fig.suptitle("Aneuploidy–immune associations are pipeline-robust but cohort-fragile", fontsize=12)
fig.tight_layout(); fig.savefig(FIG / "fig1_pipeline_vs_cohort.png", dpi=160); plt.close(fig)

# ---- Fig 2: replication metrics ----
def metrics(xk, yk):
    x, y = w[("beta", xk)][sig1], w[("beta", yk)][sig1]
    p = w[("p", yk)][sig1]
    return dict(sign=(np.sign(x) == np.sign(y)).mean(),
                rep=((np.sign(x) == np.sign(y)) & (p < .05)).mean(),
                r=stats.pearsonr(w[("beta", xk)], w[("beta", yk)])[0])
mp, mc = metrics(S1, S2), metrics(S2, S3)
fig, ax = plt.subplots(figsize=(6.6, 3.4))
labels = ["sign\nagreement", "replication rate\n(same sign & p<0.05)", "beta\ncorrelation"]
xpos = np.arange(3); wdt = .36
ax.bar(xpos - wdt/2, [mp["sign"], mp["rep"], mp["r"]], wdt, color=BLUE, label="pipeline (same donors)")
ax.bar(xpos + wdt/2, [mc["sign"], mc["rep"], mc["r"]], wdt, color=RED, label="cohort (same pipeline)")
for i, (a1, a2) in enumerate(zip([mp["sign"], mp["rep"], mp["r"]], [mc["sign"], mc["rep"], mc["r"]])):
    ax.text(i - wdt/2, a1 + .02, f"{a1:.2f}", ha="center", fontsize=8)
    ax.text(i + wdt/2, a2 + .02, f"{a2:.2f}", ha="center", fontsize=8)
ax.set_xticks(xpos); ax.set_xticklabels(labels, fontsize=9); ax.set_ylim(0, 1.15)
ax.set_ylabel("value"); ax.legend(fontsize=8)
ax.set_title("What breaks replication: the cohort, not the pipeline", fontsize=11)
fig.tight_layout(); fig.savefig(FIG / "fig2_replication_metrics.png", dpi=160); plt.close(fig)

# ---- Fig 3: effect-size distributions, incl. histology-matched ----
hm = pd.read_csv(ROOT / "results" / "histology_matched.tsv", sep="\t")
order = ["TCGA (all histologies)", "TCGA (shared histologies only)",
         "nonTCGA (all)", "nonTCGA (shared histologies only)", "nonTCGA (solid only)"]
data = [hm[hm.setting == o]["beta"].abs().dropna().values for o in order]
fig, ax = plt.subplots(figsize=(8.4, 3.8))
bp = ax.boxplot(data, patch_artist=True, widths=.6, showfliers=False)
for i, box in enumerate(bp["boxes"]):
    box.set_facecolor(BLUE if order[i].startswith("TCGA") else RED); box.set_alpha(.7)
ax.set_xticklabels([o.replace(" (", "\n(") for o in order], fontsize=8)
ax.set_ylabel("|beta| (effect size)")
ax.set_title("Effect sizes shrink in non-TCGA cohorts even within matched histologies", fontsize=11)
for i, v in enumerate(data):
    ax.text(i + 1, np.median(v) + .004, f"{np.median(v):.3f}", ha="center", fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "fig3_effect_sizes.png", dpi=160); plt.close(fig)

print("wrote", *[p.name for p in sorted(FIG.glob('*.png'))])
