"""09 — key figure: TCGA replicates itself 62%, but only 8% in an independent cohort."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
BLUE, RED, GREY = "#4C6EF5", "#E03131", "#adb5bd"

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.2))

# --- panel A: the two causes, side by side ---
labels = ["change the\nSOFTWARE\n(same patients)", "change the\nCOHORT\n(same software)"]
vals = [100, 8]
bars = ax1.bar(labels, vals, color=[BLUE, RED], alpha=.88, width=.6)
for b, v in zip(bars, vals):
    ax1.text(b.get_x()+b.get_width()/2, v+3, f"{v}%", ha="center", fontsize=13, fontweight="bold")
ax1.set_ylim(0, 118); ax1.set_ylabel("associations that replicate (%)")
ax1.set_title("What actually breaks replication", fontsize=12)
ax1.axhline(100, color=GREY, ls="--", lw=1)

# --- panel B: the control — does TCGA replicate itself? ---
labels2 = ["TCGA →\nitself", "TCGA → itself\n(prevalence\nmatched)", "TCGA →\nnon-TCGA"]
vals2 = [62, 50, 8]
err = [[62-24, 50-3, 8-0], [100-62, 100-50, 38-8]]
bars2 = ax2.bar(labels2, vals2, color=[BLUE, "#f59f00", RED], alpha=.88, width=.6)
ax2.errorbar(range(3), vals2, yerr=err, fmt="none", ecolor="#333", capsize=4, lw=1.2)
for b, v in zip(bars2, vals2):
    ax2.text(b.get_x()+b.get_width()/2, v+5, f"{v}%", ha="center", fontsize=12, fontweight="bold")
ax2.set_ylim(0, 118); ax2.set_ylabel("replication rate (%)")
ax2.set_title("Control: a cohort replicating ITSELF", fontsize=12)
ax2.text(1, 108, "only 1% of self-draws are as poor\nas a typical cross-cohort draw",
         ha="center", fontsize=8, style="italic", color="#555")
fig.suptitle("Software is innocent; the cohort is the culprit", fontsize=13)
fig.tight_layout(); fig.savefig(ROOT/"figures"/"fig4_self_vs_cross_replication.png", dpi=160)
print("wrote fig4_self_vs_cross_replication.png")
