"""09 — the metric-paradox figure: same comparison, opposite conclusions."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
RED, BLUE, GREY = "#E03131", "#4C6EF5", "#adb5bd"

metrics = ["associations still\nsignificant (FDR<0.05)", "sign\nagreement",
           "effect ratio\n(non-TCGA / TCGA)", "calibrated discrepancy\nvs re-splitting TCGA"]
vals = [0.0, 1.00, 0.54, 1.0/1.17]     # all expressed as "fraction of discovery retained"
labels = ["0/60\n(0%)", "100%", "0.54", "0.85\n(1.17x more\ndiscrepant)"]
cols = [RED, BLUE, "#f59f00", BLUE]

fig, ax = plt.subplots(figsize=(8.6, 4.0))
bars = ax.bar(metrics, vals, color=cols, alpha=.85)
for b, l in zip(bars, labels):
    ax.text(b.get_x()+b.get_width()/2, b.get_height()+.03, l, ha="center", fontsize=9)
ax.axhline(1.0, color=GREY, ls="--", lw=1)
ax.text(3.45, 1.02, "perfect\nreplication", fontsize=7, color="#666", ha="right")
ax.set_ylim(0, 1.25); ax.set_ylabel("apparent replication\n(1 = fully reproduced)")
ax.set_title("The same TCGA → non-TCGA comparison, four metrics, opposite conclusions",
             fontsize=12)
ax.text(0, .12, "'catastrophic'", ha="center", fontsize=8, color=RED, style="italic")
ax.text(3, .72, "'barely matters'", ha="center", fontsize=8, color=BLUE, style="italic")
fig.tight_layout(); fig.savefig(ROOT/"figures"/"fig4_metric_paradox.png", dpi=160)
print("wrote fig4_metric_paradox.png")
