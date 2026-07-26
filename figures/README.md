# Figures

## Schematics (BioRender)

Generated with BioRender AI. Prompted for Times New Roman at 14 pt throughout;
both render in a Times-style serif face.

| file | content |
|---|---|
| `fig0_study_design.png` | Study design: one shared genomic table feeding three arms (S1 TCGA/Pipeline A n=716, S2 same 716 donors/Pipeline B, S3 independent non-TCGA/Pipeline B n=423), with braces marking the PIPELINE and COHORT contrasts, and the 6 exposures × 10 immune readouts = 60 associations box. |
| *(pending export)* `fig5_summary.png` | Results summary: "Software is innocent; the cohort is the culprit". Left panel software-vs-cohort (94.8% agreement vs 8% replication); right panel the self-replication control bar chart (62% / 50% / 8%); bottom banner "A replication rate is meaningless without a self-replication baseline". |

**The summary schematic lives in BioRender and still needs a manual export.**
The MCP integration only returns a downloadable link at the moment the figure is
first saved, and that response is too large to receive intact, so it could not be
pulled down automatically. Open it and export as PNG into this folder as
`fig5_summary.png`:

<https://app.biorender.com/illustrations/7739648551188f0bed8e8f33?slideId=e4e76e54-470e-4574-8b51-f8dd2639a75b>

An earlier draft of the same panel is at
<https://app.biorender.com/illustrations/d091c9d618696e74a5ce552e?slideId=38357211-b188-7c00-2f64-7bebcf1bb489>
— it is superseded: its bar chart mixed the composition-controlled arms (54% / 0%,
`src/11`) with the prevalence-matched arm (50%, `src/10`), which are two different
experiments. The current version uses one consistent set (all three bars from
`src/10_self_replication.py`, n=340 per group).

Numbers in the summary schematic trace to:
- 94.8% shared variance → `results/pipeline_robustness.md`
- 62% / 50% / 8% → `results/self_replication.md`
- composition-controlled check (54% vs 0%) → `results/verify_composition_controlled.md`

## Analysis figures (matplotlib)

| file | source |
|---|---|
| `fig1_pipeline_vs_cohort.png` | `src/05_figures.py` — key scatter, β under pipeline swap vs cohort swap |
| `fig2_replication_metrics.png` | `src/05_figures.py` |
| `fig3_effect_sizes.png` | `src/05_figures.py` |
| `fig4_self_vs_cross_replication.png` | `src/09_metric_paradox_figure.py` |
