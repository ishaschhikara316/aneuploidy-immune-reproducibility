"""
01 — Build the unified analysis matrix for the cross-cohort reproducibility audit.

Produces one row per donor per RNA SOURCE, so the same donor can appear twice (once scored
from TCGA-Xena RNA, once from PCAWG RNA). That is what makes the two-arm design possible:

  ARM A (pipeline effect): TCGA donors, Xena RNA  vs  PCAWG RNA   [same donors, same exposures]
  ARM B (cohort effect):   PCAWG RNA, TCGA donors vs  non-TCGA donors  [same pipeline]

Exposures (WGD/ploidy/chromothripsis/PGA) and purity come from ONE shared donor table, so they
are defined identically everywhere; only the RNA source or the donor set varies.

Outputs data/analysis_matrix.parquet
"""
import gzip, warnings, numpy as np, pandas as pd, yaml
from pathlib import Path
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((ROOT / "config" / "panels.yaml").read_text())
PANELS = CFG["panels"]
WANT = sorted({g for v in PANELS.values() for g in v})

XENA = ROOT / "data" / "tcga" / "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"
PCAWG = ROOT / "data" / "pcawg" / "tophat_star_fpkm_uq.v2_aliquot_gl.tsv.gz"
SHEET = ROOT / "data" / "pcawg" / "pcawg_sample_sheet.tsv"
ATLAS = Path.home() / "isha-ssd" / "pd-sex-dimorphism" / "data" / "raw" / "pd_atlas.h5ad"

to_patient = lambda b: "-".join(str(b).split("-")[:3])


def zpanel(expr, genes, group):
    """within-group z-score each gene, then average -> panel score (matches prior pipelines)."""
    genes = [g for g in genes if g in expr.columns]
    if not genes:
        return pd.Series(np.nan, index=expr.index)
    o = pd.DataFrame(index=expr.index)
    for g in genes:
        o[g] = expr.groupby(group)[g].transform(
            lambda v: (v - v.mean()) / v.std(ddof=0) if v.std(ddof=0) else 0.0)
    return o.mean(axis=1)


# ---------------- exposures (shared, identical definitions) ----------------
chromo = pd.read_csv(ROOT / "data" / "chromo_per_donor.csv")
chromo["WGD"] = (chromo["ploidy"] > 2.5).astype(float)
chromo["chromo"] = chromo["chromo_high"].astype(float)
print(f"exposure table: {len(chromo)} donors "
      f"({chromo.barcode.notna().sum()} TCGA / {chromo.barcode.isna().sum()} non-TCGA)")

# ---------------- A. TCGA-Xena RNA ----------------
print("\nextracting panel genes from TCGA-Xena RNA ...")
keep_rows, header = [], None
with open(XENA) as fh:
    header = fh.readline().rstrip("\n").replace('"', "").split("\t")
    for line in fh:
        sym = line.split("\t", 1)[0].strip('"').split("|")[0]
        if sym in WANT:
            keep_rows.append(line.rstrip("\n").replace('"', "").split("\t"))
xe = pd.DataFrame(keep_rows, columns=header).set_index(header[0])
xe.index = [i.split("|")[0] for i in xe.index]
xe = xe[~xe.index.duplicated()].apply(pd.to_numeric, errors="coerce").T
xe = xe[[c.split("-")[3][:2].isdigit() and 1 <= int(c.split("-")[3][:2]) <= 9
         for c in xe.index]]                                   # tumour aliquots only
xe["barcode"] = [to_patient(i) for i in xe.index]
xe = xe.groupby("barcode").mean()
print(f"  Xena: {xe.shape[1]} panel genes x {xe.shape[0]} TCGA patients")

# ---------------- B. PCAWG RNA ----------------
print("extracting panel genes from PCAWG RNA ...")
import anndata as ad
v = ad.read_h5ad(ATLAS, backed="r").var
sym2ens = {}
for ens, sym in zip(v.index, v["feature_name"].astype(str)):
    if sym in WANT:
        sym2ens.setdefault(sym, ens.split(".")[0])
ens2sym = {e: s for s, e in sym2ens.items()}
print(f"  symbol->Ensembl resolved for {len(sym2ens)}/{len(WANT)} panel genes")

rows, hdr = [], None
with gzip.open(PCAWG, "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        gid = line.split("\t", 1)[0].split(".")[0]
        if gid in ens2sym:
            rows.append(line.rstrip("\n").split("\t"))
pe = pd.DataFrame(rows, columns=hdr).set_index("feature")
pe.index = [ens2sym[i.split(".")[0]] for i in pe.index]
pe = pe[~pe.index.duplicated()].astype(float)
pe = np.log2(pe + 1).T                                          # FPKM-UQ -> log space

sheet = pd.read_csv(SHEET, sep="\t", low_memory=False)
rna = sheet[sheet.library_strategy.eq("RNA-Seq")]
spec = rna.dcc_specimen_type.astype(str).str.lower()
rna = rna[spec.str.contains("tumour|tumor|primary|metasta|recurr", na=False)]
amap = rna.drop_duplicates("aliquot_id").set_index("aliquot_id")["donor_unique_id"].to_dict()
pe["donor_unique_id"] = [amap.get(a) for a in pe.index]
pe = pe.dropna(subset=["donor_unique_id"]).groupby("donor_unique_id").mean()
print(f"  PCAWG: {pe.shape[1]} panel genes x {pe.shape[0]} donors")

# ---------------- assemble both sources ----------------
frames = []
# Xena arm: join on TCGA barcode
a = chromo.dropna(subset=["barcode"]).merge(xe, on="barcode", how="inner")
a["rna_source"], a["cohort"] = "TCGA_Xena", "TCGA"
frames.append(a)
# PCAWG arm: join on donor_unique_id (both TCGA and non-TCGA donors)
b = chromo.merge(pe, left_on="donor_unique_id", right_index=True, how="inner")
b["rna_source"] = "PCAWG"
b["cohort"] = np.where(b["barcode"].notna(), "TCGA", "nonTCGA")
frames.append(b)

d = pd.concat(frames, ignore_index=True)

# panel scores, z-scored WITHIN histology AND within rna_source (so sources are comparable)
d["_grp"] = d["rna_source"] + "|" + d["histo"].astype(str)
for name, genes in PANELS.items():
    d[name] = zpanel(d, genes, "_grp")
d = d.drop(columns=["_grp"])

gene_cols = [g for g in WANT if g in d.columns]
d = d.drop(columns=gene_cols)                                    # keep scores, drop raw genes
d.to_parquet(ROOT / "data" / "analysis_matrix.parquet")

print("\n=== analysis matrix ===")
print(d.groupby(["rna_source", "cohort"]).size().to_string())
print(f"\npanel scores: {list(PANELS)}")
print(f"wrote data/analysis_matrix.parquet  ({d.shape[0]} rows x {d.shape[1]} cols)")
