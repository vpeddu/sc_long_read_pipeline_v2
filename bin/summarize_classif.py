import pandas as pd
import re

def remove_tx_suffix(df, col):
  #Remove version suffix from column, if present
  df[col] = df[col].str.replace(r"\.\d+$", "", regex=True)
  return df

# -----------------------------
# Load main SQANTI file
# -----------------------------
df = pd.read_csv(
    "sqanti_RulesFilter_result_classification.txt",
    sep="\t",
    header=0,
    low_memory=False
)

# Ensure FL is numeric
df["FL"] = pd.to_numeric(df["FL"], errors="coerce")

# Not needed?
#df = remove_tx_suffix(df, 'isoform')

# -----------------------------
# Load transcript duplicates file
# -----------------------------
xref = pd.read_csv('tx_dups_xref.tsv', sep="\t", names=['transcript_1', 'gene_1', 'transcript_2', 'gene_2'])

# Not needed?
#xref = remove_tx_suffix(xref, 'transcript_1')
#xref = remove_tx_suffix(xref, 'transcript_2')

# Merge sqanti classification file, and tx-duplicates
df = df.merge(
    xref[["transcript_1", "transcript_2"]],
    left_on="isoform",
    right_on="transcript_1",
    how="left"
)

# transcript_2 is highest number Ensembl tx id, so assume that is 'best'
df["isoform"] = df["transcript_2"].combine_first(df["isoform"])
df = df.drop(columns=["transcript_1", "transcript_2"])

# -----------------------------
# Load gene/transcript annotation file 
# -----------------------------
annot = pd.read_csv("transcript_xref.tsv", sep="\t")

annot = annot[["transcript_id", "gene_id"]].drop_duplicates()

# Not needed?
#annot = remove_tx_suffix(annot, 'transcript_id')

# Map transcript -> Ensembl gene_id
df = df.merge(
    annot,
    left_on="isoform",
    right_on="transcript_id",
    how="left"
)

# Overwrite gene column with Ensembl gene_id
df["associated_gene"] = df["gene_id"].combine_first(df["associated_gene"])

df = df.drop(columns=["transcript_id", "gene_id"])

# -----------------------------
# Split ENST vs others
# -----------------------------
df_enst = df[df["isoform"].str.startswith("ENST", na=False)].copy()
df_other = df[~df["isoform"].str.startswith("ENST", na=False)].copy()

# Drop rows with invalid FL
#df_enst = df_enst.dropna(subset=["FL"])

# -----------------------------
# Grouping + duplicate replacement
# -----------------------------

# Group key will be Ensembl id without suffix
df_enst["group_key"] = df_enst["isoform"].str.replace(r"-\d+$", "", regex=True)

df_enst["suffix"] = (
    df_enst["isoform"]
    .str.extract(r"-(\d+)$")[0]
    .astype(float)
)

df_enst["has_suffix"] = df_enst["suffix"].notna()
df_enst["suffix"] = df_enst["suffix"].fillna(-1)

group_sum = df_enst.groupby("group_key")["FL"].sum()

# Sort: prefer canonical transcript (no suffix), then smallest suffix
df_sorted = df_enst.sort_values(
    by=["group_key", "FL", "has_suffix", "suffix"],
    ascending=[True, False, True, True]
)

# Pick best row per group, and replace the FL column with the summed value
best_rows = df_sorted.drop_duplicates("group_key", keep="first").copy()
best_rows["FL"] = best_rows["group_key"].map(group_sum)
best_rows = best_rows.drop(columns=["group_key", "suffix", "has_suffix"])

result = pd.concat([best_rows, df_other], ignore_index=True)

print("Total rows:", len(result))
#print(
#    "Non-ENSG gene IDs remaining:",
#    (~result["associated_gene"].astype(str).str.startswith("ENSG")).sum()
#)

result.to_csv("sqanti_classif_summ.txt", sep="\t", index=False)
