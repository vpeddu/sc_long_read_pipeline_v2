import pandas as pd

df = pd.read_csv("isoform_cells.csv", header=None,
                 names=["isoform", "gene", "barcode", "count"])

iso_counts = df.groupby("isoform", as_index=False)["count"].sum()
total = iso_counts["count"].sum()
iso_counts["ratio"] = iso_counts["count"] / total

iso_counts.columns = ['pbid','count_fl','norm_fl']
iso_counts.to_csv("isoform_counts.tsv", sep='\t', index=False)
