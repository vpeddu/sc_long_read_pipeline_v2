##NOTE:  MUST run on server with AVX2 support in order to use polars
##         otherwise will abort with illegal instruction set error
##       Use fugu, iwashi, suzuki
##Run in lrseq venv (conda activate lrseq)

import argparse, sys, os
import csv, numpy as np, pandas as pd, polars as pl
import scipy.io as sio
from scipy import sparse
from gtfparse import read_gtf

#Later might want parameter to select format for novel transcripts
#Eg:
# 1: <prefix>_<n> where n is a sequential number uniquely assigned per novel transcript
# 2: <gene>_novel<n> where n is a sequential number assigned per gene/novel transcript
# 3: N<start>_<length> format where <start> is start coord from flair gtf, and <length> is ref span
#Currently formats 1 and 2 are output to x-reference file, and format 1 is used for .gtf file.
#novel_tx_fmt = 1

#Possible gene formats:
#ENSG00000100038
#19:32629000

#Possible transcript_id formats:
#ENST00000263212
#ENST00000263212-0
#CCGTGGAAGTCCGGTC_GAGCTCCTTT#233c5c35-7b9f-4900-a111-65d36214a8f4_+1of1

#Get gene symbols from features.tsv file

csv.field_size_limit(sys.maxsize)

def parse_commandline():
  default_gtf = 'flair.collapse.isoforms.gtf'
  parser=argparse.ArgumentParser()
  parser.add_argument('--gtf', '-g', help='flair gtf file', type=str, default=default_gtf, required=False)
  parser.add_argument('--tx_prefix', '-x', help='prefix for novel transcript ids', type=str, required=True)
  parser.add_argument('--features', '-f', help='features.tsv file (ensembl_id, gene_name, assay columns)', type=str, required=True)
  args=parser.parse_args()
  print(args, file=sys.stderr)
  return args

#Not currently using this function, but keeping here for now in case in future want 
#  the option for format 3  
def novel_isoform_coords(isoform, transcript_df):
  tx_row = transcript_df.filter(pl.col("transcript_id") == isoform_id)
  tx_length = tx_row[0, "end"] - tx_row[0, "start"] + 1
  tx_coords = "_".join([str(tx_row[0, "start"]), str(tx_length)])
  return(tx_coords)  
  
def write_gtf(df, out_fn):
  #gtf_df.columns (flair)
  #['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'gene_id', 'transcript_id', 'exon_number']
  #gtf_df.columns (gffcompare)
  #['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'transcript_id', 'gene_id', 'gene_name',
  # 'xloc', 'cmp_ref', 'class_code', 'tss_id', 'exon_number', 'cmp_ref_gene', 'contained_in', 'ref_gene_id']
  with open(out_fn, 'w') as f:
    for row in df.rows():
      seqname, source, feature, start, end, score, strand, frame = row[0:8]
      #Create dictionary of attribute names and values for variable part of gtf rows
      attr_dict = dict(zip(df.columns[8:], row[8:]))
      
      #Hardcode score and frame as '.' as in original flair gtf file (otherwise are None and 0 respectively)
      gtf_cols = f"{seqname}\t{source}\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t"
      gtf_attr = ""
      if feature == 'exon':
        for attr in ['gene_id', 'transcript_id', 'exon_number']:
          gtf_attr = gtf_attr + f"{attr} \"{attr_dict[attr]}\"; "
      else:
        for attr, value in attr_dict.items():
          if value:
            gtf_attr = gtf_attr + f"{attr} \"{value}\"; "
      f.write(gtf_cols + gtf_attr.rstrip() + '\n')  

args = parse_commandline()
noveltx_prefix = args.tx_prefix + "_"

features_fn = args.features
features_pd = pd.read_csv(features_fn, sep='\t', header=None)
features_pd.columns = ['ensembl_id', 'gene_name', 'assay']

#There are 10 gene names which do not map uniquely to a single ensembl id, so
#  fix these to use ensembl id as gene name, instead of gene symbol
#  TMSB15B,TBCE,MATR3,LINC01505,LINC01238,HSPA14,GOLGA8M,GGT1,CYB561D2,ARMCX5-GPRASP2
#Eg below:
#ENSG00000158427 TMSB15B Gene Expression
#ENSG00000269226 TMSB15B Gene Expression
non_unique_gene = features_pd['gene_name'].duplicated(keep=False)
features_pd.loc[non_unique_gene, 'gene_name'] = features_pd.loc[non_unique_gene, 'ensembl_id']
#Convert to polars dataframe for later joining to gtf dataframe
features_df = pl.from_pandas(features_pd).drop('assay')

#Read gtf file and filter to contain only valid genes
gtf_fn = args.gtf
out_gtf = gtf_fn[0:-4] + '.txmod.gtf'
gtf_df = read_gtf(gtf_fn)
gtf_df = gtf_df.filter(pl.col("gene_id").str.starts_with('ENSG'))

#Additionally filter for canonical chromosomes (assume len(seqname) < 6 will do this)
gtf_df = gtf_df.with_columns(pl.col("seqname").cast(pl.String))
gtf_df = gtf_df.filter(pl.col("seqname").str.len_bytes() < 6)

#For all novel transcript rows, create novel transcript ids with two potential
#  nomenclature schemes:  1/ <prefix>_xxxxx, 2/ <gene>_novelxxx
#Plotting and SQANTI cannot handle the long novel flair transcript ids and/or the special characters
transcript_df = gtf_df.filter(pl.col("feature") == "transcript").drop('exon_number')
transcript_df = transcript_df.join(features_df, left_on='gene_id', right_on='ensembl_id', how='inner')

#Split dataframe into reference transcripts and novel transcripts
ref_transcript_df = transcript_df.filter(pl.col('transcript_id').str.starts_with("ENST"))
novel_transcript_df = transcript_df.filter(~pl.col('transcript_id').str.starts_with("ENST"))

#Modify novel transcript ids
#Note: older versions of polars use 'with_row_count' instead of 'with_row_index'
#Format 1: novel_tx (<prefix>_xxxxxx)
novel_transcript_df = novel_transcript_df.with_row_index('row_idx', offset=1).with_columns(
                        pl.col('row_idx').cast(pl.String).str.zfill(6).alias('tx_num')
                        ).with_columns(
                        (pl.lit(noveltx_prefix) + pl.col("tx_num")).alias("novel_tx")
                        )

#Format 2: novel_tx_gene (<gene>_novelxxx)
novel_transcript_df = novel_transcript_df.with_columns(
                         pl.col('gene_name').rank(method='ordinal').over('gene_name').cast(pl.String)
                         .str.zfill(3).alias('gene_idx')
                         ).with_columns(
                         pl.concat_str(pl.col("gene_name"), pl.lit("_novel"), pl.col("gene_idx"))
                        .alias("novel_tx_gene")
                        )

#Drop intermediate data columns and re-merge reference and novel transcripts
#Since additional columns in novel_transcripts_df, use how="diagonal"
#  which will fill missing data with null.  Otherwise use how="vertical"
novel_transcript_df = novel_transcript_df.drop(['row_idx', 'gene_idx', 'tx_num'])

transcripts_final = pl.concat([ref_transcript_df, novel_transcript_df], how="diagonal")
transcripts_final.write_csv('transcript_xref.tsv', separator='\t')

#Create x-reference file mapping .gtf transcript ids to new transcript names 
#Rename whatever field is going to be the transcript_id, to 'tx_name'
transcript_xref = transcripts_final.rename({"novel_tx": "tx_name"})
transcript_xref = transcript_xref.with_columns(
                    pl.col('tx_name').fill_null(pl.col('transcript_id'))
                    ).select(['transcript_id', 'tx_name', 'gene_name', 'novel_tx_gene'])

#Modify transcript/gene nomenclature for gtf transcript and exon rows for each transcript
gtf_final = gtf_df.join(transcript_xref, on='transcript_id', how='inner')
gtf_final = gtf_final.with_columns(
                       gene_id = pl.col('gene_name'), transcript_id = pl.col('tx_name')
                       ).drop(['tx_name', 'gene_name', 'novel_tx_gene'])
#gtf_final.head()
write_gtf(gtf_final, out_gtf)
