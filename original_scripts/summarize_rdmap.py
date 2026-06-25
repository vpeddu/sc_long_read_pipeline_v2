import argparse, sys, os
import csv, numpy as np, pandas as pd, polars as pl
import scipy.io as sio
from scipy import sparse
from gtfparse import read_gtf

#After standard flair processing, possible transcript/isoform id formats:
#1/ ENST00000263212_ENSG00000100034
#2/ ENST00000263212-0_ENSG00000100034
#3/ CCGTGGAAGTCCGGTC_GAGCTCCTTT#233c5c35-7b9f-4900-a111-65d36214a8f4_+1of1-0_ENSG00000100038
#4/ AGGTCCGAGCTGATAA_CTTACCTGCA#0f7a20c5-1814-4911-80f9-ad99c7dfba1b_+1of1-0_chr19:32629000

#After gffcompare/sqanti, possible formats are:
#5/ ENST00000361851_MT-ATP8
#6/ TCONS_00061076_WASHC1, or TCONS_00021342_antisense_WASHC1
#7/ GCCTAATAGGTGAAGC_CCTTGTAAAGGT#b82856ba-73eb-4214-a415-e203ea3e78c0_+1of1-0_FLNA
#8/ CATGCATAGTTAGAGG_CTAAAATTTTTT#1f578e95-d76a-4324-ba9e-111e6179eb57_-1of1-0_chrX:154366000

#Get gene symbols from features.tsv file if have Ensembl gene to convert
#Get transcript start/end from flair gtf file for novel transcripts (if --rename specified)
#Convert transcript/gene names as follows
#1/ transcript: ENST00000263212, gene: PPM1F
#2/ transcript: ENST00000263212-0, gene: PPM1F
#3/ transcript: TOP3B_21962592_20198 (where 21962592 is start coord from flair gtf file, and 20198
#                  is reference span [ie end - start +1])
#4/ ignore, not a known gene
#5/ transcript: ENST00000361851, gene: MT-ATP8
#6/ transcript: TCONS_00061076, gene: WASHC1
#7/ transcript: FLNA_<start>_<span>  (same as 3/)
#8/ ignore, not a known gene

csv.field_size_limit(sys.maxsize)

def parse_commandline():
  default_fn = 'flair.collapse.combined.isoform.read.map.txt'
  default_gtf = 'flair.collapse.isoforms.gtf'
  parser=argparse.ArgumentParser()
  parser.add_argument('--sample', '-s', help='sample name (prefix for output files)', type=str, required=True)
  parser.add_argument('--read_map', '-r', help='isoform read map file', type=str, default=default_fn, required=False)
  parser.add_argument('--gtf', '-g', help='flair gtf file', type=str, default=default_gtf, required=False)
  parser.add_argument('--rename', '-n', help='rename novel isoforms', action='store_true')
  parser.add_argument('--prepend_gene', '-p', help='prepend gene name to transcript id', action='store_true')
  parser.add_argument('--matrix', '-m', help='create sparse matrix output', action='store_true')
  args=parser.parse_args()
  print(args, file=sys.stderr)
  return args

args = parse_commandline()
sample = args.sample

fn = args.read_map
gtf_fn = args.gtf
features_fn = '/mnt/ix1/Projects/M102_241107_Multiome/00_resources/features_gex.tsv'

features_df = pd.read_csv(features_fn, sep='\t', header=None)
features_df.columns = ['ensembl_id', 'gene_symbol', 'assay']
#There are 10 gene names which do not map uniquely to a single ensembl id, so
#  fix these to use ensembl id as gene name, instead of gene symbol
#  TMSB15B,TBCE,MATR3,LINC01505,LINC01238,HSPA14,GOLGA8M,GGT1,CYB561D2,ARMCX5-GPRASP2
#Eg below:
#ENSG00000158427 TMSB15B Gene Expression
#ENSG00000269226 TMSB15B Gene Expression
#ENSG00000285053 TBCE    Gene Expression
#ENSG00000284770 TBCE    Gene Expression
#ENSG00000280987 MATR3   Gene Expression
#ENSG00000015479 MATR3   Gene Expression
non_unique_gene = features_df['gene_symbol'].duplicated(keep=False)
features_df.loc[non_unique_gene, 'gene_symbol'] = features_df.loc[non_unique_gene, 'ensembl_id']

gtf_df = read_gtf(gtf_fn)
transcript_df = gtf_df.filter(pl.col("feature") == "transcript")

def osuffix(args):
  args_o = ''
  if args.rename: args_o = 'n'
  if args.prepend_gene: args_o = args_o + 'p'
  if args_o == '':
    csv_o = 'flair_id'
  else:
    csv_o = "".join(['txmod_', args_o])
  return(csv_o)

def parse_flair_isoform(isoform, features_df):
  #Split isoform by "_".  Last piece will be gene either: ENSGxxxx or gene symbol (eg TP53), or chr:coord
  #If long isoform >72 chars => novel isoform, only process if associated with known gene
  gene_name = 'NA'
  transcript_id = 'NA'
  
  iso_parts = isoform.split("_")
  if len(iso_parts) > 2 and iso_parts[-2] == 'antisense':
    gene_id = '_'.join(iso_parts[-2:])
    transcript_id = "_".join(iso_parts[0:-2])
  else:
    gene_id   = iso_parts[-1]
    transcript_id = "_".join(iso_parts[0:-1])
  
  if len(isoform) > 72:  #novel isoform formatted as BC_UMI#ONTreadname..
    iso_type = 'NA' if gene_id[0:3] == 'chr' else 'novel'
  else:
    iso_type = 'ref'   #for this purpose can consider TCONS_nnnnn as ref as well as ENSTxxx
  
  if gene_id[0:4] == 'ENSG':
    gene_name = get_gene_name(gene_id, features_df)
  else:
    gene_name = gene_id

  return(iso_type, transcript_id, gene_name)
  
def get_gene_name(ensembl_id, features_df):
  gene_row = features_df[features_df["ensembl_id"] == ensembl_id]
  #This does not always return a valid row.  Novel isoforms may not be associated with
  #  an annotated gene, in which case the gene name is <chr>:<approx_start> where approx
  #  start is in steps of 10000bp.  eg. chr4:7451000
  #Update:  Since filtering now done prior to calling this function, gene symbol should
  #         always be returned successfully
  if gene_row.shape[0] == 0:
    return(ensembl_id)
  else:
    return(gene_row["gene_symbol"].item()) 

def novel_isoform_coords(isoform, transcript_df):
  tx_row = transcript_df.filter(pl.col("transcript_id") == isoform_id)
  tx_length = tx_row[0, "end"] - tx_row[0, "start"] + 1
  tx_coords = "_".join([str(tx_row[0, "start"]), str(tx_length)])
  return(tx_coords)

isoforms = []
barcodes = []
isoform_idx_cts = []

o_suff = osuffix(args)
o_csv = ".".join(['isoform_cells', o_suff, 'csv'])
w_csv = csv.writer(open(o_csv, "w"))

#Read flair read.map file which has one line per isoform, with comma-delimited list of read names
#  associated with that isoform
with open(fn, "r") as f:
    r_csv = csv.reader(f, delimiter="\t")
    for i, row in enumerate(r_csv):
        iso_type, isoform_id, gene_name = parse_flair_isoform(row[0], features_df)
        #isoform = reformat_isoform(row[0], transcript_df, features_df)
        if iso_type == 'NA':
          continue

        if iso_type == 'novel' and args.rename:  
          #print("Novel isoform", isoform_id)       
          tx_coords = novel_isoform_coords(isoform_id, transcript_df)
          isoform = "_".join([gene_name, tx_coords]) 
        #Will we ever want to pre-pend gene name to novel transcripts if not renaming??
        #Code below assumes no
        elif iso_type == 'ref' and args.prepend_gene:
          isoform = "_".join([gene_name, isoform_id])
        else:
          isoform = isoform_id
        isoforms.append(isoform)

        isoform_idx = len(isoforms)-1
        row_barcodes = []
        row_isoforms = []

        for rdname in row[1].split(','):
          row_barcodes.append(rdname[0:16])
        barcode, isoform_ct = np.unique(row_barcodes, return_counts=True)  
        for j, cell_bc in enumerate(barcode):
          try:
            bc_idx = barcodes.index(cell_bc)
          except:
            barcodes.append(cell_bc)
            bc_idx = len(barcodes)-1

          row_isoforms.append([isoform, gene_name, cell_bc, isoform_ct[j]])
          isoform_idx_cts.append([isoform_idx, bc_idx, isoform_ct[j]])
        
        w_csv.writerows(row_isoforms)

if args.matrix:
  sc_iso = [arr[0] for arr in isoform_idx_cts]
  sc_bc = [arr[1] for arr in isoform_idx_cts]
  sc_ct = [arr[2] for arr in isoform_idx_cts]

  matrix = sparse.csr_matrix((sc_ct, (sc_iso, sc_bc)))
  sio.mmwrite(sample + '.matrix.mtx', matrix)
  #barcodes = [barcode + '-1' for barcode in barcodes] #Do this in subsequent step (when needed for Seurat load)
  np.savetxt(sample + '.barcodes.tsv', barcodes, fmt="%s")
  np.savetxt(sample + '.features.tsv', isoforms, fmt="%s")
