import argparse, sys, os
import csv, numpy as np, pandas as pd
import scipy.io as sio
from scipy import sparse

#Possible formats in x-reference files (after post-flair processing):
#1/ transcript_id: CCGTGGAAGTCCGGTC_GAGCTCCTTT#233c5c35-7b9f-4900-a111-65d36214a8f4_+1of1-0, 
#   tx_name: <disease><sample>_nnnnnn
#2/ transcript_id: ENST00000471248  (sometimes has -n suffix); tx_name = transcript_id

#gene_name (symbol) and gene_id (ensembl) are both in the x-reference file
#transcripts within intergenic regions, or associated with non-canonical chromosomes have
#  been removed from x-reference

#Possible transcript/isoform id formats in read_map file:
#1/ ENST00000349431-0_ENSG00000160087 (always -n suffix after ENST transcript)
#2/ CCGTGGAAGTCCGGTC_GAGCTCCTTT#233c5c35-7b9f-4900-a111-65d36214a8f4_+1of1-0_ENSG00000100038
#3/ AGGTCCGAGCTGATAA_CTTACCTGCA#0f7a20c5-1814-4911-80f9-ad99c7dfba1b_+1of1-0_chr19:32629000

csv.field_size_limit(sys.maxsize)

def parse_commandline():
  default_fn = 'flair.collapse.combined.isoform.read.map.txt'
  default_xref = 'transcript_xref.tsv'
  parser=argparse.ArgumentParser()
  parser.add_argument('--sample', '-s', help='sample name (prefix for output files)', type=str, required=True)
  parser.add_argument('--read_map', '-r', help='isoform read map file', type=str, default=default_fn, required=False)
  parser.add_argument('--xref', '-x', help='gtf tx cross-reference', type=str, default=default_xref, required=False)
  parser.add_argument('--prepend_gene', '-p', help='prepend gene name to transcript id', action='store_true')
  parser.add_argument('--matrix', '-m', help='create sparse matrix output', action='store_true')
  parser.add_argument('--keep_intergenic', help='keep novel isoforms whose gene tag is not an ENSG id (novel/intergenic loci) instead of dropping them', action='store_true')
  args=parser.parse_args()
  print(args, file=sys.stderr)
  return args

args = parse_commandline()
sample = args.sample

flair_fn = args.read_map
xref_fn = args.xref
w_csv = csv.writer(open('isoform_cells.csv', "w"))

def parse_flair_isoform(isoform):
  #Split isoform by "_".  If only two pieces then should be ENSTxxx_ENSGxxx
  #  If > 2 pieces (will always be 4?), then last piece could be either ENSGxxx or chr:coord
  iso_type = 'NA'
  ensembl_gene = 'NA'
  transcript_id = 'NA'
  
  iso_parts = isoform.split("_")
  if len(iso_parts) == 2:  #Ensembl transcript and gene
    iso_type = 'ref'
    ensembl_gene = iso_parts[1]
    transcript_id = iso_parts[0] 
 
  elif len(iso_parts) > 2:
    #Novel transcript associated with an Ensembl gene, or (if --keep_intergenic)
    #  a non-ENSG locus tag (eg chr19:32629000) for a novel/intergenic gene
    if iso_parts[-1][0:4] == 'ENSG' or args.keep_intergenic:
      iso_type = 'novel'
      ensembl_gene = iso_parts[-1]
      transcript_id = "_".join(iso_parts[0:-1]) #Flair novel transcript id (~= ONT read name)
  return(iso_type, transcript_id, ensembl_gene)

isoforms = []
barcodes = []
isoform_idx_cts = []

#Read transcript xref file
#Columns: ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 
#          'gene_id', 'transcript_id', 'gene_name', 'novel_tx', 'novel_tx_gene']
xref_df = pd.read_csv(xref_fn, sep='\t')
xref_df = xref_df[['transcript_id', 'gene_name', 'novel_tx', 'novel_tx_gene']]
tx_xref = xref_df.set_index('transcript_id').T.to_dict('list')

#Read flair read.map file which has one line per isoform, with comma-delimited list of read names
#  associated with that isoform
with open(flair_fn, "r") as f:
    r_csv = csv.reader(f, delimiter="\t")
    for i, row in enumerate(r_csv):
        iso_type, isoform_id, ensembl_gene = parse_flair_isoform(row[0])
        if iso_type == 'NA':
          continue

        #Transcript ids not in xreference file are either non-canonical chromosomes
        #  or in intergenic regions => ignore
        try:
          gene_name = tx_xref[isoform_id][0]
        except KeyError:
          continue

         #For novel isoforms, get new isoform name from xref  
        tx_name = tx_xref[isoform_id][1] if iso_type == 'novel' else isoform_id 
        isoform = "_".join([gene_name, tx_name]) if args.prepend_gene else tx_name       

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
