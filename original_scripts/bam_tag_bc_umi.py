#!/usr/bin/env python
import sys, os, re, pysam, csv, itertools

script_name = os.path.basename(__file__)
print("Running ", script_name)

if len(sys.argv) < 2:
  print("Usage: ", script_name, "<bam_file>")
  sys.exit(1)

bam_fn = sys.argv[1]
if bam_fn[-3:] == 'bam' and os.path.isfile(bam_fn) and os.access(bam_fn, os.R_OK):
  bam_input = pysam.AlignmentFile(bam_fn, 'rb')
else:
  print("Unable to open bam file for input:", bam_fn)
  sys.exit(1)
  
try:
  bam_outfn = bam_fn[0:-4] + '.bc_tag.bam'
  bam_output = pysam.AlignmentFile(bam_outfn, 'wb', template=bam_input)
except:
  print("Unable to open bam file for output: ", bam_outfn)
  sys.exit(1)

i=0    
for bamrd in bam_input.fetch(until_eof=True):
  i += 1
  barcode = bamrd.qname[0:16]
  umi = bamrd.qname[17:29]
  bamrd.tags = bamrd.tags + [('CB', barcode), ('UB', umi)]
  bam_output.write(bamrd)
   
bam_input.close()
bam_output.close()


