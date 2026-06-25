#!/bin/bash

FLEXIPLEX="/mnt/ix2/Experimental_tools/flexiplex/v1.01/flexiplex"

# Function to print help
print_usage()
{
	   echo "Usage: $(basename "$0") [-m]"; 
	   echo "Where -m multi-seq processing, shorter adapter sequence";
	   return
}
multiseq=0

# Parse command line options
OPTIND=1
while getopts "mh" OPT
do
  case "$OPT" in
    m) multiseq=1;;
    h) print_usage; exit 1;;
   \?) print_usage; exit 1;;
    :) echo "Option -$OPTARG requires an argument."; print_usage; exit 1;;
  esac
done

sample=$(basename $(dirname $PWD))
seq_type='3prm'

#If multi-seq demux sample, need to modify flanking adapter sequence to be shorter
if [[ $multiseq == 1 ]]; then
  flank_len=17
else
  flank_len=24
fi

if [[ ${sample: -4} == "WTX5" ]]; then
  seq_type='5prm'
fi

fastq_gz=$(ls ../${sample}*.fastq.gz)
bc_tsv=$(ls ../${sample}*.barcodes.tsv)
#sample=$(basename $fastq_gz .fastq.gz)
echo "Inputs: $fastq_gz $bc_tsv Sample: $sample"

adapter="CTACACGACGCTCTTCCGATCT"
# Parameter -d 10x3v3, is equivalent to:
# -x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ???????????? -x TTTTTTTTT -f 8 -e 2
# For MultiSeq de-mux, shorten flanking sequence, so modify to:
# -x CGACGCTCTTCCGATCT -b ???????????????? -u ???????????? -x TTTTTTTTT -f 8 -e 2
# Parameter -d 10x5v2, is equivalent to:
# -x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ?????????? -x TTTCTTATATGGG -f 8 -e 2

if [[ $flank_len -lt 24 ]]; then
  adapter=${adapter: -$flank_len}
fi

# Run flexiplex on fastq file
#gunzip -c $fastq_gz | $FLEXIPLEX -d 10x3v3 -p 20 -k $bc_tsv -n ${sample}.bc | \
#    gzip > "${sample}.bc.fastq.gz"

# Run flexiplex on fastq file
if [[ $seq_type == '3prm' ]]; then 
  gunzip -c $fastq_gz | $FLEXIPLEX -p 20 -x $adapter -b ???????????????? -u ???????????? \
                           -x TTTTTTTTT -f 8 -e 2 -k $bc_tsv -n ${sample}.bc | \
                        gzip > "${sample}.bc.fastq.gz"
else
  gunzip -c $fastq_gz | $FLEXIPLEX -p 20 -d 10x5v2 -k $bc_tsv -n ${sample}.bc | \
                        gzip > "${sample}.bc.fastq.gz"
fi

