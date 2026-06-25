#!/bin/bash

source activate /venvs2/anaconda3/envs/lrseq

CODE_DIR=/mnt/ix1/Projects/M102_241107_Multiome/00_code/LongReads
declare -A disease_code=( [Leukemia]="LEU" [CRC]="CRC" [GIM]="GIM" [Gastric]="GST" \
                          [Organoids]="GOO" )

#MUST be running from directory nomenclature:  /some/path/<sample>/Xnn_analysis for derivation
#  of sample and disease-type to be correct
#If not running from this directory structure, provide sample and tprefix directly
#  eg: sample SU968_2926A ..this will be the prefix for matrix output files
#  eg: tprefix LEU2926A ..this will be the prefix for novel transcript names
sample=$(basename $(dirname $PWD))
sample_bc="${sample#*_}"
sample_bc=${sample_bc/_WES_mrg/WES} #Shorten name to just be <sample>WES if _WES_mrg sample suffix

disease=$(echo "$PWD" | awk -F/ '{print $(NF-3)}')
tprefix="${disease_code[$disease]}${sample_bc}"
python ${CODE_DIR}/rename_gtf_txid.py --tx_prefix $tprefix 

python ${CODE_DIR}/rdmap2counts.py --sample $sample --matrix
gzip ${sample}.*.tsv ${sample}.matrix.mtx
