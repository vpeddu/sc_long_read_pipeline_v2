#!/bin/bash
# This script runs within a Singularity container where all tools are available
# Activate long_reads conda environment where featureCounts is available
source /opt/miniconda3/bin/activate long_reads

REF_DIR=/mnt/ix1/Projects/M102_241107_Multiome/00_resources

#Assume running from directory nomenclature:  /some/path/<sample>/Xnn_analysis, can derive sample as follows:
sample=$(basename $(dirname $PWD))

featureCounts \
  -a ${REF_DIR}/genes.gtf \
  -T 12 -L -t exon -g gene_name -o exon_counts.txt \
  ${sample}.umi_dd.bam
