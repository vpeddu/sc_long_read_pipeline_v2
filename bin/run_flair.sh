#!/bin/bash
# This script runs within a Singularity container where all tools are available
# Activate flair conda environment
source /opt/miniconda3/bin/activate flair

REF_DIR=/mnt/ix1/Projects/M102_241107_Multiome/00_resources

#Assume running from directory nomenclature:  /some/path/<sample>/Xnn_analysis, can derive sample as follows:
sample=$(basename $(dirname $PWD))

#DeDup directory should be A03 for standard pipeline, and A02 for merged bams
ddup_dir=$(ls -d ../A0?_deDup)

bam2Bed12 -i ${ddup_dir}/${sample}.umi_dd.bam >${sample}.umi_dd.bed

flair correct --threads 12 -q ${sample}.umi_dd.bed --gtf ${REF_DIR}/genes.gtf --genome ${REF_DIR}/genome.fa \
                --nvrna --output flair.filtered

#SG 12/5 adding quality parameter (min mapping quality of 10 to assign read to isoform)
flair collapse --threads 16 --temp_dir /mnt/ix2/TEMP -q flair.filtered_all_corrected.bed \
                 --reads ${ddup_dir}/${sample}.umi_dd.fastq.gz \
                 --genome ${REF_DIR}/genome.fa --gtf ${REF_DIR}/genes.gtf --annotation_reliant generate \
                 --generate_map --trust_ends --isoformtss --no_gtf_end_adjustment --check_splice --quality 10
