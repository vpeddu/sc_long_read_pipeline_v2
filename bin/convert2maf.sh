#!/bin/bash
# This script runs within a Singularity container where all tools are available
# Activate vep conda environment
source /opt/miniconda3/bin/activate vep

REF_DIR=/mnt/ix1/Projects/M102_241107_Multiome/00_resources
VEP_PATH="$(which vep)"
VEP_DATA_DIR="${CONDA_PREFIX}/share/ensembl-vep-data/111_GRCh38"
VCF2MAF=/mnt/ix1/Resources/tools/vcf2maf/210423-754d68a/vcf2maf.pl

#Assume running from directory nomenclature:  /some/path/<sample>/Xnn_analysis, can derive sample as follows:
sample=$(basename $(dirname $PWD))

input_vcf=$(ls ${sample}.pass.vcf)
output_maf=${input_vcf%vcf*}maf

perl $VCF2MAF --input-vcf $input_vcf --output-maf $output_maf --vep-overwrite --tumor-id SAMPLE \
   --ref-fasta ${REF_DIR}/genome.fa --ncbi-build GRCh38 --vep-path $(dirname $VEP_PATH) \
   --vep-data $VEP_DATA_DIR

sed '1d' $output_maf | cut -d$'\t' -f 9 | sort | uniq -c >maf_pass.cts.txt

#Cleanup
rm ${sample}.pass*.vcf