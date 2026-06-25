#!/bin/bash
source activate long_reads

#Assume running from directory nomenclature:  /some/path/<sample>/Xnn_analysis, can derive sample as follows:
sample=$(basename $(dirname $PWD))

REF_DIR=/mnt/ix1/Projects/M102_241107_Multiome/00_resources
GENOME_REF=${REF_DIR}/genome.fa

mkdir -p snv_debug

longshot -F --bam ${sample}.umi_dd.bam --ref $GENOME_REF --out ${sample}.vcf \
     --strand_bias_pvalue_cutoff 0.0001 --density_params 10:100:50 \
     --no_haps --variant_debug_dir snv_debug

##### Additional longshot parameters
##### --max_alignment  Use max scoring alignment algorithm rather than pair HMM forward algorithm.
##### --no_haps        Don't call HapCUT2 to phase variants.
##### --strand_bias_pvalue_cutoff (default 0.01)
##### --density_params <string>  Parameters to flag a variant as part of a "dense cluster". 
#####                            Format <n>:<l>:<gq>. If there are at least n variants within l base pairs with
#####                                genotype quality >=gq, flag variants as "dn" [default: 10:500:50]
 
#Filter to include PASS variants only
grep '^#' ${sample}.vcf > ${sample}.pass.vcf
awk '$7 == "PASS"' ${sample}.vcf >> ${sample}.pass.vcf
