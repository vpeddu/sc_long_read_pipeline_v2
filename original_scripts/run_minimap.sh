#!/bin/bash

RESOURCE_DIR='/mnt/ix1/Resources'
MINIMAP_EXE="${RESOURCE_DIR}/tools/minimap2/v2.24/minimap2"
SAMTOOLS_EXE="${RESOURCE_DIR}/tools/samtools/v1.19.2/bin/samtools"
LR_SCRIPT_DIR="${RESOURCE_DIR}/LabSoftware/software/long_reads"

#If doing isoform analysis, provide bed12 annotations for gene transcripts for better alignment around splice junctions
REF_DIR='/mnt/ix1/Projects/M102_241107_Multiome/00_resources'
JUNC_BED="${REF_DIR}/genes.junc.bed"

fastq_gz=$(ls ../A0?_flexiplex/*.bc.fastq.gz)
sample=$(basename $fastq_gz .bc.fastq.gz)
echo "Input: $fastq_gz Sample: $sample"

# Minimap2 alignment
$MINIMAP_EXE -ax splice -k14 -t 10 --secondary=no --junc-bed $JUNC_BED ${REF_DIR}/genome.fa $fastq_gz | \
    $SAMTOOLS_EXE view -Sbh -F 2048 - > ${sample}.primary.bam

#Extract bc/UMI from read name and add as tags to bam file (CB, UB)
#This is needed for htSeq (de-dup'd bam), and for Vartrix (un-dedup'd bam)
python ${LR_SCRIPT_DIR}/bam_tag_bc_umi.py ${sample}.primary.bam

$SAMTOOLS_EXE sort ${sample}.primary.bc_tag.bam >${sample}.bc_tag.st.bam
$SAMTOOLS_EXE index ${sample}.bc_tag.st.bam

$SAMTOOLS_EXE flagstat ${sample}.bc_tag.st.bam >${sample}.flagstat.txt

#Cleanup:  Remove original aligned bam if tagged bam size > original size
#  Sorted bam will be smaller than unsorted, so don't test size of that, but
#  assume it is ok if it exists and is non-zero => can delete intermediate bams
if [ -s "${sample}.primary.bam" ] && [ -s "${sample}.primary.bc_tag.bam" ] && \
   [ -s "${sample}.bc_tag.st.bam" ]; then
  f1_bytes=$(stat -c%s "${sample}.primary.bam")
  f2_bytes=$(stat -c%s "${sample}.primary.bc_tag.bam")
  if [[ "$f2_bytes" -gt "$f1_bytes" ]]; then
    rm ${sample}.primary.bam
    rm ${sample}.primary.bc_tag.bam
  fi
fi
