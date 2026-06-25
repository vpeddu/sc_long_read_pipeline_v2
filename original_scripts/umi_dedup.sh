#!/bin/bash

source activate long_reads
RESOURCE_DIR='/mnt/ix1/Resources'
CODE_DIR='/mnt/ix1/Projects/M102_241107_Multiome/00_code/LongReads'

PICARD_DIR="${RESOURCE_DIR}/tools/picard-tools/v2.23.3"
SAMTOOLS="${RESOURCE_DIR}/tools/samtools/v1.19.2/samtools"
SEQTK="${RESOURCE_DIR}/tools/seqtk/v1.4/seqtk"

LEUKEMIA_LR=/mnt/ix1/Projects/M102_241107_Multiome/Leukemia/LongReads

#Assume running from directory nomenclature:  /some/path/<sample>/Xnn_analysis, can derive sample as follows:
sample=$(basename $(dirname $PWD))

if [ -d ../A01_mergeBams ]; then
  input_bam=$(ls ../A01_mergeBams/*_mrg.bam)
  sprefix=${sample/_mrg}
else
  input_bam=$(ls ../A0?_minimap2/*.st.bam)
  sprefix=$sample
fi

echo "Input: $input_bam Sample: $sample, prefix: $sprefix"

#De-duplicate UMIs
python ${CODE_DIR}/collapse_barcodes.py -b $input_bam
sed -n 1p ${sample}.bc_umi.summary_output.tsv >${sample}.bc_umi.summary_output.st.tsv
sed '1d' ${sample}.bc_umi.summary_output.tsv | sort -k1 >>${sample}.bc_umi.summary_output.st.tsv 

#Create readname list of reads to be kept, and filter bam file accordingly
sed '1d' ${sample}.bc_umi.collapsed.txt | awk '{print $1"#"$2}' >bc_read_list.txt

java -jar ${PICARD_DIR}/picard.jar FilterSamReads \
        I=$input_bam \
        O=${sample}.umi_dd.bam \
        READ_LIST_FILE=bc_read_list.txt \
        FILTER=includeReadList \
        CREATE_INDEX=true

#Picard creates <sample>.bai as index for <sample>.bam, but some tools such as vartrix don't recognize that 
#  format => modify to <sample>.bam.bai
mv ${sample}.umi_dd.bai ${sample}.umi_dd.bam.bai

#Need de-dup'd fastq for flair collapse step (can use bed for flair correct step)
###NOTE: seqtk sometimes (server specific???) runs for many hours with no progress. Issue TBD 
#if [ -d ../A01_mergeBams ]; then
#  $SEQTK subseq ${LEUKEMIA_LR}/${sprefix}_1/A01_flexiplex/${sprefix}_1.bc.fastq.gz bc_read_list.txt | gzip - >  ${sample}.umi_dd.fastq.gz
#  $SEQTK subseq ${LEUKEMIA_LR}/${sprefix}_2/A01_flexiplex/${sprefix}_2.bc.fastq.gz bc_read_list.txt | gzip - >> ${sample}.umi_dd.fastq.gz
#else
#  $SEQTK subseq ../A01_flexiplex/${sample}.bc.fastq.gz bc_read_list.txt | gzip - > ${sample}.umi_dd.fastq.gz
#fi

#Just use samtools to convert de-dup'd bam to fastq instead
#Split into two statements since sometimes combined statement resulted in incomplete file
#$SAMTOOLS fastq ${sample}.umi_dd.bam | gzip - > ${sample}.umi_dd.fastq.gz
$SAMTOOLS fastq --threads 16 ${sample}.umi_dd.bam >${sample}.umi_dd.fastq
gzip -f ${sample}.umi_dd.fastq
