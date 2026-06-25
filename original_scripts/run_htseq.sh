#!/bin/bash
#Run from long_reads env (need env with relatively new version of htseq to get htseq-count-barcodes)
#If writing out a mtx format count matrix, need scipy also
#Ultimately make multiome venv the environment to run this in - htseq not currently installed there
source activate long_reads

REF_DIR='/mnt/ix1/Projects/M102_241107_Multiome/00_resources'
GENOME_REF="${REF_DIR}/genome.fa"
GENES_GTF="${REF_DIR}/genes.gtf"

RESOURCE_DIR='/mnt/ix1/Resources'
#HTSEQ_DIR="${RESOURCE_DIR}/tools/HTSeq/v0.5.4p5/scripts"

#Assume running from directory nomenclature:  /some/path/<sample>/Xnn_analysis, can derive sample as follows:
sample=$(basename $(dirname $PWD))

#Maybe decrease minimum alignment quality value currently at default: 10
#Decrease to 1 to effectively eliminate quality check
htseq-count-barcodes --format=bam --order=pos --stranded=no \
	--type=gene --idattr=gene_name --mode=intersection-nonempty --nonunique=all \
        --cell-barcode=CB --UMI=UB --minaqual=10 \
        -c ${sample}.bc_counts.mtx --counts_output_sparse \
        -o ${sample}.annot.bam -p BAM \
	${sample}.umi_dd.bam $GENES_GTF

#Remove id column heading in features file
sed -i '1d' ${sample}.bc_counts_features.tsv

#Add '-1' suffix to barcodes
sed -i 's/$/-1/' ${sample}.bc_counts_samples.tsv

#Gzip files and modify file names to match with cellranger type output
gzip ${sample}.bc_counts_*.tsv
gzip ${sample}.bc_counts.mtx
mv ${sample}.bc_counts_samples.tsv.gz ${sample}.barcodes.tsv.gz
mv ${sample}.bc_counts_features.tsv.gz ${sample}.features.tsv.gz
mv ${sample}.bc_counts.mtx.gz ${sample}.matrix.mtx.gz
