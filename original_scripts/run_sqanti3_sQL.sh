#!/bin/bash
source activate /venvs2/anaconda3/envs/lrseq

PROJ_DIR=/mnt/ix1/Projects/M102_241107_Multiome
CODE_DIR=${PROJ_DIR}/00_code/LongReads

SQANTI3=/home/sgrimes/shared/installs/sqanti3_iSQL #Modified version of SQANTI3 from isoSeQL author

#usage: sqanti3_qc.py [-h] --isoforms ISOFORMS --refGTF REFGTF --refFasta
#                     REFFASTA [--min_ref_len MIN_REF_LEN] [--force_id_ignore]
#                     [--fasta] [--genename]
#                     [--novel_gene_prefix NOVEL_GENE_PREFIX] [-s SITES]
#                     [-w WINDOW]
#                     [--aligner_choice {minimap2,deSALT,gmap,uLTRA}]
#                     [-x GMAP_INDEX] [--skipORF] [--orf_input ORF_INPUT]
#                     [--short_reads SHORT_READS] [--SR_bam SR_BAM]
#                     [--CAGE_peak CAGE_PEAK]
#                     [--polyA_motif_list POLYA_MOTIF_LIST]
#                     [--polyA_peak POLYA_PEAK] [--phyloP_bed PHYLOP_BED]
#                     [-e EXPRESSION] [-c COVERAGE] [-fl FL_COUNT]
#                     [--isoAnnotLite] [--gff3 GFF3] [-o OUTPUT] [-d DIR]
#                     [--saturation] [--report {html,pdf,both,skip}]
#                     [--isoform_hits]
#                     [--ratio_TSS_metric {max,mean,median,3quartile}]
#                     [-t CPUS] [-n CHUNKS] [-l {ERROR,WARNING,INFO,DEBUG}]
#                     [--is_fusion] [-v] [--tusco {human,mouse}]

#Create isoform_counts.tsv for input to SQANTI3
python ${CODE_DIR}/isoform_summ_cts.py

python ${SQANTI3}/sqanti3_qc.py -t 20 \
         --isoforms flair.collapse.isoforms.txmod.gtf \
         --refGTF ${PROJ_DIR}/00_resources/genes.gtf \
         --refFasta ${PROJ_DIR}/00_resources/genome.fa \
         --CAGE_peak ${SQANTI3}/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed \
         --polyA_motif_list ${SQANTI3}/data/polyA_motifs/mouse_and_human.polyA_motif.txt \
         --force_id_ignore \
         --fl_count isoform_counts.tsv \
         --output sqanti --dir $PWD --report skip

python ${SQANTI3}/sqanti3_filter.py rules --cpus 20 \
        --sqanti_class sqanti_classification.txt \
        --filter_gtf sqanti_corrected.gtf \
        --filter_faa sqanti_corrected.fasta \
        --dir $PWD

#Find duplicate tx structures in sqanti_corrected.genePred and write 
#  duplicates to tx_dups_xref.tsv
python ${CODE_DIR}/find_dup_isoforms.py 
 
#Read sqanti_classification.txt and output: full_length_FSM_ISM.csv
python ${CODE_DIR}/sqanti_full_length.py 

: <<'COMMENT'
COMMENT
