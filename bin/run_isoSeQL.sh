#!/bin/bash
# This script runs within a Singularity container where all tools are available
# Note: isoSeQL installation may need to be added to the container or referenced from host
source /opt/miniconda3/bin/activate sqanti3

#Set flair step number (will be A03 for _mrg, otherwise A04)
FSTEP="${1:-A03}"

PROJ_DIR=/mnt/ix1/Projects/M102_241107_Multiome/GIM
CODE_DIR="${PROJ_DIR}/../00_code/LongReads"
ISOSEQL=~/shared/installs/isoSeQL

EXP_TEMPLATE=${PROJ_DIR}/Integration/novelIsoform/A01_isoSeQL/template.exp_config.txt
SMP_TEMPLATE=${PROJ_DIR}/Integration/novelIsoform/A01_isoSeQL/template.sample_config.txt
SAMPLE_META=${PROJ_DIR}/Integration/00_samples/sample_metadata.tsv

#Assuming running from /path/to/sample_dir/analysis_dir
sample_dir=$(basename $(dirname $PWD))

#Create experiment config files for this sample
sed "s/<experiment>/$sample_dir/" $EXP_TEMPLATE > exp_config.txt

#Create sample config file for this sample
  pt_sample=$(cut -d'_' -f1-2 <<< "$sample_dir")
  head -n1 $SMP_TEMPLATE > sample_config.txt
  
  awk -F'\t' -v id="$pt_sample" '
    NR>1 && $1==id {
      printf "%s,%s_%s,%s,N/A,%s\n", id, $5, $6, $2, $8
      exit
    }
  ' "$SAMPLE_META" >> sample_config.txt

#Link flair and SQANTI3 files if necessary
if [[ "$sample_dir" == *"_mrg" ]]; then
  fstep="A03"
else
  fstep="A04"
fi

if [[ ! -e transcript_xref.tsv ]]; then
   ln -s ../${fstep}_txRename/transcript_xref.tsv .
fi

if [[ ! -e sqanti_RulesFilter_result_classification.txt ]]; then
   ln -s ../C01_SQANTI3/sqanti_RulesFilter_result_classification.txt .
fi

if [[ ! -e tx_dups_xref.tsv ]]; then
   ln -s ../C01_SQANTI3/tx_dups_xref.tsv .
fi

python ${CODE_DIR}/summarize_classif.py

python ${ISOSEQL}/isoSeQL_run.py \
    -e exp_config.txt \
    -s sample_config.txt

python ${CODE_DIR}/summarize_classif.py

#Load experiment and sample, with associated isoform details to database
python ${ISOSEQL}/isoSeQL_run.py \
       --classif sqanti_classif_summ.txt \
       --genePred ../C01_SQANTI3/sqanti_corrected.genePred \
       --sampleConfig sample_config.txt \
       --expConfig exp_config.txt \
       --db ${PROJ_DIR}/Integration/novelIsoform/A01_isoSeQL/isoforms.db
