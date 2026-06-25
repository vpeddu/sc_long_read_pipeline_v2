#!/bin/bash

CODE_DIR=/mnt/ix1/Projects/M102_241107_Multiome/00_code/LongReads
MCODE_DIR=/mnt/ix1/Projects/M102_241107_Multiome/02_metrics/LongReads

#txRename step requires avx2 support (fugu, iwashi, suzuki), if no avx2 skip that step
shopt -s extglob
if [[ "$HOSTNAME" == @(fugu|iwashi|suzuki) ]]; then
  avx2='Y'
else
  avx2='N'
fi

# Function to print help
print_usage()
{
	   echo "Usage: $(basename "$0") [-m] [step_numbers]"; 
	   echo "Where -m multi-seq processing, shorter adapter sequence";
           echo "Positional arguments are step numbers to run";
	   return
}
multiseq=0

# Parse command line options
OPTIND=1
while getopts "mh" OPT
do
  case "$OPT" in
    m) multiseq=1;;
    h) print_usage; exit 1;;
   \?) print_usage; exit 1;;
    :) echo "Option -$OPTARG requires an argument."; print_usage; exit 1;;
  esac
done

#Remove keyword argument from args, and keep any positional arguments (step #s)
if [ $multiseq -eq 1 ]; then
  param_m="-m"
  shift $((OPTIND - 1))
fi

#Assumes running from directory which is equivalent to sample name (eg SU968_2929A)
sample=$(basename $PWD)
echo "Running pipeline for sample: $sample"

#Go into flexiplex directory to ensure all cd commands can use relative link
mkdir -p A01_flexiplex
cd A01_flexiplex  

#Code to allow specifying which steps to run.  
array_contains () {
     local array="$1[@]"
     local seeking=$2
     local found=1
     for element in "${!array}"; do
         if [[ $element == "$seeking" ]]; then
             found=0
             break
         fi
     done
     return $found
 }

#Run all steps unless specific step numbers provided on command line 
#Eg: bash run_pipeline.sh 4 5 will just run flair and htSeq

if [ "$#" -eq 0 ]; then
  steps=( 1 2 3 4 5 6 7 8 )
else
  steps=( "$@" )
fi

echo "Running steps: ${steps[@]}"

if array_contains steps "1"; then
  cd ../A01_flexiplex
  echo "A01: Flexiplex barcode assignment started"
  echo "bash ${CODE_DIR}/run_flexiplex.sh $param_m >flex.out"
  bash ${CODE_DIR}/run_flexiplex.sh $param_m >flex.out
  echo "A01: Flexiplex barcode assignment completed"
fi

if array_contains steps "2"; then
  mkdir -p ../A02_minimap2
  cd ../A02_minimap2
  echo "A02: Minimap2 alignment started"
  bash ${CODE_DIR}/run_minimap.sh >mm.out
  echo "A02: Minimap2 alignment completed"
fi 

if array_contains steps "3"; then
  mkdir -p ../A03_deDup
  cd ../A03_deDup
  echo "A03: BC/UMI de-duplication started"
  bash ${CODE_DIR}/umi_dedup.sh >dd.out
  bash ${CODE_DIR}/exon_coverage.sh >cov.out
  echo "A03: BC/UMI de-duplication completed"
fi

if array_contains steps "4"; then
  mkdir -p ../A04_flair
  cd ../A04_flair
  echo "A04: flair isoform assignment started"
  bash ${CODE_DIR}/run_flair.sh >flr.out
  echo "A04: flair isoform assignment completed"

  mkdir -p ../A04_txRename
  cd ../A04_txRename
  echo "A04: flair transcript renaming started"
  ln -sf ../A04_flair/flair.collapse.isoforms.gtf .
  ln -sf ../A04_flair/flair.collapse.combined.isoform.read.map.txt .
  if [[ "$avx2" == "Y" ]]; then
    bash ${CODE_DIR}/modify_novelTx.sh >ct.out
    echo "A04: flair transcript renaming completed"
  else
    echo "A04: flair transcript renaming skipped - no avx2 support"
  fi
fi

if array_contains steps "5"; then
  echo "Array contains step5"
  mkdir -p ../A05_htSeq
  cd ../A05_htSeq
  echo "A05: Gene expression matrix creation started"
  ln -sf ../A03_deDup/${sample}.umi_dd.bam .
  ln -sf ../A03_deDup/${sample}.umi_dd.*bai .
  bash ${CODE_DIR}/run_htseq.sh >hts.out
  Rscript ${MCODE_DIR}/metrics_lr_gex.R $sample
  echo "A05: Gene expression matrix creation completed"
fi

if array_contains steps "6"; then
  mkdir -p ../B01_longshot
  cd ../B01_longshot
  echo "B01: Longshot SNV calling started"
  ln -sf ../A03_deDup/${sample}.umi_dd.bam .
  ln -sf ../A03_deDup/${sample}.umi_dd.*bai .
  bash ${CODE_DIR}/longshot_snv.sh >snv.out
  bash ${CODE_DIR}/convert2maf.sh >maf.out
  
  idt='not_idt'
  if [[ "$sample" =~ "PPC" ]]; then idt='PPC'; fi
  bash ${CODE_DIR}/filtermaf.sh $idt >filt.out
  echo "B01: Longshot SNV calling completed"
fi

if array_contains steps "7"; then
  mkdir -p ../C01_SQANTI3
  cd ../C01_SQANTI3
  echo "C01: SQANTI3 isoform classification started"
  ln -sf ../A03_txRename/flair.collapse.isoforms.txmod.gtf .
  ln -sf ../A03_txRename/transcript_xref.tsv .
  ln -sf ../A03_txRename/isoform_cells.csv

  python ${CODE_DIR}/isoform_summ_cts.py
  bash ${CODE_DIR}/run_sqanti3_sQL.sh >qc.out

  echo "C01: SQANTI3 isoform classification completed"
fi

if array_contains steps "8"; then
  mkdir -p ../C02_isoSeQL
  cd ../C02_isoSeQL
  echo "C02: isoSeQL sample isoforms load to database started"
    
  bash ${CODE_DIR}/run_isoSeQL.sh >qc.out

  echo "C02: isoSeQL sample isoforms load to database completed"
fi

echo "***All done***"
