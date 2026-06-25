#!/bin/bash
#TODO: Modify this later to have specific panel as input, might not always be IDT_pancancer

idt=${1:-notidt}
PANEL_DIR=/mnt/ix1/Projects/M102_241107_Multiome/00_resources/hybpanel

#Assume running from directory nomenclature:  /some/path/<sample>/Xnn_analysis, can derive sample as follows:
sample=$(basename $(dirname $PWD))

if [ ${idt^^} == 'PPC' ]; then
  #Filter for genes of interest (those in IDT pancancer panel)
  IDT_GENES=${PANEL_DIR}/IDT_pancancer.txt
  #Remove maf version info in row 1, and retain header from row 2
  awk 'NR==FNR {genes[$1]; next} FNR==2; $1 in genes' $IDT_GENES ${sample}.pass.maf >${sample}.pass.goi.maf
fi

#Filter for any non-synonymous coding mutation (Missense_Mutation, Nonsense_Mutation, Nonstop_Mutation)
sed '1d' ${sample}.pass.maf | awk '(NR == 1 || $9 ~ /Mutation$/)' >${sample}.ns_coding.maf