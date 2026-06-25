
process runLongshot {
    tag "B01_longshot"
    cpus 16
    memory '64 GB'

    input:
    tuple val(sample), path(bam), path(bai), path(fastq)
    path vep_data
    path genome_fasta
    path genome_fai

    
    output:
    path "*"
        publishDir "${sample}/${params.output_dir}/B01_longshot", mode: 'symlink'


    script:
    """
    /opt/miniconda3/envs/long_reads/bin/longshot -F \
        --bam ${bam} \
        --ref ${genome_fasta} \
        --out ${sample}.vcf \
         --strand_bias_pvalue_cutoff 0.0001 \
         --density_params 10:100:50 \
         --no_haps \
         --variant_debug_dir snv_debug

    grep '^#' ${sample}.vcf > ${sample}.pass.vcf
    awk '\$7 == "PASS"' ${sample}.vcf >> ${sample}.pass.vcf

    input_vcf=\$(ls ${sample}.pass.vcf)
    output_maf=\${input_vcf%vcf*}maf

    source /opt/miniconda3/bin/activate vep

#TODO cache version is hardcoded
    /opt/miniconda3/envs/vep/bin/vcf2maf.pl \
       --input-vcf \$input_vcf \
       --output-maf \$output_maf \
       --vep-overwrite \
       --tumor-id SAMPLE \
       --ref-fasta ${genome_fasta} \
       --ncbi-build GRCh38 \
       --vep-path /opt/miniconda3/envs/vep/bin/ \
       --vep-data ${vep_data} \
       --cache-version 104 
       
    sed '1d' \$output_maf | cut -d\$'\t' -f 9 | sort | uniq -c > maf_pass.cts.txt
    rm ${sample}.pass*.vcf

    idt='not_idt'
    if [[ "${sample}" == *PPC* ]]; then idt='PPC'; fi
    echo "Filter placeholder for idt=\$idt" > filt.out
    """
}
