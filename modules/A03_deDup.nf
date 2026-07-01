process runDeDup {
    tag "A03_deDup"

    publishDir path: { "${sample}/${params.output_dir}/A03_deDup" }, mode: 'symlink'

    input:
    tuple val(sample), file(minimap_bam), file(minimap_bai)
    path collapse_barcodes_py
    path ref_genes

    output:
    tuple val(sample), path("${sample}.umi_dd.bam"), path("${sample}.umi_dd.bam.bai"), path("${sample}.umi_dd.fastq.gz")

    script:
    """
 
    # UMI deduplication
    /opt/miniconda3/envs/long_reads/bin/python3 ${collapse_barcodes_py} -b ${minimap_bam}
    sed -n 1p ${sample}.bc_umi.summary_output.tsv > ${sample}.bc_umi.summary_output.st.tsv
    sed '1d' ${sample}.bc_umi.summary_output.tsv | sort -k1 >> ${sample}.bc_umi.summary_output.st.tsv

    sed '1d' ${sample}.bc_umi.collapsed.txt | awk '{print \$1"#"\$2}' > bc_read_list.txt

    /opt/miniconda3/envs/picardtools/bin/picard FilterSamReads \
        I=${minimap_bam} \
        O=${sample}.umi_dd.bam \
        READ_LIST_FILE=bc_read_list.txt \
        FILTER=includeReadList \
        CREATE_INDEX=true

    mv ${sample}.umi_dd.bai ${sample}.umi_dd.bam.bai

    /usr/bin/samtools fastq --threads ${task.cpus} ${sample}.umi_dd.bam > ${sample}.umi_dd.fastq
    /usr/bin/pigz -p ${task.cpus} -f ${sample}.umi_dd.fastq

    # Exon coverage from deduplicated UMI BAM
    /opt/miniconda3/envs/picardtools/bin/featureCounts -a ${ref_genes} \
        -T ${task.cpus} -L -t exon -g gene_name -o exon_counts.txt ${sample}.umi_dd.bam
    """
}
