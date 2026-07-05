process runHtseq {
    tag "A05_htSeq"

    publishDir path: { "${sample}/${params.output_dir}/A05_htSeq" }, mode: 'symlink'

    input:
    tuple val(sample), path(bam), path(bai), path(fastq), path(genes_gtf)
    val id_attr

    output:
    path "*"

    script:
    """

    /opt/miniconda3/envs/long_reads/bin/htseq-count-barcodes --format=bam \
        --order=pos \
        --stranded=no \
        --type=gene \
        --idattr=${id_attr} \
        --mode=intersection-nonempty \
        --nonunique=all \
        --cell-barcode=CB \
        --UMI=UB \
        --minaqual=10 \
        -c ${sample}.bc_counts.mtx \
        --counts_output_sparse \
        -o ${sample}.annot.bam -p BAM \
        ${bam} ${genes_gtf}

    sed -i '1d' ${sample}.bc_counts_features.tsv
    sed -i 's/\$/-1/' ${sample}.bc_counts_samples.tsv
    /usr/bin/pigz -p ${task.cpus} ${sample}.bc_counts_*.tsv
    /usr/bin/pigz -p ${task.cpus} ${sample}.bc_counts.mtx

    mv ${sample}.bc_counts_samples.tsv.gz ${sample}.barcodes.tsv.gz
    mv ${sample}.bc_counts_features.tsv.gz ${sample}.features.tsv.gz
    mv ${sample}.bc_counts.mtx.gz ${sample}.matrix.mtx.gz

    """
}
