process runHtseq {
    tag "A05_htSeq"

    publishDir path: { "${sample}/${params.output_dir}/A05_htSeq" }, mode: 'symlink'

    input:
    tuple val(sample), path(bam), path(bai), path(fastq), path(genes_gtf)
    val id_attr

    output:
    tuple val(sample), path("${sample}.matrix.mtx.gz"),
    path("${sample}.features.tsv.gz"),
    path("${sample}.barcodes.tsv.gz"),
    path("${sample}.annot.bam")

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

process runChrmGeneSubset {
    tag "A05_chrm_genes"

    input:
    path ref_genes_gtf

    output:
    path "chrm_genes.htseq_compat.gtf"

    script:
    // --exclude_chrm drops chrM from the flair-derived gtf entirely, so when
    // --keep_intergenic quantifies against that gtf, chrM genes would
    // otherwise get zero htseq counts. Backfill them straight from the
    // reference annotation instead. Only transcript/exon rows are kept
    // (matching the feature types flair's own gtf actually contains -- it
    // has no 'gene'-type rows), and gene_id is rewritten to hold the gene
    // SYMBOL, since that's what runTxRename repurposes gene_id to hold after
    // renaming (flair's gtf has no separate gene_name attribute) -- so a
    // single --idattr=gene_id htseq pass counts both sets of genes together.
    """
    awk -F'\t' '\$1=="chrM" && (\$3=="transcript" || \$3=="exon") {
        line = \$9
        if (match(line, /gene_name "[^"]+"/)) {
            gname = substr(line, RSTART+11, RLENGTH-12)
            sub(/gene_id "[^"]+"/, "gene_id \"" gname "\"", line)
        }
        print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\t"line
    }' ${ref_genes_gtf} > chrm_genes.htseq_compat.gtf
    """
}

process mergeChrmReferenceGenes {
    tag "A05_merge_chrm_${sample}"

    input:
    tuple val(sample), path(txmod_gtf)
    path chrm_genes_gtf

    output:
    tuple val(sample), path("${sample}.htseq_input.gtf")

    script:
    """
    cat ${txmod_gtf} ${chrm_genes_gtf} > ${sample}.htseq_input.gtf
    """
}
