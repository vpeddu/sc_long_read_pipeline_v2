process runBuildSeurat {
    tag "A06_seurat"

    publishDir path: { "${sample}/${params.output_dir}/A06_seurat" }, mode: 'symlink'

    input:
    tuple val(sample),
        path(gene_matrix, stageAs: 'gene.matrix.mtx.gz'),
        path(gene_features, stageAs: 'gene.features.tsv.gz'),
        path(gene_barcodes, stageAs: 'gene.barcodes.tsv.gz'),
        path(iso_matrix, stageAs: 'iso.matrix.mtx.gz'),
        path(iso_features, stageAs: 'iso.features.tsv.gz'),
        path(iso_barcodes, stageAs: 'iso.barcodes.tsv.gz'),
        path(sqanti_gtf),
        path(sqanti_classification),
        path(transcript_xref),
        path(sample_metrics)
    path ref_genes_gtf

    output:
    tuple val(sample), path("${sample}.seurat.rds")

    script:
    """
    # Genomic ranges for isoforms/genes aren't in the classification/matrix files
    # themselves -- pull them straight out of the transcript/gene rows of the
    # respective GTFs (POSIX awk, no gawk-only match() array extension needed).
    awk -F'\t' '\$3=="transcript" {
        line = \$9
        if (match(line, /transcript_id "[^"]+"/)) {
            val = substr(line, RSTART+15, RLENGTH-16)
        } else { val = "NA" }
        print \$1"\t"\$4"\t"\$5"\t"\$7"\t"val
    }' ${sqanti_gtf} > iso_ranges.tsv

    awk -F'\t' '\$3=="gene" {
        line = \$9
        if (match(line, /gene_name "[^"]+"/)) {
            val = substr(line, RSTART+11, RLENGTH-12)
        } else { val = "NA" }
        print \$1"\t"\$4"\t"\$5"\t"\$7"\t"val
    }' ${ref_genes_gtf} > gene_ranges.tsv

    /opt/miniconda3/envs/seurat/bin/Rscript ${projectDir}/bin/build_seurat_object.R \
        ${sample} \
        gene.matrix.mtx.gz gene.features.tsv.gz gene.barcodes.tsv.gz \
        iso.matrix.mtx.gz iso.features.tsv.gz iso.barcodes.tsv.gz \
        iso_ranges.tsv gene_ranges.tsv ${sqanti_classification} ${sample_metrics} \
        ${transcript_xref}
    """
}
