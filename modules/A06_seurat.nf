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
        path(iso_barcodes, stageAs: 'iso.barcodes.tsv.gz')

    output:
    tuple val(sample), path("${sample}.seurat.rds")

    script:
    """
    /opt/miniconda3/envs/seurat/bin/Rscript ${projectDir}/bin/build_seurat_object.R \
        ${sample} \
        gene.matrix.mtx.gz gene.features.tsv.gz gene.barcodes.tsv.gz \
        iso.matrix.mtx.gz iso.features.tsv.gz iso.barcodes.tsv.gz
    """
}
