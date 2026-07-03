process runFlair {
    //beforeScript 'export PATH=/opt/miniconda3/envs/flair/bin/:$PATH'

    tag "A04_flair"

    publishDir path: { "${sample}/${params.output_dir}/A04_flair" }, mode: 'symlink'

    input:
    tuple val(sample), path(bam), path(bai), path(fastq)
    path genome
    path gtf

    output:
    tuple val(sample),
        path("${sample}.flair.collapse.isoforms.bed"),
        path("${sample}.flair.collapse.isoforms.fa"),
        path("${sample}.flair.collapse.isoforms.gtf"),
        path("${sample}.flair.collapse.isoform.read.map.txt")

    script:
    """
    export PATH=/opt/miniconda3/envs/flair/bin/:\$PATH

    # flair transcriptome replaces correct+collapse: runs directly off the
    # sorted/indexed genome bam, and internally parallelizes isoform calling
    # by chromosome/region (see --parallelmode) instead of the single-threaded
    # collapse step that was previously the bottleneck.
    /opt/miniconda3/envs/flair/bin/flair transcriptome \
        --genomealignedbam ${bam} \
        --genome ${genome} \
        --gtf ${gtf} \
        --threads ${task.cpus} \
        --check_splice \
        --output ${sample}.flair.collapse
    """
}

process runTxRename {
    tag "A04_txRename"

    publishDir "${params.output_dir}/A04_txRename", mode: 'symlink'

    input:
    tuple val(sample),
    path(flair_collapse_isoforms_bed),
    path(flair_collapsed_isoforms_fa),
    path(flair_collapse_gtf),
    path(flair_collapse_map)

    output:
    tuple val(sample), path("${sample}.flair.collapse.isoforms.txmod.gtf"),
    path("${sample}.flair.collapse.isoform.read.map.txt"),
    path("isoform_cells.csv"),
    path("transcript_xref.tsv")

    script:
    """
    sample="${sample}"
    sample_bc="\${sample#*_}"
    sample_bc=\${sample_bc/_WES_mrg/WES}

    declare -A disease_code=( [Leukemia]="LEU" [CRC]="CRC" [GIM]="GIM" [Gastric]="GST" [Organoids]="GOO" )
    disease="GIM"
    tprefix="\${disease_code[\$disease]}\${sample_bc}"

    /opt/miniconda3/envs/long_reads/bin/python3 ${projectDir}/bin/rename_gtf_txid.py \
        --gtf ${flair_collapse_gtf} \
        --tx_prefix \$tprefix

    /opt/miniconda3/envs/long_reads/bin/python3 ${projectDir}/bin/rdmap2counts.py \
        --read_map ${flair_collapse_map} \
        --sample ${sample} \
        --matrix \
        --read_map ${flair_collapse_map} 

    /usr/bin/pigz -p ${task.cpus} ${sample}.*.tsv ${sample}.matrix.mtx
    """
}
