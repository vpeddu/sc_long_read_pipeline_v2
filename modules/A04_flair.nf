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
    #
    # --threads controls BOTH the number of concurrent region workers AND the
    # -t passed to each worker's own internal minimap2 call (not divided
    # between them) -- passing task.cpus directly oversubscribes the SLURM
    # allocation by up to task.cpus^2 threads (confirmed via Sherlock's
    # resource-usage-mismatch warning: many simultaneous `minimap2 -t <cpus>`
    # processes). Bound it so worst-case (threads * threads) stays in budget.
    flair_threads=\$(awk -v c=${task.cpus} 'BEGIN { printf "%d", sqrt(c) }')
    /opt/miniconda3/envs/flair/bin/flair transcriptome \
        --genomealignedbam ${bam} \
        --genome ${genome} \
        --gtf ${gtf} \
        --threads \$flair_threads \
        --check_splice \
        --output ${sample}.flair.collapse
    """
}

process runFlairByChrom {
    tag "A04_flair_${chrom_label}"

    publishDir path: { "${sample}/${params.output_dir}/A04_flair/chunks" }, mode: 'symlink'

    input:
    tuple val(sample), path(bam), path(bai), path(fastq), val(chrom_label), val(chroms)
    path genome
    path gtf

    output:
    tuple val(sample),
        path("${sample}.${chrom_label}.flair.collapse.isoforms.bed"),
        path("${sample}.${chrom_label}.flair.collapse.isoforms.fa"),
        path("${sample}.${chrom_label}.flair.collapse.isoforms.gtf"),
        path("${sample}.${chrom_label}.flair.collapse.isoform.read.map.txt")

    script:
    """
    export PATH=/opt/miniconda3/envs/flair/bin/:\$PATH

    /usr/bin/samtools view -@ ${task.cpus} -b ${bam} ${chroms.join(' ')} > ${chrom_label}.bam
    /usr/bin/samtools index ${chrom_label}.bam

    flair_threads=\$(awk -v c=${task.cpus} 'BEGIN { printf "%d", sqrt(c) }')
    /opt/miniconda3/envs/flair/bin/flair transcriptome \
        --genomealignedbam ${chrom_label}.bam \
        --genome ${genome} \
        --gtf ${gtf} \
        --threads \$flair_threads \
        --check_splice \
        --output ${sample}.${chrom_label}.flair.collapse || true

    # Chromosomes/contigs with no aligned reads (common for e.g. chrY in a
    # female sample, or unplaced scaffolds) make flair exit without writing
    # output; touch empty placeholders so the merge step has something to cat.
    touch ${sample}.${chrom_label}.flair.collapse.isoforms.bed
    touch ${sample}.${chrom_label}.flair.collapse.isoforms.fa
    touch ${sample}.${chrom_label}.flair.collapse.isoforms.gtf
    touch ${sample}.${chrom_label}.flair.collapse.isoform.read.map.txt
    """
}

process mergeFlairChroms {
    tag "A04_flair_merge"

    publishDir path: { "${sample}/${params.output_dir}/A04_flair" }, mode: 'symlink'

    input:
    tuple val(sample), path(isoform_beds), path(isoform_fas), path(isoform_gtfs), path(read_maps)

    output:
    tuple val(sample),
        path("${sample}.flair.collapse.isoforms.bed"),
        path("${sample}.flair.collapse.isoforms.fa"),
        path("${sample}.flair.collapse.isoforms.gtf"),
        path("${sample}.flair.collapse.isoform.read.map.txt")

    script:
    """
    # Safe to concatenate directly: chromosomes/contigs are non-overlapping,
    # reference-matched isoform IDs are locus-specific (can't repeat across
    # chromosomes), and novel isoform IDs embed a globally-unique read name.
    cat ${isoform_beds} > ${sample}.flair.collapse.isoforms.bed
    cat ${isoform_fas} > ${sample}.flair.collapse.isoforms.fa
    cat ${isoform_gtfs} > ${sample}.flair.collapse.isoforms.gtf
    cat ${read_maps} > ${sample}.flair.collapse.isoform.read.map.txt
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
