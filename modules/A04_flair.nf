process runFlair {
    //beforeScript 'export PATH=/opt/miniconda3/envs/flair/bin/:$PATH'

    tag "A04_flair"

    cpus 20
    memory '128 GB'

    input:
    tuple val(sample), path(bam), path(bai), path(fastq)
    path genome
    path gtf

    output:
    tuple val(sample), 
        path ("${sample}.flair.filtered_all_corrected.bed"), 
        path("${sample}.flair.collapse.isoforms.bed"), 
        path("${sample}.flair.collapse.isoforms.fa"), 
        path("${sample}.flair.collapse.isoforms.gtf"),
        path("${sample}.flair.collapse.isoform.read.map.txt")

    publishDir "${sample}/${params.output_dir}/A04_flair", mode: 'symlink'

    script:
    """
    export PATH=/opt/miniconda3/envs/flair/bin/:\$PATH

    /opt/miniconda3/envs/flair/bin/bam2Bed12 \
        -i ${bam} > ${sample}.umi_dd.bed

    /opt/miniconda3/envs/flair/bin/flair correct --threads ${task.cpus} \
        -q ${sample}.umi_dd.bed \
        --gtf ${gtf} \
        --genome ${genome} \
        --nvrna \
        --output ${sample}.flair.filtered

    /opt/miniconda3/envs/flair/bin/flair collapse --threads ${task.cpus} \
        -q ${sample}.flair.filtered_all_corrected.bed \
        --reads ${fastq} \
        --genome ${genome} \
        --gtf ${gtf} \
        --annotation_reliant generate \
        --generate_map \
        --trust_ends \
        --no_gtf_end_adjustment \
        --check_splice \
        --quality 10 \
        --output ${sample}.flair.collapse
    """
}

process runTxRename {
    tag "A04_txRename"

    cpus 1
    memory '24 GB'

    input:
    tuple val(sample),
    path (flair_corrected_bed),
    path(flair_collapse_isoforms_bed),
    path(flair_collapsed_isoforms_fa), 
    path(flair_collapse_gtf),
    path(flair_collapse_map)
    
    output:
    tuple val(sample), path("${sample}.flair.collapse.isoforms.txmod.gtf"), 
    path("${sample}.flair.collapse.isoform.read.map.txt"), 
    path("isoform_cells.csv"),
    path("transcript_xref.tsv")
    
    publishDir "${params.output_dir}/A04_txRename", mode: 'symlink'


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
