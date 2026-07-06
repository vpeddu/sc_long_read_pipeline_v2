process runSTARIndex {
    tag "D01_star_index"

    publishDir path: { "${params.output_dir}/D01_shortread/star_index" }, mode: 'symlink'

    input:
    path genome
    path gtf

    output:
    path "star_index"

    script:
    """
    mkdir star_index
    /opt/miniconda3/envs/sqanti3/bin/STAR --runMode genomeGenerate \
        --genomeDir star_index \
        --genomeFastaFiles ${genome} \
        --sjdbGTFfile ${gtf} \
        --sjdbOverhang ${params.star_sjdb_overhang} \
        --runThreadN ${task.cpus}
    """
}

process runSTARAlign {
    tag "D01_star_${sample}"

    publishDir path: { "${sample}/${params.output_dir}/D01_shortread/star" }, mode: 'symlink'

    input:
    tuple val(sample), path(r1), path(r2)
    path star_index

    output:
    tuple val(sample), path("${sample}.SJ.out.tab")

    script:
    """
    # --outSAMtype None: only junctions are needed here, skip alignment output.
    # --twopassMode Basic: standard best-practice for junction-discovery
    # sensitivity (re-maps using junctions discovered in a first pass).
    /opt/miniconda3/envs/sqanti3/bin/STAR --runMode alignReads \
        --genomeDir ${star_index} \
        --readFilesIn ${r1} ${r2} \
        --readFilesCommand zcat \
        --runThreadN ${task.cpus} \
        --twopassMode Basic \
        --outSAMtype None \
        --outSJtype Standard \
        --outFileNamePrefix ${sample}.
    """
}

process runShortReadQuant {
    tag "D01_quant_${sample}"

    publishDir path: { "${sample}/${params.output_dir}/D01_shortread/quant" }, mode: 'symlink'

    input:
    tuple val(sample), path(txmod_fasta), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_kallisto")

    script:
    """
    /opt/miniconda3/envs/sqanti3/bin/kallisto index -i ${sample}.kidx ${txmod_fasta}
    /opt/miniconda3/envs/sqanti3/bin/kallisto quant -i ${sample}.kidx -o ${sample}_kallisto ${r1} ${r2}
    """
}
