process runMinimap2 {
    tag "A02_minimap2"

    cpus 20
    memory '64 GB'

    publishDir path: { "${sample}/${params.output_dir}/A02_minimap2/chunks" }, mode: 'symlink'

    input:
    tuple val(sample), file(flex_fastq)
    file junc_bed
    file ref_genome

    output:
    tuple val(sample), path("*.bc_tag.st.bam"), path("*.bc_tag.st.bam.bai")

    script:
    """

  #TODO: will be faster if we pre-build the minimap2 index

    prefix=\$(basename ${flex_fastq} .bc.fastq.gz)

    echo "Running minimap2 for \${prefix} with ${task.cpus} threads"

    # -y carries flexiplex's CB:Z/UB:Z tags from the fastq header into the output BAM
    /opt/miniconda3/envs/long_reads/bin/minimap2 -ax splice -y -k14 -t ${task.cpus > 4 ? task.cpus - 4 : 1} --secondary=no --junc-bed ${junc_bed} \
        ${ref_genome} ${flex_fastq} | \
        /usr/bin/samtools view -@ 4 -Sbh -F 2048 - > \${prefix}.primary.bam

    /usr/bin/samtools sort -@ ${task.cpus} \${prefix}.primary.bam > \${prefix}.bc_tag.st.bam
    /usr/bin/samtools index \${prefix}.bc_tag.st.bam
    """
}

process mergeBam {
    tag "A02_mergeBam"

    cpus 8
    memory '32 GB'

    publishDir path: { "${sample}/${params.output_dir}/A02_minimap2" }, mode: 'symlink'

    input:
    tuple val(sample), path(chunk_bams), path(chunk_bais)

    output:
    tuple val(sample), path("${sample}.bc_tag.st.bam"), path("${sample}.bc_tag.st.bam.bai")

    script:
    """
    echo "Merging minimap2 chunk BAMs for ${sample}"

    /usr/bin/samtools merge -@ ${task.cpus} -f ${sample}.bc_tag.st.bam ${chunk_bams}
    /usr/bin/samtools index ${sample}.bc_tag.st.bam
    /usr/bin/samtools flagstat ${sample}.bc_tag.st.bam > ${sample}.flagstat.txt
    """
}
