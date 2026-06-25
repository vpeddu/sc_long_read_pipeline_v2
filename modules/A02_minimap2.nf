process runMinimap2 {
    tag "A02_minimap2"

    cpus 20
    memory '64 GB'

    input:
    tuple val(sample), file(flex_fastq)
    file junc_bed
    file ref_genome
    file bam_tag_bc_umi_py

    publishDir "${sample}/${params.output_dir}/A02_minimap2", mode: 'symlink'

    output:
    tuple val(sample), path("${sample}.bc_tag.st.bam"), path("${sample}.bc_tag.st.bam.bai")

    script:
    """

  #TODO: will be faster if we pre-build the minimap2 index

    echo "Running minimap2 for ${sample} with ${task.cpus} threads"

    /opt/miniconda3/envs/long_reads/bin/minimap2 -ax splice -k14 -t ${task.cpus > 4 ? task.cpus - 4 : 1} --secondary=no --junc-bed ${junc_bed} \
        ${ref_genome} ${flex_fastq} | \
        /usr/bin/samtools view -@ 4 -Sbh -F 2048 - > ${sample}.primary.bam

    /opt/miniconda3/envs/long_reads/bin/python3 ${bam_tag_bc_umi_py} ${sample}.primary.bam

    /usr/bin/samtools sort -@ ${task.cpus} ${sample}.primary.bc_tag.bam > ${sample}.bc_tag.st.bam
    /usr/bin/samtools index ${sample}.bc_tag.st.bam
    /usr/bin/samtools flagstat ${sample}.bc_tag.st.bam > ${sample}.flagstat.txt

  
    """
}
