process runSampleMetrics {
    tag "A07_metrics"

    publishDir path: { "${sample}/${params.output_dir}/A07_metrics" }, mode: 'symlink'

    input:
    tuple val(sample), path(raw_fastq), path(predup_bam), path(postdup_bam)

    output:
    tuple val(sample), path("${sample}.sample_metrics.tsv")

    script:
    """
    raw_reads=\$(( \$(/usr/bin/pigz -dc ${raw_fastq} | wc -l) / 4 ))
    predup_reads=\$(/usr/bin/samtools view -c -@ ${task.cpus} ${predup_bam})
    postdup_reads=\$(/usr/bin/samtools view -c -@ ${task.cpus} ${postdup_bam})

    {
        echo -e "metric\tvalue"
        echo -e "raw_reads\t\$raw_reads"
        echo -e "predup_reads\t\$predup_reads"
        echo -e "postdup_reads\t\$postdup_reads"
    } > ${sample}.sample_metrics.tsv
    """
}
