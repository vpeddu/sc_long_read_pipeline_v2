process runFlexiplex {
    tag "A01_flexiplex"
    
    cpus 20
    memory '32 GB'

    input:
    tuple val(sample), path(fastq), path(bc_tsv)
    val seq_type

    output:
    tuple val(sample), path("${sample}.bc.fastq.gz") 

    publishDir "${sample}/${params.output_dir}/A01_flexiplex", mode: 'symlink'

    script:
    """

    flank_len=${params.multiseq ? 17 : 24}
    if [[ "${seq_type}" == '5_prime' ]]; then
        seq_type='5prm'
    else
        seq_type='3prm'
    fi

    adapter="CTACACGACGCTCTTCCGATCT"
    if [[ \$flank_len -lt 24 ]]; then
        adapter=\${adapter: -\$flank_len}
    fi
/opt/miniconda3/envs/long_reads/bin/flexiplex -p ${task.cpus > 4 ? task.cpus - 8 : 1} -x "\$adapter" -b ???????????????? -u ????????????????
    if [[ "\$seq_type" == '3prm' ]]; then
        /usr/bin/pigz -dc -p 2 ${fastq} | /opt/miniconda3/envs/long_reads/bin/flexiplex -p ${task.cpus > 4 ? task.cpus - 8 : 1} -x "\$adapter" -b ???????????????? -u ???????????????? \
                           -x TTTTTTTTT -f 8 -e 2 -k ${bc_tsv} -n ${sample}.bc | \
                        /usr/bin/pigz -p 6 > ${sample}.bc.fastq.gz
    else
        /usr/bin/pigz -dc -p 2 ${fastq} |/opt/miniconda3/envs/long_reads/bin/flexiplex -p ${task.cpus > 4 ? task.cpus - 8 : 1} -d 10x5v2 -k ${bc_tsv} -n ${sample}.bc | \
                        /usr/bin/pigz -p 6 > ${sample}.bc.fastq.gz
    fi


    """
}
