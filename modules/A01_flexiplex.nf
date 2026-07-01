process splitFastq {
    tag "A01_splitFastq"

    publishDir path: { "${sample}/${params.output_dir}/A01_flexiplex/chunks" }, mode: 'symlink'

    input:
    tuple val(sample), path(fastq), path(bc_tsv)
    val n_splits

    output:
    tuple val(sample), path("${sample}.chunk*.fastq.gz"), path(bc_tsv)

    script:
    """
    echo "Splitting ${fastq} into ${n_splits} chunk(s) for ${sample}"

    # Distribute reads round-robin across N chunks (each read = 4 lines)
    /usr/bin/pigz -dc -p ${task.cpus} ${fastq} | \
        awk -v n=${n_splits} -v s="${sample}" '
            { chunk = int((NR - 1) / 4) % n
              fname = s ".chunk" chunk ".fastq"
              print > fname }'

    for f in ${sample}.chunk*.fastq; do
        /usr/bin/pigz -p ${task.cpus} -f "\$f"
    done
    """
}

process runFlexiplex {
    tag "A01_flexiplex"

    publishDir path: { "${sample}/${params.output_dir}/A01_flexiplex/chunks" }, mode: 'symlink'

    input:
    tuple val(sample), path(fastq), path(bc_tsv)
    val seq_type

    output:
    tuple val(sample), path("*.bc.fastq.gz")

    script:
    """
    prefix=\$(basename ${fastq} .fastq.gz)

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
/opt/miniconda3/envs/long_reads/bin/flexiplex -p ${task.cpus > 8 ? task.cpus - 8 : 1} -x "\$adapter" -b ???????????????? -u ????????????????
    if [[ "\$seq_type" == '3prm' ]]; then
        /usr/bin/pigz -dc -p 2 ${fastq} | /opt/miniconda3/envs/long_reads/bin/flexiplex -p ${task.cpus > 8 ? task.cpus - 8 : 1} -x "\$adapter" -b ???????????????? -u ???????????????? \
                           -x TTTTTTTTT -f 8 -e 2 -k ${bc_tsv} -n \${prefix}.bc | \
                        /usr/bin/pigz -p 6 > \${prefix}.bc.fastq.gz
    else
        /usr/bin/pigz -dc -p 2 ${fastq} |/opt/miniconda3/envs/long_reads/bin/flexiplex -p ${task.cpus > 8 ? task.cpus - 8 : 1} -d 10x5v2 -k ${bc_tsv} -n \${prefix}.bc | \
                        /usr/bin/pigz -p 6 > \${prefix}.bc.fastq.gz
    fi


    """
}
