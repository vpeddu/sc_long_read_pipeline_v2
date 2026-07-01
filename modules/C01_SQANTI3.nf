process runSQANTI3 {
    tag "C01_SQANTI3"

    publishDir path: { "${sample}/${params.output_dir}/C01_SQANTI3" }, mode: 'symlink'

    input:
    path reference_gtf
    path reference_genome
    tuple val(sample), path(txmod_gtf), path(flair_isoform_read_map), path(isoform_cells_csv), path(transcript_xref)

    output:
    tuple val(sample), path("sqanti_corrected.gtf"), path("sqanti_corrected.fasta"), path("sqanti_RulesFilter_classification.txt"), path("transcript_xref.tsv"), path("sqanti_corrected.genePred"), path("tx_dups_xref.tsv")

    script:
    """
    set +u # need this for conda activation for some reason
    
    export PYTHONPATH="/opt/miniconda3/envs/sqanti3/bin/:\$PYTHONPATH"

    source /opt/miniconda3/bin/activate sqanti3

    /opt/miniconda3/envs/sqanti3/bin/python3 ${projectDir}/bin/isoform_summ_cts.py

    /opt/miniconda3/envs/sqanti3/bin/python3 /sqanti3/sqanti3_qc.py -t ${task.cpus} \
             ${txmod_gtf} \
             ${reference_gtf} \
             ${reference_genome} \
             --cage_peak /sqanti3/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed \
             --polyA_motif_list /sqanti3/data/polyA_motifs/mouse_and_human.polyA_motif.txt \
             --force_id_ignore \
             --fl_count isoform_counts.tsv \
             --output sqanti \
             --dir \$PWD \
             --report skip \
             --skipORF

  wget --no-check-certificate https://raw.githubusercontent.com/christine-liu/SQANTI3/refs/heads/master/sqanti3_filter.py

git clone https://github.com/ConesaLab/SQANTI3
mv SQANTI3/src/ .

    /opt/miniconda3/envs/sqanti3/bin/python3 /SQANTI3/sqanti3_filter.py rules --cpus ${task.cpus} \
            --sqanti_class sqanti_classification.txt \
            --filter_gtf sqanti_corrected.gtf \
            --filter_faa sqanti_corrected.fasta \
            --dir \$PWD

    /opt/miniconda3/envs/sqanti3/bin/python3 ${projectDir}/bin/find_dup_isoforms.py
    /opt/miniconda3/envs/sqanti3/bin/python3 ${projectDir}/bin/sqanti_full_length.py
    """
}
