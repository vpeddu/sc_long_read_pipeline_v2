process runPfamDownload {
    tag "A08_pfam_download"

    publishDir path: { "${params.output_dir}/A08_orf/pfam_db" }, mode: 'symlink'

    output:
    tuple path("Pfam-A.hmm"), path("Pfam-A.hmm.h3m"), path("Pfam-A.hmm.h3i"), path("Pfam-A.hmm.h3f"), path("Pfam-A.hmm.h3p")

    script:
    // A fixed, universal reference database (unlike eg the STAR genome index,
    // which depends on the user's own --ref_dir) -- fetched once here rather
    // than baked into the container, since it's ~400MB compressed/~3.5GB
    // unpacked and would otherwise bloat every pull of the image. Cached
    // across -resume'd runs the same way runSTARIndex already is.
    """
    wget --no-check-certificate https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip Pfam-A.hmm.gz
    /opt/miniconda3/envs/orfpred/bin/hmmpress Pfam-A.hmm
    """
}

process runORFPrediction {
    tag "A08_orf_${sample}"

    publishDir path: { "${sample}/${params.output_dir}/A08_orf" }, mode: 'symlink'

    input:
    tuple val(sample), path(txmod_fasta)
    tuple path(pfam_hmm), path(pfam_h3m), path(pfam_h3i), path(pfam_h3f), path(pfam_h3p)

    output:
    tuple val(sample), path("${sample}.TD2.pep"), path("${sample}.pfam.domtblout"),
        path("${sample}.TD2.cds"), path("${sample}.TD2.gff3"), path("${sample}.TD2.bed")

    script:
    // TD2 (a from-scratch TransDecoder successor, MIT-licensed) + PSAURON
    // (its internal coding-likelihood scorer) replace the isoSeQL-pinned
    // SQANTI3 fork's own ORF prediction (--skipORF upstream), which is stuck
    // on the license-gated GeneMarkS-T and can't be enabled here. This is a
    // fully independent step -- it never touches the SQANTI3/isoSeQL branch.
    // Standard TransDecoder-style workflow: LongOrfs extracts candidate ORFs,
    // hmmsearch scores them against Pfam (using --cut_ga's built-in Pfam
    // gathering thresholds instead of a fixed e-value; hmmsearch, not
    // hmmscan, since Pfam-as-query is far faster when scanning many small
    // proteins -- it's read once instead of once per query), and Predict
    // makes the final call, retaining any candidate with a Pfam hit even if
    // its own psauron score alone wouldn't have passed.
    // TD2.Predict shells out to the bare `psauron` command rather than
    // resolving it through the same interpreter, so it must be on PATH.
    """
    export PATH=/opt/miniconda3/envs/orfpred/bin:\$PATH

    TD2.LongOrfs -t ${txmod_fasta} -O td2_out

    hmmsearch --cpu ${task.cpus} --domtblout ${sample}.pfam.domtblout --cut_ga \
        ${pfam_hmm} td2_out/longest_orfs.pep

    # TD2.Predict crashes (pandas EmptyDataError) if --retain-hmmer_hits points
    # to a domtblout with zero actual hit rows -- common whenever a sample's
    # candidate ORFs simply have no Pfam match. Only pass the flag when at
    # least one non-comment (real hit) line exists.
    retain_flag=""
    if grep -qv '^#' ${sample}.pfam.domtblout; then
        retain_flag="--retain-hmmer_hits ${sample}.pfam.domtblout"
    fi

    TD2.Predict -t ${txmod_fasta} -O td2_out \$retain_flag

    mv ${txmod_fasta}.TD2.pep ${sample}.TD2.pep
    mv ${txmod_fasta}.TD2.cds ${sample}.TD2.cds
    mv ${txmod_fasta}.TD2.gff3 ${sample}.TD2.gff3
    mv ${txmod_fasta}.TD2.bed ${sample}.TD2.bed
    """
}
