nextflow.enable.dsl=2

include { splitFastq; runFlexiplex } from './modules/A01_flexiplex.nf'
include { runMinimap2; mergeBam } from './modules/A02_minimap2.nf'
include { runDeDup } from './modules/A03_deDup.nf'
include { runFlair; runFlairByChrom; mergeFlairChroms; runTxRename } from './modules/A04_flair.nf'
include { runHtseq; runChrmGeneSubset; mergeChrmReferenceGenes } from './modules/A05_htseq.nf'
include { runBuildSeurat } from './modules/A06_seurat.nf'
include { runSampleMetrics } from './modules/A07_metrics.nf'
include { runLongshot } from './modules/B01_longshot.nf'
include { runSQANTI3 } from './modules/C01_SQANTI3.nf'
include { runIsoSeQL } from './modules/C02_isoSeQL.nf'
include { runSTARIndex; runSTARAlign; runShortReadQuant } from './modules/D01_shortread.nf'

workflow {
    params.input_dir = params.input_dir ?: '.'
    params.ref_dir = params.ref_dir 
    // Booleans/ints are read and validated/cast at each point of use below
    // rather than normalized here via reassignment: reassigning params.X inside
    // the workflow body is unreliable under the strict-syntax parser (silently
    // dropped once a conflicting config-level default exists, e.g. multiseq's
    // `false` in nextflow.config), and CLI values arrive as raw strings, where
    // Groovy's truthiness treats even the string "false" as true.
    params.seqtype = params.get('seqtype', '3_prime')
    if (!['3_prime','5_prime'].contains(params.seqtype)) {
        error("--seqtype must be either 3_prime or 5_prime")
    }
    params.avx2 = params.avx2 ?: ['fugu','iwashi','suzuki'].contains(System.getenv('HOSTNAME'))
    if ((params.minimap2_split as int) < 1) {
        error("--minimap2_split must be a positive integer")
    }

    if (params.input_list) {
        sample_data = Channel
            .fromPath(params.input_list, checkIfExists: true)
            .splitText()
            .map { it.trim() }
            .filter { it && !it.startsWith('#') }
            .map { dir_path ->
                def fastq = file("${dir_path}/*.fastq.gz")[0]
                def barcode = file("${dir_path}/*.barcodes.tsv")[0]
                def sample_name = fastq.name.replaceFirst(/\.(fastq|fq)\.gz$/, '')
                tuple(sample_name, fastq, barcode)
            }
    } else {
        // input_dir holds a single sample: pair whichever fastq and barcode file
        // are present, regardless of naming (same convention as input_list)
        sample_data = Channel
            .fromPath(params.input_dir, checkIfExists: true)
            .map { dir_path ->
                def fastq = file("${dir_path}/*.fastq.gz")[0]
                def barcode = file("${dir_path}/*.barcodes.tsv")[0]
                def sample_name = fastq.name.replaceFirst(/\.(fastq|fq)\.gz$/, '')
                tuple(sample_name, fastq, barcode)
            }
    }
    // Reference directory channel
    ref_dir = Channel.fromPath("${params.ref_dir}/*", checkIfExists: true)
    ref_junc_bed = Channel.fromPath("${params.ref_dir}/genes.junc.bed", checkIfExists: true)
    ref_genome_fa = Channel.fromPath("${params.ref_dir}/genome.fa", checkIfExists: true)
    ref_genome_fai = Channel.fromPath("${params.ref_dir}/genome.fa.fai", checkIfExists: true)
    ref_genes_gtf = Channel.fromPath("${params.ref_dir}/genes.gtf", checkIfExists: true)
    ref_features_tsv = Channel.fromPath("${params.ref_dir}/features_gex.tsv", checkIfExists: true)
    // vep_data = Channel.fromPath("/mnt/ix1/Resources/VariantAnnotation/VEP/GRCh38_v104").map { file(it, type: 'dir') }    
    // Remove the hardcoded line and use params instead
    vep_data = Channel.fromPath(params.vep_cache).map { file(it, type: 'dir') }

    //TODO: REMOVE vep_data.view { "VEP Data: ${it}" }

    // Log found samples
    sample_data.view { n, f, b -> "Processing sample: ${n}" }

    // Optional paired short-read data: headerless CSV (long_read_sample,
    // short_read_R1,short_read_R2,chemistry), '#'-comment lines skipped --
    // same headerless-list convention as --input_list. These short reads are
    // droplet single-cell (10x Chromium) libraries of the SAME cells as the
    // long-read sample, not bulk paired-end -- R1 is a cell barcode+UMI read,
    // R2 is the cDNA read. chemistry encodes both the 10x kit version
    // ('chromium' = v2, 'chromiumV3' = v3 -- salmon alevin's own barcode/UMI
    // geometry flag names, identical between a version's 3' and 5' kits) and
    // which end was sequenced (_3p/_5p suffix), since that changes the
    // expected read orientation used downstream in runShortReadQuant.
    // Long-read samples with no short-read pair still proceed (just without
    // junction correction/quant); short-read pairs whose long-read sample
    // isn't part of this run have no novel transcriptome to correct/quantify
    // against, so they're dropped here (before STAR ever runs on them).
    def VALID_SR_CHEMISTRIES = ['chromium_3p', 'chromium_5p', 'chromiumV3_3p', 'chromiumV3_5p']
    if (params.pairedshortread) {
        def sr_rows = file(params.pairedshortread, checkIfExists: true).readLines()
            .collect { it.trim() }
            .findAll { it && !it.startsWith('#') }
            .collect { it.split(',') }
        def sr_names = sr_rows.collect { it[0].trim() }
        def sr_dups = sr_names.findAll { n -> sr_names.count(n) > 1 }.unique()
        if (sr_dups) {
            error("Duplicate long-read sample name(s) in --pairedshortread: ${sr_dups}")
        }
        def sr_bad_chem = sr_rows.collect { it[3].trim() }.findAll { !VALID_SR_CHEMISTRIES.contains(it) }.unique()
        if (sr_bad_chem) {
            error("--pairedshortread chemistry column must be one of ${VALID_SR_CHEMISTRIES}; got: ${sr_bad_chem}")
        }

        def shortread_data_raw = Channel.fromList(sr_rows.collect { row ->
            tuple(row[0].trim(), file(row[1].trim(), checkIfExists: true), file(row[2].trim(), checkIfExists: true), row[3].trim())
        })
        // Plain .join() against sample_data itself, rather than a
        // .collect()-based membership filter: .collect() is a full-channel
        // barrier that can't emit until ALL of --input_list has been
        // enumerated, which stalled runSTARAlign behind that full
        // enumeration even though runSTARIndex has no such dependency and
        // finishes independently. .join() lets a matched short-read sample
        // flow through as soon as ITS OWN long-read counterpart appears in
        // sample_data. Unmatched rows on either side are silently dropped,
        // same as before (a short-read-only row has no long-read
        // transcriptome to correct/quantify against).
        shortread_data = shortread_data_raw.join(sample_data)
            .map { s, r1, r2, chem, _fastq, _barcode -> tuple(s, r1, r2, chem) }
        shortread_data.view { s, r1, r2, chem -> "Paired short-read sample: ${s} (${chem})" }
    } else {
        shortread_data = Channel.empty()
    }

    // Duplicate the sample metadata channel to allow multiple downstream consumers
    def flex_input = sample_data
    def names_input = sample_data

    // Split each sample's raw fastq into params.minimap2_split chunks before flexiplex
    split_input = splitFastq(sample_data, params.minimap2_split as int)

    // One (sample, chunk_fastq, bc_tsv) tuple per chunk -> run flexiplex in parallel
    flex_output = runFlexiplex(split_input.transpose(),
        params.seqtype
        )

    // One (sample, chunk_fastq) tuple per chunk -> run minimap2 in parallel
    minimap2_chunks = runMinimap2(flex_output,
        ref_junc_bed.first(),
        ref_genome_fa.first()
        )

    // Recombine all chunk BAMs per sample only after minimap2
    A02_minimap = mergeBam(minimap2_chunks.groupTuple())

    A03_dedup = runDeDup(A02_minimap,
        file("${projectDir}/bin/collapse_barcodes.py"),
        ref_genes_gtf
        )

    // Build the whole-genome STAR index (and align each paired short-read
    // sample against it) only when --pairedshortread is actually in use --
    // this is a ~30GB, several-GB-of-disk step, not worth paying for otherwise.
    if (params.pairedshortread) {
        star_index = runSTARIndex(ref_genome_fa, ref_genes_gtf)
        // R2 only: R1 is a cell barcode+UMI read with no genomic content --
        // see runSTARAlign's single-end rationale in modules/D01_shortread.nf.
        star_align_input = shortread_data.map { s, r1, r2, chem -> tuple(s, r2) }
        star_sj = runSTARAlign(star_align_input, star_index.first())
    } else {
        star_sj = Channel.empty()
    }

    // Pair each long-read sample's dedup bam with its optional short-read STAR
    // junctions. remainder:true is safe here specifically because shortread_data
    // was already filtered (above) to samples known to sample_data, so star_sj
    // can never contain a key absent from A03_dedup -- the only padding that
    // can occur is a long-read sample with no short-read match (intended).
    // sj_tab ?: [] uses Nextflow's "no optional file" idiom (an empty list is
    // Groovy-falsy, a bound path is truthy) instead of a sentinel file.
    flair_input = A03_dedup.join(star_sj, remainder: true)
        .map { sample, bam, bai, fastq, sj_tab -> tuple(sample, bam, bai, fastq, sj_tab ?: []) }

    // params.flair_split_by_chrom may be a raw CLI string ("true"/"false") --
    // any non-empty Groovy string (including "false") is truthy, so compare
    // the string value explicitly rather than relying on plain truthiness.
    if (!params.flair_split_by_chrom.toString().equalsIgnoreCase('false')) {
        // Group contigs for parallel per-chromosome flair runs: each "primary"
        // chromosome (no '.' in its name, e.g. chr1..chrY/chrM for GRCh38) gets
        // its own task; everything else (alt/decoy/unplaced scaffolds, which
        // typically carry very few or no aligned reads) is bundled into one
        // "other_contigs" task so it doesn't spawn hundreds of near-empty jobs.
        chrom_groups = ref_genome_fai
            .splitCsv(sep: '\t')
            .map { row -> row[0] }
            .collect()
            .flatMap { contigs ->
                def primary = contigs.findAll { !it.contains('.') }
                // chrM's read depth is vastly disproportionate to its size,
                // making flair transcript calling on it a major outlier in
                // runtime; drop it entirely here (not bundled into
                // other_contigs, which would still pay that cost) rather than
                // only filtering it out of the final output afterward.
                if (params.exclude_chrm.toString().equalsIgnoreCase('true')) {
                    primary = primary.findAll { it != 'chrM' }
                }
                def other = contigs.findAll { it.contains('.') }
                def groups = primary.collect { [it, [it]] }
                if (other) {
                    groups << ['other_contigs', other]
                }
                groups
            }

        flair_chrom_input = flair_input.combine(chrom_groups)
        A04_flair_chunks = runFlairByChrom(flair_chrom_input,
            ref_genome_fa.first(),
            ref_genes_gtf.first()
            )
        A04_flair = mergeFlairChroms(A04_flair_chunks.groupTuple())
    } else {
        A04_flair = runFlair(flair_input,
            ref_genome_fa,
            ref_genes_gtf
            )
    }


    A04_txRename_out = runTxRename(A04_flair, ref_features_tsv.first())
    A04_txRename = A04_txRename_out.renamed_gtf
    A04_iso_counts = A04_txRename_out.iso_counts

    // params.keep_intergenic may be a raw CLI string ("true"/"false") -- compare
    // the string value explicitly rather than relying on Groovy truthiness.
    if (params.keep_intergenic.toString().equalsIgnoreCase('true')) {
        // Quantify against the novel (flair-derived) transcriptome instead of the
        // static reference annotation, so novel_intergenic_NNN loci (assigned in
        // runTxRename) get their own htseq counts. Each sample's dedup bam is
        // paired with ITS OWN per-sample txmod gtf (not the shared reference gtf)
        // via a sample-keyed join. The txmod gtf only carries a "gene_id"
        // attribute (no separate "gene_name"), so --idattr switches accordingly.
        htseq_gtf_by_sample = A04_txRename.map { sample, txmod_gtf, read_map, isoform_cells, transcript_xref, txmod_fasta -> tuple(sample, txmod_gtf) }
        // --exclude_chrm drops chrM from the flair-derived gtf entirely (see
        // chrom_groups above), which would otherwise silently zero out chrM
        // gene counts here too -- backfill them from the reference annotation
        // so they aren't lost just because flair itself skipped chrM. Only
        // relevant when exclude_chrm actually removed them; when it's false,
        // flair already called chrM transcripts and txmod_gtf already has them.
        if (params.exclude_chrm.toString().equalsIgnoreCase('true')) {
            chrm_genes_gtf = runChrmGeneSubset(ref_genes_gtf)
            htseq_gtf_by_sample = mergeChrmReferenceGenes(htseq_gtf_by_sample, chrm_genes_gtf.first())
        }
        htseq_input = A03_dedup.join(htseq_gtf_by_sample)
        A05_htseq = runHtseq(htseq_input, 'gene_id')
    } else {
        htseq_input = A03_dedup.combine(ref_genes_gtf)
        A05_htseq = runHtseq(htseq_input, 'gene_name')
    }

    B01_longshot = runLongshot(A03_dedup,
    vep_data,
    ref_genome_fa,
    ref_genome_fai)

    // runSQANTI3's input signature predates the txmod fasta added to
    // A04_txRename's emit -- reshape back down to the 5-tuple it expects
    // rather than touching the (unrelated) SQANTI3 module itself.
    sqanti_input = A04_txRename.map { sample, txmod_gtf, read_map, isoform_cells, transcript_xref, txmod_fasta ->
        tuple(sample, txmod_gtf, read_map, isoform_cells, transcript_xref)
    }
    C01_sqanti3 = runSQANTI3(ref_genes_gtf,
    ref_genome_fa,
    sqanti_input)
    C02_isoSeQL = runIsoSeQL(C01_sqanti3)

    // Quantify each paired sample's short reads against its OWN novel
    // long-read-derived transcriptome (the renamed/filtered isoform fasta).
    // Plain inner join: a short-read-only sample has no transcriptome to
    // quantify against, and a long-read-only sample has no short reads.
    // The barcode whitelist is this sample's own long-read cell barcode list
    // (same one used for flexiplex demux) -- reusing it as alevin-fry's
    // --unfiltered-pl keeps cell identities consistent across both modalities.
    txmod_fasta_by_sample = A04_txRename.map { sample, txmod_gtf, read_map, isoform_cells, transcript_xref, txmod_fasta -> tuple(sample, txmod_fasta) }
    barcode_by_sample = sample_data.map { sample, fastq, barcode -> tuple(sample, barcode) }
    quant_input = txmod_fasta_by_sample.join(shortread_data).join(barcode_by_sample)
    D01_quant = runShortReadQuant(quant_input)

    // Per-sample raw/pre-dedup/post-dedup read counts, joined on sample name
    // from the raw input fastq, the pre-dedup merged bam, and the post-dedup bam.
    metrics_input = sample_data.join(A02_minimap).join(A03_dedup)
        .map { sample, fastq, barcode, bam, bai, dedup_bam, dedup_bai, dedup_fastq ->
            tuple(sample, fastq, bam, dedup_bam)
        }
    A07_metrics = runSampleMetrics(metrics_input)

    // Combine per-sample gene-level (htseq) and isoform-level (rdmap2counts, via
    // runTxRename) count matrices, SQANTI3 classification/genomic ranges, and
    // read-count metrics into one Seurat object per sample, all joined on sample
    // name. Depending on C01_sqanti3 here means Seurat object construction can't
    // start until the whole SQANTI3 branch finishes for that sample (previously
    // it only needed A04/A05, running in parallel with SQANTI3/longshot).
    seurat_input = A05_htseq.join(A04_iso_counts).join(C01_sqanti3).join(A07_metrics)
        .map { sample, gene_matrix, gene_features, gene_barcodes, annot_bam,
               iso_matrix, iso_features, iso_barcodes,
               sqanti_gtf, sqanti_fasta, sqanti_classif, sqanti_xref, sqanti_genepred, tx_dups_xref,
               sample_metrics ->
            tuple(sample, gene_matrix, gene_features, gene_barcodes,
                  iso_matrix, iso_features, iso_barcodes,
                  sqanti_gtf, sqanti_classif, sample_metrics)
        }
    A06_seurat = runBuildSeurat(seurat_input, ref_genes_gtf.first())

    //C02_isoSeQL.subscribe { dir -> println "*** Pipeline complete: ${dir}" }
}
