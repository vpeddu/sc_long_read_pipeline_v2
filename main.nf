nextflow.enable.dsl=2

include { splitFastq; runFlexiplex } from './modules/A01_flexiplex.nf'
include { runMinimap2; mergeBam } from './modules/A02_minimap2.nf'
include { runDeDup } from './modules/A03_deDup.nf'
include { runFlair; runFlairByChrom; mergeFlairChroms; runTxRename } from './modules/A04_flair.nf'
include { runHtseq } from './modules/A05_htseq.nf'
include { runBuildSeurat } from './modules/A06_seurat.nf'
include { runSampleMetrics } from './modules/A07_metrics.nf'
include { runLongshot } from './modules/B01_longshot.nf'
include { runSQANTI3 } from './modules/C01_SQANTI3.nf'
include { runIsoSeQL } from './modules/C02_isoSeQL.nf'

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
                def other = contigs.findAll { it.contains('.') }
                def groups = primary.collect { [it, [it]] }
                if (other) {
                    groups << ['other_contigs', other]
                }
                groups
            }

        flair_chrom_input = A03_dedup.combine(chrom_groups)
        A04_flair_chunks = runFlairByChrom(flair_chrom_input,
            ref_genome_fa.first(),
            ref_genes_gtf.first()
            )
        A04_flair = mergeFlairChroms(A04_flair_chunks.groupTuple())
    } else {
        A04_flair = runFlair(A03_dedup,
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
        htseq_gtf_by_sample = A04_txRename.map { sample, txmod_gtf, read_map, isoform_cells, transcript_xref -> tuple(sample, txmod_gtf) }
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

    C01_sqanti3 = runSQANTI3(ref_genes_gtf,
    ref_genome_fa,
    A04_txRename)
    C02_isoSeQL = runIsoSeQL(C01_sqanti3)

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
