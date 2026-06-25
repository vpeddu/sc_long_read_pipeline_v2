nextflow.enable.dsl=2

include { runFlexiplex } from './modules/A01_flexiplex.nf'
include { runMinimap2 } from './modules/A02_minimap2.nf'
include { runDeDup } from './modules/A03_deDup.nf'
include { runFlair; runTxRename } from './modules/A04_flair.nf'
include { runHtseq } from './modules/A05_htseq.nf'
include { runLongshot } from './modules/B01_longshot.nf'
include { runSQANTI3 } from './modules/C01_SQANTI3.nf'
include { runIsoSeQL } from './modules/C02_isoSeQL.nf'

workflow {
    params.input_dir = params.input_dir ?: '.'
    params.ref_dir = params.ref_dir 
    params.multiseq = params.get('multiseq', false) ? true : false
    params.seqtype = params.get('seqtype', '3_prime')
    if (!['3_prime','5_prime'].contains(params.seqtype)) {
        error("--seqtype must be either 3_prime or 5_prime")
    }
    params.avx2 = params.get('avx2', false) ?: ['fugu','iwashi','suzuki'].contains(System.getenv('HOSTNAME'))

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
        sample_data = Channel
            .fromPath("${params.input_dir}/*.fastq.gz", checkIfExists: true)
            .map { fastq_file ->
                def sample_name = fastq_file.name.replaceFirst(/\.(fastq|fq)\.gz$/, '')
                def barcode_file = file("${params.input_dir}/${sample_name}.barcodes.tsv")
                tuple(sample_name, fastq_file, barcode_file)
            }
    }
    // Reference directory channel
    ref_dir = Channel.fromPath("${params.ref_dir}/*", checkIfExists: true)
    ref_junc_bed = Channel.fromPath("${params.ref_dir}/genes.junc.bed", checkIfExists: true)
    ref_genome_fa = Channel.fromPath("${params.ref_dir}/genome.fa", checkIfExists: true)
    ref_genome_fai = Channel.fromPath("${params.ref_dir}/genome.fa.fai", checkIfExists: true)
    ref_genes_gtf = Channel.fromPath("${params.ref_dir}/genes.gtf", checkIfExists: true)
    // vep_data = Channel.fromPath("/mnt/ix1/Resources/VariantAnnotation/VEP/GRCh38_v104").map { file(it, type: 'dir') }    
    // Remove the hardcoded line and use params instead
    vep_data = Channel.fromPath(params.vep_cache).map { file(it, type: 'dir') }

    //TODO: REMOVE vep_data.view { "VEP Data: ${it}" }

    // Log found samples
    sample_data.view { n, f, b -> "Processing sample: ${n}" }

    // Duplicate the sample metadata channel to allow multiple downstream consumers
    def flex_input = sample_data
    def names_input = sample_data

    // Run flexiplex
    flex_output = runFlexiplex(
        sample_data,
        params.seqtype
        )
    A02_minimap = runMinimap2(flex_output,
        ref_junc_bed,
        ref_genome_fa,
        file("${projectDir}/bin/bam_tag_bc_umi.py")
        )
    A03_dedup = runDeDup(A02_minimap,
        file("${projectDir}/bin/collapse_barcodes.py"),
        ref_genes_gtf
        )
    
    A04_flair = runFlair(A03_dedup,
        ref_genome_fa,
        ref_genes_gtf
        )
        
    A04_txRename = runTxRename(A04_flair)
    
    A05_htseq = runHtseq(A03_dedup,
        ref_genes_gtf
        )

    B01_longshot = runLongshot(A03_dedup,
    vep_data,
    ref_genome_fa,
    ref_genome_fai)

    C01_sqanti3 = runSQANTI3(ref_genes_gtf,
    ref_genome_fa,
    A04_txRename)
    C02_isoSeQL = runIsoSeQL(C01_sqanti3)

    //C02_isoSeQL.subscribe { dir -> println "*** Pipeline complete: ${dir}" }
}
