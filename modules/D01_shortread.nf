process runSTARIndex {
    tag "D01_star_index"

    publishDir path: { "${params.output_dir}/D01_shortread/star_index" }, mode: 'symlink'

    input:
    path genome
    path gtf

    output:
    path "star_index"

    script:
    """
    mkdir star_index
    /opt/miniconda3/envs/sqanti3/bin/STAR --runMode genomeGenerate \
        --genomeDir star_index \
        --genomeFastaFiles ${genome} \
        --sjdbGTFfile ${gtf} \
        --sjdbOverhang ${params.star_sjdb_overhang} \
        --runThreadN ${task.cpus}
    """
}

process runSTARAlign {
    tag "D01_star_${sample}"

    publishDir path: { "${sample}/${params.output_dir}/D01_shortread/star" }, mode: 'symlink'

    input:
    tuple val(sample), path(r2)
    path star_index

    output:
    tuple val(sample), path("${sample}.SJ.out.tab")

    script:
    // Single-end on R2 only: these are droplet single-cell (10x Chromium)
    // reads, where R1 is a cell barcode+UMI read with no genomic content --
    // aligning it as a paired-end mate would corrupt junction discovery.
    // Only R2 (the actual cDNA read) carries alignable transcript sequence.
    """
    # --outSAMtype None: only junctions are needed here, skip alignment output.
    # --twopassMode Basic: standard best-practice for junction-discovery
    # sensitivity (re-maps using junctions discovered in a first pass).
    /opt/miniconda3/envs/sqanti3/bin/STAR --runMode alignReads \
        --genomeDir ${star_index} \
        --readFilesIn ${r2} \
        --readFilesCommand zcat \
        --runThreadN ${task.cpus} \
        --twopassMode Basic \
        --outSAMtype None \
        --outSJtype Standard \
        --outFileNamePrefix ${sample}.
    """
}

process runShortReadQuant {
    tag "D01_quant_${sample}"

    publishDir path: { "${sample}/${params.output_dir}/D01_shortread/quant" }, mode: 'symlink'

    input:
    tuple val(sample), path(txmod_fasta), path(r1), path(r2), val(chemistry), path(barcode_whitelist)

    output:
    tuple val(sample), path("${sample}_alevin_quant")

    script:
    // These short reads are droplet single-cell (10x Chromium) libraries, not
    // bulk paired-end: R1 is a cell barcode+UMI read, R2 is the cDNA read.
    // Quantify with salmon alevin (--rad --sketch, producing a RAD file --
    // --sketch requires --rad) + alevin-fry rather than plain `salmon quant`,
    // which would misinterpret R1 as a genuine paired-end mate. Requires
    // salmon <=1.10.2: 1.11+ deprecated `alevin` down to a stub pointing back
    // at 1.10.2, and salmon "2.x" (a full Rust rewrite) dropped it entirely
    // in favor of piscem+alevin-fry -- see bin/env_specs/salmon.explicit.txt.
    // --unfiltered-pl reuses this sample's own long-read barcode whitelist so
    // cells line up across both modalities.
    // chemistry encodes both barcode/UMI version and 3'/5' end (eg
    // 'chromiumV3_5p'): barcode/UMI geometry is identical between the 3' and
    // 5' kits of a given version, so only the trailing _3p/_5p suffix is
    // stripped for salmon alevin's own --chromium/--chromiumV3 flag; -l ISR
    // is correct for BOTH ends in the alevin-fry pipeline specifically (unlike
    // native alevin, orientation filtering happens in generate-permit-list,
    // not via libType -- see https://github.com/COMBINE-lab/alevin-fry/issues/118).
    // What DOES flip with the end is generate-permit-list's expected fragment
    // orientation: fw for 3', rc for 5' (R2 is the only mapped/biological
    // read here, and 5' R2s map to the reverse-complement strand).
    // Identity t2g map: quantify at ISOFORM resolution (flair's isoforms are
    // the unit of interest here), not the usual gene-level rollup.
    def alevin_chem = chemistry.replaceAll(/_(3p|5p)$/, '')
    def orientation = chemistry.endsWith('_5p') ? 'rc' : 'fw'
    """
    /opt/miniconda3/envs/salmon/bin/salmon index -t ${txmod_fasta} -i ${sample}.salmon_idx -k 31 -p ${task.cpus}

    awk '/^>/{id=substr(\$0,2); print id"\t"id}' ${txmod_fasta} > t2g.tsv

    /opt/miniconda3/envs/salmon/bin/salmon alevin \
        -l ISR \
        -i ${sample}.salmon_idx \
        -1 ${r1} -2 ${r2} \
        --${alevin_chem} \
        -p ${task.cpus} \
        --rad --sketch \
        -o ${sample}_map

    /opt/miniconda3/envs/salmon/bin/alevin-fry generate-permit-list \
        -d ${orientation} \
        -i ${sample}_map \
        -o ${sample}_quant \
        --unfiltered-pl ${barcode_whitelist}

    /opt/miniconda3/envs/salmon/bin/alevin-fry collate \
        -i ${sample}_quant \
        -r ${sample}_map \
        -t ${task.cpus}

    /opt/miniconda3/envs/salmon/bin/alevin-fry quant \
        -i ${sample}_quant \
        -m t2g.tsv \
        -r cr-like \
        -o ${sample}_alevin_quant \
        -t ${task.cpus} \
        --use-mtx
    """
}
