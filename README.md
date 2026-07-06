**SC_Long_Read Pipeline**

This repository contains a Nextflow-based long-read single-cell processing pipeline. It orchestrates read alignment, barcode/UMI handling, isoform collapse and classification, expression counting, variant calling, and QC using modular Nextflow workflows and helper scripts.

**Requirements**
- **Nextflow**: pipeline entry is [main.nf](main.nf#L1).
- **Java 8+** and **Singularity** (optional but recommended for reproducible runs).
- Standard bioinformatics tools used by modules (see `modules/` and `bin/` for specifics).

**Profiles & Containers**
- Use `-profile singularity` to run with Singularity (default).
- Container build script is `bin/Singularity.def`. Build with `bin/build_container.sh`
- Current container: `ghcr.io/vpeddu/sc_long_pipeline:v1.2.0`, which added a `seurat` conda env (Seurat + tximport, R 4.3) for `A06_seurat.nf` alongside the existing `long_reads`/`flair`/`vep`/`picardtools`/`isoSeQL`/`bioconductor`/`sqanti3` environments.
- `D01_shortread.nf` (`--pairedshortread`) reuses STAR and kallisto, both already present in the `sqanti3` env — no container rebuild is required for this feature.

**Pipeline Layout**
- `main.nf`: pipeline entry that composes module workflows.
- `nextflow.config`: configuration and profiles used by the pipeline.
- `modules/`: Nextflow module files for each major step (e.g., `A01_flexiplex.nf`, `A02_minimap2.nf`).
- `bin/` and `original_scripts/`: helper scripts and legacy tools invoked by workflow tasks.

**Containers & Scratch storage**
- **Container location:** the container artifacts and build helpers live in the `container/` directory (for example `container/sc_long_read_container` and `container/Dockerfile` or `container/build.sh`). Build or push these images to your preferred registry if you want cluster-wide access.
- **Default Sherlock scratch behavior:** when using the `sherlock` profile, the pipeline is expected to use `$GROUP_SCRATCH` for permanent shared files and node-local `$L_SCRATCH` for temporary per-job working files.
- **Local profile note:** the `local` profile does not set a default `workDir`. Pass `--work-dir <path>` explicitly if you want to control where work and intermediate files are written.

**Common Modules**
- `A01_flexiplex.nf`: demultiplexing/flexiplex-related steps.
- `A02_minimap2.nf`: alignment with minimap2.
- `A03_deDup.nf`: deduplication/UMI collapsing.
- `A04_flair.nf`: isoform analysis with FLAIR (`flair transcriptome`), transcript ID renaming, and per-transcript count matrix generation.
- `A05_htseq.nf`: per-cell gene expression counting.
- `A06_seurat.nf`: builds a per-sample Seurat object (`RNA` + `ISO` assays) combining gene- and isoform-level counts, SQANTI3 classification, genomic ranges, and read-count metrics (see **Outputs** below).
- `A07_metrics.nf`: raw/pre-dedup/post-dedup read counts per sample.
- `C01_SQANTI3.nf` and related: transcript classification and QC.
- `D01_shortread.nf`: optional paired short-read integration — see `--pairedshortread` below.

**Configuration options**
- Edit [nextflow.config](nextflow.config#L1) or pass overrides on the command line, e.g. `--threads 8` or `--genome GRCh38`.
- `--flair_split_by_chrom` (`true`/`false`, default `true`): when `true`, `flair transcriptome` runs once per chromosome in parallel instead of once genome-wide, which is substantially faster/lower-memory on large single-cell transcriptomes. Per-chromosome outputs are concatenated back into the usual per-sample files afterward.
- `--keep_intergenic` (`true`/`false`, default `false`): flair calls isoforms at novel loci that don't overlap any annotated gene (assigned a raw coordinate-based id instead of an Ensembl gene id). By default these are dropped from the final GTF/counts/Seurat outputs. Set `true` to keep them — each distinct intergenic locus gets a unique `novel_intergenic_NNN` label, `A05_htseq` switches to quantifying against the per-sample flair-derived transcriptome (instead of just the static reference annotation) so these loci get their own gene-level counts, and they flow through into the Seurat object like any other gene/isoform.
- `--pairedshortread <csv>`: opt-in paired short-read integration. CSV has 3 columns, **no header**, `#`-comment lines allowed: `long_read_sample,short_read_R1,short_read_R2`. The first column must exactly match a long-read sample name as derived from its fastq filename. Samples don't need to line up 1:1 — long-read samples with no row in the CSV proceed normally, and rows whose long-read sample isn't part of the current run are silently ignored (there's no long-read transcriptome to correct/quantify against). When set:
  - A one-time whole-genome STAR index is built, and each matched sample's short reads are aligned with STAR (`--twopassMode Basic`, junctions only, no BAM) to produce `SJ.out.tab`.
  - That `SJ.out.tab` is fed into `flair transcriptome` via `--junction_tab` for that sample (both the genome-wide and per-chromosome variants), improving splice-site accuracy with orthogonal short-read evidence. Unmatched samples run exactly as before (no short-read flags added).
  - The matched sample's short reads are then quantified (via `kallisto index`+`quant`, already present in the `sqanti3` env — no container changes needed) against that sample's own renamed/filtered novel transcriptome FASTA (see **Outputs** below), giving an orthogonal short-read expression estimate for the long-read-derived isoforms.
  - Tunable: `--star_sjdb_overhang` (default `100`; STAR's ideal is `read_length - 1`, override if your short reads aren't ~101bp) and `--shortread_junction_support` (default `1`, matches flair's own `--junction_support` default).

`--flair_split_by_chrom` and `--keep_intergenic` accept the literal strings `true`/`false` on the command line (e.g. `--keep_intergenic true`); pass them explicitly rather than as bare flags.

**Run examples**
Below are example commands for running the pipeline locally and on the Sherlock cluster. Replace the example paths with your input and reference locations.

Local (example using a relative path to the pipeline):

```
nextflow run ../main.nf \
  --input_dir ../path/to/input/ \
  --seqtype '3_prime' \
  --multiseq false \
  --ref_dir /path/to/resources/ \
  -resume \
  -with-tower \
  -profile "local"
```

Example `run.sh` for Sherlock (place in your run directory and edit paths as needed):

```
[vikasp1@sh04-ln04 login ~/sc_long_read_testing/sc_long_read]$ cat run.sh
nextflow run main.nf \
--input_dir /scratch/groups/hanleeji/vikas/input/ \
--seqtype '3_prime' \
--multiseq false \
--ref_dir /scratch/groups/hanleeji/vikas/reference/ \
--vep_cache /scratch/groups/hanleeji/vikas/vep_dir/GRCh38_v104/ \
--profile sherlock \
-with-tower
```

**Outputs**
- Results are written to the Nextflow work directory and the pipeline `results/` folder (or as configured in `nextflow.config`).
- Key outputs include aligned BAMs, collapsed isoform GTFs, expression count tables, SQANTI3 QC reports, and variant VCFs.
- `A06_seurat`: one `<sample>.seurat.rds` Seurat object per sample, with an `RNA` assay (gene-level counts) and an `ISO` assay (isoform-level counts). Cells are tagged with `orig.ident` and barcode-prefixed by sample so multiple samples' objects can be merged safely. Feature-level metadata includes genomic ranges (`chrom`/`start`/`end`/`strand`, both assays) and, on the `ISO` assay, SQANTI3 classification columns (`structural_category`, `associated_gene`, `length`, `exons`, coverage/coding/NMD/filter columns). Per-cell metadata includes `raw_reads`, `predup_reads`, `postdup_reads`, and `dedup_rate` (also available as a compact list in `@misc$sample_metrics`).
- `A04_txRename`: alongside the existing renamed GTF/xref/matrix outputs, always also produces `<sample>.flair.collapse.isoforms.txmod.fa` — the raw flair isoform FASTA renamed and filtered to match the same transcript set/IDs as the renamed GTF (this is "the novel long-read transcriptome" used for short-read quantification below).
- `D01_shortread` (only produced for samples matched via `--pairedshortread`): `<sample>.SJ.out.tab` (STAR-derived splice junctions fed into flair) and `<sample>_kallisto/` (standard kallisto quant output — `abundance.tsv`/`abundance.h5`/`run_info.json` — quantifying that sample's short reads against its own `txmod.fa`).

**Troubleshooting**
- If a step fails, inspect the Nextflow `work/` directory for task logs and the `trace.txt` and `report.html` files generated with `-with-trace`/`-with-report`.
- For container problems, verify Docker/Singularity is installed and the correct profile is used.

**Contributing**
- Add/modify modules in `modules/` and helper scripts in `bin/`.
- Update `SETUP_GUIDE.md` and `SETUP_CHECKLIST.md` with environment-specific instructions when adding new dependencies.

**References**
- Pipeline scripts: [main.nf](main.nf#L1) and [nextflow.config](nextflow.config#L1).
- Helper scripts: [bin/](bin/) and [original_scripts/](original_scripts/).

