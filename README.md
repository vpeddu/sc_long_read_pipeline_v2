**SC_Long_Read Pipeline**

This repository contains a Nextflow-based long-read single-cell processing pipeline. It orchestrates read alignment, barcode/UMI handling, isoform collapse and classification, expression counting, variant calling, and QC using modular Nextflow workflows and helper scripts.

**Requirements**
- **Nextflow**: pipeline entry is [main.nf](main.nf#L1).
- **Java 8+** and **Singularity** (optional but recommended for reproducible runs).
- Standard bioinformatics tools used by modules (see `modules/` and `bin/` for specifics).

**Profiles & Containers**
- Use `-profile singularity` to run with Singularity (default).
- Container build script is `bin/Singularity.def`. Build with `bin/build_container.sh`

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
- `A04_flair.nf`: isoform analysis with FLAIR.
- `C01_SQANTI3.nf` and related: transcript classification and QC.

**Configuration options**
- Edit [nextflow.config](nextflow.config#L1) or pass overrides on the command line, e.g. `--threads 8` or `--genome GRCh38`.

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

**Troubleshooting**
- If a step fails, inspect the Nextflow `work/` directory for task logs and the `trace.txt` and `report.html` files generated with `-with-trace`/`-with-report`.
- For container problems, verify Docker/Singularity is installed and the correct profile is used.

**Contributing**
- Add/modify modules in `modules/` and helper scripts in `bin/`.
- Update `SETUP_GUIDE.md` and `SETUP_CHECKLIST.md` with environment-specific instructions when adding new dependencies.

**References**
- Pipeline scripts: [main.nf](main.nf#L1) and [nextflow.config](nextflow.config#L1).
- Helper scripts: [bin/](bin/) and [original_scripts/](original_scripts/).

