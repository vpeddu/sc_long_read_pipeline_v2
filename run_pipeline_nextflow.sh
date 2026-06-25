#!/bin/bash
# Nextflow wrapper for the long-read isoform analysis pipeline

set -e

WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN_NF="$WORK_DIR/main.nf"

if [ $# -eq 0 ]; then
    cat << EOF
Usage: $(basename "$0") --input_dir /path/to/samples [OPTIONS]

Process all samples in a directory through the long-read isoform analysis pipeline.

REQUIRED:
    --input_dir DIR              Directory containing fastq.gz and barcodes.tsv files

OPTIONS:
    --multiseq                   Enable multi-seq processing (shorter adapters)
    --resume                     Resume from last checkpoint
    --dry-run                    Show what would be executed
    --container                  Use Singularity container (default: true)
    --no-container               Use local conda environments instead
    --work_dir DIR               Nextflow work directory (default: work)

EXAMPLES:
    $(basename "$0") --input_dir /path/to/samples
    $(basename "$0") --input_dir /path/to/samples --resume
    $(basename "$0") --input_dir /path/to/samples --no-container

EOF
    exit 0
fi

INPUT_DIR=""
USE_CONTAINER="true"
NF_FLAGS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --input_dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --multiseq)
            NF_FLAGS="$NF_FLAGS --multiseq true"
            shift
            ;;
        --resume)
            NF_FLAGS="$NF_FLAGS -resume"
            shift
            ;;
        --dry-run)
            NF_FLAGS="$NF_FLAGS -n"
            shift
            ;;
        --container)
            USE_CONTAINER="true"
            shift
            ;;
        --no-container)
            USE_CONTAINER="false"
            shift
            ;;
        --work_dir)
            NF_FLAGS="$NF_FLAGS -w $2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [ -z "$INPUT_DIR" ]; then
    echo "Error: --input_dir is required"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Directory not found: $INPUT_DIR"
    exit 1
fi

echo "Starting pipeline for samples in: $INPUT_DIR"
echo "Using container: $USE_CONTAINER"
echo ""

nextflow run "$MAIN_NF" \
    -Duser.country=US \
    -Duser.language=en \
    --input_dir "$INPUT_DIR" \
    --use_container "$USE_CONTAINER" \
    $NF_FLAGS
