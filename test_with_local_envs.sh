#!/bin/bash
# Test pipeline using local conda environments instead of the container

set -e

WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN_NF="$WORK_DIR/main.nf"

if [ $# -eq 0 ]; then
    echo "Usage: $(basename "$0") --input_dir /path/to/samples [--resume] [--dry-run] [additional nextflow args]"
    exit 1
fi

# Check for required local conda environments
REQUIRED_ENVS=("long_reads" "flair" "sqanti3" "vep")
MISSING_ENVS=()

for env in "${REQUIRED_ENVS[@]}"; do
    if ! conda env list | grep -q "^${env} "; then
        MISSING_ENVS+=("$env")
    fi
done

if [ ${#MISSING_ENVS[@]} -gt 0 ]; then
    echo "Error: Missing required conda environments:"
    for env in "${MISSING_ENVS[@]}"; do
        echo "  - $env"
    done
    exit 1
fi

# Run Nextflow pipeline with use_container=false
nextflow run "$MAIN_NF" \
    -Duser.country=US \
    -Duser.language=en \
    --use_container false \
    "$@"
