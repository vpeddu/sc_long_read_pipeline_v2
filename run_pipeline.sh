#!/bin/bash
# Helper script to run the Nextflow pipeline with Singularity container

set -e

WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER="$WORK_DIR/Singularity.sif"
MAIN_NF="$WORK_DIR/main.nf"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check prerequisites
if ! command -v nextflow &> /dev/null; then
    print_error "Nextflow is not installed or not in PATH"
    echo "Install Nextflow from: https://www.nextflow.io/"
    exit 1
fi

if ! command -v singularity &> /dev/null; then
    print_error "Singularity is not installed or not in PATH"
    echo "Install Singularity from: https://sylabs.io/"
    exit 1
fi

if [ ! -f "$CONTAINER" ]; then
    print_error "Container file not found: $CONTAINER"
    echo ""
    echo "Build the container first:"
    echo "  bash $WORK_DIR/build_container.sh"
    exit 1
fi

if [ ! -f "$MAIN_NF" ]; then
    print_error "Pipeline file not found: $MAIN_NF"
    exit 1
fi

# Display usage if no arguments
if [ $# -eq 0 ]; then
    cat << EOF
Usage: $(basename "$0") [OPTION]... --sample SAMPLE_NAME

Options:
    --sample SAMPLE_NAME          Sample name (required)
    --multiseq                    Enable multi-seq processing
    --avx2                        Enable AVX2 optimizations
    --code_dir PATH               Path to code directory (default: ./bin)
    --work_dir PATH               Nextflow work directory (default: ./work)
    --resume                      Resume from last checkpoint
    --dry-run                     Show what would be executed
    -h, --help                    Show this help message

Examples:
    $(basename "$0") --sample sample_001
    $(basename "$0") --sample sample_001 --multiseq
    $(basename "$0") --sample sample_001 --resume
    $(basename "$0") --sample sample_001 --dry-run

Environment:
    Container: $CONTAINER
    Pipeline: $MAIN_NF
    Work dir: ./work (can be overridden with --work_dir)

EOF
    exit 0
fi

# Parse arguments
NEXTFLOW_ARGS=""
SAMPLE_NAME=""
WORK_DIR_NEXTFLOW="work"
RESUME_FLAG=""
DRY_RUN_FLAG=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)
            SAMPLE_NAME="$2"
            NEXTFLOW_ARGS="$NEXTFLOW_ARGS --sample $2"
            shift 2
            ;;
        --multiseq)
            NEXTFLOW_ARGS="$NEXTFLOW_ARGS --multiseq"
            shift
            ;;
        --avx2)
            NEXTFLOW_ARGS="$NEXTFLOW_ARGS --avx2"
            shift
            ;;
        --code_dir)
            NEXTFLOW_ARGS="$NEXTFLOW_ARGS --code_dir $2"
            shift 2
            ;;
        --work_dir)
            WORK_DIR_NEXTFLOW="$2"
            shift 2
            ;;
        --resume)
            RESUME_FLAG="-resume"
            shift
            ;;
        --dry-run)
            DRY_RUN_FLAG="-n"
            shift
            ;;
        -h|--help)
            bash "$0"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate sample name
if [ -z "$SAMPLE_NAME" ]; then
    print_error "Sample name is required"
    echo "Use: $(basename "$0") --sample SAMPLE_NAME"
    exit 1
fi

# Create work directory if it doesn't exist
mkdir -p "$WORK_DIR_NEXTFLOW"

# Display pipeline information
echo ""
print_info "Long-Read Isoform Analysis Pipeline"
echo "=========================================="
echo "Sample name: $SAMPLE_NAME"
echo "Container: $CONTAINER"
echo "Work dir: $WORK_DIR_NEXTFLOW"
echo "Resume: ${RESUME_FLAG:-disabled}"
echo "Dry run: ${DRY_RUN_FLAG:-disabled}"
echo "=========================================="
echo ""

# Run Nextflow pipeline
print_info "Starting pipeline execution..."
echo ""

nextflow run "$MAIN_NF" \
    -w "$WORK_DIR_NEXTFLOW" \
    $RESUME_FLAG \
    $DRY_RUN_FLAG \
    $NEXTFLOW_ARGS

EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    print_info "Pipeline completed successfully!"
else
    print_error "Pipeline failed with exit code $EXIT_CODE"
fi

exit $EXIT_CODE
