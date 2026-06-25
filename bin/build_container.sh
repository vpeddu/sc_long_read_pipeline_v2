#!/bin/bash
# Quick build script for Singularity container
# This script helps build the Singularity container for the pipeline

set -e

CONTAINER_NAME="Singularity.sif"
DEF_FILE="Singularity.def"
WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=========================================="
echo "Long-Read Pipeline Singularity Build"
echo "=========================================="
echo ""
echo "Working directory: $WORK_DIR"
echo "Definition file: $DEF_FILE"
echo "Output container: $CONTAINER_NAME"
echo ""

# Check if Singularity is installed
if ! command -v singularity &> /dev/null; then
    echo "ERROR: Singularity is not installed or not in PATH"
    echo "Please install Singularity 3.0 or higher"
    exit 1
fi

echo "Singularity version:"
singularity --version
echo ""

# Check for definition file
if [ ! -f "$WORK_DIR/$DEF_FILE" ]; then
    echo "ERROR: $DEF_FILE not found in $WORK_DIR"
    exit 1
fi


# Build the container
echo "Building container..."
echo ""

# Try with fakeroot first, fall back to sudo if needed
if singularity build --fakeroot "$CONTAINER_NAME" "$DEF_FILE"; then
    BUILD_SUCCESS=true
# elif [ "$EUID" -eq 0 ] || sudo singularity build "$CONTAINER_NAME" "$DEF_FILE"; then
#     BUILD_SUCCESS=true
else
    BUILD_SUCCESS=false
fi

if [ "$BUILD_SUCCESS" = true ]; then
    echo ""
    echo "=========================================="
    echo "Build completed successfully!"
    echo "=========================================="
    echo ""
    echo "Container information:"
    ls -lh "$WORK_DIR/$CONTAINER_NAME"
    echo ""
    echo "Next steps:"
    echo "1. Test the container:"
    echo "   singularity exec $CONTAINER_NAME bash -c 'source activate long_reads && which minimap2'"
    echo ""
    echo "2. Run the pipeline:"
    echo "   nextflow run main.nf --sample <sample_name>"
    echo ""
else
    echo ""
    echo "=========================================="
    echo "Build failed!"
    echo "=========================================="
    echo ""
    echo "Possible solutions:"
    echo "1. Check that you have sufficient disk space (~10GB)"
    echo "2. Ensure Singularity is properly installed"
    echo "3. Check your internet connection"
    echo "4. Try running with sudo if fakeroot is not available"
    echo ""
    exit 1
fi
