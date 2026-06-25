#!/bin/bash
# Quick build script for Singularity container
# This script helps build the Singularity container for the pipeline

set -e

CONTAINER_NAME="Singularity.sif"
DEF_FILE="Singularity.def"
DOCKER_IMAGE_NAME="sc_long_dev_container"
DOCKERFILE="Dockerfile"
WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=========================================="
echo "Long-Read Pipeline Singularity Build"
echo "=========================================="
echo ""
echo "Working directory: $WORK_DIR"
echo "Definition file: $DEF_FILE"
echo "Output container: $CONTAINER_NAME"
echo ""

# Check for Docker/Podman first; otherwise fallback to Singularity
USE_DOCKER=0
if command -v docker &> /dev/null; then
    USE_DOCKER=1
    CONTAINER_ENGINE=docker
elif command -v podman &> /dev/null; then
    USE_DOCKER=1
    CONTAINER_ENGINE=podman
fi

if [ "$USE_DOCKER" -eq 0 ]; then
    if ! command -v singularity &> /dev/null; then
        echo "ERROR: Neither Docker/Podman nor Singularity is installed."
        echo "Please install Docker/Podman or Singularity."
        exit 1
    fi
    echo "Singularity version:"
    singularity --version
    echo ""
else
    echo "Container engine detected: $CONTAINER_ENGINE"
    $CONTAINER_ENGINE --version
    echo ""

# Check for definition file or Dockerfile
if [ ! -f "$WORK_DIR/$DEF_FILE" ] && [ ! -f "$WORK_DIR/$DOCKERFILE" ]; then
    echo "ERROR: Neither $DEF_FILE nor $DOCKERFILE found in $WORK_DIR"
    exit 1
fi

# Check disk space (need ~10GB)
AVAILABLE_SPACE=$(df "$WORK_DIR" | awk 'NR==2 {print $4}')
REQUIRED_SPACE=$((10 * 1024 * 1024))  # 10GB in KB

if [ "$AVAILABLE_SPACE" -lt "$REQUIRED_SPACE" ]; then
    echo "WARNING: Available disk space is less than 10GB"
    echo "Available: $((AVAILABLE_SPACE / 1024 / 1024))GB"
    read -p "Continue anyway? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

echo "Building container..."
echo "This may take 15-30 minutes depending on your internet connection"
echo ""

if [ "$USE_DOCKER" -eq 1 ] && [ -f "$WORK_DIR/$DOCKERFILE" ]; then
    echo "Building with $CONTAINER_ENGINE using Dockerfile"
    $CONTAINER_ENGINE build -t $DOCKER_IMAGE_NAME -f "$WORK_DIR/$DOCKERFILE" "$WORK_DIR"
    if [ $? -ne 0 ]; then
        BUILD_SUCCESS=false
    else
        BUILD_SUCCESS=true
        if command -v singularity &> /dev/null; then
            echo "Converting Docker image to Singularity SIF"
            singularity build "$CONTAINER_NAME" "docker-daemon://$DOCKER_IMAGE_NAME:latest"
        fi
    fi
else
    # Try with fakeroot first, fall back to sudo if needed
    if singularity build --fakeroot "$CONTAINER_NAME" "$DEF_FILE"; then
        BUILD_SUCCESS=true
    elif [ "$EUID" -eq 0 ] || sudo singularity build "$CONTAINER_NAME" "$DEF_FILE"; then
        BUILD_SUCCESS=true
    else
        BUILD_SUCCESS=false
    fi
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
