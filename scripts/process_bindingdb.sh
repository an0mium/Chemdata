#!/bin/bash
# Script to download and process BindingDB data

# Exit on error
set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Default values
OUTPUT_FILE="data/compounds.tsv"
CACHE_DIR="cache"
CHECKPOINT_DIR="checkpoints"
LOG_LEVEL="INFO"
LOG_FILE="logs/bindingdb.log"
PARALLEL=true
RESUME=false
TARGETS="5-HT2A,5-HT2B,5-HT2C,NMDA,D2,SERT,NET,DAT"
ACTIVITY_TYPES="Ki,IC50,EC50,Kd"
MIN_CONFIDENCE=0.7

# Function to show usage
show_help() {
    echo "Download and process BindingDB data"
    echo
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  --output FILE        Output file (default: $OUTPUT_FILE)"
    echo "  --cache-dir DIR      Cache directory (default: $CACHE_DIR)"
    echo "  --checkpoint-dir DIR Checkpoint directory (default: $CHECKPOINT_DIR)"
    echo "  --log-level LEVEL    Log level (default: $LOG_LEVEL)"
    echo "  --log-file FILE      Log file (default: $LOG_FILE)"
    echo "  --targets LIST       Target receptors (default: $TARGETS)"
    echo "  --activity-types LIST Activity types (default: $ACTIVITY_TYPES)"
    echo "  --min-confidence NUM Minimum confidence (default: $MIN_CONFIDENCE)"
    echo "  --no-parallel        Disable parallel processing"
    echo "  --resume             Resume from checkpoint"
    echo "  --help              Show this help message"
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        --cache-dir)
            CACHE_DIR="$2"
            shift 2
            ;;
        --checkpoint-dir)
            CHECKPOINT_DIR="$2"
            shift 2
            ;;
        --log-level)
            LOG_LEVEL="$2"
            shift 2
            ;;
        --log-file)
            LOG_FILE="$2"
            shift 2
            ;;
        --no-parallel)
            PARALLEL=false
            shift
            ;;
        --resume)
            RESUME=true
            shift
            ;;
        --targets)
            TARGETS="$2"
            shift 2
            ;;
        --activity-types)
            ACTIVITY_TYPES="$2"
            shift 2
            ;;
        --min-confidence)
            MIN_CONFIDENCE="$2"
            shift 2
            ;;
        --help|-h)
            show_help
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            show_help
            exit 1
            ;;
    esac
done

# Create directories
mkdir -p "$(dirname "$OUTPUT_FILE")" "$CACHE_DIR" "$CHECKPOINT_DIR" "$(dirname "$LOG_FILE")"

# Check Python version
PYTHON_VERSION=$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
REQUIRED_VERSION="3.8"

if [ "$(printf '%s\n' "$REQUIRED_VERSION" "$PYTHON_VERSION" | sort -V | head -n1)" != "$REQUIRED_VERSION" ]; then
    echo -e "${RED}Error: Python $REQUIRED_VERSION or higher is required${NC}"
    exit 1
fi

# Check if virtual environment exists
if [ ! -d ".venv" ]; then
    echo -e "${RED}Error: Virtual environment not found${NC}"
    echo "Please run setup_dev.sh first"
    exit 1
fi

# Activate virtual environment
source .venv/bin/activate

# Print configuration
echo -e "${BLUE}Starting BindingDB processing with configuration:${NC}"
echo "  Output file: $OUTPUT_FILE"
echo "  Cache directory: $CACHE_DIR"
echo "  Checkpoint directory: $CHECKPOINT_DIR"
echo "  Log level: $LOG_LEVEL"
echo "  Log file: $LOG_FILE"
echo "  Parallel processing: $PARALLEL"
echo "  Resume from checkpoint: $RESUME"
echo "  Target receptors: $TARGETS"
echo "  Activity types: $ACTIVITY_TYPES"
echo "  Minimum confidence: $MIN_CONFIDENCE"
echo

# Build command
CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli process-bindingdb"
CMD="$CMD --output $OUTPUT_FILE"
CMD="$CMD --cache-dir $CACHE_DIR"
CMD="$CMD --checkpoint-dir $CHECKPOINT_DIR"
CMD="$CMD --log-level $LOG_LEVEL"
CMD="$CMD --log-file $LOG_FILE"
CMD="$CMD --targets $TARGETS"
CMD="$CMD --activity-types $ACTIVITY_TYPES"
CMD="$CMD --min-confidence $MIN_CONFIDENCE"

if [ "$PARALLEL" = false ]; then
    CMD="$CMD --no-parallel"
fi

if [ "$RESUME" = true ]; then
    CMD="$CMD --resume"
fi

# Run processing
echo -e "${BLUE}Running command:${NC}"
echo "$CMD"
echo
exec $CMD
