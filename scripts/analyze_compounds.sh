#!/bin/bash
# Script to analyze and predict compound properties

# Exit on error
set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Default values
INPUT_FILE=""
OUTPUT_FILE=""
CACHE_DIR="cache"
CHECKPOINT_DIR="checkpoints"
LOG_DIR="logs"
LOG_LEVEL="INFO"
PARALLEL=true
RESUME=false
ANALYSIS_TYPES="activity,toxicity,abuse,bbb,nootropic"
PREDICTION_TYPES="binding,activity,toxicity,abuse"
MIN_CONFIDENCE=0.7
MODEL_DIR="models"
BATCH_SIZE=32
GPU=false

# Help message
show_help() {
    echo "Usage: $0 [options]"
    echo
    echo "Required arguments:"
    echo "  -i, --input FILE      Input TSV file with enriched compounds"
    echo "  -o, --output FILE     Output TSV file for analyzed compounds"
    echo
    echo "Optional arguments:"
    echo "  -m, --model-dir DIR   Model directory (default: $MODEL_DIR)"
    echo "  -c, --cache-dir DIR   Cache directory (default: $CACHE_DIR)"
    echo "  --checkpoint-dir DIR  Checkpoint directory (default: $CHECKPOINT_DIR)"
    echo "  --log-dir DIR         Log directory (default: $LOG_DIR)"
    echo "  -l, --log-level LVL   Log level: DEBUG, INFO, WARNING, ERROR (default: $LOG_LEVEL)"
    echo "  --analysis-types LIST Analysis types to run (default: $ANALYSIS_TYPES)"
    echo "  --prediction-types LIST Prediction types to run (default: $PREDICTION_TYPES)"
    echo "  --min-confidence N    Minimum confidence threshold (default: $MIN_CONFIDENCE)"
    echo "  -b, --batch-size N    Batch size (default: $BATCH_SIZE)"
    echo "  --no-parallel         Disable parallel processing"
    echo "  --resume              Resume from last checkpoint"
    echo "  --gpu                 Enable GPU acceleration"
    echo "  --skip-predictions    Skip ML predictions"
    echo "  --skip-analysis      Skip compound analysis"
    echo "  -h, --help           Show this help message"
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -m|--model-dir)
            MODEL_DIR="$2"
            shift 2
            ;;
        -c|--cache-dir)
            CACHE_DIR="$2"
            shift 2
            ;;
        --checkpoint-dir)
            CHECKPOINT_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        -l|--log-level)
            LOG_LEVEL="$2"
            shift 2
            ;;
        --analysis-types)
            ANALYSIS_TYPES="$2"
            shift 2
            ;;
        --prediction-types)
            PREDICTION_TYPES="$2"
            shift 2
            ;;
        --min-confidence)
            MIN_CONFIDENCE="$2"
            shift 2
            ;;
        -b|--batch-size)
            BATCH_SIZE="$2"
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
        --gpu)
            GPU=true
            shift
            ;;
        --skip-predictions)
            SKIP_PREDICTIONS="--skip-predictions"
            shift
            ;;
        --skip-analysis)
            SKIP_ANALYSIS="--skip-analysis"
            shift
            ;;
        -h|--help)
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

# Validate required arguments
if [ -z "$INPUT_FILE" ]; then
    echo -e "${RED}Error: --input argument is required${NC}"
    exit 1
fi

# Set default output file if not provided
if [ -z "$OUTPUT_FILE" ]; then
    OUTPUT_FILE="${INPUT_FILE%.*}_analyzed.tsv"
fi

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

# Create directories
mkdir -p "$(dirname "$OUTPUT_FILE")" "$MODEL_DIR" "$CACHE_DIR" "$CHECKPOINT_DIR" "$LOG_DIR"

# Activate virtual environment
source .venv/bin/activate

# Check GPU availability if requested
if [ "$GPU" = true ]; then
    if ! python -c "import torch; assert torch.cuda.is_available()"; then
        echo -e "${RED}Error: GPU requested but not available${NC}"
        exit 1
    fi
fi

# Set log file path
LOG_FILE="$LOG_DIR/analysis_$(date +%Y%m%d_%H%M%S).log"

# Build command
CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli analyze-compounds"
CMD="$CMD --input $INPUT_FILE"
CMD="$CMD --output $OUTPUT_FILE"
CMD="$CMD --model-dir $MODEL_DIR"
CMD="$CMD --cache-dir $CACHE_DIR"
CMD="$CMD --checkpoint-dir $CHECKPOINT_DIR"
CMD="$CMD --log-level $LOG_LEVEL"
CMD="$CMD --log-file $LOG_FILE"
CMD="$CMD --analysis-types $ANALYSIS_TYPES"
CMD="$CMD --prediction-types $PREDICTION_TYPES"
CMD="$CMD --min-confidence $MIN_CONFIDENCE"
CMD="$CMD --batch-size $BATCH_SIZE"

if [ "$PARALLEL" = false ]; then
    CMD="$CMD --no-parallel"
fi

if [ "$RESUME" = true ]; then
    CMD="$CMD --resume"
fi

if [ "$GPU" = true ]; then
    CMD="$CMD --gpu"
fi

if [ -n "$SKIP_PREDICTIONS" ]; then
    CMD="$CMD $SKIP_PREDICTIONS"
fi

if [ -n "$SKIP_ANALYSIS" ]; then
    CMD="$CMD $SKIP_ANALYSIS"
fi

# Print configuration
echo -e "${BLUE}Starting compound analysis with configuration:${NC}"
echo "  Input file: $INPUT_FILE"
echo "  Output file: $OUTPUT_FILE"
echo "  Model directory: $MODEL_DIR"
echo "  Cache directory: $CACHE_DIR"
echo "  Checkpoint directory: $CHECKPOINT_DIR"
echo "  Log directory: $LOG_DIR"
echo "  Log level: $LOG_LEVEL"
echo "  Log file: $LOG_FILE"
echo "  Analysis types: $ANALYSIS_TYPES"
echo "  Prediction types: $PREDICTION_TYPES"
echo "  Minimum confidence: $MIN_CONFIDENCE"
echo "  Batch size: $BATCH_SIZE"
echo "  Parallel processing: $PARALLEL"
echo "  Resume from checkpoint: $RESUME"
echo "  GPU enabled: $GPU"
echo

# Run analysis
echo -e "${BLUE}Running command:${NC}"
echo "$CMD"
echo
exec $CMD
