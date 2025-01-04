#!/bin/bash
# Script to enrich compound data with web sources

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
MODEL_DIR="models"
CACHE_DIR="cache"
CHECKPOINT_DIR="checkpoints"
LOG_DIR="logs"
LOG_LEVEL="INFO"
N_WORKERS=4
BATCH_SIZE=100
RATE_LIMIT=2
TIMEOUT=30
MAX_RETRIES=3
MIN_CONFIDENCE=0.7
PARALLEL=true
RESUME=false
SOURCES="chembl,pubchem,swiss,community,social"

# Help message
show_help() {
    echo "Usage: $0 [options]"
    echo
    echo "Required arguments:"
    echo "  -i, --input FILE      Input TSV file with compounds"
    echo "  -o, --output FILE     Output TSV file for enriched compounds"
    echo
    echo "Optional arguments:"
    echo "  -m, --model-dir DIR   Model directory (default: $MODEL_DIR)"
    echo "  -c, --cache-dir DIR   Cache directory (default: $CACHE_DIR)"
    echo "  --checkpoint-dir DIR  Checkpoint directory (default: $CHECKPOINT_DIR)"
    echo "  --log-dir DIR         Log directory (default: $LOG_DIR)"
    echo "  -l, --log-level LVL   Log level: DEBUG, INFO, WARNING, ERROR (default: $LOG_LEVEL)"
    echo "  -w, --workers N       Number of worker threads (default: $N_WORKERS)"
    echo "  -b, --batch-size N    Batch size (default: $BATCH_SIZE)"
    echo "  -r, --rate-limit N    Rate limit in requests/second (default: $RATE_LIMIT)"
    echo "  -t, --timeout N       Request timeout in seconds (default: $TIMEOUT)"
    echo "  --max-retries N       Maximum retries per request (default: $MAX_RETRIES)"
    echo "  --min-confidence N    Minimum confidence threshold (default: $MIN_CONFIDENCE)"
    echo "  --sources LIST        Comma-separated list of sources (default: $SOURCES)"
    echo "  --no-parallel         Disable parallel processing"
    echo "  --resume              Resume from last checkpoint"
    echo "  --skip-predictions    Skip ML predictions"
    echo "  --skip-web-data      Skip web data enrichment"
    echo "  --no-cache           Disable caching"
    echo "  -h, --help           Show this help message"
    echo
    echo "Environment variables:"
    echo "  REDDIT_CLIENT_ID      Reddit API client ID"
    echo "  REDDIT_CLIENT_SECRET  Reddit API client secret"
    echo "  TWITTER_BEARER_TOKEN  Twitter API bearer token"
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
        -w|--workers)
            N_WORKERS="$2"
            shift 2
            ;;
        -b|--batch-size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        -r|--rate-limit)
            RATE_LIMIT="$2"
            shift 2
            ;;
        -t|--timeout)
            TIMEOUT="$2"
            shift 2
            ;;
        --max-retries)
            MAX_RETRIES="$2"
            shift 2
            ;;
        --min-confidence)
            MIN_CONFIDENCE="$2"
            shift 2
            ;;
        --sources)
            SOURCES="$2"
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
        --skip-predictions)
            SKIP_PREDICTIONS="--skip-predictions"
            shift
            ;;
        --skip-web-data)
            SKIP_WEB_DATA="--skip-web-data"
            shift
            ;;
        --no-cache)
            NO_CACHE="--no-cache"
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
    OUTPUT_FILE="${INPUT_FILE%.*}_enriched.tsv"
fi

# Check for required environment variables
if [[ -z "$REDDIT_CLIENT_ID" ]]; then
    echo -e "${RED}Error: REDDIT_CLIENT_ID environment variable not set${NC}"
    exit 1
fi

if [[ -z "$REDDIT_CLIENT_SECRET" ]]; then
    echo -e "${RED}Error: REDDIT_CLIENT_SECRET environment variable not set${NC}"
    exit 1
fi

if [[ -z "$TWITTER_BEARER_TOKEN" ]]; then
    echo -e "${RED}Error: TWITTER_BEARER_TOKEN environment variable not set${NC}"
    exit 1
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

# Set log file path
LOG_FILE="$LOG_DIR/enrichment_$(date +%Y%m%d_%H%M%S).log"

# Build command
CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli enrich-compounds"
CMD="$CMD --input $INPUT_FILE"
CMD="$CMD --output $OUTPUT_FILE"
CMD="$CMD --model-dir $MODEL_DIR"
CMD="$CMD --cache-dir $CACHE_DIR"
CMD="$CMD --checkpoint-dir $CHECKPOINT_DIR"
CMD="$CMD --log-level $LOG_LEVEL"
CMD="$CMD --log-file $LOG_FILE"
CMD="$CMD --n-workers $N_WORKERS"
CMD="$CMD --batch-size $BATCH_SIZE"
CMD="$CMD --rate-limit $RATE_LIMIT"
CMD="$CMD --timeout $TIMEOUT"
CMD="$CMD --max-retries $MAX_RETRIES"
CMD="$CMD --min-confidence $MIN_CONFIDENCE"
CMD="$CMD --sources $SOURCES"
CMD="$CMD --reddit-id $REDDIT_CLIENT_ID"
CMD="$CMD --reddit-secret $REDDIT_CLIENT_SECRET"
CMD="$CMD --twitter-token $TWITTER_BEARER_TOKEN"

if [ "$PARALLEL" = false ]; then
    CMD="$CMD --no-parallel"
fi

if [ "$RESUME" = true ]; then
    CMD="$CMD --resume"
fi

if [ -n "$SKIP_PREDICTIONS" ]; then
    CMD="$CMD $SKIP_PREDICTIONS"
fi

if [ -n "$SKIP_WEB_DATA" ]; then
    CMD="$CMD $SKIP_WEB_DATA"
fi

if [ -n "$NO_CACHE" ]; then
    CMD="$CMD $NO_CACHE"
fi

# Print configuration
echo -e "${BLUE}Starting compound enrichment with configuration:${NC}"
echo "  Input file: $INPUT_FILE"
echo "  Output file: $OUTPUT_FILE"
echo "  Model directory: $MODEL_DIR"
echo "  Cache directory: $CACHE_DIR"
echo "  Checkpoint directory: $CHECKPOINT_DIR"
echo "  Log directory: $LOG_DIR"
echo "  Log level: $LOG_LEVEL"
echo "  Log file: $LOG_FILE"
echo "  Workers: $N_WORKERS"
echo "  Batch size: $BATCH_SIZE"
echo "  Rate limit: $RATE_LIMIT requests/second"
echo "  Timeout: $TIMEOUT seconds"
echo "  Maximum retries: $MAX_RETRIES"
echo "  Minimum confidence: $MIN_CONFIDENCE"
echo "  Data sources: $SOURCES"
echo "  Parallel processing: $PARALLEL"
echo "  Resume from checkpoint: $RESUME"
echo

# Run enrichment
echo -e "${BLUE}Running command:${NC}"
echo "$CMD"
echo
exec $CMD
