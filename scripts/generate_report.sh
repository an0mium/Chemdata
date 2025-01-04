#!/bin/bash
# Script to generate analysis reports

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
OUTPUT_DIR="reports"
REPORT_TYPES="overview,activity,predictions,safety,community"
FORMAT="html"
TEMPLATE_DIR="templates"
STATIC_DIR="static"
CACHE_DIR="cache"
LOG_DIR="logs"
LOG_LEVEL="INFO"
PARALLEL=true
INCLUDE_PLOTS=true
INCLUDE_STRUCTURES=true
INCLUDE_PREDICTIONS=true
INCLUDE_COMMUNITY=true
PLOT_FORMAT="svg"
PLOT_DPI=300
PLOT_WIDTH=800
PLOT_HEIGHT=600

# Help message
show_help() {
    echo "Usage: $0 [options]"
    echo
    echo "Required arguments:"
    echo "  -i, --input FILE      Input TSV file with analyzed compounds"
    echo
    echo "Optional arguments:"
    echo "  -o, --output-dir DIR  Output directory (default: $OUTPUT_DIR)"
    echo "  -c, --cache-dir DIR   Cache directory (default: $CACHE_DIR)"
    echo "  --log-dir DIR         Log directory (default: $LOG_DIR)"
    echo "  -l, --log-level LVL   Log level: DEBUG, INFO, WARNING, ERROR (default: $LOG_LEVEL)"
    echo "  --report-types LIST   Report types to generate (default: $REPORT_TYPES)"
    echo "  --format FMT          Report format: html, pdf, md (default: $FORMAT)"
    echo "  --template-dir DIR    Template directory (default: $TEMPLATE_DIR)"
    echo "  --static-dir DIR      Static files directory (default: $STATIC_DIR)"
    echo "  --no-plots           Disable plot generation"
    echo "  --no-structures      Disable structure rendering"
    echo "  --no-predictions     Exclude ML predictions"
    echo "  --no-community       Exclude community data"
    echo "  --plot-format FMT    Plot format: svg, png, pdf (default: $PLOT_FORMAT)"
    echo "  --plot-dpi NUM       Plot DPI for raster formats (default: $PLOT_DPI)"
    echo "  --plot-width NUM     Plot width in pixels (default: $PLOT_WIDTH)"
    echo "  --plot-height NUM    Plot height in pixels (default: $PLOT_HEIGHT)"
    echo "  --no-parallel        Disable parallel processing"
    echo "  -h, --help           Show this help message"
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -c|--cache-dir)
            CACHE_DIR="$2"
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
        --report-types)
            REPORT_TYPES="$2"
            shift 2
            ;;
        --format)
            FORMAT="$2"
            shift 2
            ;;
        --template-dir)
            TEMPLATE_DIR="$2"
            shift 2
            ;;
        --static-dir)
            STATIC_DIR="$2"
            shift 2
            ;;
        --no-plots)
            INCLUDE_PLOTS=false
            shift
            ;;
        --no-structures)
            INCLUDE_STRUCTURES=false
            shift
            ;;
        --no-predictions)
            INCLUDE_PREDICTIONS=false
            shift
            ;;
        --no-community)
            INCLUDE_COMMUNITY=false
            shift
            ;;
        --plot-format)
            PLOT_FORMAT="$2"
            shift 2
            ;;
        --plot-dpi)
            PLOT_DPI="$2"
            shift 2
            ;;
        --plot-width)
            PLOT_WIDTH="$2"
            shift 2
            ;;
        --plot-height)
            PLOT_HEIGHT="$2"
            shift 2
            ;;
        --no-parallel)
            PARALLEL=false
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
mkdir -p "$OUTPUT_DIR" "$CACHE_DIR" "$LOG_DIR"

# Activate virtual environment
source .venv/bin/activate

# Set log file path
LOG_FILE="$LOG_DIR/report_$(date +%Y%m%d_%H%M%S).log"

# Build command
CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli generate-report"
CMD="$CMD --input $INPUT_FILE"
CMD="$CMD --output-dir $OUTPUT_DIR"
CMD="$CMD --cache-dir $CACHE_DIR"
CMD="$CMD --log-level $LOG_LEVEL"
CMD="$CMD --log-file $LOG_FILE"
CMD="$CMD --report-types $REPORT_TYPES"
CMD="$CMD --format $FORMAT"
CMD="$CMD --template-dir $TEMPLATE_DIR"
CMD="$CMD --static-dir $STATIC_DIR"
CMD="$CMD --plot-format $PLOT_FORMAT"
CMD="$CMD --plot-dpi $PLOT_DPI"
CMD="$CMD --plot-width $PLOT_WIDTH"
CMD="$CMD --plot-height $PLOT_HEIGHT"

if [ "$PARALLEL" = false ]; then
    CMD="$CMD --no-parallel"
fi

if [ "$INCLUDE_PLOTS" = false ]; then
    CMD="$CMD --no-plots"
fi

if [ "$INCLUDE_STRUCTURES" = false ]; then
    CMD="$CMD --no-structures"
fi

if [ "$INCLUDE_PREDICTIONS" = false ]; then
    CMD="$CMD --no-predictions"
fi

if [ "$INCLUDE_COMMUNITY" = false ]; then
    CMD="$CMD --no-community"
fi

# Print configuration
echo -e "${BLUE}Starting report generation with configuration:${NC}"
echo "  Input file: $INPUT_FILE"
echo "  Output directory: $OUTPUT_DIR"
echo "  Cache directory: $CACHE_DIR"
echo "  Log directory: $LOG_DIR"
echo "  Log level: $LOG_LEVEL"
echo "  Log file: $LOG_FILE"
echo "  Report types: $REPORT_TYPES"
echo "  Format: $FORMAT"
echo "  Template directory: $TEMPLATE_DIR"
echo "  Static directory: $STATIC_DIR"
echo "  Include plots: $INCLUDE_PLOTS"
echo "  Include structures: $INCLUDE_STRUCTURES"
echo "  Include predictions: $INCLUDE_PREDICTIONS"
echo "  Include community data: $INCLUDE_COMMUNITY"
echo "  Plot format: $PLOT_FORMAT"
echo "  Plot DPI: $PLOT_DPI"
echo "  Plot dimensions: ${PLOT_WIDTH}x${PLOT_HEIGHT}"
echo "  Parallel processing: $PARALLEL"
echo

# Run report generation
echo -e "${BLUE}Running command:${NC}"
echo "$CMD"
echo
exec $CMD
