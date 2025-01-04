#!/bin/bash
# Script to set up and run the web application

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
OUTPUT_DIR="web/output"
HOST="localhost"
PORT=8000
DEV_MODE=false
RELOAD=false
LOG_LEVEL="INFO"
LOG_DIR="logs"
STATIC_DIR="web/static"
TEMPLATE_DIR="web/templates"

# Help message
show_help() {
    echo "Usage: $0 [options]"
    echo
    echo "Required arguments:"
    echo "  -i, --input FILE      Input TSV file with analyzed compounds"
    echo
    echo "Optional arguments:"
    echo "  -o, --output-dir DIR  Output directory (default: $OUTPUT_DIR)"
    echo "  --host HOST           Host to bind to (default: $HOST)"
    echo "  -p, --port PORT       Port to listen on (default: $PORT)"
    echo "  --static-dir DIR      Static files directory (default: $STATIC_DIR)"
    echo "  --template-dir DIR    Template directory (default: $TEMPLATE_DIR)"
    echo "  --log-dir DIR         Log directory (default: $LOG_DIR)"
    echo "  -l, --log-level LVL   Log level: DEBUG, INFO, WARNING, ERROR (default: $LOG_LEVEL)"
    echo "  --dev                Enable development mode (implies --reload)"
    echo "  --reload             Enable auto-reload on code changes"
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
        --host)
            HOST="$2"
            shift 2
            ;;
        -p|--port)
            PORT="$2"
            shift 2
            ;;
        --static-dir)
            STATIC_DIR="$2"
            shift 2
            ;;
        --template-dir)
            TEMPLATE_DIR="$2"
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
        --dev)
            DEV_MODE=true
            RELOAD=true
            shift
            ;;
        --reload)
            RELOAD=true
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
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$STATIC_DIR" "$TEMPLATE_DIR"

# Activate virtual environment
source .venv/bin/activate

# Set log file path
LOG_FILE="$LOG_DIR/webapp_$(date +%Y%m%d_%H%M%S).log"

# Build command
CMD="python -m binding_data_processor.web.app"
CMD="$CMD --input $INPUT_FILE"
CMD="$CMD --output-dir $OUTPUT_DIR"
CMD="$CMD --host $HOST"
CMD="$CMD --port $PORT"
CMD="$CMD --static-dir $STATIC_DIR"
CMD="$CMD --template-dir $TEMPLATE_DIR"
CMD="$CMD --log-level $LOG_LEVEL"
CMD="$CMD --log-file $LOG_FILE"

if [ "$DEV_MODE" = true ]; then
    CMD="$CMD --dev"
fi

if [ "$RELOAD" = true ]; then
    CMD="$CMD --reload"
fi

# Print configuration
echo -e "${BLUE}Starting web application with configuration:${NC}"
echo "  Input file: $INPUT_FILE"
echo "  Output directory: $OUTPUT_DIR"
echo "  Host: $HOST"
echo "  Port: $PORT"
echo "  Static directory: $STATIC_DIR"
echo "  Template directory: $TEMPLATE_DIR"
echo "  Log directory: $LOG_DIR"
echo "  Log level: $LOG_LEVEL"
echo "  Log file: $LOG_FILE"
echo "  Development mode: $DEV_MODE"
echo "  Auto-reload: $RELOAD"
echo

# Run web application
echo -e "${BLUE}Running command:${NC}"
echo "$CMD"
echo
exec $CMD
