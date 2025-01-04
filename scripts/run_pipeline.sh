#!/bin/bash
# Script to run the complete data processing pipeline

# Exit on error
set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Default values
OUTPUT_DIR="data"
CACHE_DIR="cache"
CHECKPOINT_DIR="checkpoints"
MODEL_DIR="models"
LOG_DIR="logs"
REPORT_DIR="reports"
LOG_LEVEL="INFO"
PARALLEL=true
RESUME=false
GPU=false
SKIP_WEB=false
SKIP_ANALYSIS=false
SKIP_REPORT=false
SKIP_WEB_APP=false
TARGETS="5-HT2A,5-HT2B,5-HT2C,NMDA,D2,SERT,NET,DAT"
ACTIVITY_TYPES="Ki,IC50,EC50,Kd"
MIN_CONFIDENCE=0.7
SOURCES="chembl,pubchem,swiss,community,social"
REPORT_TYPES="overview,activity,predictions,safety,community"
REPORT_FORMAT="html"
BATCH_SIZE=100
WEB_HOST="localhost"
WEB_PORT=8000
DEV_MODE=false

# Function to show usage
show_help() {
    echo "Run the complete data processing pipeline"
    echo
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  --output-dir DIR      Output directory (default: data)"
    echo "  --cache-dir DIR       Cache directory (default: cache)"
    echo "  --checkpoint-dir DIR  Checkpoint directory (default: checkpoints)"
    echo "  --model-dir DIR       Model directory (default: models)"
    echo "  --log-dir DIR         Log directory (default: logs)"
    echo "  --report-dir DIR      Report directory (default: reports)"
    echo "  --log-level LEVEL     Log level (default: INFO)"
    echo "  --no-parallel         Disable parallel processing"
    echo "  --resume              Resume from checkpoint"
    echo "  --gpu                 Enable GPU acceleration"
    echo "  --skip-web           Skip web data enrichment"
    echo "  --skip-analysis      Skip compound analysis"
    echo "  --skip-report        Skip report generation"
    echo "  --skip-web-app       Skip web application"
    echo "  --targets LIST        Target receptors (default: $TARGETS)"
    echo "  --activity-types LIST Activity types (default: $ACTIVITY_TYPES)"
    echo "  --min-confidence NUM  Minimum confidence (default: $MIN_CONFIDENCE)"
    echo "  --sources LIST        Data sources (default: $SOURCES)"
    echo "  --report-types LIST   Report types (default: $REPORT_TYPES)"
    echo "  --report-format FMT   Report format (default: $REPORT_FORMAT)"
    echo "  --batch-size NUM      Batch size (default: $BATCH_SIZE)"
    echo "  --web-host HOST       Web app host (default: $WEB_HOST)"
    echo "  --web-port PORT       Web app port (default: $WEB_PORT)"
    echo "  --dev                Enable development mode"
    echo "  --help               Show this help message"
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --output-dir)
            OUTPUT_DIR="$2"
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
        --model-dir)
            MODEL_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --report-dir)
            REPORT_DIR="$2"
            shift 2
            ;;
        --log-level)
            LOG_LEVEL="$2"
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
        --skip-web)
            SKIP_WEB=true
            shift
            ;;
        --skip-analysis)
            SKIP_ANALYSIS=true
            shift
            ;;
        --skip-report)
            SKIP_REPORT=true
            shift
            ;;
        --skip-web-app)
            SKIP_WEB_APP=true
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
        --sources)
            SOURCES="$2"
            shift 2
            ;;
        --report-types)
            REPORT_TYPES="$2"
            shift 2
            ;;
        --report-format)
            REPORT_FORMAT="$2"
            shift 2
            ;;
        --batch-size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --web-host)
            WEB_HOST="$2"
            shift 2
            ;;
        --web-port)
            WEB_PORT="$2"
            shift 2
            ;;
        --dev)
            DEV_MODE=true
            shift
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

# Function to check Python environment
check_env() {
    echo -e "${BLUE}Checking Python environment...${NC}"
    
    # Check if in virtual environment
    if [[ -z "$VIRTUAL_ENV" ]]; then
        echo -e "${RED}Not running in a virtual environment${NC}"
        echo "Please activate the virtual environment first:"
        echo "source .venv/bin/activate"
        exit 1
    fi
    
    # Check required packages
    python -c "import rdkit" &>/dev/null || {
        echo -e "${RED}RDKit not installed${NC}"
        exit 1
    }
    
    if [ "$GPU" = true ]; then
        python -c "import torch; assert torch.cuda.is_available()" &>/dev/null || {
            echo -e "${RED}GPU requested but PyTorch CUDA not available${NC}"
            exit 1
        }
    fi
}

# Create directories
mkdir -p "$OUTPUT_DIR" "$CACHE_DIR" "$CHECKPOINT_DIR" "$MODEL_DIR" "$LOG_DIR" "$REPORT_DIR"

# Set file paths
BINDINGDB_FILE="$OUTPUT_DIR/bindingdb.tsv"
ENRICHED_FILE="$OUTPUT_DIR/enriched.tsv"
ANALYZED_FILE="$OUTPUT_DIR/analyzed.tsv"

# Print configuration
echo -e "${BLUE}Starting pipeline with configuration:${NC}"
echo "  Output directory: $OUTPUT_DIR"
echo "  Cache directory: $CACHE_DIR"
echo "  Checkpoint directory: $CHECKPOINT_DIR"
echo "  Model directory: $MODEL_DIR"
echo "  Log directory: $LOG_DIR"
echo "  Report directory: $REPORT_DIR"
echo "  Log level: $LOG_LEVEL"
echo "  Parallel processing: $PARALLEL"
echo "  Resume from checkpoint: $RESUME"
echo "  GPU enabled: $GPU"
echo "  Skip web enrichment: $SKIP_WEB"
echo "  Skip analysis: $SKIP_ANALYSIS"
echo "  Skip report: $SKIP_REPORT"
echo "  Skip web app: $SKIP_WEB_APP"
echo "  Target receptors: $TARGETS"
echo "  Activity types: $ACTIVITY_TYPES"
echo "  Minimum confidence: $MIN_CONFIDENCE"
echo "  Data sources: $SOURCES"
echo "  Report types: $REPORT_TYPES"
echo "  Report format: $REPORT_FORMAT"
echo "  Batch size: $BATCH_SIZE"
echo "  Web host: $WEB_HOST"
echo "  Web port: $WEB_PORT"
echo "  Development mode: $DEV_MODE"
echo

# Check environment
check_env

# Step 1: Process BindingDB data
echo -e "${BLUE}Step 1: Processing BindingDB data...${NC}"
./scripts/process_bindingdb.sh \
    --output "$BINDINGDB_FILE" \
    --cache-dir "$CACHE_DIR" \
    --checkpoint-dir "$CHECKPOINT_DIR" \
    --log-level "$LOG_LEVEL" \
    --log-file "$LOG_DIR/bindingdb.log" \
    --targets "$TARGETS" \
    --activity-types "$ACTIVITY_TYPES" \
    --min-confidence "$MIN_CONFIDENCE" \
    --batch-size "$BATCH_SIZE" \
    $([ "$PARALLEL" = false ] && echo "--no-parallel") \
    $([ "$RESUME" = true ] && echo "--resume")

# Step 2: Enrich with web data (optional)
if [ "$SKIP_WEB" = false ]; then
    echo -e "${BLUE}Step 2: Enriching compounds with web data...${NC}"
    ./scripts/enrich_compounds.sh \
        --input "$BINDINGDB_FILE" \
        --output "$ENRICHED_FILE" \
        --cache-dir "$CACHE_DIR" \
        --checkpoint-dir "$CHECKPOINT_DIR" \
        --log-level "$LOG_LEVEL" \
        --log-file "$LOG_DIR/enrichment.log" \
        --sources "$SOURCES" \
        --min-confidence "$MIN_CONFIDENCE" \
        --batch-size "$BATCH_SIZE" \
        $([ "$PARALLEL" = false ] && echo "--no-parallel") \
        $([ "$RESUME" = true ] && echo "--resume")
else
    echo -e "${YELLOW}Skipping web data enrichment${NC}"
    ENRICHED_FILE="$BINDINGDB_FILE"
fi

# Step 3: Analyze compounds (optional)
if [ "$SKIP_ANALYSIS" = false ]; then
    echo -e "${BLUE}Step 3: Analyzing compounds...${NC}"
    ./scripts/analyze_compounds.sh \
        --input "$ENRICHED_FILE" \
        --output "$ANALYZED_FILE" \
        --cache-dir "$CACHE_DIR" \
        --checkpoint-dir "$CHECKPOINT_DIR" \
        --model-dir "$MODEL_DIR" \
        --log-level "$LOG_LEVEL" \
        --log-file "$LOG_DIR/analysis.log" \
        --min-confidence "$MIN_CONFIDENCE" \
        --batch-size "$BATCH_SIZE" \
        $([ "$PARALLEL" = false ] && echo "--no-parallel") \
        $([ "$RESUME" = true ] && echo "--resume") \
        $([ "$GPU" = true ] && echo "--gpu")
else
    echo -e "${YELLOW}Skipping compound analysis${NC}"
    ANALYZED_FILE="$ENRICHED_FILE"
fi

# Step 4: Generate reports (optional)
if [ "$SKIP_REPORT" = false ]; then
    echo -e "${BLUE}Step 4: Generating reports...${NC}"
    ./scripts/generate_report.sh \
        --input "$ANALYZED_FILE" \
        --output-dir "$REPORT_DIR" \
        --report-types "$REPORT_TYPES" \
        --format "$REPORT_FORMAT" \
        --log-level "$LOG_LEVEL" \
        --log-file "$LOG_DIR/report.log" \
        $([ "$PARALLEL" = false ] && echo "--no-parallel")
else
    echo -e "${YELLOW}Skipping report generation${NC}"
fi

# Step 5: Run web application (optional)
if [ "$SKIP_WEB_APP" = false ]; then
    echo -e "${BLUE}Step 5: Starting web application...${NC}"
    ./scripts/setup_and_run.sh \
        --input "$ANALYZED_FILE" \
        --host "$WEB_HOST" \
        --port "$WEB_PORT" \
        --log-level "$LOG_LEVEL" \
        --log-file "$LOG_DIR/webapp.log" \
        $([ "$DEV_MODE" = true ] && echo "--dev")
else
    echo -e "${YELLOW}Skipping web application${NC}"
fi

echo -e "${GREEN}Pipeline completed successfully!${NC}"
echo
echo "Generated files:"
echo "  BindingDB data: $BINDINGDB_FILE"
[ "$SKIP_WEB" = false ] && echo "  Enriched data: $ENRICHED_FILE"
[ "$SKIP_ANALYSIS" = false ] && echo "  Analyzed data: $ANALYZED_FILE"
[ "$SKIP_REPORT" = false ] && echo "  Reports: $REPORT_DIR/*.$REPORT_FORMAT"
[ "$SKIP_WEB_APP" = false ] && echo "  Web application: http://$WEB_HOST:$WEB_PORT"
