#!/bin/bash
set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Function to show usage
show_help() {
    echo "Run performance benchmarks"
    echo
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  --all       Run all benchmarks"
    echo "  --data      Run data processing benchmarks only"
    echo "  --ml        Run ML model benchmarks only"
    echo "  --web       Run web server benchmarks only"
    echo "  --memory    Run memory usage benchmarks only"
    echo "  --gpu       Run GPU utilization benchmarks only"
    echo "  --deps      Run dependency manager benchmarks only"
    echo "  --quick     Run quick benchmarks (smaller dataset)"
    echo "  --parallel  Run benchmarks in parallel where possible"
    echo "  --help      Show this help message"
    echo
    echo "Examples:"
    echo "  $0 --all          Run all benchmarks"
    echo "  $0 --ml --quick   Run quick ML benchmarks"
    echo "  $0 --deps         Compare pip vs uv performance"
}

# Function to measure time with microsecond precision
measure_time() {
    start_time=$(date +%s.%N)
    "$@"
    end_time=$(date +%s.%N)
    echo "$(echo "$end_time - $start_time" | bc)"
}

# Function to format time
format_time() {
    if (( $(echo "$1 < 60" | bc -l) )); then
        printf "%.2f seconds" "$1"
    elif (( $(echo "$1 < 3600" | bc -l) )); then
        mins=$(echo "$1 / 60" | bc -l)
        printf "%.2f minutes" "$mins"
    else
        hours=$(echo "$1 / 3600" | bc -l)
        printf "%.2f hours" "$hours"
    fi
}

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check environment
check_env() {
    echo -e "${BLUE}Checking environment...${NC}"
    
    # Check Python
    if ! command_exists python3; then
        echo -e "${RED}Python 3 not found${NC}"
        exit 1
    fi
    
    # Check virtual environment
    if [[ -z "${VIRTUAL_ENV}" ]]; then
        echo -e "${RED}Not in a virtual environment${NC}"
        echo "Activate virtual environment first:"
        echo "source .venv/bin/activate"
        exit 1
    fi
    
    # Check required tools
    for tool in curl wget ab time bc; do
        if ! command_exists "$tool"; then
            echo -e "${RED}Required tool not found: $tool${NC}"
            exit 1
        fi
    done
    
    # Check resource monitoring tools
    if ! command_exists top || ! command_exists ps; then
        echo -e "${YELLOW}Resource monitoring tools not found${NC}"
        echo "Some benchmarks may have limited information"
    fi
    
    echo -e "${GREEN}Environment check passed${NC}"
}

# Function to clean environment for dependency benchmarks
clean_env() {
    echo -e "${BLUE}Cleaning environment...${NC}"
    rm -rf .venv_uv/ .venv_pip/
    rm -rf ~/.cache/pip/* ~/.cache/uv/*
}

# Function to run pip benchmark
benchmark_pip() {
    echo -e "${BLUE}Running pip benchmark...${NC}"
    
    # Create virtual environment
    python3 -m venv .venv_pip
    source .venv_pip/bin/activate
    
    # Measure cold cache installation time
    pip_cold_time=$(measure_time python -m pip install -e ".[dev]")
    
    # Measure warm cache installation time
    pip_warm_time=$(measure_time python -m pip install -e ".[dev]")
    
    # Measure uninstall time
    pip_uninstall_time=$(measure_time python -m pip uninstall -y -r requirements.txt)
    
    deactivate
    
    echo "$pip_cold_time:$pip_warm_time:$pip_uninstall_time"
}

# Function to run uv benchmark
benchmark_uv() {
    echo -e "${BLUE}Running uv benchmark...${NC}"
    
    # Create virtual environment
    uv venv -p python3.12 .venv_uv
    source .venv_uv/bin/activate
    
    # Measure cold cache installation time
    uv_cold_time=$(measure_time uv pip install -e ".[dev]")
    
    # Measure warm cache installation time
    uv_warm_time=$(measure_time uv pip install -e ".[dev]")
    
    # Measure uninstall time
    uv_uninstall_time=$(measure_time uv pip uninstall -y -r requirements.txt)
    
    deactivate
    
    echo "$uv_cold_time:$uv_warm_time:$uv_uninstall_time"
}

# Function to run dependency benchmarks
benchmark_deps() {
    echo -e "${BLUE}Running dependency manager benchmarks...${NC}"
    
    # Create results directory
    mkdir -p reports/benchmark/deps
    
    # Clean environment
    clean_env
    
    # Run benchmarks
    pip_times=$(benchmark_pip)
    pip_cold_time=$(echo "$pip_times" | cut -d: -f1)
    pip_warm_time=$(echo "$pip_times" | cut -d: -f2)
    pip_uninstall_time=$(echo "$pip_times" | cut -d: -f3)
    
    uv_times=$(benchmark_uv)
    uv_cold_time=$(echo "$uv_times" | cut -d: -f1)
    uv_warm_time=$(echo "$uv_times" | cut -d: -f2)
    uv_uninstall_time=$(echo "$uv_times" | cut -d: -f3)
    
    # Calculate speedups
    cold_speedup=$(echo "scale=2; $pip_cold_time / $uv_cold_time" | bc)
    warm_speedup=$(echo "scale=2; $pip_warm_time / $uv_warm_time" | bc)
    uninstall_speedup=$(echo "scale=2; $pip_uninstall_time / $uv_uninstall_time" | bc)
    
    # Save results
    cat > reports/benchmark/deps/results.txt << EOL
pip cold install: $(format_time "$pip_cold_time")
pip warm install: $(format_time "$pip_warm_time")
pip uninstall: $(format_time "$pip_uninstall_time")

uv cold install: $(format_time "$uv_cold_time")
uv warm install: $(format_time "$uv_warm_time")
uv uninstall: $(format_time "$uv_uninstall_time")

Speedups:
Cold install: ${cold_speedup}x
Warm install: ${warm_speedup}x
Uninstall: ${uninstall_speedup}x
EOL
    
    echo -e "${GREEN}Dependency benchmarks complete${NC}"
}

# Function to create benchmark data
create_benchmark_data() {
    echo -e "${BLUE}Creating benchmark data...${NC}"
    
    mkdir -p data/benchmark
    
    if [[ "$QUICK" == "true" ]]; then
        # Create small dataset
        head -n 1000 data/raw/BindingDB_All.tsv > data/benchmark/test_data.tsv
    else
        # Create full dataset
        cp data/raw/BindingDB_All.tsv data/benchmark/test_data.tsv
    fi
    
    echo -e "${GREEN}Benchmark data created${NC}"
}

# Function to run data processing benchmarks
benchmark_data() {
    echo -e "${BLUE}Running data processing benchmarks...${NC}"
    
    # Create results directory
    mkdir -p reports/benchmark/data
    
    # Time data loading
    echo -e "${BLUE}Testing data loading...${NC}"
    /usr/bin/time -v python3 -c "
from binding_data_processor.pipeline import PipelineManager
pipeline = PipelineManager()
pipeline.load_data('data/benchmark/test_data.tsv')
" 2> reports/benchmark/data/load_time.txt
    
    # Time data processing
    echo -e "${BLUE}Testing data processing...${NC}"
    /usr/bin/time -v python3 -c "
from binding_data_processor.pipeline import PipelineManager
pipeline = PipelineManager()
pipeline.process_data('data/benchmark/test_data.tsv', 'data/benchmark/output')
" 2> reports/benchmark/data/process_time.txt
    
    echo -e "${GREEN}Data processing benchmarks complete${NC}"
}

# Function to run ML benchmarks
benchmark_ml() {
    echo -e "${BLUE}Running ML benchmarks...${NC}"
    
    # Create results directory
    mkdir -p reports/benchmark/ml
    
    # Test model loading time
    echo -e "${BLUE}Testing model loading...${NC}"
    /usr/bin/time -v python3 -c "
from binding_data_processor.processors.structure.ml.predictors.activity import ActivityPredictor
predictor = ActivityPredictor()
predictor.load_models()
" 2> reports/benchmark/ml/load_time.txt
    
    # Test inference time
    echo -e "${BLUE}Testing inference...${NC}"
    /usr/bin/time -v python3 -c "
from binding_data_processor.processors.structure.ml.predictors.activity import ActivityPredictor
from rdkit import Chem
predictor = ActivityPredictor()
mol = Chem.MolFromSmiles('CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC')
for _ in range(100):
    predictor.predict_compound(mol)
" 2> reports/benchmark/ml/inference_time.txt
    
    # Test batch inference time
    echo -e "${BLUE}Testing batch inference...${NC}"
    /usr/bin/time -v python3 -c "
from binding_data_processor.processors.structure.ml.predictors.activity import ActivityPredictor
from rdkit import Chem
import torch
predictor = ActivityPredictor()
mol = Chem.MolFromSmiles('CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC')
batch = [mol] * 100
with torch.no_grad():
    predictor.predict_batch(batch)
" 2> reports/benchmark/ml/batch_time.txt
    
    # Test GPU utilization if available
    if command_exists nvidia-smi; then
        echo -e "${BLUE}Testing GPU utilization...${NC}"
        nvidia-smi dmon -s u -c 10 > reports/benchmark/ml/gpu_util.txt &
        python3 -c "
from binding_data_processor.processors.structure.ml.predictors.activity import ActivityPredictor
from rdkit import Chem
predictor = ActivityPredictor()
mol = Chem.MolFromSmiles('CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC')
for _ in range(1000):
    predictor.predict_compound(mol)
"
        pkill -f "nvidia-smi dmon"
    fi
    
    echo -e "${GREEN}ML benchmarks complete${NC}"
}

# Function to run web benchmarks
benchmark_web() {
    echo -e "${BLUE}Running web benchmarks...${NC}"
    
    # Create results directory
    mkdir -p reports/benchmark/web
    
    # Start web server
    python3 -m web.app &
    WEB_PID=$!
    sleep 5  # Wait for server to start
    
    # Run Apache Bench tests
    echo -e "${BLUE}Testing web server performance...${NC}"
    ab -n 1000 -c 10 http://localhost:8050/ > reports/benchmark/web/ab_results.txt
    
    # Test API endpoints
    echo -e "${BLUE}Testing API endpoints...${NC}"
    ab -n 100 -c 5 http://localhost:8050/api/compounds > reports/benchmark/web/api_results.txt
    
    # Test WebSocket performance if available
    if [[ -f web/websocket.py ]]; then
        echo -e "${BLUE}Testing WebSocket performance...${NC}"
        python3 tests/websocket_benchmark.py > reports/benchmark/web/ws_results.txt
    fi
    
    # Kill web server
    kill $WEB_PID
    
    echo -e "${GREEN}Web benchmarks complete${NC}"
}

# Function to run memory benchmarks
benchmark_memory() {
    echo -e "${BLUE}Running memory benchmarks...${NC}"
    
    # Create results directory
    mkdir -p reports/benchmark/memory
    
    # Test memory usage during data processing
    echo -e "${BLUE}Testing data processing memory...${NC}"
    /usr/bin/time -v python3 -c "
from binding_data_processor.pipeline import PipelineManager
pipeline = PipelineManager()
pipeline.process_data('data/benchmark/test_data.tsv', 'data/benchmark/output')
" 2> reports/benchmark/memory/data_memory.txt
    
    # Test memory usage during ML inference
    echo -e "${BLUE}Testing ML inference memory...${NC}"
    /usr/bin/time -v python3 -c "
from binding_data_processor.processors.structure.ml.predictors.activity import ActivityPredictor
from rdkit import Chem
predictor = ActivityPredictor()
mol = Chem.MolFromSmiles('CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC')
for _ in range(100):
    predictor.predict_compound(mol)
" 2> reports/benchmark/memory/ml_memory.txt
    
    # Test memory leaks
    echo -e "${BLUE}Testing for memory leaks...${NC}"
    python3 -m memory_profiler tests/memory_test.py > reports/benchmark/memory/leaks.txt
    
    echo -e "${GREEN}Memory benchmarks complete${NC}"
}

# Function to run GPU benchmarks
benchmark_gpu() {
    if ! command_exists nvidia-smi; then
        echo -e "${YELLOW}No GPU found, skipping GPU benchmarks${NC}"
        return 0
    fi
    
    echo -e "${BLUE}Running GPU benchmarks...${NC}"
    
    # Create results directory
    mkdir -p reports/benchmark/gpu
    
    # Test GPU memory usage
    echo -e "${BLUE}Testing GPU memory usage...${NC}"
    nvidia-smi dmon -s m -c 10 > reports/benchmark/gpu/memory.txt &
    python3 -c "
from binding_data_processor.processors.structure.ml.predictors.activity import ActivityPredictor
from rdkit import Chem
predictor = ActivityPredictor()
mol = Chem.MolFromSmiles('CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC')
for _ in range(1000):
    predictor.predict_compound(mol)
"
    pkill -f "nvidia-smi dmon"
    
    # Test GPU utilization
    echo -e "${BLUE}Testing GPU utilization...${NC}"
    nvidia-smi dmon -s u -c 10 > reports/benchmark/gpu/utilization.txt &
    python3 -c "
from binding_data_processor.processors.structure.ml.predictors.activity import ActivityPredictor
from rdkit import Chem
predictor = ActivityPredictor()
mol = Chem.MolFromSmiles('CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC')
for _ in range(1000):
    predictor.predict_compound(mol)
"
    pkill -f "nvidia-smi dmon"
    
    # Test multi-GPU if available
    gpu_count=$(nvidia-smi -L | wc -l)
    if [[ $gpu_count -gt 1 ]]; then
        echo -e "${BLUE}Testing multi-GPU performance...${NC}"
        python3 tests/multi_gpu_test.py > reports/benchmark/gpu/multi_gpu.txt
    fi
    
    echo -e "${GREEN}GPU benchmarks complete${NC}"
}

# Function to generate report
generate_report() {
    echo -e "${BLUE}Generating benchmark report...${NC}"
    
    # Create report directory
    mkdir -p reports/benchmark
    
    # Generate HTML report
    python3 scripts/generate_benchmark_report.py \
        --deps reports/benchmark/deps/results.txt \
        --data reports/benchmark/data \
        --ml reports/benchmark/ml \
        --web reports/benchmark/web \
        --memory reports/benchmark/memory \
        --gpu reports/benchmark/gpu \
        --output reports/benchmark/report.html
    
    echo -e "${GREEN}Report generated: reports/benchmark/report.html${NC}"
}

# Parse command line arguments
RUN_DATA=0
RUN_ML=0
RUN_WEB=0
RUN_MEMORY=0
RUN_GPU=0
RUN_DEPS=0
QUICK=false
PARALLEL=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --all)
            RUN_DATA=1
            RUN_ML=1
            RUN_WEB=1
            RUN_MEMORY=1
            RUN_GPU=1
            RUN_DEPS=1
            shift
            ;;
        --data)
            RUN_DATA=1
            shift
            ;;
        --ml)
            RUN_ML=1
            shift
            ;;
        --web)
            RUN_WEB=1
            shift
            ;;
        --memory)
            RUN_MEMORY=1
            shift
            ;;
        --gpu)
            RUN_GPU=1
            shift
            ;;
        --deps)
            RUN_DEPS=1
            shift
            ;;
        --quick)
            QUICK=true
            shift
            ;;
        --parallel)
            PARALLEL=true
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

# If no specific options provided, run all benchmarks
if [[ $RUN_DATA -eq 0 && $RUN_ML -eq 0 && $RUN_WEB -eq 0 && \
      $RUN_MEMORY -eq 0 && $RUN_GPU -eq 0 && $RUN_DEPS -eq 0 ]]; then
    RUN_DATA=1
    RUN_ML=1
    RUN_WEB=1
    RUN_MEMORY=1
    RUN_GPU=1
    RUN_DEPS=1
fi

# Check environment
check_env

# Create benchmark data
create_benchmark_data

# Run benchmarks
if [[ "$PARALLEL" == "true" ]]; then
    # Run independent benchmarks in parallel
    pids=()
    [[ $RUN_DEPS -eq 1 ]] && benchmark_deps & pids+=($!)
    [[ $RUN_DATA -eq 1 ]] && benchmark_data & pids+=($!)
    [[ $RUN_ML -eq 1 ]] && benchmark_ml & pids+=($!)
    [[ $RUN_MEMORY -eq 1 ]] && benchmark_memory & pids+=($!)
    [[ $RUN_GPU -eq 1 ]] && benchmark_gpu & pids+=($!)
    
    # Wait for parallel benchmarks
    for pid in "${pids[@]}"; do
        wait "$pid"
    done
    
    # Run web benchmarks last (needs exclusive port)
    [[ $RUN_WEB -eq 1 ]] && benchmark_web
else
    # Run benchmarks sequentially
    [[ $RUN_DEPS -eq 1 ]] && benchmark_deps
    [[ $RUN_DATA -eq 1 ]] && benchmark_data
    [[ $RUN_ML -eq 1 ]] && benchmark_ml
    [[ $RUN_WEB -eq 1 ]] && benchmark_web
    [[ $RUN_MEMORY -eq 1 ]] && benchmark_memory
    [[ $RUN_GPU -eq 1 ]] && benchmark_gpu
fi

# Generate report
generate_report

echo -e "${GREEN}Benchmarks complete!${NC}"
echo -e "${BLUE}View results at: reports/benchmark/report.html${NC}"
