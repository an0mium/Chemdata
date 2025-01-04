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
    echo "Check development environment dependencies"
    echo
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  --minimal      Skip ML model checks"
    echo "  --no-gpu      Skip GPU checks"
    echo "  --no-net      Skip network checks"
    echo "  --help        Show this help message"
}

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check Python version
check_python() {
    echo -e "${BLUE}Checking Python...${NC}"
    
    if ! command_exists python3; then
        echo -e "${RED}Python 3 not found${NC}"
        return 1
    fi
    
    version=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
    required="3.10"
    
    if (( $(echo "$version < $required" | bc -l) )); then
        echo -e "${RED}Python $required or higher required (found $version)${NC}"
        return 1
    fi
    
    echo -e "${GREEN}Python $version found${NC}"
    return 0
}

# Function to check package manager
check_package_manager() {
    echo -e "${BLUE}Checking package manager...${NC}"
    
    if command_exists apt-get; then
        echo -e "${GREEN}apt package manager found${NC}"
        PKG_MANAGER="apt-get"
    elif command_exists brew; then
        echo -e "${GREEN}Homebrew package manager found${NC}"
        PKG_MANAGER="brew"
    elif command_exists pacman; then
        echo -e "${GREEN}pacman package manager found${NC}"
        PKG_MANAGER="pacman"
    else
        echo -e "${RED}No supported package manager found${NC}"
        return 1
    fi
    
    return 0
}

# Function to check system dependencies
check_system_deps() {
    echo -e "${BLUE}Checking system dependencies...${NC}"
    local missing=()
    
    # OpenGL dependencies
    if ! command_exists glxinfo; then
        missing+=("mesa-utils")
    fi
    
    # Graphviz
    if ! command_exists dot; then
        missing+=("graphviz")
    fi
    
    # OpenBabel
    if ! command_exists obabel; then
        missing+=("openbabel")
    fi
    
    # PostgreSQL
    if ! command_exists psql; then
        missing+=("postgresql")
    fi
    
    # Additional system tools
    for tool in curl wget git jq bc; do
        if ! command_exists "$tool"; then
            missing+=("$tool")
        fi
    done
    
    if [ ${#missing[@]} -gt 0 ]; then
        echo -e "${YELLOW}Missing system packages: ${missing[*]}${NC}"
        echo "Install with:"
        case $PKG_MANAGER in
            apt-get)
                echo "sudo apt-get install ${missing[*]}"
                ;;
            brew)
                echo "brew install ${missing[*]}"
                ;;
            pacman)
                echo "sudo pacman -S ${missing[*]}"
                ;;
        esac
        return 1
    fi
    
    echo -e "${GREEN}All system dependencies found${NC}"
    return 0
}

# Function to check Python dependencies
check_python_deps() {
    echo -e "${BLUE}Checking Python dependencies...${NC}"
    local missing=()
    
    # Core dependencies
    for pkg in numpy pandas scipy scikit-learn rdkit torch flask dash plotly requests beautifulsoup4 praw tweepy; do
        if ! python3 -c "import $pkg" 2>/dev/null; then
            missing+=("$pkg")
        fi
    done
    
    # ML dependencies
    if [[ "$MINIMAL" != "true" ]]; then
        for pkg in transformers deepchem torch_geometric; do
            if ! python3 -c "import $pkg" 2>/dev/null; then
                missing+=("$pkg")
            fi
        done
    fi
    
    if [ ${#missing[@]} -gt 0 ]; then
        echo -e "${YELLOW}Missing Python packages: ${missing[*]}${NC}"
        echo "Install with:"
        echo "uv pip install ${missing[*]}"
        return 1
    fi
    
    echo -e "${GREEN}All Python dependencies found${NC}"
    return 0
}

# Function to check development tools
check_dev_tools() {
    echo -e "${BLUE}Checking development tools...${NC}"
    local missing=()
    
    # Check git
    if ! command_exists git; then
        missing+=("git")
    fi
    
    # Check uv
    if ! command_exists uv; then
        echo -e "${YELLOW}uv not found${NC}"
        echo "Install with: curl -LsSf https://astral.sh/uv/install.sh | sh"
        return 1
    fi
    
    # Check pre-commit
    if ! command_exists pre-commit; then
        missing+=("pre-commit")
    fi
    
    # Check development packages
    for pkg in pytest black mypy flake8 isort; do
        if ! python3 -c "import $pkg" 2>/dev/null; then
            missing+=("$pkg")
        fi
    done
    
    if [ ${#missing[@]} -gt 0 ]; then
        echo -e "${YELLOW}Missing development tools: ${missing[*]}${NC}"
        return 1
    fi
    
    echo -e "${GREEN}All development tools found${NC}"
    return 0
}

# Function to check environment
check_environment() {
    echo -e "${BLUE}Checking environment...${NC}"
    
    # Check virtual environment
    if [[ -z "${VIRTUAL_ENV}" ]]; then
        echo -e "${YELLOW}Not running in a virtual environment${NC}"
        echo "Create and activate one with:"
        echo "uv venv -p python3.12 .venv"
        echo "source .venv/bin/activate"
        return 1
    fi
    
    # Check environment variables
    local missing=()
    for var in REDDIT_CLIENT_ID REDDIT_CLIENT_SECRET TWITTER_API_KEY TWITTER_API_SECRET; do
        if [[ -z "${!var}" ]]; then
            missing+=("$var")
        fi
    done
    
    if [ ${#missing[@]} -gt 0 ]; then
        echo -e "${YELLOW}Missing environment variables: ${missing[*]}${NC}"
        echo "Add them to .env file"
        return 1
    fi
    
    echo -e "${GREEN}Environment setup correctly${NC}"
    return 0
}

# Function to check data files
check_data_files() {
    echo -e "${BLUE}Checking data files...${NC}"
    
    # Check data directories
    for dir in data/{raw,processed,interim,external} models/{toxicity,abuse,activity,affinity} logs .cache; do
        if [ ! -d "$dir" ]; then
            echo -e "${YELLOW}Directory not found: $dir${NC}"
            echo "Create with: mkdir -p $dir"
            return 1
        fi
    done
    
    # Check required data files
    local missing=()
    for file in data/raw/BindingDB_All.tsv data/raw/chembl_targets.csv; do
        if [ ! -f "$file" ]; then
            missing+=("$file")
        fi
    done
    
    if [ ${#missing[@]} -gt 0 ]; then
        echo -e "${YELLOW}Missing data files: ${missing[*]}${NC}"
        echo "Run setup_dev.sh to download required data"
        return 1
    fi
    
    echo -e "${GREEN}All data files present${NC}"
    return 0
}

# Function to check ML models
check_ml_models() {
    if [[ "$MINIMAL" == "true" ]]; then
        echo -e "${YELLOW}Skipping ML model checks${NC}"
        return 0
    fi
    
    echo -e "${BLUE}Checking ML models...${NC}"
    
    # Check required model files
    local missing=()
    for model in toxicity_predictor abuse_predictor activity_predictor affinity_predictor gnn_model ensemble_model; do
        if [ ! -f "models/$model.pt" ]; then
            missing+=("$model.pt")
        fi
    done
    
    if [ ${#missing[@]} -gt 0 ]; then
        echo -e "${YELLOW}Missing ML models: ${missing[*]}${NC}"
        echo "Run setup_dev.sh without --no-ml flag to download models"
        return 1
    fi
    
    echo -e "${GREEN}All ML models present${NC}"
    return 0
}

# Function to check GPU support
check_gpu() {
    if [[ "$NO_GPU" == "true" ]]; then
        echo -e "${YELLOW}Skipping GPU checks${NC}"
        return 0
    fi
    
    echo -e "${BLUE}Checking GPU support...${NC}"
    
    # Check NVIDIA GPU
    if command_exists nvidia-smi; then
        echo -e "${GREEN}NVIDIA GPU found${NC}"
        
        # Check CUDA
        if python3 -c "import torch; print(torch.cuda.is_available())" 2>/dev/null | grep -q "True"; then
            echo -e "${GREEN}CUDA support enabled${NC}"
            
            # Check GPU memory
            local gpu_mem=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -n 1)
            if (( gpu_mem < 8000 )); then
                echo -e "${YELLOW}Warning: GPU has less than 8GB memory${NC}"
            fi
            
            return 0
        else
            echo -e "${YELLOW}CUDA support not enabled${NC}"
            echo "Install PyTorch with CUDA support"
            return 1
        fi
    else
        echo -e "${YELLOW}No NVIDIA GPU found${NC}"
        echo "GPU acceleration will not be available"
        return 1
    fi
}

# Function to check network connectivity
check_network() {
    if [[ "$NO_NET" == "true" ]]; then
        echo -e "${YELLOW}Skipping network checks${NC}"
        return 0
    fi
    
    echo -e "${BLUE}Checking network connectivity...${NC}"
    local status=0
    
    # Check API endpoints
    local endpoints=(
        "https://api.github.com"
        "https://api.twitter.com"
        "https://oauth.reddit.com"
        "https://www.bindingdb.org"
        "https://www.ebi.ac.uk"
    )
    
    for endpoint in "${endpoints[@]}"; do
        if ! curl --silent --head "$endpoint" >/dev/null; then
            echo -e "${YELLOW}Cannot reach $endpoint${NC}"
            status=1
        fi
    done
    
    if [ $status -eq 0 ]; then
        echo -e "${GREEN}All network checks passed${NC}"
    fi
    
    return $status
}

# Parse command line arguments
MINIMAL=false
NO_GPU=false
NO_NET=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --minimal)
            MINIMAL=true
            shift
            ;;
        --no-gpu)
            NO_GPU=true
            shift
            ;;
        --no-net)
            NO_NET=true
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

# Main function
main() {
    local status=0
    
    echo -e "${BLUE}Running dependency checks...${NC}"
    
    # Run all checks
    check_python || status=1
    check_package_manager || status=1
    check_system_deps || status=1
    check_python_deps || status=1
    check_dev_tools || status=1
    check_environment || status=1
    check_data_files || status=1
    check_ml_models || status=1
    check_gpu || status=1
    check_network || status=1
    
    # Print summary
    echo
    if [ $status -eq 0 ]; then
        echo -e "${GREEN}All dependency checks passed!${NC}"
    else
        echo -e "${RED}Some checks failed. See above for details.${NC}"
        echo "Run setup_dev.sh to fix missing dependencies"
    fi
    
    return $status
}

# Run main function
main
