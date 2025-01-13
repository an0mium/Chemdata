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
    echo "Set up development environment"
    echo
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  --minimal         Minimal installation (no ML models)"
    echo "  --no-deps        Skip system dependencies"
    echo "  --no-data        Skip initial data download"
    echo "  --no-ml          Skip ML model download"
    echo "  --cpu-only       Force CPU-only installation"
    echo "  --cuda-version VER CUDA version to install (default: 11.8)"
    echo "  --no-rdkit       Skip RDKit installation"
    echo "  --help           Show this help message"
    echo
    echo "Examples:"
    echo "  $0                Full installation"
    echo "  $0 --minimal      Minimal installation"
    echo "  $0 --no-ml       Skip ML model download"
}

# Function to check Python version
check_python() {
    echo -e "${BLUE}Checking Python version...${NC}"
    
    if ! command -v python3 &> /dev/null; then
        echo -e "${RED}Python 3 not found${NC}"
        exit 1
    fi
    
    version=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
    required="3.12"
    
    if (( $(echo "$version < $required" | bc -l) )); then
        echo -e "${RED}Python $required or higher required (found $version)${NC}"
        exit 1
    fi
    
    echo -e "${GREEN}Python $version found${NC}"
}

# Function to install uv
install_uv() {
    echo -e "${BLUE}Installing uv package manager...${NC}"
    
    if ! command -v uv &> /dev/null; then
        curl -LsSf https://astral.sh/uv/install.sh | sh
        
        # Add uv to PATH
        if [[ ":$PATH:" != *":$HOME/.local/bin:"* ]]; then
            if [[ -f ~/.zshrc ]]; then
                echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc
                source ~/.zshrc
            elif [[ -f ~/.bashrc ]]; then
                echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
                source ~/.bashrc
            fi
        fi
    else
        echo -e "${GREEN}uv already installed${NC}"
    fi
}

# Function to create virtual environment
create_venv() {
    echo -e "${BLUE}Creating virtual environment...${NC}"
    
    # Remove existing venv if it exists
    rm -rf .venv
    
    # Create new venv
    uv venv -p python3 .venv
    
    # Activate venv
    source .venv/bin/activate
}

# Function to install dependencies
install_deps() {
    echo -e "${BLUE}Installing dependencies...${NC}"
    
    # Install core dependencies
    uv pip install -e .
    
    # Install development dependencies
    uv pip install -e ".[dev]"
    
    if [[ "$MINIMAL" != "true" ]]; then
        # Install ML dependencies
        uv pip install -e ".[ml]"
        
        # Install optional dependencies
        uv pip install -e ".[optional]"
    fi
    
    # Install pre-commit hooks
    pre-commit install
}

# Function to install system dependencies
install_system_deps() {
    echo -e "${BLUE}Installing system dependencies...${NC}"
    
    # Install RDKit if requested
    if [[ "$NO_RDKIT" != "true" ]]; then
        echo -e "${BLUE}Installing RDKit...${NC}"
        if [[ "$OSTYPE" == "darwin"* ]]; then
            # macOS
            CONDA_SH="$HOME/miniforge3/etc/profile.d/conda.sh"
            if [ ! -f "$CONDA_SH" ]; then
                echo -e "${BLUE}Installing Miniforge...${NC}"
                curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh"
                bash Miniforge3-MacOSX-x86_64.sh -b
                rm Miniforge3-MacOSX-x86_64.sh
            fi
            source "$CONDA_SH"
            conda install -y -c conda-forge rdkit
        else
            # Linux
            uv pip install rdkit
        fi
    fi
    
    # Install PyTorch
    echo -e "${BLUE}Installing PyTorch...${NC}"
    
    # Function to get PyTorch wheel URL based on Python version and platform
    get_pytorch_wheel_url() {
        local py_version=$1
        local platform=$2
        local cuda_version=$3
        local is_cpu=$4
        
        # Get latest nightly build date
        local date=$(date +%Y%m%d)
        
        # Convert Python version to format used in wheel names (e.g., 3.12 -> cp312)
        local py_tag="cp${py_version/./}"
        
        if [[ "$platform" == "darwin"* ]]; then
            echo "https://download.pytorch.org/whl/nightly/cpu/torch-2.3.0.dev${date}-${py_tag}-none-macosx_11_0_arm64.whl"
        else
            if [[ "$is_cpu" == "true" ]]; then
                echo "https://download.pytorch.org/whl/nightly/cpu/torch-2.3.0.dev${date}-${py_tag}-none-linux_x86_64.whl"
            else
                echo "https://download.pytorch.org/whl/nightly/cu${cuda_version//.}/torch-2.3.0.dev${date}-${py_tag}-none-linux_x86_64.whl"
            fi
        fi
    }
    
    # Try installing PyTorch
    py_version=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
    
    if [[ "$CPU_ONLY" == "true" ]]; then
        echo -e "${BLUE}Installing PyTorch (CPU only)...${NC}"
        
        # First try installing from PyPI
        if ! uv pip install torch torchvision torchaudio; then
            echo -e "${YELLOW}Standard wheels not available, trying nightly build...${NC}"
            
            # Try downloading and installing nightly wheel
            wheel_url=$(get_pytorch_wheel_url "$py_version" "$OSTYPE" "" "true")
            if curl --output /dev/null --silent --head --fail "$wheel_url"; then
                echo -e "${BLUE}Downloading PyTorch wheel from $wheel_url${NC}"
                curl -L -o torch.whl "$wheel_url"
                if ! uv pip install torch.whl; then
                    echo -e "${RED}Failed to install PyTorch wheel${NC}"
                    rm torch.whl
                    
                    echo -e "${YELLOW}Building from source...${NC}"
                    # Install build dependencies
                    uv pip install cmake ninja
                    
                    # Clone and build PyTorch
                    git clone --recursive https://github.com/pytorch/pytorch
                    cd pytorch
                    export CMAKE_PREFIX_PATH=${CONDA_PREFIX:-"$(dirname $(which conda))/../"}
                    export BUILD_TEST=0
                    python setup.py install
                    cd ..
                    rm -rf pytorch
                else
                    rm torch.whl
                fi
            else
                echo -e "${YELLOW}Nightly wheel not available, building from source...${NC}"
                # Same build process as above
                uv pip install cmake ninja
                git clone --recursive https://github.com/pytorch/pytorch
                cd pytorch
                export CMAKE_PREFIX_PATH=${CONDA_PREFIX:-"$(dirname $(which conda))/../"}
                export BUILD_TEST=0
                python setup.py install
                cd ..
                rm -rf pytorch
            fi
        fi
        
        # Install torchvision and torchaudio
        echo -e "${BLUE}Installing torchvision and torchaudio...${NC}"
        uv pip install --pre torchvision torchaudio --index-url https://download.pytorch.org/whl/nightly/cpu
        
    else
        echo -e "${BLUE}Installing PyTorch with CUDA $CUDA_VERSION...${NC}"
        if [[ "$OSTYPE" == "darwin"* ]]; then
            # macOS (no CUDA support)
            if ! uv pip install torch torchvision torchaudio; then
                wheel_url=$(get_pytorch_wheel_url "$py_version" "$OSTYPE" "" "true")
                if curl --output /dev/null --silent --head --fail "$wheel_url"; then
                    echo -e "${BLUE}Downloading PyTorch wheel from $wheel_url${NC}"
                    curl -L -o torch.whl "$wheel_url"
                    if ! uv pip install torch.whl; then
                        echo -e "${RED}Failed to install PyTorch wheel${NC}"
                        rm torch.whl
                        
                        echo -e "${YELLOW}Building from source...${NC}"
                        uv pip install cmake ninja
                        git clone --recursive https://github.com/pytorch/pytorch
                        cd pytorch
                        export CMAKE_PREFIX_PATH=${CONDA_PREFIX:-"$(dirname $(which conda))/../"}
                        export BUILD_TEST=0
                        python setup.py install
                        cd ..
                        rm -rf pytorch
                    else
                        rm torch.whl
                    fi
                else
                    echo -e "${YELLOW}Nightly wheel not available, building from source...${NC}"
                    uv pip install cmake ninja
                    git clone --recursive https://github.com/pytorch/pytorch
                    cd pytorch
                    export CMAKE_PREFIX_PATH=${CONDA_PREFIX:-"$(dirname $(which conda))/../"}
                    export BUILD_TEST=0
                    python setup.py install
                    cd ..
                    rm -rf pytorch
                fi
            fi
            
            # Install torchvision and torchaudio
            uv pip install --pre torchvision torchaudio --index-url https://download.pytorch.org/whl/nightly/cpu
            
        else
            # Linux with CUDA
            if ! uv pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu${CUDA_VERSION//.}; then
                wheel_url=$(get_pytorch_wheel_url "$py_version" "$OSTYPE" "$CUDA_VERSION" "false")
                if curl --output /dev/null --silent --head --fail "$wheel_url"; then
                    echo -e "${BLUE}Downloading PyTorch wheel from $wheel_url${NC}"
                    curl -L -o torch.whl "$wheel_url"
                    if ! uv pip install torch.whl; then
                        echo -e "${RED}Failed to install PyTorch wheel${NC}"
                        rm torch.whl
                        
                        echo -e "${YELLOW}Building from source with CUDA support...${NC}"
                        uv pip install cmake ninja
                        git clone --recursive https://github.com/pytorch/pytorch
                        cd pytorch
                        export CMAKE_PREFIX_PATH=${CONDA_PREFIX:-"$(dirname $(which conda))/../"}
                        export BUILD_TEST=0
                        export TORCH_CUDA_ARCH_LIST="6.0 6.1 7.0 7.5 8.0 8.6"
                        export USE_CUDA=1
                        python setup.py install
                        cd ..
                        rm -rf pytorch
                    else
                        rm torch.whl
                    fi
                else
                    echo -e "${YELLOW}Nightly wheel not available, building from source with CUDA support...${NC}"
                    uv pip install cmake ninja
                    git clone --recursive https://github.com/pytorch/pytorch
                    cd pytorch
                    export CMAKE_PREFIX_PATH=${CONDA_PREFIX:-"$(dirname $(which conda))/../"}
                    export BUILD_TEST=0
                    export TORCH_CUDA_ARCH_LIST="6.0 6.1 7.0 7.5 8.0 8.6"
                    export USE_CUDA=1
                    python setup.py install
                    cd ..
                    rm -rf pytorch
                fi
            fi
            
            # Install torchvision and torchaudio
            uv pip install --pre torchvision torchaudio --index-url https://download.pytorch.org/whl/nightly/cu${CUDA_VERSION//.}
        fi
    fi
    
    # Make special deps script executable
    chmod +x scripts/install_special_deps.sh
    
    # Run special deps installation
    ./scripts/install_special_deps.sh
}

# Function to create config files
create_config() {
    echo -e "${BLUE}Creating configuration files...${NC}"
    
    # Create config directory
    mkdir -p ~/.config/chemdata
    
    # Create default config if it doesn't exist
    if [ ! -f ~/.config/chemdata/config.json ]; then
        cat > ~/.config/chemdata/config.json << EOL
{
    "data_dir": "$(pwd)/data",
    "model_dir": "$(pwd)/models",
    "cache_dir": "$(pwd)/.cache",
    "log_dir": "$(pwd)/logs",
    "settings": {
        "workers": 4,
        "batch_size": 100,
        "use_gpu": ${CPU_ONLY:=false},
        "cache_ttl": 3600
    },
    "ml": {
        "model_batch_size": 32,
        "prediction_confidence_threshold": 0.7,
        "enable_predictions": ${MINIMAL:=false}
    },
    "web": {
        "host": "localhost",
        "port": 8050,
        "debug": false,
        "cache_type": "filesystem"
    },
    "api": {
        "pubchem_rate_limit": 5,
        "swiss_rate_limit": 2,
        "community_rate_limit": 1,
        "social_rate_limit": 0.2
    },
    "credentials": {}
}
EOL
    fi
    
    # Create .env file if it doesn't exist
    if [ ! -f .env ]; then
        cat > .env << EOL
# API Credentials
REDDIT_CLIENT_ID=
REDDIT_CLIENT_SECRET=
TWITTER_API_KEY=
TWITTER_API_SECRET=

# Directories
DATA_DIR=./data
MODEL_DIR=./models
CACHE_DIR=./.cache
LOG_DIR=./logs

# Processing Settings
WORKERS=4
BATCH_SIZE=100
USE_GPU=${CPU_ONLY:=false}

# ML Settings
MODEL_BATCH_SIZE=32
PREDICTION_CONFIDENCE_THRESHOLD=0.7
ENABLE_PREDICTIONS=${MINIMAL:=false}

# Web Settings
FLASK_DEBUG=false
FLASK_HOST=localhost
FLASK_PORT=8050
EOL
    fi
}

# Function to create project directories
create_dirs() {
    echo -e "${BLUE}Creating project directories...${NC}"
    
    # Create data directories
    mkdir -p data/{raw,processed,interim,external}
    
    # Create model directories
    mkdir -p models/{toxicity,abuse,activity,affinity}
    
    # Create other directories
    mkdir -p logs
    mkdir -p .cache
    mkdir -p reports/{coverage,test-results,profiling}
    mkdir -p web/{static,templates}
}

# Function to download ML models
download_models() {
    if [[ "$SKIP_ML" == "true" || "$MINIMAL" == "true" ]]; then
        echo -e "${YELLOW}Skipping ML model download${NC}"
        return
    fi
    
    echo -e "${BLUE}Downloading ML models...${NC}"
    
    # Create models directory
    mkdir -p models
    
    # Download models from cloud storage
    models=(
        "toxicity_predictor.pt"
        "abuse_predictor.pt"
        "activity_predictor.pt"
        "affinity_predictor.pt"
        "gnn_model.pt"
        "ensemble_model.pt"
    )
    
    for model in "${models[@]}"; do
        if [ ! -f "models/$model" ]; then
            echo -e "${BLUE}Downloading $model...${NC}"
            curl -L "https://storage.googleapis.com/chemdata-models/$model" \
                -o "models/$model"
        else
            echo -e "${GREEN}$model already exists${NC}"
        fi
    done
}

# Function to download initial data
download_initial_data() {
    if [[ "$SKIP_DATA" == "true" ]]; then
        echo -e "${YELLOW}Skipping initial data download${NC}"
        return
    fi
    
    echo -e "${BLUE}Downloading initial data...${NC}"
    
    # Create data directory
    mkdir -p data/raw
    
    # Download BindingDB data
    if [ ! -f data/raw/BindingDB_All.tsv ]; then
        echo -e "${BLUE}Downloading BindingDB data...${NC}"
        curl -L "https://www.bindingdb.org/bind/downloads/BindingDB_All.tsv.zip" \
            -o data/raw/BindingDB_All.tsv.zip
        unzip data/raw/BindingDB_All.tsv.zip -d data/raw
        rm data/raw/BindingDB_All.tsv.zip
    else
        echo -e "${GREEN}BindingDB data already exists${NC}"
    fi
    
    # Download ChEMBL data
    if [ ! -f data/raw/chembl_targets.csv ]; then
        echo -e "${BLUE}Downloading ChEMBL data...${NC}"
        curl -L "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_target_dictionary.txt" \
            -o data/raw/chembl_targets.csv
    else
        echo -e "${GREEN}ChEMBL data already exists${NC}"
    fi
    
    # Validate downloaded data
    echo -e "${BLUE}Validating downloaded data...${NC}"
    
    if [ ! -f data/raw/BindingDB_All.tsv ]; then
        echo -e "${RED}BindingDB data not found${NC}"
        exit 1
    fi
    
    if [ ! -f data/raw/chembl_targets.csv ]; then
        echo -e "${RED}ChEMBL data not found${NC}"
        exit 1
    fi
}

# Function to initialize git
init_git() {
    echo -e "${BLUE}Initializing git repository...${NC}"
    
    if [ ! -d .git ]; then
        git init
        
        # Create .gitignore if it doesn't exist
        if [ ! -f .gitignore ]; then
            cat > .gitignore << EOL
# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
build/
develop-eggs/
dist/
downloads/
eggs/
.eggs/
lib/
lib64/
parts/
sdist/
var/
wheels/
*.egg-info/
.installed.cfg
*.egg

# Virtual Environment
.env
.venv/
venv/
ENV/

# IDE
.idea/
.vscode/
*.swp
*.swo

# Data
data/raw/
data/processed/
data/interim/
data/external/
*.tsv
*.csv
*.json
!package.json
!tsconfig.json

# Models
models/*/
*.pt
*.pth
*.h5
*.ckpt

# Logs
logs/
*.log

# Cache
.cache/
.pytest_cache/
.mypy_cache/
.ruff_cache/
.coverage
htmlcov/
EOL
        fi
        
        # Set up git hooks
        pre-commit install
    fi
}

# Parse command line arguments
MINIMAL=false
SKIP_DEPS=false
SKIP_DATA=false
SKIP_ML=false
CPU_ONLY=false
NO_RDKIT=false
CUDA_VERSION="11.8"

while [[ $# -gt 0 ]]; do
    case $1 in
        --minimal)
            MINIMAL=true
            shift
            ;;
        --no-deps)
            SKIP_DEPS=true
            shift
            ;;
        --no-data)
            SKIP_DATA=true
            shift
            ;;
        --no-ml)
            SKIP_ML=true
            shift
            ;;
        --cpu-only)
            CPU_ONLY=true
            shift
            ;;
        --cuda-version)
            CUDA_VERSION="$2"
            shift 2
            ;;
        --no-rdkit)
            NO_RDKIT=true
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

# Run setup steps
echo -e "${BLUE}Setting up development environment...${NC}"

# Check Python version
check_python

# Install uv
install_uv

# Create virtual environment
create_venv

# Install system dependencies if not skipped
if [[ "$SKIP_DEPS" != "true" ]]; then
    install_system_deps
fi

# Install dependencies
install_deps

# Create config files
create_config

# Create project directories
create_dirs

# Download initial data
download_initial_data

# Download ML models
download_models

# Initialize git
init_git

# Print success message
echo -e "\n${GREEN}Development environment setup complete!${NC}"

# Print next steps
echo -e "\n${BLUE}Next steps:${NC}"
echo "1. Add API credentials to .env file"
echo "2. Run tests: pytest"
echo "3. Start development server: python -m binding_data_processor.web.app"
echo "4. Visit http://localhost:8050 in your browser"

# Print warnings based on installation type
if [[ "$MINIMAL" == "true" ]]; then
    echo -e "\n${YELLOW}Note: Minimal installation completed${NC}"
    echo "Some features requiring ML models will not be available"
    echo "Run without --minimal flag for full installation"
fi

if [[ "$CPU_ONLY" == "true" ]]; then
    echo -e "\n${YELLOW}Note: CPU-only installation completed${NC}"
    echo "GPU acceleration will not be available"
    echo "Run without --cpu-only flag to enable GPU support"
fi
