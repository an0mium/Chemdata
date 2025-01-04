#!/bin/bash
# Script to set up the project environment and dependencies

# Exit on error
set -e

# Default values
PYTHON_VERSION="3.8"
VENV_DIR="venv"
INSTALL_CUDA=false
INSTALL_RDKIT=true
INSTALL_SPECIAL=true
DEV_MODE=false
FORCE=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --python-version)
            PYTHON_VERSION="$2"
            shift 2
            ;;
        --venv-dir)
            VENV_DIR="$2"
            shift 2
            ;;
        --cuda)
            INSTALL_CUDA=true
            shift
            ;;
        --no-rdkit)
            INSTALL_RDKIT=false
            shift
            ;;
        --no-special)
            INSTALL_SPECIAL=false
            shift
            ;;
        --dev)
            DEV_MODE=true
            shift
            ;;
        --force)
            FORCE=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if virtual environment exists
if [ -d "$VENV_DIR" ] && [ "$FORCE" = false ]; then
    echo "Virtual environment already exists at $VENV_DIR"
    echo "Use --force to recreate it"
    exit 1
fi

# Check Python version
CURRENT_PYTHON_VERSION=$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
if [ "$(printf '%s\n' "$PYTHON_VERSION" "$CURRENT_PYTHON_VERSION" | sort -V | head -n1)" != "$PYTHON_VERSION" ]; then
    echo "Error: Python $PYTHON_VERSION or higher is required (found $CURRENT_PYTHON_VERSION)"
    exit 1
fi

# Create project directories
echo "Creating project directories..."
mkdir -p data/{raw,processed,interim,external}
mkdir -p models/{activity,toxicity,abuse,bbb}
mkdir -p reports/{figures,tables,html,pdf}
mkdir -p logs
mkdir -p cache
mkdir -p checkpoints
mkdir -p notebooks
mkdir -p docs
mkdir -p tests/data

# Create/recreate virtual environment
echo "Creating virtual environment..."
if [ -d "$VENV_DIR" ]; then
    rm -rf "$VENV_DIR"
fi
python3 -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip setuptools wheel

# Install base dependencies
echo "Installing base dependencies..."
pip install -r requirements.txt

# Install development dependencies
if [ "$DEV_MODE" = true ]; then
    echo "Installing development dependencies..."
    pip install -r requirements-dev.txt
fi

# Install RDKit
if [ "$INSTALL_RDKIT" = true ]; then
    echo "Installing RDKit..."
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS
        brew install rdkit
    else
        # Linux
        conda install -c conda-forge rdkit
    fi
fi

# Install CUDA if requested
if [ "$INSTALL_CUDA" = true ]; then
    echo "Installing CUDA dependencies..."
    pip install torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu116
fi

# Install special dependencies
if [ "$INSTALL_SPECIAL" = true ]; then
    echo "Installing special dependencies..."
    ./scripts/install_special_deps.sh
fi

# Set up pre-commit hooks if in dev mode
if [ "$DEV_MODE" = true ]; then
    echo "Setting up pre-commit hooks..."
    pip install pre-commit
    pre-commit install
fi

# Create .env file if it doesn't exist
if [ ! -f .env ]; then
    echo "Creating .env file..."
    cat > .env << EOL
# Project settings
PROJECT_ROOT=$(pwd)
PYTHONPATH=\${PROJECT_ROOT}
DATA_DIR=\${PROJECT_ROOT}/data
MODEL_DIR=\${PROJECT_ROOT}/models
CACHE_DIR=\${PROJECT_ROOT}/cache
LOG_DIR=\${PROJECT_ROOT}/logs

# API keys (replace with your keys)
CHEMBL_API_KEY=
PUBCHEM_API_KEY=
SWISS_API_KEY=
REDDIT_CLIENT_ID=
REDDIT_CLIENT_SECRET=
TWITTER_API_KEY=
TWITTER_API_SECRET=

# ML settings
USE_GPU=false
BATCH_SIZE=32
NUM_WORKERS=4

# Web settings
HOST=localhost
PORT=8000
DEBUG=true
EOL
fi

# Initialize git if not already
if [ ! -d .git ]; then
    echo "Initializing git repository..."
    git init
    git add .
    git commit -m "Initial commit"
fi

# Print success message
echo
echo "Project setup completed successfully!"
echo
echo "Project structure:"
echo "  data/            - Data files"
echo "    raw/          - Raw data"
echo "    processed/    - Processed data"
echo "    interim/      - Intermediate data"
echo "    external/     - External data"
echo "  models/         - ML models"
echo "  reports/        - Generated reports"
echo "  logs/          - Log files"
echo "  cache/         - Cache files"
echo "  checkpoints/   - Checkpoints"
echo "  notebooks/     - Jupyter notebooks"
echo "  docs/          - Documentation"
echo "  tests/         - Test files"
echo
echo "Virtual environment: $VENV_DIR"
echo "Python version: $CURRENT_PYTHON_VERSION"
echo
echo "Next steps:"
echo "1. Edit .env file to add your API keys"
echo "2. Run the pipeline:"
echo "   ./scripts/run_pipeline.sh"
echo
echo "For development:"
echo "1. Install dev dependencies:"
echo "   ./scripts/setup_project.sh --dev"
echo "2. Run tests:"
echo "   pytest"
echo "3. Start coding!"
