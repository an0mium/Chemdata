#!/bin/bash
# Script to install special dependencies and development tools

# Exit on error
set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Default values
INSTALL_RDKIT=true
INSTALL_OPENBABEL=true
INSTALL_INDIGO=true
INSTALL_CHEMAXON=false  # Requires license
INSTALL_CUDA=false
INSTALL_DEV=false
FORCE=false
VENV_DIR=".venv"

# Function to show usage
show_help() {
    echo "Install special dependencies and development tools"
    echo
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  --no-rdkit         Skip RDKit installation"
    echo "  --no-openbabel    Skip OpenBabel installation"
    echo "  --no-indigo       Skip Indigo installation"
    echo "  --chemaxon        Install ChemAxon tools (requires license)"
    echo "  --cuda            Install CUDA support"
    echo "  --dev             Install development tools"
    echo "  --force           Force reinstallation"
    echo "  --venv-dir DIR    Virtual environment directory (default: $VENV_DIR)"
    echo "  --help            Show this help message"
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --no-rdkit)
            INSTALL_RDKIT=false
            shift
            ;;
        --no-openbabel)
            INSTALL_OPENBABEL=false
            shift
            ;;
        --no-indigo)
            INSTALL_INDIGO=false
            shift
            ;;
        --chemaxon)
            INSTALL_CHEMAXON=true
            shift
            ;;
        --cuda)
            INSTALL_CUDA=true
            shift
            ;;
        --dev)
            INSTALL_DEV=true
            shift
            ;;
        --force)
            FORCE=true
            shift
            ;;
        --venv-dir)
            VENV_DIR="$2"
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

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Detect OS and architecture
if [[ "$OSTYPE" == "darwin"* ]]; then
    OS="macos"
    if [[ $(uname -m) == 'arm64' ]]; then
        ARCH="arm64"
    else
        ARCH="x86_64"
    fi
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    OS="linux"
    ARCH=$(uname -m)
else
    echo -e "${RED}Unsupported OS: $OSTYPE${NC}"
    exit 1
fi

# Create temporary directory
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT

# Function to install system packages
install_system_packages() {
    echo -e "${BLUE}Installing system packages...${NC}"
    
    if [ "$OS" = "macos" ]; then
        if ! command_exists brew; then
            echo -e "${YELLOW}Installing Homebrew...${NC}"
            /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        fi
        
        brew install \
            cmake \
            boost \
            openbabel \
            rdkit \
            graphviz \
            cairo \
            jpeg \
            giflib \
            libpng \
            libtiff \
            libxml2 \
            libxslt \
            libffi \
            openssl \
            zlib \
            xz \
            ncurses \
            readline \
            sqlite \
            freetype \
            pkg-config
            
    elif [ "$OS" = "linux" ]; then
        if command_exists apt-get; then
            # Ubuntu/Debian
            sudo apt-get update
            sudo apt-get install -y \
                build-essential \
                cmake \
                libboost-all-dev \
                libopenbabel-dev \
                librdkit-dev \
                python3-dev \
                python3-pip \
                python3-venv \
                graphviz \
                libgraphviz-dev \
                pkg-config \
                libcairo2-dev \
                libjpeg-dev \
                libgif-dev \
                libpng-dev \
                libtiff-dev \
                libxml2-dev \
                libxslt-dev \
                libffi-dev \
                libssl-dev \
                zlib1g-dev \
                libbz2-dev \
                liblzma-dev \
                libncurses5-dev \
                libreadline-dev \
                libsqlite3-dev \
                libfreetype6-dev
                
        elif command_exists yum; then
            # RHEL/CentOS/Fedora
            sudo yum groupinstall -y "Development Tools"
            sudo yum install -y \
                cmake \
                boost-devel \
                openbabel-devel \
                rdkit-devel \
                python3-devel \
                python3-pip \
                graphviz \
                graphviz-devel \
                pkgconfig \
                cairo-devel \
                libjpeg-turbo-devel \
                giflib-devel \
                libpng-devel \
                libtiff-devel \
                libxml2-devel \
                libxslt-devel \
                libffi-devel \
                openssl-devel \
                zlib-devel \
                bzip2-devel \
                xz-devel \
                ncurses-devel \
                readline-devel \
                sqlite-devel \
                freetype-devel
        fi
    fi
}

# Function to install uv
install_uv() {
    echo -e "${BLUE}Installing uv package manager...${NC}"
    
    if ! command_exists uv; then
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

# Function to install cheminformatics tools
install_chem_tools() {
    echo -e "${BLUE}Installing cheminformatics tools...${NC}"
    
    # Install RDKit
    if [ "$INSTALL_RDKIT" = true ]; then
        echo -e "${BLUE}Installing RDKit...${NC}"
        if [ "$OS" = "macos" ]; then
            brew install rdkit
        else
            uv pip install rdkit
        fi
    fi
    
    # Install OpenBabel
    if [ "$INSTALL_OPENBABEL" = true ]; then
        echo -e "${BLUE}Installing OpenBabel...${NC}"
        if [ "$OS" = "macos" ]; then
            brew install open-babel
        else
            sudo apt-get install -y openbabel
        fi
        uv pip install openbabel-wheel
    fi
    
    # Install Indigo
    if [ "$INSTALL_INDIGO" = true ]; then
        echo -e "${BLUE}Installing Indigo...${NC}"
        uv pip install epam.indigo
    fi
    
    # Install ChemAxon tools
    if [ "$INSTALL_CHEMAXON" = true ]; then
        echo -e "${YELLOW}ChemAxon installation requires a license${NC}"
        echo "Please download and install manually from https://chemaxon.com/"
    fi
}

# Function to install ML dependencies
install_ml_deps() {
    echo -e "${BLUE}Installing ML dependencies...${NC}"
    
    # Install PyTorch
    if [ "$INSTALL_CUDA" = true ]; then
        if [ "$OS" = "macos" ]; then
            echo -e "${YELLOW}CUDA not supported on macOS, installing CPU-only PyTorch...${NC}"
            uv pip install torch torchvision torchaudio
        else
            echo -e "${BLUE}Installing PyTorch with CUDA support...${NC}"
            uv pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
        fi
    else
        echo -e "${BLUE}Installing CPU-only PyTorch...${NC}"
        uv pip install torch torchvision torchaudio
    fi
    
    # Install other ML libraries
    uv pip install \
        "dgl>=1.0.0" \
        "dgllife>=0.3.0" \
        "deepchem>=2.7.0" \
        "mordred>=1.2.0" \
        "descriptastorus>=2.6.0" \
        "scikit-learn>=1.0.0" \
        "xgboost>=1.7.0" \
        "lightgbm>=4.0.0" \
        "optuna>=3.0.0" \
        "ray[tune]>=2.0.0"
}

# Function to install visualization dependencies
install_viz_deps() {
    echo -e "${BLUE}Installing visualization dependencies...${NC}"
    uv pip install \
        "plotly>=5.0.0" \
        "dash>=2.0.0" \
        "dash-bio>=1.0.0" \
        "dash-bootstrap-components>=1.0.0" \
        "networkx>=3.0" \
        "graphviz>=0.20.0" \
        "py3Dmol>=2.0.0" \
        "nglview>=3.0.0"
}

# Function to install web scraping dependencies
install_web_deps() {
    echo -e "${BLUE}Installing web scraping dependencies...${NC}"
    uv pip install \
        "selenium>=4.0.0" \
        "beautifulsoup4>=4.10.0" \
        "requests>=2.28.0" \
        "aiohttp>=3.8.0" \
        "playwright>=1.30.0"
}

# Function to install development tools
install_dev_tools() {
    if [ "$INSTALL_DEV" = true ]; then
        echo -e "${BLUE}Installing development tools...${NC}"
        
        # Install pre-commit hooks
        uv pip install pre-commit
        pre-commit install
        
        # Install code quality tools
        uv pip install \
            "black>=23.0.0" \
            "ruff>=0.1.0" \
            "mypy>=1.0.0" \
            "bandit>=1.7.0" \
            "safety>=2.0.0"
        
        # Install testing tools
        uv pip install \
            "pytest>=7.0.0" \
            "pytest-cov>=4.0.0" \
            "pytest-mock>=3.10.0" \
            "pytest-xdist>=3.0.0" \
            "pytest-timeout>=2.1.0" \
            "pytest-randomly>=3.12.0"
        
        # Install documentation tools
        uv pip install \
            "sphinx>=6.0.0" \
            "sphinx-rtd-theme>=1.2.0" \
            "nbsphinx>=0.9.0" \
            "jupyter>=1.0.0"
    fi
}

# Main installation process
echo -e "${BLUE}Starting installation of special dependencies...${NC}"

# Install system packages
install_system_packages

# Install uv
install_uv

# Install cheminformatics tools
install_chem_tools

# Install ML dependencies
install_ml_deps

# Install visualization dependencies
install_viz_deps

# Install web scraping dependencies
install_web_deps

# Install development tools
install_dev_tools

# Print success message
echo -e "\n${GREEN}Special dependencies installed successfully!${NC}"
echo
echo "Installed components:"
echo "  RDKit: $([ "$INSTALL_RDKIT" = true ] && echo "Yes" || echo "No")"
echo "  OpenBabel: $([ "$INSTALL_OPENBABEL" = true ] && echo "Yes" || echo "No")"
echo "  Indigo: $([ "$INSTALL_INDIGO" = true ] && echo "Yes" || echo "No")"
echo "  CUDA: $([ "$INSTALL_CUDA" = true ] && echo "Yes" || echo "No")"
echo "  Development tools: $([ "$INSTALL_DEV" = true ] && echo "Yes" || echo "No")"
echo
echo "Additional components:"
echo "  ML libraries"
echo "  Visualization libraries"
echo "  Web scraping tools"
echo
echo -e "${BLUE}Next steps:${NC}"
echo "1. Run the pipeline:"
echo "   ./scripts/run_pipeline.sh"
