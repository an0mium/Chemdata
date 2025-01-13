#!/bin/bash
# Script to install special dependencies and development tools

# Exit on error
set -e

# Detect OS
case "$OSTYPE" in
    darwin*)  OS="macos" ;;
    linux*)   OS="linux" ;;
    *)        echo "Unsupported OS: $OSTYPE"; exit 1 ;;
esac

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

# Function to detect CPU architecture
get_cpu_arch() {
    if [[ "$(uname -m)" == "arm64" ]]; then
        echo "arm64"
    else
        echo "x86_64"
    fi
}

# Function to set platform-specific compiler flags
set_compiler_flags() {
    local arch=$(get_cpu_arch)
    if [[ "$OSTYPE" == "darwin"* && "$arch" == "arm64" ]]; then
        # Apple Silicon specific flags
        export CC=/usr/bin/clang
        export CXX=/usr/bin/clang++
        export ARCHFLAGS="-arch arm64"
        
        # Base flags
        export CFLAGS="-O3 -arch arm64"
        export CXXFLAGS="-O3 -arch arm64"
        export LDFLAGS=""
        
        # Include paths
        export CFLAGS="$CFLAGS -I/opt/homebrew/include"
        export CXXFLAGS="$CXXFLAGS -I/opt/homebrew/include"
        export LDFLAGS="$LDFLAGS -L/opt/homebrew/lib"
        
        # Python paths
        local PYTHON_INCLUDE=$(python3-config --includes)
        export CFLAGS="$CFLAGS $PYTHON_INCLUDE"
        export CXXFLAGS="$CXXFLAGS $PYTHON_INCLUDE"
        
        # OpenMP support
        export CPPFLAGS="-Xpreprocessor -fopenmp"
        export CFLAGS="$CFLAGS -I/opt/homebrew/opt/libomp/include"
        export CXXFLAGS="$CXXFLAGS -I/opt/homebrew/opt/libomp/include"
        export LDFLAGS="$LDFLAGS -L/opt/homebrew/opt/libomp/lib -lomp"
        
        # ARM64 optimizations
        export CFLAGS="$CFLAGS -mcpu=apple-a14 -mtune=native"
        export CXXFLAGS="$CXXFLAGS -mcpu=apple-a14 -mtune=native"
        
        # Disable x86 instructions
        export CFLAGS="$CFLAGS -mno-avx -mno-avx2 -mno-sse4.2"
        export CXXFLAGS="$CXXFLAGS -mno-avx -mno-avx2 -mno-sse4.2"
        
        # Warning controls
        export CFLAGS="$CFLAGS -Wno-deprecated-declarations -Wno-unreachable-code -Wno-unused-function"
        export CXXFLAGS="$CXXFLAGS -Wno-deprecated-declarations -Wno-unreachable-code -Wno-unused-function"
        
        # Python compatibility
        export CFLAGS="$CFLAGS -DPY_SSIZE_T_CLEAN -DCYTHON_COMPILING_IN_CPYTHON"
        export CXXFLAGS="$CXXFLAGS -DPY_SSIZE_T_CLEAN -DCYTHON_COMPILING_IN_CPYTHON"
        
        # Additional ARM64 flags
        export CFLAGS="$CFLAGS -DAPPLE_ARM64=1"
        export CXXFLAGS="$CXXFLAGS -DAPPLE_ARM64=1"
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        # Linux flags
        export CFLAGS="-O3"
        export CXXFLAGS="-O3"
        if [[ "$(get_cpu_arch)" == "x86_64" ]]; then
            # Enable AVX2 only on x86_64
            export CFLAGS="$CFLAGS -march=native"
            export CXXFLAGS="$CXXFLAGS -march=native"
        fi
    fi
}

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
            open-babel \
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
            pkg-config \
            libomp
            
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
                libfreetype6-dev \
                libomp-dev
                
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
                freetype-devel \
                libomp-devel
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

# Function to install C extension packages
install_c_extensions() {
    echo -e "${BLUE}Installing packages with C extensions...${NC}"
    
    # Set compiler flags with additional compatibility flags
    set_compiler_flags
    
    # Get Python include paths and add them directly
    PYTHON_INCLUDES=$(python3-config --includes)
    export CFLAGS="$CFLAGS $PYTHON_INCLUDES"
    export CXXFLAGS="$CXXFLAGS $PYTHON_INCLUDES"
    
    # Add Python 3.10+ compatibility flags
    export CFLAGS="$CFLAGS -DPY_SSIZE_T_CLEAN -DCYTHON_COMPILING_IN_CPYTHON -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION"
    export CXXFLAGS="$CXXFLAGS -DPY_SSIZE_T_CLEAN -DCYTHON_COMPILING_IN_CPYTHON -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION"
    
    # Disable deprecation warnings and enable additional compatibility flags
    export CFLAGS="$CFLAGS -Wno-deprecated-declarations -Wno-unreachable-code -Wno-unused-function -Wno-#warnings -Wno-error=implicit-function-declaration"
    export CXXFLAGS="$CXXFLAGS -Wno-deprecated-declarations -Wno-unreachable-code -Wno-unused-function -Wno-#warnings -Wno-error=implicit-function-declaration"
    
    # Install base dependencies with specific versions
    pip install --no-deps 'setuptools>=65.0.0' 'wheel>=0.38.0'
    
    # Install numpy first as it's required by many packages
    pip install 'numpy>=1.24.0,<2.0.0'
    
    # Define package versions compatible with Python 3.10+ on ARM64
    PACKAGES=(
        "cymem>=2.0.7,<3.0.0"
        "preshed>=3.0.8,<4.0.0"
        "murmurhash>=1.0.9,<2.0.0"
        "wasabi>=1.1.2,<2.0.0"
        "srsly>=2.4.7,<3.0.0"
        "plac>=1.4.0,<2.0.0"
        "thinc>=8.1.10,<9.0.0"
        "blis>=0.7.11,<0.8.0"
    )

    # Add Python 3.10+ and ARM64 specific flags
    if python3 -c "import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)" 2>/dev/null; then
        # Python 3.10+ flags
        export CFLAGS="$CFLAGS -DPY_SSIZE_T_CLEAN -DCYTHON_COMPILING_IN_CPYTHON -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION"
        export CXXFLAGS="$CXXFLAGS -DPY_SSIZE_T_CLEAN -DCYTHON_COMPILING_IN_CPYTHON -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION"
        # Unicode handling flags
        export CFLAGS="$CFLAGS -DPY_UNICODE_WIDE -DPyUnicode_GET_SIZE=PyUnicode_GET_LENGTH -DPyUnicode_WSTR_LENGTH=PyUnicode_GET_LENGTH"
        export CXXFLAGS="$CXXFLAGS -DPY_UNICODE_WIDE -DPyUnicode_GET_SIZE=PyUnicode_GET_LENGTH -DPyUnicode_WSTR_LENGTH=PyUnicode_GET_LENGTH"
        # ARM64 specific flags
        if [[ "$(uname -m)" == "arm64" ]]; then
            export ARCHFLAGS="-arch arm64"
            export CFLAGS="$CFLAGS -arch arm64"
            export CXXFLAGS="$CXXFLAGS -arch arm64"
            # Disable AVX instructions
            export CFLAGS="$CFLAGS -mno-avx -mno-avx2"
            export CXXFLAGS="$CXXFLAGS -mno-avx -mno-avx2"
        fi
    fi
    
    # Function to patch source before building
    patch_source() {
        local package=$1
        local src_dir=$2
        
        # Apply common patches for all packages
        find "$src_dir" -type f -name "*.c" -o -name "*.cpp" -exec sed -i.bak '
            # Reference counting fixes
            s/++Py_REFCNT(\([^)]*\))/Py_INCREF(\1)/g
            s/--Py_REFCNT(\([^)]*\))/Py_DECREF(\1)/g
            s/Py_REFCNT(\([^)]*\))++/Py_INCREF(\1)/g
            s/Py_REFCNT(\([^)]*\))--/Py_DECREF(\1)/g
            # Python 3.10+ compatibility fixes
            s/PyCode_New([^)]*)/PyCode_NewEmpty("", "", 0)/g
            s/_PyGen_Send/PyIter_Send/g
            s/PyUnicode_GET_SIZE/PyUnicode_GET_LENGTH/g
            s/PyUnicode_WSTR_LENGTH/PyUnicode_GET_LENGTH/g
            s/PyUnicode_AsUnicode/PyUnicode_AsUTF8/g
            # Remove deprecated tp_print
            s/[^>]tp_print/tp_repr/g
            # Fix other deprecated APIs
            s/PyInt_/PyLong_/g
            s/PyString_/PyBytes_/g
            s/Py_TPFLAGS_HAVE_WEAKREFS/0/g
            s/Py_TPFLAGS_HAVE_ITER/0/g
        ' {} +
        
        # Additional package-specific patches
        case "$package" in
            preshed*|thinc*)
                # Fix memory management
                find "$src_dir" -type f -name "*.c" -o -name "*.cpp" -exec sed -i.bak '
                    s/PyMem_Malloc/PyMem_RawMalloc/g
                    s/PyMem_Realloc/PyMem_RawRealloc/g
                    s/PyMem_Free/PyMem_RawFree/g
                ' {} +
                ;;
            blis*)
                # Remove AVX instructions for ARM64
                if [[ "$(uname -m)" == "arm64" ]]; then
                    find "$src_dir" -type f -name "*.c" -exec sed -i.bak 's/__AVX__/0/g' {} +
                fi
                ;;
        esac
    }
    
    # Try installing packages with patching if needed
    for package in "${PACKAGES[@]}"; do
        echo "Installing $package..."
        
        # Try pre-built wheel first
        if pip install --only-binary :all: "$package" 2>/dev/null; then
            echo "Installed pre-built wheel for $package"
            continue
        fi
        
        echo "Pre-built wheel not available for $package, building from source..."
        
        # Create temp directory for source
        pkg_temp_dir=$(mktemp -d)
        pkg_name=$(echo "$package" | cut -d= -f1)
        
        # Download source
        pip download --no-binary :all: "$package" -d "$pkg_temp_dir"
        
        # Extract source
        cd "$pkg_temp_dir"
        tar xf "$pkg_name"*.tar.gz || unzip "$pkg_name"*.zip
        cd "$pkg_name"*
        
        # Patch source files
        patch_source "$package" "."
        
        # Build and install
        ARCHFLAGS="-arch $(get_cpu_arch)" pip install --no-deps --no-binary :all: .
        
        # Clean up
        cd -
        rm -rf "$pkg_temp_dir"
    done
    
    # Install spacy with compatible version
    pip install --no-deps 'spacy>=3.4.1,<3.5.0'
    
    # Install scispacy and dependencies
    pip install --no-deps 'scispacy==0.4.0'
    pip install --no-deps 'scipy<1.11' 'requests>=2.0.0,<3.0.0' 'conllu' 'numpy' 'joblib' 'scikit-learn>=0.20.3' 'pysbd'
    
    # Install atproto dependency
    pip install --no-deps 'atproto>=0.0.20'
    
    # Install en-core-sci-lg model
    pip install --no-deps https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_core_sci_lg-0.5.1.tar.gz
    
    echo "All packages installed successfully"
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
            sudo apt-get install -y libopenbabel-dev
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

# Function to install nmslib
install_nmslib() {
    echo -e "${BLUE}Installing nmslib...${NC}"
    
    # Install nmslib-metabrainz version 2.1.3 or later
    if pip install --no-deps "nmslib-metabrainz>=2.1.3"; then
        echo -e "${GREEN}Successfully installed nmslib-metabrainz${NC}"
        return 0
    fi
    
    echo "MetaBrainz installation failed, trying conda..."
    if command_exists conda; then
        echo "Attempting to install nmslib via conda..."
        if conda install -y -c conda-forge nmslib; then
            echo "Successfully installed nmslib via conda"
            return 0
        fi
    fi
    
    echo "Trying original nmslib..."
    if [[ "$OSTYPE" == "darwin"* && "$(get_cpu_arch)" == "arm64" ]]; then
        echo "On Apple Silicon, installing with optimized flags..."
        # Set compiler flags
        set_compiler_flags
        # Try installing with optimized flags
        if pip install --no-deps --no-binary :all: nmslib; then
            echo "Successfully installed nmslib from source"
            return 0
        fi
        echo "Failed to install nmslib. Please use alternative nearest neighbor implementations like scikit-learn or annoy"
        return 1
    else
        if pip install --no-deps --no-cache-dir nmslib; then
            echo "Successfully installed nmslib"
            return 0
        fi
        echo "Failed to install nmslib"
    return 1
    fi
}

# Function to install ML dependencies
install_ml_deps() {
    echo -e "${BLUE}Installing ML dependencies...${NC}"
    
    # Install nmslib first
    install_nmslib
    
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

# Install C extensions with special handling
install_c_extensions

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
