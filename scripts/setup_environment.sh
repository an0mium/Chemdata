#!/bin/bash
set -e

echo "Setting up development environment..."

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check if a Python package is installed
package_installed() {
    local package_name=$1
    python -c "import ${package_name%%[<>=]*}" 2>/dev/null
    return $?
}

# Function to check package version compatibility
check_package_version() {
    local package=$1
    local min_version=$2
    python -c "import pkg_resources; pkg_resources.require('$package>=$min_version')" 2>/dev/null
    return $?
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
        export CFLAGS="-O3 -arch arm64 -I/opt/homebrew/include"
        export CXXFLAGS="-O3 -arch arm64 -I/opt/homebrew/include"
        export LDFLAGS="-L/opt/homebrew/lib"
        # Add OpenMP support
        export CPPFLAGS="-Xpreprocessor -fopenmp"
        export CFLAGS="$CFLAGS -I/opt/homebrew/opt/libomp/include"
        export CXXFLAGS="$CXXFLAGS -I/opt/homebrew/opt/libomp/include"
        export LDFLAGS="$LDFLAGS -L/opt/homebrew/opt/libomp/lib -lomp"
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        # Linux flags
        export CFLAGS="-O3"
        export CXXFLAGS="-O3"
    fi
}

# Function to install packages using uv or pip with platform-specific options
install_packages() {
    local packages=("$@")
    local install_cmd="pip install"
    local arch=$(get_cpu_arch)
    
    # Set platform-specific compiler flags
    set_compiler_flags
    
    if command_exists uv; then
        install_cmd="uv pip install"
    fi
    
    # Platform-specific installation handling
    if [[ "$OSTYPE" == "darwin"* && "$arch" == "arm64" ]]; then
        # Try binary wheel first
        if ! $install_cmd --only-binary :all: "${packages[@]}" 2>/dev/null; then
            # If binary fails, try source with optimized flags
            echo "Binary wheel not available, attempting source build..."
            $install_cmd --no-binary :all: "${packages[@]}" || return 1
        fi
    else
        # Default installation
        $install_cmd "${packages[@]}" || return 1
    fi
}

# Function to build and install package from source with enhanced error handling
build_from_source() {
    local package_name=$1
    local git_url=$2
    local branch=${3:-main}
    local fallback_url=${4:-""}
    local requirements=${5:-""}

    echo "Building $package_name from source..."
    temp_dir=$(mktemp -d)
    
    # Set compiler flags before building
    set_compiler_flags
    
    if git clone "$git_url" "$temp_dir"; then
        cd "$temp_dir"
        if git checkout "$branch"; then
            if [ -n "$requirements" ]; then
                if ! pip install -r "$requirements"; then
                    echo "Failed to install requirements for $package_name"
                    cd - >/dev/null
                    rm -rf "$temp_dir"
                    return 1
                fi
            fi
            
            # Try to build with optimized flags first
            if ! pip install -e .; then
                echo "Failed to install $package_name from source"
                if [ -n "$fallback_url" ]; then
                    echo "Trying fallback installation for $package_name..."
                    if ! pip install "$fallback_url"; then
                        echo "Fallback installation failed for $package_name"
                        cd - >/dev/null
                        rm -rf "$temp_dir"
                        return 1
                    fi
                else
                    cd - >/dev/null
                    rm -rf "$temp_dir"
                    return 1
                fi
            fi
        else
            echo "Failed to checkout branch $branch for $package_name"
            cd - >/dev/null
            rm -rf "$temp_dir"
            return 1
        fi
        cd - >/dev/null
    else
        echo "Failed to clone $package_name repository"
        if [ -n "$fallback_url" ]; then
            echo "Trying fallback installation for $package_name..."
            if ! pip install "$fallback_url"; then
                echo "Fallback installation failed for $package_name"
                rm -rf "$temp_dir"
                return 1
            fi
        else
            rm -rf "$temp_dir"
            return 1
        fi
    fi
    rm -rf "$temp_dir"
    return 0
}

# Function to get Python version
get_python_version() {
    python -c 'import sys; v = sys.version_info; print(f"{v.major}.{v.minor}")'
}

# Function to check if Python version is >= 3.11
needs_legacy_handling() {
    python -c 'import sys; v = sys.version_info; exit(0 if v.major == 3 and v.minor >= 11 or v.major > 3 else 1)'
}

# Check if requirements.txt exists
if [ ! -f "requirements.txt" ]; then
    echo "Error: requirements.txt not found"
    exit 1
fi

# Check if conda is available and use it if present
if command_exists conda; then
    echo "Conda found. Checking environment..."
    if [ -z "$CONDA_DEFAULT_ENV" ] || [ "$CONDA_DEFAULT_ENV" != "chemdata_311" ]; then
        if conda env list | grep -q "chemdata_311"; then
            echo "Environment 'chemdata_311' exists. Activating..."
        else
            echo "Creating environment 'chemdata_311' with Python 3.11..."
            conda create -y -n chemdata_311 python=3.11
        fi
        eval "$(conda shell.bash hook)"
        conda activate chemdata_311
    else
        echo "Using conda environment: chemdata_311"
    fi
else
    echo "Conda not found. Setting up with venv..."
    # Create virtual environment if it doesn't exist
    if [ ! -d "venv" ]; then
        echo "Creating virtual environment..."
        python3 -m venv venv
    fi

    # Activate virtual environment
    source venv/bin/activate

    # Check if uv is available, install if not
    if ! command_exists uv; then
        echo "uv not found. Installing uv..."
        curl -LsSf https://astral.sh/uv/install.sh | sh
    fi
fi

# Create necessary directories
echo "Creating project directories..."
mkdir -p data
mkdir -p cache
mkdir -p output
mkdir -p logs

# Install platform-specific packages
echo "Installing platform-specific packages..."
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    echo "Detected macOS..."
    if ! command_exists brew; then
        echo "Homebrew not found. Please install Homebrew first."
        exit 1
    fi
    
    # Install required build tools and libraries
    echo "Installing build dependencies..."
    brew install cmake open-babel libomp
    
    # Additional Apple Silicon specific setup
    if [[ "$(get_cpu_arch)" == "arm64" ]]; then
        echo "Configuring for Apple Silicon..."
        # Set platform-specific compiler flags
        set_compiler_flags
        
        # Install Rosetta 2 if needed (for x86_64 binaries)
        if ! pkgutil --pkg-info com.apple.pkg.RosettaUpdateAuto >/dev/null 2>&1; then
            echo "Installing Rosetta 2..."
            softwareupdate --install-rosetta --agree-to-license
        fi
    fi
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    # Linux
    echo "Detected Linux..."
    sudo apt-get update
    sudo apt-get install -y \
        build-essential \
        cmake \
        libopenbabel-dev \
        libomp-dev \
        python3-dev \
        libffi-dev
fi

# Update basic tools
echo "Updating basic tools..."
install_packages pip setuptools wheel

# Install core dependencies using conda first
echo "Installing core dependencies..."
if command_exists conda; then
    echo "Installing core packages with conda..."
    conda install -y -c conda-forge numpy pandas scipy pyarrow biopython
else
    echo "Installing core packages with pip..."
    install_packages numpy pandas scipy pyarrow
fi

# Install chemistry tools
echo "Installing chemistry tools..."
if command_exists conda; then
    echo "Installing chemistry packages with conda..."
    # Install core chemistry packages from conda-forge
    conda install -y -c conda-forge \
        rdkit \
        openbabel \
        mordred \
        molvs

    # Verify and install specific versions if needed
    echo "Verifying package versions..."
    if ! check_package_version "mordred" "1.2.0"; then
        echo "Installing specific version of mordred..."
        install_packages "mordred==1.2.0"
    fi
    if ! check_package_version "molvs" "0.1.1"; then
        echo "Installing specific version of molvs..."
        install_packages "molvs==0.1.1"
    fi
else
    echo "Installing chemistry packages with pip..."
    install_packages "rdkit>=2023.3.1" "mordred==1.2.0" "molvs==0.1.1"
fi

# Function to handle package-specific installation requirements
get_package_info() {
    local package_name=$1
    local info=()
    
    case "$package_name" in
        "pytorch")
            if [[ "$OSTYPE" == "darwin"* ]]; then
                info+=("pytorch/pytorch" "pytorch" "torch" "main")
            else
                info+=("pytorch/pytorch" "pytorch,nvidia" "torch" "main")
            fi
            ;;
        "torchvision")
            info+=("pytorch/vision" "pytorch" "torchvision" "main")
            ;;
        "torchaudio")
            info+=("pytorch/audio" "pytorch" "torchaudio" "main")
            ;;
        "pytorch-lightning")
            info+=("Lightning-AI/lightning" "conda-forge" "pytorch-lightning" "main")
            ;;
        "scikit-learn")
            info+=("scikit-learn/scikit-learn" "conda-forge" "scikit-learn" "main")
            ;;
        "xgboost")
            info+=("dmlc/xgboost" "conda-forge" "xgboost" "master")
            ;;
        "lightgbm")
            info+=("microsoft/LightGBM" "conda-forge" "lightgbm" "master")
            ;;
        "optuna")
            info+=("optuna/optuna" "conda-forge" "optuna" "master")
            ;;
        "deepchem")
            info+=("deepchem/deepchem" "conda-forge" "deepchem" "master")
            ;;
        "transformers")
            info+=("huggingface/transformers" "conda-forge" "transformers" "main")
            ;;
        *)
            info+=("" "conda-forge" "$package_name" "main")
            ;;
    esac
    
    echo "${info[@]}"
}

# Function to try installing packages through different methods
try_install_package() {
    local package_name=$1
    local base_name="${package_name%%[<>=]*}"
    
    # Skip if package is already installed and version check passes
    if package_installed "$base_name"; then
        if [[ "$package_name" == *"<"* || "$package_name" == *">"* || "$package_name" == *"="* ]]; then
            version_req="${package_name#$base_name}"
            if check_package_version "$base_name" "${version_req#[<>=]}"; then
                echo "Package $package_name is already installed with compatible version"
                return 0
            fi
        else
            echo "Package $base_name is already installed"
            return 0
        fi
    fi
    read -r repo channels pkg branch <<< "$(get_package_info "$package_name")"
    
    echo "Attempting to install $package_name..."
    
    # Try conda first if available
    if command_exists conda; then
        echo "Trying conda install..."
        
        # Handle package-specific dependencies
        case "$package_name" in
            "pytorch"|"torchvision"|"torchaudio")
                if [[ "$OSTYPE" == "darwin"* ]]; then
                    # macOS - no CUDA needed
                    if conda install -y -c pytorch "$pkg" 2>/dev/null; then
                        return 0
                    fi
                else
                    # Linux - try with CUDA
                    if conda install -y -c pytorch -c nvidia "$pkg" pytorch-cuda=11.8 2>/dev/null; then
                        return 0
                    fi
                fi
                ;;
            "deepchem")
                # Install with RDKit dependency
                if conda install -y -c conda-forge deepchem rdkit 2>/dev/null; then
                    return 0
                fi
                ;;
            *)
                # Try conda-forge first
                if conda install -y -c conda-forge "$pkg" 2>/dev/null; then
                    return 0
                fi
                
                # Try specified channels
                IFS=',' read -ra channel_array <<< "$channels"
                for channel in "${channel_array[@]}"; do
                    if conda install -y -c "$channel" "$pkg" 2>/dev/null; then
                        return 0
                    fi
                done
                ;;
        esac
    fi
    
    # Try uv/pip with any package-specific options
    echo "Trying pip install..."
    case "$package_name" in
        "pytorch"|"torchvision"|"torchaudio")
            if [[ "$OSTYPE" == "darwin"* ]]; then
                # macOS - no CUDA
                if install_packages "$pkg"; then
                    return 0
                fi
            else
                # Linux - try with CUDA
                if install_packages "$pkg[cuda]"; then
                    return 0
                fi
            fi
            ;;
        *)
            if install_packages "$pkg"; then
                return 0
            fi
            ;;
    esac
    
    # Try building from source if repo provided
    if [ -n "$repo" ]; then
        echo "Trying to build from source..."
        local git_url="https://github.com/$repo.git"
        
        # Handle package-specific build requirements
        case "$package_name" in
            "pytorch")
                # PyTorch needs special build setup
                if build_from_source "$pkg" "$git_url" "$branch" "" "requirements.txt"; then
                    return 0
                fi
                ;;
            *)
                if build_from_source "$pkg" "$git_url" "$branch"; then
                    return 0
                fi
                ;;
        esac
    fi
    
    echo "Failed to install $package_name through any method"
    return 1
}

# Install special dependencies with platform-specific handling
echo "Installing special dependencies..."
if ! bash scripts/install_special_deps.sh; then
    echo "Error: Failed to install special dependencies"
    exit 1
fi

# Install machine learning packages
echo "Installing machine learning packages..."

# Install PyTorch ecosystem
echo "Installing PyTorch ecosystem..."
try_install_package "pytorch"
try_install_package "torchvision"
try_install_package "torchaudio"
try_install_package "pytorch-lightning"

# Install other ML packages
echo "Installing other ML packages..."
try_install_package "scikit-learn"
try_install_package "xgboost"
try_install_package "lightgbm"
try_install_package "optuna"
try_install_package "deepchem"

# Install TensorFlow dependencies in the correct order
echo "Installing TensorFlow dependencies..."
install_packages "tensorflow"
install_packages "tf-keras<3.0.0"  # Use older version for compatibility
try_install_package "transformers"
try_install_package "sentence_transformers"

# Install web and API packages
echo "Installing web and API packages..."
install_packages requests aiohttp beautifulsoup4 lxml selenium webdriver_manager chembl_webresource_client
install_packages praw tweepy crawl4ai playwright tf-playwright-stealth

# Install Playwright browsers
echo "Installing Playwright browsers..."
if command_exists playwright; then
    playwright install
else
    echo "Warning: playwright not found after installation. Browser automation may not work correctly."
fi

# Handle scholarly installation
echo "Installing scholarly..."
if needs_legacy_handling; then
    build_from_source "scholarly" "https://github.com/scholarly-python-package/scholarly.git" "main" "scholarly>=1.7.0"
else
    install_packages scholarly
fi

install_packages atproto patent-client

# Install web interface packages
echo "Installing web interface packages..."
install_packages dash dash-bio dash-bootstrap-components plotly flask flask-restful flask-sqlalchemy
install_packages flask-cors flask-caching fastapi uvicorn gunicorn starlette

# Install pydantic with specific version for compatibility
echo "Installing pydantic with specific version..."
install_packages "pydantic>=2.7.1,<3.0.0"

# Install NLP packages
echo "Installing NLP packages..."
if command_exists conda; then
    echo "Installing NLP packages with conda-forge..."
    conda install -y -c conda-forge nltk gensim
else
    echo "Installing NLP packages with pip..."
    install_packages nltk gensim
fi

install_packages conllu pysbd

# Install database packages
echo "Installing database packages..."
install_packages sqlalchemy psycopg2-binary alembic

# Install development tools
echo "Installing development tools..."
if command_exists conda; then
    echo "Installing development tools with conda-forge..."
    conda install -y -c conda-forge black isort flake8 mypy bandit jupyterlab ipython
else
    echo "Installing development tools with pip..."
    install_packages black isort flake8 mypy bandit jupyterlab ipython
fi

# Install pre-commit with --user to ensure it's available globally
echo "Installing pre-commit..."
pip install --user pre-commit

# Initialize pre-commit hooks with error handling
echo "Setting up pre-commit hooks..."
if command_exists pre-commit; then
    pre-commit install || echo "Warning: Failed to initialize pre-commit hooks"
else
    echo "Warning: pre-commit not found after installation. Git hooks will not be set up."
fi

# Install testing packages
echo "Installing testing packages..."
if command_exists conda; then
    echo "Installing testing packages with conda-forge..."
    conda install -y -c conda-forge pytest pytest-cov pytest-asyncio pytest-mock responses factory-boy
else
    echo "Installing testing packages with pip..."
    install_packages pytest pytest-cov pytest-asyncio pytest-mock responses factory-boy
fi

# Install documentation packages
echo "Installing documentation packages..."
install_packages sphinx sphinx-rtd-theme sphinx-autodoc-typehints nbsphinx jupyter myst-parser ipykernel

# Handle epo-ops-client installation
echo "Installing epo-ops-client..."
pip install python-epo-ops-client==4.0.0 || pip install python-epo-ops-client==3.1.4

# Install remaining dependencies (excluding already installed packages)
echo "Installing remaining dependencies..."
while read -r package; do
    if ! echo "$package" | grep -qE "^#|^$|epo-ops-client|chemdataextractor|scispacy|en-core-sci-lg|rdkit|torch-geometric|scholarly"; then
        install_packages "$package"
    fi
done < requirements.txt

# Create example configuration file if it doesn't exist
if [ ! -f ".env" ]; then
    echo "Creating example .env file..."
    cat > .env << EOL
# API Keys
REDDIT_CLIENT_ID=your_client_id
REDDIT_CLIENT_SECRET=your_client_secret
TWITTER_API_KEY=your_api_key
TWITTER_API_SECRET=your_api_secret
ESPACENET_KEY=your_key
ESPACENET_SECRET=your_secret

# Database Configuration
DB_HOST=localhost
DB_PORT=5432
DB_NAME=chemdata
DB_USER=your_username
DB_PASSWORD=your_password

# Cache Configuration
CACHE_DIR=cache
CACHE_TIMEOUT=3600

# Web Configuration
WEB_HOST=0.0.0.0
WEB_PORT=8000
DEBUG=true
EOL
    echo "Created example .env file. Please update with your actual credentials."
fi

# Initialize git hooks
if [ -d ".git" ]; then
    echo "Setting up git hooks..."
    cp scripts/pre-commit .git/hooks/
    chmod +x .git/hooks/pre-commit
fi

# Install the package in development mode
echo "Installing binding_data_processor in development mode..."
pip install -e .

echo "Environment setup complete!"
echo "Next steps:"
echo "1. Update .env file with your credentials"
echo "2. Run './scripts/export_all_compounds.sh' to start the export process"

# Run export scripts with error handling
echo "Running export scripts..."
if [ -f "scripts/export_all_compounds.sh" ]; then
    if bash scripts/export_all_compounds.sh; then
        echo "Successfully ran export_all_compounds.sh"
    else
        echo "Warning: export_all_compounds.sh failed but continuing..."
    fi
else
    echo "Warning: export_all_compounds.sh not found, skipping..."
fi

# Run export_compounds.py last
if [ -f "scripts/export_compounds.py" ]; then
    if python scripts/export_compounds.py; then
        echo "Successfully ran export_compounds.py"
    else
        echo "Warning: export_compounds.py failed"
    fi
else
    echo "Warning: export_compounds.py not found, skipping..."
fi
