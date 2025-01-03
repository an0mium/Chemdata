#!/bin/bash
set -e

# This script checks if all required dependencies are properly installed

# Ensure we're in the project root
cd "$(dirname "$0")/.."

# Check if uv is installed
if ! command -v uv &> /dev/null; then
    echo "❌ uv is not installed"
    echo "Install with: curl -LsSf https://astral.sh/uv/install.sh | sh"
    exit 1
fi

# Check virtual environment
if [ ! -d ".venv" ]; then
    echo "❌ Virtual environment not found"
    echo "Create with: uv venv -p python3.12 .venv"
    exit 1
fi

# Activate virtual environment
source .venv/bin/activate

# Function to check Python package
check_package() {
    local pkg=$1
    local min_version=$2
    
    if python3 -c "import $pkg" 2>/dev/null; then
        version=$(python3 -c "import $pkg; print($pkg.__version__ if hasattr($pkg, '__version__') else 'unknown')" 2>/dev/null)
        if [ -n "$min_version" ]; then
            if python3 -c "from packaging import version; exit(0 if version.parse('$version') >= version.parse('$min_version') else 1)" 2>/dev/null; then
                echo "✅ $pkg $version (>= $min_version)"
            else
                echo "⚠️  $pkg $version (< $min_version required)"
            fi
        else
            echo "✅ $pkg $version"
        fi
    else
        echo "❌ $pkg not found"
    fi
}

echo "Checking Python version..."
python_version=$(python3 --version | cut -d' ' -f2)
if python3 -c "from packaging import version; exit(0 if version.parse('$python_version') >= version.parse('3.12') else 1)"; then
    echo "✅ Python $python_version"
else
    echo "❌ Python $python_version (>= 3.12 required)"
    exit 1
fi

echo -e "\nChecking core dependencies..."
check_package "rdkit" "2024.3.1"
check_package "pandas" "2.2.0"
check_package "flask" "2.0.0"
check_package "requests" "2.31.0"
check_package "bs4" "4.12.0"  # beautifulsoup4
check_package "pubchempy" "1.0.4"
check_package "tqdm" "4.66.0"
check_package "pypdf" "3.0.0"

echo -e "\nChecking web dependencies..."
check_package "aiohttp" "3.8.0"
check_package "jinja2" "3.0.0"
check_package "werkzeug" "2.0.0"

echo -e "\nChecking development dependencies..."
check_package "pytest" "7.0.0"
check_package "black" "23.0.0"
check_package "isort" "5.0.0"
check_package "mypy" "1.0.0"
check_package "flake8" "6.0.0"

echo -e "\nChecking optional dependencies..."
check_package "ipywidgets" "8.1.0"
check_package "openpyxl" "3.1.2"

echo -e "\nChecking dependency conflicts..."
uv pip check

echo -e "\nChecking for updates..."
uv pip list --outdated

echo -e "\nChecking virtual environment..."
if [ -f ".venv/pyvenv.cfg" ]; then
    echo "✅ Virtual environment configuration found"
    cat .venv/pyvenv.cfg
else
    echo "❌ Virtual environment configuration not found"
fi

echo -e "\nChecking pre-commit hooks..."
if [ -f ".git/hooks/pre-commit" ]; then
    echo "✅ Pre-commit hooks installed"
else
    echo "⚠️  Pre-commit hooks not installed"
    echo "Install with: pre-commit install"
fi

echo -e "\nAll checks complete!"
