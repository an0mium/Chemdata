#!/bin/bash

# Function to show usage
show_usage() {
    echo "Usage: source ./scripts/manage_venv.sh [3.11|3.12]"
    echo "  3.11 - Activate Python 3.11 environment (for most packages)"
    echo "  3.12 - Activate Python 3.12 environment (for PyTorch and ML packages)"
    echo ""
    echo "Example:"
    echo "  source ./scripts/manage_venv.sh 3.11  # For general development"
    echo "  source ./scripts/manage_venv.sh 3.12  # For ML/PyTorch work"
}

# Check if the script is being sourced
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "Error: This script must be sourced"
    echo "Usage: source ./scripts/manage_venv.sh [3.11|3.12]"
    exit 1
fi

# Check for argument
if [ -z "$1" ]; then
    show_usage
    return 1
fi

# Function to deactivate current environment if one is active
deactivate_current() {
    if [ ! -z "$VIRTUAL_ENV" ]; then
        deactivate
    fi
}

# Activate the appropriate virtual environment
case "$1" in
    "3.11")
        if [ -d "venv-3.11" ]; then
            deactivate_current
            source venv-3.11/bin/activate
            echo "Activated Python 3.11 environment"
            echo "This environment contains:"
            echo "- Web frameworks (Flask, FastAPI)"
            echo "- Database tools (SQLAlchemy)"
            echo "- Development tools (black, flake8)"
            echo "- Documentation tools (Sphinx)"
            echo "- Testing tools (pytest)"
        else
            echo "Error: Python 3.11 environment not found"
            echo "Run ./scripts/setup_environment.sh first"
            return 1
        fi
        ;;
    "3.12")
        if [ -d "venv-3.12" ]; then
            deactivate_current
            source venv-3.12/bin/activate
            echo "Activated Python 3.12 environment"
            echo "This environment contains:"
            echo "- PyTorch and related packages"
            echo "- Deep learning frameworks"
            echo "- ML tools (transformers, pytorch-lightning)"
        else
            echo "Error: Python 3.12 environment not found"
            echo "Run ./scripts/setup_environment.sh first"
            return 1
        fi
        ;;
    *)
        show_usage
        return 1
        ;;
esac

# Print current Python version and environment
echo ""
echo "Current environment:"
echo "Python: $(python --version)"
echo "Environment: $VIRTUAL_ENV"
