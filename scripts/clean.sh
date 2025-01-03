#!/bin/bash
set -e

# This script cleans up uv-related files and caches

# Ensure we're in the project root
cd "$(dirname "$0")/.."

# Function to confirm action
confirm() {
    read -r -p "${1:-Are you sure?} [y/N] " response
    case "$response" in
        [yY][eE][sS]|[yY]) 
            true
            ;;
        *)
            false
            ;;
    esac
}

# Parse command line arguments
ALL=false
VENV=false
CACHE=false
BUILD=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --all)
            ALL=true
            shift
            ;;
        --venv)
            VENV=true
            shift
            ;;
        --cache)
            CACHE=true
            shift
            ;;
        --build)
            BUILD=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# If no specific options, ask what to clean
if [ "$ALL" = false ] && [ "$VENV" = false ] && [ "$CACHE" = false ] && [ "$BUILD" = false ]; then
    echo "What would you like to clean?"
    echo "1) Everything (--all)"
    echo "2) Virtual environment only (--venv)"
    echo "3) Cache only (--cache)"
    echo "4) Build artifacts only (--build)"
    read -r -p "Enter option number [1-4]: " option
    case "$option" in
        1) ALL=true ;;
        2) VENV=true ;;
        3) CACHE=true ;;
        4) BUILD=true ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
fi

# Clean virtual environment
if [ "$ALL" = true ] || [ "$VENV" = true ]; then
    if [ ! -d ".venv" ] || confirm "Remove virtual environment?"; then
        echo "Removing virtual environment..."
        rm -rf .venv
    fi
fi

# Clean uv cache
if [ "$ALL" = true ] || [ "$CACHE" = true ]; then
    if [ ! -d ".uv" ] || confirm "Clear uv cache?"; then
        echo "Clearing uv cache..."
        rm -rf .uv
        rm -rf .uv-cache
    fi
fi

# Clean build artifacts
if [ "$ALL" = true ] || [ "$BUILD" = true ]; then
    if confirm "Remove build artifacts?"; then
        echo "Removing build artifacts..."
        rm -rf build/
        rm -rf dist/
        rm -rf *.egg-info/
        rm -rf .eggs/
        rm -rf .pytest_cache/
        rm -rf .mypy_cache/
        rm -rf .coverage
        rm -rf htmlcov/
        find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
    fi
fi

echo "Clean complete!"
