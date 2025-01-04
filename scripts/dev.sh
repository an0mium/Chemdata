#!/bin/bash
set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

# Function to show usage
show_help() {
    echo "Development utility script"
    echo
    echo "Usage: $0 [command] [options]"
    echo
    echo "Commands:"
    echo "  test        Run tests"
    echo "  lint        Run linting"
    echo "  format      Format code"
    echo "  check       Run all checks (test + lint)"
    echo "  clean       Clean build artifacts"
    echo "  deps        Update dependencies"
    echo "  serve       Start development server"
    echo "  build       Build package"
    echo
    echo "Options:"
    echo "  --coverage  Generate coverage report (with test)"
    echo "  --watch     Watch for changes (with test/serve)"
    echo "  --all       Clean everything (with clean)"
    echo "  --check     Check only, don't format (with format)"
    echo
    echo "Examples:"
    echo "  $0 test --coverage     Run tests with coverage"
    echo "  $0 serve --watch       Start dev server with auto-reload"
    echo "  $0 format             Format code"
}

# Function to run tests
run_tests() {
    echo -e "${BLUE}Running tests...${NC}"
    
    # Base command
    cmd="python -m pytest"
    
    # Add coverage if requested
    if [[ "$*" == *"--coverage"* ]]; then
        cmd="$cmd --cov=binding_data_processor --cov-report=html"
    fi
    
    # Add watch if requested
    if [[ "$*" == *"--watch"* ]]; then
        cmd="$cmd --watch"
    fi
    
    # Run tests
    $cmd
    
    # Show coverage report if generated
    if [[ "$*" == *"--coverage"* ]]; then
        echo -e "${BLUE}Coverage report generated in htmlcov/index.html${NC}"
    fi
}

# Function to run linting
run_lint() {
    echo -e "${BLUE}Running linting...${NC}"
    
    echo "Running black..."
    python -m black . --check
    
    echo "Running isort..."
    python -m isort . --check-only
    
    echo "Running flake8..."
    python -m flake8 .
    
    echo "Running mypy..."
    python -m mypy .
}

# Function to format code
format_code() {
    if [[ "$*" == *"--check"* ]]; then
        echo -e "${BLUE}Checking code format...${NC}"
        python -m black . --check
        python -m isort . --check-only
    else
        echo -e "${BLUE}Formatting code...${NC}"
        python -m black .
        python -m isort .
    fi
}

# Function to clean build artifacts
clean_build() {
    echo -e "${BLUE}Cleaning build artifacts...${NC}"
    
    # Always clean these
    rm -rf build/ dist/ *.egg-info/
    find . -type d -name __pycache__ -exec rm -rf {} +
    find . -type f -name "*.pyc" -delete
    find . -type f -name "*.pyo" -delete
    find . -type f -name "*.pyd" -delete
    find . -type f -name ".coverage" -delete
    find . -type d -name "htmlcov" -exec rm -rf {} +
    find . -type f -name ".coverage.*" -delete
    
    # Clean everything if requested
    if [[ "$*" == *"--all"* ]]; then
        echo -e "${BLUE}Cleaning everything...${NC}"
        rm -rf .venv/
        rm -rf .pytest_cache/
        rm -rf .mypy_cache/
        rm -rf .tox/
    fi
}

# Function to update dependencies
update_deps() {
    echo -e "${BLUE}Updating dependencies...${NC}"
    
    # Update core dependencies
    uv pip install -U -e .
    
    # Update development dependencies
    uv pip install -U -e ".[dev]"
    
    # Update pre-commit hooks
    pre-commit autoupdate
}

# Function to start development server
start_server() {
    echo -e "${BLUE}Starting development server...${NC}"
    
    if [[ "$*" == *"--watch"* ]]; then
        # Use uvicorn with auto-reload
        python -m uvicorn web.app:app --reload --port 8050
    else
        # Start normally
        python -m web.app
    fi
}

# Function to build package
build_package() {
    echo -e "${BLUE}Building package...${NC}"
    
    # Clean first
    clean_build
    
    # Build package
    python -m build
    
    echo -e "${GREEN}Build complete! Artifacts in dist/${NC}"
}

# Main command processing
case "$1" in
    test)
        shift
        run_tests "$@"
        ;;
    lint)
        run_lint
        ;;
    format)
        shift
        format_code "$@"
        ;;
    check)
        run_tests
        run_lint
        ;;
    clean)
        shift
        clean_build "$@"
        ;;
    deps)
        update_deps
        ;;
    serve)
        shift
        start_server "$@"
        ;;
    build)
        build_package
        ;;
    help|-h|--help)
        show_help
        ;;
    *)
        echo -e "${RED}Unknown command: $1${NC}"
        show_help
        exit 1
        ;;
esac
