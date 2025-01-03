#!/bin/bash
set -e

# This script runs tests using uv's fast dependency resolution and test caching

# Ensure we're in the project root
cd "$(dirname "$0")/.."

# Check if uv is installed
if ! command -v uv &> /dev/null; then
    echo "uv is not installed. Installing now..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
fi

# Create and activate virtual environment if needed
if [ ! -d ".venv" ]; then
    echo "Creating virtual environment..."
    uv venv -p python3.12 .venv
fi
source .venv/bin/activate

# Install test dependencies if needed
if [ ! -d ".venv/lib/python3.12/site-packages/pytest" ]; then
    echo "Installing test dependencies..."
    uv pip install -e ".[dev]"
fi

# Parse command line arguments
PYTEST_ARGS=""
COVERAGE=false
WATCH=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --coverage)
            COVERAGE=true
            shift
            ;;
        --watch)
            WATCH=true
            shift
            ;;
        *)
            PYTEST_ARGS="$PYTEST_ARGS $1"
            shift
            ;;
    esac
done

# Function to run tests
run_tests() {
    if [ "$COVERAGE" = true ]; then
        pytest --cov=chemdata --cov-report=html --cov-report=term $PYTEST_ARGS
    else
        pytest $PYTEST_ARGS
    fi
}

if [ "$WATCH" = true ]; then
    # Watch for changes and run tests
    echo "Watching for changes..."
    watchmedo shell-command \
        --patterns="*.py" \
        --recursive \
        --command="clear; echo 'Running tests...'; $(declare -f run_tests); run_tests" \
        .
else
    # Run tests once
    run_tests
fi

# If coverage was generated, show the report location
if [ "$COVERAGE" = true ]; then
    echo "Coverage report generated in htmlcov/index.html"
fi
