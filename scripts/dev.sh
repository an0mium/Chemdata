#!/bin/bash
set -e

# This script sets up and runs the development environment using uv

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

# Install dependencies if needed (using uv's fast dependency resolution)
if [ ! -d ".venv/lib/python3.12/site-packages/chemdata" ]; then
    echo "Installing dependencies..."
    uv pip install -e ".[dev]"
fi

# Run pre-commit install if not already installed
if [ ! -f ".git/hooks/pre-commit" ]; then
    echo "Setting up pre-commit hooks..."
    pre-commit install
fi

# Start development processes in parallel
echo "Starting development environment..."

# Function to cleanup background processes on exit
cleanup() {
    echo "Cleaning up..."
    kill $(jobs -p) 2>/dev/null
}
trap cleanup EXIT

# Start Flask development server
echo "Starting web server..."
python -m chemdata.web &

# Start file watcher for auto-reloading
echo "Starting file watcher..."
watchmedo auto-restart \
    --directory=./chemdata \
    --pattern="*.py" \
    --recursive \
    python -m chemdata.web &

# Show logs
echo "Development environment is running!"
echo "Web interface: http://localhost:5000"
echo "Press Ctrl+C to stop"

# Wait for any process to exit
wait -n

# Exit with status of process that exited first
exit $?
