#!/bin/bash
set -e

# This script updates dependencies using uv's dependency resolution

# Ensure we're in the project root
cd "$(dirname "$0")/.."

# Check if uv is installed
if ! command -v uv &> /dev/null; then
    echo "uv is not installed. Installing now..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
fi

# Activate virtual environment if it exists, create it if it doesn't
if [ ! -d ".venv" ]; then
    echo "Creating virtual environment..."
    uv venv -p python3.12 .venv
fi
source .venv/bin/activate

# Update all dependencies to their latest compatible versions
echo "Updating dependencies..."
uv pip compile pyproject.toml --upgrade-all -o requirements.txt

# Install updated dependencies
echo "Installing updated dependencies..."
uv pip install -e ".[dev]"

# Show outdated packages
echo "Checking for outdated packages..."
uv pip list --outdated

# Run tests to ensure everything still works
echo "Running tests..."
pytest

echo "Dependency update complete!"
echo "Review the changes in requirements.txt and commit if everything looks good."
