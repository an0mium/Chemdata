#!/bin/bash
set -e

# This script demonstrates how to use the ChemData package with uv

# Ensure we're in the project root
cd "$(dirname "$0")/.."

# Check if uv is installed
if ! command -v uv &> /dev/null; then
    echo "uv is not installed. Installing now..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
fi

# Create a new virtual environment if it doesn't exist
if [ ! -d ".venv" ]; then
    echo "Creating virtual environment..."
    uv venv -p python3.12 .venv
fi

# Activate virtual environment
source .venv/bin/activate

# Install dependencies if needed
if [ ! -d ".venv/lib/python3.12/site-packages/chemdata" ]; then
    echo "Installing dependencies..."
    uv pip install -e ".[dev]"
fi

# Example 1: Search for 5-HT2A receptor ligands
echo "Example 1: Searching for 5-HT2A receptor ligands..."
python -m chemdata.cli --target "5-HT2A" --output data/5ht2a_ligands.tsv

# Example 2: Process a list of compounds
echo "Example 2: Processing compound list..."
echo -e "LSD\nDMT\nPsilocin" > data/example_compounds.txt
python -m chemdata.cli --compounds data/example_compounds.txt --output data/example_results.tsv

# Example 3: Start the web interface
echo "Example 3: Starting web interface..."
echo "Visit http://localhost:5000 in your browser"
python -m chemdata.web

# Note: The web interface will keep running until you press Ctrl+C
