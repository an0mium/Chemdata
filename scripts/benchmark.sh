#!/bin/bash
set -e

# This script benchmarks uv's performance against pip

# Ensure we're in the project root
cd "$(dirname "$0")/.."

# Function to time command execution
time_cmd() {
    local start=$(date +%s.%N)
    "$@"
    local end=$(date +%s.%N)
    echo "$(echo "$end - $start" | bc) seconds"
}

# Function to clean environment
clean_env() {
    rm -rf .venv
    rm -rf .uv
    rm -rf .uv-cache
}

echo "Benchmarking dependency installation performance..."
echo "------------------------------------------------"

# Test uv
echo "Testing uv..."
clean_env

echo "1. Creating virtual environment with uv..."
time_cmd uv venv -p python3.12 .venv
source .venv/bin/activate

echo "2. Installing dependencies with uv..."
time_cmd uv pip install -e ".[dev]"

echo "3. Installing dependencies again (cached) with uv..."
time_cmd uv pip install -e ".[dev]"

deactivate
clean_env

# Test pip
echo -e "\nTesting pip..."
clean_env

echo "1. Creating virtual environment with venv..."
time_cmd python3.12 -m venv .venv
source .venv/bin/activate

echo "2. Installing dependencies with pip..."
time_cmd pip install -e ".[dev]"

echo "3. Installing dependencies again (cached) with pip..."
time_cmd pip install -e ".[dev]"

deactivate
clean_env

echo -e "\nBenchmarking dependency resolution performance..."
echo "------------------------------------------------"

# Test uv pip compile
echo "Testing uv pip compile..."
time_cmd uv pip compile pyproject.toml -o requirements.txt

# Test pip-compile
echo -e "\nTesting pip-compile..."
time_cmd pip-compile pyproject.toml -o requirements.txt

echo -e "\nBenchmarking complete!"
