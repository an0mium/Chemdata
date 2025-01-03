#!/bin/bash
set -e

# This script runs the complete data collection pipeline

# Ensure we're in the project root
cd "$(dirname "$0")/.."

# Function to check required environment variables
check_env_vars() {
    local missing=()
    if [ -z "$OPENAI_API_KEY" ]; then missing+=("OPENAI_API_KEY"); fi
    if [ -z "$SERP_API_KEY" ]; then missing+=("SERP_API_KEY"); fi
    
    if [ ${#missing[@]} -ne 0 ]; then
        echo "❌ Missing required environment variables:"
        printf '%s\n' "${missing[@]}"
        exit 1
    fi
}

# Function to show step header
show_step() {
    echo -e "\n=== $1 ===\n"
}

# Check environment variables
show_step "Checking environment variables"
check_env_vars

# Check dependencies
show_step "Checking dependencies"
./scripts/check_deps.sh || {
    echo "❌ Dependency check failed"
    exit 1
}

# Download BindingDB data
show_step "Downloading BindingDB data"
./scripts/get_bindingdb.sh || {
    echo "❌ Failed to download BindingDB data"
    exit 1
}

# Collect 5-HT2 receptor ligands
show_step "Collecting 5-HT2 receptor ligands"
python3 -c "
from binding_data_processor import BindingDataProcessor
from api_client import PubMedClient
processor = BindingDataProcessor(PubMedClient())
processor.save_5ht2_ligands(
    'data/5ht2_ligands.tsv',
    llm_api_key='$OPENAI_API_KEY'
)" || {
    echo "❌ Failed to collect 5-HT2 receptor ligands"
    exit 1
}

# Process and enrich the data
show_step "Processing and enriching data"
python3 -c "
from binding_data_processor import BindingDataProcessor
from api_client import PubMedClient
processor = BindingDataProcessor(PubMedClient())
processor.process_file(
    'data/5ht2_ligands.tsv',
    'data/5ht2_ligands_enriched.tsv'
)" || {
    echo "❌ Failed to process and enrich data"
    exit 1
}

# Start web interface
show_step "Starting web interface"
echo "Starting web server at http://localhost:5000"
python3 -m chemdata.web

echo -e "\n✅ Data collection pipeline complete!"
echo "Enriched data available in: data/5ht2_ligands_enriched.tsv"
echo "Web interface running at: http://localhost:5000"
