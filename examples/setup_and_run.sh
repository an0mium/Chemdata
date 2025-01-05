#!/bin/bash
set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Print step header
print_step() {
    echo -e "\n${YELLOW}=== $1 ===${NC}\n"
}

# Check if command exists
check_command() {
    if ! command -v $1 &> /dev/null; then
        echo -e "${RED}Error: $1 is required but not installed.${NC}"
        exit 1
    fi
}

# Main setup
main() {
    print_step "Checking prerequisites"
    check_command python3
    check_command pip
    check_command virtualenv

    print_step "Creating virtual environment"
    if [ -d "venv" ]; then
        echo "Removing existing virtual environment..."
        rm -rf venv
    fi
    python3 -m virtualenv venv

    print_step "Activating virtual environment"
    source venv/bin/activate

    print_step "Installing dependencies"
    pip install -r requirements-dev.txt
    pip install -e .

    print_step "Running tests"
    pytest tests/

    print_step "Running example enrichment"
    # Process example compounds with Swiss predictions only (no web data)
    enrich-compounds \
        --skip-web-data \
        --verbose \
        data/example_compounds.json \
        enriched_compounds.json

    print_step "Setup complete!"
    echo -e "${GREEN}The examples have been set up successfully!${NC}"
    echo
    echo "To activate the virtual environment:"
    echo "  source venv/bin/activate"
    echo
    echo "To run the enrichment script:"
    echo "  enrich-compounds data/example_compounds.json output.json"
    echo
    echo "Available options:"
    echo "  --reddit-id ID        Reddit client ID"
    echo "  --reddit-secret KEY   Reddit client secret"
    echo "  --twitter-token TOKEN Twitter bearer token"
    echo "  --skip-predictions    Skip Swiss predictions"
    echo "  --skip-web-data       Skip web data collection"
    echo "  --no-cache           Disable caching"
    echo "  --workers N          Number of worker threads (default: 4)"
    echo "  --batch-size N       Batch size for processing (default: 10)"
    echo "  --verbose            Enable debug logging"
    echo
    echo "Example with all data sources:"
    echo "  enrich-compounds \\"
    echo "    --reddit-id YOUR_REDDIT_ID \\"
    echo "    --reddit-secret YOUR_REDDIT_SECRET \\"
    echo "    --twitter-token YOUR_TWITTER_TOKEN \\"
    echo "    data/example_compounds.json \\"
    echo "    enriched_compounds.json"
}

# Run main function
main
