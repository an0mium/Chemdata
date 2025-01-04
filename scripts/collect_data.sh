#!/bin/bash
set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Function to show usage
show_help() {
    echo "Collect chemical compound data from various sources"
    echo
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  --bindingdb    Collect BindingDB data"
    echo "  --web          Collect web source data"
    echo "  --social       Collect social media data"
    echo "  --all          Collect all data (default)"
    echo "  --update       Update existing data"
    echo "  --force        Force data collection even if recent"
    echo "  --help         Show this help message"
    echo
    echo "Examples:"
    echo "  $0 --all          Collect all data"
    echo "  $0 --bindingdb    Collect BindingDB data only"
    echo "  $0 --web --update Update web source data"
}

# Function to check if in virtual environment
check_venv() {
    if [[ -z "${VIRTUAL_ENV}" ]]; then
        echo -e "${RED}Error: Not in a virtual environment${NC}"
        echo "Activate your virtual environment first:"
        echo "source .venv/bin/activate"
        exit 1
    fi
}

# Function to check API credentials
check_credentials() {
    local missing=()
    
    # Check environment variables
    for var in REDDIT_CLIENT_ID REDDIT_CLIENT_SECRET TWITTER_API_KEY TWITTER_API_SECRET; do
        if [[ -z "${!var}" ]]; then
            missing+=("$var")
        fi
    done
    
    if [ ${#missing[@]} -gt 0 ]; then
        echo -e "${RED}Error: Missing API credentials${NC}"
        echo "The following environment variables are not set:"
        printf '%s\n' "${missing[@]}"
        echo "Add them to .env file or export them"
        exit 1
    fi
}

# Function to check if data is recent
check_data_age() {
    local data_file=$1
    local max_age=$2  # in days
    
    if [[ -f "$data_file" ]]; then
        local file_age=$(( ($(date +%s) - $(date -r "$data_file" +%s)) / 86400 ))
        if [[ $file_age -lt $max_age ]]; then
            return 0  # Data is recent
        fi
    fi
    return 1  # Data is old or doesn't exist
}

# Function to collect BindingDB data
collect_bindingdb() {
    local force=$1
    local update=$2
    
    echo -e "${BLUE}Collecting BindingDB data...${NC}"
    
    # Check if data is recent (less than 7 days old)
    if [[ "$force" != "true" ]] && check_data_age "data/raw/bindingdb.tsv" 7; then
        if [[ "$update" != "true" ]]; then
            echo -e "${YELLOW}BindingDB data is recent. Use --force to collect anyway.${NC}"
            return 0
        fi
    fi
    
    # Run data collection
    python -m binding_data_processor.main collect-bindingdb \
        --output data/raw/bindingdb.tsv \
        $([ "$update" == "true" ] && echo "--update")
}

# Function to collect web source data
collect_web() {
    local force=$1
    local update=$2
    
    echo -e "${BLUE}Collecting web source data...${NC}"
    
    # Check if data is recent (less than 1 day old)
    if [[ "$force" != "true" ]] && check_data_age "data/raw/web_sources.json" 1; then
        if [[ "$update" != "true" ]]; then
            echo -e "${YELLOW}Web source data is recent. Use --force to collect anyway.${NC}"
            return 0
        fi
    fi
    
    # Run data collection
    python -m binding_data_processor.main collect-web \
        --output data/raw/web_sources.json \
        $([ "$update" == "true" ] && echo "--update")
}

# Function to collect social media data
collect_social() {
    local force=$1
    local update=$2
    
    echo -e "${BLUE}Collecting social media data...${NC}"
    
    # Check if data is recent (less than 1 day old)
    if [[ "$force" != "true" ]] && check_data_age "data/raw/social_media.json" 1; then
        if [[ "$update" != "true" ]]; then
            echo -e "${YELLOW}Social media data is recent. Use --force to collect anyway.${NC}"
            return 0
        fi
    fi
    
    # Run data collection
    python -m binding_data_processor.main collect-social \
        --output data/raw/social_media.json \
        $([ "$update" == "true" ] && echo "--update")
}

# Parse command line arguments
COLLECT_BINDINGDB=0
COLLECT_WEB=0
COLLECT_SOCIAL=0
FORCE=false
UPDATE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --bindingdb)
            COLLECT_BINDINGDB=1
            shift
            ;;
        --web)
            COLLECT_WEB=1
            shift
            ;;
        --social)
            COLLECT_SOCIAL=1
            shift
            ;;
        --all)
            COLLECT_BINDINGDB=1
            COLLECT_WEB=1
            COLLECT_SOCIAL=1
            shift
            ;;
        --force)
            FORCE=true
            shift
            ;;
        --update)
            UPDATE=true
            shift
            ;;
        --help|-h)
            show_help
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            show_help
            exit 1
            ;;
    esac
done

# If no specific options provided, collect all
if [[ $COLLECT_BINDINGDB -eq 0 && $COLLECT_WEB -eq 0 && $COLLECT_SOCIAL -eq 0 ]]; then
    COLLECT_BINDINGDB=1
    COLLECT_WEB=1
    COLLECT_SOCIAL=1
fi

# Check environment
check_venv
check_credentials

# Create data directories if they don't exist
mkdir -p data/{raw,processed,interim,external}

# Collect data based on options
if [[ $COLLECT_BINDINGDB -eq 1 ]]; then
    collect_bindingdb "$FORCE" "$UPDATE"
fi

if [[ $COLLECT_WEB -eq 1 ]]; then
    collect_web "$FORCE" "$UPDATE"
fi

if [[ $COLLECT_SOCIAL -eq 1 ]]; then
    collect_social "$FORCE" "$UPDATE"
fi

echo -e "${GREEN}Data collection complete!${NC}"

# Print next steps
echo -e "\n${BLUE}Next steps:${NC}"
echo "1. Process collected data: python -m binding_data_processor.main process"
echo "2. Update web interface: python -m web.app"
echo "3. Review collected data in data/raw/"
