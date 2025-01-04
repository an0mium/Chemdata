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
    echo "Clean development artifacts"
    echo
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  --all       Clean everything (including venv and caches)"
    echo "  --venv      Clean virtual environment only"
    echo "  --cache     Clean cache files only"
    echo "  --build     Clean build artifacts only"
    echo "  --data      Clean data files only"
    echo "  --models    Clean model files only"
    echo "  --logs      Clean log files only"
    echo "  --temp      Clean temporary files only"
    echo "  --reports   Clean report files only"
    echo "  --help      Show this help message"
    echo
    echo "Examples:"
    echo "  $0 --all          Clean everything"
    echo "  $0 --cache        Clean cache files only"
    echo "  $0 --build --data Clean build artifacts and data files"
}

# Function to get directory size
get_size() {
    if [[ -d "$1" ]]; then
        du -sh "$1" 2>/dev/null | cut -f1
    else
        echo "0B"
    fi
}

# Function to clean virtual environment
clean_venv() {
    echo -e "${BLUE}Cleaning virtual environment...${NC}"
    local size=$(get_size ".venv")
    rm -rf .venv/
    rm -rf .venv_uv/
    rm -rf .venv_pip/
    echo -e "${GREEN}Cleaned virtual environment (freed ${size})${NC}"
}

# Function to clean cache files
clean_cache() {
    echo -e "${BLUE}Cleaning cache files...${NC}"
    local total_size=0
    
    # Python cache
    find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
    find . -type f -name "*.pyc" -delete
    find . -type f -name "*.pyo" -delete
    find . -type f -name "*.pyd" -delete
    
    # Tool cache
    for cache in .pytest_cache .mypy_cache .ruff_cache .coverage htmlcov .tox .cache; do
        if [[ -d "$cache" ]]; then
            local size=$(get_size "$cache")
            total_size="$((total_size + $(echo $size | tr -d 'GMK')))"
            rm -rf "$cache"
            echo -e "${GREEN}Cleaned $cache (freed ${size})${NC}"
        fi
    done
    
    # Project cache
    for cache in web/.webassets-cache models/cache data/cache; do
        if [[ -d "$cache" ]]; then
            local size=$(get_size "$cache")
            total_size="$((total_size + $(echo $size | tr -d 'GMK')))"
            rm -rf "$cache"
            echo -e "${GREEN}Cleaned $cache (freed ${size})${NC}"
        fi
    done
    
    echo -e "${GREEN}Total cache space freed: ${total_size}MB${NC}"
}

# Function to clean build artifacts
clean_build() {
    echo -e "${BLUE}Cleaning build artifacts...${NC}"
    local total_size=0
    
    for dir in build dist *.egg-info .eggs pip-wheel-metadata; do
        if [[ -d "$dir" ]]; then
            local size=$(get_size "$dir")
            total_size="$((total_size + $(echo $size | tr -d 'GMK')))"
            rm -rf "$dir"
            echo -e "${GREEN}Cleaned $dir (freed ${size})${NC}"
        fi
    done
    
    echo -e "${GREEN}Total build space freed: ${total_size}MB${NC}"
}

# Function to clean data files
clean_data() {
    echo -e "${BLUE}Cleaning data files...${NC}"
    local total_size=0
    
    # Data directories
    for dir in data/raw data/processed data/interim data/external; do
        if [[ -d "$dir" ]]; then
            local size=$(get_size "$dir")
            total_size="$((total_size + $(echo $size | tr -d 'GMK')))"
            rm -rf "$dir"/*
            echo -e "${GREEN}Cleaned $dir (freed ${size})${NC}"
        fi
    done
    
    # Data files
    find . -type f -name "*.tsv" -delete
    find . -type f -name "*.csv" -delete
    find . -type f -name "*.json" ! -name "package.json" ! -name "tsconfig.json" -delete
    
    echo -e "${GREEN}Total data space freed: ${total_size}MB${NC}"
}

# Function to clean model files
clean_models() {
    echo -e "${BLUE}Cleaning model files...${NC}"
    local total_size=0
    
    # Model directories
    for dir in models/toxicity models/abuse models/activity models/affinity; do
        if [[ -d "$dir" ]]; then
            local size=$(get_size "$dir")
            total_size="$((total_size + $(echo $size | tr -d 'GMK')))"
            rm -rf "$dir"/*
            echo -e "${GREEN}Cleaned $dir (freed ${size})${NC}"
        fi
    done
    
    # Model files
    find . -type f -name "*.pt" -delete
    find . -type f -name "*.pth" -delete
    find . -type f -name "*.h5" -delete
    find . -type f -name "*.ckpt" -delete
    
    echo -e "${GREEN}Total model space freed: ${total_size}MB${NC}"
}

# Function to clean log files
clean_logs() {
    echo -e "${BLUE}Cleaning log files...${NC}"
    local total_size=0
    
    # Log directories
    for dir in logs reports/profiling; do
        if [[ -d "$dir" ]]; then
            local size=$(get_size "$dir")
            total_size="$((total_size + $(echo $size | tr -d 'GMK')))"
            rm -rf "$dir"/*
            echo -e "${GREEN}Cleaned $dir (freed ${size})${NC}"
        fi
    done
    
    # Log files
    find . -type f -name "*.log" -delete
    
    echo -e "${GREEN}Total log space freed: ${total_size}MB${NC}"
}

# Function to clean temporary files
clean_temp() {
    echo -e "${BLUE}Cleaning temporary files...${NC}"
    
    # Temp directories
    rm -rf temp/
    rm -rf tmp/
    
    # Editor files
    find . -type f -name "*.swp" -delete
    find . -type f -name "*.swo" -delete
    find . -type f -name "*~" -delete
    find . -type f -name "*.bak" -delete
    
    # System files
    find . -type f -name ".DS_Store" -delete
    
    echo -e "${GREEN}Temporary files cleaned${NC}"
}

# Function to clean report files
clean_reports() {
    echo -e "${BLUE}Cleaning report files...${NC}"
    local total_size=0
    
    for dir in reports/coverage reports/test-results reports/profiling; do
        if [[ -d "$dir" ]]; then
            local size=$(get_size "$dir")
            total_size="$((total_size + $(echo $size | tr -d 'GMK')))"
            rm -rf "$dir"/*
            echo -e "${GREEN}Cleaned $dir (freed ${size})${NC}"
        fi
    done
    
    echo -e "${GREEN}Total report space freed: ${total_size}MB${NC}"
}

# Parse command line arguments
CLEAN_ALL=0
CLEAN_VENV=0
CLEAN_CACHE=0
CLEAN_BUILD=0
CLEAN_DATA=0
CLEAN_MODELS=0
CLEAN_LOGS=0
CLEAN_TEMP=0
CLEAN_REPORTS=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --all)
            CLEAN_ALL=1
            shift
            ;;
        --venv)
            CLEAN_VENV=1
            shift
            ;;
        --cache)
            CLEAN_CACHE=1
            shift
            ;;
        --build)
            CLEAN_BUILD=1
            shift
            ;;
        --data)
            CLEAN_DATA=1
            shift
            ;;
        --models)
            CLEAN_MODELS=1
            shift
            ;;
        --logs)
            CLEAN_LOGS=1
            shift
            ;;
        --temp)
            CLEAN_TEMP=1
            shift
            ;;
        --reports)
            CLEAN_REPORTS=1
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

# If no specific options provided, show help
if [[ $CLEAN_ALL -eq 0 && $CLEAN_VENV -eq 0 && $CLEAN_CACHE -eq 0 && $CLEAN_BUILD -eq 0 && \
      $CLEAN_DATA -eq 0 && $CLEAN_MODELS -eq 0 && $CLEAN_LOGS -eq 0 && $CLEAN_TEMP -eq 0 && \
      $CLEAN_REPORTS -eq 0 ]]; then
    show_help
    exit 1
fi

# Confirm if cleaning everything
if [[ $CLEAN_ALL -eq 1 ]]; then
    echo -e "${YELLOW}WARNING: This will remove all development artifacts, including:${NC}"
    echo "- Virtual environment"
    echo "- Cache files"
    echo "- Build artifacts"
    echo "- Data files"
    echo "- Model files"
    echo "- Log files"
    echo "- Temporary files"
    echo "- Report files"
    echo
    read -p "Are you sure? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${BLUE}Operation cancelled${NC}"
        exit 0
    fi
fi

# Clean based on options
if [[ $CLEAN_ALL -eq 1 || $CLEAN_VENV -eq 1 ]]; then clean_venv; fi
if [[ $CLEAN_ALL -eq 1 || $CLEAN_CACHE -eq 1 ]]; then clean_cache; fi
if [[ $CLEAN_ALL -eq 1 || $CLEAN_BUILD -eq 1 ]]; then clean_build; fi
if [[ $CLEAN_ALL -eq 1 || $CLEAN_DATA -eq 1 ]]; then clean_data; fi
if [[ $CLEAN_ALL -eq 1 || $CLEAN_MODELS -eq 1 ]]; then clean_models; fi
if [[ $CLEAN_ALL -eq 1 || $CLEAN_LOGS -eq 1 ]]; then clean_logs; fi
if [[ $CLEAN_ALL -eq 1 || $CLEAN_TEMP -eq 1 ]]; then clean_temp; fi
if [[ $CLEAN_ALL -eq 1 || $CLEAN_REPORTS -eq 1 ]]; then clean_reports; fi

echo -e "${GREEN}Cleaning complete!${NC}"

# Recreate necessary directories
mkdir -p data/{raw,processed,interim,external}
mkdir -p models/{toxicity,abuse,activity,affinity}
mkdir -p logs
mkdir -p reports/{coverage,test-results,profiling}
mkdir -p .cache

echo -e "${BLUE}Recreated necessary directories${NC}"
