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
    echo "Update project dependencies"
    echo
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  --check     Check for updates without installing"
    echo "  --dev       Update development dependencies only"
    echo "  --prod      Update production dependencies only"
    echo "  --models    Update ML models only"
    echo "  --system    Update system dependencies only"
    echo "  --hooks     Update pre-commit hooks only"
    echo "  --all       Update all dependencies (default)"
    echo "  --minimal   Skip ML-related updates"
    echo "  --help      Show this help message"
    echo
    echo "Examples:"
    echo "  $0 --check      Check for available updates"
    echo "  $0 --dev        Update development dependencies"
    echo "  $0 --prod       Update production dependencies"
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

# Function to check for uv
check_uv() {
    if ! command -v uv &> /dev/null; then
        echo -e "${RED}Error: uv not found${NC}"
        echo "Install uv first:"
        echo "curl -LsSf https://astral.sh/uv/install.sh | sh"
        exit 1
    fi
}

# Function to backup files
backup_files() {
    echo -e "${BLUE}Creating backups...${NC}"
    
    # Create backup directory with timestamp
    local backup_dir="backups/$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$backup_dir"
    
    # Backup requirements files
    for file in requirements*.txt; do
        if [[ -f "$file" ]]; then
            cp "$file" "$backup_dir/"
        fi
    done
    
    # Backup lock files
    if [[ -f "poetry.lock" ]]; then
        cp poetry.lock "$backup_dir/"
    fi
    
    # Backup model files if updating models
    if [[ $UPDATE_MODELS -eq 1 && ! "$CHECK_ONLY" == "true" ]]; then
        if [[ -d "models" ]]; then
            cp -r models "$backup_dir/"
        fi
    fi
    
    echo -e "${GREEN}Backups created in $backup_dir${NC}"
}

# Function to restore from backup
restore_from_backup() {
    local backup_dir="$1"
    echo -e "${YELLOW}Restoring from backup $backup_dir...${NC}"
    
    # Restore requirements files
    for file in "$backup_dir"/requirements*.txt; do
        if [[ -f "$file" ]]; then
            cp "$file" .
        fi
    done
    
    # Restore lock files
    if [[ -f "$backup_dir/poetry.lock" ]]; then
        cp "$backup_dir/poetry.lock" .
    fi
    
    # Restore model files
    if [[ -d "$backup_dir/models" ]]; then
        rm -rf models
        cp -r "$backup_dir/models" .
    fi
    
    echo -e "${GREEN}Restore complete${NC}"
}

# Function to update production dependencies
update_prod() {
    echo -e "${BLUE}Updating production dependencies...${NC}"
    
    if [[ "$CHECK_ONLY" == "true" ]]; then
        echo "Checking for updates..."
        uv pip list --outdated
    else
        # Update core dependencies
        uv pip install -U -e .
        
        # Update lockfile
        uv pip freeze > requirements.txt
        
        echo -e "${GREEN}Production dependencies updated${NC}"
    fi
}

# Function to update development dependencies
update_dev() {
    echo -e "${BLUE}Updating development dependencies...${NC}"
    
    if [[ "$CHECK_ONLY" == "true" ]]; then
        echo "Checking for updates..."
        uv pip list --outdated
    else
        # Update development dependencies
        uv pip install -U -e ".[dev]"
        
        # Update pre-commit hooks
        pre-commit autoupdate
        
        # Update development requirements
        uv pip freeze > requirements-dev.txt
        
        echo -e "${GREEN}Development dependencies updated${NC}"
    fi
}

# Function to update ML models
update_models() {
    if [[ "$MINIMAL" == "true" ]]; then
        echo -e "${YELLOW}Skipping ML model updates (minimal mode)${NC}"
        return 0
    fi
    
    echo -e "${BLUE}Updating ML models...${NC}"
    
    # Create models directory if it doesn't exist
    mkdir -p models/{toxicity,abuse,activity,affinity}
    
    if [[ "$CHECK_ONLY" == "true" ]]; then
        # Check for model updates
        echo -e "${BLUE}Checking for model updates...${NC}"
        for model in toxicity_predictor abuse_predictor activity_predictor affinity_predictor gnn_model ensemble_model; do
            if curl --silent --head "https://storage.googleapis.com/chemdata-models/$model.pt" | grep -q "200 OK"; then
                local_hash=""
                if [[ -f "models/$model.pt" ]]; then
                    local_hash=$(sha256sum "models/$model.pt" | cut -d' ' -f1)
                fi
                remote_hash=$(curl -s "https://storage.googleapis.com/chemdata-models/$model.pt.sha256")
                
                if [[ "$local_hash" != "$remote_hash" ]]; then
                    echo -e "${YELLOW}Update available for $model${NC}"
                else
                    echo -e "${GREEN}$model is up to date${NC}"
                fi
            else
                echo -e "${RED}Could not check $model${NC}"
            fi
        done
    else
        # Download/update models in parallel
        echo -e "${BLUE}Downloading models...${NC}"
        local pids=()
        for model in toxicity_predictor abuse_predictor activity_predictor affinity_predictor gnn_model ensemble_model; do
            (
                if curl -L "https://storage.googleapis.com/chemdata-models/$model.pt" -o "models/$model.pt.tmp"; then
                    mv "models/$model.pt.tmp" "models/$model.pt"
                    echo -e "${GREEN}Updated $model${NC}"
                else
                    echo -e "${RED}Failed to update $model${NC}"
                    rm -f "models/$model.pt.tmp"
                fi
            ) &
            pids+=($!)
        done
        
        # Wait for all downloads to complete
        for pid in "${pids[@]}"; do
            wait "$pid"
        done
    fi
}

# Function to update system dependencies
update_system() {
    echo -e "${BLUE}Updating system dependencies...${NC}"
    
    # Determine package manager
    if command -v apt-get &> /dev/null; then
        PKG_MANAGER="apt-get"
        if [[ "$CHECK_ONLY" == "true" ]]; then
            sudo apt-get update
            apt list --upgradable
        else
            sudo apt-get update
            sudo apt-get upgrade -y
        fi
    elif command -v brew &> /dev/null; then
        PKG_MANAGER="brew"
        if [[ "$CHECK_ONLY" == "true" ]]; then
            brew update
            brew outdated
        else
            brew update
            brew upgrade
        fi
    else
        echo -e "${RED}No supported package manager found${NC}"
        return 1
    fi
    
    # Update system packages
    local packages=(
        "cmake"
        "boost"
        "openbabel"
        "rdkit"
        "graphviz"
        "postgresql"
    )
    
    if [[ "$CHECK_ONLY" == "true" ]]; then
        echo -e "${BLUE}Checking system packages...${NC}"
        for package in "${packages[@]}"; do
            case $PKG_MANAGER in
                apt-get)
                    apt list --upgradable | grep -q "^$package/" && echo "$package needs update"
                    ;;
                brew)
                    brew outdated | grep -q "^$package" && echo "$package needs update"
                    ;;
            esac
        done
    else
        echo -e "${BLUE}Updating system packages...${NC}"
        case $PKG_MANAGER in
            apt-get)
                sudo apt-get install -y "${packages[@]}"
                ;;
            brew)
                brew install "${packages[@]}"
                ;;
        esac
    fi
}

# Function to update pre-commit hooks
update_hooks() {
    echo -e "${BLUE}Updating pre-commit hooks...${NC}"
    
    if [[ "$CHECK_ONLY" == "true" ]]; then
        pre-commit autoupdate --dry-run
    else
        pre-commit autoupdate
        pre-commit clean
        pre-commit install
        pre-commit install --hook-type pre-push
    fi
}

# Parse command line arguments
CHECK_ONLY=false
UPDATE_DEV=0
UPDATE_PROD=0
UPDATE_MODELS=0
UPDATE_SYSTEM=0
UPDATE_HOOKS=0
MINIMAL=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --check)
            CHECK_ONLY=true
            shift
            ;;
        --dev)
            UPDATE_DEV=1
            shift
            ;;
        --prod)
            UPDATE_PROD=1
            shift
            ;;
        --models)
            UPDATE_MODELS=1
            shift
            ;;
        --system)
            UPDATE_SYSTEM=1
            shift
            ;;
        --hooks)
            UPDATE_HOOKS=1
            shift
            ;;
        --all)
            UPDATE_DEV=1
            UPDATE_PROD=1
            UPDATE_MODELS=1
            UPDATE_SYSTEM=1
            UPDATE_HOOKS=1
            shift
            ;;
        --minimal)
            MINIMAL=true
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

# If no specific options provided, update all
if [[ $UPDATE_DEV -eq 0 && $UPDATE_PROD -eq 0 && $UPDATE_MODELS -eq 0 && \
      $UPDATE_SYSTEM -eq 0 && $UPDATE_HOOKS -eq 0 ]]; then
    UPDATE_DEV=1
    UPDATE_PROD=1
    UPDATE_MODELS=1
    UPDATE_SYSTEM=1
    UPDATE_HOOKS=1
fi

# Check environment
check_venv
check_uv

# Create backups if not just checking
if [[ "$CHECK_ONLY" != "true" ]]; then
    backup_files
fi

# Trap errors
trap 'echo -e "${RED}Error occurred. Rolling back...${NC}"; restore_from_backup "backups/$(ls -t backups | head -n1)"' ERR

# Run updates based on options
if [[ $UPDATE_SYSTEM -eq 1 ]]; then update_system; fi
if [[ $UPDATE_PROD -eq 1 ]]; then update_prod; fi
if [[ $UPDATE_DEV -eq 1 ]]; then update_dev; fi
if [[ $UPDATE_MODELS -eq 1 ]]; then update_models; fi
if [[ $UPDATE_HOOKS -eq 1 ]]; then update_hooks; fi

# Remove trap
trap - ERR

echo -e "${GREEN}Update complete!${NC}"

# Print next steps
if [[ "$CHECK_ONLY" == "true" ]]; then
    echo -e "\n${BLUE}To update dependencies, run:${NC}"
    echo "$0 --all"
else
    echo -e "\n${BLUE}Next steps:${NC}"
    echo "1. Review changes in requirements files"
    echo "2. Run tests: ./scripts/test.sh"
    echo "3. Start development server: ./scripts/dev.sh serve"
    echo "4. Commit changes if everything works"
fi
