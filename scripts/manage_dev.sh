#!/bin/bash
# Script to manage development environment and tools

# Exit on error
set -e

# Default values
VENV_DIR=".venv"
CONFIG_DIR=".config"
CACHE_DIR=".cache"
LOG_DIR="logs"
DEV_TOOLS="black,flake8,mypy,isort,bandit,pytest"
LINT=false
FORMAT=false
TYPE_CHECK=false
SECURITY=false
CLEAN=false
BACKUP=false
RESTORE=false
FIX=false
FORCE=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --venv-dir)
            VENV_DIR="$2"
            shift 2
            ;;
        --config-dir)
            CONFIG_DIR="$2"
            shift 2
            ;;
        --cache-dir)
            CACHE_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --dev-tools)
            DEV_TOOLS="$2"
            shift 2
            ;;
        --lint)
            LINT=true
            shift
            ;;
        --format)
            FORMAT=true
            shift
            ;;
        --type-check)
            TYPE_CHECK=true
            shift
            ;;
        --security)
            SECURITY=true
            shift
            ;;
        --clean)
            CLEAN=true
            shift
            ;;
        --backup)
            BACKUP=true
            shift
            ;;
        --restore)
            RESTORE=true
            shift
            ;;
        --fix)
            FIX=true
            shift
            ;;
        --force)
            FORCE=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to create directory structure
create_dirs() {
    echo "Creating directory structure..."
    
    # Development directories
    mkdir -p "$CONFIG_DIR"/{black,flake8,mypy,isort,bandit}
    mkdir -p "$CACHE_DIR"
    mkdir -p "$LOG_DIR"
}

# Function to install development tools
install_dev_tools() {
    echo "Installing development tools..."
    
    # Ensure pip is up to date
    pip install --upgrade pip setuptools wheel
    
    # Install tools
    for tool in ${DEV_TOOLS//,/ }; do
        echo "Installing $tool..."
        case "$tool" in
            black)
                pip install black
                ;;
            flake8)
                pip install flake8 flake8-docstrings flake8-bugbear flake8-comprehensions flake8-simplify
                ;;
            mypy)
                pip install mypy types-all
                ;;
            isort)
                pip install isort
                ;;
            bandit)
                pip install bandit safety
                ;;
            pytest)
                pip install pytest pytest-cov pytest-mock pytest-xdist pytest-timeout pytest-randomly pytest-sugar
                ;;
            *)
                echo "Unknown tool: $tool"
                ;;
        esac
    done
}

# Function to run linting
run_lint() {
    echo "Running linting..."
    
    # Run flake8
    if command_exists flake8; then
        echo "Running flake8..."
        if [ "$FIX" = true ]; then
            flake8 binding_data_processor tests scripts || true
        else
            flake8 binding_data_processor tests scripts
        fi
    fi
    
    # Run isort check
    if command_exists isort; then
        echo "Running isort..."
        if [ "$FIX" = true ]; then
            isort binding_data_processor tests scripts
        else
            isort --check-only binding_data_processor tests scripts
        fi
    fi
}

# Function to run formatting
run_format() {
    echo "Running formatting..."
    
    # Run black
    if command_exists black; then
        echo "Running black..."
        if [ "$FIX" = true ]; then
            black binding_data_processor tests scripts
        else
            black --check binding_data_processor tests scripts
        fi
    fi
}

# Function to run type checking
run_type_check() {
    echo "Running type checking..."
    
    # Run mypy
    if command_exists mypy; then
        echo "Running mypy..."
        if [ "$FIX" = true ]; then
            mypy binding_data_processor tests scripts || true
        else
            mypy binding_data_processor tests scripts
        fi
    fi
}

# Function to run security checks
run_security() {
    echo "Running security checks..."
    
    # Run bandit
    if command_exists bandit; then
        echo "Running bandit..."
        if [ "$FIX" = true ]; then
            bandit -r binding_data_processor tests scripts || true
        else
            bandit -r binding_data_processor tests scripts
        fi
    fi
    
    # Run safety
    if command_exists safety; then
        echo "Running safety..."
        safety check
    fi
}

# Function to clean development files
clean_dev() {
    echo "Cleaning development files..."
    
    # Clean cache directories
    rm -rf "$CACHE_DIR"/*
    rm -rf .mypy_cache/
    rm -rf .pytest_cache/
    
    # Clean Python cache
    find . -type d -name "__pycache__" -exec rm -rf {} +
    find . -type f -name "*.pyc" -delete
    find . -type f -name "*.pyo" -delete
    find . -type f -name "*.pyd" -delete
    
    # Clean editor files
    find . -type f -name ".*.swp" -delete
    find . -type f -name ".DS_Store" -delete
}

# Function to backup development files
backup_dev() {
    echo "Backing up development files..."
    
    # Create backup directory
    BACKUP_DIR="backups/dev_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$CONFIG_DIR" "$BACKUP_DIR/"
    cp -r "$CACHE_DIR" "$BACKUP_DIR/"
    
    # Copy config files
    cp -f .flake8 "$BACKUP_DIR/" 2>/dev/null || true
    cp -f .mypy.ini "$BACKUP_DIR/" 2>/dev/null || true
    cp -f pyproject.toml "$BACKUP_DIR/" 2>/dev/null || true
    cp -f setup.cfg "$BACKUP_DIR/" 2>/dev/null || true
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore development files
restore_dev() {
    echo "Restoring development files..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/dev_*.tar.gz 2>/dev/null | head -n1)
    if [ -z "$LATEST_BACKUP" ]; then
        echo "No backup found"
        exit 1
    fi
    
    # Extract backup
    echo "Extracting backup: $LATEST_BACKUP"
    tar -xzf "$LATEST_BACKUP"
    BACKUP_DIR="${LATEST_BACKUP%.tar.gz}"
    
    # Restore directories
    if [ "$FORCE" = true ]; then
        rm -rf "$CONFIG_DIR" "$CACHE_DIR"
    fi
    
    cp -r "$BACKUP_DIR/config" ./
    cp -r "$BACKUP_DIR/cache" ./
    
    # Restore config files
    cp -f "$BACKUP_DIR/.flake8" . 2>/dev/null || true
    cp -f "$BACKUP_DIR/.mypy.ini" . 2>/dev/null || true
    cp -f "$BACKUP_DIR/pyproject.toml" . 2>/dev/null || true
    cp -f "$BACKUP_DIR/setup.cfg" . 2>/dev/null || true
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Development files restored from: $BACKUP_DIR"
}

# Function to show development statistics
show_stats() {
    echo "Development statistics:"
    echo
    
    echo "Python files:"
    echo "  Source: $(find binding_data_processor -type f -name "*.py" | wc -l) files"
    echo "  Tests: $(find tests -type f -name "test_*.py" | wc -l) files"
    echo "  Scripts: $(find scripts -type f -name "*.py" -o -name "*.sh" | wc -l) files"
    echo
    
    echo "Code quality:"
    if command_exists flake8; then
        echo "  Flake8: $(flake8 binding_data_processor tests scripts 2>&1 | wc -l) issues"
    fi
    if command_exists black; then
        echo "  Black: $(black --check binding_data_processor tests scripts 2>&1 | wc -l) files would be reformatted"
    fi
    if command_exists mypy; then
        echo "  Mypy: $(mypy binding_data_processor tests scripts 2>&1 | wc -l) type issues"
    fi
    if command_exists bandit; then
        echo "  Bandit: $(bandit -r binding_data_processor tests scripts 2>&1 | grep "Issue" | wc -l) security issues"
    fi
    echo
    
    echo "Development files:"
    echo "  Config: $(find "$CONFIG_DIR" -type f | wc -l) files"
    echo "  Cache: $(find "$CACHE_DIR" -type f | wc -l) files"
    echo "  Logs: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing development environment..."

# Create directory structure
create_dirs

# Install development tools
install_dev_tools

# Run linting if requested
if [ "$LINT" = true ]; then
    run_lint
fi

# Run formatting if requested
if [ "$FORMAT" = true ]; then
    run_format
fi

# Run type checking if requested
if [ "$TYPE_CHECK" = true ]; then
    run_type_check
fi

# Run security checks if requested
if [ "$SECURITY" = true ]; then
    run_security
fi

# Clean development files if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean development files? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_dev
    fi
fi

# Backup development files if requested
if [ "$BACKUP" = true ]; then
    backup_dev
fi

# Restore development files if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore development files? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_dev
    fi
fi

# Show statistics
show_stats

echo
echo "Development management completed successfully!"
