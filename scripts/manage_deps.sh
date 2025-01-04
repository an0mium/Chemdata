#!/bin/bash
# Script to manage project dependencies and requirements

# Exit on error
set -e

# Default values
VENV_DIR=".venv"
CACHE_DIR=".cache"
LOG_DIR="logs"
REQUIREMENTS_FILE="requirements.txt"
DEV_REQUIREMENTS="requirements-dev.txt"
INSTALL=false
UPDATE=false
CHECK=false
CLEAN=false
BACKUP=false
RESTORE=false
DEV_MODE=false
FORCE=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --venv-dir)
            VENV_DIR="$2"
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
        --requirements)
            REQUIREMENTS_FILE="$2"
            shift 2
            ;;
        --dev-requirements)
            DEV_REQUIREMENTS="$2"
            shift 2
            ;;
        --install)
            INSTALL=true
            shift
            ;;
        --update)
            UPDATE=true
            shift
            ;;
        --check)
            CHECK=true
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
        --dev)
            DEV_MODE=true
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
    
    # Create directories
    mkdir -p "$CACHE_DIR"
    mkdir -p "$LOG_DIR"
}

# Function to install dependencies
install_deps() {
    echo "Installing dependencies..."
    
    # Ensure pip is up to date
    pip install --upgrade pip setuptools wheel
    
    # Install base requirements
    if [ -f "$REQUIREMENTS_FILE" ]; then
        echo "Installing base requirements..."
        pip install -r "$REQUIREMENTS_FILE"
    fi
    
    # Install development requirements
    if [ "$DEV_MODE" = true ] && [ -f "$DEV_REQUIREMENTS" ]; then
        echo "Installing development requirements..."
        pip install -r "$DEV_REQUIREMENTS"
    fi
    
    # Install project in editable mode
    pip install -e .
}

# Function to update dependencies
update_deps() {
    echo "Updating dependencies..."
    
    # Update base requirements
    if [ -f "$REQUIREMENTS_FILE" ]; then
        echo "Updating base requirements..."
        pip install -r "$REQUIREMENTS_FILE" --upgrade
        
        # Update requirements file
        pip freeze > "$REQUIREMENTS_FILE"
    fi
    
    # Update development requirements
    if [ "$DEV_MODE" = true ] && [ -f "$DEV_REQUIREMENTS" ]; then
        echo "Updating development requirements..."
        pip install -r "$DEV_REQUIREMENTS" --upgrade
        
        # Update requirements file
        pip freeze > "$DEV_REQUIREMENTS"
    fi
}

# Function to check dependencies
check_deps() {
    echo "Checking dependencies..."
    
    # Check base requirements
    if [ -f "$REQUIREMENTS_FILE" ]; then
        echo "Checking base requirements..."
        pip check
        
        # Check for outdated packages
        pip list --outdated
        
        # Check for security vulnerabilities
        if command_exists safety; then
            safety check
        fi
    fi
    
    # Check development requirements
    if [ "$DEV_MODE" = true ] && [ -f "$DEV_REQUIREMENTS" ]; then
        echo "Checking development requirements..."
        pip check -r "$DEV_REQUIREMENTS"
    fi
}

# Function to clean dependencies
clean_deps() {
    echo "Cleaning dependencies..."
    
    # Clean pip cache
    pip cache purge
    
    # Clean project cache
    rm -rf "$CACHE_DIR"/*
    rm -rf .eggs/
    rm -rf *.egg-info/
    rm -rf build/
    rm -rf dist/
    
    # Clean virtual environment
    if [ "$FORCE" = true ]; then
        rm -rf "$VENV_DIR"
    fi
}

# Function to backup dependencies
backup_deps() {
    echo "Backing up dependencies..."
    
    # Create backup directory
    BACKUP_DIR="backups/deps_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy requirements files
    cp -f "$REQUIREMENTS_FILE" "$BACKUP_DIR/" 2>/dev/null || true
    cp -f "$DEV_REQUIREMENTS" "$BACKUP_DIR/" 2>/dev/null || true
    
    # Export current environment
    pip freeze > "$BACKUP_DIR/environment.txt"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore dependencies
restore_deps() {
    echo "Restoring dependencies..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/deps_*.tar.gz 2>/dev/null | head -n1)
    if [ -z "$LATEST_BACKUP" ]; then
        echo "No backup found"
        exit 1
    fi
    
    # Extract backup
    echo "Extracting backup: $LATEST_BACKUP"
    tar -xzf "$LATEST_BACKUP"
    BACKUP_DIR="${LATEST_BACKUP%.tar.gz}"
    
    # Restore requirements files
    if [ "$FORCE" = true ]; then
        cp -f "$BACKUP_DIR/requirements.txt" . 2>/dev/null || true
        cp -f "$BACKUP_DIR/requirements-dev.txt" . 2>/dev/null || true
    fi
    
    # Restore environment
    if [ -f "$BACKUP_DIR/environment.txt" ]; then
        pip install -r "$BACKUP_DIR/environment.txt"
    fi
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Dependencies restored from: $BACKUP_DIR"
}

# Function to show dependency statistics
show_stats() {
    echo "Dependency statistics:"
    echo
    
    echo "Package counts:"
    echo "  Base: $(pip freeze | wc -l) packages"
    if [ "$DEV_MODE" = true ]; then
        echo "  Development: $(pip freeze -r "$DEV_REQUIREMENTS" 2>/dev/null | wc -l) packages"
    fi
    echo
    
    echo "Storage:"
    echo "  Cache: $(du -sh "$CACHE_DIR" | cut -f1)"
    if [ -d "$VENV_DIR" ]; then
        echo "  Virtual environment: $(du -sh "$VENV_DIR" | cut -f1)"
    fi
    echo
    
    echo "Updates available:"
    pip list --outdated | tail -n +3
    echo
    
    if command_exists safety; then
        echo "Security issues:"
        safety check || true
    fi
}

# Main process
echo "Managing dependencies..."

# Create directory structure
create_dirs

# Install dependencies if requested
if [ "$INSTALL" = true ]; then
    install_deps
fi

# Update dependencies if requested
if [ "$UPDATE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to update dependencies? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        update_deps
    fi
fi

# Check dependencies if requested
if [ "$CHECK" = true ]; then
    check_deps
fi

# Clean dependencies if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean dependencies? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_deps
    fi
fi

# Backup dependencies if requested
if [ "$BACKUP" = true ]; then
    backup_deps
fi

# Restore dependencies if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore dependencies? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_deps
    fi
fi

# Show statistics
show_stats

echo
echo "Dependency management completed successfully!"
