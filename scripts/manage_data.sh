#!/bin/bash
# Script to manage project data files and directories

# Exit on error
set -e

# Default values
DATA_DIR="data"
CACHE_DIR="cache"
CHECKPOINT_DIR="checkpoints"
MODEL_DIR="models"
REPORT_DIR="reports"
LOG_DIR="logs"
CLEAN_CACHE=false
CLEAN_CHECKPOINTS=false
CLEAN_REPORTS=false
CLEAN_LOGS=false
CLEAN_ALL=false
BACKUP=false
RESTORE=false
COMPRESS=false
FORCE=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --cache-dir)
            CACHE_DIR="$2"
            shift 2
            ;;
        --checkpoint-dir)
            CHECKPOINT_DIR="$2"
            shift 2
            ;;
        --model-dir)
            MODEL_DIR="$2"
            shift 2
            ;;
        --report-dir)
            REPORT_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --clean-cache)
            CLEAN_CACHE=true
            shift
            ;;
        --clean-checkpoints)
            CLEAN_CHECKPOINTS=true
            shift
            ;;
        --clean-reports)
            CLEAN_REPORTS=true
            shift
            ;;
        --clean-logs)
            CLEAN_LOGS=true
            shift
            ;;
        --clean-all)
            CLEAN_ALL=true
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
        --compress)
            COMPRESS=true
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

# Function to create directory structure
create_dirs() {
    echo "Creating directory structure..."
    
    # Data directories
    mkdir -p "$DATA_DIR"/{raw,processed,interim,external}
    
    # Model directories
    mkdir -p "$MODEL_DIR"/{activity,toxicity,abuse,bbb}
    
    # Report directories
    mkdir -p "$REPORT_DIR"/{figures,tables,html,pdf}
    
    # Other directories
    mkdir -p "$CACHE_DIR"
    mkdir -p "$CHECKPOINT_DIR"
    mkdir -p "$LOG_DIR"
}

# Function to clean directories
clean_dirs() {
    if [ "$CLEAN_ALL" = true ] || [ "$CLEAN_CACHE" = true ]; then
        echo "Cleaning cache directory..."
        rm -rf "$CACHE_DIR"/*
    fi
    
    if [ "$CLEAN_ALL" = true ] || [ "$CLEAN_CHECKPOINTS" = true ]; then
        echo "Cleaning checkpoints directory..."
        rm -rf "$CHECKPOINT_DIR"/*
    fi
    
    if [ "$CLEAN_ALL" = true ] || [ "$CLEAN_REPORTS" = true ]; then
        echo "Cleaning reports directory..."
        rm -rf "$REPORT_DIR"/*
    fi
    
    if [ "$CLEAN_ALL" = true ] || [ "$CLEAN_LOGS" = true ]; then
        echo "Cleaning logs directory..."
        rm -rf "$LOG_DIR"/*
    fi
}

# Function to backup data
backup_data() {
    echo "Backing up data..."
    
    # Create backup directory
    BACKUP_DIR="backups/$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$DATA_DIR" "$BACKUP_DIR/"
    cp -r "$MODEL_DIR" "$BACKUP_DIR/"
    cp -r "$REPORT_DIR" "$BACKUP_DIR/"
    
    if [ "$COMPRESS" = true ]; then
        echo "Compressing backup..."
        tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
        rm -rf "$BACKUP_DIR"
        echo "Backup saved to: $BACKUP_DIR.tar.gz"
    else
        echo "Backup saved to: $BACKUP_DIR"
    fi
}

# Function to restore data
restore_data() {
    echo "Restoring data..."
    
    # Find latest backup
    if [ "$COMPRESS" = true ]; then
        LATEST_BACKUP=$(ls -t backups/*.tar.gz 2>/dev/null | head -n1)
        if [ -z "$LATEST_BACKUP" ]; then
            echo "No backup found"
            exit 1
        fi
        
        # Extract backup
        echo "Extracting backup: $LATEST_BACKUP"
        tar -xzf "$LATEST_BACKUP"
        BACKUP_DIR="${LATEST_BACKUP%.tar.gz}"
    else
        LATEST_BACKUP=$(ls -td backups/*/ 2>/dev/null | head -n1)
        if [ -z "$LATEST_BACKUP" ]; then
            echo "No backup found"
            exit 1
        fi
        BACKUP_DIR="${LATEST_BACKUP%/}"
    fi
    
    # Restore directories
    if [ "$FORCE" = true ]; then
        rm -rf "$DATA_DIR" "$MODEL_DIR" "$REPORT_DIR"
    fi
    
    cp -r "$BACKUP_DIR/data" ./
    cp -r "$BACKUP_DIR/models" ./
    cp -r "$BACKUP_DIR/reports" ./
    
    if [ "$COMPRESS" = true ]; then
        rm -rf "$BACKUP_DIR"
    fi
    
    echo "Data restored from: $BACKUP_DIR"
}

# Function to show data statistics
show_stats() {
    echo "Data statistics:"
    echo
    echo "Data directory ($DATA_DIR):"
    echo "  Raw data: $(find "$DATA_DIR/raw" -type f | wc -l) files ($(du -sh "$DATA_DIR/raw" | cut -f1))"
    echo "  Processed data: $(find "$DATA_DIR/processed" -type f | wc -l) files ($(du -sh "$DATA_DIR/processed" | cut -f1))"
    echo "  Interim data: $(find "$DATA_DIR/interim" -type f | wc -l) files ($(du -sh "$DATA_DIR/interim" | cut -f1))"
    echo "  External data: $(find "$DATA_DIR/external" -type f | wc -l) files ($(du -sh "$DATA_DIR/external" | cut -f1))"
    echo
    echo "Model directory ($MODEL_DIR):"
    echo "  Activity models: $(find "$MODEL_DIR/activity" -type f | wc -l) files ($(du -sh "$MODEL_DIR/activity" | cut -f1))"
    echo "  Toxicity models: $(find "$MODEL_DIR/toxicity" -type f | wc -l) files ($(du -sh "$MODEL_DIR/toxicity" | cut -f1))"
    echo "  Abuse models: $(find "$MODEL_DIR/abuse" -type f | wc -l) files ($(du -sh "$MODEL_DIR/abuse" | cut -f1))"
    echo "  BBB models: $(find "$MODEL_DIR/bbb" -type f | wc -l) files ($(du -sh "$MODEL_DIR/bbb" | cut -f1))"
    echo
    echo "Report directory ($REPORT_DIR):"
    echo "  Figures: $(find "$REPORT_DIR/figures" -type f | wc -l) files ($(du -sh "$REPORT_DIR/figures" | cut -f1))"
    echo "  Tables: $(find "$REPORT_DIR/tables" -type f | wc -l) files ($(du -sh "$REPORT_DIR/tables" | cut -f1))"
    echo "  HTML reports: $(find "$REPORT_DIR/html" -type f | wc -l) files ($(du -sh "$REPORT_DIR/html" | cut -f1))"
    echo "  PDF reports: $(find "$REPORT_DIR/pdf" -type f | wc -l) files ($(du -sh "$REPORT_DIR/pdf" | cut -f1))"
    echo
    echo "Other directories:"
    echo "  Cache: $(find "$CACHE_DIR" -type f | wc -l) files ($(du -sh "$CACHE_DIR" | cut -f1))"
    echo "  Checkpoints: $(find "$CHECKPOINT_DIR" -type f | wc -l) files ($(du -sh "$CHECKPOINT_DIR" | cut -f1))"
    echo "  Logs: $(find "$LOG_DIR" -type f | wc -l) files ($(du -sh "$LOG_DIR" | cut -f1))"
}

# Main process
echo "Managing project data..."

# Create directory structure
create_dirs

# Clean directories if requested
if [ "$CLEAN_ALL" = true ] || [ "$CLEAN_CACHE" = true ] || [ "$CLEAN_CHECKPOINTS" = true ] || [ "$CLEAN_REPORTS" = true ] || [ "$CLEAN_LOGS" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean directories? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_dirs
    fi
fi

# Backup data if requested
if [ "$BACKUP" = true ]; then
    backup_data
fi

# Restore data if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore data? This will overwrite existing data. [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_data
    fi
fi

# Show statistics
show_stats

echo
echo "Data management completed successfully!"
