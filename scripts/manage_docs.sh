#!/bin/bash
# Script to manage project documentation and reports

# Exit on error
set -e

# Default values
DOCS_DIR="docs"
REPORT_DIR="reports"
NOTEBOOK_DIR="notebooks"
OUTPUT_DIR="site"
FORMAT="html"
DOC_TYPES="api,models,data,web,cli"
BUILD=false
SERVE=false
CLEAN=false
BACKUP=false
RESTORE=false
COMPRESS=false
DEV_MODE=false
FORCE=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --docs-dir)
            DOCS_DIR="$2"
            shift 2
            ;;
        --report-dir)
            REPORT_DIR="$2"
            shift 2
            ;;
        --notebook-dir)
            NOTEBOOK_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --format)
            FORMAT="$2"
            shift 2
            ;;
        --doc-types)
            DOC_TYPES="$2"
            shift 2
            ;;
        --build)
            BUILD=true
            shift
            ;;
        --serve)
            SERVE=true
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
        --compress)
            COMPRESS=true
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

# Function to create directory structure
create_dirs() {
    echo "Creating directory structure..."
    
    # Documentation directories
    mkdir -p "$DOCS_DIR"/{api,models,data,web,cli}/{_static,_templates}
    
    # Report directories
    mkdir -p "$REPORT_DIR"/{figures,tables,html,pdf}
    
    # Notebook directories
    mkdir -p "$NOTEBOOK_DIR"/{examples,tutorials,analysis}
    
    # Output directory
    mkdir -p "$OUTPUT_DIR"
}

# Function to build documentation
build_docs() {
    echo "Building documentation..."
    
    for type in ${DOC_TYPES//,/ }; do
        echo "Building $type documentation..."
        
        # Build command
        CMD="sphinx-build"
        if [ "$DEV_MODE" = true ]; then
            CMD="$CMD -W -n"  # Warnings as errors, nitpicky mode
        fi
        
        CMD="$CMD -b $FORMAT"  # Output format
        CMD="$CMD $DOCS_DIR/$type"  # Source directory
        CMD="$CMD $OUTPUT_DIR/$type"  # Output directory
        
        # Run command
        echo "Running: $CMD"
        $CMD
    done
    
    # Build notebooks if they exist
    if [ -d "$NOTEBOOK_DIR" ]; then
        echo "Building notebooks..."
        jupyter nbconvert \
            --to html \
            --output-dir "$OUTPUT_DIR/notebooks" \
            "$NOTEBOOK_DIR"/**/*.ipynb
    fi
}

# Function to serve documentation
serve_docs() {
    echo "Serving documentation..."
    
    # Build documentation first
    build_docs
    
    # Start server
    if command -v python3 &>/dev/null; then
        echo "Starting server at http://localhost:8000"
        (cd "$OUTPUT_DIR" && python3 -m http.server 8000)
    else
        echo "Error: Python 3 not found"
        exit 1
    fi
}

# Function to clean documentation
clean_docs() {
    echo "Cleaning documentation..."
    
    # Clean output directory
    rm -rf "$OUTPUT_DIR"/*
    
    # Clean build artifacts
    find "$DOCS_DIR" -type d -name "_build" -exec rm -rf {} +
    find "$DOCS_DIR" -type d -name ".doctrees" -exec rm -rf {} +
    find "$DOCS_DIR" -type f -name "*.pyc" -delete
    find "$DOCS_DIR" -type f -name ".DS_Store" -delete
    
    # Clean notebook artifacts
    find "$NOTEBOOK_DIR" -type f -name ".ipynb_checkpoints" -delete
}

# Function to backup documentation
backup_docs() {
    echo "Backing up documentation..."
    
    # Create backup directory
    BACKUP_DIR="backups/docs_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$DOCS_DIR" "$BACKUP_DIR/"
    cp -r "$REPORT_DIR" "$BACKUP_DIR/"
    cp -r "$NOTEBOOK_DIR" "$BACKUP_DIR/"
    
    if [ "$COMPRESS" = true ]; then
        echo "Compressing backup..."
        tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
        rm -rf "$BACKUP_DIR"
        echo "Backup saved to: $BACKUP_DIR.tar.gz"
    else
        echo "Backup saved to: $BACKUP_DIR"
    fi
}

# Function to restore documentation
restore_docs() {
    echo "Restoring documentation..."
    
    # Find latest backup
    if [ "$COMPRESS" = true ]; then
        LATEST_BACKUP=$(ls -t backups/docs_*.tar.gz 2>/dev/null | head -n1)
        if [ -n "$LATEST_BACKUP" ]; then
            echo "Extracting backup: $LATEST_BACKUP"
            tar -xzf "$LATEST_BACKUP"
            BACKUP_DIR="${LATEST_BACKUP%.tar.gz}"
        fi
    else
        LATEST_BACKUP=$(ls -td backups/docs_*/ 2>/dev/null | head -n1)
        if [ -n "$LATEST_BACKUP" ]; then
            BACKUP_DIR="${LATEST_BACKUP%/}"
        fi
    fi
    
    if [ -z "$BACKUP_DIR" ]; then
        echo "No backup found"
        exit 1
    fi
    
    # Restore directories
    if [ "$FORCE" = true ]; then
        rm -rf "$DOCS_DIR" "$REPORT_DIR" "$NOTEBOOK_DIR"
    fi
    
    cp -r "$BACKUP_DIR/docs" ./
    cp -r "$BACKUP_DIR/reports" ./
    cp -r "$BACKUP_DIR/notebooks" ./
    
    if [ "$COMPRESS" = true ]; then
        rm -rf "$BACKUP_DIR"
    fi
    
    echo "Documentation restored from: $BACKUP_DIR"
}

# Function to show documentation statistics
show_stats() {
    echo "Documentation statistics:"
    echo
    
    echo "Documentation files:"
    for type in ${DOC_TYPES//,/ }; do
        echo "  $type: $(find "$DOCS_DIR/$type" -type f -name "*.rst" | wc -l) files"
        echo "    Source: $(find "$DOCS_DIR/$type" -type f -name "*.rst" -exec cat {} \; | wc -l) lines"
        if [ -d "$OUTPUT_DIR/$type" ]; then
            echo "    Built: $(du -sh "$OUTPUT_DIR/$type" | cut -f1)"
        fi
    done
    echo
    
    echo "Reports:"
    echo "  Figures: $(find "$REPORT_DIR/figures" -type f | wc -l) files ($(du -sh "$REPORT_DIR/figures" | cut -f1))"
    echo "  Tables: $(find "$REPORT_DIR/tables" -type f | wc -l) files ($(du -sh "$REPORT_DIR/tables" | cut -f1))"
    echo "  HTML: $(find "$REPORT_DIR/html" -type f | wc -l) files ($(du -sh "$REPORT_DIR/html" | cut -f1))"
    echo "  PDF: $(find "$REPORT_DIR/pdf" -type f | wc -l) files ($(du -sh "$REPORT_DIR/pdf" | cut -f1))"
    echo
    
    echo "Notebooks:"
    echo "  Examples: $(find "$NOTEBOOK_DIR/examples" -type f -name "*.ipynb" | wc -l) notebooks"
    echo "  Tutorials: $(find "$NOTEBOOK_DIR/tutorials" -type f -name "*.ipynb" | wc -l) notebooks"
    echo "  Analysis: $(find "$NOTEBOOK_DIR/analysis" -type f -name "*.ipynb" | wc -l) notebooks"
}

# Main process
echo "Managing documentation..."

# Create directory structure
create_dirs

# Build documentation if requested
if [ "$BUILD" = true ]; then
    build_docs
fi

# Serve documentation if requested
if [ "$SERVE" = true ]; then
    serve_docs
fi

# Clean documentation if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean documentation? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_docs
    fi
fi

# Backup documentation if requested
if [ "$BACKUP" = true ]; then
    backup_docs
fi

# Restore documentation if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore documentation? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_docs
    fi
fi

# Show statistics
show_stats

echo
echo "Documentation management completed successfully!"
