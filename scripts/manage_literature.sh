#!/bin/bash
# Script to manage literature data and references

# Exit on error
set -e

# Default values
DATA_DIR="data"
LITERATURE_DIR="literature"
LOG_DIR="logs"
SOURCE_TYPES="pubmed,chembl,patents,articles,books"
CONTENT_TYPES="abstracts,fulltext,references,citations,metadata"
MIN_CONFIDENCE=0.7
HARVEST=false
VALIDATE=false
CLEAN=false
BACKUP=false
RESTORE=false
FORCE=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --literature-dir)
            LITERATURE_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --source-types)
            SOURCE_TYPES="$2"
            shift 2
            ;;
        --content-types)
            CONTENT_TYPES="$2"
            shift 2
            ;;
        --min-confidence)
            MIN_CONFIDENCE="$2"
            shift 2
            ;;
        --harvest)
            HARVEST=true
            shift
            ;;
        --validate)
            VALIDATE=true
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
    
    # Literature directories
    for source in ${SOURCE_TYPES//,/ }; do
        for content in ${CONTENT_TYPES//,/ }; do
            mkdir -p "$LITERATURE_DIR/$source/$content"/{raw,processed,extracted}
        done
    done
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to harvest literature data
harvest_data() {
    echo "Harvesting literature data..."
    
    for source in ${SOURCE_TYPES//,/ }; do
        for content in ${CONTENT_TYPES//,/ }; do
            echo "Harvesting $content from $source..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli harvest-literature"
            CMD="$CMD --data-dir $DATA_DIR"
            CMD="$CMD --literature-dir $LITERATURE_DIR/$source/$content"
            CMD="$CMD --source $source"
            CMD="$CMD --content $content"
            CMD="$CMD --min-confidence $MIN_CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if harvesting fails
            
            # Generate harvest report
            echo "Generating harvest report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-literature-report \
                --literature-dir "$LITERATURE_DIR/$source/$content" \
                --output "$LITERATURE_DIR/$source/$content/raw/report.html"
        done
    done
}

# Function to validate literature data
validate_data() {
    echo "Validating literature data..."
    
    for source in ${SOURCE_TYPES//,/ }; do
        for content in ${CONTENT_TYPES//,/ }; do
            echo "Validating $content from $source..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli validate-literature"
            CMD="$CMD --literature-dir $LITERATURE_DIR/$source/$content"
            CMD="$CMD --source $source"
            CMD="$CMD --content $content"
            CMD="$CMD --min-confidence $MIN_CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if validation fails
            
            # Generate validation report
            echo "Generating validation report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-validation-report \
                --literature-dir "$LITERATURE_DIR/$source/$content" \
                --output "$LITERATURE_DIR/$source/$content/extracted/validation.html"
        done
    done
}

# Function to clean literature data
clean_literature() {
    echo "Cleaning literature data..."
    
    # Clean literature directories
    rm -rf "$LITERATURE_DIR"/*
}

# Function to backup literature data
backup_literature() {
    echo "Backing up literature data..."
    
    # Create backup directory
    BACKUP_DIR="backups/literature_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$LITERATURE_DIR" "$BACKUP_DIR/"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore literature data
restore_literature() {
    echo "Restoring literature data..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/literature_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$LITERATURE_DIR"
    fi
    
    cp -r "$BACKUP_DIR/literature" ./
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Literature data restored from: $BACKUP_DIR"
}

# Function to show literature statistics
show_stats() {
    echo "Literature statistics:"
    echo
    
    echo "Literature files:"
    for source in ${SOURCE_TYPES//,/ }; do
        echo "$source data:"
        for content in ${CONTENT_TYPES//,/ }; do
            echo "  $content data:"
            
            # Raw data
            echo "    Raw data:"
            echo "      Files: $(find "$LITERATURE_DIR/$source/$content/raw" -type f | wc -l) files"
            echo "      Size: $(du -sh "$LITERATURE_DIR/$source/$content/raw" | cut -f1)"
            
            # Count entries by content type
            case "$content" in
                abstracts)
                    if [ -f "$LITERATURE_DIR/$source/$content/raw/abstracts.json" ]; then
                        echo "      Abstracts: $(jq '.abstracts | length' "$LITERATURE_DIR/$source/$content/raw/abstracts.json")"
                        echo "      Papers: $(jq '.papers | length' "$LITERATURE_DIR/$source/$content/raw/abstracts.json")"
                    fi
                    ;;
                fulltext)
                    if [ -f "$LITERATURE_DIR/$source/$content/raw/fulltext.json" ]; then
                        echo "      Papers: $(jq '.papers | length' "$LITERATURE_DIR/$source/$content/raw/fulltext.json")"
                        echo "      Sections: $(jq '.sections | length' "$LITERATURE_DIR/$source/$content/raw/fulltext.json")"
                    fi
                    ;;
                references)
                    if [ -f "$LITERATURE_DIR/$source/$content/raw/references.json" ]; then
                        echo "      References: $(jq '.references | length' "$LITERATURE_DIR/$source/$content/raw/references.json")"
                        echo "      Papers: $(jq '.papers | length' "$LITERATURE_DIR/$source/$content/raw/references.json")"
                    fi
                    ;;
                citations)
                    if [ -f "$LITERATURE_DIR/$source/$content/raw/citations.json" ]; then
                        echo "      Citations: $(jq '.citations | length' "$LITERATURE_DIR/$source/$content/raw/citations.json")"
                        echo "      Papers: $(jq '.papers | length' "$LITERATURE_DIR/$source/$content/raw/citations.json")"
                    fi
                    ;;
                metadata)
                    if [ -f "$LITERATURE_DIR/$source/$content/raw/metadata.json" ]; then
                        echo "      Papers: $(jq '.papers | length' "$LITERATURE_DIR/$source/$content/raw/metadata.json")"
                        echo "      Fields: $(jq '.fields | length' "$LITERATURE_DIR/$source/$content/raw/metadata.json")"
                    fi
                    ;;
            esac
            
            # Processed data
            echo "    Processed data:"
            echo "      Files: $(find "$LITERATURE_DIR/$source/$content/processed" -type f | wc -l) files"
            echo "      Size: $(du -sh "$LITERATURE_DIR/$source/$content/processed" | cut -f1)"
            
            # Extracted data
            echo "    Extracted data:"
            echo "      Files: $(find "$LITERATURE_DIR/$source/$content/extracted" -type f | wc -l) files"
            echo "      Size: $(du -sh "$LITERATURE_DIR/$source/$content/extracted" | cut -f1)"
            
            # Show reports
            if [ -f "$LITERATURE_DIR/$source/$content/raw/report.html" ]; then
                echo "    Harvest report: $LITERATURE_DIR/$source/$content/raw/report.html"
            fi
            if [ -f "$LITERATURE_DIR/$source/$content/extracted/validation.html" ]; then
                echo "    Validation report: $LITERATURE_DIR/$source/$content/extracted/validation.html"
            fi
            echo
        done
    done
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing literature data..."

# Create directory structure
create_dirs

# Harvest data if requested
if [ "$HARVEST" = true ]; then
    harvest_data
fi

# Validate data if requested
if [ "$VALIDATE" = true ]; then
    validate_data
fi

# Clean data if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean literature data? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_literature
    fi
fi

# Backup data if requested
if [ "$BACKUP" = true ]; then
    backup_literature
fi

# Restore data if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore literature data? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_literature
    fi
fi

# Show statistics
show_stats

echo
echo "Literature management completed successfully!"
