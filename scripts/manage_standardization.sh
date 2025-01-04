#!/bin/bash
# Script to manage data standardization and normalization

# Exit on error
set -e

# Default values
DATA_DIR="data"
STANDARD_DIR="standardized"
LOG_DIR="logs"
STANDARD_TYPES="structure,activity,names,units,identifiers"
STANDARD_FORMATS="tsv,json,sdf"
STRICT=false
STANDARDIZE=false
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
        --standard-dir)
            STANDARD_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --standard-types)
            STANDARD_TYPES="$2"
            shift 2
            ;;
        --standard-formats)
            STANDARD_FORMATS="$2"
            shift 2
            ;;
        --strict)
            STRICT=true
            shift
            ;;
        --standardize)
            STANDARDIZE=true
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
    
    # Standardization directories
    for type in ${STANDARD_TYPES//,/ }; do
        for format in ${STANDARD_FORMATS//,/ }; do
            mkdir -p "$STANDARD_DIR/$type/$format"/{input,output,mapping}
        done
    done
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to run standardization
run_standardization() {
    echo "Running standardization..."
    
    for type in ${STANDARD_TYPES//,/ }; do
        for format in ${STANDARD_FORMATS//,/ }; do
            echo "Standardizing $type data in $format format..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli standardize-data"
            CMD="$CMD --data-dir $DATA_DIR"
            CMD="$CMD --standard-dir $STANDARD_DIR/$type/$format"
            CMD="$CMD --standard-type $type"
            CMD="$CMD --format $format"
            CMD="$CMD --log-dir $LOG_DIR"
            
            if [ "$STRICT" = true ]; then
                CMD="$CMD --strict"
            fi
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if standardization fails
            
            # Generate mapping report
            echo "Generating mapping report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-mapping-report \
                --standard-dir "$STANDARD_DIR/$type/$format" \
                --output "$STANDARD_DIR/$type/$format/mapping/report.html"
        done
    done
}

# Function to clean standardization
clean_standardization() {
    echo "Cleaning standardization..."
    
    # Clean standardization directories
    rm -rf "$STANDARD_DIR"/*
}

# Function to backup standardization
backup_standardization() {
    echo "Backing up standardization..."
    
    # Create backup directory
    BACKUP_DIR="backups/standard_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$STANDARD_DIR" "$BACKUP_DIR/"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore standardization
restore_standardization() {
    echo "Restoring standardization..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/standard_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$STANDARD_DIR"
    fi
    
    cp -r "$BACKUP_DIR/standardized" ./
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Standardization restored from: $BACKUP_DIR"
}

# Function to show standardization statistics
show_stats() {
    echo "Standardization statistics:"
    echo
    
    echo "Standardization files:"
    for type in ${STANDARD_TYPES//,/ }; do
        echo "$type standardization:"
        for format in ${STANDARD_FORMATS//,/ }; do
            echo "  $format format:"
            
            # Input files
            echo "    Input:"
            echo "      Files: $(find "$STANDARD_DIR/$type/$format/input" -type f | wc -l) files"
            echo "      Size: $(du -sh "$STANDARD_DIR/$type/$format/input" | cut -f1)"
            
            # Output files
            echo "    Output:"
            echo "      Files: $(find "$STANDARD_DIR/$type/$format/output" -type f | wc -l) files"
            echo "      Size: $(du -sh "$STANDARD_DIR/$type/$format/output" | cut -f1)"
            
            # Mapping files
            echo "    Mapping:"
            echo "      Files: $(find "$STANDARD_DIR/$type/$format/mapping" -type f | wc -l) files"
            
            # Count standardized entries
            case "$type" in
                structure)
                    echo "    Standardized structures: $(find "$STANDARD_DIR/$type/$format/output" -name "*.sdf" -exec grep -c '$$$$' {} \; 2>/dev/null | awk '{s+=$1} END {print s}' || echo 0)"
                    ;;
                activity)
                    echo "    Standardized activities: $(find "$STANDARD_DIR/$type/$format/output" -name "*.json" -exec jq '.activities | length' {} \; 2>/dev/null | awk '{s+=$1} END {print s}' || echo 0)"
                    ;;
                names)
                    echo "    Standardized names: $(find "$STANDARD_DIR/$type/$format/output" -name "*.json" -exec jq '.names | length' {} \; 2>/dev/null | awk '{s+=$1} END {print s}' || echo 0)"
                    ;;
                units)
                    echo "    Standardized units: $(find "$STANDARD_DIR/$type/$format/output" -name "*.json" -exec jq '.units | length' {} \; 2>/dev/null | awk '{s+=$1} END {print s}' || echo 0)"
                    ;;
                identifiers)
                    echo "    Standardized identifiers: $(find "$STANDARD_DIR/$type/$format/output" -name "*.json" -exec jq '.identifiers | length' {} \; 2>/dev/null | awk '{s+=$1} END {print s}' || echo 0)"
                    ;;
            esac
            
            # Show mapping report
            if [ -f "$STANDARD_DIR/$type/$format/mapping/report.html" ]; then
                echo "    Report: $STANDARD_DIR/$type/$format/mapping/report.html"
            fi
            echo
        done
    done
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing standardization..."

# Create directory structure
create_dirs

# Run standardization if requested
if [ "$STANDARDIZE" = true ]; then
    run_standardization
fi

# Clean standardization if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean standardization? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_standardization
    fi
fi

# Backup standardization if requested
if [ "$BACKUP" = true ]; then
    backup_standardization
fi

# Restore standardization if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore standardization? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_standardization
    fi
fi

# Show statistics
show_stats

echo
echo "Standardization management completed successfully!"
