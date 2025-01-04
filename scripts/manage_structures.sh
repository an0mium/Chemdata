#!/bin/bash
# Script to manage chemical structure data and formats

# Exit on error
set -e

# Default values
DATA_DIR="data"
STRUCTURE_DIR="structures"
LOG_DIR="logs"
FORMAT_TYPES="smiles,sdf,mol,inchi,cdxml,pdb"
OPERATION_TYPES="convert,standardize,validate,analyze,depict"
MIN_CONFIDENCE=0.7
PROCESS=false
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
        --structure-dir)
            STRUCTURE_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --format-types)
            FORMAT_TYPES="$2"
            shift 2
            ;;
        --operation-types)
            OPERATION_TYPES="$2"
            shift 2
            ;;
        --min-confidence)
            MIN_CONFIDENCE="$2"
            shift 2
            ;;
        --process)
            PROCESS=true
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
    
    # Structure directories
    for format in ${FORMAT_TYPES//,/ }; do
        for operation in ${OPERATION_TYPES//,/ }; do
            mkdir -p "$STRUCTURE_DIR/$format/$operation"/{input,output,reports}
        done
    done
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to process structures
process_structures() {
    echo "Processing structures..."
    
    for format in ${FORMAT_TYPES//,/ }; do
        for operation in ${OPERATION_TYPES//,/ }; do
            echo "Performing $operation on $format structures..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli process-structures"
            CMD="$CMD --data-dir $DATA_DIR"
            CMD="$CMD --structure-dir $STRUCTURE_DIR/$format/$operation"
            CMD="$CMD --format $format"
            CMD="$CMD --operation $operation"
            CMD="$CMD --min-confidence $MIN_CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if processing fails
            
            # Generate processing report
            echo "Generating processing report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-structure-report \
                --structure-dir "$STRUCTURE_DIR/$format/$operation" \
                --output "$STRUCTURE_DIR/$format/$operation/reports/report.html"
        done
    done
}

# Function to validate structures
validate_structures() {
    echo "Validating structures..."
    
    for format in ${FORMAT_TYPES//,/ }; do
        for operation in ${OPERATION_TYPES//,/ }; do
            echo "Validating $operation results for $format structures..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli validate-structures"
            CMD="$CMD --structure-dir $STRUCTURE_DIR/$format/$operation"
            CMD="$CMD --format $format"
            CMD="$CMD --operation $operation"
            CMD="$CMD --min-confidence $MIN_CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if validation fails
            
            # Generate validation report
            echo "Generating validation report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-validation-report \
                --structure-dir "$STRUCTURE_DIR/$format/$operation" \
                --output "$STRUCTURE_DIR/$format/$operation/reports/validation.html"
        done
    done
}

# Function to clean structures
clean_structures() {
    echo "Cleaning structures..."
    
    # Clean structure directories
    rm -rf "$STRUCTURE_DIR"/*
}

# Function to backup structures
backup_structures() {
    echo "Backing up structures..."
    
    # Create backup directory
    BACKUP_DIR="backups/structures_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$STRUCTURE_DIR" "$BACKUP_DIR/"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore structures
restore_structures() {
    echo "Restoring structures..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/structures_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$STRUCTURE_DIR"
    fi
    
    cp -r "$BACKUP_DIR/structures" ./
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Structures restored from: $BACKUP_DIR"
}

# Function to show structure statistics
show_stats() {
    echo "Structure statistics:"
    echo
    
    echo "Structure files:"
    for format in ${FORMAT_TYPES//,/ }; do
        echo "$format structures:"
        for operation in ${OPERATION_TYPES//,/ }; do
            echo "  $operation operation:"
            
            # Input files
            echo "    Input:"
            echo "      Files: $(find "$STRUCTURE_DIR/$format/$operation/input" -type f | wc -l) files"
            echo "      Size: $(du -sh "$STRUCTURE_DIR/$format/$operation/input" | cut -f1)"
            
            # Count structures by operation
            case "$operation" in
                convert)
                    if [ -f "$STRUCTURE_DIR/$format/$operation/output/conversion.json" ]; then
                        echo "      Converted: $(jq '.converted | length' "$STRUCTURE_DIR/$format/$operation/output/conversion.json")"
                        echo "      Failed: $(jq '.failed | length' "$STRUCTURE_DIR/$format/$operation/output/conversion.json")"
                    fi
                    ;;
                standardize)
                    if [ -f "$STRUCTURE_DIR/$format/$operation/output/standardization.json" ]; then
                        echo "      Standardized: $(jq '.standardized | length' "$STRUCTURE_DIR/$format/$operation/output/standardization.json")"
                        echo "      Modified: $(jq '.modified | length' "$STRUCTURE_DIR/$format/$operation/output/standardization.json")"
                    fi
                    ;;
                validate)
                    if [ -f "$STRUCTURE_DIR/$format/$operation/output/validation.json" ]; then
                        echo "      Valid: $(jq '.valid | length' "$STRUCTURE_DIR/$format/$operation/output/validation.json")"
                        echo "      Invalid: $(jq '.invalid | length' "$STRUCTURE_DIR/$format/$operation/output/validation.json")"
                    fi
                    ;;
                analyze)
                    if [ -f "$STRUCTURE_DIR/$format/$operation/output/analysis.json" ]; then
                        echo "      Properties: $(jq '.properties | length' "$STRUCTURE_DIR/$format/$operation/output/analysis.json")"
                        echo "      Features: $(jq '.features | length' "$STRUCTURE_DIR/$format/$operation/output/analysis.json")"
                    fi
                    ;;
                depict)
                    if [ -f "$STRUCTURE_DIR/$format/$operation/output/depiction.json" ]; then
                        echo "      2D: $(jq '.2d | length' "$STRUCTURE_DIR/$format/$operation/output/depiction.json")"
                        echo "      3D: $(jq '.3d | length' "$STRUCTURE_DIR/$format/$operation/output/depiction.json")"
                    fi
                    ;;
            esac
            
            # Output files
            echo "    Output:"
            echo "      Files: $(find "$STRUCTURE_DIR/$format/$operation/output" -type f | wc -l) files"
            echo "      Size: $(du -sh "$STRUCTURE_DIR/$format/$operation/output" | cut -f1)"
            
            # Show reports
            if [ -f "$STRUCTURE_DIR/$format/$operation/reports/report.html" ]; then
                echo "    Processing report: $STRUCTURE_DIR/$format/$operation/reports/report.html"
            fi
            if [ -f "$STRUCTURE_DIR/$format/$operation/reports/validation.html" ]; then
                echo "    Validation report: $STRUCTURE_DIR/$format/$operation/reports/validation.html"
            fi
            echo
        done
    done
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing structures..."

# Create directory structure
create_dirs

# Process structures if requested
if [ "$PROCESS" = true ]; then
    process_structures
fi

# Validate structures if requested
if [ "$VALIDATE" = true ]; then
    validate_structures
fi

# Clean structures if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean structures? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_structures
    fi
fi

# Backup structures if requested
if [ "$BACKUP" = true ]; then
    backup_structures
fi

# Restore structures if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore structures? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_structures
    fi
fi

# Show statistics
show_stats

echo
echo "Structure management completed successfully!"
