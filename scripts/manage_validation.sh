#!/bin/bash
# Script to manage data validation and quality checks

# Exit on error
set -e

# Default values
DATA_DIR="data"
VALIDATION_DIR="validation"
LOG_DIR="logs"
VALIDATION_TYPES="structure,activity,predictions,community,safety"
VALIDATION_LEVELS="error,warning,info"
MIN_CONFIDENCE=0.7
STRICT=false
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
        --validation-dir)
            VALIDATION_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --validation-types)
            VALIDATION_TYPES="$2"
            shift 2
            ;;
        --validation-levels)
            VALIDATION_LEVELS="$2"
            shift 2
            ;;
        --min-confidence)
            MIN_CONFIDENCE="$2"
            shift 2
            ;;
        --strict)
            STRICT=true
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
    
    # Validation directories
    for type in ${VALIDATION_TYPES//,/ }; do
        mkdir -p "$VALIDATION_DIR/$type"/{reports,issues,fixes}
    done
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to run validation
run_validation() {
    echo "Running validation..."
    
    for type in ${VALIDATION_TYPES//,/ }; do
        echo "Validating $type data..."
        
        # Build command
        CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli validate-data"
        CMD="$CMD --data-dir $DATA_DIR"
        CMD="$CMD --validation-dir $VALIDATION_DIR/$type"
        CMD="$CMD --validation-type $type"
        CMD="$CMD --validation-levels $VALIDATION_LEVELS"
        CMD="$CMD --min-confidence $MIN_CONFIDENCE"
        CMD="$CMD --log-dir $LOG_DIR"
        
        if [ "$STRICT" = true ]; then
            CMD="$CMD --strict"
        fi
        
        # Run command
        echo "Running: $CMD"
        $CMD || true  # Continue even if validation fails
        
        # Generate validation report
        echo "Generating validation report..."
        python -m binding_data_processor.processors.psychopharm.predictors.cli generate-validation-report \
            --validation-dir "$VALIDATION_DIR/$type" \
            --output "$VALIDATION_DIR/$type/reports/report.html"
    done
}

# Function to clean validation
clean_validation() {
    echo "Cleaning validation..."
    
    # Clean validation directories
    rm -rf "$VALIDATION_DIR"/*
}

# Function to backup validation
backup_validation() {
    echo "Backing up validation..."
    
    # Create backup directory
    BACKUP_DIR="backups/validation_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$VALIDATION_DIR" "$BACKUP_DIR/"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore validation
restore_validation() {
    echo "Restoring validation..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/validation_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$VALIDATION_DIR"
    fi
    
    cp -r "$BACKUP_DIR/validation" ./
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Validation restored from: $BACKUP_DIR"
}

# Function to show validation statistics
show_stats() {
    echo "Validation statistics:"
    echo
    
    echo "Validation files:"
    for type in ${VALIDATION_TYPES//,/ }; do
        echo "$type validation:"
        
        # Count issues by level
        for level in ${VALIDATION_LEVELS//,/ }; do
            echo "  $level issues:"
            echo "    Files: $(find "$VALIDATION_DIR/$type/issues" -name "*_$level.json" | wc -l) files"
            
            # Count total issues
            total_issues=0
            for file in "$VALIDATION_DIR/$type/issues"/*"_$level.json"; do
                if [ -f "$file" ]; then
                    issues=$(jq '.issues | length' "$file" 2>/dev/null || echo 0)
                    total_issues=$((total_issues + issues))
                fi
            done
            echo "    Total: $total_issues issues"
        done
        
        # Count fixes
        echo "  Fixes:"
        echo "    Files: $(find "$VALIDATION_DIR/$type/fixes" -type f | wc -l) files"
        echo "    Applied: $(find "$VALIDATION_DIR/$type/fixes" -name "*_applied.json" | wc -l) fixes"
        
        # Show report status
        if [ -f "$VALIDATION_DIR/$type/reports/report.html" ]; then
            echo "  Report: $VALIDATION_DIR/$type/reports/report.html"
        fi
        echo
    done
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing validation..."

# Create directory structure
create_dirs

# Run validation if requested
if [ "$VALIDATE" = true ]; then
    run_validation
fi

# Clean validation if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean validation? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_validation
    fi
fi

# Backup validation if requested
if [ "$BACKUP" = true ]; then
    backup_validation
fi

# Restore validation if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore validation? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_validation
    fi
fi

# Show statistics
show_stats

echo
echo "Validation management completed successfully!"
