#!/bin/bash
# Script to manage data exports and file formats

# Exit on error
set -e

# Default values
DATA_DIR="data"
EXPORT_DIR="exports"
LOG_DIR="logs"
EXPORT_TYPES="tsv,csv,json,sdf,xlsx"
EXPORT_MODES="full,filtered,summary"
EXPORT_FIELDS="all,minimal,custom"
FILTER_FILE="filters.json"
BATCH_SIZE=1000
INCLUDE_PREDICTIONS=false
INCLUDE_COMMUNITY=false
INCLUDE_SAFETY=false
EXPORT=false
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
        --export-dir)
            EXPORT_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --export-types)
            EXPORT_TYPES="$2"
            shift 2
            ;;
        --export-modes)
            EXPORT_MODES="$2"
            shift 2
            ;;
        --export-fields)
            EXPORT_FIELDS="$2"
            shift 2
            ;;
        --filter-file)
            FILTER_FILE="$2"
            shift 2
            ;;
        --batch-size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --include-predictions)
            INCLUDE_PREDICTIONS=true
            shift
            ;;
        --include-community)
            INCLUDE_COMMUNITY=true
            shift
            ;;
        --include-safety)
            INCLUDE_SAFETY=true
            shift
            ;;
        --export)
            EXPORT=true
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
    
    # Export directories
    for type in ${EXPORT_TYPES//,/ }; do
        for mode in ${EXPORT_MODES//,/ }; do
            for fields in ${EXPORT_FIELDS//,/ }; do
                mkdir -p "$EXPORT_DIR/$type/$mode/$fields"/{data,reports}
            done
        done
    done
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to export data
export_data() {
    echo "Exporting data..."
    
    for type in ${EXPORT_TYPES//,/ }; do
        for mode in ${EXPORT_MODES//,/ }; do
            for fields in ${EXPORT_FIELDS//,/ }; do
                echo "Performing $mode export with $fields fields in $type format..."
                
                # Build command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli export-data"
                CMD="$CMD --data-dir $DATA_DIR"
                CMD="$CMD --export-dir $EXPORT_DIR/$type/$mode/$fields"
                CMD="$CMD --export-type $type"
                CMD="$CMD --mode $mode"
                CMD="$CMD --fields $fields"
                CMD="$CMD --batch-size $BATCH_SIZE"
                CMD="$CMD --log-dir $LOG_DIR"
                
                if [ -n "$FILTER_FILE" ]; then
                    CMD="$CMD --filter-file $FILTER_FILE"
                fi
                
                if [ "$INCLUDE_PREDICTIONS" = true ]; then
                    CMD="$CMD --include-predictions"
                fi
                
                if [ "$INCLUDE_COMMUNITY" = true ]; then
                    CMD="$CMD --include-community"
                fi
                
                if [ "$INCLUDE_SAFETY" = true ]; then
                    CMD="$CMD --include-safety"
                fi
                
                # Run command
                echo "Running: $CMD"
                $CMD || true  # Continue even if export fails
                
                # Generate export report
                echo "Generating export report..."
                python -m binding_data_processor.processors.psychopharm.predictors.cli generate-export-report \
                    --export-dir "$EXPORT_DIR/$type/$mode/$fields" \
                    --output "$EXPORT_DIR/$type/$mode/$fields/reports/report.html"
            done
        done
    done
}

# Function to validate exports
validate_exports() {
    echo "Validating exports..."
    
    for type in ${EXPORT_TYPES//,/ }; do
        for mode in ${EXPORT_MODES//,/ }; do
            for fields in ${EXPORT_FIELDS//,/ }; do
                echo "Validating $mode export with $fields fields in $type format..."
                
                # Build command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli validate-export"
                CMD="$CMD --export-dir $EXPORT_DIR/$type/$mode/$fields"
                CMD="$CMD --export-type $type"
                CMD="$CMD --mode $mode"
                CMD="$CMD --fields $fields"
                CMD="$CMD --log-dir $LOG_DIR"
                
                # Run command
                echo "Running: $CMD"
                $CMD || true  # Continue even if validation fails
                
                # Generate validation report
                echo "Generating validation report..."
                python -m binding_data_processor.processors.psychopharm.predictors.cli generate-validation-report \
                    --export-dir "$EXPORT_DIR/$type/$mode/$fields" \
                    --output "$EXPORT_DIR/$type/$mode/$fields/reports/validation.html"
            done
        done
    done
}

# Function to clean exports
clean_exports() {
    echo "Cleaning exports..."
    
    # Clean export directories
    rm -rf "$EXPORT_DIR"/*
}

# Function to backup exports
backup_exports() {
    echo "Backing up exports..."
    
    # Create backup directory
    BACKUP_DIR="backups/exports_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$EXPORT_DIR" "$BACKUP_DIR/"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore exports
restore_exports() {
    echo "Restoring exports..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/exports_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$EXPORT_DIR"
    fi
    
    cp -r "$BACKUP_DIR/exports" ./
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Exports restored from: $BACKUP_DIR"
}

# Function to show export statistics
show_stats() {
    echo "Export statistics:"
    echo
    
    echo "Export files:"
    for type in ${EXPORT_TYPES//,/ }; do
        echo "$type exports:"
        for mode in ${EXPORT_MODES//,/ }; do
            for fields in ${EXPORT_FIELDS//,/ }; do
                echo "  $mode mode ($fields fields):"
                
                # Data files
                echo "    Data:"
                echo "      Files: $(find "$EXPORT_DIR/$type/$mode/$fields/data" -type f | wc -l) files"
                echo "      Size: $(du -sh "$EXPORT_DIR/$type/$mode/$fields/data" | cut -f1)"
                
                # Count entries by type
                case "$type" in
                    tsv|csv)
                        if [ -f "$EXPORT_DIR/$type/$mode/$fields/data/compounds.$type" ]; then
                            echo "      Rows: $(($(wc -l < "$EXPORT_DIR/$type/$mode/$fields/data/compounds.$type") - 1))"
                            echo "      Columns: $(head -n1 "$EXPORT_DIR/$type/$mode/$fields/data/compounds.$type" | tr '\t,' '\n' | wc -l)"
                        fi
                        ;;
                    json)
                        if [ -f "$EXPORT_DIR/$type/$mode/$fields/data/compounds.json" ]; then
                            echo "      Compounds: $(jq '.compounds | length' "$EXPORT_DIR/$type/$mode/$fields/data/compounds.json")"
                            echo "      Fields: $(jq '.fields | length' "$EXPORT_DIR/$type/$mode/$fields/data/compounds.json")"
                        fi
                        ;;
                    sdf)
                        if [ -f "$EXPORT_DIR/$type/$mode/$fields/data/compounds.sdf" ]; then
                            echo "      Molecules: $(grep -c '$$$$' "$EXPORT_DIR/$type/$mode/$fields/data/compounds.sdf")"
                            echo "      Properties: $(grep -c '<' "$EXPORT_DIR/$type/$mode/$fields/data/compounds.sdf")"
                        fi
                        ;;
                    xlsx)
                        if [ -f "$EXPORT_DIR/$type/$mode/$fields/data/compounds.xlsx" ]; then
                            echo "      Sheets: $(python -c "import openpyxl; print(len(openpyxl.load_workbook('$EXPORT_DIR/$type/$mode/$fields/data/compounds.xlsx').sheetnames))")"
                            echo "      Rows: $(python -c "import openpyxl; print(openpyxl.load_workbook('$EXPORT_DIR/$type/$mode/$fields/data/compounds.xlsx').active.max_row - 1)")"
                        fi
                        ;;
                esac
                
                # Show reports
                if [ -f "$EXPORT_DIR/$type/$mode/$fields/reports/report.html" ]; then
                    echo "    Export report: $EXPORT_DIR/$type/$mode/$fields/reports/report.html"
                fi
                if [ -f "$EXPORT_DIR/$type/$mode/$fields/reports/validation.html" ]; then
                    echo "    Validation report: $EXPORT_DIR/$type/$mode/$fields/reports/validation.html"
                fi
                echo
            done
        done
    done
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing exports..."

# Create directory structure
create_dirs

# Export data if requested
if [ "$EXPORT" = true ]; then
    export_data
fi

# Validate exports if requested
if [ "$VALIDATE" = true ]; then
    validate_exports
fi

# Clean exports if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean exports? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_exports
    fi
fi

# Backup exports if requested
if [ "$BACKUP" = true ]; then
    backup_exports
fi

# Restore exports if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore exports? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_exports
    fi
fi

# Show statistics
show_stats

echo
echo "Export management completed successfully!"
