#!/bin/bash
# Script to manage chemical property data and calculations

# Exit on error
set -e

# Default values
DATA_DIR="data"
PROPERTY_DIR="properties"
LOG_DIR="logs"
PROPERTY_TYPES="physical,chemical,quantum,descriptors,fingerprints"
CALCULATION_TYPES="predict,compute,estimate,measure,validate"
MIN_CONFIDENCE=0.7
CALCULATE=false
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
        --property-dir)
            PROPERTY_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --property-types)
            PROPERTY_TYPES="$2"
            shift 2
            ;;
        --calculation-types)
            CALCULATION_TYPES="$2"
            shift 2
            ;;
        --min-confidence)
            MIN_CONFIDENCE="$2"
            shift 2
            ;;
        --calculate)
            CALCULATE=true
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
    
    # Property directories
    for property in ${PROPERTY_TYPES//,/ }; do
        for calculation in ${CALCULATION_TYPES//,/ }; do
            mkdir -p "$PROPERTY_DIR/$property/$calculation"/{input,output,reports}
        done
    done
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to calculate properties
calculate_properties() {
    echo "Calculating properties..."
    
    for property in ${PROPERTY_TYPES//,/ }; do
        for calculation in ${CALCULATION_TYPES//,/ }; do
            echo "Performing $calculation for $property properties..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli calculate-properties"
            CMD="$CMD --data-dir $DATA_DIR"
            CMD="$CMD --property-dir $PROPERTY_DIR/$property/$calculation"
            CMD="$CMD --property $property"
            CMD="$CMD --calculation $calculation"
            CMD="$CMD --min-confidence $MIN_CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if calculation fails
            
            # Generate calculation report
            echo "Generating calculation report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-property-report \
                --property-dir "$PROPERTY_DIR/$property/$calculation" \
                --output "$PROPERTY_DIR/$property/$calculation/reports/report.html"
        done
    done
}

# Function to validate properties
validate_properties() {
    echo "Validating properties..."
    
    for property in ${PROPERTY_TYPES//,/ }; do
        for calculation in ${CALCULATION_TYPES//,/ }; do
            echo "Validating $calculation results for $property properties..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli validate-properties"
            CMD="$CMD --property-dir $PROPERTY_DIR/$property/$calculation"
            CMD="$CMD --property $property"
            CMD="$CMD --calculation $calculation"
            CMD="$CMD --min-confidence $MIN_CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if validation fails
            
            # Generate validation report
            echo "Generating validation report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-validation-report \
                --property-dir "$PROPERTY_DIR/$property/$calculation" \
                --output "$PROPERTY_DIR/$property/$calculation/reports/validation.html"
        done
    done
}

# Function to clean properties
clean_properties() {
    echo "Cleaning properties..."
    
    # Clean property directories
    rm -rf "$PROPERTY_DIR"/*
}

# Function to backup properties
backup_properties() {
    echo "Backing up properties..."
    
    # Create backup directory
    BACKUP_DIR="backups/properties_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$PROPERTY_DIR" "$BACKUP_DIR/"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore properties
restore_properties() {
    echo "Restoring properties..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/properties_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$PROPERTY_DIR"
    fi
    
    cp -r "$BACKUP_DIR/properties" ./
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Properties restored from: $BACKUP_DIR"
}

# Function to show property statistics
show_stats() {
    echo "Property statistics:"
    echo
    
    echo "Property files:"
    for property in ${PROPERTY_TYPES//,/ }; do
        echo "$property properties:"
        for calculation in ${CALCULATION_TYPES//,/ }; do
            echo "  $calculation calculation:"
            
            # Input files
            echo "    Input:"
            echo "      Files: $(find "$PROPERTY_DIR/$property/$calculation/input" -type f | wc -l) files"
            echo "      Size: $(du -sh "$PROPERTY_DIR/$property/$calculation/input" | cut -f1)"
            
            # Count properties by type
            case "$property" in
                physical)
                    if [ -f "$PROPERTY_DIR/$property/$calculation/output/physical.json" ]; then
                        echo "      Properties: $(jq '.properties | length' "$PROPERTY_DIR/$property/$calculation/output/physical.json")"
                        echo "      Compounds: $(jq '.compounds | length' "$PROPERTY_DIR/$property/$calculation/output/physical.json")"
                    fi
                    ;;
                chemical)
                    if [ -f "$PROPERTY_DIR/$property/$calculation/output/chemical.json" ]; then
                        echo "      Properties: $(jq '.properties | length' "$PROPERTY_DIR/$property/$calculation/output/chemical.json")"
                        echo "      Compounds: $(jq '.compounds | length' "$PROPERTY_DIR/$property/$calculation/output/chemical.json")"
                    fi
                    ;;
                quantum)
                    if [ -f "$PROPERTY_DIR/$property/$calculation/output/quantum.json" ]; then
                        echo "      Properties: $(jq '.properties | length' "$PROPERTY_DIR/$property/$calculation/output/quantum.json")"
                        echo "      Compounds: $(jq '.compounds | length' "$PROPERTY_DIR/$property/$calculation/output/quantum.json")"
                    fi
                    ;;
                descriptors)
                    if [ -f "$PROPERTY_DIR/$property/$calculation/output/descriptors.json" ]; then
                        echo "      Descriptors: $(jq '.descriptors | length' "$PROPERTY_DIR/$property/$calculation/output/descriptors.json")"
                        echo "      Compounds: $(jq '.compounds | length' "$PROPERTY_DIR/$property/$calculation/output/descriptors.json")"
                    fi
                    ;;
                fingerprints)
                    if [ -f "$PROPERTY_DIR/$property/$calculation/output/fingerprints.json" ]; then
                        echo "      Types: $(jq '.types | length' "$PROPERTY_DIR/$property/$calculation/output/fingerprints.json")"
                        echo "      Compounds: $(jq '.compounds | length' "$PROPERTY_DIR/$property/$calculation/output/fingerprints.json")"
                    fi
                    ;;
            esac
            
            # Output files
            echo "    Output:"
            echo "      Files: $(find "$PROPERTY_DIR/$property/$calculation/output" -type f | wc -l) files"
            echo "      Size: $(du -sh "$PROPERTY_DIR/$property/$calculation/output" | cut -f1)"
            
            # Show reports
            if [ -f "$PROPERTY_DIR/$property/$calculation/reports/report.html" ]; then
                echo "    Calculation report: $PROPERTY_DIR/$property/$calculation/reports/report.html"
            fi
            if [ -f "$PROPERTY_DIR/$property/$calculation/reports/validation.html" ]; then
                echo "    Validation report: $PROPERTY_DIR/$property/$calculation/reports/validation.html"
            fi
            echo
        done
    done
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing properties..."

# Create directory structure
create_dirs

# Calculate properties if requested
if [ "$CALCULATE" = true ]; then
    calculate_properties
fi

# Validate properties if requested
if [ "$VALIDATE" = true ]; then
    validate_properties
fi

# Clean properties if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean properties? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_properties
    fi
fi

# Backup properties if requested
if [ "$BACKUP" = true ]; then
    backup_properties
fi

# Restore properties if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore properties? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_properties
    fi
fi

# Show statistics
show_stats

echo
echo "Property management completed successfully!"
