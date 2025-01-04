#!/bin/bash
# Script to manage regulatory data and legal status

# Exit on error
set -e

# Default values
DATA_DIR="data"
REGULATORY_DIR="regulatory"
LOG_DIR="logs"
JURISDICTION_TYPES="international,national,state,local"
STATUS_TYPES="controlled,scheduled,restricted,approved,investigational"
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
        --regulatory-dir)
            REGULATORY_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --jurisdiction-types)
            JURISDICTION_TYPES="$2"
            shift 2
            ;;
        --status-types)
            STATUS_TYPES="$2"
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
    
    # Regulatory directories
    for jurisdiction in ${JURISDICTION_TYPES//,/ }; do
        for status in ${STATUS_TYPES//,/ }; do
            mkdir -p "$REGULATORY_DIR/$jurisdiction/$status"/{data,updates,reports}
        done
    done
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to harvest regulatory data
harvest_data() {
    echo "Harvesting regulatory data..."
    
    for jurisdiction in ${JURISDICTION_TYPES//,/ }; do
        for status in ${STATUS_TYPES//,/ }; do
            echo "Harvesting $status data for $jurisdiction jurisdiction..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli harvest-regulatory"
            CMD="$CMD --data-dir $DATA_DIR"
            CMD="$CMD --regulatory-dir $REGULATORY_DIR/$jurisdiction/$status"
            CMD="$CMD --jurisdiction $jurisdiction"
            CMD="$CMD --status $status"
            CMD="$CMD --min-confidence $MIN_CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if harvesting fails
            
            # Generate harvest report
            echo "Generating harvest report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-regulatory-report \
                --regulatory-dir "$REGULATORY_DIR/$jurisdiction/$status" \
                --output "$REGULATORY_DIR/$jurisdiction/$status/reports/report.html"
        done
    done
}

# Function to validate regulatory data
validate_data() {
    echo "Validating regulatory data..."
    
    for jurisdiction in ${JURISDICTION_TYPES//,/ }; do
        for status in ${STATUS_TYPES//,/ }; do
            echo "Validating $status data for $jurisdiction jurisdiction..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli validate-regulatory"
            CMD="$CMD --regulatory-dir $REGULATORY_DIR/$jurisdiction/$status"
            CMD="$CMD --jurisdiction $jurisdiction"
            CMD="$CMD --status $status"
            CMD="$CMD --min-confidence $MIN_CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if validation fails
            
            # Generate validation report
            echo "Generating validation report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-validation-report \
                --regulatory-dir "$REGULATORY_DIR/$jurisdiction/$status" \
                --output "$REGULATORY_DIR/$jurisdiction/$status/reports/validation.html"
        done
    done
}

# Function to clean regulatory data
clean_regulatory() {
    echo "Cleaning regulatory data..."
    
    # Clean regulatory directories
    rm -rf "$REGULATORY_DIR"/*
}

# Function to backup regulatory data
backup_regulatory() {
    echo "Backing up regulatory data..."
    
    # Create backup directory
    BACKUP_DIR="backups/regulatory_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$REGULATORY_DIR" "$BACKUP_DIR/"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore regulatory data
restore_regulatory() {
    echo "Restoring regulatory data..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/regulatory_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$REGULATORY_DIR"
    fi
    
    cp -r "$BACKUP_DIR/regulatory" ./
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Regulatory data restored from: $BACKUP_DIR"
}

# Function to show regulatory statistics
show_stats() {
    echo "Regulatory statistics:"
    echo
    
    echo "Regulatory files:"
    for jurisdiction in ${JURISDICTION_TYPES//,/ }; do
        echo "$jurisdiction jurisdiction:"
        for status in ${STATUS_TYPES//,/ }; do
            echo "  $status status:"
            
            # Data files
            echo "    Data:"
            echo "      Files: $(find "$REGULATORY_DIR/$jurisdiction/$status/data" -type f | wc -l) files"
            echo "      Size: $(du -sh "$REGULATORY_DIR/$jurisdiction/$status/data" | cut -f1)"
            
            # Count entries by status
            case "$status" in
                controlled)
                    if [ -f "$REGULATORY_DIR/$jurisdiction/$status/data/controlled.json" ]; then
                        echo "      Substances: $(jq '.substances | length' "$REGULATORY_DIR/$jurisdiction/$status/data/controlled.json")"
                        echo "      Schedules: $(jq '.schedules | length' "$REGULATORY_DIR/$jurisdiction/$status/data/controlled.json")"
                    fi
                    ;;
                scheduled)
                    if [ -f "$REGULATORY_DIR/$jurisdiction/$status/data/scheduled.json" ]; then
                        echo "      Substances: $(jq '.substances | length' "$REGULATORY_DIR/$jurisdiction/$status/data/scheduled.json")"
                        echo "      Categories: $(jq '.categories | length' "$REGULATORY_DIR/$jurisdiction/$status/data/scheduled.json")"
                    fi
                    ;;
                restricted)
                    if [ -f "$REGULATORY_DIR/$jurisdiction/$status/data/restricted.json" ]; then
                        echo "      Substances: $(jq '.substances | length' "$REGULATORY_DIR/$jurisdiction/$status/data/restricted.json")"
                        echo "      Restrictions: $(jq '.restrictions | length' "$REGULATORY_DIR/$jurisdiction/$status/data/restricted.json")"
                    fi
                    ;;
                approved)
                    if [ -f "$REGULATORY_DIR/$jurisdiction/$status/data/approved.json" ]; then
                        echo "      Substances: $(jq '.substances | length' "$REGULATORY_DIR/$jurisdiction/$status/data/approved.json")"
                        echo "      Indications: $(jq '.indications | length' "$REGULATORY_DIR/$jurisdiction/$status/data/approved.json")"
                    fi
                    ;;
                investigational)
                    if [ -f "$REGULATORY_DIR/$jurisdiction/$status/data/investigational.json" ]; then
                        echo "      Substances: $(jq '.substances | length' "$REGULATORY_DIR/$jurisdiction/$status/data/investigational.json")"
                        echo "      Trials: $(jq '.trials | length' "$REGULATORY_DIR/$jurisdiction/$status/data/investigational.json")"
                    fi
                    ;;
            esac
            
            # Updates
            echo "    Updates:"
            echo "      Files: $(find "$REGULATORY_DIR/$jurisdiction/$status/updates" -type f | wc -l) files"
            echo "      Size: $(du -sh "$REGULATORY_DIR/$jurisdiction/$status/updates" | cut -f1)"
            
            # Show reports
            if [ -f "$REGULATORY_DIR/$jurisdiction/$status/reports/report.html" ]; then
                echo "    Harvest report: $REGULATORY_DIR/$jurisdiction/$status/reports/report.html"
            fi
            if [ -f "$REGULATORY_DIR/$jurisdiction/$status/reports/validation.html" ]; then
                echo "    Validation report: $REGULATORY_DIR/$jurisdiction/$status/reports/validation.html"
            fi
            echo
        done
    done
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing regulatory data..."

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
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean regulatory data? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_regulatory
    fi
fi

# Backup data if requested
if [ "$BACKUP" = true ]; then
    backup_regulatory
fi

# Restore data if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore regulatory data? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_regulatory
    fi
fi

# Show statistics
show_stats

echo
echo "Regulatory management completed successfully!"
