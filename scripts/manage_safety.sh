#!/bin/bash
# Script to manage safety checks and risk assessments

# Exit on error
set -e

# Default values
DATA_DIR="data"
SAFETY_DIR="safety"
LOG_DIR="logs"
SAFETY_TYPES="toxicity,abuse,interactions,warnings,regulatory"
SAFETY_LEVELS="high,medium,low"
MIN_CONFIDENCE=0.7
CHECK=false
CLEAN=false
BACKUP=false
RESTORE=false
STRICT=false
FORCE=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --safety-dir)
            SAFETY_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --safety-types)
            SAFETY_TYPES="$2"
            shift 2
            ;;
        --safety-levels)
            SAFETY_LEVELS="$2"
            shift 2
            ;;
        --min-confidence)
            MIN_CONFIDENCE="$2"
            shift 2
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
        --strict)
            STRICT=true
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
    
    # Safety directories
    for type in ${SAFETY_TYPES//,/ }; do
        for level in ${SAFETY_LEVELS//,/ }; do
            mkdir -p "$SAFETY_DIR/$type/$level"/{checks,alerts,reports}
        done
    done
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to run safety checks
run_checks() {
    echo "Running safety checks..."
    
    for type in ${SAFETY_TYPES//,/ }; do
        for level in ${SAFETY_LEVELS//,/ }; do
            echo "Running $level $type safety checks..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli check-safety"
            CMD="$CMD --data-dir $DATA_DIR"
            CMD="$CMD --safety-dir $SAFETY_DIR/$type/$level"
            CMD="$CMD --safety-type $type"
            CMD="$CMD --safety-level $level"
            CMD="$CMD --min-confidence $MIN_CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            if [ "$STRICT" = true ]; then
                CMD="$CMD --strict"
            fi
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if checks fail
            
            # Generate safety report
            echo "Generating safety report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-safety-report \
                --safety-dir "$SAFETY_DIR/$type/$level" \
                --output "$SAFETY_DIR/$type/$level/reports/report.html"
        done
    done
}

# Function to clean safety checks
clean_safety() {
    echo "Cleaning safety checks..."
    
    # Clean safety directories
    rm -rf "$SAFETY_DIR"/*
}

# Function to backup safety checks
backup_safety() {
    echo "Backing up safety checks..."
    
    # Create backup directory
    BACKUP_DIR="backups/safety_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$SAFETY_DIR" "$BACKUP_DIR/"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore safety checks
restore_safety() {
    echo "Restoring safety checks..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/safety_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$SAFETY_DIR"
    fi
    
    cp -r "$BACKUP_DIR/safety" ./
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Safety checks restored from: $BACKUP_DIR"
}

# Function to show safety statistics
show_stats() {
    echo "Safety statistics:"
    echo
    
    echo "Safety files:"
    for type in ${SAFETY_TYPES//,/ }; do
        echo "$type safety checks:"
        for level in ${SAFETY_LEVELS//,/ }; do
            echo "  $level level:"
            
            # Check files
            echo "    Checks:"
            echo "      Files: $(find "$SAFETY_DIR/$type/$level/checks" -type f | wc -l) files"
            echo "      Size: $(du -sh "$SAFETY_DIR/$type/$level/checks" | cut -f1)"
            
            # Alert files
            echo "    Alerts:"
            echo "      Files: $(find "$SAFETY_DIR/$type/$level/alerts" -type f | wc -l) files"
            echo "      Size: $(du -sh "$SAFETY_DIR/$type/$level/alerts" | cut -f1)"
            
            # Count safety issues
            case "$type" in
                toxicity)
                    echo "    Toxicity:"
                    if [ -f "$SAFETY_DIR/$type/$level/checks/toxicity.json" ]; then
                        echo "      Risks: $(jq '.risks | length' "$SAFETY_DIR/$type/$level/checks/toxicity.json")"
                        echo "      Compounds: $(jq '.compounds | length' "$SAFETY_DIR/$type/$level/checks/toxicity.json")"
                    fi
                    ;;
                abuse)
                    echo "    Abuse potential:"
                    if [ -f "$SAFETY_DIR/$type/$level/checks/abuse.json" ]; then
                        echo "      Risks: $(jq '.risks | length' "$SAFETY_DIR/$type/$level/checks/abuse.json")"
                        echo "      Compounds: $(jq '.compounds | length' "$SAFETY_DIR/$type/$level/checks/abuse.json")"
                    fi
                    ;;
                interactions)
                    echo "    Interactions:"
                    if [ -f "$SAFETY_DIR/$type/$level/checks/interactions.json" ]; then
                        echo "      Pairs: $(jq '.pairs | length' "$SAFETY_DIR/$type/$level/checks/interactions.json")"
                        echo "      Severity: $(jq '.severity' "$SAFETY_DIR/$type/$level/checks/interactions.json")"
                    fi
                    ;;
                warnings)
                    echo "    Warnings:"
                    if [ -f "$SAFETY_DIR/$type/$level/checks/warnings.json" ]; then
                        echo "      Alerts: $(jq '.alerts | length' "$SAFETY_DIR/$type/$level/checks/warnings.json")"
                        echo "      Categories: $(jq '.categories | length' "$SAFETY_DIR/$type/$level/checks/warnings.json")"
                    fi
                    ;;
                regulatory)
                    echo "    Regulatory:"
                    if [ -f "$SAFETY_DIR/$type/$level/checks/regulatory.json" ]; then
                        echo "      Status: $(jq '.status | length' "$SAFETY_DIR/$type/$level/checks/regulatory.json")"
                        echo "      Jurisdictions: $(jq '.jurisdictions | length' "$SAFETY_DIR/$type/$level/checks/regulatory.json")"
                    fi
                    ;;
            esac
            
            # Show report
            if [ -f "$SAFETY_DIR/$type/$level/reports/report.html" ]; then
                echo "    Report: $SAFETY_DIR/$type/$level/reports/report.html"
            fi
            echo
        done
    done
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing safety checks..."

# Create directory structure
create_dirs

# Run safety checks if requested
if [ "$CHECK" = true ]; then
    run_checks
fi

# Clean safety checks if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean safety checks? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_safety
    fi
fi

# Backup safety checks if requested
if [ "$BACKUP" = true ]; then
    backup_safety
fi

# Restore safety checks if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore safety checks? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_safety
    fi
fi

# Show statistics
show_stats

echo
echo "Safety management completed successfully!"
