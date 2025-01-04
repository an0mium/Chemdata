#!/bin/bash
# Script to manage community data and user reports

# Exit on error
set -e

# Default values
DATA_DIR="data"
COMMUNITY_DIR="community"
LOG_DIR="logs"
SOURCE_TYPES="reddit,twitter,discord,bluesky,erowid,psychonautwiki,tripsit"
DATA_TYPES="reports,discussions,reviews,alerts,questions"
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
        --community-dir)
            COMMUNITY_DIR="$2"
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
        --data-types)
            DATA_TYPES="$2"
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
    
    # Community directories
    for source in ${SOURCE_TYPES//,/ }; do
        for type in ${DATA_TYPES//,/ }; do
            mkdir -p "$COMMUNITY_DIR/$source/$type"/{raw,processed,validated}
        done
    done
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to harvest community data
harvest_data() {
    echo "Harvesting community data..."
    
    for source in ${SOURCE_TYPES//,/ }; do
        for type in ${DATA_TYPES//,/ }; do
            echo "Harvesting $type data from $source..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli harvest-community"
            CMD="$CMD --data-dir $DATA_DIR"
            CMD="$CMD --community-dir $COMMUNITY_DIR/$source/$type"
            CMD="$CMD --source $source"
            CMD="$CMD --type $type"
            CMD="$CMD --min-confidence $MIN_CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if harvesting fails
            
            # Generate harvest report
            echo "Generating harvest report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-harvest-report \
                --community-dir "$COMMUNITY_DIR/$source/$type" \
                --output "$COMMUNITY_DIR/$source/$type/raw/report.html"
        done
    done
}

# Function to validate community data
validate_data() {
    echo "Validating community data..."
    
    for source in ${SOURCE_TYPES//,/ }; do
        for type in ${DATA_TYPES//,/ }; do
            echo "Validating $type data from $source..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli validate-community"
            CMD="$CMD --community-dir $COMMUNITY_DIR/$source/$type"
            CMD="$CMD --source $source"
            CMD="$CMD --type $type"
            CMD="$CMD --min-confidence $MIN_CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if validation fails
            
            # Generate validation report
            echo "Generating validation report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-validation-report \
                --community-dir "$COMMUNITY_DIR/$source/$type" \
                --output "$COMMUNITY_DIR/$source/$type/validated/report.html"
        done
    done
}

# Function to clean community data
clean_community() {
    echo "Cleaning community data..."
    
    # Clean community directories
    rm -rf "$COMMUNITY_DIR"/*
}

# Function to backup community data
backup_community() {
    echo "Backing up community data..."
    
    # Create backup directory
    BACKUP_DIR="backups/community_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$COMMUNITY_DIR" "$BACKUP_DIR/"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore community data
restore_community() {
    echo "Restoring community data..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/community_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$COMMUNITY_DIR"
    fi
    
    cp -r "$BACKUP_DIR/community" ./
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Community data restored from: $BACKUP_DIR"
}

# Function to show community statistics
show_stats() {
    echo "Community statistics:"
    echo
    
    echo "Community files:"
    for source in ${SOURCE_TYPES//,/ }; do
        echo "$source data:"
        for type in ${DATA_TYPES//,/ }; do
            echo "  $type data:"
            
            # Raw data
            echo "    Raw data:"
            echo "      Files: $(find "$COMMUNITY_DIR/$source/$type/raw" -type f | wc -l) files"
            echo "      Size: $(du -sh "$COMMUNITY_DIR/$source/$type/raw" | cut -f1)"
            
            # Count entries by type
            case "$type" in
                reports)
                    if [ -f "$COMMUNITY_DIR/$source/$type/raw/reports.json" ]; then
                        echo "      Reports: $(jq '.reports | length' "$COMMUNITY_DIR/$source/$type/raw/reports.json")"
                        echo "      Users: $(jq '.users | length' "$COMMUNITY_DIR/$source/$type/raw/reports.json")"
                    fi
                    ;;
                discussions)
                    if [ -f "$COMMUNITY_DIR/$source/$type/raw/discussions.json" ]; then
                        echo "      Threads: $(jq '.threads | length' "$COMMUNITY_DIR/$source/$type/raw/discussions.json")"
                        echo "      Posts: $(jq '.posts | length' "$COMMUNITY_DIR/$source/$type/raw/discussions.json")"
                    fi
                    ;;
                reviews)
                    if [ -f "$COMMUNITY_DIR/$source/$type/raw/reviews.json" ]; then
                        echo "      Reviews: $(jq '.reviews | length' "$COMMUNITY_DIR/$source/$type/raw/reviews.json")"
                        echo "      Ratings: $(jq '.ratings | length' "$COMMUNITY_DIR/$source/$type/raw/reviews.json")"
                    fi
                    ;;
                alerts)
                    if [ -f "$COMMUNITY_DIR/$source/$type/raw/alerts.json" ]; then
                        echo "      Alerts: $(jq '.alerts | length' "$COMMUNITY_DIR/$source/$type/raw/alerts.json")"
                        echo "      Severity: $(jq '.severity | length' "$COMMUNITY_DIR/$source/$type/raw/alerts.json")"
                    fi
                    ;;
                questions)
                    if [ -f "$COMMUNITY_DIR/$source/$type/raw/questions.json" ]; then
                        echo "      Questions: $(jq '.questions | length' "$COMMUNITY_DIR/$source/$type/raw/questions.json")"
                        echo "      Answers: $(jq '.answers | length' "$COMMUNITY_DIR/$source/$type/raw/questions.json")"
                    fi
                    ;;
            esac
            
            # Processed data
            echo "    Processed data:"
            echo "      Files: $(find "$COMMUNITY_DIR/$source/$type/processed" -type f | wc -l) files"
            echo "      Size: $(du -sh "$COMMUNITY_DIR/$source/$type/processed" | cut -f1)"
            
            # Validated data
            echo "    Validated data:"
            echo "      Files: $(find "$COMMUNITY_DIR/$source/$type/validated" -type f | wc -l) files"
            echo "      Size: $(du -sh "$COMMUNITY_DIR/$source/$type/validated" | cut -f1)"
            
            # Show reports
            if [ -f "$COMMUNITY_DIR/$source/$type/raw/report.html" ]; then
                echo "    Harvest report: $COMMUNITY_DIR/$source/$type/raw/report.html"
            fi
            if [ -f "$COMMUNITY_DIR/$source/$type/validated/report.html" ]; then
                echo "    Validation report: $COMMUNITY_DIR/$source/$type/validated/report.html"
            fi
            echo
        done
    done
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing community data..."

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
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean community data? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_community
    fi
fi

# Backup data if requested
if [ "$BACKUP" = true ]; then
    backup_community
fi

# Restore data if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore community data? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_community
    fi
fi

# Show statistics
show_stats

echo
echo "Community management completed successfully!"
