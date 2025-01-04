#!/bin/bash
# Script to manage data analysis and insights

# Exit on error
set -e

# Default values
DATA_DIR="data"
ANALYSIS_DIR="analysis"
LOG_DIR="logs"
ANALYSIS_TYPES="structure,activity,predictions,community,safety"
ANALYSIS_METHODS="statistics,clustering,similarity,trends,patterns"
MIN_CONFIDENCE=0.7
ANALYZE=false
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
        --analysis-dir)
            ANALYSIS_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --analysis-types)
            ANALYSIS_TYPES="$2"
            shift 2
            ;;
        --analysis-methods)
            ANALYSIS_METHODS="$2"
            shift 2
            ;;
        --min-confidence)
            MIN_CONFIDENCE="$2"
            shift 2
            ;;
        --analyze)
            ANALYZE=true
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
    
    # Analysis directories
    for type in ${ANALYSIS_TYPES//,/ }; do
        for method in ${ANALYSIS_METHODS//,/ }; do
            mkdir -p "$ANALYSIS_DIR/$type/$method"/{data,plots,reports}
        done
    done
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to run analysis
run_analysis() {
    echo "Running analysis..."
    
    for type in ${ANALYSIS_TYPES//,/ }; do
        for method in ${ANALYSIS_METHODS//,/ }; do
            echo "Running $method analysis on $type data..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli analyze-data"
            CMD="$CMD --data-dir $DATA_DIR"
            CMD="$CMD --analysis-dir $ANALYSIS_DIR/$type/$method"
            CMD="$CMD --analysis-type $type"
            CMD="$CMD --analysis-method $method"
            CMD="$CMD --min-confidence $MIN_CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if analysis fails
            
            # Generate analysis report
            echo "Generating analysis report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-analysis-report \
                --analysis-dir "$ANALYSIS_DIR/$type/$method" \
                --output "$ANALYSIS_DIR/$type/$method/reports/report.html"
        done
    done
}

# Function to clean analysis
clean_analysis() {
    echo "Cleaning analysis..."
    
    # Clean analysis directories
    rm -rf "$ANALYSIS_DIR"/*
}

# Function to backup analysis
backup_analysis() {
    echo "Backing up analysis..."
    
    # Create backup directory
    BACKUP_DIR="backups/analysis_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$ANALYSIS_DIR" "$BACKUP_DIR/"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore analysis
restore_analysis() {
    echo "Restoring analysis..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/analysis_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$ANALYSIS_DIR"
    fi
    
    cp -r "$BACKUP_DIR/analysis" ./
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Analysis restored from: $BACKUP_DIR"
}

# Function to show analysis statistics
show_stats() {
    echo "Analysis statistics:"
    echo
    
    echo "Analysis files:"
    for type in ${ANALYSIS_TYPES//,/ }; do
        echo "$type analysis:"
        for method in ${ANALYSIS_METHODS//,/ }; do
            echo "  $method analysis:"
            
            # Data files
            echo "    Data:"
            echo "      Files: $(find "$ANALYSIS_DIR/$type/$method/data" -type f | wc -l) files"
            echo "      Size: $(du -sh "$ANALYSIS_DIR/$type/$method/data" | cut -f1)"
            
            # Plot files
            echo "    Plots:"
            echo "      Files: $(find "$ANALYSIS_DIR/$type/$method/plots" -type f | wc -l) files"
            echo "      Size: $(du -sh "$ANALYSIS_DIR/$type/$method/plots" | cut -f1)"
            
            # Count analysis results
            case "$method" in
                statistics)
                    echo "    Statistics:"
                    if [ -f "$ANALYSIS_DIR/$type/$method/data/stats.json" ]; then
                        echo "      Metrics: $(jq '.metrics | length' "$ANALYSIS_DIR/$type/$method/data/stats.json")"
                        echo "      Tests: $(jq '.tests | length' "$ANALYSIS_DIR/$type/$method/data/stats.json")"
                    fi
                    ;;
                clustering)
                    echo "    Clustering:"
                    if [ -f "$ANALYSIS_DIR/$type/$method/data/clusters.json" ]; then
                        echo "      Clusters: $(jq '.clusters | length' "$ANALYSIS_DIR/$type/$method/data/clusters.json")"
                        echo "      Members: $(jq '.clusters[].members | length' "$ANALYSIS_DIR/$type/$method/data/clusters.json" | awk '{s+=$1} END {print s}')"
                    fi
                    ;;
                similarity)
                    echo "    Similarity:"
                    if [ -f "$ANALYSIS_DIR/$type/$method/data/similarity.json" ]; then
                        echo "      Pairs: $(jq '.pairs | length' "$ANALYSIS_DIR/$type/$method/data/similarity.json")"
                        echo "      Threshold: $(jq '.threshold' "$ANALYSIS_DIR/$type/$method/data/similarity.json")"
                    fi
                    ;;
                trends)
                    echo "    Trends:"
                    if [ -f "$ANALYSIS_DIR/$type/$method/data/trends.json" ]; then
                        echo "      Patterns: $(jq '.patterns | length' "$ANALYSIS_DIR/$type/$method/data/trends.json")"
                        echo "      Timespan: $(jq '.timespan' "$ANALYSIS_DIR/$type/$method/data/trends.json")"
                    fi
                    ;;
                patterns)
                    echo "    Patterns:"
                    if [ -f "$ANALYSIS_DIR/$type/$method/data/patterns.json" ]; then
                        echo "      Rules: $(jq '.rules | length' "$ANALYSIS_DIR/$type/$method/data/patterns.json")"
                        echo "      Support: $(jq '.min_support' "$ANALYSIS_DIR/$type/$method/data/patterns.json")"
                    fi
                    ;;
            esac
            
            # Show report
            if [ -f "$ANALYSIS_DIR/$type/$method/reports/report.html" ]; then
                echo "    Report: $ANALYSIS_DIR/$type/$method/reports/report.html"
            fi
            echo
        done
    done
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing analysis..."

# Create directory structure
create_dirs

# Run analysis if requested
if [ "$ANALYZE" = true ]; then
    run_analysis
fi

# Clean analysis if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean analysis? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_analysis
    fi
fi

# Backup analysis if requested
if [ "$BACKUP" = true ]; then
    backup_analysis
fi

# Restore analysis if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore analysis? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_analysis
    fi
fi

# Show statistics
show_stats

echo
echo "Analysis management completed successfully!"
