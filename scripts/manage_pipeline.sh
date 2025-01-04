#!/bin/bash
# Script to manage data processing pipeline and workflow

# Exit on error
set -e

# Default values
DATA_DIR="data"
PIPELINE_DIR="pipeline"
CACHE_DIR="cache/pipeline"
LOG_DIR="logs"
STAGE_TYPES="collect,process,analyze,predict,validate,export"
PIPELINE_MODES="bindingdb,community,literature,regulatory,development,production,testing"
BATCH_SIZE=100
MIN_CONFIDENCE=0.7
PARALLEL=4
TIMEOUT=3600
RETRIES=3
RUN=false
MONITOR=false
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
        --pipeline-dir)
            PIPELINE_DIR="$2"
            shift 2
            ;;
        --cache-dir)
            CACHE_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --stage-types)
            STAGE_TYPES="$2"
            shift 2
            ;;
        --pipeline-modes)
            PIPELINE_MODES="$2"
            shift 2
            ;;
        --batch-size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --min-confidence)
            MIN_CONFIDENCE="$2"
            shift 2
            ;;
        --parallel)
            PARALLEL="$2"
            shift 2
            ;;
        --timeout)
            TIMEOUT="$2"
            shift 2
            ;;
        --retries)
            RETRIES="$2"
            shift 2
            ;;
        --run)
            RUN=true
            shift
            ;;
        --monitor)
            MONITOR=true
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
    
    # Pipeline directories
    for mode in ${PIPELINE_MODES//,/ }; do
        for stage in ${STAGE_TYPES//,/ }; do
            mkdir -p "$PIPELINE_DIR/$mode/$stage"/{input,output,checkpoints,reports}
        done
    done
    
    # Cache directory
    mkdir -p "$CACHE_DIR"
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to run pipeline
run_pipeline() {
    echo "Running pipeline..."
    
    for mode in ${PIPELINE_MODES//,/ }; do
        echo "Running pipeline in $mode mode..."
        
        for stage in ${STAGE_TYPES//,/ }; do
            echo "Running $stage stage..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli run-pipeline"
            CMD="$CMD --data-dir $DATA_DIR"
            CMD="$CMD --pipeline-dir $PIPELINE_DIR/$mode/$stage"
            CMD="$CMD --cache-dir $CACHE_DIR"
            CMD="$CMD --mode $mode"
            CMD="$CMD --stage $stage"
            CMD="$CMD --batch-size $BATCH_SIZE"
            CMD="$CMD --min-confidence $MIN_CONFIDENCE"
            CMD="$CMD --parallel $PARALLEL"
            CMD="$CMD --timeout $TIMEOUT"
            CMD="$CMD --retries $RETRIES"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if stage fails
            
            # Generate stage report
            echo "Generating stage report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-pipeline-report \
                --pipeline-dir "$PIPELINE_DIR/$mode/$stage" \
                --output "$PIPELINE_DIR/$mode/$stage/reports/stage.html"
        done
        
        # Generate pipeline report
        echo "Generating pipeline report..."
        python -m binding_data_processor.processors.psychopharm.predictors.cli generate-pipeline-report \
            --pipeline-dir "$PIPELINE_DIR" \
            --mode "$mode" \
            --output "$PIPELINE_DIR/pipeline_$mode.html"
    done
}

# Function to monitor pipeline
monitor_pipeline() {
    echo "Monitoring pipeline..."
    
    # Build command
    CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli monitor-pipeline"
    CMD="$CMD --pipeline-dir $PIPELINE_DIR"
    CMD="$CMD --log-dir $LOG_DIR"
    
    # Run command
    echo "Running: $CMD"
    $CMD
}

# Function to validate pipeline
validate_pipeline() {
    echo "Validating pipeline..."
    
    for mode in ${PIPELINE_MODES//,/ }; do
        for stage in ${STAGE_TYPES//,/ }; do
            echo "Validating $stage stage in $mode mode..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli validate-pipeline"
            CMD="$CMD --pipeline-dir $PIPELINE_DIR/$mode/$stage"
            CMD="$CMD --mode $mode"
            CMD="$CMD --stage $stage"
            CMD="$CMD --min-confidence $MIN_CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if validation fails
            
            # Generate validation report
            echo "Generating validation report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-validation-report \
                --pipeline-dir "$PIPELINE_DIR/$mode/$stage" \
                --output "$PIPELINE_DIR/$mode/$stage/reports/validation.html"
        done
    done
}

# Function to clean pipeline
clean_pipeline() {
    echo "Cleaning pipeline..."
    
    # Clean pipeline directories
    rm -rf "$PIPELINE_DIR"/*
    
    # Clean cache
    rm -rf "$CACHE_DIR"/*
}

# Function to backup pipeline
backup_pipeline() {
    echo "Backing up pipeline..."
    
    # Create backup directory
    BACKUP_DIR="backups/pipeline_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$PIPELINE_DIR" "$BACKUP_DIR/"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore pipeline
restore_pipeline() {
    echo "Restoring pipeline..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/pipeline_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$PIPELINE_DIR"
    fi
    
    cp -r "$BACKUP_DIR/pipeline" ./
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Pipeline restored from: $BACKUP_DIR"
}

# Function to show pipeline statistics
show_stats() {
    echo "Pipeline statistics:"
    echo
    
    echo "Pipeline files:"
    for mode in ${PIPELINE_MODES//,/ }; do
        echo "$mode mode:"
        for stage in ${STAGE_TYPES//,/ }; do
            echo "  $stage stage:"
            
            # Input files
            echo "    Input:"
            echo "      Files: $(find "$PIPELINE_DIR/$mode/$stage/input" -type f | wc -l) files"
            echo "      Size: $(du -sh "$PIPELINE_DIR/$mode/$stage/input" | cut -f1)"
            
            # Output files
            echo "    Output:"
            echo "      Files: $(find "$PIPELINE_DIR/$mode/$stage/output" -type f | wc -l) files"
            echo "      Size: $(du -sh "$PIPELINE_DIR/$mode/$stage/output" | cut -f1)"
            
            # Checkpoint files
            echo "    Checkpoints:"
            echo "      Files: $(find "$PIPELINE_DIR/$mode/$stage/checkpoints" -type f | wc -l) files"
            echo "      Size: $(du -sh "$PIPELINE_DIR/$mode/$stage/checkpoints" | cut -f1)"
            
            # Count entries by stage
            case "$stage" in
                collect)
                    if [ -f "$PIPELINE_DIR/$mode/$stage/output/collection.json" ]; then
                        echo "    Collection:"
                        echo "      Sources: $(jq '.sources | length' "$PIPELINE_DIR/$mode/$stage/output/collection.json")"
                        echo "      Compounds: $(jq '.compounds | length' "$PIPELINE_DIR/$mode/$stage/output/collection.json")"
                    fi
                    ;;
                process)
                    if [ -f "$PIPELINE_DIR/$mode/$stage/output/processing.json" ]; then
                        echo "    Processing:"
                        echo "      Processed: $(jq '.processed | length' "$PIPELINE_DIR/$mode/$stage/output/processing.json")"
                        echo "      Failed: $(jq '.failed | length' "$PIPELINE_DIR/$mode/$stage/output/processing.json")"
                    fi
                    ;;
                analyze)
                    if [ -f "$PIPELINE_DIR/$mode/$stage/output/analysis.json" ]; then
                        echo "    Analysis:"
                        echo "      Features: $(jq '.features | length' "$PIPELINE_DIR/$mode/$stage/output/analysis.json")"
                        echo "      Patterns: $(jq '.patterns | length' "$PIPELINE_DIR/$mode/$stage/output/analysis.json")"
                    fi
                    ;;
                predict)
                    if [ -f "$PIPELINE_DIR/$mode/$stage/output/predictions.json" ]; then
                        echo "    Predictions:"
                        echo "      Models: $(jq '.models | length' "$PIPELINE_DIR/$mode/$stage/output/predictions.json")"
                        echo "      Predictions: $(jq '.predictions | length' "$PIPELINE_DIR/$mode/$stage/output/predictions.json")"
                        if [ "$mode" = "bindingdb" ]; then
                            echo "      With activity: $(jq '.predictions[] | select(.activities != null) | .id' "$PIPELINE_DIR/$mode/$stage/output/predictions.json" | wc -l)"
                        fi
                    fi
                    ;;
                validate)
                    if [ -f "$PIPELINE_DIR/$mode/$stage/output/validation.json" ]; then
                        echo "    Validation:"
                        echo "      Tests: $(jq '.tests | length' "$PIPELINE_DIR/$mode/$stage/output/validation.json")"
                        echo "      Score: $(jq '.score' "$PIPELINE_DIR/$mode/$stage/output/validation.json")"
                    fi
                    ;;
                export)
                    if [ -f "$PIPELINE_DIR/$mode/$stage/output/export.json" ]; then
                        echo "    Export:"
                        echo "      Formats: $(jq '.formats | length' "$PIPELINE_DIR/$mode/$stage/output/export.json")"
                        echo "      Files: $(jq '.files | length' "$PIPELINE_DIR/$mode/$stage/output/export.json")"
                    fi
                    ;;
            esac
            
            # Show reports
            if [ -f "$PIPELINE_DIR/$mode/$stage/reports/stage.html" ]; then
                echo "    Stage report: $PIPELINE_DIR/$mode/$stage/reports/stage.html"
            fi
            if [ -f "$PIPELINE_DIR/$mode/$stage/reports/validation.html" ]; then
                echo "    Validation report: $PIPELINE_DIR/$mode/$stage/reports/validation.html"
            fi
            echo
        done
        
        # Show pipeline report
        if [ -f "$PIPELINE_DIR/pipeline_$mode.html" ]; then
            echo "  Pipeline report: $PIPELINE_DIR/pipeline_$mode.html"
        fi
        echo
    done
    
    echo "Cache:"
    echo "  Size: $(du -sh "$CACHE_DIR" | cut -f1)"
    echo "  Files: $(find "$CACHE_DIR" -type f | wc -l) files"
    echo
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing pipeline..."

# Create directory structure
create_dirs

# Run pipeline if requested
if [ "$RUN" = true ]; then
    run_pipeline
fi

# Monitor pipeline if requested
if [ "$MONITOR" = true ]; then
    monitor_pipeline
fi

# Validate pipeline if requested
if [ "$VALIDATE" = true ]; then
    validate_pipeline
fi

# Clean pipeline if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean pipeline? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_pipeline
    fi
fi

# Backup pipeline if requested
if [ "$BACKUP" = true ]; then
    backup_pipeline
fi

# Restore pipeline if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore pipeline? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_pipeline
    fi
fi

# Show statistics
show_stats

echo
echo "Pipeline management completed successfully!"
