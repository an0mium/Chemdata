#!/bin/bash
# Script to manage data visualizations and plots

# Exit on error
set -e

# Default values
DATA_DIR="data"
PLOT_DIR="plots"
LOG_DIR="logs"
PLOT_TYPES="structure,activity,predictions,community,safety"
PLOT_FORMATS="png,svg,html,pdf"
PLOT_THEME="light"
RESOLUTION="1200x800"
DPI=300
INTERACTIVE=false
PLOT=false
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
        --plot-dir)
            PLOT_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --plot-types)
            PLOT_TYPES="$2"
            shift 2
            ;;
        --plot-formats)
            PLOT_FORMATS="$2"
            shift 2
            ;;
        --plot-theme)
            PLOT_THEME="$2"
            shift 2
            ;;
        --resolution)
            RESOLUTION="$2"
            shift 2
            ;;
        --dpi)
            DPI="$2"
            shift 2
            ;;
        --interactive)
            INTERACTIVE=true
            shift
            ;;
        --plot)
            PLOT=true
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
    
    # Plot directories
    for type in ${PLOT_TYPES//,/ }; do
        for format in ${PLOT_FORMATS//,/ }; do
            mkdir -p "$PLOT_DIR/$type/$format"
        done
    done
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to generate plots
generate_plots() {
    echo "Generating plots..."
    
    for type in ${PLOT_TYPES//,/ }; do
        for format in ${PLOT_FORMATS//,/ }; do
            echo "Generating $type plots in $format format..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli generate-plots"
            CMD="$CMD --data-dir $DATA_DIR"
            CMD="$CMD --plot-dir $PLOT_DIR/$type/$format"
            CMD="$CMD --plot-type $type"
            CMD="$CMD --format $format"
            CMD="$CMD --theme $PLOT_THEME"
            CMD="$CMD --resolution $RESOLUTION"
            CMD="$CMD --dpi $DPI"
            CMD="$CMD --log-dir $LOG_DIR"
            
            if [ "$INTERACTIVE" = true ] && [ "$format" = "html" ]; then
                CMD="$CMD --interactive"
            fi
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if plotting fails
        done
    done
}

# Function to clean plots
clean_plots() {
    echo "Cleaning plots..."
    
    # Clean plot directories
    rm -rf "$PLOT_DIR"/*
}

# Function to backup plots
backup_plots() {
    echo "Backing up plots..."
    
    # Create backup directory
    BACKUP_DIR="backups/plots_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$PLOT_DIR" "$BACKUP_DIR/"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore plots
restore_plots() {
    echo "Restoring plots..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/plots_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$PLOT_DIR"
    fi
    
    cp -r "$BACKUP_DIR/plots" ./
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Plots restored from: $BACKUP_DIR"
}

# Function to show plot statistics
show_stats() {
    echo "Plot statistics:"
    echo
    
    echo "Plot files:"
    for type in ${PLOT_TYPES//,/ }; do
        echo "$type plots:"
        for format in ${PLOT_FORMATS//,/ }; do
            echo "  $format format:"
            echo "    Files: $(find "$PLOT_DIR/$type/$format" -type f | wc -l) files"
            echo "    Size: $(du -sh "$PLOT_DIR/$type/$format" | cut -f1)"
            
            # Count by plot subtype
            case "$type" in
                structure)
                    echo "    2D structures: $(find "$PLOT_DIR/$type/$format" -name "*_2d.$format" | wc -l)"
                    echo "    3D structures: $(find "$PLOT_DIR/$type/$format" -name "*_3d.$format" | wc -l)"
                    ;;
                activity)
                    echo "    Binding plots: $(find "$PLOT_DIR/$type/$format" -name "*_binding.$format" | wc -l)"
                    echo "    SAR plots: $(find "$PLOT_DIR/$type/$format" -name "*_sar.$format" | wc -l)"
                    ;;
                predictions)
                    echo "    Confidence plots: $(find "$PLOT_DIR/$type/$format" -name "*_confidence.$format" | wc -l)"
                    echo "    Distribution plots: $(find "$PLOT_DIR/$type/$format" -name "*_distribution.$format" | wc -l)"
                    ;;
                community)
                    echo "    Activity reports: $(find "$PLOT_DIR/$type/$format" -name "*_activity.$format" | wc -l)"
                    echo "    Experience reports: $(find "$PLOT_DIR/$type/$format" -name "*_experience.$format" | wc -l)"
                    ;;
                safety)
                    echo "    Risk plots: $(find "$PLOT_DIR/$type/$format" -name "*_risk.$format" | wc -l)"
                    echo "    Warning plots: $(find "$PLOT_DIR/$type/$format" -name "*_warning.$format" | wc -l)"
                    ;;
            esac
        done
        echo
    done
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing visualizations..."

# Create directory structure
create_dirs

# Generate plots if requested
if [ "$PLOT" = true ]; then
    generate_plots
fi

# Clean plots if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean plots? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_plots
    fi
fi

# Backup plots if requested
if [ "$BACKUP" = true ]; then
    backup_plots
fi

# Restore plots if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore plots? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_plots
    fi
fi

# Show statistics
show_stats

echo
echo "Visualization management completed successfully!"
