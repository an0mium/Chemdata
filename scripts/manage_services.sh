#!/bin/bash
# Script to manage external services, microservices, and API integrations

# Exit on error
set -e

# Default values
DATA_DIR="data"
CONFIG_DIR="config"
SERVICE_DIR="services"
CACHE_DIR="cache/services"
LOG_DIR="logs"
SERVICE_TYPES="bindingdb,chembl,pubchem,swiss,community,social"
SERVICE_MODES="development,production,testing"
PORT_BASE=8000
HOST="localhost"
WORKERS=4
CACHE_TTL=3600
TIMEOUT=30
RETRIES=3
START=false
STOP=false
MONITOR=false
TEST=false
CHECK=false
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
        --config-dir)
            CONFIG_DIR="$2"
            shift 2
            ;;
        --service-dir)
            SERVICE_DIR="$2"
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
        --service-types)
            SERVICE_TYPES="$2"
            shift 2
            ;;
        --service-modes)
            SERVICE_MODES="$2"
            shift 2
            ;;
        --port-base)
            PORT_BASE="$2"
            shift 2
            ;;
        --host)
            HOST="$2"
            shift 2
            ;;
        --workers)
            WORKERS="$2"
            shift 2
            ;;
        --cache-ttl)
            CACHE_TTL="$2"
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
        --start)
            START=true
            shift
            ;;
        --stop)
            STOP=true
            shift
            ;;
        --monitor)
            MONITOR=true
            shift
            ;;
        --test)
            TEST=true
            shift
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
    
    # Service directories
    for service in ${SERVICE_TYPES//,/ }; do
        for mode in ${SERVICE_MODES//,/ }; do
            mkdir -p "$SERVICE_DIR/$service/$mode"/{config,data,cache,logs,reports}
        done
    done
    
    # Cache directory
    mkdir -p "$CACHE_DIR"
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to start services
start_services() {
    echo "Starting services..."
    
    # Track ports
    local port=$PORT_BASE
    
    for mode in ${SERVICE_MODES//,/ }; do
        echo "Starting services in $mode mode..."
        
        for service in ${SERVICE_TYPES//,/ }; do
            echo "Starting $service service..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli start-service"
            CMD="$CMD --data-dir $DATA_DIR"
            CMD="$CMD --config-dir $CONFIG_DIR"
            CMD="$CMD --service-dir $SERVICE_DIR/$service/$mode"
            CMD="$CMD --cache-dir $CACHE_DIR"
            CMD="$CMD --service $service"
            CMD="$CMD --mode $mode"
            CMD="$CMD --host $HOST"
            CMD="$CMD --port $port"
            CMD="$CMD --workers $WORKERS"
            CMD="$CMD --cache-ttl $CACHE_TTL"
            CMD="$CMD --timeout $TIMEOUT"
            CMD="$CMD --retries $RETRIES"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if service fails
            
            # Generate service report
            echo "Generating service report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-service-report \
                --service-dir "$SERVICE_DIR/$service/$mode" \
                --output "$SERVICE_DIR/$service/$mode/reports/service.html"
            
            # Increment port
            ((port++))
        done
    done
}

# Function to stop services
stop_services() {
    echo "Stopping services..."
    
    for mode in ${SERVICE_MODES//,/ }; do
        echo "Stopping services in $mode mode..."
        
        for service in ${SERVICE_TYPES//,/ }; do
            echo "Stopping $service service..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli stop-service"
            CMD="$CMD --service-dir $SERVICE_DIR/$service/$mode"
            CMD="$CMD --service $service"
            CMD="$CMD --mode $mode"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if stop fails
        done
    done
}

# Function to monitor services
monitor_services() {
    echo "Monitoring services..."
    
    # Build command
    CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli monitor-services"
    CMD="$CMD --service-dir $SERVICE_DIR"
    CMD="$CMD --log-dir $LOG_DIR"
    
    # Run command
    echo "Running: $CMD"
    $CMD
}

# Function to test services
test_services() {
    echo "Testing services..."
    
    for mode in ${SERVICE_MODES//,/ }; do
        echo "Testing services in $mode mode..."
        
        for service in ${SERVICE_TYPES//,/ }; do
            echo "Testing $service service..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli test-service"
            CMD="$CMD --service-dir $SERVICE_DIR/$service/$mode"
            CMD="$CMD --service $service"
            CMD="$CMD --mode $mode"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if tests fail
            
            # Generate test report
            echo "Generating test report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-test-report \
                --service-dir "$SERVICE_DIR/$service/$mode" \
                --output "$SERVICE_DIR/$service/$mode/reports/test.html"
        done
    done
}

# Function to check services
check_services() {
    echo "Checking services..."
    
    for mode in ${SERVICE_MODES//,/ }; do
        echo "Checking services in $mode mode..."
        
        for service in ${SERVICE_TYPES//,/ }; do
            echo "Checking $service service..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli check-service"
            CMD="$CMD --service-dir $SERVICE_DIR/$service/$mode"
            CMD="$CMD --service $service"
            CMD="$CMD --mode $mode"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if check fails
            
            # Generate check report
            echo "Generating check report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-check-report \
                --service-dir "$SERVICE_DIR/$service/$mode" \
                --output "$SERVICE_DIR/$service/$mode/reports/check.html"
        done
    done
}

# Function to clean services
clean_services() {
    echo "Cleaning services..."
    
    # Clean service directories
    rm -rf "$SERVICE_DIR"/*
    
    # Clean cache
    rm -rf "$CACHE_DIR"/*
}

# Function to backup services
backup_services() {
    echo "Backing up services..."
    
    # Create backup directory
    BACKUP_DIR="backups/services_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$SERVICE_DIR" "$BACKUP_DIR/"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore services
restore_services() {
    echo "Restoring services..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/services_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$SERVICE_DIR"
    fi
    
    cp -r "$BACKUP_DIR/services" ./
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Services restored from: $BACKUP_DIR"
}

# Function to show service statistics
show_stats() {
    echo "Service statistics:"
    echo
    
    echo "Service files:"
    for mode in ${SERVICE_MODES//,/ }; do
        echo "$mode mode:"
        for service in ${SERVICE_TYPES//,/ }; do
            echo "  $service service:"
            
            # Config files
            echo "    Config:"
            echo "      Files: $(find "$SERVICE_DIR/$service/$mode/config" -type f | wc -l) files"
            echo "      Size: $(du -sh "$SERVICE_DIR/$service/$mode/config" | cut -f1)"
            
            # Data files
            echo "    Data:"
            echo "      Files: $(find "$SERVICE_DIR/$service/$mode/data" -type f | wc -l) files"
            echo "      Size: $(du -sh "$SERVICE_DIR/$service/$mode/data" | cut -f1)"
            
            # Cache files
            echo "    Cache:"
            echo "      Files: $(find "$SERVICE_DIR/$service/$mode/cache" -type f | wc -l) files"
            echo "      Size: $(du -sh "$SERVICE_DIR/$service/$mode/cache" | cut -f1)"
            
            # Count entries by service
            case "$service" in
                bindingdb)
                    if [ -f "$SERVICE_DIR/$service/$mode/data/compounds.json" ]; then
                        echo "    Compounds:"
                        echo "      Total: $(jq '.compounds | length' "$SERVICE_DIR/$service/$mode/data/compounds.json")"
                        echo "      With activity: $(jq '.compounds[] | select(.activities != null) | .id' "$SERVICE_DIR/$service/$mode/data/compounds.json" | wc -l)"
                    fi
                    if [ -f "$SERVICE_DIR/$service/$mode/logs/metrics.json" ]; then
                        echo "    Metrics:"
                        echo "      Requests: $(jq '.requests' "$SERVICE_DIR/$service/$mode/logs/metrics.json")"
                        echo "      Compounds: $(jq '.compounds' "$SERVICE_DIR/$service/$mode/logs/metrics.json")"
                    fi
                    ;;
                chembl)
                    if [ -f "$SERVICE_DIR/$service/$mode/data/molecules.json" ]; then
                        echo "    Molecules:"
                        echo "      Total: $(jq '.molecules | length' "$SERVICE_DIR/$service/$mode/data/molecules.json")"
                        echo "      With assays: $(jq '.molecules[] | select(.assays != null) | .id' "$SERVICE_DIR/$service/$mode/data/molecules.json" | wc -l)"
                    fi
                    if [ -f "$SERVICE_DIR/$service/$mode/logs/metrics.json" ]; then
                        echo "    Metrics:"
                        echo "      Queries: $(jq '.queries' "$SERVICE_DIR/$service/$mode/logs/metrics.json")"
                        echo "      Results: $(jq '.results' "$SERVICE_DIR/$service/$mode/logs/metrics.json")"
                    fi
                    ;;
                pubchem)
                    if [ -f "$SERVICE_DIR/$service/$mode/data/substances.json" ]; then
                        echo "    Substances:"
                        echo "      Total: $(jq '.substances | length' "$SERVICE_DIR/$service/$mode/data/substances.json")"
                        echo "      With bioactivity: $(jq '.substances[] | select(.bioactivity != null) | .id' "$SERVICE_DIR/$service/$mode/data/substances.json" | wc -l)"
                    fi
                    if [ -f "$SERVICE_DIR/$service/$mode/logs/metrics.json" ]; then
                        echo "    Metrics:"
                        echo "      Searches: $(jq '.searches' "$SERVICE_DIR/$service/$mode/logs/metrics.json")"
                        echo "      Compounds: $(jq '.compounds' "$SERVICE_DIR/$service/$mode/logs/metrics.json")"
                    fi
                    ;;
                swiss)
                    if [ -f "$SERVICE_DIR/$service/$mode/data/predictions.json" ]; then
                        echo "    Predictions:"
                        echo "      Total: $(jq '.predictions | length' "$SERVICE_DIR/$service/$mode/data/predictions.json")"
                        echo "      High confidence: $(jq '.predictions[] | select(.confidence > 0.7) | .id' "$SERVICE_DIR/$service/$mode/data/predictions.json" | wc -l)"
                    fi
                    if [ -f "$SERVICE_DIR/$service/$mode/logs/metrics.json" ]; then
                        echo "    Metrics:"
                        echo "      Predictions: $(jq '.predictions' "$SERVICE_DIR/$service/$mode/logs/metrics.json")"
                        echo "      Properties: $(jq '.properties' "$SERVICE_DIR/$service/$mode/logs/metrics.json")"
                    fi
                    ;;
                community)
                    if [ -f "$SERVICE_DIR/$service/$mode/data/reports.json" ]; then
                        echo "    Reports:"
                        echo "      Total: $(jq '.reports | length' "$SERVICE_DIR/$service/$mode/data/reports.json")"
                        echo "      With effects: $(jq '.reports[] | select(.effects != null) | .id' "$SERVICE_DIR/$service/$mode/data/reports.json" | wc -l)"
                    fi
                    if [ -f "$SERVICE_DIR/$service/$mode/logs/metrics.json" ]; then
                        echo "    Metrics:"
                        echo "      Reports: $(jq '.reports' "$SERVICE_DIR/$service/$mode/logs/metrics.json")"
                        echo "      Sources: $(jq '.sources' "$SERVICE_DIR/$service/$mode/logs/metrics.json")"
                    fi
                    ;;
                social)
                    if [ -f "$SERVICE_DIR/$service/$mode/data/posts.json" ]; then
                        echo "    Posts:"
                        echo "      Total: $(jq '.posts | length' "$SERVICE_DIR/$service/$mode/data/posts.json")"
                        echo "      With compounds: $(jq '.posts[] | select(.compounds != null) | .id' "$SERVICE_DIR/$service/$mode/data/posts.json" | wc -l)"
                    fi
                    if [ -f "$SERVICE_DIR/$service/$mode/logs/metrics.json" ]; then
                        echo "    Metrics:"
                        echo "      Posts: $(jq '.posts' "$SERVICE_DIR/$service/$mode/logs/metrics.json")"
                        echo "      Mentions: $(jq '.mentions' "$SERVICE_DIR/$service/$mode/logs/metrics.json")"
                    fi
                    ;;
            esac
            
            # Show reports
            if [ -f "$SERVICE_DIR/$service/$mode/reports/service.html" ]; then
                echo "    Service report: $SERVICE_DIR/$service/$mode/reports/service.html"
            fi
            if [ -f "$SERVICE_DIR/$service/$mode/reports/test.html" ]; then
                echo "    Test report: $SERVICE_DIR/$service/$mode/reports/test.html"
            fi
            if [ -f "$SERVICE_DIR/$service/$mode/reports/check.html" ]; then
                echo "    Check report: $SERVICE_DIR/$service/$mode/reports/check.html"
            fi
            echo
        done
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
echo "Managing services..."

# Create directory structure
create_dirs

# Start services if requested
if [ "$START" = true ]; then
    start_services
fi

# Stop services if requested
if [ "$STOP" = true ]; then
    stop_services
fi

# Monitor services if requested
if [ "$MONITOR" = true ]; then
    monitor_services
fi

# Test services if requested
if [ "$TEST" = true ]; then
    test_services
fi

# Check services if requested
if [ "$CHECK" = true ]; then
    check_services
fi

# Clean services if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean services? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_services
    fi
fi

# Backup services if requested
if [ "$BACKUP" = true ]; then
    backup_services
fi

# Restore services if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore services? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_services
    fi
fi

# Show statistics
show_stats

echo
echo "Service management completed successfully!"
