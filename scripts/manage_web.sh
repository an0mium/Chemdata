#!/bin/bash
# Script to manage web interface, visualization components, and server

# Exit on error
set -e

# Default values
DATA_DIR="data/web"
WEB_DIR="web"
STATIC_DIR="web/static"
TEMPLATE_DIR="web/templates"
CACHE_DIR="cache/web"
LOG_DIR="logs"
COMPONENT_TYPES="frontend,backend,api,docs,analysis,export,input,visualization"
DASHBOARD_TYPES="main,detail,list,plot"
WEB_COMPONENTS="dashboard,visualization,analysis,export"
WEB_MODES="development,production,testing"
SERVER_MODES="development,production,testing"
DATASET_TYPES="bindingdb,chembl,pubchem,swiss,community,literature"
PORT=3000
HOST="localhost"
WORKERS=4
TIMEOUT=30
START=false
STOP=false
RELOAD=false
BUILD=false
TEST=false
DEPLOY=false
CLEAN=false
BACKUP=false
RESTORE=false
DEBUG=false
FORCE=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --web-dir)
            WEB_DIR="$2"
            shift 2
            ;;
        --static-dir)
            STATIC_DIR="$2"
            shift 2
            ;;
        --template-dir)
            TEMPLATE_DIR="$2"
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
        --component-types)
            COMPONENT_TYPES="$2"
            shift 2
            ;;
        --dashboard-types)
            DASHBOARD_TYPES="$2"
            shift 2
            ;;
        --web-components)
            WEB_COMPONENTS="$2"
            shift 2
            ;;
        --web-modes)
            WEB_MODES="$2"
            shift 2
            ;;
        --server-modes)
            SERVER_MODES="$2"
            shift 2
            ;;
        --dataset-types)
            DATASET_TYPES="$2"
            shift 2
            ;;
        --port)
            PORT="$2"
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
        --timeout)
            TIMEOUT="$2"
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
        --reload)
            RELOAD=true
            shift
            ;;
        --build)
            BUILD=true
            shift
            ;;
        --test)
            TEST=true
            shift
            ;;
        --deploy)
            DEPLOY=true
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
        --debug)
            DEBUG=true
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

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to create directory structure
create_dirs() {
    echo "Creating directory structure..."
    
    # Web directories
    mkdir -p "$WEB_DIR"
    
    # Static directories
    mkdir -p "$STATIC_DIR"/{css,js,img,fonts}
    
    # Template directories
    mkdir -p "$TEMPLATE_DIR"/{components,layouts,pages}
    
    # Component directories
    for component in ${COMPONENT_TYPES//,/ }; do
        for mode in ${SERVER_MODES//,/ }; do
            mkdir -p "$WEB_DIR/$component/$mode"/{src,dist,static,reports}
            mkdir -p "$WEB_DIR/components/$component"/{templates,static/{css,js}}
        done
    done
    
    # Dashboard directories
    for type in ${DASHBOARD_TYPES//,/ }; do
        mkdir -p "$WEB_DIR/dashboard/$type"/{templates,static/{css,js}}
    done
    
    # Web component directories
    for mode in ${WEB_MODES//,/ }; do
        for component in ${WEB_COMPONENTS//,/ }; do
            mkdir -p "$STATIC_DIR/$mode/$component"/{css,js,img}
            mkdir -p "$TEMPLATE_DIR/$mode/$component"
        done
    done
    
    # Data directories
    for type in ${DATASET_TYPES//,/ }; do
        mkdir -p "$DATA_DIR/$type"/{raw,processed,cache}
    done
    
    # Cache directory
    mkdir -p "$CACHE_DIR"
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to build web assets
build_assets() {
    echo "Building web assets..."
    
    # Build components
    for component in ${COMPONENT_TYPES//,/ }; do
        for mode in ${SERVER_MODES//,/ }; do
            echo "Building $component in $mode mode..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli build-web"
            CMD="$CMD --data-dir $DATA_DIR"
            CMD="$CMD --web-dir $WEB_DIR/$component/$mode"
            CMD="$CMD --component $component"
            CMD="$CMD --mode $mode"
            CMD="$CMD --log-dir $LOG_DIR"
            
            if [ "$DEBUG" = true ]; then
                CMD="$CMD --debug"
            fi
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if build fails
            
            # Generate build report
            echo "Generating build report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-build-report \
                --web-dir "$WEB_DIR/$component/$mode" \
                --output "$WEB_DIR/$component/$mode/reports/build.html"
        done
    done
    
    # Build web components
    for mode in ${WEB_MODES//,/ }; do
        echo "Building assets for $mode mode..."
        
        for component in ${WEB_COMPONENTS//,/ }; do
            echo "Building $component component..."
            
            # Build command
            CMD="python -m binding_data_processor.web.cli build-assets"
            CMD="$CMD --static-dir $STATIC_DIR/$mode/$component"
            CMD="$CMD --template-dir $TEMPLATE_DIR/$mode/$component"
            CMD="$CMD --mode $mode"
            CMD="$CMD --component $component"
            CMD="$CMD --log-dir $LOG_DIR"
            
            if [ "$DEBUG" = true ]; then
                CMD="$CMD --debug"
            fi
            
            # Run build
            echo "Running: $CMD"
            $CMD || true  # Continue even if build fails
        done
    done
    
    # Build CSS
    echo "Building CSS..."
    for component in ${COMPONENT_TYPES//,/ }; do
        if [ -d "$WEB_DIR/components/$component/static/css" ]; then
            cat "$WEB_DIR/components/$component/static/css"/*.css > "$STATIC_DIR/css/$component.css"
        fi
    done
    for type in ${DASHBOARD_TYPES//,/ }; do
        if [ -d "$WEB_DIR/dashboard/$type/static/css" ]; then
            cat "$WEB_DIR/dashboard/$type/static/css"/*.css > "$STATIC_DIR/css/$type.css"
        fi
    done
    
    # Build JavaScript
    echo "Building JavaScript..."
    for component in ${COMPONENT_TYPES//,/ }; do
        if [ -d "$WEB_DIR/components/$component/static/js" ]; then
            cat "$WEB_DIR/components/$component/static/js"/*.js > "$STATIC_DIR/js/$component.js"
        fi
    done
    for type in ${DASHBOARD_TYPES//,/ }; do
        if [ -d "$WEB_DIR/dashboard/$type/static/js" ]; then
            cat "$WEB_DIR/dashboard/$type/static/js"/*.js > "$STATIC_DIR/js/$type.js"
        fi
    done
}

# Function to start web server
start_server() {
    echo "Starting web server..."
    
    # Build server command
    CMD="python -m binding_data_processor.web.cli start-server"
    CMD="$CMD --data-dir $DATA_DIR"
    CMD="$CMD --static-dir $STATIC_DIR"
    CMD="$CMD --template-dir $TEMPLATE_DIR"
    CMD="$CMD --cache-dir $CACHE_DIR"
    CMD="$CMD --host $HOST"
    CMD="$CMD --port $PORT"
    CMD="$CMD --workers $WORKERS"
    CMD="$CMD --timeout $TIMEOUT"
    CMD="$CMD --log-dir $LOG_DIR"
    
    if [ "$DEBUG" = true ]; then
        CMD="$CMD --debug"
    fi
    
    # Run server
    echo "Running: $CMD"
    $CMD
}

# Function to stop web server
stop_server() {
    echo "Stopping web server..."
    
    # Build stop command
    CMD="python -m binding_data_processor.web.cli stop-server"
    CMD="$CMD --log-dir $LOG_DIR"
    
    # Run command
    echo "Running: $CMD"
    $CMD
}

# Function to reload web server
reload_server() {
    echo "Reloading web server..."
    
    # Build reload command
    CMD="python -m binding_data_processor.web.cli reload-server"
    CMD="$CMD --log-dir $LOG_DIR"
    
    # Run command
    echo "Running: $CMD"
    $CMD
}

# Function to test web interface
test_web() {
    echo "Testing web interface..."
    
    # Test components
    for component in ${COMPONENT_TYPES//,/ }; do
        for mode in ${SERVER_MODES//,/ }; do
            echo "Testing $component in $mode mode..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli test-web"
            CMD="$CMD --web-dir $WEB_DIR/$component/$mode"
            CMD="$CMD --component $component"
            CMD="$CMD --mode $mode"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run command
            echo "Running: $CMD"
            $CMD || true  # Continue even if tests fail
            
            # Generate test report
            echo "Generating test report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-test-report \
                --web-dir "$WEB_DIR/$component/$mode" \
                --output "$WEB_DIR/$component/$mode/reports/test.html"
        done
    done
    
    # Test dashboards
    for type in ${DASHBOARD_TYPES//,/ }; do
        echo "Testing $type dashboard..."
        
        # Build command
        CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli test-web"
        CMD="$CMD --web-dir $WEB_DIR"
        CMD="$CMD --dashboard-type $type"
        CMD="$CMD --log-dir $LOG_DIR"
        
        # Run command
        echo "Running: $CMD"
        $CMD || true  # Continue even if tests fail
    done
    
    # Generate combined test report
    echo "Generating combined test report..."
    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-web-test-report \
        --web-dir "$WEB_DIR" \
        --output "$WEB_DIR/test_report.html"
}

# Function to deploy web application
deploy_app() {
    echo "Deploying web application..."
    
    # Build deploy command
    CMD="python -m binding_data_processor.web.cli deploy-app"
    CMD="$CMD --data-dir $DATA_DIR"
    CMD="$CMD --static-dir $STATIC_DIR"
    CMD="$CMD --template-dir $TEMPLATE_DIR"
    CMD="$CMD --cache-dir $CACHE_DIR"
    CMD="$CMD --log-dir $LOG_DIR"
    
    if [ "$DEBUG" = true ]; then
        CMD="$CMD --debug"
    fi
    
    # Run deploy
    echo "Running: $CMD"
    $CMD
}

# Function to clean web files
clean_web() {
    echo "Cleaning web files..."
    
    # Clean web directories
    rm -rf "$WEB_DIR"/*
    
    # Clean static directories
    rm -rf "$STATIC_DIR"/*
    
    # Clean template directories
    rm -rf "$TEMPLATE_DIR"/*
    
    # Clean data directories
    rm -rf "$DATA_DIR"/*
    
    # Clean cache
    rm -rf "$CACHE_DIR"/*
}

# Function to backup web files
backup_web() {
    echo "Backing up web files..."
    
    # Create backup directory
    BACKUP_DIR="backups/web_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$WEB_DIR" "$BACKUP_DIR/"
    cp -r "$STATIC_DIR" "$BACKUP_DIR/"
    cp -r "$TEMPLATE_DIR" "$BACKUP_DIR/"
    cp -r "$DATA_DIR" "$BACKUP_DIR/"
    
    # Create encrypted backup if gpg is available
    if command_exists gpg; then
        echo "Encrypting backup..."
        tar -czf - "$BACKUP_DIR" | gpg --symmetric --output "$BACKUP_DIR.tar.gz.gpg"
        rm -rf "$BACKUP_DIR"
        echo "Encrypted backup saved to: $BACKUP_DIR.tar.gz.gpg"
    else
        # Create unencrypted backup
        tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
        rm -rf "$BACKUP_DIR"
        echo "Backup saved to: $BACKUP_DIR.tar.gz"
    fi
}

# Function to restore web files
restore_web() {
    echo "Restoring web files..."
    
    # Find latest backup
    local LATEST_BACKUP=""
    if command_exists gpg; then
        LATEST_BACKUP=$(ls -t backups/web_*.tar.gz.gpg 2>/dev/null | head -n1)
        if [ -n "$LATEST_BACKUP" ]; then
            echo "Decrypting backup: $LATEST_BACKUP"
            gpg --decrypt "$LATEST_BACKUP" | tar -xz
            BACKUP_DIR="${LATEST_BACKUP%.tar.gz.gpg}"
        fi
    fi
    
    if [ -z "$LATEST_BACKUP" ]; then
        LATEST_BACKUP=$(ls -t backups/web_*.tar.gz 2>/dev/null | head -n1)
        if [ -n "$LATEST_BACKUP" ]; then
            echo "Extracting backup: $LATEST_BACKUP"
            tar -xzf "$LATEST_BACKUP"
            BACKUP_DIR="${LATEST_BACKUP%.tar.gz}"
        fi
    fi
    
    if [ -z "$BACKUP_DIR" ]; then
        echo "No backup found"
        exit 1
    fi
    
    # Restore directories
    if [ "$FORCE" = true ]; then
        rm -rf "$WEB_DIR" "$STATIC_DIR" "$TEMPLATE_DIR" "$DATA_DIR"
    fi
    
    cp -r "$BACKUP_DIR/web" ./
    cp -r "$BACKUP_DIR/static" ./web/
    cp -r "$BACKUP_DIR/templates" ./web/
    cp -r "$BACKUP_DIR/data/web" ./data/
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Web files restored from: $BACKUP_DIR"
}

# Function to show web statistics
show_stats() {
    echo "Web statistics:"
    echo
    
    echo "Component files:"
    for component in ${COMPONENT_TYPES//,/ }; do
        echo "$component component:"
        
        # Show mode-specific files
        for mode in ${SERVER_MODES//,/ }; do
            echo "  $mode mode:"
            
            # Source files
            echo "    Source:"
            echo "      Files: $(find "$WEB_DIR/$component/$mode/src" -type f | wc -l) files"
            echo "      Size: $(du -sh "$WEB_DIR/$component/$mode/src" | cut -f1)"
            
            # Distribution files
            echo "    Distribution:"
            echo "      Files: $(find "$WEB_DIR/$component/$mode/dist" -type f | wc -l) files"
            echo "      Size: $(du -sh "$WEB_DIR/$component/$mode/dist" | cut -f1)"
            
            # Show reports
            if [ -f "$WEB_DIR/$component/$mode/reports/build.html" ]; then
                echo "    Build report: $WEB_DIR/$component/$mode/reports/build.html"
            fi
            if [ -f "$WEB_DIR/$component/$mode/reports/test.html" ]; then
                echo "    Test report: $WEB_DIR/$component/$mode/reports/test.html"
            fi
        done
        
        # Show component-specific files
        echo "  Component files:"
        
        # Template files
        echo "    Templates:"
        echo "      Files: $(find "$WEB_DIR/components/$component/templates" -type f | wc -l) files"
        echo "      Size: $(du -sh "$WEB_DIR/components/$component/templates" | cut -f1)"
        
        # Static files
        if [ -d "$WEB_DIR/components/$component/static" ]; then
            echo "    Static:"
            echo "      CSS: $(find "$WEB_DIR/components/$component/static/css" -type f | wc -l) files"
            echo "      JavaScript: $(find "$WEB_DIR/components/$component/static/js" -type f | wc -l) files"
            echo "      Size: $(du -sh "$WEB_DIR/components/$component/static" | cut -f1)"
        fi
        echo
    done
    
    echo "Dashboard files:"
    for type in ${DASHBOARD_TYPES//,/ }; do
        echo "$type dashboard:"
        
        # Template files
        echo "  Templates:"
        echo "    Files: $(find "$WEB_DIR/dashboard/$type/templates" -type f | wc -l) files"
        echo "    Size: $(du -sh "$WEB_DIR/dashboard/$type/templates" | cut -f1)"
        
        # Static files
        if [ -d "$WEB_DIR/dashboard/$type/static" ]; then
            echo "  Static:"
            echo "    CSS: $(find "$WEB_DIR/dashboard/$type/static/css" -type f | wc -l) files"
            echo "    JavaScript: $(find "$WEB_DIR/dashboard/$type/static/js" -type f | wc -l) files"
            echo "    Size: $(du -sh "$WEB_DIR/dashboard/$type/static" | cut -f1)"
        fi
        echo
    done
    
    echo "Web component files:"
    for mode in ${WEB_MODES//,/ }; do
        echo "$mode mode:"
        for component in ${WEB_COMPONENTS//,/ }; do
            echo "  $component component:"
            
            # CSS files
            echo "    CSS files:"
            echo "      Files: $(find "$STATIC_DIR/$mode/$component/css" -type f | wc -l) files"
            echo "      Size: $(du -sh "$STATIC_DIR/$mode/$component/css" | cut -f1)"
            
            # JavaScript files
            echo "    JavaScript files:"
            echo "      Files: $(find "$STATIC_DIR/$mode/$component/js" -type f | wc -l) files"
            echo "      Size: $(du -sh "$STATIC_DIR/$mode/$component/js" | cut -f1)"
            
            # Image files
            echo "    Image files:"
            echo "      Files: $(find "$STATIC_DIR/$mode/$component/img" -type f | wc -l) files"
            echo "      Size: $(du -sh "$STATIC_DIR/$mode/$component/img" | cut -f1)"
            echo
        done
    done
    
    echo "Static files:"
    echo "  CSS: $(find "$STATIC_DIR/css" -type f | wc -l) files"
    echo "  JavaScript: $(find "$STATIC_DIR/js" -type f | wc -l) files"
    echo "  Images: $(find "$STATIC_DIR/img" -type f | wc -l) files"
    echo "  Fonts: $(find "$STATIC_DIR/fonts" -type f | wc -l) files"
    echo "  Size: $(du -sh "$STATIC_DIR" | cut -f1)"
    echo
    
    echo "Template files:"
    echo "  Components: $(find "$TEMPLATE_DIR/components" -type f | wc -l) files"
    echo "  Layouts: $(find "$TEMPLATE_DIR/layouts" -type f | wc -l) files"
    echo "  Pages: $(find "$TEMPLATE_DIR/pages" -type f | wc -l) files"
    echo "  Size: $(du -sh "$TEMPLATE_DIR" | cut -f1)"
    echo
    
    echo "Data files:"
    for type in ${DATASET_TYPES//,/ }; do
        echo "$type dataset:"
        
        # Raw data
        echo "  Raw data:"
        echo "    Files: $(find "$DATA_DIR/$type/raw" -type f | wc -l) files"
        echo "    Size: $(du -sh "$DATA_DIR/$type/raw" | cut -f1)"
        
        # Processed data
        echo "  Processed data:"
        echo "    Files: $(find "$DATA_DIR/$type/processed" -type f | wc -l) files"
        echo "    Size: $(du -sh "$DATA_DIR/$type/processed" | cut -f1)"
        
        # Cache data
        echo "  Cache data:"
        echo "    Files: $(find "$DATA_DIR/$type/cache" -type f | wc -l) files"
        echo "    Size: $(du -sh "$DATA_DIR/$type/cache" | cut -f1)"
        echo
    done
    
    echo "Cache:"
    echo "  Size: $(du -sh "$CACHE_DIR" | cut -f1)"
    echo "  Files: $(find "$CACHE_DIR" -type f | wc -l) files"
    echo
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
    
    # Show test report
    if [ -f "$WEB_DIR/test_report.html" ]; then
        echo
        echo "Test report: $WEB_DIR/test_report.html"
    fi
}

# Main process
echo "Managing web interface..."

# Create directory structure
create_dirs

# Build assets if requested
if [ "$BUILD" = true ]; then
    build_assets
fi

# Start server if requested
if [ "$START" = true ]; then
    start_server
fi

# Stop server if requested
if [ "$STOP" = true ]; then
    stop_server
fi

# Reload server if requested
if [ "$RELOAD" = true ]; then
    reload_server
fi

# Test web interface if requested
if [ "$TEST" = true ]; then
    test_web
fi

# Deploy application if requested
if [ "$DEPLOY" = true ]; then
    deploy_app
fi

# Clean web files if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean web files? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_web
    fi
fi

# Backup web files if requested
if [ "$BACKUP" = true ]; then
    backup_web
fi

# Restore web files if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore web files? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_web
    fi
fi

# Show statistics
show_stats

echo
echo "Web management completed successfully!"
