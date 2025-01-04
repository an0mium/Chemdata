#!/bin/bash
# Script to manage API integrations, endpoints, documentation, and versioning

# Exit on error
set -e

# Default values
CONFIG_FILE=".env"
API_DIR="api"
DOCS_DIR="docs/api"
SPEC_DIR="specs/api"
CACHE_DIR="cache/api"
LOG_DIR="logs"
API_TYPES="rest,graphql,websocket,chembl,pubchem,swiss,reddit,twitter,discord,bluesky"
API_VERSIONS="v1,v2,v3"
API_MODES="development,production,testing"
PORT_BASE=9000
HOST="localhost"
WORKERS=4
CACHE_TTL=3600
TIMEOUT=30
RETRIES=3
START=false
STOP=false
MONITOR=false
TEST=false
DOCS=false
CLEAN=false
BACKUP=false
RESTORE=false
FORCE=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --api-dir)
            API_DIR="$2"
            shift 2
            ;;
        --docs-dir)
            DOCS_DIR="$2"
            shift 2
            ;;
        --spec-dir)
            SPEC_DIR="$2"
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
        --api-types)
            API_TYPES="$2"
            shift 2
            ;;
        --api-versions)
            API_VERSIONS="$2"
            shift 2
            ;;
        --api-modes)
            API_MODES="$2"
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
        --docs)
            DOCS=true
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

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to get API key from config
get_api_key() {
    local api="$1"
    local key_name="${api^^}_API_KEY"
    grep "^$key_name=" "$CONFIG_FILE" 2>/dev/null | cut -d'=' -f2 || echo ""
}

# Function to set API key in config
set_api_key() {
    local api="$1"
    local key="$2"
    local key_name="${api^^}_API_KEY"
    
    # Create config file if it doesn't exist
    touch "$CONFIG_FILE"
    
    # Update or add key
    if grep -q "^$key_name=" "$CONFIG_FILE"; then
        sed -i "s|^$key_name=.*|$key_name=$key|" "$CONFIG_FILE"
    else
        echo "$key_name=$key" >> "$CONFIG_FILE"
    fi
}

# Function to create directory structure
create_dirs() {
    echo "Creating directory structure..."
    
    # API directories
    for version in ${API_VERSIONS//,/ }; do
        for mode in ${API_MODES//,/ }; do
            for type in ${API_TYPES//,/ }; do
                mkdir -p "$API_DIR/$version/$mode/$type"/{endpoints,schemas,tests,docs}
            done
        done
    done
    
    # Documentation directories
    for version in ${API_VERSIONS//,/ }; do
        mkdir -p "$DOCS_DIR/$version"/{reference,guides,examples}
    done
    
    # Specification directories
    for version in ${API_VERSIONS//,/ }; do
        mkdir -p "$SPEC_DIR/$version"/{openapi,graphql,asyncapi}
    done
    
    # Cache directory
    mkdir -p "$CACHE_DIR"
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to start APIs
start_apis() {
    echo "Starting APIs..."
    
    # Track ports
    local port=$PORT_BASE
    
    for version in ${API_VERSIONS//,/ }; do
        for mode in ${API_MODES//,/ }; do
            for type in ${API_TYPES//,/ }; do
                echo "Starting $type API $version in $mode mode..."
                
                # Get API key if external service
                local key=""
                case "$type" in
                    chembl|pubchem|swiss|reddit|twitter|discord|bluesky)
                        key=$(get_api_key "$type")
                        if [ -z "$key" ]; then
                            echo "Skipping $type API (no key configured)"
                            continue
                        fi
                        ;;
                esac
                
                # Build command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli start-api"
                CMD="$CMD --api-dir $API_DIR/$version/$mode/$type"
                CMD="$CMD --cache-dir $CACHE_DIR"
                CMD="$CMD --version $version"
                CMD="$CMD --mode $mode"
                CMD="$CMD --type $type"
                CMD="$CMD --host $HOST"
                CMD="$CMD --port $port"
                CMD="$CMD --workers $WORKERS"
                CMD="$CMD --cache-ttl $CACHE_TTL"
                CMD="$CMD --timeout $TIMEOUT"
                CMD="$CMD --retries $RETRIES"
                CMD="$CMD --log-dir $LOG_DIR"
                
                if [ -n "$key" ]; then
                    CMD="$CMD --key $key"
                fi
                
                # Run command
                echo "Running: $CMD"
                $CMD || true  # Continue even if API fails
                
                # Generate API documentation
                if [ "$DOCS" = true ]; then
                    echo "Generating API documentation..."
                    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-api-docs \
                        --api-dir "$API_DIR/$version/$mode/$type" \
                        --docs-dir "$DOCS_DIR/$version" \
                        --output "$DOCS_DIR/$version/reference/$type.html"
                fi
                
                # Increment port
                ((port++))
            done
        done
    done
}

# Function to stop APIs
stop_apis() {
    echo "Stopping APIs..."
    
    for version in ${API_VERSIONS//,/ }; do
        for mode in ${API_MODES//,/ }; do
            for type in ${API_TYPES//,/ }; do
                echo "Stopping $type API $version in $mode mode..."
                
                # Build command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli stop-api"
                CMD="$CMD --api-dir $API_DIR/$version/$mode/$type"
                CMD="$CMD --version $version"
                CMD="$CMD --mode $mode"
                CMD="$CMD --type $type"
                CMD="$CMD --log-dir $LOG_DIR"
                
                # Run command
                echo "Running: $CMD"
                $CMD || true  # Continue even if stop fails
            done
        done
    done
}

# Function to monitor APIs
monitor_apis() {
    echo "Monitoring APIs..."
    
    # Build command
    CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli monitor-apis"
    CMD="$CMD --api-dir $API_DIR"
    CMD="$CMD --log-dir $LOG_DIR"
    
    # Run command
    echo "Running: $CMD"
    $CMD
}

# Function to test APIs
test_apis() {
    echo "Testing APIs..."
    
    for version in ${API_VERSIONS//,/ }; do
        for mode in ${API_MODES//,/ }; do
            for type in ${API_TYPES//,/ }; do
                echo "Testing $type API $version in $mode mode..."
                
                # Get API key if external service
                local key=""
                case "$type" in
                    chembl|pubchem|swiss|reddit|twitter|discord|bluesky)
                        key=$(get_api_key "$type")
                        if [ -z "$key" ]; then
                            echo "Skipping $type API (no key configured)"
                            continue
                        fi
                        ;;
                esac
                
                # Build command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli test-api"
                CMD="$CMD --api-dir $API_DIR/$version/$mode/$type"
                CMD="$CMD --version $version"
                CMD="$CMD --mode $mode"
                CMD="$CMD --type $type"
                CMD="$CMD --log-dir $LOG_DIR"
                
                if [ -n "$key" ]; then
                    CMD="$CMD --key $key"
                fi
                
                # Run command
                echo "Running: $CMD"
                $CMD || true  # Continue even if tests fail
                
                # Generate test report
                echo "Generating test report..."
                python -m binding_data_processor.processors.psychopharm.predictors.cli generate-test-report \
                    --api-dir "$API_DIR/$version/$mode/$type" \
                    --output "$API_DIR/$version/$mode/$type/tests/report.html"
            done
        done
    done
}

# Function to generate API documentation
generate_docs() {
    echo "Generating API documentation..."
    
    for version in ${API_VERSIONS//,/ }; do
        echo "Generating documentation for API $version..."
        
        for type in ${API_TYPES//,/ }; do
            echo "Generating $type API documentation..."
            
            # Generate API specifications
            echo "Generating API specifications..."
            case "$type" in
                rest)
                    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-openapi-spec \
                        --api-dir "$API_DIR/$version/production/$type" \
                        --output "$SPEC_DIR/$version/openapi/spec.yaml"
                    ;;
                graphql)
                    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-graphql-schema \
                        --api-dir "$API_DIR/$version/production/$type" \
                        --output "$SPEC_DIR/$version/graphql/schema.graphql"
                    ;;
                websocket)
                    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-asyncapi-spec \
                        --api-dir "$API_DIR/$version/production/$type" \
                        --output "$SPEC_DIR/$version/asyncapi/spec.yaml"
                    ;;
                *)
                    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-api-spec \
                        --api-dir "$API_DIR/$version/production/$type" \
                        --output "$SPEC_DIR/$version/$type/spec.yaml"
                    ;;
            esac
            
            # Generate reference documentation
            echo "Generating reference documentation..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-api-reference \
                --api-dir "$API_DIR/$version/production/$type" \
                --spec-dir "$SPEC_DIR/$version" \
                --output "$DOCS_DIR/$version/reference/$type.html"
            
            # Generate guides
            echo "Generating guides..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-api-guides \
                --api-dir "$API_DIR/$version/production/$type" \
                --output "$DOCS_DIR/$version/guides/$type"
            
            # Generate examples
            echo "Generating examples..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-api-examples \
                --api-dir "$API_DIR/$version/production/$type" \
                --output "$DOCS_DIR/$version/examples/$type"
        done
        
        # Generate combined documentation
        echo "Generating combined documentation..."
        python -m binding_data_processor.processors.psychopharm.predictors.cli generate-api-docs \
            --api-dir "$API_DIR/$version" \
            --spec-dir "$SPEC_DIR/$version" \
            --docs-dir "$DOCS_DIR/$version" \
            --output "$DOCS_DIR/$version/index.html"
    done
}

# Function to clean APIs
clean_apis() {
    echo "Cleaning APIs..."
    
    # Clean API directories
    rm -rf "$API_DIR"/*
    
    # Clean documentation
    rm -rf "$DOCS_DIR"/*
    
    # Clean specifications
    rm -rf "$SPEC_DIR"/*
    
    # Clean cache
    rm -rf "$CACHE_DIR"/*
}

# Function to backup APIs
backup_apis() {
    echo "Backing up APIs..."
    
    # Create backup directory
    BACKUP_DIR="backups/apis_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories and config
    cp -r "$API_DIR" "$BACKUP_DIR/"
    cp -r "$DOCS_DIR" "$BACKUP_DIR/"
    cp -r "$SPEC_DIR" "$BACKUP_DIR/"
    cp "$CONFIG_FILE" "$BACKUP_DIR/"
    
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

# Function to restore APIs
restore_apis() {
    echo "Restoring APIs..."
    
    # Find latest backup
    local LATEST_BACKUP=""
    if command_exists gpg; then
        LATEST_BACKUP=$(ls -t backups/apis_*.tar.gz.gpg 2>/dev/null | head -n1)
        if [ -n "$LATEST_BACKUP" ]; then
            echo "Decrypting backup: $LATEST_BACKUP"
            gpg --decrypt "$LATEST_BACKUP" | tar -xz
            BACKUP_DIR="${LATEST_BACKUP%.tar.gz.gpg}"
        fi
    fi
    
    if [ -z "$LATEST_BACKUP" ]; then
        LATEST_BACKUP=$(ls -t backups/apis_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$API_DIR" "$DOCS_DIR" "$SPEC_DIR"
        cp "$BACKUP_DIR/.env" "$CONFIG_FILE"
    else
        # Merge configs
        while IFS='=' read -r key value; do
            if [ -n "$key" ]; then
                set_api_key "${key%_API_KEY}" "$value"
            fi
        done < "$BACKUP_DIR/.env"
    fi
    
    cp -r "$BACKUP_DIR/api" ./
    cp -r "$BACKUP_DIR/docs/api" ./docs/
    cp -r "$BACKUP_DIR/specs/api" ./specs/
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "APIs restored from: $BACKUP_DIR"
}

# Function to configure APIs
configure_apis() {
    echo "Configuring APIs..."
    echo
    echo "Please enter your API keys (press Enter to skip):"
    echo
    
    for api in ${API_TYPES//,/ }; do
        # Only configure external services
        case "$api" in
            chembl|pubchem|swiss|reddit|twitter|discord|bluesky)
                # Get current key
                local current_key=$(get_api_key "$api")
                
                # Prompt for key
                local prompt="$api API key"
                if [ -n "$current_key" ]; then
                    prompt+=" (current: ${current_key:0:4}...${current_key: -4})"
                fi
                prompt+=": "
                
                read -p "$prompt" key
                
                # Update key if provided
                if [ -n "$key" ]; then
                    set_api_key "$api" "$key"
                    echo "Updated $api API key"
                fi
                
                echo
                ;;
        esac
    done
}

# Function to show API statistics
show_stats() {
    echo "API statistics:"
    echo
    
    echo "API files:"
    for version in ${API_VERSIONS//,/ }; do
        echo "$version API:"
        for mode in ${API_MODES//,/ }; do
            echo "  $mode mode:"
            for type in ${API_TYPES//,/ }; do
                echo "    $type API:"
                
                # Check key for external services
                case "$type" in
                    chembl|pubchem|swiss|reddit|twitter|discord|bluesky)
                        local key=$(get_api_key "$type")
                        if [ -n "$key" ]; then
                            echo "      Key: ${key:0:4}...${key: -4}"
                        else
                            echo "      Key: Not configured"
                        fi
                        ;;
                esac
                
                # Endpoint files
                echo "      Endpoints:"
                echo "        Files: $(find "$API_DIR/$version/$mode/$type/endpoints" -type f | wc -l) files"
                echo "        Size: $(du -sh "$API_DIR/$version/$mode/$type/endpoints" | cut -f1)"
                
                # Schema files
                echo "      Schemas:"
                echo "        Files: $(find "$API_DIR/$version/$mode/$type/schemas" -type f | wc -l) files"
                echo "        Size: $(du -sh "$API_DIR/$version/$mode/$type/schemas" | cut -f1)"
                
                # Count metrics by type
                case "$type" in
                    rest)
                        if [ -f "$API_DIR/$version/$mode/$type/endpoints/metrics.json" ]; then
                            echo "      Metrics:"
                            echo "        Routes: $(jq '.routes' "$API_DIR/$version/$mode/$type/endpoints/metrics.json")"
                            echo "        Controllers: $(jq '.controllers' "$API_DIR/$version/$mode/$type/endpoints/metrics.json")"
                        fi
                        ;;
                    graphql)
                        if [ -f "$API_DIR/$version/$mode/$type/endpoints/metrics.json" ]; then
                            echo "      Metrics:"
                            echo "        Types: $(jq '.types' "$API_DIR/$version/$mode/$type/endpoints/metrics.json")"
                            echo "        Resolvers: $(jq '.resolvers' "$API_DIR/$version/$mode/$type/endpoints/metrics.json")"
                        fi
                        ;;
                    websocket)
                        if [ -f "$API_DIR/$version/$mode/$type/endpoints/metrics.json" ]; then
                            echo "      Metrics:"
                            echo "        Channels: $(jq '.channels' "$API_DIR/$version/$mode/$type/endpoints/metrics.json")"
                            echo "        Handlers: $(jq '.handlers' "$API_DIR/$version/$mode/$type/endpoints/metrics.json")"
                        fi
                        ;;
                    chembl)
                        if [ -f "$API_DIR/$version/$mode/$type/endpoints/metrics.json" ]; then
                            echo "      Metrics:"
                            echo "        Queries: $(jq '.queries' "$API_DIR/$version/$mode/$type/endpoints/metrics.json")"
                            echo "        Results: $(jq '.results' "$API_DIR/$version/$mode/$type/endpoints/metrics.json")"
                        fi
                        ;;
                    pubchem)
                        if [ -f "$API_DIR/$version/$mode/$type/endpoints/metrics.json" ]; then
                            echo "      Metrics:"
                            echo "        Searches: $(jq '.searches' "$API_DIR/$version/$mode/$type/endpoints/metrics.json")"
                            echo "        Compounds: $(jq '.compounds' "$API_DIR/$version/$mode/$type/endpoints/metrics.json")"
                        fi
                        ;;
                    swiss)
                        if [ -f "$API_DIR/$version/$mode/$type/endpoints/metrics.json" ]; then
                            echo "      Metrics:"
                            echo "        Predictions: $(jq '.predictions' "$API_DIR/$version/$mode/$type/endpoints/metrics.json")"
                            echo "        Properties: $(jq '.properties' "$API_DIR/$version/$mode/$type/endpoints/metrics.json")"
                        fi
                        ;;
                    reddit|twitter|discord|bluesky)
                        if [ -f "$API_DIR/$version/$mode/$type/endpoints/metrics.json" ]; then
                            echo "      Metrics:"
                            echo "        Posts: $(jq '.posts' "$API_DIR/$version/$mode/$type/endpoints/metrics.json")"
                            echo "        Mentions: $(jq '.mentions' "$API_DIR/$version/$mode/$type/endpoints/metrics.json")"
                        fi
                        ;;
                esac
                
                # Show reports
                if [ -f "$API_DIR/$version/$mode/$type/tests/report.html" ]; then
                    echo "      Test report: $API_DIR/$version/$mode/$type/tests/report.html"
                fi
                echo
            done
        done
    done
    
    echo "Documentation files:"
    for version in ${API_VERSIONS//,/ }; do
        echo "$version documentation:"
        
        # Reference files
        echo "  Reference:"
        echo "    Files: $(find "$DOCS_DIR/$version/reference" -type f | wc -l) files"
        echo "    Size: $(du -sh "$DOCS_DIR/$version/reference" | cut -f1)"
        
        # Guide files
        echo "  Guides:"
        echo "    Files: $(find "$DOCS_DIR/$version/guides" -type f | wc -l) files"
        echo "    Size: $(du -sh "$DOCS_DIR/$version/guides" | cut -f1)"
        
        # Example files
        echo "  Examples:"
        echo "    Files: $(find "$DOCS_DIR/$version/examples" -type f | wc -l) files"
        echo "    Size: $(du -sh "$DOCS_DIR/$version/examples" | cut -f1)"
        echo
    done
    
    echo "Specification files:"
    for version in ${API_VERSIONS//,/ }; do
        echo "$version specifications:"
        
        # OpenAPI specs
        echo "  OpenAPI:"
        echo "    Files: $(find "$SPEC_DIR/$version/openapi" -type f | wc -l) files"
        echo "    Size: $(du -sh "$SPEC_DIR/$version/openapi" | cut -f1)"
        
        # GraphQL schemas
        echo "  GraphQL:"
        echo "    Files: $(find "$SPEC_DIR/$version/graphql" -type f | wc -l) files"
        echo "    Size: $(du -sh "$SPEC_DIR/$version/graphql" | cut -f1)"
        
        # AsyncAPI specs
        echo "  AsyncAPI:"
        echo "    Files: $(find "$SPEC_DIR/$version/asyncapi" -type f | wc -l) files"
        echo "    Size: $(du -sh "$SPEC_DIR/$version/asyncapi" | cut -f1)"
        
        # External API specs
        for type in chembl pubchem swiss reddit twitter discord bluesky; do
            if [ -d "$SPEC_DIR/$version/$type" ]; then
                echo "  $type:"
                echo "    Files: $(find "$SPEC_DIR/$version/$type" -type f | wc -l) files"
                echo "    Size: $(du -sh "$SPEC_DIR/$version/$type" | cut -f1)"
            fi
        done
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
echo "Managing APIs..."

# Create directory structure
create_dirs

# Configure APIs if no keys exist
if [ ! -f "$CONFIG_FILE" ]; then
    configure_apis
fi

# Start APIs if requested
if [ "$START" = true ]; then
    start_apis
fi

# Stop APIs if requested
if [ "$STOP" = true ]; then
    stop_apis
fi

# Monitor APIs if requested
if [ "$MONITOR" = true ]; then
    monitor_apis
fi

# Test APIs if requested
if [ "$TEST" = true ]; then
    test_apis
fi

# Generate documentation if requested
if [ "$DOCS" = true ]; then
    generate_docs
fi

# Clean APIs if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean APIs? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_apis
    fi
fi

# Backup APIs if requested
if [ "$BACKUP" = true ]; then
    backup_apis
fi

# Restore APIs if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore APIs? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_apis
    fi
fi

# Show statistics
show_stats

echo
echo "API management completed successfully!"
