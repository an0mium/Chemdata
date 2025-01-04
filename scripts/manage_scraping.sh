#!/bin/bash
# Script to manage web scraping, data extraction, and enrichment

# Exit on error
set -e

# Default values
CONFIG_FILE=".env"
DATA_DIR="data/web"
SCRAPE_DIR="scraping"
CACHE_DIR="cache/web"
LOG_DIR="logs"
SOURCE_TYPES="reddit,twitter,discord,bluesky,erowid,psychonautwiki,tripsit,pubmed,chembl,pubchem,swiss,patents"
CONTENT_TYPES="posts,comments,reports,articles,papers,structures,properties"
SCRAPE_TYPES="compounds,activity,safety,community"
SCRAPER_MODES="development,production,testing"
SEARCH_TERMS_FILE="search_terms.txt"
PROXY_FILE="proxies.txt"
BATCH_SIZE=100
RATE_LIMIT=2
DELAY=2
TIMEOUT=3600
RETRIES=3
PARALLEL=4
SCRAPE=false
MONITOR=false
TEST=false
VALIDATE=false
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
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --scrape-dir)
            SCRAPE_DIR="$2"
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
        --source-types)
            SOURCE_TYPES="$2"
            shift 2
            ;;
        --content-types)
            CONTENT_TYPES="$2"
            shift 2
            ;;
        --scrape-types)
            SCRAPE_TYPES="$2"
            shift 2
            ;;
        --scraper-modes)
            SCRAPER_MODES="$2"
            shift 2
            ;;
        --search-terms)
            SEARCH_TERMS_FILE="$2"
            shift 2
            ;;
        --proxies)
            PROXY_FILE="$2"
            shift 2
            ;;
        --batch-size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --rate-limit)
            RATE_LIMIT="$2"
            shift 2
            ;;
        --delay)
            DELAY="$2"
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
        --parallel)
            PARALLEL="$2"
            shift 2
            ;;
        --scrape)
            SCRAPE=true
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
    
    # Data directories
    for source in ${SOURCE_TYPES//,/ }; do
        for type in ${SCRAPE_TYPES//,/ }; do
            mkdir -p "$DATA_DIR/$source/$type"
        done
    done
    
    # Scraping directories
    for source in ${SOURCE_TYPES//,/ }; do
        for mode in ${SCRAPER_MODES//,/ }; do
            for content in ${CONTENT_TYPES//,/ }; do
                mkdir -p "$SCRAPE_DIR/$source/$mode/$content"/{raw,processed,enriched,validated,reports}
            done
        done
    done
    
    # Cache directory
    mkdir -p "$CACHE_DIR"
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to run scraping
run_scraping() {
    echo "Running scraping..."
    
    for mode in ${SCRAPER_MODES//,/ }; do
        echo "Running scrapers in $mode mode..."
        
        for source in ${SOURCE_TYPES//,/ }; do
            # Get API key if external service
            local key=""
            case "$source" in
                chembl|pubchem|swiss|reddit|twitter|discord|bluesky)
                    key=$(get_api_key "$source")
                    if [ -z "$key" ]; then
                        echo "Skipping $source scraper (no key configured)"
                        continue
                    fi
                    ;;
            esac
            
            # Run content-based scraping
            for content in ${CONTENT_TYPES//,/ }; do
                echo "Scraping $content from $source..."
                
                # Build command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli run-scraping"
                CMD="$CMD --data-dir $DATA_DIR"
                CMD="$CMD --scrape-dir $SCRAPE_DIR/$source/$mode/$content"
                CMD="$CMD --cache-dir $CACHE_DIR"
                CMD="$CMD --source $source"
                CMD="$CMD --mode $mode"
                CMD="$CMD --content $content"
                CMD="$CMD --search-terms $SEARCH_TERMS_FILE"
                CMD="$CMD --proxies $PROXY_FILE"
                CMD="$CMD --batch-size $BATCH_SIZE"
                CMD="$CMD --rate-limit $RATE_LIMIT"
                CMD="$CMD --delay $DELAY"
                CMD="$CMD --timeout $TIMEOUT"
                CMD="$CMD --retries $RETRIES"
                CMD="$CMD --parallel $PARALLEL"
                CMD="$CMD --log-dir $LOG_DIR"
                
                if [ -n "$key" ]; then
                    CMD="$CMD --key $key"
                fi
                
                # Run command
                echo "Running: $CMD"
                $CMD || true  # Continue even if scraping fails
                
                # Generate scraping report
                echo "Generating scraping report..."
                python -m binding_data_processor.processors.psychopharm.predictors.cli generate-scraping-report \
                    --scrape-dir "$SCRAPE_DIR/$source/$mode/$content" \
                    --output "$SCRAPE_DIR/$source/$mode/$content/reports/scraper.html"
            done
            
            # Run type-based scraping
            for type in ${SCRAPE_TYPES//,/ }; do
                echo "Scraping $type data from $source..."
                
                # Build command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli scrape-data"
                CMD="$CMD --source $source"
                CMD="$CMD --mode $mode"
                CMD="$CMD --type $type"
                CMD="$CMD --data-dir $DATA_DIR"
                CMD="$CMD --cache-dir $CACHE_DIR"
                CMD="$CMD --log-dir $LOG_DIR"
                CMD="$CMD --rate-limit $RATE_LIMIT"
                CMD="$CMD --timeout $TIMEOUT"
                CMD="$CMD --retries $RETRIES"
                CMD="$CMD --parallel $PARALLEL"
                
                if [ -n "$key" ]; then
                    CMD="$CMD --key $key"
                fi
                
                # Run command
                echo "Running: $CMD"
                $CMD || true  # Continue even if scraping fails
            done
        done
    done
}

# Function to monitor scraping
monitor_scraping() {
    echo "Monitoring scraping..."
    
    # Build command
    CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli monitor-scraping"
    CMD="$CMD --data-dir $DATA_DIR"
    CMD="$CMD --scrape-dir $SCRAPE_DIR"
    CMD="$CMD --log-dir $LOG_DIR"
    
    # Run command
    echo "Running: $CMD"
    $CMD
}

# Function to test scraping
test_scraping() {
    echo "Testing scraping..."
    
    for mode in ${SCRAPER_MODES//,/ }; do
        echo "Testing scrapers in $mode mode..."
        
        for source in ${SOURCE_TYPES//,/ }; do
            # Get API key if external service
            local key=""
            case "$source" in
                chembl|pubchem|swiss|reddit|twitter|discord|bluesky)
                    key=$(get_api_key "$source")
                    if [ -z "$key" ]; then
                        echo "Skipping $source scraper (no key configured)"
                        continue
                    fi
                    ;;
            esac
            
            # Test content-based scraping
            for content in ${CONTENT_TYPES//,/ }; do
                echo "Testing $content scraping from $source..."
                
                # Build command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli test-scraping"
                CMD="$CMD --scrape-dir $SCRAPE_DIR/$source/$mode/$content"
                CMD="$CMD --source $source"
                CMD="$CMD --mode $mode"
                CMD="$CMD --content $content"
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
                    --scrape-dir "$SCRAPE_DIR/$source/$mode/$content" \
                    --output "$SCRAPE_DIR/$source/$mode/$content/reports/test.html"
            done
            
            # Test type-based scraping
            for type in ${SCRAPE_TYPES//,/ }; do
                echo "Testing $type data scraping from $source..."
                
                # Build command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli test-data"
                CMD="$CMD --source $source"
                CMD="$CMD --mode $mode"
                CMD="$CMD --type $type"
                CMD="$CMD --data-dir $DATA_DIR"
                CMD="$CMD --log-dir $LOG_DIR"
                
                if [ -n "$key" ]; then
                    CMD="$CMD --key $key"
                fi
                
                # Run command
                echo "Running: $CMD"
                $CMD || true  # Continue even if tests fail
            done
        done
    done
}

# Function to validate scraped data
validate_data() {
    echo "Validating scraped data..."
    
    for mode in ${SCRAPER_MODES//,/ }; do
        echo "Validating data in $mode mode..."
        
        for source in ${SOURCE_TYPES//,/ }; do
            # Validate content-based data
            for content in ${CONTENT_TYPES//,/ }; do
                echo "Validating $content from $source..."
                
                # Build command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli validate-scraping"
                CMD="$CMD --scrape-dir $SCRAPE_DIR/$source/$mode/$content"
                CMD="$CMD --source $source"
                CMD="$CMD --mode $mode"
                CMD="$CMD --content $content"
                CMD="$CMD --log-dir $LOG_DIR"
                
                # Run command
                echo "Running: $CMD"
                $CMD || true  # Continue even if validation fails
                
                # Generate validation report
                echo "Generating validation report..."
                python -m binding_data_processor.processors.psychopharm.predictors.cli generate-validation-report \
                    --scrape-dir "$SCRAPE_DIR/$source/$mode/$content" \
                    --output "$SCRAPE_DIR/$source/$mode/$content/reports/validation.html"
            done
            
            # Validate type-based data
            for type in ${SCRAPE_TYPES//,/ }; do
                echo "Validating $type data from $source..."
                
                # Build command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli validate-data"
                CMD="$CMD --source $source"
                CMD="$CMD --mode $mode"
                CMD="$CMD --type $type"
                CMD="$CMD --data-dir $DATA_DIR"
                CMD="$CMD --log-dir $LOG_DIR"
                
                # Run command
                echo "Running: $CMD"
                $CMD || true  # Continue even if validation fails
            done
        done
    done
}

# Function to clean scraped data
clean_data() {
    echo "Cleaning scraped data..."
    
    # Clean data directories
    rm -rf "$DATA_DIR"/*
    
    # Clean scraping directories
    rm -rf "$SCRAPE_DIR"/*
    
    # Clean cache
    rm -rf "$CACHE_DIR"/*
}

# Function to backup scraped data
backup_data() {
    echo "Backing up scraped data..."
    
    # Create backup directory
    BACKUP_DIR="backups/web_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories and config
    cp -r "$DATA_DIR" "$BACKUP_DIR/"
    cp -r "$SCRAPE_DIR" "$BACKUP_DIR/"
    cp "$CONFIG_FILE" "$BACKUP_DIR/"
    cp "$SEARCH_TERMS_FILE" "$BACKUP_DIR/"
    cp "$PROXY_FILE" "$BACKUP_DIR/"
    
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

# Function to restore scraped data
restore_data() {
    echo "Restoring scraped data..."
    
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
        rm -rf "$DATA_DIR" "$SCRAPE_DIR"
        cp "$BACKUP_DIR/.env" "$CONFIG_FILE"
        cp "$BACKUP_DIR/search_terms.txt" "$SEARCH_TERMS_FILE"
        cp "$BACKUP_DIR/proxies.txt" "$PROXY_FILE"
    else
        # Merge configs
        while IFS='=' read -r key value; do
            if [ -n "$key" ]; then
                set_api_key "${key%_API_KEY}" "$value"
            fi
        done < "$BACKUP_DIR/.env"
    fi
    
    cp -r "$BACKUP_DIR/web" ./data/
    cp -r "$BACKUP_DIR/scraping" ./
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Data restored from: $BACKUP_DIR"
}

# Function to configure scraping
configure_scraping() {
    echo "Configuring scraping..."
    echo
    echo "Please enter your API keys (press Enter to skip):"
    echo
    
    for source in ${SOURCE_TYPES//,/ }; do
        # Only configure external services
        case "$source" in
            chembl|pubchem|swiss|reddit|twitter|discord|bluesky)
                # Get current key
                local current_key=$(get_api_key "$source")
                
                # Prompt for key
                local prompt="$source API key"
                if [ -n "$current_key" ]; then
                    prompt+=" (current: ${current_key:0:4}...${current_key: -4})"
                fi
                prompt+=": "
                
                read -p "$prompt" key
                
                # Update key if provided
                if [ -n "$key" ]; then
                    set_api_key "$source" "$key"
                    echo "Updated $source API key"
                fi
                
                echo
                ;;
        esac
    done
}

# Function to show scraping statistics
show_stats() {
    echo "Scraping statistics:"
    echo
    
    echo "Data files:"
    for mode in ${SCRAPER_MODES//,/ }; do
        echo "$mode mode:"
        for source in ${SOURCE_TYPES//,/ }; do
            echo "  $source data:"
            
            # Check key for external services
            case "$source" in
                chembl|pubchem|swiss|reddit|twitter|discord|bluesky)
                    local key=$(get_api_key "$source")
                    if [ -n "$key" ]; then
                        echo "    Key: ${key:0:4}...${key: -4}"
                    else
                        echo "    Key: Not configured"
                    fi
                    ;;
            esac
            
            # Show type-based stats
            for type in ${SCRAPE_TYPES//,/ }; do
                echo "    $type data:"
                echo "      Files: $(find "$DATA_DIR/$source/$type" -type f | wc -l) files"
                echo "      Size: $(du -sh "$DATA_DIR/$source/$type" | cut -f1)"
                
                # Count entries by type
                case "$type" in
                    compounds)
                        if [ -f "$DATA_DIR/$source/$type/compounds.json" ]; then
                            echo "      Compounds: $(jq '.compounds | length' "$DATA_DIR/$source/$type/compounds.json")"
                        fi
                        ;;
                    activity)
                        if [ -f "$DATA_DIR/$source/$type/activities.json" ]; then
                            echo "      Activities: $(jq '.activities | length' "$DATA_DIR/$source/$type/activities.json")"
                        fi
                        ;;
                    safety)
                        if [ -f "$DATA_DIR/$source/$type/reports.json" ]; then
                            echo "      Reports: $(jq '.reports | length' "$DATA_DIR/$source/$type/reports.json")"
                        fi
                        ;;
                    community)
                        if [ -f "$DATA_DIR/$source/$type/posts.json" ]; then
                            echo "      Posts: $(jq '.posts | length' "$DATA_DIR/$source/$type/posts.json")"
                        fi
                        ;;
                esac
            done
            
            # Show content-based stats
            for content in ${CONTENT_TYPES//,/ }; do
                echo "    $content data:"
                
                # Raw data
                echo "      Raw data:"
                echo "        Files: $(find "$SCRAPE_DIR/$source/$mode/$content/raw" -type f | wc -l) files"
                echo "        Size: $(du -sh "$SCRAPE_DIR/$source/$mode/$content/raw" | cut -f1)"
                
                # Count entries by content type
                case "$content" in
                    posts)
                        if [ -f "$SCRAPE_DIR/$source/$mode/$content/raw/posts.json" ]; then
                            echo "        Posts: $(jq '.posts | length' "$SCRAPE_DIR/$source/$mode/$content/raw/posts.json")"
                            echo "        Authors: $(jq '.authors | length' "$SCRAPE_DIR/$source/$mode/$content/raw/posts.json")"
                        fi
                        ;;
                    comments)
                        if [ -f "$SCRAPE_DIR/$source/$mode/$content/raw/comments.json" ]; then
                            echo "        Comments: $(jq '.comments | length' "$SCRAPE_DIR/$source/$mode/$content/raw/comments.json")"
                            echo "        Threads: $(jq '.threads | length' "$SCRAPE_DIR/$source/$mode/$content/raw/comments.json")"
                        fi
                        ;;
                    reports)
                        if [ -f "$SCRAPE_DIR/$source/$mode/$content/raw/reports.json" ]; then
                            echo "        Reports: $(jq '.reports | length' "$SCRAPE_DIR/$source/$mode/$content/raw/reports.json")"
                            echo "        Substances: $(jq '.substances | length' "$SCRAPE_DIR/$source/$mode/$content/raw/reports.json")"
                        fi
                        ;;
                    articles)
                        if [ -f "$SCRAPE_DIR/$source/$mode/$content/raw/articles.json" ]; then
                            echo "        Articles: $(jq '.articles | length' "$SCRAPE_DIR/$source/$mode/$content/raw/articles.json")"
                            echo "        Sources: $(jq '.sources | length' "$SCRAPE_DIR/$source/$mode/$content/raw/articles.json")"
                        fi
                        ;;
                    papers)
                        if [ -f "$SCRAPE_DIR/$source/$mode/$content/raw/papers.json" ]; then
                            echo "        Papers: $(jq '.papers | length' "$SCRAPE_DIR/$source/$mode/$content/raw/papers.json")"
                            echo "        Citations: $(jq '.citations | length' "$SCRAPE_DIR/$source/$mode/$content/raw/papers.json")"
                        fi
                        ;;
                    structures)
                        if [ -f "$SCRAPE_DIR/$source/$mode/$content/raw/structures.json" ]; then
                            echo "        Structures: $(jq '.structures | length' "$SCRAPE_DIR/$source/$mode/$content/raw/structures.json")"
                            echo "        Formats: $(jq '.formats | length' "$SCRAPE_DIR/$source/$mode/$content/raw/structures.json")"
                        fi
                        ;;
                    properties)
                        if [ -f "$SCRAPE_DIR/$source/$mode/$content/raw/properties.json" ]; then
                            echo "        Properties: $(jq '.properties | length' "$SCRAPE_DIR/$source/$mode/$content/raw/properties.json")"
                            echo "        Compounds: $(jq '.compounds | length' "$SCRAPE_DIR/$source/$mode/$content/raw/properties.json")"
                        fi
                        ;;
                esac
                
                # Processed data
                echo "      Processed data:"
                echo "        Files: $(find "$SCRAPE_DIR/$source/$mode/$content/processed" -type f | wc -l) files"
                echo "        Size: $(du -sh "$SCRAPE_DIR/$source/$mode/$content/processed" | cut -f1)"
                
                # Enriched data
                echo "      Enriched data:"
                echo "        Files: $(find "$SCRAPE_DIR/$source/$mode/$content/enriched" -type f | wc -l) files"
                echo "        Size: $(du -sh "$SCRAPE_DIR/$source/$mode/$content/enriched" | cut -f1)"
                
                # Validated data
                echo "      Validated data:"
                echo "        Files: $(find "$SCRAPE_DIR/$source/$mode/$content/validated" -type f | wc -l) files"
                echo "        Size: $(du -sh "$SCRAPE_DIR/$source/$mode/$content/validated" | cut -f1)"
                
                # Show reports
                if [ -f "$SCRAPE_DIR/$source/$mode/$content/reports/scraper.html" ]; then
                    echo "      Scraping report: $SCRAPE_DIR/$source/$mode/$content/reports/scraper.html"
                fi
                if [ -f "$SCRAPE_DIR/$source/$mode/$content/reports/test.html" ]; then
                    echo "      Test report: $SCRAPE_DIR/$source/$mode/$content/reports/test.html"
                fi
                if [ -f "$SCRAPE_DIR/$source/$mode/$content/reports/validation.html" ]; then
                    echo "      Validation report: $SCRAPE_DIR/$source/$mode/$content/reports/validation.html"
                fi
            done
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
echo "Managing web scraping..."

# Create directory structure
create_dirs

# Configure scraping if no keys exist
if [ ! -f "$CONFIG_FILE" ]; then
    configure_scraping
fi

# Run scraping if requested
if [ "$SCRAPE" = true ]; then
    run_scraping
fi

# Monitor scraping if requested
if [ "$MONITOR" = true ]; then
    monitor_scraping
fi

# Test scraping if requested
if [ "$TEST" = true ]; then
    test_scraping
fi

# Validate data if requested
if [ "$VALIDATE" = true ]; then
    validate_data
fi

# Clean data if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean scraped data? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_data
    fi
fi

# Backup data if requested
if [ "$BACKUP" = true ]; then
    backup_data
fi

# Restore data if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore scraped data? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_data
    fi
fi

# Show statistics
show_stats

echo
echo "Web scraping management completed successfully!"
