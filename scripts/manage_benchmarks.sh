#!/bin/bash
# Script to manage performance benchmarks, profiling, and analysis

# Exit on error
set -e

# Default values
DATA_DIR="data/benchmarks"
BENCHMARK_DIR="benchmarks"
PROFILE_DIR="profiles"
REPORT_DIR="reports/benchmarks"
CACHE_DIR="cache/benchmarks"
LOG_DIR="logs"
BENCHMARK_TYPES="pipeline,ml,web,api,services,scraping,processing,analysis,prediction,validation,export,enrichment,social,community,structure"
BENCHMARK_MODES="speed,memory,cpu,gpu,development,production,testing,accuracy,reliability"
PROFILE_TYPES="cprofile,memory_profiler,line_profiler,memory,cpu,io,network,database,api,ml,web"
ML_METRICS="accuracy,precision,recall,f1,auc,mae,mse,rmse"
WEB_METRICS="success_rate,latency,throughput,error_rate,coverage"
DATA_METRICS="completeness,accuracy,consistency,timeliness"
BATCH_SIZES="100,500,1000,5000"
PARALLEL_LEVELS="1,2,4,8"
BENCHMARK_PATTERN="benchmark_*.py"
ITERATIONS=3
WARMUP=1
TIMEOUT=3600
RETRIES=3
RUN=false
MONITOR=false
COMPARE=false
PROFILE=false
REPORT=false
CLEAN=false
BACKUP=false
RESTORE=false
GPU=false
FORCE=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --benchmark-dir)
            BENCHMARK_DIR="$2"
            shift 2
            ;;
        --profile-dir)
            PROFILE_DIR="$2"
            shift 2
            ;;
        --report-dir)
            REPORT_DIR="$2"
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
        --benchmark-types)
            BENCHMARK_TYPES="$2"
            shift 2
            ;;
        --benchmark-modes)
            BENCHMARK_MODES="$2"
            shift 2
            ;;
        --profile-types)
            PROFILE_TYPES="$2"
            shift 2
            ;;
        --ml-metrics)
            ML_METRICS="$2"
            shift 2
            ;;
        --web-metrics)
            WEB_METRICS="$2"
            shift 2
            ;;
        --data-metrics)
            DATA_METRICS="$2"
            shift 2
            ;;
        --batch-sizes)
            BATCH_SIZES="$2"
            shift 2
            ;;
        --parallel-levels)
            PARALLEL_LEVELS="$2"
            shift 2
            ;;
        --pattern)
            BENCHMARK_PATTERN="$2"
            shift 2
            ;;
        --iterations)
            ITERATIONS="$2"
            shift 2
            ;;
        --warmup)
            WARMUP="$2"
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
        --compare)
            COMPARE=true
            shift
            ;;
        --profile)
            PROFILE=true
            shift
            ;;
        --report)
            REPORT=true
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
        --gpu)
            GPU=true
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
    
    # Benchmark directories
    for type in ${BENCHMARK_TYPES//,/ }; do
        for mode in ${BENCHMARK_MODES//,/ }; do
            mkdir -p "$BENCHMARK_DIR/$type/$mode"/{data,profiles,reports,metrics}
            mkdir -p "$DATA_DIR/$type/$mode"/{raw,processed,reports,metrics}
        done
    done
    
    # Profile directories
    for type in ${PROFILE_TYPES//,/ }; do
        mkdir -p "$PROFILE_DIR/$type"/{snapshots,flamegraphs,reports,metrics}
    done
    
    # Report directories
    mkdir -p "$REPORT_DIR"/{html,json,profiles,comparisons,metrics}
    
    # Metrics directories
    mkdir -p "$REPORT_DIR/metrics"/{ml,web,data,performance}
    
    # Cache directory
    mkdir -p "$CACHE_DIR"
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to run ML benchmarks
run_ml_benchmarks() {
    echo "Running ML benchmarks..."
    
    for mode in ${BENCHMARK_MODES//,/ }; do
        echo "Running ML benchmarks in $mode mode..."
        
        # Benchmark different ML components
        for component in binding activity toxicity abuse psychoactive nootropic; do
            echo "Benchmarking $component predictions..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli benchmark-ml"
            CMD="$CMD --data-dir $DATA_DIR/ml/$mode"
            CMD="$CMD --component $component"
            CMD="$CMD --mode $mode"
            CMD="$CMD --metrics $ML_METRICS"
            CMD="$CMD --iterations $ITERATIONS"
            CMD="$CMD --log-dir $LOG_DIR"
            
            if [ "$GPU" = true ]; then
                CMD="$CMD --gpu"
            fi
            
            # Run benchmark
            echo "Running: $CMD"
            $CMD || true
            
            # Generate ML report
            echo "Generating ML report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-ml-report \
                --data-dir "$DATA_DIR/ml/$mode" \
                --component "$component" \
                --output "$DATA_DIR/ml/$mode/reports/${component}_report.html"
        done
    done
}

# Function to run web benchmarks
run_web_benchmarks() {
    echo "Running web benchmarks..."
    
    for mode in ${BENCHMARK_MODES//,/ }; do
        echo "Running web benchmarks in $mode mode..."
        
        # Benchmark different web components
        for component in scraping enrichment social community structure; do
            echo "Benchmarking $component processing..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli benchmark-web"
            CMD="$CMD --data-dir $DATA_DIR/web/$mode"
            CMD="$CMD --component $component"
            CMD="$CMD --mode $mode"
            CMD="$CMD --metrics $WEB_METRICS"
            CMD="$CMD --iterations $ITERATIONS"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run benchmark
            echo "Running: $CMD"
            $CMD || true
            
            # Generate web report
            echo "Generating web report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-web-report \
                --data-dir "$DATA_DIR/web/$mode" \
                --component "$component" \
                --output "$DATA_DIR/web/$mode/reports/${component}_report.html"
        done
    done
}

# Function to run data quality benchmarks
run_data_benchmarks() {
    echo "Running data quality benchmarks..."
    
    for mode in ${BENCHMARK_MODES//,/ }; do
        echo "Running data benchmarks in $mode mode..."
        
        # Benchmark different data components
        for component in bindingdb chembl pubchem swiss community social; do
            echo "Benchmarking $component data quality..."
            
            # Build command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli benchmark-data"
            CMD="$CMD --data-dir $DATA_DIR/data/$mode"
            CMD="$CMD --component $component"
            CMD="$CMD --mode $mode"
            CMD="$CMD --metrics $DATA_METRICS"
            CMD="$CMD --iterations $ITERATIONS"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run benchmark
            echo "Running: $CMD"
            $CMD || true
            
            # Generate data report
            echo "Generating data report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-data-report \
                --data-dir "$DATA_DIR/data/$mode" \
                --component "$component" \
                --output "$DATA_DIR/data/$mode/reports/${component}_report.html"
        done
    done
}

# Function to run benchmarks
run_benchmarks() {
    echo "Running benchmarks..."
    
    # Run ML benchmarks
    run_ml_benchmarks
    
    # Run web benchmarks
    run_web_benchmarks
    
    # Run data benchmarks
    run_data_benchmarks
    
    # Run standard benchmarks
    for mode in ${BENCHMARK_MODES//,/ }; do
        echo "Running benchmarks in $mode mode..."
        
        for type in ${BENCHMARK_TYPES//,/ }; do
            echo "Running $type benchmarks..."
            
            # Build pytest command for benchmark discovery
            PYTEST_CMD="pytest"
            PYTEST_CMD="$PYTEST_CMD $BENCHMARK_DIR/$type/$BENCHMARK_PATTERN"
            PYTEST_CMD="$PYTEST_CMD --benchmark-only"
            PYTEST_CMD="$PYTEST_CMD --benchmark-json=$REPORT_DIR/json/$type-$mode.json"
            
            if [ "$GPU" = true ]; then
                PYTEST_CMD="$PYTEST_CMD --gpu"
            fi
            
            # Run pytest benchmarks
            echo "Running pytest benchmarks: $PYTEST_CMD"
            $PYTEST_CMD || true
            
            # Run batch size benchmarks
            for size in ${BATCH_SIZES//,/ }; do
                echo "Testing batch size $size..."
                
                # Run parallel level benchmarks
                for level in ${PARALLEL_LEVELS//,/ }; do
                    echo "Testing parallel level $level..."
                    
                    # Run iterations
                    for ((i=1; i<=$ITERATIONS; i++)); do
                        echo "Running iteration $i..."
                        
                        # Warmup if enabled
                        if [ "$WARMUP" = true ]; then
                            echo "Running warmup..."
                            
                            # Build warmup command
                            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli run-benchmark"
                            CMD="$CMD --data-dir $DATA_DIR/$type/$mode"
                            CMD="$CMD --benchmark-dir $BENCHMARK_DIR/$type/$mode"
                            CMD="$CMD --cache-dir $CACHE_DIR"
                            CMD="$CMD --type $type"
                            CMD="$CMD --mode $mode"
                            CMD="$CMD --batch-size $size"
                            CMD="$CMD --parallel $level"
                            CMD="$CMD --warmup true"
                            CMD="$CMD --timeout $TIMEOUT"
                            CMD="$CMD --retries $RETRIES"
                            CMD="$CMD --log-dir $LOG_DIR"
                            
                            if [ "$GPU" = true ]; then
                                CMD="$CMD --gpu"
                            fi
                            
                            # Run warmup
                            echo "Running: $CMD"
                            $CMD || true
                        fi
                        
                        # Build benchmark command
                        CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli run-benchmark"
                        CMD="$CMD --data-dir $DATA_DIR/$type/$mode"
                        CMD="$CMD --benchmark-dir $BENCHMARK_DIR/$type/$mode"
                        CMD="$CMD --cache-dir $CACHE_DIR"
                        CMD="$CMD --type $type"
                        CMD="$CMD --mode $mode"
                        CMD="$CMD --batch-size $size"
                        CMD="$CMD --parallel $level"
                        CMD="$CMD --warmup false"
                        CMD="$CMD --timeout $TIMEOUT"
                        CMD="$CMD --retries $RETRIES"
                        CMD="$CMD --log-dir $LOG_DIR"
                        
                        if [ "$GPU" = true ]; then
                            CMD="$CMD --gpu"
                        fi
                        
                        # Run benchmark
                        echo "Running: $CMD"
                        $CMD || true
                    done
                    
                    # Generate benchmark report
                    echo "Generating benchmark report..."
                    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-benchmark-report \
                        --data-dir "$DATA_DIR/$type/$mode" \
                        --benchmark-dir "$BENCHMARK_DIR/$type/$mode" \
                        --batch-size "$size" \
                        --parallel "$level" \
                        --output "$DATA_DIR/$type/$mode/reports/benchmark_${size}_${level}.html"
                done
            done
        done
    done
    
    # Generate combined reports
    echo "Generating combined reports..."
    python -m pytest_benchmark.cli compare \
        "$REPORT_DIR/json"/*.json \
        --csv "$REPORT_DIR/report.csv" \
        --histogram "$REPORT_DIR/html/histogram.svg" \
        --sort name
}

# Function to monitor benchmarks
monitor_benchmarks() {
    echo "Monitoring benchmarks..."
    
    # Build command
    CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli monitor-benchmarks"
    CMD="$CMD --data-dir $DATA_DIR"
    CMD="$CMD --benchmark-dir $BENCHMARK_DIR"
    CMD="$CMD --log-dir $LOG_DIR"
    
    # Run command
    echo "Running: $CMD"
    $CMD
}

# Function to run profiling
run_profiling() {
    echo "Running profiling..."
    
    for type in ${PROFILE_TYPES//,/ }; do
        echo "Running $type profiling..."
        
        # Build profile command
        CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli run-profile"
        CMD="$CMD --profile-dir $PROFILE_DIR/$type"
        CMD="$CMD --type $type"
        CMD="$CMD --timeout $TIMEOUT"
        CMD="$CMD --log-dir $LOG_DIR"
        
        # Run profiling
        echo "Running: $CMD"
        $CMD || true
        
        # Generate profile report
        echo "Generating profile report..."
        python -m binding_data_processor.processors.psychopharm.predictors.cli generate-profile-report \
            --profile-dir "$PROFILE_DIR/$type" \
            --output "$PROFILE_DIR/$type/reports/profile.html"
        
        # Generate flamegraph
        if command_exists flamegraph; then
            echo "Generating flamegraph..."
            flamegraph.pl "$PROFILE_DIR/$type/snapshots/profile.folded" > "$PROFILE_DIR/$type/flamegraphs/profile.svg"
        fi
    done
}

# Function to compare benchmarks
compare_benchmarks() {
    echo "Comparing benchmarks..."
    
    # Compare ML benchmarks
    echo "Comparing ML benchmarks..."
    python -m binding_data_processor.processors.psychopharm.predictors.cli compare-ml-benchmarks \
        --data-dir "$DATA_DIR/ml" \
        --output "$REPORT_DIR/comparisons/ml_comparison.html"
    
    # Compare web benchmarks
    echo "Comparing web benchmarks..."
    python -m binding_data_processor.processors.psychopharm.predictors.cli compare-web-benchmarks \
        --data-dir "$DATA_DIR/web" \
        --output "$REPORT_DIR/comparisons/web_comparison.html"
    
    # Compare data benchmarks
    echo "Comparing data benchmarks..."
    python -m binding_data_processor.processors.psychopharm.predictors.cli compare-data-benchmarks \
        --data-dir "$DATA_DIR/data" \
        --output "$REPORT_DIR/comparisons/data_comparison.html"
    
    # Compare standard benchmarks
    for type in ${BENCHMARK_TYPES//,/ }; do
        for mode in ${BENCHMARK_MODES//,/ }; do
            echo "Comparing $mode benchmarks for $type..."
            
            # Find benchmark files
            BENCHMARK_FILES=()
            if [ -f "$REPORT_DIR/json/$type-$mode.json" ]; then
                BENCHMARK_FILES+=("$REPORT_DIR/json/$type-$mode.json")
            fi
            if [ -f "$BENCHMARK_DIR/$type/$mode/data/benchmarks.json" ]; then
                BENCHMARK_FILES+=("$BENCHMARK_DIR/$type/$mode/data/benchmarks.json")
            fi
            
            if [ ${#BENCHMARK_FILES[@]} -gt 0 ]; then
                # Compare benchmarks
                python -m pytest_benchmark.cli compare \
                    "${BENCHMARK_FILES[@]}" \
                    --csv "$REPORT_DIR/comparisons/$type-$mode.csv" \
                    --histogram "$REPORT_DIR/html/$type-$mode-comparison.svg" \
                    --sort name
                
                # Build comparison command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli compare-benchmarks"
                CMD="$CMD --data-dir $DATA_DIR/$type/$mode"
                CMD="$CMD --benchmark-dir $BENCHMARK_DIR/$type/$mode"
                CMD="$CMD --type $type"
                CMD="$CMD --mode $mode"
                CMD="$CMD --log-dir $LOG_DIR"
                
                # Run command
                echo "Running: $CMD"
                $CMD || true
                
                # Generate comparison report
                echo "Generating comparison report..."
                python -m binding_data_processor.processors.psychopharm.predictors.cli generate-comparison-report \
                    --data-dir "$DATA_DIR/$type/$mode" \
                    --benchmark-dir "$BENCHMARK_DIR/$type/$mode" \
                    --output "$DATA_DIR/$type/$mode/reports/comparison.html"
            else
                echo "No benchmark files found for $type $mode"
            fi
        done
    done
    
    # Generate combined comparison
    echo "Generating combined comparison..."
    python -m pytest_benchmark.cli compare \
        "$REPORT_DIR/json"/*.json \
        --csv "$REPORT_DIR/comparisons/combined.csv" \
        --histogram "$REPORT_DIR/html/combined-comparison.svg" \
        --sort name
}

# Function to generate reports
generate_reports() {
    echo "Generating reports..."
    
    # Generate ML summary
    echo "Generating ML summary..."
    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-ml-summary \
        --data-dir "$DATA_DIR/ml" \
        --output "$REPORT_DIR/metrics/ml/summary.html"
    
    # Generate web summary
    echo "Generating web summary..."
    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-web-summary \
        --data-dir "$DATA_DIR/web" \
        --output "$REPORT_DIR/metrics/web/summary.html"
    
    # Generate data summary
    echo "Generating data summary..."
    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-data-summary \
        --data-dir "$DATA_DIR/data" \
        --output "$REPORT_DIR/metrics/data/summary.html"
    
    # Generate benchmark summary
    echo "Generating benchmark summary..."
    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-benchmark-summary \
        --data-dir "$DATA_DIR" \
        --benchmark-dir "$BENCHMARK_DIR" \
        --output "$DATA_DIR/benchmark_summary.html"
    
    # Generate profile summary
    echo "Generating profile summary..."
    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-profile-summary \
        --profile-dir "$PROFILE_DIR" \
        --output "$PROFILE_DIR/profile_summary.html"
    
    # Generate combined report
    echo "Generating combined report..."
    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-performance-report \
        --data-dir "$DATA_DIR" \
        --benchmark-dir "$BENCHMARK_DIR" \
        --profile-dir "$PROFILE_DIR" \
        --output "performance_report.html"
}

# Function to clean benchmarks
clean_benchmarks() {
    echo "Cleaning benchmarks..."
    
    # Clean benchmark directories
    rm -rf "$BENCHMARK_DIR"/*
    rm -rf "$DATA_DIR"/*
    
    # Clean profile directories
    rm -rf "$PROFILE_DIR"/*
    
    # Clean report directories
    rm -rf "$REPORT_DIR"/*
    
    # Clean cache
    rm -rf "$CACHE_DIR"/*
    
    # Clean pytest cache
    find . -type d -name ".pytest_cache" -exec rm -rf {} +
    find . -type d -name "__pycache__" -exec rm -rf {} +
}

# Function to backup benchmarks
backup_benchmarks() {
    echo "Backing up benchmarks..."
    
    # Create backup directory
    BACKUP_DIR="backups/benchmarks_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$BENCHMARK_DIR" "$BACKUP_DIR/"
    cp -r "$DATA_DIR" "$BACKUP_DIR/"
    cp -r "$PROFILE_DIR" "$BACKUP_DIR/"
    cp -r "$REPORT_DIR" "$BACKUP_DIR/"
    
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

# Function to restore benchmarks
restore_benchmarks() {
    echo "Restoring benchmarks..."
    
    # Find latest backup
    local LATEST_BACKUP=""
    if command_exists gpg; then
        LATEST_BACKUP=$(ls -t backups/benchmarks_*.tar.gz.gpg 2>/dev/null | head -n1)
        if [ -n "$LATEST_BACKUP" ]; then
            echo "Decrypting backup: $LATEST_BACKUP"
            gpg --decrypt "$LATEST_BACKUP" | tar -xz
            BACKUP_DIR="${LATEST_BACKUP%.tar.gz.gpg}"
        fi
    fi
    
    if [ -z "$LATEST_BACKUP" ]; then
        LATEST_BACKUP=$(ls -t backups/benchmarks_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$BENCHMARK_DIR" "$DATA_DIR" "$PROFILE_DIR" "$REPORT_DIR"
    fi
    
    cp -r "$BACKUP_DIR/benchmarks" ./
    cp -r "$BACKUP_DIR/data/benchmarks" ./data/
    cp -r "$BACKUP_DIR/profiles" ./
    cp -r "$BACKUP_DIR/reports/benchmarks" ./reports/
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Benchmarks restored from: $BACKUP_DIR"
}

# Function to show benchmark statistics
show_stats() {
    echo "Benchmark statistics:"
    echo
    
    # Show ML metrics
    echo "ML metrics:"
    for mode in ${BENCHMARK_MODES//,/ }; do
        echo "$mode mode:"
        for metric in ${ML_METRICS//,/ }; do
            if [ -f "$REPORT_DIR/metrics/ml/${mode}_${metric}.json" ]; then
                echo "  $metric: $(jq -r '.value' "$REPORT_DIR/metrics/ml/${mode}_${metric}.json")"
            fi
        done
    done
    echo
    
    # Show web metrics
    echo "Web metrics:"
    for mode in ${BENCHMARK_MODES//,/ }; do
        echo "$mode mode:"
        for metric in ${WEB_METRICS//,/ }; do
            if [ -f "$REPORT_DIR/metrics/web/${mode}_${metric}.json" ]; then
                echo "  $metric: $(jq -r '.value' "$REPORT_DIR/metrics/web/${mode}_${metric}.json")"
            fi
        done
    done
    echo
    
    # Show data metrics
    echo "Data metrics:"
    for mode in ${BENCHMARK_MODES//,/ }; do
        echo "$mode mode:"
        for metric in ${DATA_METRICS//,/ }; do
            if [ -f "$REPORT_DIR/metrics/data/${mode}_${metric}.json" ]; then
                echo "  $metric: $(jq -r '.value' "$REPORT_DIR/metrics/data/${mode}_${metric}.json")"
            fi
        done
    done
    echo
    
    # Show benchmark files
    echo "Benchmark files:"
    for mode in ${BENCHMARK_MODES//,/ }; do
        echo "$mode mode:"
        for type in ${BENCHMARK_TYPES//,/ }; do
            echo "  $type benchmarks:"
            echo "    Test files: $(find "$BENCHMARK_DIR/$type" -type f -name "$BENCHMARK_PATTERN" | wc -l) files"
            
            # Raw data
            echo "    Raw data:"
            echo "      Files: $(find "$DATA_DIR/$type/$mode/raw" -type f | wc -l) files"
            echo "      Size: $(du -sh "$DATA_DIR/$type/$mode/raw" | cut -f1)"
            
            # Processed data
            echo "    Processed data:"
            echo "      Files: $(find "$DATA_DIR/$type/$mode/processed" -type f | wc -l) files"
            echo "      Size: $(du -sh "$DATA_DIR/$type/$mode/processed" | cut -f1)"
            
            # Show batch size stats
            for size in ${BATCH_SIZES//,/ }; do
                echo "    Batch size $size:"
                
                # Show parallel level stats
                for level in ${PARALLEL_LEVELS//,/ }; do
                    if [ -f "$DATA_DIR/$type/$mode/processed/benchmark_${size}_${level}.json" ]; then
                        echo "      Parallel level $level:"
                        echo "        Runtime: $(jq '.runtime' "$DATA_DIR/$type/$mode/processed/benchmark_${size}_${level}.json") seconds"
                        echo "        Memory: $(jq '.memory' "$DATA_DIR/$type/$mode/processed/benchmark_${size}_${level}.json") MB"
                        echo "        Throughput: $(jq '.throughput' "$DATA_DIR/$type/$mode/processed/benchmark_${size}_${level}.json") items/sec"
                    fi
                done
            done
            
            # Show pytest results
            if [ -f "$REPORT_DIR/json/$type-$mode.json" ]; then
                echo "    Pytest benchmarks:"
                echo "      Tests: $(jq '.benchmarks | length' "$REPORT_DIR/json/$type-$mode.json") benchmarks"
                echo "      Last run: $(jq -r '.datetime' "$REPORT_DIR/json/$type-$mode.json")"
            fi
            
            # Show reports
            if [ -f "$DATA_DIR/$type/$mode/reports/benchmark.html" ]; then
                echo "    Benchmark report: $DATA_DIR/$type/$mode/reports/benchmark.html"
            fi
            if [ -f "$DATA_DIR/$type/$mode/reports/comparison.html" ]; then
                echo "    Comparison report: $DATA_DIR/$type/$mode/reports/comparison.html"
            fi
            echo
        done
    done
    
    # Show profile files
    echo "Profile files:"
    for type in ${PROFILE_TYPES//,/ }; do
        echo "$type profiles:"
        
        # Snapshot files
        echo "  Snapshots:"
        echo "    Files: $(find "$PROFILE_DIR/$type/snapshots" -type f | wc -l) files"
        echo "    Size: $(du -sh "$PROFILE_DIR/$type/snapshots" | cut -f1)"
        
        # Flamegraph files
        echo "  Flamegraphs:"
        echo "    Files: $(find "$PROFILE_DIR/$type/flamegraphs" -type f | wc -l) files"
        echo "    Size: $(du -sh "$PROFILE_DIR/$type/flamegraphs" | cut -f1)"
        
        # Show reports
        if [ -f "$PROFILE_DIR/$type/reports/profile.html" ]; then
            echo "  Profile report: $PROFILE_DIR/$type/reports/profile.html"
        fi
        echo
    done
    
    # Show reports
    echo "Reports:"
    echo "  JSON: $(find "$REPORT_DIR/json" -type f -name "*.json" | wc -l) files"
    echo "  HTML: $(find "$REPORT_DIR/html" -type f | wc -l) files"
    echo "  Comparisons: $(find "$REPORT_DIR/comparisons" -type f | wc -l) files"
    echo "  Profiles: $(find "$REPORT_DIR/profiles" -type f | wc -l) files"
    echo "  Metrics: $(find "$REPORT_DIR/metrics" -type f | wc -l) files"
    echo
    
    # Show summary reports
    if [ -f "$DATA_DIR/benchmark_summary.html" ]; then
        echo "Benchmark summary: $DATA_DIR/benchmark_summary.html"
    fi
    if [ -f "$PROFILE_DIR/profile_summary.html" ]; then
        echo "Profile summary: $PROFILE_DIR/profile_summary.html"
    fi
    if [ -f "performance_report.html" ]; then
        echo "Performance report: performance_report.html"
    fi
    
    if [ -f "$REPORT_DIR/comparisons/combined.csv" ]; then
        echo
        echo "Latest comparison:"
        cat "$REPORT_DIR/comparisons/combined.csv"
    fi
    echo
    
    echo "Cache:"
    echo "  Size: $(du -sh "$CACHE_DIR" | cut -f1)"
    echo "  Files: $(find "$CACHE_DIR" -type f | wc -l) files"
    echo
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing benchmarks..."

# Create directory structure
create_dirs

# Run benchmarks if requested
if [ "$RUN" = true ]; then
    run_benchmarks
fi

# Monitor benchmarks if requested
if [ "$MONITOR" = true ]; then
    monitor_benchmarks
fi

# Run profiling if requested
if [ "$PROFILE" = true ]; then
    run_profiling
fi

# Compare benchmarks if requested
if [ "$COMPARE" = true ]; then
    compare_benchmarks
fi

# Generate reports if requested
if [ "$REPORT" = true ]; then
    generate_reports
fi

# Clean benchmarks if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean benchmarks? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_benchmarks
    fi
fi

# Backup benchmarks if requested
if [ "$BACKUP" = true ]; then
    backup_benchmarks
fi

# Restore benchmarks if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore benchmarks? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_benchmarks
    fi
fi

# Show statistics
show_stats

echo
echo "Benchmark management completed successfully!"
