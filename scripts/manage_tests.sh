#!/bin/bash
# Script to manage project tests and test data

# Exit on error
set -e

# Default values
TEST_DIR="tests"
TEST_DATA_DIR="tests/data"
COVERAGE_DIR="coverage"
LOG_DIR="logs"
TEST_TYPES="unit,integration,functional,performance"
TEST_PATTERN="test_*.py"
RUN=false
COVERAGE=false
BENCHMARK=false
CLEAN=false
BACKUP=false
RESTORE=false
PARALLEL=true
VERBOSE=false
FORCE=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --test-dir)
            TEST_DIR="$2"
            shift 2
            ;;
        --test-data-dir)
            TEST_DATA_DIR="$2"
            shift 2
            ;;
        --coverage-dir)
            COVERAGE_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --test-types)
            TEST_TYPES="$2"
            shift 2
            ;;
        --pattern)
            TEST_PATTERN="$2"
            shift 2
            ;;
        --run)
            RUN=true
            shift
            ;;
        --coverage)
            COVERAGE=true
            shift
            ;;
        --benchmark)
            BENCHMARK=true
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
        --no-parallel)
            PARALLEL=false
            shift
            ;;
        --verbose)
            VERBOSE=true
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
    
    # Test directories
    for type in ${TEST_TYPES//,/ }; do
        mkdir -p "$TEST_DIR/$type"
    done
    
    # Test data directory
    mkdir -p "$TEST_DATA_DIR"
    
    # Coverage directory
    mkdir -p "$COVERAGE_DIR"
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to run tests
run_tests() {
    echo "Running tests..."
    
    # Build base command
    CMD="pytest"
    if [ "$VERBOSE" = true ]; then
        CMD="$CMD -v"
    fi
    if [ "$PARALLEL" = true ]; then
        CMD="$CMD -n auto"
    fi
    if [ "$COVERAGE" = true ]; then
        CMD="$CMD --cov=binding_data_processor --cov-report=html:$COVERAGE_DIR"
    fi
    
    # Run tests for each type
    for type in ${TEST_TYPES//,/ }; do
        echo "Running $type tests..."
        
        # Build test command
        TEST_CMD="$CMD $TEST_DIR/$type/$TEST_PATTERN"
        if [ "$BENCHMARK" = true ] && [ "$type" = "performance" ]; then
            TEST_CMD="$TEST_CMD --benchmark-only --benchmark-json=$COVERAGE_DIR/benchmarks.json"
        fi
        
        # Run command
        echo "Running: $TEST_CMD"
        $TEST_CMD || true  # Continue even if tests fail
    done
    
    # Generate combined coverage report
    if [ "$COVERAGE" = true ]; then
        coverage combine
        coverage html -d "$COVERAGE_DIR"
        coverage report
    fi
    
    # Generate benchmark report
    if [ "$BENCHMARK" = true ]; then
        python -m pytest_benchmark.cli compare "$COVERAGE_DIR/benchmarks.json"
    fi
}

# Function to clean test files
clean_tests() {
    echo "Cleaning test files..."
    
    # Clean coverage files
    rm -rf "$COVERAGE_DIR"/*
    rm -f .coverage*
    
    # Clean pytest cache
    find "$TEST_DIR" -type d -name "__pycache__" -exec rm -rf {} +
    find "$TEST_DIR" -type d -name ".pytest_cache" -exec rm -rf {} +
    
    # Clean test data
    rm -rf "$TEST_DATA_DIR"/*
}

# Function to backup test data
backup_tests() {
    echo "Backing up test data..."
    
    # Create backup directory
    BACKUP_DIR="backups/tests_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$TEST_DATA_DIR" "$BACKUP_DIR/"
    cp -r "$COVERAGE_DIR" "$BACKUP_DIR/"
    
    # Create archive
    tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
    rm -rf "$BACKUP_DIR"
    
    echo "Backup saved to: $BACKUP_DIR.tar.gz"
}

# Function to restore test data
restore_tests() {
    echo "Restoring test data..."
    
    # Find latest backup
    LATEST_BACKUP=$(ls -t backups/tests_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$TEST_DATA_DIR" "$COVERAGE_DIR"
    fi
    
    cp -r "$BACKUP_DIR/data" "$TEST_DATA_DIR"
    cp -r "$BACKUP_DIR/coverage" "$COVERAGE_DIR"
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Test data restored from: $BACKUP_DIR"
}

# Function to show test statistics
show_stats() {
    echo "Test statistics:"
    echo
    
    echo "Test files:"
    for type in ${TEST_TYPES//,/ }; do
        echo "  $type tests: $(find "$TEST_DIR/$type" -type f -name "$TEST_PATTERN" | wc -l) files"
        echo "    Test cases: $(grep -r "def test_" "$TEST_DIR/$type" | wc -l)"
        echo "    Assertions: $(grep -r "assert" "$TEST_DIR/$type" | wc -l)"
    done
    echo
    
    echo "Test data:"
    echo "  Files: $(find "$TEST_DATA_DIR" -type f | wc -l) files"
    echo "  Size: $(du -sh "$TEST_DATA_DIR" | cut -f1)"
    echo
    
    if [ -d "$COVERAGE_DIR" ]; then
        echo "Coverage:"
        if [ -f "$COVERAGE_DIR/index.html" ]; then
            echo "  Report: $COVERAGE_DIR/index.html"
        fi
        if [ -f "$COVERAGE_DIR/benchmarks.json" ]; then
            echo "  Benchmarks: $COVERAGE_DIR/benchmarks.json"
        fi
    fi
}

# Main process
echo "Managing tests..."

# Create directory structure
create_dirs

# Run tests if requested
if [ "$RUN" = true ]; then
    run_tests
fi

# Clean test files if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean test files? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_tests
    fi
fi

# Backup test data if requested
if [ "$BACKUP" = true ]; then
    backup_tests
fi

# Restore test data if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore test data? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_tests
    fi
fi

# Show statistics
show_stats

echo
echo "Test management completed successfully!"
