#!/bin/bash
# Run tests and generate coverage reports

# Set environment variables
export PYTHONPATH="."
export PYTHONWARNINGS="ignore::DeprecationWarning"

# Create directories if they don't exist
mkdir -p reports/coverage
mkdir -p reports/test-results

# Clean up previous reports
rm -rf reports/coverage/*
rm -rf reports/test-results/*

# Run tests with coverage
echo "Running tests with coverage..."
pytest \
    --verbose \
    --cov=binding_data_processor \
    --cov=web_enrichment \
    --cov=web \
    --cov-report=html:reports/coverage \
    --cov-report=xml:reports/coverage/coverage.xml \
    --cov-report=term \
    --html=reports/test-results/report.html \
    --self-contained-html \
    --durations=10 \
    tests/

# Check test status
status=$?
if [ $status -eq 0 ]; then
    echo "All tests passed!"
else
    echo "Some tests failed!"
    exit $status
fi

# Run type checking
echo "Running type checking..."
mypy \
    binding_data_processor \
    web_enrichment \
    web \
    tests \
    --ignore-missing-imports \
    --html-report reports/mypy

# Run linting
echo "Running linting..."
flake8 \
    binding_data_processor \
    web_enrichment \
    web \
    tests \
    --max-line-length=100 \
    --ignore=E203,W503 \
    --statistics \
    --count

# Run code formatting check
echo "Checking code formatting..."
black \
    binding_data_processor \
    web_enrichment \
    web \
    tests \
    --check \
    --diff

# Run import sorting check
echo "Checking import sorting..."
isort \
    binding_data_processor \
    web_enrichment \
    web \
    tests \
    --check-only \
    --diff

# Print coverage report
echo "Coverage Report:"
coverage report

# Print summary
echo
echo "Test Results: reports/test-results/report.html"
echo "Coverage Report: reports/coverage/index.html"
echo "Type Check Report: reports/mypy/index.html"
