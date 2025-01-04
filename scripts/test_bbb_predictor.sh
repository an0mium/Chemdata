#!/bin/bash
# Script to run BBB predictor tests and generate coverage report

# Set up environment
export PYTHONPATH=.

# Create directories
mkdir -p reports/coverage
mkdir -p reports/test_results

# Run tests with coverage
echo "Running BBB predictor tests with coverage..."
pytest binding_data_processor/processors/psychopharm/predictors/bbb/tests/ \
    --cov=binding_data_processor.processors.psychopharm.predictors.bbb \
    --cov-report=html:reports/coverage \
    --cov-report=term \
    --verbose \
    --junit-xml=reports/test_results/bbb_predictor_tests.xml

# Run integration test with example compounds
echo -e "\nRunning integration test with example compounds..."
./examples/scripts/run_bbb_prediction.sh

# Display coverage report location
echo -e "\nTest Results:"
echo "=============="
echo "Coverage report: reports/coverage/index.html"
echo "Test results: reports/test_results/bbb_predictor_tests.xml"

# Check if any tests failed
if [ $? -eq 0 ]; then
    echo -e "\nAll tests passed successfully!"
else
    echo -e "\nSome tests failed. Please check the test results for details."
    exit 1
fi
