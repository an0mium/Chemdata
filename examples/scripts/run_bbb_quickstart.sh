#!/bin/bash
# Script to run BBB permeability prediction quickstart example

# Set up environment
export PYTHONPATH=.

# Create necessary directories
mkdir -p models/bbb
mkdir -p cache
mkdir -p output

# Print header
echo "BBB Permeability Prediction Quickstart"
echo "====================================="
echo

# Check if example compounds exist
if [ ! -f "examples/data/example_compounds.tsv" ]; then
    echo "Error: Example compounds file not found!"
    echo "Please ensure examples/data/example_compounds.tsv exists."
    exit 1
fi

# Check if models are installed
if [ ! -f "models/bbb/fingerprint_model.pt" ]; then
    echo "Installing BBB predictor models..."
    cd binding_data_processor/processors/psychopharm/predictors/bbb
    ./setup.py --base-dir ../../../../../../
    cd ../../../../../
    echo
fi

# Run quickstart example
echo "Running BBB prediction example..."
echo "--------------------------------"
./examples/scripts/bbb_quickstart.py

# Check exit status
if [ $? -eq 0 ]; then
    echo
    echo "Success! BBB predictions have been generated."
    echo "Check output/bbb_predictions.tsv for results."
    
    # Show preview of results if they exist
    if [ -f "output/bbb_predictions.tsv" ]; then
        echo
        echo "Preview of results:"
        echo "-----------------"
        head -n 5 output/bbb_predictions.tsv
    fi
else
    echo
    echo "Error: BBB prediction failed. Please check the error messages above."
    exit 1
fi
