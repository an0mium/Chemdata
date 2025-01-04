#!/bin/bash
# Script to demonstrate BBB permeability prediction with validation compounds

# Set up environment
export PYTHONPATH=.

# Create directories
mkdir -p models/bbb
mkdir -p cache
mkdir -p output

# Run prediction script
echo "Running BBB permeability predictions..."
python scripts/predict_bbb.py \
    examples/data/example_compounds.tsv \
    --output-file output/bbb_predictions.tsv \
    --model-dir models/bbb \
    --cache-dir cache \
    --log-level DEBUG

# Analyze results
echo -e "\nResults saved to output/bbb_predictions.tsv"
echo -e "\nAnalysis Summary:"
echo "==================="
echo "1. CNS-Active Compounds (Expected: High BBB Permeability)"
echo "   - Caffeine, Amphetamine, Ketamine, etc."
echo "   These compounds should show high BBB permeability scores"
echo -e "\n2. Peripherally Selective Compounds (Expected: Low BBB Permeability)"
echo "   - Loperamide, Diphenoxylate, Domperidone, etc."
echo "   These compounds should show low BBB permeability scores"
echo -e "\nValidation:"
echo "==========="
echo "The predictions should correctly classify:"
echo "- CNS-active drugs as BBB permeable"
echo "- Peripherally selective drugs as BBB impermeable"
echo "This validates the predictor's ability to distinguish between:"
echo "- Compounds that can cross the BBB"
echo "- Compounds with similar receptor binding but low BBB permeability"

# Display first few predictions
echo -e "\nSample Predictions:"
echo "=================="
head -n 5 output/bbb_predictions.tsv

# Make results directory readable
chmod -R 755 output/

echo -e "\nDone! Full results available in output/bbb_predictions.tsv"
