#!/bin/bash
set -e

echo "Starting automated compound harvesting and export process..."

# Create required directories
mkdir -p data cache output logs/{error,validation,warnings,sources}

# Setup Python environment
if [ -d "venv-3.11" ]; then
    source venv-3.11/bin/activate
elif [ -d "venv-3.12" ]; then
    source venv-3.12/bin/activate
elif [ -d "venv" ]; then
    source venv/bin/activate
else
    echo "Creating virtual environment..."
    python3 -m venv venv
    source venv/bin/activate
fi

# Check for conda and use if available
if command -v conda &> /dev/null; then
    echo "Conda detected. You can optionally use conda for RDKit installation:"
    echo "conda install -c conda-forge rdkit"
fi

# Install/update dependencies
echo "Installing/updating dependencies..."
pip install --upgrade pip setuptools wheel
pip install -r requirements.txt

# Check for uv and use it if available
if command -v uv &> /dev/null; then
    echo "Using uv for faster package installation..."
    uv pip install pandas openpyxl rdkit requests beautifulsoup4 lxml aiohttp
else
    echo "Installing additional dependencies..."
    pip install pandas openpyxl requests beautifulsoup4 lxml aiohttp
fi

# Check for RDKit with detailed error handling
echo "Checking RDKit installation..."
python3 -c "import rdkit" &> /dev/null || {
    echo "Error: RDKit not found. Please install RDKit first."
    echo -e "\nInstallation options:"
    echo "1. Build from source (recommended):"
    echo "   git clone https://github.com/rdkit/rdkit.git"
    echo "   cd rdkit && mkdir build && cd build"
    echo "   cmake .. && make -j4 && make install"
    echo -e "\n2. Use conda (alternative):"
    echo "   conda install -c conda-forge rdkit"
    echo -e "\n3. Build with specific options:"
    echo "   cmake -DRDK_BUILD_PYTHON_WRAPPERS=ON -DRDK_BUILD_COMPRESSED_SUPPLIERS=ON .."
    echo -e "\nFor detailed installation instructions, visit:"
    echo "https://www.rdkit.org/docs/Install.html"
    exit 1
}

# Check if custom compounds file exists
if [ ! -f "data/custom_compounds.tsv" ]; then
    echo "Creating template for custom compounds..."
    echo -e "Name\tCAS_Number\tSMILES\tCompound_Type\tNotes" > data/custom_compounds.tsv
    echo "Created template at data/custom_compounds.tsv"
    echo "Add your custom compounds to this file before running the export"
    echo -e "\nExample compound categories that can be added:"
    echo "1. Receptor-Based Compounds:"
    echo "   - 5-HT2 agonists (e.g., DOI, 2C-B, DOM)"
    echo "   - NMDA antagonists (e.g., Ketamine, Memantine, MK-801)"
    echo "   - Anti-addictive agents (e.g., Ibogaine, 18-MC, Naltrexone)"
    echo "   - Physical enhancement compounds"
    echo "   - Longevity enhancement compounds"
    echo "   - Documented recreational/nootropic compounds"
    echo -e "\n2. Proteins/Peptides:"
    echo "   - Follistatin-288/315"
    echo "   - alpha-Klotho"
    echo "   - Myoglobin/Hemoglobin"
    echo "   - Profilin"
    echo "   - Apolipoprotein E/A-I Milano"
    echo "   - Ferritin"
    echo "   - Tubulin (alpha/beta/gamma/delta/epsilon)"
    echo "   - Actin/Troponin/Myosin"
    echo -e "\n3. Basic Biomolecules:"
    echo "   - Creatinine/Creatine"
    echo "   - ATP/ADP/AMP"
    echo "   - NAD+/NADH"
    echo "   - NADP+/NADPH"
    echo "   - FAD/CoA/Acetyl-CoA"
    echo -e "\n4. Example Format:"
    echo "   Name           CAS          SMILES                          Type        Notes"
    echo "   Ketamine      6740-88-1    CN1CCCCC1=O                     NMDA_antag  Used in medicine"
    echo "   Ibogaine      83-74-9      CC1=C2CC3C4CC...               Anti_add    Plant alkaloid"
fi

# Clear previous log files
rm -f logs/{error,validation,warnings,sources}/*

# Run the export script with comprehensive logging
echo "Running automated compound harvesting and export..."
{
    python3 scripts/export_compounds.py 2>&1 | tee -a logs/error/export.log
} || {
    echo "Export script failed. Check logs/error/export.log for details."
    exit 1
}

# Check exit status and provide detailed output
if [ $? -eq 0 ]; then
    echo "Export completed successfully"
    echo -e "\nOutput files:"
    
    # List generated files with details
    echo "Generated files:"
    ls -lh output/compounds_*
    
    # Print summary if files exist
    latest_tsv=$(ls -t output/compounds_*.tsv 2>/dev/null | head -n1)
    if [ -n "$latest_tsv" ]; then
        echo -e "\nExport summary:"
        total_compounds=$(( $(wc -l < "$latest_tsv") - 1 ))
        echo "Total compounds: $total_compounds"
        echo "File sizes:"
        ls -lh output/compounds_* | awk '{print $9 ": " $5}'
        
        # Print compound type breakdown
        if command -v awk &> /dev/null; then
            echo -e "\nCompound type breakdown:"
            awk -F'\t' 'NR>1 {count[$4]++} END {
                printf "\n%-30s %s\n", "Type", "Count"
                printf "%-30s %s\n", "--------------------", "-----"
                for (type in count) printf "%-30s %d\n", type, count[type]
            }' "$latest_tsv" | sort -rn -k2
        fi
        
        # Print validation statistics
        echo -e "\nValidation statistics:"
        structures=$(awk -F'\t' 'NR>1 && $3!="" {count++} END {print count}' "$latest_tsv")
        cas_numbers=$(awk -F'\t' 'NR>1 && $2!="" {count++} END {print count}' "$latest_tsv")
        complete=$(awk -F'\t' 'NR>1 && $2!="" && $3!="" && $4!="" {count++} END {print count}' "$latest_tsv")
        
        printf "%-30s %5d (%d%%)\n" "Compounds with structures:" $structures $(( structures * 100 / total_compounds ))
        printf "%-30s %5d (%d%%)\n" "Compounds with CAS numbers:" $cas_numbers $(( cas_numbers * 100 / total_compounds ))
        printf "%-30s %5d (%d%%)\n" "Compounds with complete data:" $complete $(( complete * 100 / total_compounds ))
    fi
    
    echo -e "\nOutput locations:"
    echo "- TSV format: output/compounds_*.tsv"
    echo "- JSON format: output/compounds_*.json"
    echo "- Excel format: output/compounds_*.xlsx"
    
    # Print validation summary if available
    if [ -f "logs/validation/summary.log" ]; then
        echo -e "\nValidation summary:"
        tail -n 5 logs/validation/summary.log
    fi
    
    # Print any warnings if they exist
    if [ -f "logs/warnings/export.log" ]; then
        echo -e "\nWarnings during export:"
        cat logs/warnings/export.log
    fi
    
    echo -e "\nNext steps:"
    echo "1. Review the exported files in the output directory"
    echo "2. Check logs directory for detailed processing information"
    echo "3. Add any missing compounds to data/custom_compounds.tsv"
    echo "4. Run this script again to update exports with new data"
    
    # Print data source summary if available
    if [ -f "logs/sources/summary.log" ]; then
        echo -e "\nData source summary:"
        cat logs/sources/summary.log
    fi
else
    echo "Export failed. Check logs for details."
    if [ -f "logs/error/export.log" ]; then
        echo -e "\nLast few error messages:"
        tail -n 10 logs/error/export.log
    fi
    exit 1
fi

# Deactivate virtual environment
deactivate
