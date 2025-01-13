# Compound Export Guide

This guide explains how to use the compound export functionality to generate comprehensive TSV and Excel files containing compound data from multiple sources.

## Overview

The export functionality:
1. Collects compounds from BindingDB matching specific criteria (5-HT2 agonists, NMDA antagonists, etc.)
2. Adds known compounds not in BindingDB
3. Includes protein/peptide data from UniProt
4. Includes basic biomolecule data
5. Enriches all data with additional information from:
   - PubChem
   - ChEMBL
   - Patents
   - Scientific literature
   - Community sources
   - Social media mentions

## Prerequisites

1. Python 3.8 or higher
2. Virtual environment (recommended)
3. Required Python packages:
   - pandas
   - openpyxl
   - rdkit
   - All packages listed in requirements.txt

## Setup

1. Clone the repository:
```bash
git clone https://github.com/yourusername/chemdata.git
cd chemdata
```

2. Run the setup script:
```bash
./scripts/export_all_compounds.sh
```

This will:
- Create necessary directories
- Set up a virtual environment
- Install dependencies
- Create a template for custom compounds

## Adding Custom Compounds

You can add custom compounds in two ways:

### 1. Edit the Custom Compounds Template

The script creates a template file at `data/custom_compounds.tsv` with these columns:
- Name: Common name of the compound
- CAS_Number: CAS registry number
- SMILES: SMILES notation of the structure
- Compound_Type: One of [PSYCHOACTIVE, NOOTROPIC, BIOMOLECULE, OTHER]

Example:
```tsv
Name    CAS_Number      SMILES                                  Compound_Type
DOI     83619-26-1     CC(NC)C(C1=CC(=C(C=C1)OC)I)C           PSYCHOACTIVE
```

### 2. Modify the Source Code

For programmatically defined compounds, edit the `CUSTOM_COMPOUNDS` dictionary in `scripts/export_compounds.py`:

```python
CUSTOM_COMPOUNDS = {
    "5-HT2_AGONISTS": [
        {
            "name": "DOI",
            "cas": "83619-26-1",
            "smiles": "CC(NC)C(C1=CC(=C(C=C1)OC)I)C"
        },
        # Add more compounds...
    ]
}
```

## Running the Export

1. Simply run:
```bash
./scripts/export_all_compounds.sh
```

2. The script will:
   - Download data from BindingDB
   - Process custom compounds
   - Enrich data from multiple sources
   - Export results in TSV and Excel formats

## Output Files

The script generates two files in the `output` directory:

1. `compounds_YYYYMMDD_HHMMSS.tsv`: Tab-separated values file containing:
   - Basic compound information (name, CAS, SMILES, etc.)
   - Physical properties (MW, LogP, TPSA, etc.)
   - Target information
   - Binding data
   - References (DOIs, PubMed IDs, patents)
   - Predictions (activity, toxicity, abuse potential)

2. `compounds_YYYYMMDD_HHMMSS.xlsx`: Excel version of the same data, useful for:
   - Easier viewing and filtering
   - Data analysis
   - Sharing with collaborators

## Understanding the Output

The output files contain these columns:

### Basic Information
- Name: Common name of the compound
- CAS_Number: CAS registry number
- SMILES: Standardized SMILES notation
- InChI: International Chemical Identifier
- InChI_Key: Hashed InChI for easy lookup
- Compound_Type: Classification of the compound

### Physical Properties
- Molecular_Weight: Exact molecular weight
- LogP: Calculated partition coefficient
- TPSA: Topological polar surface area
- HBD: Number of hydrogen bond donors
- HBA: Number of hydrogen bond acceptors
- Rotatable_Bonds: Number of rotatable bonds

### Target Information
- Target_Name: Name of the primary target
- Target_Organism: Source organism
- Affinity_Type: Type of binding measurement (Ki, IC50, etc.)
- Affinity_Value: Numerical value in nM
- Confidence: Confidence score for the data

### References
- DOIs: Digital Object Identifiers for publications
- PubMed_IDs: PubMed reference IDs
- Patents: Related patent numbers

### Predictions
- Predicted_Activity: Predicted biological activity
- Predicted_Toxicity: Predicted toxicity risks
- Predicted_Abuse: Predicted abuse potential

## Troubleshooting

1. If the script fails to download data:
   - Check your internet connection
   - Verify the BindingDB URLs are accessible
   - Check the logs directory for error messages

2. If compound processing fails:
   - Verify SMILES strings are valid
   - Check CAS numbers are correctly formatted
   - Look for error messages in the logs

3. For enrichment failures:
   - Check API access to external services
   - Verify rate limits haven't been exceeded
   - Review the logs for specific error messages

## Best Practices

1. Always validate custom compound data:
   - Use correct SMILES notation
   - Verify CAS numbers
   - Choose appropriate compound types

2. Keep the custom compounds file organized:
   - Group similar compounds together
   - Add comments for clarity
   - Document data sources

3. Regular maintenance:
   - Update custom compound data
   - Check for new versions of BindingDB data
   - Verify external service access

## Advanced Usage

### Filtering Output

You can filter the output using pandas:

```python
import pandas as pd

# Read the TSV file
df = pd.read_csv("output/compounds_20230615_120000.tsv", sep="\t")

# Filter for specific compound types
psychoactive = df[df["Compound_Type"] == "PSYCHOACTIVE"]

# Filter by affinity
strong_binders = df[df["Affinity_Value"] < 100]  # Less than 100 nM

# Export filtered results
strong_binders.to_csv("strong_binders.tsv", sep="\t", index=False)
```

### Adding New Data Sources

To add a new data source:

1. Create a new client in `binding_data_processor/web_enrichment/clients/`
2. Add the client to the `clients` dictionary in `export_compounds.py`
3. Update the enrichment logic in `enrich_compound_data()`

### Customizing Output Format

To modify the output format:

1. Update the row dictionary in `export_results()`
2. Add new columns as needed
3. Modify the export format in `df.to_csv()` or `df.to_excel()`
