# Chemdata

Chemical data processing and analysis tools with a focus on receptor ligands.

## Features

- Comprehensive data collection from multiple sources:
  - BindingDB
  - PubChem
  - ChEMBL
  - PubMed
  - Wikipedia
  - Patents
  - Community sources (PsychonautWiki, Erowid)
- Support for 5-HT2 and NMDA receptor ligands
- Progress indicators for long-running operations
- Detailed binding data analysis
- Chemical structure validation
- Web scraping and data enrichment

## Setup

1. Clone the repository
2. Install dependencies:
```bash
pip install pandas rdkit requests beautifulsoup4 pubchempy tqdm
```

## Data Files

Large data files are not included in the repository. They should be placed in the `data/` directory:

- BindingDB_All.tsv (Download from [BindingDB](https://bindingdb.org/bind/chemsearch/marvin/Download.jsp))
- Other CSV files will be generated during processing

## Usage

```python
from binding_data_processor import BindingDataProcessor
from api_client import PubMedClient

# Initialize processor
processor = BindingDataProcessor(PubMedClient())

# Generate input CSV with 5-HT2 ligands
processor.save_5ht2_ligands("data/5ht2_ligands.tsv")

# Process ligands with full data enrichment
processor.process_file("data/5ht2_ligands.tsv", "data/enriched_5ht2_ligands.tsv")
```

## Progress Indicators

The script provides detailed progress information for:
- Reading BindingDB data
- Processing compounds
- Fetching data from external sources
- Enriching compound information
- Saving results

## Output Format

The enriched TSV file includes:
- Basic compound information (name, SMILES, InChI)
- Common names (ranked by search results)
- Binding data (up to 12 targets)
- Activity types
- Legal status
- Pharmacology
- References and URLs

