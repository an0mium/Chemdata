# Web Enrichment Examples

This directory contains example files and scripts demonstrating how to use the web enrichment functionality.

## Files

- `data/example_compounds.json`: Example compounds in JSON format with BBB permeability categories
- `data/example_compounds.tsv`: Same compounds in TSV format with category comments
- `scripts/enrich_compounds.py`: Script to enrich compounds with web data
- `scripts/search_patents.py`: Script to search and analyze patents across multiple sources

## Using the Enrichment Script

The `enrich_compounds.py` script demonstrates how to:
1. Configure web enrichment
2. Create a manager instance
3. Load compounds from various sources
4. Enrich compounds with web data
5. Export enriched data

### Input Formats

The script supports two input formats:

#### JSON Format
```json
{
    "compounds": [
        {
            "category": "CNS-Active (BBB Permeable)",
            "compounds": [
                {
                    "name": "Caffeine",
                    "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                    "cas_number": "58-08-2"
                }
            ]
        }
    ]
}
```

#### TSV Format
```
name    smiles    cas_number
# CNS-Active (BBB Permeable)
Caffeine    CN1C=NC2=C1C(=O)N(C(=O)N2C)C    58-08-2
```

### Usage

Basic usage:
```bash
python examples/scripts/enrich_compounds.py input.json output.json
```

With options:
```bash
python examples/scripts/enrich_compounds.py \
    --reddit-id YOUR_REDDIT_ID \
    --reddit-secret YOUR_REDDIT_SECRET \
    --twitter-token YOUR_TWITTER_TOKEN \
    --skip-predictions \
    --skip-web-data \
    --no-cache \
    --workers 4 \
    --batch-size 10 \
    --verbose \
    input.json output.json
```

### Options

- `input_path`: Path to input file (JSON or TSV)
- `output_path`: Path to output JSON file
- `--reddit-id`: Reddit client ID for social data
- `--reddit-secret`: Reddit client secret for social data
- `--twitter-token`: Twitter bearer token for social data
- `--skip-predictions`: Skip Swiss predictions
- `--skip-web-data`: Skip web data collection
- `--no-cache`: Disable caching
- `--workers`: Number of worker threads (default: 4)
- `--batch-size`: Batch size for processing (default: 10)
- `--verbose`: Enable debug logging

## Using the Patent Search Script

The `search_patents.py` script demonstrates how to:
1. Search patents across Google Patents, USPTO, and Espacenet
2. Get comprehensive patent data including:
   - Family information from INPADOC
   - Legal status from Espacenet
   - Citations with relevance scoring
   - Extracted chemical compounds
3. Export results in both JSON and TSV formats
4. Generate detailed statistics

### Basic Usage

Search by text query:
```bash
python examples/scripts/search_patents.py \
    --query "5-HT2A antagonist" \
    --output results/
```

Search by chemical structure:
```bash
python examples/scripts/search_patents.py \
    --query "serotonin antagonist" \
    --structure "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC3=CC=C(C=C3)F" \
    --output results/
```

### Advanced Options

- `--query`: Search query (required)
- `--structure`: Chemical structure (SMILES/InChI)
- `--date-from`: Start date (YYYY-MM-DD)
- `--date-to`: End date (YYYY-MM-DD)
- `--classification`: Patent classification code
- `--max-results`: Maximum number of results (default: 100)
- `--output`: Output directory (default: results/)
- `--espacenet-key`: Espacenet API key
- `--uspto-key`: USPTO API key

### Example Commands

1. Search with date range:
```bash
python examples/scripts/search_patents.py \
    --query "NMDA antagonist" \
    --date-from 2020-01-01 \
    --date-to 2024-01-01 \
    --output results/
```

2. Search with classification:
```bash
python examples/scripts/search_patents.py \
    --query "nootropic" \
    --classification "A61K31" \
    --output results/
```

3. Full search with all data sources:
```bash
python examples/scripts/search_patents.py \
    --query "5-HT2A antagonist" \
    --structure "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC3=CC=C(C=C3)F" \
    --date-from 2020-01-01 \
    --classification "A61K31" \
    --espacenet-key YOUR_ESPACENET_KEY \
    --uspto-key YOUR_USPTO_KEY \
    --output results/
```

### Output Format

The script generates two output files:

#### 1. patents.json
Detailed JSON with all patent data:
```json
{
    "title": "Example Patent",
    "abstract": "...",
    "inventors": ["..."],
    "assignee": "...",
    "filing_date": "2024-01-01",
    "publication_date": "2024-01-05",
    "patent_number": "US12345678",
    "compounds": ["..."],
    "effects": ["..."],
    "mechanisms": ["..."],
    "safety_notes": ["..."],
    "family": {
        "family_id": "12345",
        "members": ["US12345678", "EP12345678"],
        "priority_date": "2023-12-01",
        "countries": ["US", "EP"]
    },
    "citations": [
        {
            "patent_number": "US87654321",
            "title": "...",
            "filing_date": "2023-01-01",
            "relevance": 0.8,
            "citation_type": "backward"
        }
    ],
    "classifications": [
        {
            "system": "IPC",
            "code": "A61K31",
            "description": "...",
            "level": "subclass"
        }
    ],
    "legal_status": {
        "status": "granted",
        "date": "2024-01-05",
        "country": "US",
        "description": "Patent granted"
    },
    "metadata": {
        "sources": ["google_patents", "espacenet", "uspto"],
        "scrape_date": "2024-01-10T12:00:00"
    }
}
```

#### 2. patents_summary.tsv
Tabular summary with key fields:
```
patent_number    title    assignee    filing_date    family_size    countries    compounds    legal_status
US12345678    Example Patent    Company Inc    2024-01-01    2    US,EP    Compound A,Compound B    granted
```

### Statistics

The script outputs detailed statistics including:
- Total patents found
- Source distribution
- Family coverage and country distribution
- Legal status distribution
- Compound extraction metrics

### Caching

The script uses a cache directory at `~/.cache/binding_data_processor` to store:
- HTTP responses
- Intermediate results
- Patent family data

Use `--no-cache` to disable caching and force fresh data retrieval.

### Error Handling

The script includes robust error handling:
- Validates input parameters
- Handles API rate limits and retries
- Reports detailed error messages
- Tracks failed requests
- Provides metrics on success/failure rates

### Metrics

The script outputs detailed metrics including:
- Number of processed patents
- Number of failed requests
- Per-source metrics (requests, cache hits, errors)
- Data coverage statistics
