# Web Enrichment Examples

This directory contains example files and scripts demonstrating how to use the web enrichment functionality.

## Files

- `data/example_compounds.json`: Example compounds in JSON format with BBB permeability categories
- `data/example_compounds.tsv`: Same compounds in TSV format with category comments
- `scripts/enrich_compounds.py`: Script to enrich compounds with web data

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

### Example Commands

1. Process JSON file with all data sources:
```bash
python examples/scripts/enrich_compounds.py \
    --reddit-id YOUR_REDDIT_ID \
    --reddit-secret YOUR_REDDIT_SECRET \
    --twitter-token YOUR_TWITTER_TOKEN \
    examples/data/example_compounds.json \
    enriched_compounds.json
```

2. Process TSV file with only Swiss predictions:
```bash
python examples/scripts/enrich_compounds.py \
    --skip-web-data \
    examples/data/example_compounds.tsv \
    enriched_compounds.json
```

3. Process with custom batch size and workers:
```bash
python examples/scripts/enrich_compounds.py \
    --workers 8 \
    --batch-size 20 \
    examples/data/example_compounds.json \
    enriched_compounds.json
```

### Output Format

The script outputs enriched data in JSON format:
```json
{
    "compounds": [
        {
            "category": "CNS-Active (BBB Permeable)",
            "compounds": [
                {
                    "name": "Caffeine",
                    "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                    "cas_number": "58-08-2",
                    "swiss_data": {
                        "targets": [
                            {
                                "target": "Adenosine A2a receptor",
                                "probability": 0.95
                            }
                        ]
                    },
                    "community_data": {
                        "reports": [
                            {
                                "source": "PsychonautWiki",
                                "text": "Example report"
                            }
                        ]
                    },
                    "social_data": {
                        "posts": [
                            {
                                "platform": "Reddit",
                                "text": "Example post"
                            }
                        ]
                    },
                    "enrichment_metadata": {
                        "timestamp": "2024-01-01T12:00:00",
                        "category": "CNS-Active (BBB Permeable)",
                        "sources": [
                            "swiss",
                            "community",
                            "social"
                        ]
                    }
                }
            ]
        }
    ]
}
```

### Caching

The script uses a cache directory at `~/.cache/binding_data_processor` to store:
- HTTP responses
- Model files
- Intermediate results

Use `--no-cache` to disable caching and force fresh data retrieval.

### Error Handling

The script includes robust error handling:
- Validates input files and formats
- Handles API rate limits and retries
- Reports detailed error messages
- Tracks failed compounds
- Provides metrics on success/failure rates

### Metrics

The script outputs detailed metrics including:
- Number of processed compounds
- Number of failed compounds
- Per-client metrics (requests, cache hits, errors)
