# Chemical Data Collector

A robust tool for collecting and aggregating chemical compound data from multiple sources, with a focus on psychoactive compounds, NPS (novel psychoactive substances), research chemicals, and neurotransmitters.

## Features

- **Multi-source Data Collection**: Aggregates data from PubChem, BindingDB, and other authoritative sources
- **Robust API Handling**: 
  - Circuit breaker pattern for API resilience
  - Rate limiting and automatic retries
  - Response caching to minimize API calls
- **Efficient Processing**:
  - Batch processing for large datasets
  - Parallel execution using ThreadPoolExecutor
  - Progress tracking with tqdm
- **Data Validation**:
  - CAS number format and checksum validation
  - Chemical structure validation using RDKit
  - Comprehensive data field validation
- **Flexible Output**:
  - CSV export with configurable fields
  - Detailed processing statistics
  - Cache usage analytics

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/chemical-data-collector.git
cd chemical-data-collector
```

2. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Command Line Interface

Process a CSV file containing compound data:

```bash
python cli.py -i input.csv -o output.csv
```

Additional options:
```bash
python cli.py --help
```

```
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input CSV file containing compound data
  -o OUTPUT, --output OUTPUT
                        Output CSV file for processed data
  -s SOURCES [SOURCES ...], --sources SOURCES [SOURCES ...]
                        Data sources to use (default: all enabled sources)
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Logging level
  --clear-cache         Clear cache before processing
```

### Input CSV Format

The input CSV should contain at least one of these identifier columns:
- `cas`: CAS Registry Number
- `name`: Common or generic name
- `smiles`: SMILES notation
- `inchi_key`: InChIKey

Example input.csv:
```csv
cas,name,smiles
50-78-2,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
58-08-2,Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
```

### Python API

```python
from data_processor import DataProcessor

# Initialize processor
processor = DataProcessor()

# Process file
stats = processor.process_file(
    input_file="input.csv",
    output_file="output.csv",
    sources=['pubchem', 'bindingdb']
)

# Print statistics
print(f"Processed {stats['total_compounds']} compounds")
print(f"Success rate: {stats['success_rate']:.1f}%")
```

## Architecture

The project is organized into several modules:

- `models.py`: Core data structures and validation logic
- `api_client.py`: API interaction with circuit breaker pattern
- `cache_manager.py`: Caching system for API responses
- `data_processor.py`: Batch processing and parallel execution
- `logger.py`: Logging configuration and utilities
- `cli.py`: Command-line interface
- `config.py`: Configuration settings

## Improvements

1. **Code Organization**:
   - Modular architecture with clear separation of concerns
   - Type hints and comprehensive docstrings
   - Consistent code style with flake8 and black

2. **Performance**:
   - Caching system to reduce API calls
   - Parallel processing for better throughput
   - Batch operations for efficient data handling

3. **Reliability**:
   - Circuit breaker pattern for API resilience
   - Comprehensive error handling
   - Data validation at multiple levels

4. **Monitoring**:
   - Detailed logging system
   - Processing statistics
   - Cache analytics

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Run tests: `pytest`
5. Submit a pull request

## License

MIT License - see LICENSE file for details
