# ChemData

A comprehensive pipeline for processing and analyzing chemical compound data, with a focus on psychopharmacological compounds. The system combines data from BindingDB with web-enriched information, patent data, and machine learning predictions.

## Features

### Data Processing
- **BindingDB Integration**
  - Automated data processing
  - Structure validation and standardization
  - Property calculation
  - Binding data analysis

- **Additional Data Sources**
  - ChEMBL API integration
  - PubChem data harvesting
  - Swiss* services (SwissTargetPrediction, SwissADME)
  - Patent database search and analysis

- **Community Data**
  - PsychonautWiki API integration
  - Erowid experience reports
  - TripSit factsheets
  - Reddit discussions (r/researchchemicals, r/nootropics)
  - Twitter mentions and trends
  - Bluesky integration

### Machine Learning
- **Binding Predictions**
  - Target-specific models
  - Cross-target interactions
  - Binding site prediction
  - Uncertainty estimation

- **Activity Classification**
  - Mechanism of action
  - Effect classification
  - Duration prediction
  - Structure-activity relationships

- **BBB Permeability Prediction**
  - Core fingerprint-based prediction
  - Transporter analysis (P-gp, BCRP)
  - ML model integration
  - Web data enrichment
  - Comprehensive validation suite

- **Safety Assessment**
  - Toxicity prediction
  - Drug interaction risks
  - Side effect profiles
  - Abuse potential

### Data Enrichment
- **Structure Analysis**
  - 2D/3D conformer generation
  - Pharmacophore detection
  - Similarity search
  - Substructure analysis

- **Property Calculation**
  - Physicochemical properties
  - Drug-likeness scores
  - ADMET predictions
  - Blood-brain barrier penetration

- **Literature Mining**
  - PubMed integration
  - Patent analysis
  - Citation tracking
  - Regulatory status

### Web Interface
- **Compound Browser**
  - Advanced search and filtering
  - Structure visualization
  - Activity data display
  - Prediction visualization

- **Detail Views**
  - Chemical properties
  - Binding profiles
  - Safety information
  - Community data
  - Patent references

- **Export System**
  - Flexible column selection
  - Custom filtering
  - Multiple formats
  - Batch processing

## Requirements

- Python 3.8 or higher
- Docker and Docker Compose
- RDKit
- OpenBabel
- PyTorch (optional, for ML features)
- PostgreSQL
- Redis

## Quick Start

1. Clone the repository:
```bash
git clone https://github.com/yourusername/chemdata.git
cd chemdata
```

2. Create and activate a virtual environment:
```bash
python -m venv .venv
source .venv/bin/activate  # Linux/macOS
# or
.venv\Scripts\activate  # Windows
```

3. Install dependencies:
```bash
./scripts/setup_dev.sh
```

4. Set up environment variables:
```bash
cp .env.example .env
# Edit .env with your configuration
```

5. Start the services:
```bash
docker-compose up -d
```

6. Run the pipeline:
```bash
./scripts/run_pipeline.sh
```

7. Access the web interface:
```
http://localhost:8000
```

## Development

### Setup Development Environment

```bash
# Start development container
docker-compose up -d dev

# Enter development shell
docker-compose exec dev bash

# Install development dependencies
./scripts/setup_dev.sh --dev
```

### Running Tests

```bash
# Run all tests
docker-compose run --rm test

# Run specific tests
docker-compose run --rm test pytest path/to/test.py

# Run tests with coverage
docker-compose run --rm test pytest --cov=binding_data_processor
```

### Code Quality

```bash
# Run linters
pre-commit run --all-files

# Run type checking
mypy binding_data_processor

# Run security checks
bandit -r binding_data_processor
```

### Building Documentation

```bash
# Build documentation
cd docs
make html
```

## Usage Examples

### Command Line Interface

Process compounds from BindingDB:
```bash
python -m binding_data_processor.cli process-compounds \
    --input bindingdb.tsv \
    --output results/ \
    --enable-ml \
    --enable-web \
    --enable-social
```

### Web Application

Run the web interface:
```bash
streamlit run examples/web_app/app.py
```

### Python API

```python
from binding_data_processor.pipeline import ProcessingPipeline
from binding_data_processor.pipeline.config import ProcessingConfig

# Create pipeline
pipeline = ProcessingPipeline(
    config=ProcessingConfig(
        use_ml_predictions=True,
        use_web_enrichment=True,
        use_social_monitoring=True,
    )
)

# Process compounds
compounds = pipeline.process_compounds(
    input_file="bindingdb.tsv",
    output_dir="results/",
)

# Use BBB predictor
from binding_data_processor.processors.psychopharm.predictors.bbb import (
    BBBPredictorWebEnriched
)

predictor = BBBPredictorWebEnriched(
    model_dir="models/bbb",
    cache_dir="cache",
)

result = predictor.predict(compound)
print(f"BBB Class: {result.value}")
print(f"Confidence: {result.confidence:.2f}")
print("\nSupporting Data:")
for key, value in result.supporting_data.items():
    print(f"  {key}: {value}")
```

### Data Processing Scripts

```bash
# Process BindingDB data
./scripts/process_bindingdb.sh \
    --input data/raw/BindingDB_All.tsv \
    --output data/processed/compounds.tsv \
    --workers 4 \
    --batch-size 100

# Enrich compounds
./scripts/enrich_compounds.sh \
    --input data/processed/compounds.tsv \
    --output data/enriched/compounds.tsv \
    --workers 4 \
    --batch-size 100 \
    --rate-limit 2 \
    --sources "chembl,pubchem,swiss,community,social"

# Analyze compounds
./scripts/analyze_compounds.sh \
    --input data/enriched/compounds.tsv \
    --output data/analyzed/compounds.tsv \
    --patent-search \
    --structure-analysis \
    --property-calculation

# Generate report
./scripts/generate_report.sh \
    --input data/analyzed/compounds.tsv \
    --output-dir reports \
    --format html \
    --include-plots
```

## Project Structure

```
binding_data_processor/
├── data_sources/          # Data source integrations
│   ├── bindingdb.py      # BindingDB processing
│   ├── chembl.py         # ChEMBL API client
│   └── pubchem.py        # PubChem integration
├── models/               # Data models and ML
│   ├── compound/        # Compound data models
│   └── psychopharm/     # Psychopharm models
├── pipeline/            # Processing pipeline
│   ├── base.py         # Pipeline coordination
│   ├── ml.py          # ML predictions
│   └── web.py         # Web enrichment
├── processors/         # Data processors
│   ├── structure/     # Structure processing
│   ├── patent/       # Patent analysis
│   └── psychopharm/  # Psychopharm analysis
├── web_enrichment/    # Web data enrichment
│   ├── manager.py    # Enrichment coordination
│   ├── swiss/       # Swiss tools integration
│   └── community/   # Community data sources
└── web/             # Web interface
    ├── api/        # REST API endpoints
    ├── components/ # UI components
    └── pages/      # Web pages
```

## Configuration

The application can be configured through environment variables or a `.env` file:

```bash
# Data directories
CHEMDATA_DATA_DIR=./data
CHEMDATA_CACHE_DIR=./cache
CHEMDATA_LOG_DIR=./logs
CHEMDATA_OUTPUT_DIR=./output
CHEMDATA_MODEL_DIR=./models

# API credentials
REDDIT_CLIENT_ID=your_client_id
REDDIT_CLIENT_SECRET=your_client_secret
TWITTER_API_KEY=your_api_key
TWITTER_API_SECRET=your_api_secret

# Database
POSTGRES_USER=chemdata
POSTGRES_PASSWORD=chemdata
POSTGRES_DB=chemdata
POSTGRES_HOST=postgres

# Redis
REDIS_HOST=redis
REDIS_PORT=6379

# Web server
FLASK_APP=binding_data_processor.web.app
FLASK_ENV=development
FLASK_DEBUG=1
```

## Docker Services

- `web`: Web application and API
- `worker`: Background task worker
- `redis`: Cache and message broker
- `postgres`: Database
- `dev`: Development environment
- `test`: Test runner

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- BindingDB for providing compound data
- ChEMBL for their comprehensive API
- RDKit team for cheminformatics tools
- Open source community for various libraries used
- PsychonautWiki and Erowid for community data
- Swiss Institute of Bioinformatics for web services
- Patent offices for making data publicly accessible

## Citation

If you use this software in your research, please cite:

```bibtex
@software{chemdata2024,
  author = {anomium},
  title = {ChemData: A Comprehensive Pipeline for Psychoactive Compound Analysis},
  year = {2024},
  url = {https://github.com/anomium/chemdata}
}
