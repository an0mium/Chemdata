# BBB Permeability Prediction Package

This package provides a comprehensive set of tools for predicting blood-brain barrier (BBB) permeability of chemical compounds. It combines machine learning models, transporter analysis, and web data enrichment to provide accurate predictions with supporting evidence.

## Package Structure

```
bbb/
├── __init__.py          # Package exports
├── base.py             # Core BBB prediction functionality
├── predictors.py       # BBB prediction with transporter analysis
├── integration.py      # BBB prediction with integrated ML models
├── enrichment.py       # BBB prediction with web data enrichment
└── tests/             # Test suite
    └── test_bbb.py    # Unit tests
```

## Features

1. Core BBB Prediction
   - Molecular fingerprint analysis
   - Descriptor-based prediction
   - Enhanced feature extraction
   - Confidence scoring

2. Transporter Analysis
   - P-glycoprotein (P-gp) substrate prediction
   - BCRP substrate prediction
   - Other BBB-related transporters
   - Receptor-mediated transport

3. ML Model Integration
   - Abuse potential prediction
   - Toxicity prediction
   - Receptor binding prediction
   - Psychoactive effects prediction
   - Nootropic activity prediction

4. Web Data Enrichment
   - ChEMBL data integration
   - PubChem data integration
   - Swiss* services integration
   - Community data sources
   - Social media monitoring
   - LLM-based data extraction

## Validation Dataset

The package includes a validation dataset (`examples/data/example_compounds.tsv`) with two groups of compounds:

1. CNS-Active (BBB Permeable)
   - Caffeine
   - Amphetamine
   - Ketamine
   - Morphine
   - LSD
   - Fluoxetine
   - Diazepam
   - Nicotine

2. Peripherally Selective (Low BBB Permeability)
   - Loperamide
   - Diphenoxylate
   - Domperidone
   - Fexofenadine
   - N-Methylnaltrexone
   - Cetirizine
   - Ranitidine
   - Ondansetron
   - Butorphanol
   - Naloxegol

This dataset is designed to validate the predictor's ability to distinguish between:
- Compounds that readily cross the BBB (CNS-active)
- Compounds with similar receptor binding but low BBB permeability (peripherally selective)

## Quick Start

### Basic Usage

```python
from binding_data_processor.processors.psychopharm.predictors.bbb import (
    BBBPredictorWebEnriched
)

# Initialize predictor
predictor = BBBPredictorWebEnriched(
    model_dir="models/bbb",
    cache_dir="cache",
)

# Make predictions
result = predictor.predict(compound)
print(f"BBB Class: {result.value}")
print(f"Confidence: {result.confidence:.2f}")

# Export predictions
predictor.export_predictions(
    "predictions.tsv",
    include_supporting_data=True,
    include_web_data=True,
)
```


### Advanced Usage

```python
# Custom transporter configuration
transporters = {
    "efflux": {"p_glycoprotein", "bcrp"},
    "uptake": {"lat1", "mct1"},
}

# Custom web clients
web_clients = {
    "chembl": ChemblClient(cache_dir="cache"),
    "pubchem": PubchemClient(cache_dir="cache"),
}

# Initialize with custom configuration
predictor = BBBPredictorWebEnriched(
    model_dir="models/bbb",
    cache_dir="cache",
    transporters=transporters,
    web_clients=web_clients,
)

# Make predictions with detailed output
result = predictor.predict(compound)
print(f"BBB Class: {result.value}")
print(f"Confidence: {result.confidence:.2f}")
print("\nSupporting Data:")
for key, value in result.supporting_data.items():
    print(f"  {key}: {value}")

# Export with custom columns
predictor.export_predictions(
    "predictions.tsv",
    columns=[
        "compound_name",
        "permeability_class",
        "confidence",
        "transporter_data",
        "receptor_data",
    ],
    include_supporting_data=True,
    include_web_data=True,
)
```


### Running Example Script

The package includes a script to run predictions on the validation dataset:

```bash
# Make script executable
chmod +x examples/scripts/run_bbb_prediction.sh

# Run predictions
./examples/scripts/run_bbb_prediction.sh
```

This will:
1. Create necessary directories (models, cache, output)
2. Run predictions on example compounds
3. Export results to output/bbb_predictions.tsv
4. Display a summary of the results

## Development

### Setting Up Development Environment

1. Create virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # Linux/macOS
# or
venv\Scripts\activate  # Windows
```

2. Install development dependencies:
```bash
pip install -r requirements-dev.txt
```

3. Install pre-commit hooks:
```bash
pre-commit install
```

### Running Tests

```bash
# Run all tests
pytest binding_data_processor/processors/psychopharm/predictors/bbb/tests/

# Run specific test class
pytest binding_data_processor/processors/psychopharm/predictors/bbb/tests/test_bbb.py::TestBBBPredictorWebEnriched

# Run with coverage
pytest --cov=binding_data_processor.processors.psychopharm.predictors.bbb tests/
```

### Code Style

- Follow PEP 8 guidelines
- Use type hints
- Write docstrings in Google format
- Keep functions focused and under 50 lines
- Add tests for new functionality

### Adding New Features

1. Create feature branch:
```bash
git checkout -b feature/your-feature-name
```

2. Implement feature:
   - Add tests first (TDD)
   - Implement feature
   - Update documentation
   - Run tests and linting

3. Submit pull request:
   - Clear description of changes
   - Link to related issues
   - Test coverage report
   - Documentation updates

## Contributing

See [CONTRIBUTING.md](../../../../../CONTRIBUTING.md) for detailed guidelines on:
- Code style
- Testing requirements
- Documentation standards
- Pull request process

## License

This project is licensed under the MIT License. See [LICENSE](../../../../../LICENSE) for details.
