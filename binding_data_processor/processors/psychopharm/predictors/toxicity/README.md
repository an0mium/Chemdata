# Toxicity Prediction Package

This package provides a comprehensive set of tools for predicting toxicological properties of chemical compounds. It combines multiple prediction models to assess various aspects of toxicity, including mechanism of action, organ-specific effects, and safety concerns.

## Package Structure

```
toxicity/
├── __init__.py          # Package exports
├── base.py             # Core toxicity prediction
├── integration.py      # Enhanced prediction with ensemble models
├── enrichment.py       # Web data enrichment
└── tests/             # Test suite
    ├── test_toxicity.py       # Base predictor tests
    └── test_pipeline.py       # Integration tests
```

## Features

1. Core Toxicity Prediction
   - Toxicity class prediction
   - Mechanism of action prediction
   - Organ-specific toxicity assessment
   - Safety concern evaluation
   - Confidence scoring

2. Toxicity Mechanisms
   - Cellular mechanisms:
     - Oxidative stress
     - Mitochondrial toxicity
     - DNA damage
     - Protein adducts
     - Lipid peroxidation
   - Molecular mechanisms:
     - Reactive metabolites
     - Free radical formation
     - Protein crosslinking
     - Ion channel effects
   - Systemic mechanisms:
     - Immune activation
     - Inflammation
     - Organ failure
     - Metabolic disruption

3. Organ-Specific Toxicity
   - Liver toxicity:
     - Hepatocellular damage
     - Cholestasis
     - Steatosis
     - Enzyme elevation
   - Kidney toxicity:
     - Tubular damage
     - Glomerular damage
     - Filtration impairment
   - Cardiac toxicity:
     - Arrhythmia
     - QT prolongation
     - Structural changes
   - Neurotoxicity:
     - Cognitive impairment
     - Seizure risk
     - Behavioral changes

4. Safety Assessment
   - Acute toxicity:
     - LD50 prediction
     - Immediate effects
     - Overdose risk
   - Chronic toxicity:
     - Carcinogenicity
     - Mutagenicity
     - Reproductive toxicity
   - Special concerns:
     - Drug interactions
     - Contraindications
     - Vulnerable populations

5. Ensemble Models
   - Class prediction ensemble
   - Mechanism prediction ensembles
   - Organ toxicity ensembles
   - Safety prediction ensembles
   - Model versioning and persistence

6. Web Data Enrichment
   - Literature analysis
   - Clinical trial data
   - Safety databases
   - Case reports
   - Regulatory documents

## Validation Dataset

The package includes a validation dataset with compounds of known toxicity:

1. High Toxicity
   - Strychnine
   - Tetrodotoxin
   - Ricin
   - Botulinum toxin
   - Amatoxins

2. Moderate Toxicity
   - Acetaminophen
   - Ethanol
   - Caffeine
   - Nicotine
   - Aspirin

3. Low Toxicity
   - Vitamin C
   - Glycine
   - Glucose
   - Melatonin
   - L-Theanine

## Quick Start

### Basic Usage

```python
from binding_data_processor.processors.psychopharm.predictors.toxicity import (
    ToxicityPredictor
)

# Initialize predictor
predictor = ToxicityPredictor(
    model_dir="models/toxicity",
    cache_dir="cache",
)

# Make predictions
result = predictor.predict(compound)
print(f"Toxicity class: {result.value}")
print(f"Confidence: {result.confidence:.2f}")

# Get detailed predictions
mechanisms = predictor.predict_mechanisms(compound)
organ_effects = predictor.predict_organ_toxicity(compound)
safety_concerns = predictor.predict_safety_concerns(compound)

# Get prediction statistics
stats = predictor.get_prediction_statistics()
```

### Advanced Usage

```python
# Custom toxicity mechanisms
mechanisms = {
    "cellular": {"oxidative_stress", "dna_damage", "membrane_disruption"},
    "molecular": {"protein_binding", "enzyme_inhibition"},
    "systemic": {"inflammation", "organ_failure"},
}

# Custom organ toxicity types
organ_toxicity = {
    "liver": {"hepatotoxicity", "enzyme_elevation"},
    "kidney": {"nephrotoxicity", "filtration_impairment"},
    "heart": {"cardiotoxicity", "arrhythmia"},
}

# Custom safety concerns
safety_concerns = {
    "acute": {"immediate_effects", "overdose_risk"},
    "chronic": {"long_term_effects", "cumulative_toxicity"},
    "special": {"drug_interactions", "contraindications"},
}

# Initialize with custom configuration
predictor = ToxicityPredictor(
    model_dir="models/toxicity",
    cache_dir="cache",
    toxicity_mechanisms=mechanisms,
    organ_toxicity=organ_toxicity,
    safety_concerns=safety_concerns,
)

# Make predictions
result = predictor.predict(compound)
print(f"Toxicity class: {result.value}")
print(f"Confidence: {result.confidence:.2f}")

# Analyze mechanisms
for category, mechanisms in result.supporting_data['mechanisms'].items():
    print(f"\n{category.title()} Mechanisms:")
    for mechanism, data in mechanisms.items():
        print(f"  {mechanism}: score={data['score']:.2f}, confidence={data['confidence']:.2f}")

# Analyze organ toxicity
for organ, toxicities in result.supporting_data['organs'].items():
    print(f"\n{organ.title()} Toxicity:")
    for toxicity, data in toxicities.items():
        print(f"  {toxicity}: severity={data['severity']:.2f}, confidence={data['confidence']:.2f}")

# Analyze safety concerns
for category, concerns in result.supporting_data['concerns'].items():
    print(f"\n{category.title()} Concerns:")
    for concern, data in concerns.items():
        print(f"  {concern}: risk={data['risk']:.2f}, confidence={data['confidence']:.2f}")
```

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
pytest binding_data_processor/processors/psychopharm/predictors/toxicity/tests/

# Run specific test class
pytest binding_data_processor/processors/psychopharm/predictors/toxicity/tests/test_toxicity.py::TestToxicityPredictor

# Run with coverage
pytest --cov=binding_data_processor.processors.psychopharm.predictors.toxicity tests/
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
