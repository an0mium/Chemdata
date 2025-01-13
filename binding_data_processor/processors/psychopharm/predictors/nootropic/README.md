# Nootropic Effects Prediction Package

This package provides a comprehensive set of tools for predicting nootropic effects of chemical compounds. It combines ensemble machine learning models, BBB permeability analysis, and web data enrichment to provide accurate predictions with supporting evidence.

## Package Structure

```
nootropic/
├── __init__.py          # Package exports
├── base.py             # Core predictor class
├── features.py         # Feature extraction utilities
├── model_loading.py    # Model loading and management
├── prediction.py       # Prediction utilities
└── tests/             # Test suite
    └── test_nootropic.py  # Comprehensive tests
```

## Features

1. Core Nootropic Prediction
   - Mechanism of action prediction with ensemble models
   - Cognitive domain effects with confidence scoring
   - Side effect profiling with risk assessment
   - Comprehensive feature extraction
   - Model persistence and versioning

2. Feature Extraction
   - Molecular fingerprints
   - Pharmacophore features
   - Binding site features
   - Literature-derived features
   - Community data features
   - Enhanced descriptors

3. BBB Integration
   - BBB permeability prediction
   - Effect scaling by BBB permeability
   - Transport mechanism analysis
   - Confidence adjustment
   - Supporting evidence integration

4. Cognitive Domains
   - Memory enhancement
   - Attention improvement
   - Learning facilitation
   - Executive function
   - Processing speed
   - Mental clarity
   - Neuroplasticity effects

5. Side Effect Analysis
   - Physical side effects
   - Cognitive side effects
   - Tolerance development
   - Drug interactions
   - Safety profiling
   - Risk assessment
   - Confidence scoring

6. Model Ensembles
   - Mechanism prediction ensemble
   - Effect prediction ensembles
   - Side effect prediction ensembles
   - Confidence estimation
   - Model versioning
   - Feature importance analysis

7. Web Data Integration
   - Literature analysis
   - Clinical trial data
   - Community reports
   - Safety databases
   - Patent analysis
   - Evidence weighting

## Validation Dataset

The package includes a validation dataset with known nootropic compounds:

1. Strong Nootropics
   - Piracetam
   - Aniracetam
   - Modafinil
   - Noopept
   - Phenylpiracetam

2. Mild Nootropics
   - Caffeine
   - L-Theanine
   - Bacopa monnieri
   - Ginkgo biloba
   - Lion's mane

3. Control Compounds
   - Aspirin
   - Ibuprofen
   - Acetaminophen
   - Diphenhydramine
   - Cetirizine

## Quick Start

### Basic Usage

```python
from binding_data_processor.processors.psychopharm.predictors.nootropic import NootropicPredictor

# Initialize predictor
predictor = NootropicPredictor(
    model_dir="models/nootropic",
    cache_dir="cache",
)

# Make predictions
result = predictor.predict(compound)
print(f"Mechanism: {result.value}")
print(f"Confidence: {result.confidence:.2f}")

# Get supporting data
print("\nBBB Data:")
print(f"  Permeability: {result.supporting_data['bbb_prediction']['value']}")
print(f"  Confidence: {result.supporting_data['bbb_prediction']['confidence']:.2f}")

# Analyze effects
for domain, effects in result.supporting_data['effects'].items():
    print(f"\n{domain.title()} Effects:")
    for effect, data in effects.items():
        print(f"  {effect}: score={data['score']:.2f}, confidence={data['confidence']:.2f}")

# Analyze side effects
for category, effects in result.supporting_data['side_effects'].items():
    print(f"\n{category.title()} Side Effects:")
    for effect, data in effects.items():
        print(f"  {effect}: risk={data['risk']:.2f}, confidence={data['confidence']:.2f}")
```

### Advanced Usage

```python
# Custom cognitive domains
domains = {
    "memory": {"working_memory", "long_term_memory", "recall_speed"},
    "attention": {"sustained_attention", "divided_attention", "focus"},
    "learning": {"skill_acquisition", "knowledge_retention"},
}

# Custom side effects
side_effects = {
    "physical": {"headache", "insomnia", "fatigue"},
    "cognitive": {"brain_fog", "anxiety", "mood_changes"},
    "tolerance": {"acute_tolerance", "withdrawal", "dependence"},
}

# Initialize with custom configuration
predictor = NootropicPredictor(
    model_dir="models/nootropic",
    cache_dir="cache",
    cognitive_domains=domains,
    side_effects=side_effects,
    device="cuda" if torch.cuda.is_available() else "cpu",
)

# Get feature importance
importances = predictor.get_feature_importance()
for feature_type, scores in importances.items():
    print(f"\n{feature_type} Feature Importance:")
    for feature, score in scores.items():
        print(f"  {feature}: {score:.3f}")

# Retrain models
metrics = predictor.retrain(
    compounds=training_compounds,
    labels=mechanism_labels,
    effects={"memory": {"working_memory": scores}},
    mechanisms={"cholinergic": scores},
)
print("\nTraining Metrics:")
for metric, value in metrics.items():
    print(f"  {metric}: {value:.3f}")
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
pytest binding_data_processor/processors/psychopharm/predictors/nootropic/tests/

# Run specific test class
pytest binding_data_processor/processors/psychopharm/predictors/nootropic/tests/test_nootropic.py::TestNootropicPredictor

# Run with coverage
pytest --cov=binding_data_processor.processors.psychopharm.predictors.nootropic tests/
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
