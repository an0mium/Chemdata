# Model Consolidation Plan

## Overview

Current issues:
1. Duplicate model definitions between models/ and models/compound/
2. Psychopharm functionality not integrated
3. Analysis code spread across modules
4. Validation logic duplicated

## Step 1: Directory Structure

```
binding_data_processor/models/compound/
├── base/
│   ├── __init__.py
│   ├── core.py          # Core data model
│   ├── mixins.py        # Shared functionality
│   └── types.py         # Type definitions
├── ml/
│   ├── __init__.py
│   ├── predictors.py    # ML models
│   ├── features.py      # Feature extraction
│   └── training.py      # Model training
├── enrichment/
│   ├── __init__.py
│   ├── web.py          # Web data
│   ├── community.py    # Community data
│   └── social.py       # Social media data
├── analysis/
│   ├── __init__.py
│   ├── binding.py      # Binding analysis
│   ├── activity.py     # Activity analysis
│   ├── safety.py       # Safety analysis
│   └── properties.py   # Property analysis
└── export/
    ├── __init__.py
    ├── formats.py      # Export formats
    └── validation.py   # Export validation
```

## Step 2: Class Hierarchy

```python
# Base Classes
class CompoundBase:
    """Base compound data model."""
    
class ValidationMixin:
    """Validation functionality."""
    
class SerializationMixin:
    """Serialization functionality."""

# ML Classes
class MLCompound(CompoundBase):
    """ML-enabled compound."""
    
class PredictorMixin:
    """Prediction functionality."""

# Enrichment Classes
class EnrichedCompound(MLCompound):
    """Web-enriched compound."""
    
class WebDataMixin:
    """Web data functionality."""

# Analysis Classes
class AnalyzedCompound(EnrichedCompound):
    """Analyzed compound."""
    
class AnalysisMixin:
    """Analysis functionality."""

# Export Classes
class ExportableCompound(AnalyzedCompound):
    """Export-ready compound."""
    
class ExportMixin:
    """Export functionality."""
```

## Step 3: Migration Steps

### Day 1: Setup
1. Create new directory structure
```bash
mkdir -p binding_data_processor/models/compound/{base,ml,enrichment,analysis,export}
touch binding_data_processor/models/compound/{base,ml,enrichment,analysis,export}/__init__.py
```

2. Move existing files
```bash
# Move base files
mv binding_data_processor/models/compound.py binding_data_processor/models/compound/base/core.py
mv binding_data_processor/models/mixins.py binding_data_processor/models/compound/base/mixins.py
mv binding_data_processor/models/types.py binding_data_processor/models/compound/base/types.py

# Move ML files
mv binding_data_processor/models/compound_ml.py binding_data_processor/models/compound/ml/predictors.py

# Move enrichment files
mv binding_data_processor/models/compound_enrichment.py binding_data_processor/models/compound/enrichment/web.py

# Move analysis files
mv binding_data_processor/models/compound_analysis.py binding_data_processor/models/compound/analysis/base.py

# Move export files
mv binding_data_processor/models/compound_export.py binding_data_processor/models/compound/export/formats.py
```

### Day 2: Base Classes
1. Consolidate base functionality
```python
# In base/core.py
class CompoundBase:
    """Base compound data model."""
    def __init__(self):
        self.identifiers = {}
        self.properties = {}
        self.metadata = {}
```

2. Update mixins
```python
# In base/mixins.py
class ValidationMixin:
    """Validation functionality."""
    def validate(self):
        pass

class SerializationMixin:
    """Serialization functionality."""
    def to_dict(self):
        pass
```

### Day 3: ML Integration
1. Consolidate ML functionality
```python
# In ml/predictors.py
class MLCompound(CompoundBase, PredictorMixin):
    """ML-enabled compound."""
    def predict(self, model_name: str):
        pass
```

2. Add feature extraction
```python
# In ml/features.py
class FeatureExtractor:
    """Feature extraction."""
    def extract_features(self, compound: CompoundBase):
        pass
```

### Day 4: Enrichment Integration
1. Consolidate enrichment functionality
```python
# In enrichment/web.py
class EnrichedCompound(MLCompound, WebDataMixin):
    """Web-enriched compound."""
    def enrich(self):
        pass
```

2. Add community data
```python
# In enrichment/community.py
class CommunityDataMixin:
    """Community data functionality."""
    def get_community_data(self):
        pass
```

### Day 5: Analysis Integration
1. Consolidate analysis functionality
```python
# In analysis/base.py
class AnalyzedCompound(EnrichedCompound, AnalysisMixin):
    """Analyzed compound."""
    def analyze(self):
        pass
```

2. Add specific analyses
```python
# In analysis/binding.py
class BindingAnalyzer:
    """Binding analysis."""
    def analyze_binding(self, compound: AnalyzedCompound):
        pass
```

## Step 4: Testing

### Test Structure
```
tests/models/compound/
├── base/
│   ├── test_core.py
│   ├── test_mixins.py
│   └── test_types.py
├── ml/
│   ├── test_predictors.py
│   └── test_features.py
├── enrichment/
│   ├── test_web.py
│   └── test_community.py
├── analysis/
│   ├── test_binding.py
│   └── test_activity.py
└── export/
    ├── test_formats.py
    └── test_validation.py
```

### Test Cases
1. Base functionality
```python
def test_compound_initialization():
    compound = CompoundBase()
    assert compound.identifiers == {}
    assert compound.properties == {}
```

2. ML functionality
```python
def test_prediction():
    compound = MLCompound()
    result = compound.predict("binding")
    assert isinstance(result, dict)
```

3. Integration tests
```python
def test_full_pipeline():
    compound = ExportableCompound()
    compound.predict("binding")
    compound.enrich()
    compound.analyze()
    result = compound.export()
    assert isinstance(result, str)
```

## Success Criteria

### Code Quality
- [ ] No duplicate implementations
- [ ] Clear inheritance hierarchy
- [ ] Comprehensive docstrings
- [ ] Type hints

### Test Coverage
- [ ] Unit tests for all classes
- [ ] Integration tests
- [ ] Edge cases covered
- [ ] 90%+ coverage

### Documentation
- [ ] API documentation
- [ ] Usage examples
- [ ] Migration guide
- [ ] Architecture docs

## Commands

### Setup
```bash
# Create structure
./scripts/setup_models.sh

# Run tests
pytest tests/models/compound/

# Check coverage
pytest --cov=binding_data_processor/models/compound/

# Build docs
cd docs && make html
```

### Development
```bash
# Run specific tests
pytest tests/models/compound/base/test_core.py -v

# Run linters
flake8 binding_data_processor/models/compound/
mypy binding_data_processor/models/compound/

# Format code
black binding_data_processor/models/compound/
```

## Notes

1. Keep backward compatibility during migration
2. Update imports gradually
3. Run tests after each change
4. Update documentation as you go
5. Monitor performance impacts
