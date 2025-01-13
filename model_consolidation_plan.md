# Model Consolidation Plan

## Overview

Current issues:
1. Duplicate model definitions between models/ and models/compound/
2. Psychopharm functionality not integrated
3. Analysis code spread across modules
4. Validation logic duplicated

## Step 1: Directory Structure

```
binding_data_processor/models/quantum/
├── base/
│   ├── __init__.py
│   ├── core.py          # Core quantum model
│   ├── electronic.py    # Electronic structure
│   └── criticality.py   # Phase transitions
├── analysis/
│   ├── __init__.py
│   ├── electronic.py    # Electronic analysis
│   ├── phase.py        # Phase analysis
│   └── scaling.py      # Scaling analysis
├── properties/
│   ├── __init__.py
│   ├── observables.py  # Quantum observables
│   └── correlation.py  # Correlation functions
└── tests/
    ├── test_core.py
    ├── test_electronic.py
    └── test_criticality.py

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

### Quantum Models
```python
# Base Classes
class QuantumBase:
    """Base quantum model."""
    def __init__(self):
        self.wavefunction = None
        self.density_matrix = None
        self.observables = {}

class ElectronicStructure(QuantumBase):
    """Electronic structure model."""
    def calculate_electronic_structure(self):
        pass

class CriticalityAnalysis(QuantumBase):
    """Phase transition analysis."""
    def analyze_criticality(self):
        pass

# Analysis Classes
class ElectronicAnalyzer:
    """Electronic structure analysis."""
    def analyze_electronic_structure(self, compound: ElectronicStructure):
        pass

class PhaseAnalyzer:
    """Phase transition analysis."""
    def analyze_phase_transition(self, compound: CriticalityAnalysis):
        pass

class ScalingAnalyzer:
    """Critical point scaling analysis."""
    def analyze_scaling(self, compound: CriticalityAnalysis):
        pass
```

### Compound Models

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

### Day 0: Quantum Setup
1. Create quantum directory structure
```bash
mkdir -p binding_data_processor/models/quantum/{base,analysis,properties,tests}
touch binding_data_processor/models/quantum/{base,analysis,properties,tests}/__init__.py
```

2. Create quantum base files
```bash
# Create base files
touch binding_data_processor/models/quantum/base/{core,electronic,criticality}.py
touch binding_data_processor/models/quantum/analysis/{electronic,phase,scaling}.py
touch binding_data_processor/models/quantum/properties/{observables,correlation}.py
```

3. Create quantum test files
```bash
touch binding_data_processor/models/quantum/tests/{test_core,test_electronic,test_criticality}.py
```

### Day 1: Setup

### Quantum Integration
1. Implement base classes
```python
# In quantum/base/core.py
class QuantumBase:
    """Base quantum model."""
    def __init__(self):
        self.wavefunction = None
        self.density_matrix = None
        self.observables = {}
```

2. Add electronic structure
```python
# In quantum/base/electronic.py
class ElectronicStructure(QuantumBase):
    """Electronic structure model."""
    def calculate_electronic_structure(self):
        pass
```

3. Add criticality analysis
```python
# In quantum/base/criticality.py
class CriticalityAnalysis(QuantumBase):
    """Phase transition analysis."""
    def analyze_criticality(self):
        pass
```
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

### Quantum Analysis
1. Implement electronic analysis
```python
# In quantum/analysis/electronic.py
class ElectronicAnalyzer:
    """Electronic structure analysis."""
    def analyze_electronic_structure(self, compound: ElectronicStructure):
        pass
```

2. Add phase analysis
```python
# In quantum/analysis/phase.py
class PhaseAnalyzer:
    """Phase transition analysis."""
    def analyze_phase_transition(self, compound: CriticalityAnalysis):
        pass
```

3. Add scaling analysis
```python
# In quantum/analysis/scaling.py
class ScalingAnalyzer:
    """Critical point scaling analysis."""
    def analyze_scaling(self, compound: CriticalityAnalysis):
        pass
```
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

### Quantum Tests
```
tests/models/quantum/
├── base/
│   ├── test_core.py
│   ├── test_electronic.py
│   └── test_criticality.py
├── analysis/
│   ├── test_electronic.py
│   ├── test_phase.py
│   └── test_scaling.py
└── properties/
    ├── test_observables.py
    └── test_correlation.py
```

### Test Cases
1. Quantum base functionality
```python
def test_quantum_initialization():
    quantum = QuantumBase()
    assert quantum.wavefunction is None
    assert quantum.density_matrix is None
    assert quantum.observables == {}
```

2. Electronic structure
```python
def test_electronic_structure():
    structure = ElectronicStructure()
    result = structure.calculate_electronic_structure()
    assert isinstance(result, dict)
```

3. Criticality analysis
```python
def test_criticality():
    analysis = CriticalityAnalysis()
    result = analysis.analyze_criticality()
    assert isinstance(result, dict)
```

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

### Quantum Integration
- [ ] Complete quantum model implementation
- [ ] Electronic structure analysis
- [ ] Phase transition detection
- [ ] Critical point analysis
- [ ] Scaling behavior analysis

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

### Quantum Setup
```bash
# Create quantum structure
./scripts/setup_quantum.sh

# Run quantum tests
pytest tests/models/quantum/

# Check quantum coverage
pytest --cov=binding_data_processor/models/quantum/
```

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
