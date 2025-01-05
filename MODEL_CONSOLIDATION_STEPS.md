# Model Consolidation Steps

## Overview

The codebase currently has duplicate model implementations:
1. models/compound/ - Core compound models
2. models/psychopharm/ - Psychopharm-specific models

These need to be consolidated into a single, cohesive model structure.

## Current Structure

### Compound Models
```
models/compound/
├── base.py          # Base compound model
├── ml.py           # ML functionality
├── enrichment.py   # Web enrichment
├── analysis.py     # Analysis tools
└── export.py       # Export features
```

### Psychopharm Models
```
models/psychopharm/
├── base.py         # Base psychopharm model
├── binding.py      # Binding predictions
├── activity.py     # Activity analysis
├── safety.py       # Safety assessment
└── enrichment.py   # Web enrichment
```

## Target Structure

```
models/compound/
├── base/
│   ├── __init__.py
│   ├── core.py        # Core model
│   ├── mixins.py      # Shared mixins
│   └── types.py       # Type definitions
├── analysis/
│   ├── __init__.py
│   ├── binding.py     # Binding analysis
│   ├── activity.py    # Activity analysis
│   ├── safety.py      # Safety analysis
│   └── properties.py  # Property analysis
├── ml/
│   ├── __init__.py
│   ├── predictors.py  # ML models
│   ├── features.py    # Feature extraction
│   └── training.py    # Model training
├── enrichment/
│   ├── __init__.py
│   ├── web.py        # Web enrichment
│   ├── community.py  # Community data
│   └── social.py     # Social data
└── export/
    ├── __init__.py
    ├── formats.py    # Export formats
    └── validation.py # Export validation
```

## Step-by-Step Plan

### 1. Create New Structure
```bash
# Create directories
mkdir -p models/compound/{base,analysis,ml,enrichment,export}

# Create __init__.py files
touch models/compound/{base,analysis,ml,enrichment,export}/__init__.py
```

### 2. Move Base Models
```python
# In models/compound/base/core.py
from typing import Dict, List, Optional

class CompoundBase:
    """Base compound model with core functionality."""
    def __init__(self):
        self.data = {}
        self.predictions = {}
        self.analysis = {}
```

### 3. Add Mixins
```python
# In models/compound/base/mixins.py
class MLMixin:
    """ML functionality mixin."""
    def predict(self):
        pass

class WebMixin:
    """Web enrichment mixin."""
    def enrich(self):
        pass
```

### 4. Merge Analysis
```python
# In models/compound/analysis/binding.py
class BindingAnalysis:
    """Binding affinity analysis."""
    def analyze_binding(self):
        pass

# In models/compound/analysis/activity.py
class ActivityAnalysis:
    """Activity analysis."""
    def analyze_activity(self):
        pass
```

### 5. Consolidate ML
```python
# In models/compound/ml/predictors.py
class BindingPredictor:
    """Binding affinity prediction."""
    def predict(self):
        pass

class ActivityPredictor:
    """Activity prediction."""
    def predict(self):
        pass
```

### 6. Merge Enrichment
```python
# In models/compound/enrichment/web.py
class WebEnrichment:
    """Web data enrichment."""
    def enrich(self):
        pass

# In models/compound/enrichment/community.py
class CommunityEnrichment:
    """Community data enrichment."""
    def enrich(self):
        pass
```

### 7. Update Imports
```python
# Update all imports to use new structure
from models.compound.base.core import CompoundBase
from models.compound.analysis.binding import BindingAnalysis
from models.compound.ml.predictors import BindingPredictor
```

### 8. Add Tests
```python
# In tests/models/compound/test_core.py
def test_compound_base():
    compound = CompoundBase()
    assert compound.data == {}

# In tests/models/compound/test_analysis.py
def test_binding_analysis():
    analysis = BindingAnalysis()
    result = analysis.analyze_binding()
    assert result is not None
```

## Validation Steps

### 1. Code Quality
- [ ] Run linters
- [ ] Run type checks
- [ ] Run tests
- [ ] Check coverage

### 2. Functionality
- [ ] Test core features
- [ ] Test ML models
- [ ] Test analysis
- [ ] Test enrichment

### 3. Integration
- [ ] Test pipeline
- [ ] Test web app
- [ ] Test exports
- [ ] Test imports

## Success Criteria

### 1. Code Structure
- Single source of truth for models
- Clear separation of concerns
- No duplicate code
- Type safety

### 2. Functionality
- All existing features preserved
- All tests passing
- No regressions
- Full coverage

### 3. Documentation
- Updated docstrings
- Updated README
- Updated examples
- Updated guides

## Next Steps

1. Create new directory structure
2. Move and merge models
3. Update imports
4. Add tests
5. Validate functionality
6. Update documentation

## Timeline

### Day 1
- Create structure
- Move base models
- Add mixins

### Day 2
- Merge analysis
- Consolidate ML
- Merge enrichment

### Day 3
- Update imports
- Add tests
- Validate

### Day 4
- Update docs
- Final testing
- Deploy
