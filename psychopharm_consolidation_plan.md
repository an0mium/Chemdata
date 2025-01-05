# Psychopharm Model Consolidation Plan

## Overview

Current State:
- Core psychopharm models in models/psychopharm/
- Legacy models in models/compound/
- Duplicate functionality needs merging
- ML features need integration

## Step 1: Base Model Consolidation

### Target: models/psychopharm/base.py
1. Merge from:
   - models/compound_base.py
   - models/mixins.py
   - models/compound/base/core.py

```python
# models/psychopharm/base.py
from typing import Dict, Any, Optional
from dataclasses import dataclass, field

@dataclass
class PsychopharmBase:
    """Base class for psychopharmacological compounds."""
    
    identifiers: Dict[str, str] = field(default_factory=dict)
    properties: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    # Core validation
    def validate(self) -> bool:
        """Validate compound data."""
        return all([
            self._validate_identifiers(),
            self._validate_properties(),
            self._validate_metadata()
        ])
    
    # Enhanced serialization
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary with all data."""
        return {
            "identifiers": self.identifiers,
            "properties": self.properties,
            "metadata": self.metadata,
            "predictions": self.get_all_predictions(),
            "web_data": self.get_all_web_data()
        }
```

## Step 2: ML Integration

### Target: models/psychopharm/binding.py
1. Merge from:
   - models/compound_ml.py
   - models/compound/ml/predictors.py

```python
# models/psychopharm/binding.py
from typing import Dict, Tuple, Optional
import numpy as np

class BindingPredictor:
    """Enhanced binding prediction with uncertainty."""
    
    def predict_binding(self, 
                       receptor: str,
                       include_uncertainty: bool = True
                      ) -> Tuple[float, float]:
        """Predict binding affinity with uncertainty."""
        prediction = self._base_prediction(receptor)
        uncertainty = self._estimate_uncertainty(receptor)
        return prediction, uncertainty
    
    def predict_ensemble(self,
                        receptor: str,
                        n_models: int = 5
                       ) -> Dict[str, Any]:
        """Get ensemble predictions."""
        predictions = []
        for _ in range(n_models):
            pred = self._single_model_prediction(receptor)
            predictions.append(pred)
        
        return {
            "mean": np.mean(predictions),
            "std": np.std(predictions),
            "individual": predictions
        }
```

## Step 3: Web Enrichment

### Target: models/psychopharm/enrichment.py
1. Merge from:
   - models/compound_enrichment.py
   - models/compound/enrichment/web.py

```python
# models/psychopharm/enrichment.py
from typing import Dict, Any, Optional
import aiohttp

class WebEnrichment:
    """Enhanced web data enrichment."""
    
    async def enrich_compound(self,
                            sources: Optional[List[str]] = None
                           ) -> Dict[str, Any]:
        """Enrich compound with web data."""
        sources = sources or ["community", "social", "literature"]
        
        async with aiohttp.ClientSession() as session:
            tasks = []
            for source in sources:
                task = self._fetch_source_data(session, source)
                tasks.append(task)
            
            results = await asyncio.gather(*tasks)
            
        return self._merge_results(results)
```

## Step 4: Analysis Integration

### Target: models/psychopharm/analysis.py
1. Merge from:
   - models/compound_analysis.py
   - models/compound/analysis/*.py

```python
# models/psychopharm/analysis.py
from typing import Dict, Any, Optional

class CompoundAnalyzer:
    """Enhanced compound analysis."""
    
    def analyze_compound(self,
                        analysis_types: Optional[List[str]] = None
                       ) -> Dict[str, Any]:
        """Comprehensive compound analysis."""
        analysis_types = analysis_types or [
            "binding",
            "activity",
            "safety",
            "properties"
        ]
        
        results = {}
        for analysis_type in analysis_types:
            method = getattr(self, f"_analyze_{analysis_type}")
            results[analysis_type] = method()
            
        return results
```

## Step 5: Testing

### Test Structure
```
tests/models/psychopharm/
├── test_base.py
├── test_binding.py
├── test_enrichment.py
└── test_analysis.py
```

### Example Test Cases
```python
# tests/models/psychopharm/test_binding.py
def test_binding_prediction_with_uncertainty():
    """Test binding prediction with uncertainty."""
    compound = PsychopharmCompound()
    prediction, uncertainty = compound.predict_binding("5HT2A")
    
    assert 0 <= prediction <= 1
    assert 0 <= uncertainty <= 1

def test_ensemble_prediction():
    """Test ensemble prediction."""
    compound = PsychopharmCompound()
    result = compound.predict_ensemble("5HT2A", n_models=5)
    
    assert "mean" in result
    assert "std" in result
    assert len(result["individual"]) == 5
```

## Migration Steps

1. Create New Structure
```bash
mkdir -p models/psychopharm/{base,binding,enrichment,analysis}
touch models/psychopharm/{base,binding,enrichment,analysis}/__init__.py
```

2. Move Files
```bash
mv models/compound_base.py models/psychopharm/base.py
mv models/compound_ml.py models/psychopharm/binding.py
mv models/compound_enrichment.py models/psychopharm/enrichment.py
mv models/compound_analysis.py models/psychopharm/analysis.py
```

3. Update Imports
```bash
find . -type f -name "*.py" -exec sed -i '' \
    's/from models.compound/from models.psychopharm/g' {} +
```

4. Run Tests
```bash
pytest tests/models/psychopharm/
```

## Success Criteria

1. Code Quality
- [ ] No duplicate code
- [ ] Clear inheritance
- [ ] Type hints
- [ ] Docstrings

2. Functionality
- [ ] All features preserved
- [ ] ML integration
- [ ] Web enrichment
- [ ] Analysis tools

3. Testing
- [ ] Unit tests pass
- [ ] Integration tests pass
- [ ] 90%+ coverage
- [ ] Edge cases covered

## Next Steps

1. Execute migration plan
2. Update documentation
3. Run full test suite
4. Clean up old files
