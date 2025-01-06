# Model Consolidation Steps

## Overview

The goal is to merge PsychopharmBase features into CompoundData while maintaining clean inheritance and avoiding duplication.

## 1. Core Data Model (CompoundData)

### Current Features ✓
- Basic identifiers
- Chemical properties
- Database IDs
- Common names
- Binding & Activity data
- ML predictions
- Web enrichment
- Analysis features

### Features to Add

1. From PsychopharmBase
- Enhanced Social Media Data
```python
social_media_data: Dict[str, List[Dict]] = {
    "reddit": [],    # [{subreddit, title, url, score, date}]
    "twitter": [],   # [{user, text, url, date}]
    "bluesky": [],   # [{user, text, url, date}]
    "discord": [],   # [{server, channel, text, date}]
}
```

2. Enhanced Patent Data
```python
patent_data: Dict[str, Dict] = {
    "numbers": Set[str],
    "titles": Dict[str, str],      # number -> title
    "abstracts": Dict[str, str],   # number -> abstract
    "claims": Dict[str, List[str]], # number -> claims
    "citations": Dict[str, List[str]] # number -> cited by
}
```

3. Enhanced Literature Data
```python
literature_data: Dict[str, Dict] = {
    "pubmed_ids": Set[str],
    "titles": Dict[str, str],       # pmid -> title
    "abstracts": Dict[str, str],    # pmid -> abstract
    "citations": Dict[str, List[str]], # pmid -> cited by
    "keywords": Dict[str, List[str]]  # pmid -> keywords
}
```

4. Enhanced Community Data
```python
community_data: Dict[str, Dict] = {
    "psychonaut": {
        "url": Optional[str],
        "data": Dict
    },
    "erowid": {
        "url": Optional[str],
        "data": Dict
    },
    "tripsit": {
        "url": Optional[str],
        "data": Dict
    }
}
```

5. From EnhancedCompound
```python
# Timing data
onset_times: Dict[str, TimeRange] = field(default_factory=dict)  # Route -> Time
duration_times: Dict[str, TimeRange] = field(default_factory=dict)  # Route -> Time

# Scoring data
effect_scores: Dict[str, EffectScore] = field(default_factory=dict)  # Effect -> Score
risk_scores: Dict[str, RiskScore] = field(default_factory=dict)  # Risk -> Score
receptor_bindings: Dict[str, ReceptorBinding] = field(default_factory=dict)

# Analysis data
analysis_data: Dict[str, Dict] = {
    "binding_profiles": Dict[str, ReceptorBinding],
    "binding_confidence": Dict[str, float],
    "binding_sources": Dict[str, List[str]],
    "activity_profiles": Dict[str, EffectScore],
    "activity_confidence": Dict[str, float],
    "activity_sources": Dict[str, List[str]],
    "risk_profiles": Dict[str, RiskScore],
    "risk_confidence": Dict[str, float],
    "risk_sources": Dict[str, List[str]],
    "dose_ranges": Dict[str, DoseRange],
    "time_ranges": Dict[str, TimeRange],
    "dose_confidence": Dict[str, float]
}
```

6. From Mixins
```python
# Base functionality
BaseMixin:
    - Core validation
    - Serialization
    - Type checking
    - Error handling

# Web enrichment
WebEnrichmentMixin:
    - Web scraping
    - Data extraction
    - Rate limiting
    - Error recovery

# ML predictions
PredictionsMixin:
    - Feature extraction
    - Model integration
    - Prediction caching
    - Uncertainty estimation

# Analysis capabilities
AnalysisMixin:
    - Data analysis
    - Visualization
    - Statistical tests
    - Report generation
```

## 2. Implementation Steps

1. Integrate Mixins
- Add mixin classes to CompoundData
- Ensure proper initialization order
- Handle method conflicts
- Maintain clean inheritance

2. Update Type Definitions
- Add new type definitions to compound/types.py
- Ensure consistent type usage across modules

2. Update Core Model
- Add new fields to CompoundData
- Add validation methods for new fields
- Add helper methods for data access

3. Update Analysis Methods
- Add methods for analyzing new data types
- Add methods for merging analysis results
- Add methods for exporting analysis

4. Update Web Enrichment
- Add methods for enriching new data types
- Add methods for merging enrichment results
- Add methods for validating enrichment

5. Update ML Integration
- Add methods for feature extraction from new data
- Add methods for prediction using new features
- Add methods for ensemble predictions

## 3. Migration Steps

1. Core Data Migration
```python
# In compound/base/core.py
class CompoundData:
    # Add new fields
    social_media_data: Dict[str, List[Dict]] = field(default_factory=dict)
    patent_data: Dict[str, Dict] = field(default_factory=dict)
    literature_data: Dict[str, Dict] = field(default_factory=dict)
    community_data: Dict[str, Dict] = field(default_factory=dict)
    analysis_data: Dict[str, Dict] = field(default_factory=dict)
```

2. Validation Migration
```python
def _validate(self):
    """Validate all compound data."""
    errors = []
    
    # Base validation
    errors.extend(self._validate_base_data())
    
    # Chemical validation
    if self.smiles:
        is_valid, error = validate_smiles(self.smiles)
        if not is_valid:
            errors.append(error)
            
    if self.inchi:
        is_valid, error = validate_inchi(self.inchi)
        if not is_valid:
            errors.append(error)
            
    if self.cas_number:
        is_valid, error = validate_cas_number(self.cas_number)
        if not is_valid:
            errors.append(error)
    
    # Property validation
    if self.properties:
        for prop, value in self.properties.items():
            if prop == "molecular_weight":
                is_valid, error = validate_molecular_weight(value)
                if not is_valid:
                    errors.append(error)
            elif prop == "logp":
                is_valid, error = validate_logp(value)
                if not is_valid:
                    errors.append(error)
    
    # Enrichment validation
    errors.extend(self._validate_social_data())
    errors.extend(self._validate_patent_data())
    errors.extend(self._validate_literature_data())
    errors.extend(self._validate_community_data())
    
    # Analysis validation
    errors.extend(self._validate_analysis_data())
    
    # Timing validation
    for route, time_range in self.onset_times.items():
        if not isinstance(time_range, TimeRange):
            errors.append(f"Invalid onset time range for route {route}")
            
    for route, time_range in self.duration_times.items():
        if not isinstance(time_range, TimeRange):
            errors.append(f"Invalid duration time range for route {route}")
    
    # Scoring validation
    for effect, score in self.effect_scores.items():
        if not isinstance(score, EffectScore):
            errors.append(f"Invalid effect score for {effect}")
            
    for risk, score in self.risk_scores.items():
        if not isinstance(score, RiskScore):
            errors.append(f"Invalid risk score for {risk}")
            
    for receptor, binding in self.receptor_bindings.items():
        if not isinstance(binding, ReceptorBinding):
            errors.append(f"Invalid receptor binding for {receptor}")
    
    # Prediction validation
    if self.predictions:
        for pred_type, pred in self.predictions.items():
            if pred.confidence is not None:
                is_valid, error = validate_confidence(pred.confidence)
                if not is_valid:
                    errors.append(f"Invalid confidence for {pred_type}: {error}")
                    
            if pred.probability is not None:
                is_valid, error = validate_probability(pred.probability)
                if not is_valid:
                    errors.append(f"Invalid probability for {pred_type}: {error}")
    
    # Raise if any errors
    if errors:
        raise ValidationError("\n".join(errors))
```

3. Additional Validation Functions
```python
def validate_time_range(time_range: TimeRange) -> Tuple[bool, Optional[str]]:
    """Validate time range values."""
    if not isinstance(time_range, TimeRange):
        return False, "Invalid time range type"
        
    if time_range.min < 0:
        return False, "Minimum time cannot be negative"
        
    if time_range.max <= time_range.min:
        return False, "Maximum time must be greater than minimum"
        
    return True, None

def validate_effect_score(score: EffectScore) -> Tuple[bool, Optional[str]]:
    """Validate effect score values."""
    if not isinstance(score, EffectScore):
        return False, "Invalid effect score type"
        
    if not 0 <= score.value <= 1:
        return False, "Score value must be between 0 and 1"
        
    if score.confidence is not None:
        is_valid, error = validate_confidence(score.confidence)
        if not is_valid:
            return False, error
            
    return True, None

def validate_risk_score(score: RiskScore) -> Tuple[bool, Optional[str]]:
    """Validate risk score values."""
    if not isinstance(score, RiskScore):
        return False, "Invalid risk score type"
        
    if not 0 <= score.value <= 1:
        return False, "Score value must be between 0 and 1"
        
    if score.confidence is not None:
        is_valid, error = validate_confidence(score.confidence)
        if not is_valid:
            return False, error
            
    return True, None

def validate_receptor_binding(binding: ReceptorBinding) -> Tuple[bool, Optional[str]]:
    """Validate receptor binding values."""
    if not isinstance(binding, ReceptorBinding):
        return False, "Invalid receptor binding type"
        
    if binding.affinity is not None:
        is_valid, error = validate_affinity_value(binding.affinity)
        if not is_valid:
            return False, error
            
    if binding.affinity_type is not None:
        is_valid, error = validate_affinity_type(binding.affinity_type)
        if not is_valid:
            return False, error
            
    if binding.affinity_unit is not None:
        is_valid, error = validate_affinity_unit(binding.affinity_unit)
        if not is_valid:
            return False, error
            
    return True, None
```

3. Analysis Migration
```python
def analyze_social_data(self):
    """Analyze social media mentions."""
    pass

def analyze_patent_data(self):
    """Analyze patent information."""
    pass

def analyze_literature_data(self):
    """Analyze scientific literature."""
    pass

def analyze_community_data(self):
    """Analyze community reports."""
    pass
```

4. Enrichment Migration
```python
def enrich_social_data(self):
    """Enrich with social media data."""
    pass

def enrich_patent_data(self):
    """Enrich with patent data."""
    pass

def enrich_literature_data(self):
    """Enrich with literature data."""
    pass

def enrich_community_data(self):
    """Enrich with community data."""
    pass
```

## 4. Testing Steps

1. Unit Tests
- Add tests for new fields
- Add tests for validation
- Add tests for analysis
- Add tests for enrichment

2. Integration Tests
- Add tests for data merging
- Add tests for full pipeline
- Add tests for web interface

3. Performance Tests
- Add tests for large datasets
- Add tests for concurrent access
- Add tests for memory usage

## 5. Documentation Steps

1. API Documentation
- Document new fields
- Document new methods
- Document validation rules
- Document analysis features

2. Usage Examples
- Add examples for new features
- Add examples for analysis
- Add examples for enrichment
- Add examples for export

3. Migration Guide
- Document migration process
- Document breaking changes
- Document new capabilities
- Document best practices

## Success Criteria

1. Code Quality
- All tests passing
- No duplicate code
- Clear documentation
- Type hints complete

2. Functionality
- All features preserved
- Enhanced capabilities
- Good performance
- Clean interfaces

3. Documentation
- Updated API docs
- Clear examples
- Migration guide
- Best practices

4. Testing
- High coverage
- Integration tests
- Performance tests
- Error cases


# Model Consolidation Steps

## Overview

The codebase currently has several model implementations that need to be consolidated:
1. Legacy models in deprecated/ directory
2. Enhanced models in compound/ directory
3. Psychopharm-specific models in psychopharm/ directory

## Current Status

### Core Infrastructure
1. Directory Structure ✓
   ```
   binding_data_processor/models/compound/
   ├── base/
   │   ├── core.py        # Core model
   │   ├── mixins.py      # Shared mixins
   │   ├── validation.py  # Validation logic
   │   └── types.py       # Type definitions ✓
   ├── analysis/
   │   ├── binding/       # Binding analysis
   │   ├── activity/      # Activity analysis
   │   ├── safety/        # Safety analysis
   │   └── properties/    # Property analysis
   ├── ml/
   │   ├── features.py    # Feature extraction
   │   ├── training.py    # Model training
   │   ├── predictors.py  # Model predictors
   │   └── ensemble.py    # Ensemble models
   ├── enrichment/
   │   ├── web.py        # Web enrichment
   │   ├── community.py  # Community data
   │   └── social.py     # Social data
   └── export/
       ├── formats.py    # Export formats
       └── validation.py # Export validation
   ```

2. Type System Enhancement ✓
   - Merged compound types from both versions
   - Added comprehensive legal status options
   - Extended psychoactive classifications
   - Enhanced nootropic mechanisms
   - Added BBB permeability levels
   - Preserved all data structures and type aliases

3. Legacy Code Migration ✓
   ```
   binding_data_processor/models/compound/
   ├── base/
   │   ├── core.py        # Core model
   │   ├── mixins.py      # Shared mixins
   │   ├── validation.py  # Validation logic
   │   └── types.py       # Type definitions
   ├── analysis/
   │   ├── binding/       # Binding analysis
   │   ├── activity/      # Activity analysis
   │   ├── safety/        # Safety analysis
   │   └── properties/    # Property analysis
   ├── ml/
   │   ├── features.py    # Feature extraction
   │   ├── training.py    # Model training
   │   ├── predictors.py  # Model predictors
   │   └── ensemble.py    # Ensemble models
   ├── enrichment/
   │   ├── web.py        # Web enrichment
   │   ├── community.py  # Community data
   │   └── social.py     # Social data
   └── export/
       ├── formats.py    # Export formats
       └── validation.py # Export validation
   ```

2. Legacy Code Migration ✓
   - Created deprecated/ directory structure
   - Moved legacy files to appropriate subdirectories
   - Preserved backup files in backups/
   - Maintained documentation of migrated code

3. Test Framework Enhancement ✓
   - Enhanced nootropic prediction test suite
   - Added comprehensive test compounds:
     ```python
     # Known nootropics with strong evidence
     piracetam = CompoundData(
         name="Piracetam",
         smiles="O=C1N(CC(=O)N(CC1)CC)CC",
         cas_number="7491-74-9",
     )
     aniracetam = CompoundData(
         name="Aniracetam", 
         smiles="O=C1N(C(=O)CC(N1)c1ccccc1)CC",
         cas_number="72432-10-1",
     )
     modafinil = CompoundData(
         name="Modafinil",
         smiles="CC(=O)C(CS(=O)(=O)C)NC(=O)C",
         cas_number="68693-11-8",
     )
     ```
   - Added test coverage for:
     - Basic prediction functionality
     - Mechanism prediction
     - Cognitive effects
     - Side effects
     - Safety analysis
     - Literature analysis
     - Model persistence
     - Error handling

## Model Architecture

### Core Models
1. CompoundData
   ```python
   class CompoundData:
       """Core compound data model."""
       name: str
       smiles: str
       cas_number: str
       properties: Dict[str, Any]
       predictions: Dict[str, Any]
       web_data: Dict[str, Any]
       analysis: Dict[str, Any]
       
       def validate(self) -> ValidationResult:
           """Validate compound data."""
           pass
   ```

2. TargetData
   ```python
   class TargetData:
       """Target binding data model."""
       name: str
       type: str
       organism: str
       affinity: float
       conditions: Dict[str, Any]
       
       def merge(self, other: "TargetData") -> "TargetData":
           """Merge target data."""
           pass
   ```

3. PredictionData
   ```python
   class PredictionData:
       """ML prediction data model."""
       model: str
       value: float
       confidence: float
       metadata: Dict[str, Any]
       
       def to_dict(self) -> Dict[str, Any]:
           """Convert to dictionary."""
           pass
   ```

### Type System
```python
class CompoundType(Enum):
    SMALL_MOLECULE = "small_molecule"
    NATURAL_PRODUCT = "natural_product"
    PHARMACEUTICAL = "pharmaceutical"

class PsychoactiveClass(Enum):
    PSYCHEDELIC = "psychedelic"
    DISSOCIATIVE = "dissociative"
    STIMULANT = "stimulant"
    DEPRESSANT = "depressant"
    NOOTROPIC = "nootropic"

class NootropicMechanism(Enum):
    MEMORY_ENHANCEMENT = "memory_enhancement"
    FOCUS_IMPROVEMENT = "focus_improvement"
    NEUROPROTECTION = "neuroprotection"
    NEUROPLASTICITY = "neuroplasticity"
    CHOLINERGIC = "cholinergic"
    GLUTAMATERGIC = "glutamatergic"
    DOPAMINERGIC = "dopaminergic"
    SEROTONERGIC = "serotonergic"
    NOOTROPIC_SYNERGY = "nootropic_synergy"
    COGNITIVE_MODULATION = "cognitive_modulation"
    BRAIN_METABOLISM = "brain_metabolism"
```

### Class Hierarchy
```python
# Base class with core functionality
class CompoundBase(ValidationMixin, SerializationMixin):
    """Base compound class."""
    name: str
    smiles: str
    cas_number: str
    properties: Dict[str, Any]

# ML-enabled compound
class MLCompound(CompoundBase):
    """ML-enabled compound functionality."""
    predictions: Dict[str, PredictionResult]
    features: Dict[str, Any]
    history: List[Dict[str, Any]]

# Web-enriched compound
class EnrichedCompound(MLCompound):
    """Web-enriched compound functionality."""
    web_data: WebData
    literature_data: LiteratureData
    social_data: Dict[str, Any]
    patent_data: List[Dict[str, Any]]

# Analyzed compound
class AnalyzedCompound(EnrichedCompound):
    """Analyzed compound functionality."""
    analysis: Dict[str, Any]
    metrics: Dict[str, float]
    visualizations: Dict[str, Any]
```

### Analysis Models
1. BBBPrediction ✓
   ```python
   class BBBPrediction:
       """BBB permeability prediction."""
       permeability: float
       confidence: float
       transporters: Dict[str, float]
       mechanisms: Dict[str, float]
       
       def get_primary_mechanisms(self) -> List[str]:
           """Get primary transport mechanisms."""
           pass
   ```

2. ToxicityPrediction ✓
   ```python
   class ToxicityPrediction:
       """Toxicity prediction."""
       toxicity_score: float
       confidence: float
       mechanisms: Dict[str, float]
       organ_effects: Dict[str, float]
       
       def get_risk_level(self) -> str:
           """Get overall risk level."""
           pass
   ```

3. AbusePrediction ✓
   ```python
   class AbusePrediction:
       """Abuse potential prediction."""
       abuse_score: float
       confidence: float
       reward_pathways: Dict[str, float]
       tolerance_profile: Dict[str, float]
       
       def get_risk_category(self) -> str:
           """Get abuse risk category."""
           pass
   ```

### Web Models
1. WebData ✓
   ```python
   class WebData:
       """Web-enriched compound data."""
       community_reports: List[Dict[str, Any]]
       literature_data: List[Dict[str, Any]]
       social_data: Dict[str, Any]
       patent_data: List[Dict[str, Any]]
       
       def merge(self, other: "WebData") -> "WebData":
           """Merge web data."""
           pass
   ```

2. LiteratureData ✓
   ```python
   class LiteratureData:
       """Literature analysis data."""
       papers: List[Dict[str, Any]]
       findings: List[Dict[str, Any]]
       mechanisms: Dict[str, List[str]]
       safety_data: Dict[str, Any]
       
       def summarize(self) -> Dict[str, Any]:
           """Summarize literature data."""
           pass
   ```

3. PatentData ✓
   ```python
   class PatentData:
       """Patent analysis data."""
       patents: List[Dict[str, Any]]
       compounds: List[Dict[str, Any]]
       activities: Dict[str, Any]
       synthesis: Dict[str, Any]
       
       def extract_compounds(self) -> List[CompoundData]:
           """Extract compounds from patents."""
           pass
   ```

### ML Models
1. EnsemblePredictor
   ```python
   class EnsemblePredictor:
       """Ensemble ML predictor."""
       predictors: List[BasePredictor]
       weights: Dict[str, float]
       
       def predict(self, compound: CompoundData) -> PredictionResult:
           """Get ensemble prediction."""
           pass
           
       def update_weights(self, metrics: Dict[str, float]):
           """Update predictor weights."""
           pass
   ```

2. UncertaintyEstimator
   ```python
   class UncertaintyEstimator:
       """Prediction uncertainty estimator."""
       def estimate(self, predictions: List[PredictionResult]) -> float:
           """Estimate prediction uncertainty."""
           pass
           
       def get_confidence_interval(self, 
           predictions: List[PredictionResult]
       ) -> Tuple[float, float]:
           """Get confidence interval."""
           pass
   ```

3. FeatureImportance
   ```python
   class FeatureImportance:
       """Feature importance analyzer."""
       def analyze(self, model: BasePredictor) -> Dict[str, float]:
           """Analyze feature importance."""
           pass
           
       def explain_prediction(self,
           model: BasePredictor,
           compound: CompoundData
       ) -> Dict[str, float]:
           """Explain specific prediction."""
           pass
   ```

### Integration Models
1. CompoundCollection
   ```python
   class CompoundCollection:
       """Collection of compounds."""
       compounds: List[EnrichedCompound]
       metadata: Dict[str, Any]
       
       def filter(self, criteria: Dict[str, Any]) -> "CompoundCollection":
           """Filter compounds."""
           pass
           
       def sort(self, key: str, reverse: bool = False) -> "CompoundCollection":
           """Sort compounds."""
           pass
           
       def export(self, format: str) -> bytes:
           """Export collection."""
           pass
   ```

2. AnalysisResult
   ```python
   class AnalysisResult:
       """Compound analysis result."""
       compound: EnrichedCompound
       predictions: Dict[str, PredictionResult]
       metrics: Dict[str, float]
       visualizations: Dict[str, Any]
       
       def to_html(self) -> str:
           """Convert to HTML."""
           pass
           
       def to_json(self) -> Dict[str, Any]:
           """Convert to JSON."""
           pass
   ```

3. ValidationResult
   ```python
   class ValidationResult:
       """Validation result."""
       is_valid: bool
       errors: List[str]
       warnings: List[str]
       suggestions: List[str]
       
       def to_dict(self) -> Dict[str, Any]:
           """Convert to dictionary."""
           pass
   ```

## Next Steps

### Day 1: Core Integration
1. BBB Integration
   - Integrate BBB predictor with nootropic analysis
   - Add BBB-specific test cases
   - Validate predictions against known data
   - Add integration tests

2. Model Enhancement
   - Add ensemble models for:
     - Nootropic prediction
     - BBB prediction
     - Toxicity prediction
     - Abuse potential prediction
   - Implement cross-validation
   - Add feature importance analysis
   - Add uncertainty estimation

3. Data Validation
   - Add schema validation for:
     - Input compounds
     - Prediction results
     - Analysis results
     - Web data
   - Implement data quality checks
   - Add consistency validation
   - Improve error handling

### Day 2: Web Enhancement
1. LLM Integration
   - Add SciBERT for entity extraction
   - Add PubMedBert for relevance scoring
   - Add text classification
   - Add relationship extraction

2. Web Scraping
   - Add rate limiting
   - Add proxy support
   - Add error recovery
   - Add data validation

3. Monitoring
   - Implement performance monitoring
   - Add prediction tracking
   - Enhance error logging
   - Add usage analytics

### Day 3: Infrastructure
1. Code Structure
   - Move remaining legacy code to deprecated/
   - Update import statements
   - Fix circular dependencies
   - Add missing __init__ files

2. Documentation
   - Update docstrings
   - Add type hints
   - Add examples
   - Update README files

3. Testing
   - Add missing test cases
   - Improve test coverage
   - Add integration tests
   - Add performance tests

## Success Metrics

### Code Quality
- [ ] All tests passing
- [ ] >90% test coverage
- [ ] No circular imports
- [ ] Clean architecture
- [ ] All files under 700 lines
- [ ] No duplicate code
- [ ] Clear inheritance
- [ ] Type hints complete

### Documentation
- [ ] Complete docstrings
- [ ] Up-to-date READMEs
- [ ] Clear examples
- [ ] Good API docs
- [ ] Architecture docs
- [ ] Usage guides

### Performance
- [ ] Fast prediction times (<500ms)
- [ ] Efficient memory usage (<2GB)
- [ ] Good scalability
- [ ] Reliable caching
- [ ] Error recovery
- [ ] Monitoring

### Usability
- [ ] Clear interfaces
- [ ] Good error messages
- [ ] Helpful documentation
- [ ] Easy deployment
- [ ] Intuitive API
- [ ] Good UX

## Notes
1. Keep models focused and small
2. Use clear inheritance patterns
3. Add comprehensive validation
4. Handle errors gracefully
5. Monitor performance carefully
6. Document everything thoroughly
