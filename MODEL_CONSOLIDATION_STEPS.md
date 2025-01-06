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
