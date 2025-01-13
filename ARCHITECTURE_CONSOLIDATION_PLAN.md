# Architecture Consolidation Plan

## Overview

This document outlines a comprehensive plan to consolidate and enhance the architecture of the binding_data_processor system, addressing current issues while preserving and optimizing valuable functionality.

## Implementation Phases

### Phase 1: Core Framework (Weeks 1-2)

#### 1. Enhanced Dynamics Analysis
```python
@dataclass
class DynamicsConfig:
    """Configuration for dynamics analysis."""
    # Analysis settings
    use_normal_modes: bool = True
    use_correlations: bool = True
    use_domain_motions: bool = True
    use_contacts: bool = True
    use_alphafold_confidence: bool = True
    
    # Analysis parameters
    n_modes: int = 10
    contact_cutoff: float = 8.0
    min_plddt: float = 70.0
    interface_cutoff: float = 5.0
    
    # Domain analysis
    analyze_hinge_regions: bool = True
    analyze_interfaces: bool = True
    analyze_correlations: bool = True
    
    # Performance
    cache_results: bool = True
    n_jobs: int = 1

class DynamicsAnalyzer(BaseAnalyzer):
    """Enhanced protein dynamics analyzer."""
    
    def __init__(self, config: Optional[DynamicsConfig] = None):
        self.config = config or DynamicsConfig()
        self.logger = logging.getLogger(self.__class__.__name__)
        self._cache = {}
        
    def analyze(self, structure: Structure) -> Dict[str, Any]:
        """Analyze protein dynamics."""
        try:
            results = {}
            
            # Basic dynamics analysis
            coords = self.get_coordinates()
            flexibility = self.analyze_flexibility()
            
            # Normal modes analysis
            if self.config.use_normal_modes:
                modes = self.analyze_normal_modes(
                    n_modes=self.config.n_modes
                )
                results["modes"] = modes
                
                # Calculate correlations
                if self.config.use_correlations:
                    correlations = self._get_correlation_matrix()
                    results["correlations"] = correlations
            
            # Domain motion analysis
            if self.config.use_domain_motions:
                domains = self._get_domain_assignments()
                domain_results = analyze_domain_motions(
                    structure,
                    domains,
                    self.config
                )
                results["domains"] = domain_results
                
                # Analyze interfaces and hinges
                if self.config.analyze_interfaces:
                    interfaces = self._analyze_interfaces(
                        structure,
                        domains,
                        domain_results["motions"]
                    )
                    results["interfaces"] = interfaces
                    
                if self.config.analyze_hinge_regions:
                    hinges = self._identify_hinge_regions(
                        structure,
                        domains,
                        domain_results["motions"]
                    )
                    results["hinges"] = hinges
            
            # Weight by AlphaFold confidence
            if self.config.use_alphafold_confidence:
                results = self._apply_confidence_weights(
                    results,
                    min_plddt=self.config.min_plddt
                )
                
            return results
            
        except Exception as e:
            self.logger.error(f"Analysis failed: {str(e)}")
            return {}
            
    def analyze_site_dynamics(
        self,
        site_residues: List[int],
        include_context: bool = True
    ) -> Dict[str, Any]:
        """Analyze binding site dynamics."""
        try:
            # Get site-specific dynamics
            site_dynamics = self._analyze_site_specific(site_residues)
            
            # Get flexibility profile
            flexibility = self.analyze_flexibility()
            site_flexibility = {
                res: flexibility["scores"].get(res, 0.0)
                for res in site_residues
            }
            
            # Get contact network
            contacts = self.analyze_contacts()
            site_contacts = self._get_site_contacts(
                site_residues,
                contacts["network"]
            )
            
            # Calculate metrics
            metrics = {
                "avg_flexibility": float(np.mean(list(site_flexibility.values()))),
                "contact_density": site_contacts["density"],
                "surface_exposure": site_contacts["exposure"]
            }
            
            results = {
                "dynamics": site_dynamics,
                "flexibility": site_flexibility,
                "contacts": site_contacts,
                "metrics": metrics
            }
            
            # Add domain context
            if include_context:
                results["context"] = self._analyze_domain_context(site_residues)
                
            return results
            
        except Exception as e:
            self.logger.error(f"Site analysis failed: {str(e)}")
            return {}
```

### Phase 2: Integration (Weeks 3-4)

#### 1. Enhanced Toxicity Integration
```python
@dataclass
class ToxicityConfig:
    """Configuration for toxicity analysis."""
    # Analysis settings
    use_dynamics: bool = True
    use_binding_site: bool = True
    use_domain_context: bool = True
    use_confidence_weighting: bool = True
    use_alphafold_confidence: bool = True
    use_enhanced_scoring: bool = True
    
    # Confidence thresholds
    min_confidence: float = 70.0
    min_plddt: float = 70.0
    min_binding_score: float = 0.5
    
    # Risk weights
    flexibility_weight: float = 1.0
    accessibility_weight: float = 1.0
    domain_motion_weight: float = 0.5
    interface_weight: float = 0.5
    
    # Scoring thresholds
    high_flexibility_threshold: float = 1.5
    high_accessibility_threshold: float = 0.8
    high_motion_threshold: float = 2.0
    high_interface_threshold: float = 5.0

class ToxicityAnalyzer(BaseAnalyzer):
    """Enhanced toxicity analyzer with dynamics integration."""
    
    def __init__(self, config: Optional[ToxicityConfig] = None):
        super().__init__()
        self.config = config or ToxicityConfig()
        self.dynamics = DynamicsAnalyzer()
        self.bbb = BBBPredictorEnhanced()
        
    def analyze_binding_site(
        self,
        structure: Structure,
        site_residues: List[int]
    ) -> Dict[str, Any]:
        """Analyze binding site with dynamics integration."""
        try:
            results = {}
            
            # Get dynamics analysis if enabled
            if self.config.use_dynamics:
                dynamics = self.dynamics.analyze_site_dynamics(
                    site_residues,
                    include_context=self.config.use_domain_context
                )
                results["dynamics"] = dynamics
                
                # Calculate risk factors from dynamics
                risk_factors = self._calculate_dynamic_risks(
                    dynamics,
                    site_residues
                )
                results["risk_factors"] = risk_factors
                
            # Get binding site properties
            if self.config.use_binding_site:
                binding = self._analyze_binding_site(
                    structure,
                    site_residues,
                    include_dynamics="dynamics" in results
                )
                results["binding"] = binding
                
                # Update risk factors with binding properties
                if "risk_factors" in results:
                    binding_risks = self._calculate_binding_risks(binding)
                    results["risk_factors"].update(binding_risks)
                    
            # Weight results by confidence if enabled
            if self.config.use_confidence_weighting:
                results = self._apply_confidence_weights(
                    results,
                    min_confidence=self.config.min_confidence
                )
                
            # Calculate overall toxicity score
            results["toxicity_score"] = self._calculate_toxicity_score(
                results.get("risk_factors", {}),
                weights={
                    "flexibility": self.config.flexibility_weight,
                    "accessibility": self.config.accessibility_weight,
                    "domain_motion": self.config.domain_motion_weight,
                    "interface": self.config.interface_weight
                }
            )
            
            return results
            
        except Exception as e:
            self.logger.error(f"Toxicity analysis failed: {str(e)}")
            return {}
            
    def _calculate_dynamic_risks(
        self,
        dynamics: Dict[str, Any],
        site_residues: List[int]
    ) -> Dict[str, float]:
        """Calculate risk factors from dynamics analysis."""
        risks = {}
        
        # Assess flexibility risk
        if "flexibility" in dynamics:
            flex_values = [
                dynamics["flexibility"].get(res, 0.0)
                for res in site_residues
            ]
            risks["flexibility"] = self._score_flexibility(
                mean=np.mean(flex_values),
                std=np.std(flex_values),
                threshold=self.config.high_flexibility_threshold
            )
            
        # Assess accessibility risk
        if "contacts" in dynamics:
            risks["accessibility"] = self._score_accessibility(
                exposure=dynamics["contacts"]["surface_exposure"],
                density=dynamics["contacts"]["contact_density"],
                threshold=self.config.high_accessibility_threshold
            )
            
        # Assess domain motion risk
        if "domain_context" in dynamics:
            context = dynamics["domain_context"]
            if "motions" in context:
                risks["domain_motion"] = self._score_domain_motion(
                    amplitudes=[
                        motion["variance"][0]
                        for motion in context["motions"].values()
                    ],
                    correlations=context.get("correlations", {}),
                    threshold=self.config.high_motion_threshold
                )
                
        # Assess interface risk
        if "contacts" in dynamics:
            risks["interface"] = self._score_interface(
                dynamics["contacts"],
                threshold=self.config.high_interface_threshold
            )
            
        return risks
        
    def _calculate_toxicity_score(
        self,
        risks: Dict[str, float],
        weights: Dict[str, float]
    ) -> float:
        """Calculate overall toxicity score."""
        if not risks:
            return 0.0
            
        weighted_sum = 0.0
        total_weight = 0.0
        
        for factor, score in risks.items():
            weight = weights.get(factor, 1.0)
            weighted_sum += score * weight
            total_weight += weight
            
        return weighted_sum / total_weight if total_weight > 0 else 0.0
```

### Phase 3: Testing & Validation (Weeks 5-6)

#### 1. Unit Tests
```python
def test_dynamics_analyzer():
    """Test comprehensive dynamics analysis."""
    config = DynamicsConfig(
        use_normal_modes=True,
        use_domain_motions=True,
        analyze_hinge_regions=True
    )
    analyzer = DynamicsAnalyzer(config)
    structure = load_test_structure()
    
    results = analyzer.analyze(structure)
    assert "modes" in results
    assert "domains" in results
    assert "hinges" in results
    
def test_toxicity_integration():
    """Test toxicity prediction with dynamics."""
    config = ToxicityConfig(
        use_dynamics=True,
        use_binding_site=True,
        use_domain_context=True
    )
    analyzer = ToxicityAnalyzer(config)
    structure = load_test_structure()
    site_residues = [10, 11, 12, 13]
    
    results = analyzer.analyze_binding_site(
        structure,
        site_residues
    )
    
    assert "dynamics" in results
    assert "risk_factors" in results
    assert "toxicity_score" in results
    assert 0.0 <= results["toxicity_score"] <= 1.0
```

## Success Metrics

1. Code Quality
- Zero circular imports
- >90% test coverage
- Comprehensive documentation
- Type hints and validation

2. Performance
- Cache hit rate >80%
- Analysis time <30s
- Memory usage <2GB
- Efficient data structures

3. Integration
- Clean interfaces
- Proper error handling
- Consistent patterns
- Good separation of concerns

4. Maintainability
- Clear documentation
- Modular design
- Easy configuration
- Simple deployment

## Next Steps

1. Infrastructure Setup
- Set up new directory structure
- Create base classes
- Implement core interfaces
- Add initial tests

2. Migration Planning
- Identify critical paths
- Plan phased rollout
- Create backup strategy
- Document procedures

3. Development
- Start with core infrastructure
- Move to domain logic
- Add integration layer
- Implement testing

4. Validation
- Run comprehensive tests
- Validate performance
- Check compatibility
- Verify functionality

5. Deployment
- Stage changes
- Monitor metrics
- Gather feedback
- Make adjustments

## Conclusion

This consolidation plan provides a comprehensive approach to:
1. Fix current architectural issues
2. Improve code organization
3. Enhance performance
4. Increase maintainability
5. Better developer experience

The phased implementation ensures minimal disruption while achieving significant improvements in code quality, performance, and maintainability.


# Architecture Consolidation Plan

## Detailed Implementation Analysis

### Toxicity Prediction Implementations

1. **psychopharm/predictors/toxicity.py**
   - Comprehensive ensemble-based implementation
   - Rich feature set including:
     * Cellular mechanisms (oxidative stress, mitochondrial toxicity, etc.)
     * Molecular mechanisms (reactive metabolites, protein binding, etc.)
     * Systemic effects (immune activation, organ failure, etc.)
   - Sophisticated model management with versioning
   - Extensive prediction history tracking
   - Comprehensive validation and metrics

2. **psychopharm/predictors/toxicity_enhanced.py**
   - Extends base toxicity predictor
   - Adds BBB permeability integration
   - Implements model calibration
   - Enhanced uncertainty estimation
   - Improved feature importance analysis
   - Better model versioning

3. **structure/ml/toxicity.py**
   - Simpler ML-based implementation
   - Focused on structural features
   - Basic model lifecycle management
   - Limited to core toxicity endpoints
   - Minimal infrastructure integration

### Key Architectural Differences

1. **Model Management**
   - psychopharm: Sophisticated ensemble management with versioning
   - structure/ml: Basic single model per endpoint
   - Enhanced: Adds calibration and uncertainty

2. **Feature Engineering**
   - psychopharm: Rich domain-specific features
   - structure/ml: Structure-based features only
   - Enhanced: Integrates BBB features

3. **Prediction Workflow**
   - psychopharm: Complex multi-stage prediction
   - structure/ml: Direct endpoint prediction
   - Enhanced: Adds confidence calibration

4. **Infrastructure Integration**
   - psychopharm: Deep integration with monitoring/caching
   - structure/ml: Minimal infrastructure usage
   - Enhanced: Adds cross-domain integration

## Structure ML Module Analysis

The structure/ml module and its parallel implementations should be consolidated following this detailed analysis and integration plan:

### 1. Activity Prediction Components

#### Current Parallel Implementations
- `processors/structure/ml/activity/` - Core activity prediction
- `processors/psychopharm/predictors/` - Psychopharmacology-specific predictors
- `models/compound/ml/predictors.py` - General compound predictors
- `pipeline/ml/` - ML pipeline implementations

#### Integration Strategy
1. Create unified activity prediction framework:
```
binding_data_processor/
  └── processors/
      └── structure/
          └── ml/
              └── activity/
                  ├── base.py         # Base prediction interfaces
                  ├── generic/        # General activity predictors
                  ├── psychopharm/    # Psychopharm-specific (moved from processors/psychopharm)
                  └── specialized/    # Domain-specific predictors
```

2. Standardize prediction interfaces:
- Define common base classes in `activity/base.py`
- Implement shared utilities and metrics
- Standardize model input/output formats
- Unify configuration handling

3. Migrate domain-specific logic:
- Move psychopharm predictors while preserving domain expertise
- Maintain backwards compatibility during transition
- Add proper abstraction layers for specialization

### 2. Protein Analysis Components

#### Current Parallel Implementations
- `processors/structure/ml/activity/protein/` - Core protein analysis
- `processors/structure/ml/activity/binding/structure/` - Binding site analysis
- `models/compound/analysis/` - Compound-protein interaction analysis

#### Integration Strategy
1. Consolidate protein analysis modules:
```
binding_data_processor/
  └── processors/
      └── structure/
          └── ml/
              └── protein/
                  ├── analysis/       # Core analysis capabilities
                  ├── binding/        # Binding site analysis
                  ├── dynamics/       # Protein dynamics
                  └── prediction/     # Structure prediction
```

2. Unify analysis interfaces:
- Create consistent API for protein analysis
- Standardize data structures and formats
- Implement shared validation and utilities

3. Enhance integration points:
- Connect binding analysis with dynamics
- Link structure prediction to analysis
- Enable compound-protein analysis workflows

### 3. ML Model Components

#### Current Parallel Implementations
- `processors/structure/ml/models/` - Core ML models
- `models/compound/ml/` - Compound-specific models
- `pipeline/ml/core/` - Pipeline models
- Multiple model implementations across predictors

#### Integration Strategy
1. Create centralized model framework:
```
binding_data_processor/
  └── processors/
      └── structure/
          └── ml/
              └── models/
                  ├── base/          # Base model interfaces
                  ├── architectures/ # Neural network architectures
                  ├── ensemble/      # Ensemble methods
                  └── specialized/   # Domain-specific models
```

2. Standardize model interfaces:
- Define common model lifecycle methods
- Implement shared training utilities
- Standardize serialization formats
- Unify evaluation metrics

3. Enhance model capabilities:
- Centralize feature engineering
- Implement model registry system
- Add model versioning support
- Enable model composition

### 4. Infrastructure Components

#### Current Parallel Implementations
- `processors/structure/ml/activity/binding/structure/cache.py`
- `processors/structure/ml/activity/binding/structure/monitoring.py`
- `infrastructure/cache/`
- `infrastructure/monitoring/`

#### Integration Strategy
1. Consolidate infrastructure services:
```
binding_data_processor/
  └── infrastructure/
      ├── cache/         # Unified caching system
      ├── monitoring/    # Centralized monitoring
      ├── validation/    # Shared validation
      └── utils/         # Common utilities
```

2. Standardize service interfaces:
- Create consistent API for infrastructure services
- Implement proper dependency injection
- Add configuration management
- Enable service discovery

## Implementation Phases

### Phase 1: Core Framework (Weeks 1-2)

#### 1. Unified Analysis Framework
```python
class BaseAnalyzer(ABC):
    """Abstract base class for all analyzers."""
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()
        self.logger = logging.getLogger(self.__class__.__name__)
        self._cache = {}
        
    @abstractmethod
    def analyze(self, data: Any) -> Dict[str, Any]:
        """Core analysis method."""
        pass
        
    def _validate_input(self, data: Any) -> bool:
        """Validate input data."""
        pass
        
    def _format_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Format analysis results."""
        pass
```

#### 2. Enhanced Configuration Management
```python
@dataclass
class AnalysisConfig:
    """Base configuration for analysis."""
    cache_results: bool = True
    validate_input: bool = True
    include_metadata: bool = True
    log_level: int = logging.INFO

@dataclass 
class DynamicsConfig(AnalysisConfig):
    """Configuration for dynamics analysis."""
    use_normal_modes: bool = True
    use_correlations: bool = True
    use_domain_motions: bool = True
    n_modes: int = 10
    contact_cutoff: float = 8.0
```

#### 3. Integrated Caching System
```python
class AnalysisCache:
    """Multi-level caching for analysis results."""
    
    def __init__(self):
        self.memory_cache = {}
        self.file_cache = FileCache()
        
    async def get(self, key: str, domain: str) -> Optional[Dict[str, Any]]:
        """Get cached results with fallback."""
        # Try memory first
        if result := self.memory_cache.get(key):
            return result
            
        # Try file cache
        if result := await self.file_cache.get(key):
            self.memory_cache[key] = result
            return result
            
        return None
```

#### 4. Result Validation Framework
```python
class AnalysisValidator:
    """Validate analysis results."""
    
    def validate_structure(self, structure: Structure) -> bool:
        """Validate protein structure."""
        if not isinstance(structure, Structure):
            return False
        if not structure.get_list():
            return False
        return True
        
    def validate_results(self, results: Dict[str, Any]) -> bool:
        """Validate analysis results."""
        required = {"basic", "detailed", "metadata"}
        return all(k in results for k in required)
```

### Phase 2: Domain Integration (Weeks 3-4)

#### 1. Comprehensive Dynamics Analysis
```python
@dataclass
class DynamicsConfig:
    """Configuration for dynamics analysis."""
    use_normal_modes: bool = True
    use_correlations: bool = True
    use_domain_motions: bool = True
    use_contacts: bool = True
    use_alphafold_confidence: bool = True
    n_modes: int = 10
    contact_cutoff: float = 8.0
    min_plddt: float = 70.0

class DynamicsAnalyzer(BaseAnalyzer):
    """Enhanced protein dynamics analyzer."""
    
    def analyze(self, structure: Structure) -> Dict[str, Any]:
        """Analyze protein dynamics."""
        try:
            # Get structure coordinates
            coords = self.get_coordinates()
            
            # Calculate normal modes
            modes = None
            if self.config.use_normal_modes:
                modes = self.analyze_normal_modes(
                    n_modes=self.config.n_modes
                )
                
            # Calculate flexibility profile
            flexibility = self.analyze_flexibility()
            
            # Calculate correlations
            correlations = None
            if self.config.use_correlations and modes is not None:
                correlations = self._get_correlation_matrix()
                
            # Analyze domain motions
            domain_motions = None
            if self.config.use_domain_motions:
                domain_motions = self.analyze_domain_motions(
                    modes=modes,
                    include_correlations=True
                )
                
            # Weight by AlphaFold confidence
            if self.config.use_alphafold_confidence:
                flexibility = self.apply_confidence_weights(
                    flexibility,
                    min_plddt=self.config.min_plddt
                )
                
            results = {
                "flexibility": flexibility,
                "correlations": correlations,
                "modes": modes,
                "domain_motions": domain_motions
            }
            
            # Add contact analysis
            if self.config.use_contacts:
                results["contacts"] = self.analyze_contacts(
                    cutoff=self.config.contact_cutoff
                )
                
            return results
            
        except Exception as e:
            self.logger.error(f"Analysis failed: {str(e)}")
            return {}
            
    def analyze_site_dynamics(
        self,
        site_residues: List[int],
        include_context: bool = True
    ) -> Dict[str, Any]:
        """Analyze binding site dynamics."""
        try:
            # Get site-specific dynamics
            site_dynamics = self._analyze_site_specific(site_residues)
            
            # Get flexibility profile
            flexibility = self.analyze_flexibility()
            site_flexibility = {
                res: flexibility["scores"].get(res, 0.0)
                for res in site_residues
            }
            
            # Get contact network
            contacts = self.analyze_contacts()
            site_contacts = self._get_site_contacts(
                site_residues,
                contacts["network"]
            )
            
            # Calculate metrics
            metrics = {
                "avg_flexibility": float(np.mean(list(site_flexibility.values()))),
                "contact_density": site_contacts["density"],
                "surface_exposure": site_contacts["exposure"]
            }
            
            results = {
                "dynamics": site_dynamics,
                "flexibility": site_flexibility,
                "contacts": site_contacts,
                "metrics": metrics
            }
            
            # Add domain context
            if include_context:
                results["context"] = self._analyze_domain_context(site_residues)
                
            return results
            
        except Exception as e:
            self.logger.error(f"Site analysis failed: {str(e)}")
            return {}
```
            
    def analyze_domain_motions(
        self,
        structure: Structure,
        domains: Dict[str, List[int]],
        include_correlations: bool = True,
    ) -> Dict[str, Any]:
        """Analyze domain motions using PCA."""
        # Get CA coordinates
        coords = self._get_ca_coordinates(structure)
        
        # Calculate domain centers and motions
        centers = {}
        motions = {}
        for domain_id, residues in domains.items():
            # Get domain coordinates
            domain_coords = coords[residues]
            
            # Calculate center
            center = np.mean(domain_coords, axis=0)
            centers[domain_id] = center
            
            # Calculate principal components
            U, S, Vt = np.linalg.svd(
                domain_coords - center
            )
            
            motions[domain_id] = {
                "axes": Vt,
                "variance": S**2,
                "center": center,
            }
            
        # Calculate correlations if requested
        correlations = None
        if include_correlations:
            correlations = self._calculate_domain_correlations(
                coords,
                domains,
                centers,
            )
            
        return {
            "centers": centers,
            "motions": motions,
            "correlations": correlations,
        }
```

#### 2. Toxicity Integration with Dynamics
```python
class ToxicityAnalyzer(BaseAnalyzer):
    """Enhanced toxicity analyzer with dynamics."""
    
    def __init__(self):
        super().__init__()
        self.dynamics = DynamicsAnalyzer()
        self.bbb = BBBPredictorEnhanced()
        
    def analyze_binding_site(
        self,
        structure: Structure,
        site_residues: List[int],
    ) -> Dict[str, Any]:
        """Analyze binding site with dynamics."""
        # Get dynamics analysis
        dynamics = self.dynamics.analyze(structure)
        
        # Get domain context
        domain_context = self._get_domain_context(
            site_residues,
            dynamics["domains"],
        )
        
        # Calculate site flexibility
        site_flexibility = self._calculate_site_flexibility(
            site_residues,
            dynamics["flexibility"],
            domain_context,
        )
        
        # Assess binding site accessibility
        accessibility = self._assess_accessibility(
            site_residues,
            dynamics["domains"]["interfaces"],
            site_flexibility,
        )
        
        # Calculate risk factors
        risk_factors = self._calculate_risk_factors(
            site_flexibility,
            accessibility,
            domain_context,
        )
        
        return {
            "dynamics": {
                "flexibility": site_flexibility,
                "domain_context": domain_context,
                "accessibility": accessibility,
            },
            "risk_factors": risk_factors,
        }
        
    def _calculate_risk_factors(
        self,
        flexibility: Dict[int, float],
        accessibility: float,
        domain_context: Dict[str, Any],
    ) -> Dict[str, float]:
        """Calculate toxicity risk factors from dynamics."""
        risks = {}
        
        # Assess flexibility risk
        avg_flex = np.mean(list(flexibility.values()))
        risks["flexibility"] = self._score_flexibility_risk(avg_flex)
        
        # Assess accessibility risk
        risks["accessibility"] = self._score_accessibility_risk(
            accessibility
        )
        
        # Assess domain motion risk
        if "motions" in domain_context:
            motion_amplitude = np.mean([
                m["variance"][0]  # Principal component
                for m in domain_context["motions"].values()
            ])
            risks["domain_motion"] = self._score_motion_risk(
                motion_amplitude
            )
            
        # Assess interface risk
        if "interfaces" in domain_context:
            interface_size = len(domain_context["interfaces"])
            risks["interface"] = self._score_interface_risk(
                interface_size
            )
            
        return risks
```

#### 2. Toxicity Integration
```python
class ToxicityAnalyzer(BaseAnalyzer):
    """Enhanced toxicity analyzer with dynamics."""
    
    def __init__(self):
        super().__init__()
        self.dynamics = DynamicsAnalyzer()
        self.bbb = BBBPredictorEnhanced()
        
    def analyze_binding_site(
        self,
        structure: Structure,
        site_residues: List[int]
    ) -> Dict[str, Any]:
        """Analyze binding site with dynamics."""
        # Get site dynamics
        dynamics = self.dynamics.analyze_site_dynamics(
            structure,
            site_residues
        )
        
        # Adjust predictions based on dynamics
        flexibility = dynamics["flexibility"]
        accessibility = self._calculate_accessibility(
            dynamics["contacts"],
            dynamics["surface"]
        )
        
        return {
            "dynamics": dynamics,
            "accessibility": accessibility,
            "risk_factors": self._assess_risks(
                flexibility,
                accessibility
            )
        }
```

#### 3. Infrastructure Integration
```python
class AnalysisPipeline:
    """Integrated analysis pipeline."""
    
    def __init__(self):
        self.cache = AnalysisCache()
        self.validator = AnalysisValidator()
        self.dynamics = DynamicsAnalyzer()
        self.toxicity = ToxicityAnalyzer()
        
    async def analyze(
        self,
        structure: Structure,
        analysis_type: str
    ) -> Dict[str, Any]:
        """Run integrated analysis."""
        # Validate input
        if not self.validator.validate_structure(structure):
            raise ValueError("Invalid structure")
            
        # Check cache
        cache_key = self._get_cache_key(structure, analysis_type)
        if cached := await self.cache.get(cache_key):
            return cached
            
        # Run analysis
        if analysis_type == "dynamics":
            results = self.dynamics.analyze(structure)
        elif analysis_type == "toxicity":
            results = self.toxicity.analyze(structure)
        else:
            raise ValueError(f"Unknown analysis type: {analysis_type}")
            
        # Validate and cache results
        if self.validator.validate_results(results):
            await self.cache.set(cache_key, results)
            return results
        else:
            raise ValueError("Invalid analysis results")
```

### Phase 3: Testing & Validation (Weeks 5-6)

#### 1. Unit Tests
```python
def test_dynamics_analyzer():
    """Test dynamics analysis."""
    structure = load_test_structure()
    analyzer = DynamicsAnalyzer()
    
    results = analyzer.analyze(structure)
    assert "flexibility" in results
    assert "modes" in results
    assert results["metadata"]["timestamp"]
    
def test_integrated_pipeline():
    """Test analysis pipeline."""
    structure = load_test_structure()
    pipeline = AnalysisPipeline()
    
    results = pipeline.analyze(structure, "dynamics")
    assert results["flexibility"]["mean_b_factor"] > 0
    assert len(results["modes"]) == 10
```

#### 2. Integration Tests
```python
def test_toxicity_with_dynamics():
    """Test toxicity prediction with dynamics."""
    structure = load_test_structure()
    analyzer = ToxicityAnalyzer()
    
    site_residues = [10, 11, 12, 13]
    results = analyzer.analyze_binding_site(
        structure,
        site_residues
    )
    
    assert "dynamics" in results
    assert "accessibility" in results
    assert "risk_factors" in results
```

### Phase 4: Documentation & Deployment (Weeks 7-8)

1. API Documentation
```python
class DynamicsAnalyzer:
    """Analyze protein dynamics using multiple methods.
    
    This analyzer provides:
    - Flexibility analysis using B-factors
    - Normal mode analysis using elastic network model
    - Domain motion analysis
    - Contact network analysis
    - Site-specific dynamics
    
    The analysis can be configured via DynamicsConfig to enable/disable
    specific analysis types and adjust parameters.
    
    Example:
        analyzer = DynamicsAnalyzer()
        results = analyzer.analyze(structure)
        flexibility = results["flexibility"]
        modes = results["modes"]
    """
```

2. Usage Examples
```python
# Basic dynamics analysis
analyzer = DynamicsAnalyzer()
results = analyzer.analyze(structure)

# Site-specific analysis
site_results = analyzer.analyze_site_dynamics(
    structure,
    site_residues=[10, 11, 12, 13]
)

# Integrated toxicity prediction
pipeline = AnalysisPipeline()
results = pipeline.analyze(structure, "toxicity")
```

3. Performance Optimization
```python
class DynamicsCache:
    """Optimized caching for dynamics analysis."""
    
    def __init__(self):
        self.normal_modes = {}
        self.contact_maps = {}
        
    def cache_modes(self, structure_id: str, modes: np.ndarray):
        """Cache normal modes for reuse."""
        self.normal_modes[structure_id] = modes
        
    def cache_contacts(self, structure_id: str, contacts: np.ndarray):
        """Cache contact maps for reuse."""
        self.contact_maps[structure_id] = contacts
```

### Success Metrics

1. Code Quality
- Zero circular imports
- >90% test coverage
- Comprehensive documentation
- Type hints and validation

2. Performance
- Cache hit rate >80%
- Analysis time <30s
- Memory usage <2GB
- Efficient data structures

3. Integration
- Clean interfaces
- Proper error handling
- Consistent patterns
- Good separation of concerns

4. Maintainability
- Clear documentation
- Modular design
- Easy configuration
- Simple deployment

### Phase 2: Migration
1. Move psychopharm predictors to new structure
2. Consolidate protein analysis modules
3. Migrate ML models to central framework
4. Update infrastructure usage

### Phase 3: Enhancement
1. Implement advanced model capabilities
2. Add comprehensive monitoring
3. Enhance integration testing
4. Optimize performance

### Phase 4: Cleanup
1. Remove deprecated implementations
2. Update documentation
3. Verify backwards compatibility
4. Validate consolidated architecture

## Benefits

1. **Reduced Duplication**
- Consolidated code base
- Unified interfaces
- Shared utilities
- Centralized configuration

2. **Enhanced Maintainability**
- Clear architecture
- Standardized patterns
- Better documentation
- Easier testing

3. **Improved Functionality**
- Better integration
- Enhanced features
- Optimized performance
- Easier extensions

4. **Better Developer Experience**
- Clear structure
- Consistent patterns
- Better tooling
- Easier onboarding

## Risks and Mitigation

1. **Breaking Changes**
- Maintain compatibility layers
- Version major changes
- Provide migration guides
- Update gradually

2. **Performance Impact**
- Profile consolidated code
- Optimize critical paths
- Cache effectively
- Monitor metrics

3. **Integration Issues**
- Comprehensive testing
- Staged rollout
- Feature flags
- Rollback plans

## Success Metrics

1. **Code Quality**
- Reduced code duplication
- Improved test coverage
- Better documentation
- Cleaner architecture

2. **Performance**
- Response times
- Resource usage
- Cache effectiveness
- Training speed

3. **Developer Productivity**
- Development velocity
- Bug resolution time
- Code review efficiency
- Onboarding time

## Next Steps

1. Create detailed implementation plan for each phase
2. Set up metrics collection for baseline
3. Begin with core framework implementation
4. Create proof-of-concept for key integrations

# Architecture Consolidation Plan

## Overview

This document outlines a comprehensive plan to consolidate and enhance the architecture of the binding_data_processor system, addressing current issues while preserving and optimizing valuable functionality.

## Current State Analysis

After analyzing all surface area and volume calculation implementations across the codebase, here is the consolidated architecture:

Solvent Accessible Surface Area (SASA) calculations:

Consolidated in surface/base.py and surface.py using calculate_sasa() from utils.py
Purpose: Calculate protein solvent accessibility
Convex Hull calculations (kept separate):

binding/structure/surface.py - For protein shape analysis
pockets.py - For binding pocket volume calculation
geometry.py - For general geometric analysis (volume, surface area, shape descriptors)
Each serves a distinct purpose and should remain separate
Small Molecule calculations:

properties/base.py using RDKit's ComputeMolSurf
Purpose: Calculate small molecule surface area
Kept separate as it uses specialized RDKit functionality
The consolidation improves the codebase by:

Standardizing SASA calculations across protein solvent accessibility code
Maintaining specialized implementations where appropriate
Preserving clear separation of concerns between different types of surface/volume analysis
Keeping efficient, purpose-built implementations for each use case

Looking at the code organization, I see that:

The code is currently split across:
prediction.py (main file)
prediction/prediction.py (implementation details)
prediction/init.py (exports)
The functionality can be better organized into:
Core prediction interface
Structure prediction implementations (AlphaFold, local model)
Property prediction
Analysis integration
Metrics calculation
A better structure would be:
protein/
  prediction/
    __init__.py           # Exports main interface
    base.py              # Abstract base classes
    alphafold/
      __init__.py        # AlphaFold implementation
      api.py             # API client
      local.py           # Local model
    local/
      __init__.py        # Local model implementation  
      features.py        # Feature extraction
      model.py           # Model definition
    properties/
      __init__.py        # Property prediction
      stability.py       # Stability prediction
      flexibility.py     # Flexibility prediction
      function.py        # Function prediction
    metrics/
      __init__.py        # Metrics calculation
      confidence.py      # Confidence scores
      validation.py      # Validation metrics




### 1. Key Insights from Implementation

#### Configuration Centralization Success
- Moving ProteinAnalysisConfig to core/config.py resolved circular imports
- Demonstrates effectiveness of centralized configuration pattern
- Provides template for handling other config classes
- Enables better dependency management
- Facilitates configuration validation and type checking

#### Infrastructure Patterns
- Common infrastructure needs across web, pipeline, and processors
- Centralization reduces duplication and standardizes interfaces
- Core infrastructure should be independent of specific domains
- Shared utilities improve consistency and maintainability
- Common patterns emerge for error handling, logging, and monitoring

#### Validation Framework Opportunities
- Multiple validation approaches can be unified
- Common validation patterns emerge across modules
- Centralized validation improves consistency
- Shared schemas reduce duplication
- Standardized error handling improves debugging

### 2. Implementation Priorities

1. Core Configuration
- Move all config classes to core/config.py
- Use dataclasses for type safety
- Add validation methods
- Implement config inheritance where appropriate
- Add comprehensive documentation

2. Infrastructure Layer
- Create unified base classes
- Implement consistent interfaces
- Add proper error handling
- Standardize logging patterns
- Centralize monitoring

3. Validation Framework
- Create core validation utilities
- Define standard schemas
- Implement reusable validators
- Add comprehensive testing
- Document validation patterns

### 3. Specific Import Chain Resolution

Current import chain causing errors in export_compounds.py:
```python
scripts/export_compounds.py
└── binding_data_processor.data_sources.bindingdb_enhanced
    └── binding_data_processor.pipeline.base
        └── binding_data_processor.processors.psychopharm
            └── binding_data_processor.processors.structure
                └── binding_data_processor.processors.structure.ml.activity.protein.analysis.config
```

Resolution Steps:
1. Move ProteinAnalysisConfig to core/config.py (Completed)
2. Update import paths to use absolute imports:
```python
# In export_compounds.py
from binding_data_processor.core.config import ProteinAnalysisConfig
from binding_data_processor.data_sources.bindingdb_enhanced import BindingDBSourceEnhanced
```

3. Break circular dependencies:
- Move shared interfaces to core module
- Use dependency injection for configuration
- Implement proper layering:
  * Core (configs, interfaces)
  * Domain (analysis, processing)
  * Infrastructure (pipeline, web)
  * Application (scripts, CLI)


# Architecture Consolidation Plan

## Overview

This document outlines a comprehensive plan to consolidate and enhance the architecture of the binding_data_processor system, addressing current issues while preserving and optimizing valuable functionality.

## Current State Analysis

### 1. Immediate Issues

#### Import Chain Problem
Current import chain causing errors:
```python
# binding/base.py
from ..protein import ProteinStructureAnalyzer

# protein/__init__.py
from .analysis import ProteinStructureAnalyzer

# protein/analysis/base.py
from .config import ProteinAnalysisConfig  # Fails here
```

Solution:
1. Move ProteinAnalysisConfig to a separate core config module
2. Update imports to use absolute paths
3. Break circular dependencies

### 2. Core Infrastructure Issues

#### Cache Management
- Duplicate caching logic between pipeline and web components
- Inconsistent cache interfaces
- Different implementation patterns
- Lack of proper cache invalidation

#### Monitoring and Metrics
- Scattered monitoring logic
- Different metric collection approaches
- Inconsistent reporting
- Mixed monitoring implementations

#### Logging
- Inconsistent logging patterns
- Different log levels and formats
- Duplicate logging code
- Lack of structured logging

#### Error Handling
- Inconsistent error handling
- Duplicate error types
- Mixed error hierarchies
- Unclear error recovery paths

#### Type System
- Duplicate type definitions
- Inconsistent type usage
- Mixed type hierarchies
- Lack of proper validation

### 2. Module-Specific Issues

#### Protein Structure Analysis
```
binding_data_processor/processors/structure/ml/activity/
├── protein/
│ ├── analysis/
│ │ ├── __init__.py
│ │ ├── base.py (imports ProteinAnalysisConfig)
│ │ ├── config.py (missing ProteinAnalysisConfig)
│ │ └── surface.py
│ ├── __init__.py
│ └── prediction.py
└── binding/
└── base.py (imports ProteinStructureAnalyzer)
```

Current Issues:
- Circular imports between analysis and config modules
- Missing configuration class definitions
- Unclear separation between analysis and prediction
- Mixed concerns in protein analysis modules

#### Domain Logic
- Duplicate functionality across modules
- Inconsistent interfaces and patterns
- Scattered infrastructure code
- Complex dependencies
- Mixed concerns

## Proposed Architecture

### 1. Core Infrastructure Layer

```
infrastructure/
├── cache/           # Unified caching system
│   ├── base.py
│   ├── memory.py
│   ├── file.py
│   └── redis.py
├── monitoring/      # System monitoring and metrics
│   ├── base.py
│   ├── metrics.py
│   └── alerts.py
├── logging/         # Centralized logging
│   ├── base.py
│   └── formatters.py
├── validation/      # Data validation framework
│   ├── base.py
│   └── schemas.py
├── errors/          # Error handling
│   ├── base.py
│   └── handlers.py
├── config/          # Configuration management
│   ├── base.py
│   └── validation.py
├── types/           # Core type definitions
│   ├── base.py
│   └── validators.py
└── testing/         # Test utilities and fixtures
    ├── base.py
    └── fixtures.py
```

### 2. Domain Layer

#### Analysis Module
```
analysis/
├── structure/       # Structure analysis
├── binding/         # Binding site analysis
├── activity/        # Activity prediction
├── properties/      # Property calculation
└── validation/      # Analysis validation
```

#### Machine Learning Module
```
ml/
├── models/          # Model definitions
├── training/        # Training pipelines
├── prediction/      # Prediction pipelines
├── evaluation/      # Model evaluation
└── optimization/    # Model optimization
```

#### Processing Module
```
processing/
├── pipeline/        # Pipeline framework
├── transforms/      # Data transformations
├── validation/      # Process validation
└── enrichment/      # Data enrichment
```

### 3. Integration Layer

#### Clients Module
```
clients/
├── base/           # Base client interfaces
├── web/            # Web service clients
├── database/       # Database clients
└── external/       # External API clients
```

#### Web Module
```
web/
├── api/            # API endpoints
├── components/     # UI components
├── templates/      # Page templates
└── static/         # Static assets
```

## Implementation Plan

### Phase 1: Fix Import Chain (Week 1)

1. Create Core Config Module
```python
# binding_data_processor/core/config.py
from dataclasses import dataclass, field
from typing import Dict, Optional
import torch

@dataclass
class ProteinAnalysisConfig:
    """Configuration for protein structure analysis."""
    # Model configuration
    model_path: Optional[str] = None
    device: str = "cuda" if torch.cuda.is_available() else "cpu"
    cache_dir: Optional[str] = None
    
    # Analysis settings
    analyze_pockets: bool = True
    analyze_dynamics: bool = True
    analyze_conservation: bool = True
    analyze_interfaces: bool = True
    analyze_quality: bool = True
    
    # Performance
    cache_results: bool = True
    n_jobs: int = 1
```

2. Update Import Structure
```python
# binding_data_processor/processors/structure/ml/activity/protein/analysis/base.py
from binding_data_processor.core.config import ProteinAnalysisConfig

class ProteinStructureAnalyzer:
    def __init__(self, config: ProteinAnalysisConfig):
        self.config = config
```

3. Fix Circular Dependencies
```python
# binding_data_processor/processors/structure/ml/activity/protein/__init__.py
from .analysis.base import ProteinStructureAnalyzer
from binding_data_processor.core.config import ProteinAnalysisConfig

__all__ = ['ProteinStructureAnalyzer', 'ProteinAnalysisConfig']
```

### Phase 2: Core Infrastructure (Weeks 2-3)

1. Configuration Management
```python
@dataclass
class ProteinAnalysisConfig:
    """Configuration for protein structure analysis."""
    # Analysis settings
    analyze_pockets: bool = True
    analyze_dynamics: bool = True
    analyze_conservation: bool = True
    analyze_interfaces: bool = True
    analyze_quality: bool = True

    # Thresholds
    pocket_min_volume: float = 100.0
    interface_cutoff: float = 5.0
    clash_cutoff: float = 0.4

    # Performance
    cache_results: bool = True
    n_jobs: int = 1
```

2. Error Handling
```python
class StructureAnalyzer:
    """Base class for structure analysis."""
    def __init__(self, config: ProteinAnalysisConfig):
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)

    def analyze(self, structure: Structure) -> Dict[str, Any]:
        """Template method for structure analysis."""
        try:
            self._validate_structure(structure)
            properties = self._analyze_structure(structure)
            self._validate_results(properties)
            return properties
        except Exception as e:
            self.logger.error(f"Analysis failed: {str(e)}")
            raise
```

3. Validation
```python
def _validate_structure(self, structure: Structure) -> None:
    """Validate input structure."""
    if not isinstance(structure, Structure):
        raise ValueError("Invalid structure type")
    if not structure.get_list():
        raise ValueError("Empty structure")

def _validate_results(self, results: Dict[str, Any]) -> None:
    """Validate analysis results."""
    required_keys = {"basic", "surface", "quality"}
    if not all(k in results for k in required_keys):
        raise ValueError("Missing required result properties")
```

### Phase 2: Domain Logic Migration (Weeks 3-4)

1. Create Base Classes
```python
class BaseStructureAnalyzer:
    """Base implementation of structure analyzer."""
    def __init__(self, config: ProteinAnalysisConfig):
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)

    def analyze(self, structure: Structure) -> Dict[str, Any]:
        """Template method for structure analysis."""
        try:
            self._validate_structure(structure)
            properties = self._analyze_structure(structure)
            self._validate_results(properties)
            return properties
        except Exception as e:
            self.logger.error(f"Analysis failed: {str(e)}")
            raise
```

2. Implement Interfaces
```python
class StructureAnalyzer(Protocol):
    def analyze(self, structure: Structure) -> Dict[str, Any]: ...
    def get_properties(self, structure: Structure) -> Dict[str, Any]: ...
```

### Phase 3: Integration and Testing (Weeks 5-6)

1. Unit Tests
```python
def test_protein_analysis():
    config = ProteinAnalysisConfig()
    analyzer = ProteinStructureAnalyzer(config)
    structure = load_test_structure()
    results = analyzer.analyze(structure)
    assert "basic" in results
    assert "surface" in results
    assert results["quality"]["clash_score"] < config.clash_cutoff
```

2. Integration Tests
```python
def test_binding_prediction():
    structure = load_test_structure()
    ligand = load_test_ligand()
    predictor = BindingPredictor()
    result = predictor.predict_binding(structure, ligand)
    assert result.score > 0.5
    assert len(result.interactions) > 0
```

## Success Metrics

### 1. Code Quality
- Zero circular imports
- Clear module boundaries
- Consistent error handling
- >90% test coverage
- Comprehensive documentation

### 2. Performance
- Reduced memory usage
- Faster processing times
- Efficient caching
- Optimized algorithms
- Better resource utilization

### 3. Maintainability
- Clear separation of concerns
- Consistent patterns
- Well-documented interfaces
- Easy to extend
- Reusable components

### 4. Developer Experience
- Faster onboarding
- Better tooling
- Clear documentation
- Easier debugging
- Streamlined workflows

## Dependencies

### 1. Core Dependencies
- Python 3.7+
- BioPython
- NumPy
- RDKit
- PyTorch (optional)

### 2. Development Tools
- pytest
- mypy
- black
- isort
- flake8

### 3. Infrastructure
- Redis (optional for caching)
- PostgreSQL
- Docker
- Kubernetes (optional)

## Risk Management

### 1. Breaking Changes
- Careful interface design
- Comprehensive testing
- Phased rollout
- Clear documentation
- Migration guides
- Rollback plans

### 2. Performance Impact
- Benchmark critical paths
- Profile hot spots
- Monitor metrics
- Optimize bottlenecks
- Performance testing
- Load testing

### 3. Integration Issues
- Clear interfaces
- Good error handling
- Proper validation
- Integration testing
- Monitoring
- Alerting

## Next Steps

1. Infrastructure Setup
- Set up new directory structure
- Create base classes
- Implement core interfaces
- Add initial tests

2. Migration Planning
- Identify critical paths
- Plan phased rollout
- Create backup strategy
- Document procedures

3. Development
- Start with core infrastructure
- Move to domain logic
- Add integration layer
- Implement testing

4. Validation
- Run comprehensive tests
- Validate performance
- Check compatibility
- Verify functionality

5. Deployment
- Stage changes
- Monitor metrics
- Gather feedback
- Make adjustments

## Conclusion

This consolidation plan provides a comprehensive approach to:
1. Fix current architectural issues
2. Improve code organization
3. Enhance performance
4. Increase maintainability
5. Better developer experience

The phased implementation ensures minimal disruption while achieving significant improvements in code quality, performance, and maintainability.
