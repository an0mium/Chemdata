"""Chemical compound data models.

This package provides comprehensive data models for representing chemical compounds
and their associated data, with support for:

Core Data:
- Basic identifiers (name, SMILES, InChI, CAS)
- Chemical properties (MW, LogP, TPSA, etc.)
- Database IDs (ChEMBL, PubChem, DrugBank, BindingDB)
- Common names and search results
- Legal status and scheduling

Binding & Activity:
- BindingDB data integration
- Target predictions and profiles
- Activity classifications
- SAR analysis
- Pharmacophore detection
- Receptor interactions

ML Predictions:
- Feature management and caching
- Model version tracking
- Binding affinity predictions
- Toxicity predictions
- Abuse potential assessment
- Activity predictions
- Blood-brain barrier penetration
- Psychoactive effects
- Nootropic activity

Web-Enriched Data:
- Patent information and analysis
- Literature references and citations
- Community experience reports
- Safety profiles and warnings
- Dosage information
- Administration routes
- Duration statistics
- Drug combinations
- Risk factors
- Regulatory status
- Clinical trial status
- Research status

Analysis Capabilities:
- Binding profile analysis
- Activity pattern detection
- Safety assessment
- SAR analysis
- Mechanism prediction
- Structure validation
- Property calculation
- Risk assessment

Export Features:
- Flexible TSV export
- Configurable JSON export
- Detailed report generation
- Custom column selection
- Data validation
- Format conversion

The models support both flexible data structures for dynamic content and
structured fields for core data to ensure consistency and enable efficient
filtering and analysis. The main interface is through a clean inheritance chain:

BaseCompound -> EnrichedCompound -> MLCompound -> AnalyzedCompound

For simpler use cases, the base classes can be used independently.
"""

import warnings

# Import new modular models
from .compound import (
    Compound,
    BaseCompound,
    EnrichedCompound,
    MLCompound,
    AnalyzedCompound,
    CompoundType,
    LegalStatus,
    PsychoactiveClass,
    NootropicMechanism,
    BBBPermeability,
    BindingType,
    ActivityType,
    RiskLevel,
    TargetData,
    ValidationError,
)

# Import mixins for custom implementations
from .enrichment import WebEnrichmentMixin
from .predictions import PredictionsMixin
from .analysis import AnalysisMixin

# Version reflecting combined functionality
__version__ = "2.1.0"

# Export all public classes and types
__all__ = [
    # Main interface
    'Compound',
    
    # Base classes
    'BaseCompound',
    'EnrichedCompound',
    'MLCompound',
    'AnalyzedCompound',
    
    # Types
    'CompoundType',
    'LegalStatus',
    'PsychoactiveClass',
    'NootropicMechanism',
    'BBBPermeability',
    'BindingType',
    'ActivityType',
    'RiskLevel',
    'TargetData',
    
    # Mixins
    'WebEnrichmentMixin',
    'PredictionsMixin',
    'AnalysisMixin',
    
    # Exceptions
    'ValidationError',
]

# Deprecation warnings for old files
def _warn_deprecated(old_name: str, new_name: str) -> None:
    warnings.warn(
        f"{old_name} is deprecated and will be removed in a future version. "
        f"Use {new_name} instead.",
        DeprecationWarning,
        stacklevel=2
    )

class _DeprecatedModule:
    def __init__(self, old_name: str, new_name: str):
        self._old_name = old_name
        self._new_name = new_name
        
    def __getattr__(self, name):
        _warn_deprecated(self._old_name, self._new_name)
        return getattr(Compound, name)

# Deprecated imports with warnings
compound = _DeprecatedModule('compound.py', 'compound/__init__.py')
compound_base = _DeprecatedModule('compound_base.py', 'compound/base.py')
compound_ml = _DeprecatedModule('compound_ml.py', 'compound/ml.py')
compound_enrichment = _DeprecatedModule('compound_enrichment.py', 'compound/enrichment.py')
compound_analysis = _DeprecatedModule('compound_analysis.py', 'compound/analysis.py')
compound_export = _DeprecatedModule('compound_export.py', 'compound/export.py')
