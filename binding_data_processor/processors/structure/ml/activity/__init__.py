"""Activity prediction and analysis module.

This module provides functionality for:
1. Binding affinity prediction
2. Target prediction 
3. Activity classification
4. Structure-activity relationships
5. Binding site analysis
6. Protein structure prediction and analysis
7. AlphaFold integration
8. Structure-based drug design
"""

from typing import Dict, List, Optional, Set, Tuple, Union, Any

from .base import ActivityPredictor, ActivityPredictorConfig
from .binding import (
    BindingPredictor,
    BindingSitePredictor,
    BindingModePredictor,
)
from .target import TargetPredictor
from .interaction import InteractionPredictor
from .protein import (
    ProteinPredictor,
    ProteinStructureAnalyzer,
)
from .protein.analysis.dynamics import (
    analyze_dynamics,
    analyze_flexibility,
    analyze_normal_modes,
    analyze_contacts,
)
from .alphafold import (
    AlphaFoldPredictor,
    AlphaFoldConfig,
    AlphaFoldResult,
)

__all__ = [
    # Base classes
    "ActivityPredictor",
    "ActivityPredictorConfig",
    # Binding prediction
    "BindingPredictor",
    "BindingSitePredictor",
    "BindingModePredictor",
    # Target prediction
    "TargetPredictor",
    # Interaction prediction
    "InteractionPredictor",
    # Protein prediction and analysis
    "ProteinPredictor",
    "ProteinStructureAnalyzer",
    # Protein dynamics analysis
    "analyze_dynamics",
    "analyze_flexibility",
    "analyze_normal_modes",
    "analyze_contacts",
    # AlphaFold integration
    "AlphaFoldPredictor",
    "AlphaFoldConfig",
    "AlphaFoldResult",
]
