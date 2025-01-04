"""Analysis manager for pipeline.

This package provides analysis functionality for:
1. Binding analysis
2. Activity analysis
3. Safety analysis
4. SAR analysis
5. Property analysis

The main interface is the AnalysisManager class which coordinates
all analysis components.
"""

from .base import AnalysisManager, AnalysisConfig, AnalysisStats
from .binding import BindingAnalyzer
from .activity import ActivityAnalyzer
from .safety import SafetyAnalyzer
from .sar import SARAnalyzer
from .properties import PropertyAnalyzer

__all__ = [
    # Main interface
    "AnalysisManager",
    "AnalysisConfig",
    "AnalysisStats",
    
    # Component analyzers
    "BindingAnalyzer",
    "ActivityAnalyzer",
    "SafetyAnalyzer",
    "SARAnalyzer",
    "PropertyAnalyzer",
]

# Version info
__version__ = "1.0.0"
