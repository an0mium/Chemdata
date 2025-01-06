"""Compound analysis package.

This package provides comprehensive analysis capabilities for compound data:

Core Analysis Classes:
- CompoundAnalysis: Main class combining all analysis capabilities
- BindingAnalysisMixin: Analysis of binding properties and patterns
- ActivityAnalysisMixin: Analysis of activity and mechanisms
- SafetyAnalysisMixin: Analysis of safety risks and alerts
- PropertyAnalysisMixin: Analysis of physicochemical properties
- SARAnalysisMixin: Analysis of structure-activity relationships

Usage:
    from binding_data_processor.models.compound.analysis import CompoundAnalysis

    # Create analysis object
    analysis = CompoundAnalysis(...)

    # Run all analyses
    results = analysis.analyze_all()

    # Get summary
    summary = analysis.get_analysis_summary()

    # Get key findings
    findings = analysis.get_key_findings()

    # Get optimization suggestions
    suggestions = analysis.get_optimization_suggestions()
"""

from .base import AnalyzedCompound
from .binding_analysis import BindingAnalyzer
from .activity_analysis import ActivityAnalyzer
from .safety_analysis import SafetyAnalyzer
from .property_analysis import PropertyAnalyzer
from .sar_analysis import SARAnalyzer

__all__ = [
    'AnalyzedCompound',
    'BindingAnalyzer',
    'ActivityAnalyzer',
    'SafetyAnalyzer',
    'PropertyAnalyzer',
    'SARAnalyzer',
]
