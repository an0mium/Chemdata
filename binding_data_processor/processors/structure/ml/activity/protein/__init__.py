"""Protein structure analysis and prediction module."""

from binding_data_processor.core.config import ProteinAnalysisConfig
from .prediction import ProteinPredictor
from .analysis import ProteinStructureAnalyzer
from .analysis.dynamics import analyze_dynamics, analyze_flexibility, analyze_normal_modes, analyze_contacts
from .analysis.surface import analyze_surface

__all__ = [
    "ProteinAnalysisConfig",
    "ProteinPredictor",
    "ProteinStructureAnalyzer",
    "analyze_dynamics",
    "analyze_flexibility",
    "analyze_normal_modes",
    "analyze_contacts",
    "analyze_surface",
]
