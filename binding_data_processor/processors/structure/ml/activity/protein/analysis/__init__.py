"""Protein structure analysis package."""

# Core functionality
from .base import ProteinStructureAnalyzer
from binding_data_processor.core.config import ProteinAnalysisConfig

# Analysis modules
from .surface import analyze_surface
from .dynamics import analyze_dynamics, analyze_flexibility, analyze_normal_modes, analyze_contacts
from .quality import calculate_quality_metrics, analyze_clashes, analyze_rotamers
from .pockets import find_binding_pockets, analyze_pocket_properties
from .interfaces import analyze_interfaces
from .depth import calculate_residue_depth

# Utility functions
from .utils import (
    get_atom_radius,
    calculate_center_of_mass,
    calculate_radius_of_gyration,
    get_residue_property,
    calculate_surface_exposure,
    get_secondary_structure,
    calculate_interface_area,
    calculate_cavity_volume,
    get_atom_coords,
    get_backbone_atoms,
    get_sequence_from_structure,
    find_contacts,
    calculate_sasa,
    get_ca_atoms,
)

__all__ = [
    "ProteinStructureAnalyzer",
    "ProteinAnalysisConfig",
    "analyze_surface",
    "analyze_dynamics",
    "analyze_flexibility",
    "analyze_normal_modes",
    "analyze_contacts",
    "calculate_quality_metrics",
    "analyze_clashes",
    "analyze_rotamers",
    "find_binding_pockets",
    "analyze_pocket_properties",
    "analyze_interfaces",
    "get_atom_radius",
    "calculate_center_of_mass",
    "calculate_radius_of_gyration",
    "get_residue_property",
    "calculate_surface_exposure",
    "calculate_residue_depth",
    "get_secondary_structure",
    "calculate_interface_area",
    "calculate_cavity_volume",
    "get_atom_coords",
    "get_backbone_atoms",
    "get_sequence_from_structure",
    "find_contacts",
    "calculate_sasa",
    "get_ca_atoms",
]
