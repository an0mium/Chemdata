"""Constants module for structure analysis system.

This module provides:
1. System constants
2. Configuration values
3. Standard parameters
4. Default settings
"""

from enum import Enum, auto
from typing import Dict, List, Set, FrozenSet
from pathlib import Path

# System paths
DEFAULT_CACHE_DIR = Path("cache/structure")
DEFAULT_METRICS_DIR = Path("metrics/structure")
DEFAULT_LOG_DIR = Path("logs/structure")


# Component names
class ComponentName(str, Enum):
    """Component name constants."""

    ANALYZER = "analyzer"
    INTEGRATOR = "integrator"
    COORDINATOR = "coordinator"
    FACTORY = "factory"
    SYSTEM = "system"


# Status values
class Status(str, Enum):
    """Status constants."""

    HEALTHY = "healthy"
    DEGRADED = "degraded"
    FAILED = "failed"


# Operation types
class Operation(str, Enum):
    """Operation type constants."""

    ANALYZE = "analyze"
    INTEGRATE = "integrate"
    COORDINATE = "coordinate"
    CREATE = "create"
    PROCESS = "process"


# Structure types
class StructureType(str, Enum):
    """Structure type constants."""

    PROTEIN = "protein"
    LIGAND = "ligand"
    COMPLEX = "complex"


# Analysis parameters
ANALYSIS_PARAMS = {
    # Distance thresholds (Angstroms)
    "contact_distance": 4.5,
    "hydrogen_bond_distance": 3.5,
    "hydrophobic_distance": 4.0,
    "ionic_distance": 4.0,
    "metal_coordination_distance": 3.0,
    "pi_stacking_distance": 5.0,
    "water_bridge_distance": 3.5,
    # Angle thresholds (degrees)
    "hydrogen_bond_angle": 63.0,
    "pi_stacking_angle": 30.0,
    # Energy thresholds (kcal/mol)
    "interaction_energy": -0.5,
    "stability_energy": -1.0,
    # Surface parameters
    "probe_radius": 1.4,
    "surface_density": 3.0,
    # Pocket parameters
    "min_pocket_volume": 100.0,
    "max_pocket_depth": 20.0,
    # Dynamics parameters
    "simulation_time": 1000,  # ps
    "time_step": 2.0,  # fs
    "temperature": 300.0,  # K
    # Analysis settings
    "num_modes": 10,
    "cutoff_distance": 8.0,
    "contact_threshold": 0.5,
}

# Atom types
ATOM_TYPES: FrozenSet[str] = frozenset(
    [
        "C",
        "N",
        "O",
        "S",
        "P",
        "H",
        "F",
        "Cl",
        "Br",
        "I",
        "Na",
        "K",
        "Mg",
        "Ca",
        "Zn",
        "Fe",
    ]
)

# Residue types
RESIDUE_TYPES: FrozenSet[str] = frozenset(
    [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]
)

# Secondary structure types
SECONDARY_STRUCTURE: FrozenSet[str] = frozenset(
    [
        "helix",
        "sheet",
        "loop",
        "turn",
        "coil",
    ]
)

# Interaction types
INTERACTION_TYPES: FrozenSet[str] = frozenset(
    [
        "hydrogen_bond",
        "hydrophobic",
        "ionic",
        "metal_coordination",
        "pi_stacking",
        "water_bridge",
    ]
)

# Chemical properties
CHEMICAL_PROPERTIES = {
    # Atomic radii (Angstroms)
    "atomic_radii": {
        "H": 1.20,
        "C": 1.70,
        "N": 1.55,
        "O": 1.52,
        "F": 1.47,
        "P": 1.80,
        "S": 1.80,
        "Cl": 1.75,
        "Br": 1.85,
        "I": 1.98,
    },
    # VDW radii (Angstroms)
    "vdw_radii": {
        "H": 1.09,
        "C": 1.70,
        "N": 1.55,
        "O": 1.52,
        "F": 1.47,
        "P": 1.80,
        "S": 1.80,
        "Cl": 1.75,
        "Br": 1.85,
        "I": 1.98,
    },
    # Electronegativity (Pauling scale)
    "electronegativity": {
        "H": 2.20,
        "C": 2.55,
        "N": 3.04,
        "O": 3.44,
        "F": 3.98,
        "P": 2.19,
        "S": 2.58,
        "Cl": 3.16,
        "Br": 2.96,
        "I": 2.66,
    },
}

# Performance settings
PERFORMANCE_SETTINGS = {
    # Cache settings
    "cache_size": 1000,
    "cache_ttl": 3600,  # seconds
    # Batch settings
    "batch_size": 100,
    "max_batch_size": 1000,
    # Threading settings
    "num_threads": 4,
    "max_threads": 16,
    # Memory settings
    "memory_limit": 1024 * 1024 * 1024,  # 1 GB
    "min_free_memory": 100 * 1024 * 1024,  # 100 MB
}

# Validation settings
VALIDATION_SETTINGS = {
    # Structure validation
    "min_atoms": 10,
    "max_atoms": 100000,
    "min_residues": 1,
    "max_residues": 10000,
    "min_chains": 1,
    "max_chains": 100,
    # Quality thresholds
    "max_clash_score": 50.0,
    "min_rotamer_score": 0.8,
    "max_rama_outliers": 0.05,
    # Geometry limits
    "min_bond_length": 0.5,  # Angstroms
    "max_bond_length": 4.0,
    "min_bond_angle": 60.0,  # degrees
    "max_bond_angle": 180.0,
}


# File formats
class FileFormat(str, Enum):
    """File format constants."""

    PDB = "pdb"
    CIF = "cif"
    MOL2 = "mol2"
    SDF = "sdf"
    JSON = "json"
    CSV = "csv"


# Error codes
class ErrorCode(str, Enum):
    """Error code constants."""

    VALIDATION_ERROR = "validation_error"
    PROCESSING_ERROR = "processing_error"
    SYSTEM_ERROR = "system_error"
    IO_ERROR = "io_error"
    RESOURCE_ERROR = "resource_error"


# Default values
DEFAULT_VALUES = {
    # Analysis defaults
    "resolution": 1.0,  # Angstroms
    "probe_radius": 1.4,  # Angstroms
    "surface_density": 3.0,  # points/A^2
    # Performance defaults
    "cache_size": 1000,
    "batch_size": 100,
    "num_threads": 4,
    # Validation defaults
    "clash_threshold": 0.4,  # Angstroms
    "rama_z_score": -2.0,
    "max_b_factor": 100.0,
}

# Unit conversions
UNIT_CONVERSIONS = {
    # Distance
    "angstrom_to_nm": 0.1,
    "nm_to_angstrom": 10.0,
    # Energy
    "kcal_to_kj": 4.184,
    "kj_to_kcal": 0.239,
    # Time
    "ps_to_ns": 0.001,
    "ns_to_ps": 1000.0,
}


# Logging levels
class LogLevel(str, Enum):
    """Log level constants."""

    DEBUG = "debug"
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


# Monitoring settings
MONITORING_SETTINGS = {
    # Metric collection
    "collect_interval": 60,  # seconds
    "max_metrics": 10000,
    "metric_ttl": 86400,  # 24 hours
    # Health checks
    "health_check_interval": 300,  # 5 minutes
    "max_failures": 3,
    "recovery_time": 600,  # 10 minutes
    # Alerts
    "alert_cooldown": 3600,  # 1 hour
    "max_alerts": 100,
}

# System limits
SYSTEM_LIMITS = {
    # Resource limits
    "max_memory": 8 * 1024 * 1024 * 1024,  # 8 GB
    "max_cpu_percent": 90,
    "max_disk_percent": 95,
    # Operation limits
    "max_retries": 3,
    "timeout": 3600,  # 1 hour
    "max_batch_operations": 1000,
    # Rate limits
    "max_requests_per_second": 10,
    "max_concurrent_operations": 4,
}
