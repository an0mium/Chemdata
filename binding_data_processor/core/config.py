"""Core configuration classes for protein structure and dynamics analysis."""

from dataclasses import dataclass, field
from typing import Dict, Optional, Any
import torch


@dataclass
class ProteinAnalysisConfig:
    """Configuration for protein structure analysis.

    Attributes:
        # Model Configuration
        model_path: Path to model weights
        device: Device to run computations on
        cache_dir: Directory for caching results
        dssp_binary: Path to DSSP executable

        # AlphaFold Integration
        use_alphafold: Whether to use AlphaFold predictions
        alphafold_config: AlphaFold-specific configuration
        use_alphafold_api: Whether to use AlphaFold API
        alphafold_api_base_url: Base URL for AlphaFold API
        alphafold_api_timeout: Timeout for API calls

        # Analysis Configuration
        analyze_pockets: Whether to analyze binding pockets
        analyze_dynamics: Whether to analyze protein dynamics
        analyze_conservation: Whether to analyze sequence conservation
        analyze_interfaces: Whether to analyze protein interfaces
        analyze_quality: Whether to analyze structure quality

        # Feature Configuration
        use_pssm: Whether to use PSSM features
        use_dssp: Whether to use DSSP features
        use_surface_features: Whether to use surface features
        use_conservation: Whether to use conservation features
        use_dynamics: Whether to use dynamics features
        use_quality: Whether to use quality metrics
        use_surface: Whether to use surface analysis
        use_pockets: Whether to use pocket detection
        use_interfaces: Whether to use interface analysis

        # Performance Settings
        batch_size: Batch size for parallel processing
        num_workers: Number of worker processes
        n_jobs: Number of parallel jobs
        cache_results: Whether to cache results

        # Analysis Thresholds
        min_pocket_volume: Minimum volume for binding pockets
        pocket_max_depth: Maximum pocket depth
        cavity_min_depth: Minimum cavity depth
        interface_min_area: Minimum interface area
        interface_distance_cutoff: Distance cutoff for interfaces
        surface_probe_radius: Probe radius for surface calculations
        clash_distance_cutoff: Distance cutoff for clashes
        hbond_distance_cutoff: Distance cutoff for H-bonds
        hbond_angle_cutoff: Angle cutoff for H-bonds
        hydrophobic_contact_cutoff: Distance cutoff for hydrophobic contacts
        residue_exposure_threshold: SASA threshold for exposed residues
        contact_distance: Distance cutoff for residue contacts
        clash_overlap: Overlap distance for atomic clashes

        # Quality Thresholds
        quality_thresholds: Thresholds for structure quality metrics

        # Confidence Thresholds
        confidence_thresholds: Thresholds for confidence levels

        # Atom Properties
        atom_radii: Van der Waals radii for atoms

        # Residue Properties
        residue_properties: Physical properties for amino acids

        # Domain Analysis Configuration
        include_domain_analysis: Whether to include domain motion analysis
        domain_analysis_params: Parameters for domain analysis

        # Dynamics Analysis Configuration
        flexibility_cutoff: B-factor cutoff for flexibility analysis
        contact_cutoff: Distance cutoff for contact analysis
        n_modes: Number of normal modes to calculate
        logging_level: Logging level for analysis
    """

    # Model Configuration
    model_path: Optional[str] = None
    device: str = "cuda" if torch.cuda.is_available() else "cpu"
    cache_dir: Optional[str] = None
    dssp_binary: Optional[str] = None

    # AlphaFold Integration
    use_alphafold: bool = True
    alphafold_config: Optional["AlphaFoldConfig"] = None
    use_alphafold_api: bool = True
    alphafold_api_base_url: str = "https://alphafold.ebi.ac.uk/api"
    alphafold_api_timeout: int = 300
    use_alphafold_pae: bool = True  # Use predicted aligned error for dynamics
    use_alphafold_confidence: bool = True  # Use pLDDT confidence scores
    min_plddt_score: float = 70.0  # Minimum pLDDT score to consider
    pae_power: float = 1.0  # Power for PAE-based contact weighting

    # Analysis Configuration
    analyze_pockets: bool = True
    analyze_dynamics: bool = True
    analyze_conservation: bool = True
    analyze_interfaces: bool = True
    analyze_quality: bool = True

    # Feature Configuration
    use_pssm: bool = True
    use_dssp: bool = True
    use_surface_features: bool = True
    use_conservation: bool = True
    use_dynamics: bool = True
    use_quality: bool = True
    use_surface: bool = True
    use_pockets: bool = True
    use_interfaces: bool = True

    # Performance Settings
    batch_size: int = 100
    num_workers: int = 4
    n_jobs: int = 1
    cache_results: bool = True

    # Analysis Thresholds
    min_pocket_volume: float = 100.0
    pocket_max_depth: float = 20.0
    cavity_min_depth: float = 2.0
    interface_min_area: float = 100.0
    interface_distance_cutoff: float = 5.0
    surface_probe_radius: float = 1.4
    clash_distance_cutoff: float = 0.4
    hbond_distance_cutoff: float = 3.5
    hbond_angle_cutoff: float = 30.0
    hydrophobic_contact_cutoff: float = 4.5
    residue_exposure_threshold: float = 2.8
    contact_distance: float = 5.0
    clash_overlap: float = 0.4

    # Quality Thresholds
    quality_thresholds: Dict[str, float] = field(
        default_factory=lambda: {
            "clash_score": 20.0,
            "rama_favored": 0.98,
            "rotamer_outliers": 0.01,
            "cbeta_deviations": 0.01,
        }
    )

    # Confidence Thresholds
    confidence_thresholds: Dict[str, float] = field(
        default_factory=lambda: {
            "high": 90.0,
            "medium": 70.0,
            "low": 50.0,
        }
    )

    # Atom Properties
    atom_radii: Dict[str, float] = field(
        default_factory=lambda: {
            "C": 1.7,
            "N": 1.55,
            "O": 1.52,
            "S": 1.8,
            "P": 1.8,
            "H": 1.2,
            "F": 1.47,
            "Cl": 1.75,
            "Br": 1.85,
            "I": 1.98,
        }
    )

    # Residue Properties
    residue_properties: Dict[str, Dict[str, float]] = field(
        default_factory=lambda: {
            "hydrophobicity": {  # Kyte-Doolittle scale
                "ILE": 4.5,
                "VAL": 4.2,
                "LEU": 3.8,
                "PHE": 2.8,
                "CYS": 2.5,
                "MET": 1.9,
                "ALA": 1.8,
                "GLY": -0.4,
                "THR": -0.7,
                "SER": -0.8,
                "TRP": -0.9,
                "TYR": -1.3,
                "PRO": -1.6,
                "HIS": -3.2,
                "GLU": -3.5,
                "GLN": -3.5,
                "ASP": -3.5,
                "ASN": -3.5,
                "LYS": -3.9,
                "ARG": -4.5,
            },
            "charge": {  # Net charge at pH 7
                "ARG": 1,
                "LYS": 1,
                "ASP": -1,
                "GLU": -1,
                "HIS": 0.5,
            },
            "volume": {  # Residue volumes in Å³
                "ALA": 88.6,
                "ARG": 173.4,
                "ASN": 114.1,
                "ASP": 111.1,
                "CYS": 108.5,
                "GLN": 143.8,
                "GLU": 138.4,
                "GLY": 60.1,
                "HIS": 153.2,
                "ILE": 166.7,
                "LEU": 166.7,
                "LYS": 168.6,
                "MET": 162.9,
                "PHE": 189.9,
                "PRO": 112.7,
                "SER": 89.0,
                "THR": 116.1,
                "TRP": 227.8,
                "TYR": 193.6,
                "VAL": 140.0,
            },
        }
    )

    # Domain Analysis Configuration
    include_domain_analysis: bool = True
    domain_analysis_params: Dict[str, Any] = field(
        default_factory=lambda: {
            "use_structure_based": True,
            "use_sequence_based": False,
            "min_domain_size": 30,
            "contact_cutoff": 8.0,
            "interface_cutoff": 5.0,
            "hinge_score_cutoff": 0.5,
        }
    )

    # Dynamics Analysis Configuration
    flexibility_cutoff: float = 1.0
    contact_cutoff: float = 8.0
    n_modes: int = 10
    logging_level: str = "INFO"

    def __post_init__(self):
        """Validate configuration after initialization."""
        self._validate_config()

    def _validate_config(self):
        """Validate configuration parameters."""
        self._validate_thresholds()
        self._validate_dynamics_params()
        self._validate_domain_params()

    def _validate_thresholds(self):
        """Validate threshold parameters."""
        if self.min_pocket_volume <= 0:
            raise ValueError("min_pocket_volume must be positive")
        if self.pocket_max_depth <= 0:
            raise ValueError("pocket_max_depth must be positive")
        if self.cavity_min_depth <= 0:
            raise ValueError("cavity_min_depth must be positive")
        if self.interface_min_area <= 0:
            raise ValueError("interface_min_area must be positive")
        if self.interface_distance_cutoff <= 0:
            raise ValueError("interface_distance_cutoff must be positive")
        if self.surface_probe_radius <= 0:
            raise ValueError("surface_probe_radius must be positive")
        if self.clash_distance_cutoff <= 0:
            raise ValueError("clash_distance_cutoff must be positive")
        if self.hbond_distance_cutoff <= 0:
            raise ValueError("hbond_distance_cutoff must be positive")
        if not 0 <= self.hbond_angle_cutoff <= 180:
            raise ValueError("hbond_angle_cutoff must be between 0 and 180")
        if self.hydrophobic_contact_cutoff <= 0:
            raise ValueError("hydrophobic_contact_cutoff must be positive")
        if self.residue_exposure_threshold <= 0:
            raise ValueError("residue_exposure_threshold must be positive")
        if self.contact_distance <= 0:
            raise ValueError("contact_distance must be positive")
        if self.clash_overlap <= 0:
            raise ValueError("clash_overlap must be positive")

    def _validate_dynamics_params(self):
        """Validate dynamics analysis parameters."""
        if self.flexibility_cutoff <= 0:
            raise ValueError("flexibility_cutoff must be positive")
        if self.contact_cutoff <= 0:
            raise ValueError("contact_cutoff must be positive")
        if self.n_modes < 1:
            raise ValueError("n_modes must be at least 1")

    def _validate_domain_params(self):
        """Validate domain analysis parameters."""
        if self.domain_analysis_params["min_domain_size"] < 1:
            raise ValueError("min_domain_size must be at least 1")
        if self.domain_analysis_params["contact_cutoff"] <= 0:
            raise ValueError("domain contact_cutoff must be positive")
        if self.domain_analysis_params["interface_cutoff"] <= 0:
            raise ValueError("interface_cutoff must be positive")
        if not 0 <= self.domain_analysis_params["hinge_score_cutoff"] <= 1:
            raise ValueError("hinge_score_cutoff must be between 0 and 1")


@dataclass
class AlphaFoldConfig:
    """Configuration for AlphaFold integration.

    Attributes:
        model_preset: Preset model configuration ("monomer" or "multimer")
        num_recycle: Number of recycle iterations
        use_templates: Whether to use templates
        max_templates: Maximum number of templates to use
        use_amber: Whether to use AMBER relaxation
        use_pae: Whether to use predicted aligned error
        pae_power: Power for PAE-based contact weighting
        min_plddt: Minimum pLDDT score to consider
        confidence_bins: pLDDT score bins for confidence levels
    """

    model_preset: str = "monomer"
    num_recycle: int = 3
    use_templates: bool = True
    max_templates: int = 4
    use_amber: bool = True
    use_pae: bool = True
    pae_power: float = 1.0
    min_plddt: float = 70.0
    confidence_bins: Dict[str, float] = field(
        default_factory=lambda: {
            "very_high": 90.0,
            "high": 70.0,
            "medium": 50.0,
            "low": 0.0,
        }
    )

    def __post_init__(self):
        """Validate configuration after initialization."""
        if self.model_preset not in ["monomer", "multimer"]:
            raise ValueError("model_preset must be 'monomer' or 'multimer'")
        if self.num_recycle < 1:
            raise ValueError("num_recycle must be at least 1")
        if self.max_templates < 0:
            raise ValueError("max_templates must be non-negative")
        if not 0 <= self.pae_power <= 2:
            raise ValueError("pae_power must be between 0 and 2")
        if not 0 <= self.min_plddt <= 100:
            raise ValueError("min_plddt must be between 0 and 100")
        self._validate_confidence_bins()

    def _validate_confidence_bins(self):
        """Validate confidence bin thresholds."""
        prev_value = 100.0
        for key in ["very_high", "high", "medium", "low"]:
            if key not in self.confidence_bins:
                raise ValueError(f"Missing confidence bin: {key}")
            value = self.confidence_bins[key]
            if not 0 <= value <= 100:
                raise ValueError(f"Confidence threshold {key} must be between 0 and 100")
            if value > prev_value:
                raise ValueError("Confidence thresholds must be monotonically decreasing")
            prev_value = value


@dataclass
class DynamicsAnalysisConfig:
    """Configuration for protein dynamics analysis.

    Attributes:
        protein_config: Base protein analysis configuration
        use_bfactors: Whether to use B-factors for flexibility analysis
        use_normal_modes: Whether to use normal modes analysis
        use_contacts: Whether to use contact network analysis
        use_domain_motions: Whether to analyze domain motions
        use_correlations: Whether to analyze residue correlations
        advanced_params: Additional parameters for advanced analysis
    """

    protein_config: ProteinAnalysisConfig = field(default_factory=lambda: ProteinAnalysisConfig(alphafold_config=AlphaFoldConfig()))
    use_bfactors: bool = True
    use_normal_modes: bool = True
    use_contacts: bool = True
    use_domain_motions: bool = True
    use_correlations: bool = True
    advanced_params: Dict[str, Any] = field(
        default_factory=lambda: {
            "mode_cutoff": 0.1,
            "correlation_cutoff": 0.5,
            "flexibility_window": 5,
            "contact_weight": 1.0,
            "use_weighted_modes": True,
            "use_pae_weights": True,  # Weight contacts by PAE confidence
            "pae_power": 1.0,  # Power for PAE-based weights
            "min_plddt": 70.0,  # Minimum pLDDT score for residue inclusion
            "plddt_weight": True,  # Weight by pLDDT confidence scores
        }
    )

    def __post_init__(self):
        """Validate configuration after initialization."""
        self._validate_config()

    def _validate_config(self):
        """Validate configuration parameters."""
        if not 0 <= self.advanced_params["mode_cutoff"] <= 1:
            raise ValueError("mode_cutoff must be between 0 and 1")
        if not 0 <= self.advanced_params["correlation_cutoff"] <= 1:
            raise ValueError("correlation_cutoff must be between 0 and 1")
        if self.advanced_params["flexibility_window"] < 1:
            raise ValueError("flexibility_window must be at least 1")
        if self.advanced_params["contact_weight"] <= 0:
            raise ValueError("contact_weight must be positive")
