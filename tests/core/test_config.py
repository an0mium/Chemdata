"""Tests for core configuration classes."""

import pytest
from binding_data_processor.core.config import (
    ProteinAnalysisConfig,
    DynamicsAnalysisConfig,
    AlphaFoldConfig,
)


def test_protein_analysis_config_defaults():
    """Test default values for ProteinAnalysisConfig."""
    config = ProteinAnalysisConfig()

    # Test model configuration defaults
    assert config.model_path is None
    assert config.cache_dir is None
    assert config.dssp_binary is None

    # Test AlphaFold integration defaults
    assert config.use_alphafold is True
    assert isinstance(config.alphafold_config, AlphaFoldConfig)
    assert config.use_alphafold_api is True
    assert config.alphafold_api_base_url == "https://alphafold.ebi.ac.uk/api"
    assert config.alphafold_api_timeout == 300
    assert config.use_alphafold_pae is True
    assert config.use_alphafold_confidence is True
    assert config.min_plddt_score == 70.0
    assert config.pae_power == 1.0

    # Test analysis configuration defaults
    assert config.analyze_pockets is True
    assert config.analyze_dynamics is True
    assert config.analyze_conservation is True
    assert config.analyze_interfaces is True
    assert config.analyze_quality is True

    # Test feature configuration defaults
    assert config.use_pssm is True
    assert config.use_dssp is True
    assert config.use_surface_features is True
    assert config.use_conservation is True
    assert config.use_dynamics is True
    assert config.use_quality is True
    assert config.use_surface is True
    assert config.use_pockets is True
    assert config.use_interfaces is True

    # Test performance settings defaults
    assert config.batch_size == 100
    assert config.num_workers == 4
    assert config.n_jobs == 1
    assert config.cache_results is True

    # Test analysis thresholds defaults
    assert config.min_pocket_volume == 100.0
    assert config.pocket_max_depth == 20.0
    assert config.cavity_min_depth == 2.0
    assert config.interface_min_area == 100.0
    assert config.interface_distance_cutoff == 5.0
    assert config.surface_probe_radius == 1.4
    assert config.clash_distance_cutoff == 0.4
    assert config.hbond_distance_cutoff == 3.5
    assert config.hbond_angle_cutoff == 30.0
    assert config.hydrophobic_contact_cutoff == 4.5
    assert config.residue_exposure_threshold == 2.8
    assert config.contact_distance == 5.0
    assert config.clash_overlap == 0.4

    # Test domain analysis defaults
    assert config.include_domain_analysis is True
    assert config.domain_analysis_params["use_structure_based"] is True
    assert config.domain_analysis_params["use_sequence_based"] is False
    assert config.domain_analysis_params["min_domain_size"] == 30
    assert config.domain_analysis_params["contact_cutoff"] == 8.0
    assert config.domain_analysis_params["interface_cutoff"] == 5.0
    assert config.domain_analysis_params["hinge_score_cutoff"] == 0.5

    # Test dynamics analysis defaults
    assert config.flexibility_cutoff == 1.0
    assert config.contact_cutoff == 8.0
    assert config.n_modes == 10
    assert config.logging_level == "INFO"


def test_protein_analysis_config_validation():
    """Test validation of ProteinAnalysisConfig parameters."""
    # Test invalid flexibility cutoff
    with pytest.raises(ValueError, match="flexibility_cutoff must be positive"):
        ProteinAnalysisConfig(flexibility_cutoff=0)

    # Test invalid contact cutoff
    with pytest.raises(ValueError, match="contact_cutoff must be positive"):
        ProteinAnalysisConfig(contact_cutoff=-1)

    # Test invalid n_modes
    with pytest.raises(ValueError, match="n_modes must be at least 1"):
        ProteinAnalysisConfig(n_modes=0)

    # Test invalid domain parameters
    with pytest.raises(ValueError, match="min_domain_size must be at least 1"):
        ProteinAnalysisConfig(domain_analysis_params={"min_domain_size": 0})

    with pytest.raises(ValueError, match="domain contact_cutoff must be positive"):
        ProteinAnalysisConfig(domain_analysis_params={"contact_cutoff": -1})

    with pytest.raises(ValueError, match="interface_cutoff must be positive"):
        ProteinAnalysisConfig(domain_analysis_params={"interface_cutoff": 0})

    with pytest.raises(ValueError, match="hinge_score_cutoff must be between 0 and 1"):
        ProteinAnalysisConfig(domain_analysis_params={"hinge_score_cutoff": 1.5})


def test_alphafold_config_defaults():
    """Test default values for AlphaFoldConfig."""
    config = AlphaFoldConfig()

    # Test basic settings
    assert config.model_preset == "monomer"
    assert config.num_recycle == 3
    assert config.use_templates is True
    assert config.max_templates == 4
    assert config.use_amber is True
    assert config.use_pae is True
    assert config.pae_power == 1.0
    assert config.min_plddt == 70.0

    # Test confidence bins
    assert config.confidence_bins["very_high"] == 90.0
    assert config.confidence_bins["high"] == 70.0
    assert config.confidence_bins["medium"] == 50.0
    assert config.confidence_bins["low"] == 0.0


def test_alphafold_config_validation():
    """Test validation of AlphaFoldConfig parameters."""
    # Test invalid model preset
    with pytest.raises(ValueError, match="model_preset must be 'monomer' or 'multimer'"):
        AlphaFoldConfig(model_preset="invalid")

    # Test invalid num_recycle
    with pytest.raises(ValueError, match="num_recycle must be at least 1"):
        AlphaFoldConfig(num_recycle=0)

    # Test invalid max_templates
    with pytest.raises(ValueError, match="max_templates must be non-negative"):
        AlphaFoldConfig(max_templates=-1)

    # Test invalid pae_power
    with pytest.raises(ValueError, match="pae_power must be between 0 and 2"):
        AlphaFoldConfig(pae_power=2.5)

    # Test invalid min_plddt
    with pytest.raises(ValueError, match="min_plddt must be between 0 and 100"):
        AlphaFoldConfig(min_plddt=150)

    # Test invalid confidence bins
    with pytest.raises(ValueError, match="Missing confidence bin: very_high"):
        AlphaFoldConfig(confidence_bins={"high": 70.0})

    with pytest.raises(ValueError, match="Confidence threshold very_high must be between 0 and 100"):
        AlphaFoldConfig(confidence_bins={"very_high": 150.0, "high": 70.0, "medium": 50.0, "low": 0.0})

    with pytest.raises(ValueError, match="Confidence thresholds must be monotonically decreasing"):
        AlphaFoldConfig(confidence_bins={"very_high": 70.0, "high": 90.0, "medium": 50.0, "low": 0.0})


def test_dynamics_analysis_config_defaults():
    """Test default values for DynamicsAnalysisConfig."""
    config = DynamicsAnalysisConfig()

    # Test basic configuration defaults
    assert isinstance(config.protein_config, ProteinAnalysisConfig)
    assert config.use_bfactors is True
    assert config.use_normal_modes is True
    assert config.use_contacts is True
    assert config.use_domain_motions is True
    assert config.use_correlations is True

    # Test advanced parameters defaults
    assert config.advanced_params["mode_cutoff"] == 0.1
    assert config.advanced_params["correlation_cutoff"] == 0.5
    assert config.advanced_params["flexibility_window"] == 5
    assert config.advanced_params["contact_weight"] == 1.0
    assert config.advanced_params["use_weighted_modes"] is True
    assert config.advanced_params["use_pae_weights"] is True
    assert config.advanced_params["pae_power"] == 1.0
    assert config.advanced_params["min_plddt"] == 70.0
    assert config.advanced_params["plddt_weight"] is True


def test_dynamics_analysis_config_validation():
    """Test validation of DynamicsAnalysisConfig parameters."""
    # Test invalid mode cutoff
    with pytest.raises(ValueError, match="mode_cutoff must be between 0 and 1"):
        DynamicsAnalysisConfig(advanced_params={"mode_cutoff": 1.5})

    # Test invalid correlation cutoff
    with pytest.raises(ValueError, match="correlation_cutoff must be between 0 and 1"):
        DynamicsAnalysisConfig(advanced_params={"correlation_cutoff": -0.1})

    # Test invalid flexibility window
    with pytest.raises(ValueError, match="flexibility_window must be at least 1"):
        DynamicsAnalysisConfig(advanced_params={"flexibility_window": 0})

    # Test invalid contact weight
    with pytest.raises(ValueError, match="contact_weight must be positive"):
        DynamicsAnalysisConfig(advanced_params={"contact_weight": -1.0})


def test_protein_analysis_config_residue_properties():
    """Test residue properties in ProteinAnalysisConfig."""
    config = ProteinAnalysisConfig()

    # Test hydrophobicity values
    assert config.residue_properties["hydrophobicity"]["ILE"] == 4.5
    assert config.residue_properties["hydrophobicity"]["GLY"] == -0.4
    assert config.residue_properties["hydrophobicity"]["ARG"] == -4.5

    # Test charge values
    assert config.residue_properties["charge"]["ARG"] == 1
    assert config.residue_properties["charge"]["ASP"] == -1
    assert config.residue_properties["charge"]["HIS"] == 0.5

    # Test volume values
    assert config.residue_properties["volume"]["ALA"] == 88.6
    assert config.residue_properties["volume"]["TRP"] == 227.8
    assert config.residue_properties["volume"]["VAL"] == 140.0


def test_protein_analysis_config_atom_radii():
    """Test atom radii in ProteinAnalysisConfig."""
    config = ProteinAnalysisConfig()

    # Test common atom radii
    assert config.atom_radii["C"] == 1.7
    assert config.atom_radii["N"] == 1.55
    assert config.atom_radii["O"] == 1.52
    assert config.atom_radii["H"] == 1.2

    # Test halogen radii
    assert config.atom_radii["F"] == 1.47
    assert config.atom_radii["Cl"] == 1.75
    assert config.atom_radii["Br"] == 1.85
    assert config.atom_radii["I"] == 1.98


def test_protein_analysis_config_quality_thresholds():
    """Test quality thresholds in ProteinAnalysisConfig."""
    config = ProteinAnalysisConfig()

    # Test quality metric thresholds
    assert config.quality_thresholds["clash_score"] == 20.0
    assert config.quality_thresholds["rama_favored"] == 0.98
    assert config.quality_thresholds["rotamer_outliers"] == 0.01
    assert config.quality_thresholds["cbeta_deviations"] == 0.01


def test_protein_analysis_config_confidence_thresholds():
    """Test confidence thresholds in ProteinAnalysisConfig."""
    config = ProteinAnalysisConfig()

    # Test confidence level thresholds
    assert config.confidence_thresholds["high"] == 90.0
    assert config.confidence_thresholds["medium"] == 70.0
    assert config.confidence_thresholds["low"] == 50.0
