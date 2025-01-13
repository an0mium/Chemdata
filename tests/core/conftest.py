"""Test fixtures for core module tests."""

import pytest
from binding_data_processor.core.config import ProteinAnalysisConfig, DynamicsAnalysisConfig


@pytest.fixture
def default_protein_config():
    """Return a default ProteinAnalysisConfig instance."""
    return ProteinAnalysisConfig()


@pytest.fixture
def custom_protein_config():
    """Return a ProteinAnalysisConfig instance with custom settings."""
    return ProteinAnalysisConfig(
        model_path="/path/to/model",
        cache_dir="/path/to/cache",
        dssp_binary="/usr/local/bin/dssp",
        batch_size=200,
        num_workers=8,
        n_jobs=4,
        min_pocket_volume=150.0,
        pocket_max_depth=25.0,
        cavity_min_depth=3.0,
        interface_min_area=150.0,
        domain_analysis_params={
            "use_structure_based": True,
            "use_sequence_based": True,
            "min_domain_size": 50,
            "contact_cutoff": 10.0,
            "interface_cutoff": 6.0,
            "hinge_score_cutoff": 0.7,
        },
    )


@pytest.fixture
def default_dynamics_config():
    """Return a default DynamicsAnalysisConfig instance."""
    return DynamicsAnalysisConfig()


@pytest.fixture
def custom_dynamics_config():
    """Return a DynamicsAnalysisConfig instance with custom settings."""
    return DynamicsAnalysisConfig(
        use_bfactors=False,
        use_normal_modes=True,
        use_contacts=True,
        use_domain_motions=True,
        use_correlations=True,
        advanced_params={
            "mode_cutoff": 0.2,
            "correlation_cutoff": 0.7,
            "flexibility_window": 7,
            "contact_weight": 1.5,
            "use_weighted_modes": False,
        },
    )
