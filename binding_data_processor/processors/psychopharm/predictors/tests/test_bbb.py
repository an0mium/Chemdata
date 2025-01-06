"""Tests for BBB permeability prediction."""

import logging
import pytest
from pathlib import Path
from typing import List, Dict

from .....models.core import CompoundData
from ..bbb import BBBPermeabilityPredictor
from ..bbb_enhanced import BBBPermeabilityPredictorEnhanced


@pytest.fixture
def test_compounds() -> List[CompoundData]:
    """Create test compounds."""
    # Known BBB permeable compounds
    caffeine = CompoundData(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    morphine = CompoundData(
        name="Morphine",
        smiles="CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
        cas_number="57-27-2",
    )

    # Known BBB impermeable compounds
    sucrose = CompoundData(
        name="Sucrose",
        smiles="C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O",
        cas_number="57-50-1",
    )
    insulin = CompoundData(
        name="Insulin",
        smiles="CC(C)CC1=CC=C(C=C1)CC(C(=O)O)N",  # Simplified structure
        cas_number="11061-68-0",
    )

    return [caffeine, morphine, sucrose, insulin]


@pytest.fixture
def bbb_predictor(tmp_path: Path) -> BBBPermeabilityPredictor:
    """Create BBB permeability predictor."""
    return BBBPermeabilityPredictor(
        model_dir=str(tmp_path / "models"),
        cache_dir=str(tmp_path / "cache"),
        log_level=logging.DEBUG,
    )


@pytest.fixture
def enhanced_predictor(tmp_path: Path) -> BBBPermeabilityPredictorEnhanced:
    """Create enhanced BBB permeability predictor."""
    return BBBPermeabilityPredictorEnhanced(
        model_dir=str(tmp_path / "models"),
        cache_dir=str(tmp_path / "cache"),
        log_level=logging.DEBUG,
    )


def test_predictor_initialization(bbb_predictor: BBBPermeabilityPredictor):
    """Test predictor initialization."""
    assert bbb_predictor.ml_pipeline is not None
    assert bbb_predictor.feature_types == ["morgan", "maccs", "rdkit"]
    assert bbb_predictor.permeability_types == [
        "PASSIVE_DIFFUSION",
        "ACTIVE_TRANSPORT",
        "EFFLUX_TRANSPORT",
        "CARRIER_MEDIATED",
    ]


def test_basic_prediction(
    bbb_predictor: BBBPermeabilityPredictor,
    test_compounds: List[CompoundData],
):
    """Test basic BBB permeability prediction."""
    # Train predictor
    bbb_predictor.train(test_compounds[:2], [1.0, 0.9])  # High permeability
    bbb_predictor.train(test_compounds[2:], [0.1, 0.1])  # Low permeability

    # Test predictions
    for compound in test_compounds:
        score, confidence = bbb_predictor.predict(compound)
        assert isinstance(score, float)
        assert isinstance(confidence, float)
        assert 0 <= score <= 1
        assert 0 <= confidence <= 1


def test_mechanism_prediction(
    bbb_predictor: BBBPermeabilityPredictor,
    test_compounds: List[CompoundData],
):
    """Test permeability mechanism prediction."""
    # Train predictor
    bbb_predictor.train(test_compounds[:2], [1.0, 0.9])
    bbb_predictor.train(test_compounds[2:], [0.1, 0.1])

    # Test mechanism predictions
    for compound in test_compounds:
        mechanisms = bbb_predictor.predict_mechanisms(compound)
        assert isinstance(mechanisms, Dict)
        assert all(m in bbb_predictor.permeability_types for m in mechanisms)
        assert all(isinstance(v, float) for v in mechanisms.values())
        assert all(0 <= v <= 1 for v in mechanisms.values())


def test_enhanced_prediction(
    enhanced_predictor: BBBPermeabilityPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test enhanced BBB permeability prediction."""
    # Train predictor
    enhanced_predictor.train(test_compounds[:2], [1.0, 0.9])
    enhanced_predictor.train(test_compounds[2:], [0.1, 0.1])

    # Test predictions with literature data
    for compound in test_compounds:
        prediction = enhanced_predictor.predict_with_evidence(compound)
        assert isinstance(prediction, Dict)
        assert "permeability_score" in prediction
        assert "confidence" in prediction
        assert "mechanisms" in prediction
        assert "literature_evidence" in prediction
        assert "experimental_data" in prediction
        assert "property_analysis" in prediction
        assert isinstance(prediction["permeability_score"], float)
        assert isinstance(prediction["confidence"], float)
        assert isinstance(prediction["mechanisms"], Dict)
        assert isinstance(prediction["literature_evidence"], List)
        assert isinstance(prediction["experimental_data"], Dict)
        assert isinstance(prediction["property_analysis"], Dict)


def test_model_persistence(
    bbb_predictor: BBBPermeabilityPredictor,
    test_compounds: List[CompoundData],
    tmp_path: Path,
):
    """Test model persistence."""
    # Train predictor
    bbb_predictor.train(test_compounds[:2], [1.0, 0.9])
    bbb_predictor.train(test_compounds[2:], [0.1, 0.1])

    # Save models
    save_dir = tmp_path / "saved_models"
    bbb_predictor.save_models(save_dir)

    # Create new predictor and load models
    new_predictor = BBBPermeabilityPredictor(
        model_dir=str(save_dir),
        cache_dir=str(tmp_path / "new_cache"),
    )

    # Compare predictions
    for compound in test_compounds:
        score1, conf1 = bbb_predictor.predict(compound)
        score2, conf2 = new_predictor.predict(compound)
        assert abs(score1 - score2) < 1e-6
        assert abs(conf1 - conf2) < 1e-6

        mech1 = bbb_predictor.predict_mechanisms(compound)
        mech2 = new_predictor.predict_mechanisms(compound)
        assert set(mech1.keys()) == set(mech2.keys())
        assert all(abs(mech1[k] - mech2[k]) < 1e-6 for k in mech1)


def test_error_handling(bbb_predictor: BBBPermeabilityPredictor):
    """Test error handling."""
    # Test with invalid compound
    invalid_compound = CompoundData(
        name="Invalid",
        smiles="INVALID",
        cas_number="000-00-0",
    )

    with pytest.raises(ValueError):
        bbb_predictor.predict(invalid_compound)

    with pytest.raises(ValueError):
        bbb_predictor.predict_mechanisms(invalid_compound)

    # Test with untrained model
    valid_compound = CompoundData(
        name="Valid",
        smiles="CCO",
        cas_number="64-17-5",
    )

    with pytest.raises(RuntimeError):
        bbb_predictor.predict(valid_compound)

    with pytest.raises(RuntimeError):
        bbb_predictor.predict_mechanisms(valid_compound)


def test_prediction_metrics(
    bbb_predictor: BBBPermeabilityPredictor,
    test_compounds: List[CompoundData],
):
    """Test prediction metrics."""
    # Train predictor
    metrics = bbb_predictor.train(
        test_compounds[:2],
        [1.0, 0.9],
        test_compounds[2:],
        [0.1, 0.1],
    )

    # Check metrics
    assert isinstance(metrics, Dict)
    assert "mse" in metrics
    assert "rmse" in metrics
    assert "mae" in metrics
    assert "r2" in metrics
    assert all(isinstance(v, float) for v in metrics.values())
    assert all(v >= 0 for v in metrics.values())


def test_enhanced_metrics(
    enhanced_predictor: BBBPermeabilityPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test enhanced prediction metrics."""
    # Train predictor
    metrics = enhanced_predictor.train(
        test_compounds[:2],
        [1.0, 0.9],
        test_compounds[2:],
        [0.1, 0.1],
    )

    # Check metrics
    assert isinstance(metrics, Dict)
    assert "ml_metrics" in metrics
    assert "literature_metrics" in metrics
    assert "experimental_metrics" in metrics
    assert "combined_metrics" in metrics
    assert all(isinstance(v, Dict) for v in metrics.values())
    assert all(all(isinstance(v, float) for v in m.values()) for m in metrics.values())
    assert all(all(v >= 0 for v in m.values()) for m in metrics.values())


def test_property_analysis(
    enhanced_predictor: BBBPermeabilityPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test property analysis for BBB permeability."""
    # Train predictor
    enhanced_predictor.train(test_compounds[:2], [1.0, 0.9])
    enhanced_predictor.train(test_compounds[2:], [0.1, 0.1])

    # Test property analysis
    for compound in test_compounds:
        analysis = enhanced_predictor.analyze_properties(compound)
        assert isinstance(analysis, Dict)
        assert "molecular_weight" in analysis
        assert "logp" in analysis
        assert "hbd" in analysis
        assert "hba" in analysis
        assert "tpsa" in analysis
        assert "rotatable_bonds" in analysis
        assert "rule_violations" in analysis
        assert all(isinstance(v, (float, int, list)) for v in analysis.values())
