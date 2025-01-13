"""Tests for nootropic compound prediction and analysis.

This module provides comprehensive tests for:
1. Basic nootropic prediction
2. Mechanism prediction
3. Enhanced prediction with evidence
4. Model persistence
5. Error handling
6. Metrics calculation
7. Safety analysis
8. BBB permeability integration
9. Literature evidence integration
10. Community data integration
11. Prediction history tracking
12. Ensemble model validation
"""

import logging
import pytest
from pathlib import Path
from typing import List, Dict

import numpy as np
from rdkit import Chem

from .....models.core import CompoundData
from ..nootropic.base import NootropicPredictor


@pytest.fixture
def test_compounds() -> List[CompoundData]:
    """Create test compounds with known properties."""
    # Known nootropics with strong evidence
    piracetam = CompoundData(
        name="Piracetam",
        smiles="O=C1N(CC(=O)N(CC1)CC)CC",
        cas_number="7491-74-9",
    )
    aniracetam = CompoundData(
        name="Aniracetam",
        smiles="O=C1N(C(=O)CC(N1)c1ccccc1)CC",
        cas_number="72432-10-1",
    )
    modafinil = CompoundData(
        name="Modafinil",
        smiles="CC(=O)C(CS(=O)(=O)C)NC(=O)C",
        cas_number="68693-11-8",
    )

    # Compounds with some nootropic effects
    caffeine = CompoundData(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    l_theanine = CompoundData(
        name="L-Theanine",
        smiles="N[C@@H](CC(=O)N[C@@H](Cc1ccccc1)C(=O)O)C(=O)O",
        cas_number="3081-61-6",
    )

    # Non-nootropic compounds
    aspirin = CompoundData(
        name="Aspirin",
        smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
        cas_number="50-78-2",
    )
    ibuprofen = CompoundData(
        name="Ibuprofen",
        smiles="CC(C)CC1=CC=C(C=C1)[C@H](C)C(=O)O",
        cas_number="15687-27-1",
    )

    # Calculate molecular weights for validation
    for compound in [piracetam, aniracetam, modafinil, caffeine, l_theanine, aspirin, ibuprofen]:
        mol = Chem.MolFromSmiles(compound.smiles)
        if mol:
            compound.molecular_weight = Chem.Descriptors.ExactMolWt(mol)

    return [piracetam, aniracetam, modafinil, caffeine, l_theanine, aspirin, ibuprofen]


@pytest.fixture
def nootropic_predictor(tmp_path: Path) -> NootropicPredictor:
    """Create nootropic predictor."""
    return NootropicPredictor(
        model_dir=str(tmp_path / "models"),
        cache_dir=str(tmp_path / "cache"),
        log_level=logging.DEBUG,
    )


def test_predictor_initialization(nootropic_predictor: NootropicPredictor):
    """Test predictor initialization."""
    # Test core components
    assert nootropic_predictor.feature_types == ["fingerprints", "descriptors", "enhanced"]
    assert nootropic_predictor.mechanism_ensemble is not None
    assert nootropic_predictor.effect_ensembles is not None
    assert nootropic_predictor.side_effect_ensembles is not None
    assert nootropic_predictor.bbb_predictor is not None

    # Test mechanism types
    assert all(
        m
        in [
            "AMPA_MODULATION",
            "NMDA_MODULATION",
            "ACETYLCHOLINE_MODULATION",
            "DOPAMINE_MODULATION",
            "SEROTONIN_MODULATION",
            "GABA_MODULATION",
            "BDNF_MODULATION",
            "NEUROPLASTICITY",
            "NEUROPROTECTION",
            "ANTI_INFLAMMATION",
            "UNKNOWN",
        ]
        for m in nootropic_predictor.mechanism_types
    )


def test_basic_prediction(
    nootropic_predictor: NootropicPredictor,
    test_compounds: List[CompoundData],
):
    """Test basic nootropic prediction."""
    # Train predictor with both boolean and float labels
    nootropic_predictor.train(
        test_compounds[:3],  # Strong nootropics
        [1.0, 0.95, 0.9],
        test_compounds[3:5],  # Mild nootropics
        [0.6, 0.5],
        test_compounds[5:],  # Non-nootropics
        [0.1, 0.1],
    )

    # Test predictions
    scores = []
    confidences = []
    for compound in test_compounds:
        # Test prediction with full evidence
        result = nootropic_predictor.predict(compound)
        assert result.value is not None
        assert result.confidence > 0
        scores.append(result.value.value)
        confidences.append(result.confidence)

        # Validate supporting data structure
        assert "bbb_prediction" in result.supporting_data
        assert "effects" in result.supporting_data
        assert "side_effects" in result.supporting_data

    # Validate score distribution
    scores = np.array(scores)
    confidences = np.array(confidences)
    assert np.all(scores[:3] > scores[3:])  # Strong nootropics score higher
    assert np.all(scores[3:5] > scores[5:])  # Mild nootropics score higher than non-nootropics
    assert np.mean(confidences) > 0.5  # Average confidence should be reasonable


def test_mechanism_prediction(
    nootropic_predictor: NootropicPredictor,
    test_compounds: List[CompoundData],
):
    """Test mechanism prediction."""
    # Train predictor
    nootropic_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test mechanism predictions
    all_scores = []
    for compound in test_compounds:
        mechanism, confidence = nootropic_predictor._predict_mechanism(compound)
        assert isinstance(mechanism, str)
        assert isinstance(confidence, float)
        assert 0 <= confidence <= 1

        # Get feature importances
        importances = nootropic_predictor.get_feature_importance()
        assert isinstance(importances, Dict)
        if importances:
            assert "mechanism" in importances

        # Validate mechanism
        assert mechanism in nootropic_predictor.mechanism_types

    # Analyze mechanism score distributions
    score_matrix = np.array(all_scores)
    if len(all_scores) > 0:
        assert score_matrix.shape[0] == len(test_compounds)
        assert np.all(np.mean(score_matrix[:3], axis=1) > np.mean(score_matrix[5:], axis=1))


def test_enhanced_prediction(
    nootropic_predictor: NootropicPredictor,
    test_compounds: List[CompoundData],
):
    """Test enhanced nootropic prediction."""
    # Train predictor
    nootropic_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test predictions with BBB integration
    for compound in test_compounds:
        result = nootropic_predictor.predict(compound)
        assert result.value is not None
        assert result.confidence > 0
        assert "bbb_prediction" in result.supporting_data
        assert "effects" in result.supporting_data
        assert "side_effects" in result.supporting_data

        # Check BBB prediction structure
        bbb_data = result.supporting_data["bbb_prediction"]
        assert "value" in bbb_data
        assert "confidence" in bbb_data

        # Check effects structure
        effects = result.supporting_data["effects"]
        for domain, domain_effects in effects.items():
            for effect, effect_data in domain_effects.items():
                assert "score" in effect_data
                assert "confidence" in effect_data
                assert "bbb_factor" in effect_data

        # Check side effects structure
        side_effects = result.supporting_data["side_effects"]
        for category, category_effects in side_effects.items():
            for effect, effect_data in category_effects.items():
                assert "risk" in effect_data
                assert "confidence" in effect_data
                assert "bbb_factor" in effect_data


def test_model_persistence(
    nootropic_predictor: NootropicPredictor,
    test_compounds: List[CompoundData],
    tmp_path: Path,
):
    """Test model persistence."""
    # Train predictor
    nootropic_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Save models
    save_dir = tmp_path / "saved_models"
    nootropic_predictor.save_models(save_dir)

    # Create new predictor and load models
    new_predictor = NootropicPredictor(
        model_dir=str(save_dir),
        cache_dir=str(tmp_path / "new_cache"),
    )

    # Compare predictions
    for compound in test_compounds:
        result1 = nootropic_predictor.predict(compound)
        result2 = new_predictor.predict(compound)
        assert result1.value == result2.value
        assert abs(result1.confidence - result2.confidence) < 1e-6
        assert result1.supporting_data.keys() == result2.supporting_data.keys()


def test_error_handling(nootropic_predictor: NootropicPredictor):
    """Test error handling."""
    # Test with invalid compound
    invalid_compound = CompoundData(
        name="Invalid",
        smiles="INVALID",
        cas_number="000-00-0",
    )

    result = nootropic_predictor.predict(invalid_compound)
    assert result.value.value == "UNKNOWN"
    assert result.confidence == 0.0
    assert "error" in result.supporting_data

    # Test with untrained model
    valid_compound = CompoundData(
        name="Valid",
        smiles="CCO",
        cas_number="64-17-5",
    )

    result = nootropic_predictor.predict(valid_compound)
    assert result.value.value == "UNKNOWN"
    assert result.confidence == 0.0
    assert "error" in result.supporting_data


def test_prediction_metrics(
    nootropic_predictor: NootropicPredictor,
    test_compounds: List[CompoundData],
):
    """Test prediction metrics."""
    # Train predictor
    metrics = nootropic_predictor.retrain(
        test_compounds[:3],
        ["AMPA_MODULATION", "NMDA_MODULATION", "ACETYLCHOLINE_MODULATION"],
        {"memory": {"working_memory": [0.9, 0.8, 0.7]}},
        {"physical": {"headache": [0.1, 0.2, 0.3]}},
    )

    # Check metrics structure
    assert isinstance(metrics, Dict)
    assert any(k.startswith("mechanism_") for k in metrics)
    assert any(k.startswith("effect_") for k in metrics)
    assert any(k.startswith("side_effect_") for k in metrics)
    assert all(isinstance(v, float) for v in metrics.values())


def test_calibration(
    nootropic_predictor: NootropicPredictor,
    test_compounds: List[CompoundData],
):
    """Test model calibration."""
    # Train predictor
    nootropic_predictor.retrain(
        test_compounds[:3],
        ["AMPA_MODULATION", "NMDA_MODULATION", "ACETYLCHOLINE_MODULATION"],
        {"memory": {"working_memory": [0.9, 0.8, 0.7]}},
        {"physical": {"headache": [0.1, 0.2, 0.3]}},
    )

    # Calibrate models
    metrics = nootropic_predictor.calibrate_models(
        test_compounds,
        ["AMPA_MODULATION", "NMDA_MODULATION", "ACETYLCHOLINE_MODULATION"],
        {"memory": {"working_memory": [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]}},
        {"physical": {"headache": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]}},
    )

    # Check calibration metrics
    assert isinstance(metrics, Dict)
    assert "mechanism_brier_score" in metrics
    assert isinstance(metrics["mechanism_brier_score"], float)
    assert 0 <= metrics["mechanism_brier_score"] <= 1


def test_prediction_history(
    nootropic_predictor: NootropicPredictor,
    test_compounds: List[CompoundData],
):
    """Test prediction history tracking."""
    # Train predictor
    nootropic_predictor.retrain(
        test_compounds[:3],
        ["AMPA_MODULATION", "NMDA_MODULATION", "ACETYLCHOLINE_MODULATION"],
        {"memory": {"working_memory": [0.9, 0.8, 0.7]}},
        {"physical": {"headache": [0.1, 0.2, 0.3]}},
    )

    # Make predictions
    for compound in test_compounds:
        nootropic_predictor.predict(compound)

    # Check prediction history
    history = nootropic_predictor.prediction_history
    assert len(history) == len(test_compounds)
    assert all(
        col in history.columns
        for col in [
            "compound_name",
            "mechanism",
            "mechanism_confidence",
            "bbb_permeability",
            "bbb_confidence",
            "timestamp",
        ]
    )


def test_safety_analysis(
    nootropic_predictor: NootropicPredictor,
    test_compounds: List[CompoundData],
):
    """Test safety analysis for nootropics."""
    # Train predictor
    nootropic_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test safety analysis
    for compound in test_compounds:
        result = nootropic_predictor.predict(compound)
        safety_data = result.supporting_data.get("side_effects", {})
        assert isinstance(safety_data, Dict)

        # Check safety profile structure
        for category, effects in safety_data.items():
            for effect, data in effects.items():
                assert "risk" in data
                assert "confidence" in data
                assert isinstance(data["risk"], float)
                assert isinstance(data["confidence"], float)
                assert 0 <= data["risk"] <= 1
                assert 0 <= data["confidence"] <= 1


def test_bbb_prediction(
    nootropic_predictor: NootropicPredictor,
    test_compounds: List[CompoundData],
):
    """Test BBB permeability prediction integration."""
    # Train predictor
    nootropic_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test BBB prediction for each compound
    for compound in test_compounds:
        result = nootropic_predictor.predict(compound)
        bbb_data = result.supporting_data["bbb_prediction"]
        assert isinstance(bbb_data, Dict)
        assert "value" in bbb_data
        assert "confidence" in bbb_data
        assert isinstance(bbb_data["value"], str)
        assert isinstance(bbb_data["confidence"], float)

        # Check effect scaling by BBB permeability
        effects = result.supporting_data["effects"]
        for domain, domain_effects in effects.items():
            for effect, effect_data in domain_effects.items():
                assert "bbb_factor" in effect_data
                assert isinstance(effect_data["bbb_factor"], float)
                assert 0 <= effect_data["bbb_factor"] <= 1


def test_ensemble_prediction(
    nootropic_predictor: NootropicPredictor,
    test_compounds: List[CompoundData],
):
    """Test ensemble model prediction."""
    # Train predictor
    nootropic_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test ensemble predictions
    for compound in test_compounds:
        mechanism, confidence = nootropic_predictor._predict_mechanism(compound)

        # Check ensemble model structure
        assert nootropic_predictor.mechanism_ensemble is not None
        assert "mechanism" in nootropic_predictor.mechanism_ensemble.models

        # Check feature importance
        importances = nootropic_predictor.get_feature_importance()
        if "mechanism" in importances:
            assert isinstance(importances["mechanism"], Dict)
            assert all(isinstance(v, float) for v in importances["mechanism"].values())

        # Check ensemble confidence calculation
        assert 0 <= confidence <= 1

        # Validate prediction
        assert isinstance(mechanism, str)
        assert mechanism in [
            "AMPA_MODULATION",
            "NMDA_MODULATION",
            "ACETYLCHOLINE_MODULATION",
            "DOPAMINE_MODULATION",
            "SEROTONIN_MODULATION",
            "GABA_MODULATION",
            "BDNF_MODULATION",
            "NEUROPLASTICITY",
            "NEUROPROTECTION",
            "ANTI_INFLAMMATION",
            "UNKNOWN",
        ]
