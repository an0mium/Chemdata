"""Tests for nootropic compound prediction and analysis.

This module provides comprehensive tests for:
1. Basic nootropic prediction
2. Mechanism prediction
3. Enhanced prediction with evidence
4. Model persistence
5. Error handling
6. Metrics calculation
7. Safety analysis
"""

import logging
import pytest
from pathlib import Path
from typing import List, Dict

import numpy as np
from rdkit import Chem

from .....models.core import CompoundData
from ..nootropic import NootropicPredictor
from ..nootropic_enhanced import NootropicPredictorEnhanced


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


@pytest.fixture
def enhanced_predictor(tmp_path: Path) -> NootropicPredictorEnhanced:
    """Create enhanced nootropic predictor."""
    return NootropicPredictorEnhanced(
        model_dir=str(tmp_path / "models"),
        cache_dir=str(tmp_path / "cache"),
        log_level=logging.DEBUG,
    )


def test_predictor_initialization(nootropic_predictor: NootropicPredictor):
    """Test predictor initialization."""
    assert nootropic_predictor.ml_pipeline is not None
    assert nootropic_predictor.feature_types == ["morgan", "maccs", "rdkit"]
    assert nootropic_predictor.mechanism_types == [
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
    ]


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
        # Test float score prediction
        score, confidence = nootropic_predictor.predict(compound)
        assert isinstance(score, float)
        assert isinstance(confidence, float)
        assert 0 <= score <= 1
        assert 0 <= confidence <= 1
        scores.append(score)
        confidences.append(confidence)

        # Test boolean prediction
        is_nootropic = nootropic_predictor.predict_binary(compound, threshold=0.5)
        assert isinstance(is_nootropic, bool)

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
        mechanisms = nootropic_predictor.predict_mechanisms(compound)
        assert isinstance(mechanisms, Dict)
        assert all(m in nootropic_predictor.mechanism_types for m in mechanisms)
        assert all(isinstance(v, float) for v in mechanisms.values())
        assert all(0 <= v <= 1 for v in mechanisms.values())
        all_scores.append(list(mechanisms.values()))

        # Test mechanism thresholds
        primary = nootropic_predictor.get_primary_mechanisms(compound, threshold=0.5)
        assert isinstance(primary, List)
        assert all(m in nootropic_predictor.mechanism_types for m in primary)
        assert len(primary) <= len(mechanisms)

    # Analyze mechanism score distributions
    score_matrix = np.array(all_scores)
    assert score_matrix.shape == (len(test_compounds), len(nootropic_predictor.mechanism_types))
    assert np.all(np.mean(score_matrix[:3], axis=1) > np.mean(score_matrix[5:], axis=1))


def test_enhanced_prediction(
    enhanced_predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test enhanced nootropic prediction."""
    # Train predictor
    enhanced_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test predictions with evidence
    for compound in test_compounds:
        prediction = enhanced_predictor.predict_with_evidence(compound)
        assert isinstance(prediction, Dict)
        assert "nootropic_score" in prediction
        assert "confidence" in prediction
        assert "mechanisms" in prediction
        assert "literature_evidence" in prediction
        assert "community_data" in prediction
        assert "safety_profile" in prediction
        assert isinstance(prediction["nootropic_score"], float)
        assert isinstance(prediction["confidence"], float)
        assert isinstance(prediction["mechanisms"], Dict)
        assert isinstance(prediction["literature_evidence"], List)
        assert isinstance(prediction["community_data"], Dict)
        assert isinstance(prediction["safety_profile"], Dict)

        # Check literature evidence structure
        for evidence in prediction["literature_evidence"]:
            assert "source" in evidence
            assert "title" in evidence
            assert "year" in evidence
            assert "relevance" in evidence
            assert "findings" in evidence

        # Check community data structure
        assert "total_reports" in prediction["community_data"]
        assert "average_rating" in prediction["community_data"]
        assert "reported_effects" in prediction["community_data"]
        assert "reported_mechanisms" in prediction["community_data"]


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
        score1, conf1 = nootropic_predictor.predict(compound)
        score2, conf2 = new_predictor.predict(compound)
        assert abs(score1 - score2) < 1e-6
        assert abs(conf1 - conf2) < 1e-6

        mech1 = nootropic_predictor.predict_mechanisms(compound)
        mech2 = new_predictor.predict_mechanisms(compound)
        assert set(mech1.keys()) == set(mech2.keys())
        assert all(abs(mech1[k] - mech2[k]) < 1e-6 for k in mech1)


def test_error_handling(nootropic_predictor: NootropicPredictor):
    """Test error handling."""
    # Test with invalid compound
    invalid_compound = CompoundData(
        name="Invalid",
        smiles="INVALID",
        cas_number="000-00-0",
    )

    with pytest.raises(ValueError):
        nootropic_predictor.predict(invalid_compound)

    with pytest.raises(ValueError):
        nootropic_predictor.predict_mechanisms(invalid_compound)

    # Test with untrained model
    valid_compound = CompoundData(
        name="Valid",
        smiles="CCO",
        cas_number="64-17-5",
    )

    with pytest.raises(RuntimeError):
        nootropic_predictor.predict(valid_compound)

    with pytest.raises(RuntimeError):
        nootropic_predictor.predict_mechanisms(valid_compound)


def test_prediction_metrics(
    nootropic_predictor: NootropicPredictor,
    test_compounds: List[CompoundData],
):
    """Test prediction metrics."""
    # Train predictor
    metrics = nootropic_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Check regression metrics
    assert isinstance(metrics, Dict)
    assert "mse" in metrics
    assert "rmse" in metrics
    assert "mae" in metrics
    assert "r2" in metrics
    assert all(isinstance(v, float) for v in metrics.values())
    assert all(v >= 0 for v in metrics.values())

    # Check classification metrics
    binary_metrics = nootropic_predictor.get_binary_metrics(threshold=0.5)
    assert "accuracy" in binary_metrics
    assert "precision" in binary_metrics
    assert "recall" in binary_metrics
    assert "f1" in binary_metrics
    assert all(isinstance(v, float) for v in binary_metrics.values())
    assert all(0 <= v <= 1 for v in binary_metrics.values())


def test_enhanced_metrics(
    enhanced_predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test enhanced prediction metrics."""
    # Train predictor
    metrics = enhanced_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Check metrics
    assert isinstance(metrics, Dict)
    assert "ml_metrics" in metrics
    assert "literature_metrics" in metrics
    assert "community_metrics" in metrics
    assert "combined_metrics" in metrics
    assert all(isinstance(v, Dict) for v in metrics.values())
    assert all(all(isinstance(v, float) for v in m.values()) for m in metrics.values())
    assert all(all(v >= 0 for v in m.values()) for m in metrics.values())


def test_safety_analysis(
    enhanced_predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test safety analysis for nootropics."""
    # Train predictor
    enhanced_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test safety analysis
    for compound in test_compounds:
        analysis = enhanced_predictor.analyze_safety(compound)
        assert isinstance(analysis, Dict)
        assert "toxicity_score" in analysis
        assert "addiction_potential" in analysis
        assert "side_effects" in analysis
        assert "interactions" in analysis
        assert "contraindications" in analysis
        assert "long_term_risks" in analysis
        assert "tolerance_profile" in analysis
        assert "withdrawal_profile" in analysis
        assert isinstance(analysis["toxicity_score"], float)
        assert isinstance(analysis["addiction_potential"], float)
        assert isinstance(analysis["side_effects"], List)
        assert isinstance(analysis["interactions"], Dict)
        assert isinstance(analysis["contraindications"], List)
        assert isinstance(analysis["long_term_risks"], List)
        assert isinstance(analysis["tolerance_profile"], Dict)
        assert isinstance(analysis["withdrawal_profile"], Dict)


def test_bbb_prediction(
    enhanced_predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test BBB permeability prediction integration."""
    # Train predictor
    enhanced_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test BBB prediction for each compound
    for compound in test_compounds:
        bbb_result = enhanced_predictor.predict_bbb(compound)
        assert isinstance(bbb_result, Dict)
        assert "value" in bbb_result
        assert "confidence" in bbb_result
        assert "supporting_data" in bbb_result
        assert isinstance(bbb_result["value"], str)
        assert isinstance(bbb_result["confidence"], float)
        assert isinstance(bbb_result["supporting_data"], Dict)

        # Check supporting data structure
        supporting_data = bbb_result["supporting_data"]
        assert "transporters" in supporting_data
        assert "mechanisms" in supporting_data
        assert "literature_evidence" in supporting_data
        assert isinstance(supporting_data["transporters"], Dict)
        assert isinstance(supporting_data["mechanisms"], Dict)
        assert isinstance(supporting_data["literature_evidence"], List)

        # Validate prediction impact
        prediction = enhanced_predictor.predict_with_evidence(compound)
        assert "bbb_prediction" in prediction
        assert prediction["bbb_prediction"] == bbb_result

        # Check effect scaling by BBB permeability
        for domain, effects in prediction["cognitive"].items():
            for effect, data in effects.items():
                assert "bbb_adjusted_score" in data
                assert isinstance(data["bbb_adjusted_score"], float)
                assert data["bbb_adjusted_score"] <= data["score"]

        # Check confidence adjustment
        assert "bbb_confidence_factor" in prediction
        assert isinstance(prediction["bbb_confidence_factor"], float)
        assert 0 <= prediction["bbb_confidence_factor"] <= 1


def test_literature_analysis(
    enhanced_predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test literature analysis for nootropics."""
    # Train predictor
    enhanced_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test literature analysis
    for compound in test_compounds:
        analysis = enhanced_predictor.analyze_literature(compound)
        assert isinstance(analysis, Dict)
        assert "total_papers" in analysis
        assert "publication_trend" in analysis
        assert "key_findings" in analysis
        assert "mechanism_evidence" in analysis
        assert "safety_evidence" in analysis
        assert "clinical_trials" in analysis
        assert isinstance(analysis["total_papers"], int)
        assert isinstance(analysis["publication_trend"], Dict)
        assert isinstance(analysis["key_findings"], List)
        assert isinstance(analysis["mechanism_evidence"], Dict)
        assert isinstance(analysis["safety_evidence"], Dict)
        assert isinstance(analysis["clinical_trials"], List)
