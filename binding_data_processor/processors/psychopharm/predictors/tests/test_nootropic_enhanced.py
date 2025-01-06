"""Tests for enhanced nootropic prediction and analysis.

This module provides comprehensive tests for:
1. Basic nootropic prediction 
2. Mechanism prediction
3. Enhanced prediction with evidence
4. Cognitive domain prediction
5. Side effects prediction
6. Model persistence
7. Error handling
8. Metrics calculation
9. Safety analysis
10. Literature analysis
11. Prediction history tracking
"""

import logging
import pytest
from pathlib import Path
from typing import Dict, List, Set

import numpy as np
import pandas as pd
from rdkit import Chem

from .....models.core import CompoundData
from .....models.psychopharm import NootropicMechanism
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

    # Additional nootropics
    noopept = CompoundData(
        name="Noopept",
        smiles="CC(C)[C@H](NC(=O)[C@H](Cc1ccccc1)n1ncc2ccccc21)C(=O)N",
        cas_number="157115-85-0",
    )
    phenylpiracetam = CompoundData(
        name="Phenylpiracetam",
        smiles="O=C1N(CC(=O)N(CC1)CC(c1ccccc1))CC",
        cas_number="77472-70-9",
    )
    semax = CompoundData(
        name="Semax",
        smiles=(
            "CC(C)[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](CC(=O)N)"
            "NC(=O)CNC(=O)[C@H](N)CC(=O)N)C(=O)N[C@@H](CC(=O)N)C(=O)N"
        ),
        cas_number="80714-61-0",
    )

    # Compounds with anxiolytic effects
    phenibut = CompoundData(
        name="Phenibut",
        smiles="O=C(O)CCNC(c1ccccc1)H",
        cas_number="1078-21-3",
    )
    selank = CompoundData(
        name="Selank",
        smiles=(
            "CC(C)[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](CC(=O)N)"
            "NC(=O)CNC(=O)[C@H](N)CC(=O)N)C(=O)N[C@@H](CC(=O)N)C(=O)N"
        ),
        cas_number="129954-34-3",
    )

    # Compounds with neuroprotective effects
    cerebrolysin = CompoundData(
        name="Cerebrolysin",
        smiles=(
            "CC(C)[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](CC(=O)N)"
            "NC(=O)CNC(=O)[C@H](N)CC(=O)N)C(=O)N[C@@H](CC(=O)N)C(=O)N"
        ),
        cas_number="85166-18-9",
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
    compounds = [
        piracetam,
        aniracetam,
        modafinil,
        caffeine,
        l_theanine,
        noopept,
        phenylpiracetam,
        semax,
        phenibut,
        selank,
        cerebrolysin,
        aspirin,
        ibuprofen,
    ]
    for compound in compounds:
        mol = Chem.MolFromSmiles(compound.smiles)
        if mol:
            compound.molecular_weight = Chem.Descriptors.ExactMolWt(mol)

    return compounds


@pytest.fixture
def test_cognitive_domains() -> Dict[str, Set[str]]:
    """Create test cognitive domains."""
    return {
        "memory": {"working_memory", "long_term_memory", "recall_speed"},
        "attention": {"sustained_attention", "divided_attention", "focus"},
        "learning": {"skill_acquisition", "knowledge_retention", "comprehension"},
        "processing": {"processing_speed", "mental_clarity", "cognitive_flexibility"},
        "executive": {"planning", "decision_making", "impulse_control"},
    }


@pytest.fixture
def test_side_effects() -> Dict[str, Set[str]]:
    """Create test side effects."""
    return {
        "physical": {"headache", "nausea", "insomnia", "fatigue"},
        "cognitive": {"brain_fog", "anxiety", "irritability", "mood_changes"},
        "tolerance": {"rapid_tolerance", "withdrawal", "dependence"},
        "interactions": {"drug_interactions", "supplement_interactions"},
    }


@pytest.fixture
def predictor(
    tmp_path: Path,
    test_cognitive_domains: Dict[str, Set[str]],
    test_side_effects: Dict[str, Set[str]],
) -> NootropicPredictorEnhanced:
    """Create enhanced nootropic predictor."""
    return NootropicPredictorEnhanced(
        model_dir=str(tmp_path / "models"),
        cache_dir=str(tmp_path / "cache"),
        log_level=logging.DEBUG,
        cognitive_domains=test_cognitive_domains,
        side_effects=test_side_effects,
    )


def test_predictor_initialization(
    predictor: NootropicPredictorEnhanced,
    test_cognitive_domains: Dict[str, Set[str]],
    test_side_effects: Dict[str, Set[str]],
):
    """Test predictor initialization."""
    assert predictor.ml_pipeline is not None
    assert predictor.feature_types == ["morgan", "maccs", "rdkit"]
    assert predictor.mechanism_types == [
        "MEMORY_ENHANCEMENT",
        "FOCUS_IMPROVEMENT",
        "NEUROPROTECTION",
        "NEUROPLASTICITY",
        "CHOLINERGIC",
        "GLUTAMATERGIC",
        "DOPAMINERGIC",
        "SEROTONERGIC",
        "NOOTROPIC_SYNERGY",
        "COGNITIVE_MODULATION",
        "BRAIN_METABOLISM",
    ]

    # Check cognitive domains
    assert predictor.cognitive_domains == test_cognitive_domains
    assert predictor.side_effects == test_side_effects
    assert predictor.mechanism_ensemble is not None

    # Check ensemble initialization
    for domain, effects in test_cognitive_domains.items():
        assert domain in predictor.effect_ensembles
        for effect in effects:
            assert effect in predictor.effect_ensembles[domain]

    for category, effects in test_side_effects.items():
        assert category in predictor.side_effect_ensembles
        for effect in effects:
            assert effect in predictor.side_effect_ensembles[category]


def test_basic_prediction(
    predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test basic nootropic prediction."""
    # Train predictor with both boolean and float labels
    predictor.train(
        test_compounds[:3],  # Strong nootropics (piracetam, aniracetam, modafinil)
        [1.0, 0.95, 0.9],
        test_compounds[3:8],  # Additional nootropics (caffeine through semax)
        [0.7, 0.6, 0.8, 0.7, 0.8],
        test_compounds[8:],  # Anxiolytics, neuroprotectives, and non-nootropics
        [0.4, 0.3, 0.2, 0.1, 0.1],
    )

    # Test predictions
    scores = []
    confidences = []
    for compound in test_compounds:
        # Test float score prediction
        score, confidence = predictor.predict(compound)
        assert isinstance(score, float)
        assert isinstance(confidence, float)
        assert 0 <= score <= 1
        assert 0 <= confidence <= 1
        scores.append(score)
        confidences.append(confidence)

        # Test boolean prediction
        is_nootropic = predictor.predict_binary(compound, threshold=0.5)
        assert isinstance(is_nootropic, bool)

    # Validate score distribution
    scores = np.array(scores)
    confidences = np.array(confidences)
    assert np.all(scores[:3] > scores[3:])  # Strong nootropics score higher
    assert np.all(scores[3:5] > scores[5:])  # Mild nootropics score higher than non-nootropics
    assert np.mean(confidences) > 0.5  # Average confidence should be reasonable


def test_mechanism_prediction(
    predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test mechanism prediction."""
    # Train predictor
    predictor.train(
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
        mechanisms = predictor.predict_mechanisms(compound)
        assert isinstance(mechanisms, Dict)
        assert all(m in predictor.mechanism_types for m in mechanisms)
        assert all(isinstance(v, float) for v in mechanisms.values())
        assert all(0 <= v <= 1 for v in mechanisms.values())
        all_scores.append(list(mechanisms.values()))

        # Test mechanism thresholds
        primary = predictor.get_primary_mechanisms(compound, threshold=0.5)
        assert isinstance(primary, List)
        assert all(m in predictor.mechanism_types for m in primary)
        assert len(primary) <= len(mechanisms)

    # Analyze mechanism score distributions
    score_matrix = np.array(all_scores)
    assert score_matrix.shape == (len(test_compounds), len(predictor.mechanism_types))
    assert np.all(np.mean(score_matrix[:3], axis=1) > np.mean(score_matrix[5:], axis=1))


def test_cognitive_effects(
    predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test cognitive effects prediction."""
    compound = test_compounds[0]
    effect_predictions = predictor._predict_effects(compound)

    # Check prediction structure
    for domain, effects in effect_predictions.items():
        assert domain in predictor.cognitive_domains
        for effect, pred in effects.items():
            assert effect in predictor.cognitive_domains[domain]
            assert "score" in pred
            assert "confidence" in pred
            assert isinstance(pred["score"], float)
            assert isinstance(pred["confidence"], float)
            assert 0 <= pred["confidence"] <= 1


def test_side_effect_prediction(
    predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test side effects prediction."""
    compound = test_compounds[0]
    side_effect_predictions = predictor._predict_side_effects(compound)

    # Check prediction structure
    for category, effects in side_effect_predictions.items():
        assert category in predictor.side_effects
        for effect, pred in effects.items():
            assert effect in predictor.side_effects[category]
            assert "risk" in pred
            assert "confidence" in pred
            assert isinstance(pred["risk"], float)
            assert isinstance(pred["confidence"], float)
            assert 0 <= pred["confidence"] <= 1


def test_enhanced_prediction(
    predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test enhanced prediction pipeline."""
    compound = test_compounds[0]
    result = predictor.predict(compound)

    # Check result structure
    assert result.value in NootropicMechanism
    assert isinstance(result.confidence, float)
    assert 0 <= result.confidence <= 1

    # Check supporting data
    assert "effects" in result.supporting_data
    assert "side_effects" in result.supporting_data
    assert "literature_evidence" in result.supporting_data
    assert "safety_profile" in result.supporting_data

    # Check prediction history
    assert not predictor.prediction_history.empty
    assert "compound_name" in predictor.prediction_history.columns
    assert "mechanism" in predictor.prediction_history.columns
    assert "mechanism_confidence" in predictor.prediction_history.columns
    assert "timestamp" in predictor.prediction_history.columns


def test_model_persistence(
    predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
    tmp_path: Path,
):
    """Test model persistence."""
    # Train predictor
    predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Save models
    save_dir = tmp_path / "saved_models"
    predictor.save_models(save_dir)

    # Create new predictor and load models
    new_predictor = NootropicPredictorEnhanced(
        model_dir=str(save_dir),
        cache_dir=str(tmp_path / "new_cache"),
    )

    # Compare predictions
    for compound in test_compounds:
        result1 = predictor.predict(compound)
        result2 = new_predictor.predict(compound)

        assert result1.value == result2.value
        assert abs(result1.confidence - result2.confidence) < 1e-6

        mech1 = predictor.predict_mechanisms(compound)
        mech2 = new_predictor.predict_mechanisms(compound)
        assert set(mech1.keys()) == set(mech2.keys())
        assert all(abs(mech1[k] - mech2[k]) < 1e-6 for k in mech1)


def test_error_handling(predictor: NootropicPredictorEnhanced):
    """Test error handling."""
    # Test with invalid compound
    invalid_compound = CompoundData(
        name="Invalid",
        smiles="INVALID",
        cas_number="000-00-0",
    )
    result = predictor.predict(invalid_compound)

    assert result.value == NootropicMechanism.UNKNOWN
    assert result.confidence == 0.0
    assert "error" in result.supporting_data

    # Test with untrained model
    valid_compound = CompoundData(
        name="Valid",
        smiles="CCO",
        cas_number="64-17-5",
    )

    with pytest.raises(RuntimeError):
        predictor.predict_mechanisms(valid_compound)


def test_prediction_metrics(
    predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test prediction metrics."""
    # Train predictor
    metrics = predictor.train(
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

    # Check binary metrics
    binary_metrics = predictor.get_binary_metrics(threshold=0.5)
    assert "accuracy" in binary_metrics
    assert "precision" in binary_metrics
    assert "recall" in binary_metrics
    assert "f1" in binary_metrics
    assert all(isinstance(v, float) for v in binary_metrics.values())
    assert all(0 <= v <= 1 for v in binary_metrics.values())


def test_prediction_history(
    predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test prediction history tracking."""
    # Make predictions
    for compound in test_compounds:
        predictor.predict(compound)

    # Check history
    history = predictor.prediction_history
    assert len(history) > 0
    assert all(
        col in history.columns
        for col in [
            "compound_name",
            "mechanism",
            "mechanism_confidence",
            "timestamp",
        ]
    )

    # Check timestamps are ordered
    timestamps = pd.to_datetime(history["timestamp"])
    assert (timestamps.diff()[1:] >= pd.Timedelta(0)).all()


def test_safety_analysis(
    predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test safety analysis."""
    # Train predictor
    predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test safety analysis
    for compound in test_compounds:
        analysis = predictor.analyze_safety(compound)
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


def test_uncertainty_estimation(
    predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test uncertainty estimation in predictions."""
    # Train predictor
    predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test uncertainty metrics
    for compound in test_compounds:
        mechanism, confidence = predictor._predict_mechanism(compound)

        # Check confidence calculation components
        assert hasattr(predictor.mechanism_ensemble.models["mechanism"], "estimators_")
        estimators = predictor.mechanism_ensemble.models["mechanism"].estimators_

        # Get individual predictions
        predictions = []
        probabilities = []
        for _, model in estimators:
            pred = model.predict([compound])[0]
            prob = model.predict_proba([compound])[0]
            predictions.append(pred)
            probabilities.append(prob)

        # Verify disagreement calculation
        predictions = np.array(predictions)
        agreement = np.sum(predictions == mechanism) / len(predictions)
        assert 0 <= agreement <= 1

        # Verify probability variance
        probabilities = np.array(probabilities)
        variance = np.var(probabilities, axis=0)
        assert np.all(variance >= 0)

        # Verify confidence bounds
        assert 0 <= confidence <= 1


def test_confidence_calibration(
    predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test confidence calibration."""
    # Train predictor
    predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test calibration
    memory_labels = [NootropicMechanism.MEMORY_ENHANCEMENT] * 3
    focus_labels = [NootropicMechanism.FOCUS_IMPROVEMENT] * 2
    mechanism_labels = memory_labels + focus_labels
    effect_data = {
        "memory": {
            "working_memory": [0.9, 0.8, 0.7, 0.4, 0.3],
            "long_term_memory": [0.8, 0.7, 0.6, 0.3, 0.2],
        }
    }
    side_effect_data = {
        "physical": {
            "headache": [0.1, 0.2, 0.3, 0.4, 0.5],
            "nausea": [0.2, 0.3, 0.4, 0.5, 0.6],
        }
    }

    # Calibrate models
    metrics = predictor.calibrate_models(
        test_compounds[:5],
        mechanism_labels,
        effect_data,
        side_effect_data,
    )

    # Check calibration metrics
    assert "mechanism_brier_score" in metrics
    assert metrics["mechanism_brier_score"] >= 0
    assert metrics["mechanism_brier_score"] <= 1

    # Verify calibrated predictions
    for compound in test_compounds:
        result = predictor.predict(compound)
        assert 0 <= result.confidence <= 1

        # Check effect confidences
        for domain, effects in result.supporting_data["effects"].items():
            for effect, data in effects.items():
                assert 0 <= data["confidence"] <= 1

        # Check side effect confidences
        for category, effects in result.supporting_data["side_effects"].items():
            for effect, data in effects.items():
                assert 0 <= data["confidence"] <= 1


def test_feature_importance(
    predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test feature importance tracking."""
    # Train predictor
    predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Get feature importances
    importances = predictor.get_feature_importance()

    # Check mechanism importances
    assert "mechanism" in importances
    mechanism_importances = importances["mechanism"]
    for model_name, importance in mechanism_importances.items():
        assert isinstance(importance, np.ndarray)
        assert len(importance) > 0
        assert np.all(importance >= 0)
        assert np.isclose(np.sum(importance), 1.0, rtol=1e-5)

    # Check effect importances
    for domain, effects in predictor.cognitive_domains.items():
        for effect in effects:
            key = f"effect_{domain}_{effect}"
            if key in importances:
                effect_importance = importances[key]["importance"]
                assert isinstance(effect_importance, np.ndarray)
                assert len(effect_importance) > 0
                assert np.all(effect_importance >= 0)
                assert np.isclose(np.sum(effect_importance), 1.0, rtol=1e-5)

    # Check side effect importances
    for category, effects in predictor.side_effects.items():
        for effect in effects:
            key = f"side_effect_{category}_{effect}"
            if key in importances:
                effect_importance = importances[key]["importance"]
                assert isinstance(effect_importance, np.ndarray)
                assert len(effect_importance) > 0
                assert np.all(effect_importance >= 0)
                assert np.isclose(np.sum(effect_importance), 1.0, rtol=1e-5)


def test_bbb_integration(
    predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test BBB prediction integration."""
    # Train predictor
    predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test BBB integration
    for compound in test_compounds:
        result = predictor.predict(compound)

        # Check BBB data in result
        assert "bbb_prediction" in result.supporting_data
        bbb_data = result.supporting_data["bbb_prediction"]
        assert "value" in bbb_data
        assert "confidence" in bbb_data
        assert "supporting_data" in bbb_data
        assert isinstance(bbb_data["value"], str)
        assert isinstance(bbb_data["confidence"], float)
        assert isinstance(bbb_data["supporting_data"], dict)

        # Check BBB impact on predictions
        for domain, effects in result.supporting_data["effects"].items():
            for effect, data in effects.items():
                # Effect scores should be scaled by BBB permeability
                assert data["score"] <= data["bbb_factor"]
                assert 0 <= data["confidence"] <= 1

        # Check BBB impact on side effects
        for category, effects in result.supporting_data["side_effects"].items():
            for effect, data in effects.items():
                # CNS effects should be scaled by BBB permeability
                if category in ["cognitive", "neurological"]:
                    assert data["risk"] <= data["bbb_factor"]
                assert 0 <= data["confidence"] <= 1


def test_literature_analysis(
    predictor: NootropicPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test literature analysis."""
    # Train predictor
    predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test literature analysis
    for compound in test_compounds:
        analysis = predictor.analyze_literature(compound)
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
