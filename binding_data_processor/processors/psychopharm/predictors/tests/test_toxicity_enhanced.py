"""Tests for enhanced toxicity prediction and analysis.

This module provides comprehensive tests for:
1. Basic toxicity prediction 
2. Mechanism prediction
3. Enhanced prediction with evidence
4. Organ toxicity prediction
5. Safety concern prediction
6. Model persistence
7. Error handling
8. Metrics calculation
9. BBB integration
10. Prediction history tracking
"""

import logging
import pytest
from pathlib import Path
from typing import Dict, List, Set

import numpy as np
import pandas as pd
from rdkit import Chem

from .....models.core import CompoundData
from .....models.psychopharm import ToxicityClass
from ..toxicity_enhanced import ToxicityPredictorEnhanced


@pytest.fixture
def test_compounds() -> List[CompoundData]:
    """Create test compounds with known properties."""
    # Known toxic compounds
    acetaminophen = CompoundData(
        name="Acetaminophen",
        smiles="CC(=O)NC1=CC=C(O)C=C1",
        cas_number="103-90-2",
    )
    ethanol = CompoundData(
        name="Ethanol",
        smiles="CCO",
        cas_number="64-17-5",
    )
    caffeine = CompoundData(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )

    # Compounds with liver toxicity
    isoniazid = CompoundData(
        name="Isoniazid",
        smiles="NNC(=O)C1=CN=CC=C1",
        cas_number="54-85-3",
    )
    rifampicin = CompoundData(
        name="Rifampicin",
        smiles=(
            "CC1COC2=C(C3=C(C=C2C(=O)C(=C/C4=CC=C(NC(=O)C)C(O)=C4)/C)" "C(=O)C5=C(O)C=CC=C5O3)" "C6=C(O)C=C(N(C)C)C=C16"
        ),
        cas_number="13292-46-1",
    )

    # Compounds with kidney toxicity
    cisplatin = CompoundData(
        name="Cisplatin",
        smiles="N.N.Cl[Pt]Cl",
        cas_number="15663-27-1",
    )
    gentamicin = CompoundData(
        name="Gentamicin",
        smiles="CC(C1CCC(N)C(OC2C(N)CC(N)C(OC3OCC(C)(O)C(NC)C3O)C2O)O1)NC",
        cas_number="1403-66-3",
    )

    # Compounds with neurotoxicity
    lead = CompoundData(
        name="Lead",
        smiles="[Pb]",
        cas_number="7439-92-1",
    )
    mercury = CompoundData(
        name="Mercury",
        smiles="[Hg]",
        cas_number="7439-97-6",
    )

    # Non-toxic compounds
    glucose = CompoundData(
        name="Glucose",
        smiles="C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O",
        cas_number="50-99-7",
    )
    glycine = CompoundData(
        name="Glycine",
        smiles="NCC(=O)O",
        cas_number="56-40-6",
    )

    # Calculate molecular weights for validation
    compounds = [
        acetaminophen,
        ethanol,
        caffeine,
        isoniazid,
        rifampicin,
        cisplatin,
        gentamicin,
        lead,
        mercury,
        glucose,
        glycine,
    ]
    for compound in compounds:
        mol = Chem.MolFromSmiles(compound.smiles)
        if mol:
            compound.molecular_weight = Chem.Descriptors.ExactMolWt(mol)

    return compounds


@pytest.fixture
def test_toxicity_mechanisms() -> Dict[str, Set[str]]:
    """Create test toxicity mechanisms."""
    return {
        "cellular": {"oxidative_stress", "mitochondrial_toxicity", "dna_damage"},
        "molecular": {"reactive_metabolites", "protein_binding", "enzyme_inhibition"},
        "systemic": {"inflammation", "apoptosis", "organ_failure"},
    }


@pytest.fixture
def test_organ_toxicity() -> Dict[str, Set[str]]:
    """Create test organ toxicity categories."""
    return {
        "liver": {"hepatocellular_damage", "cholestasis", "steatosis"},
        "kidney": {"tubular_damage", "glomerular_damage", "crystal_formation"},
        "brain": {"neurotoxicity", "cognitive_impairment", "seizure_risk"},
        "heart": {"arrhythmia", "contractility_changes", "qt_prolongation"},
    }


@pytest.fixture
def test_safety_concerns() -> Dict[str, Set[str]]:
    """Create test safety concerns."""
    return {
        "acute": {"ld50", "immediate_effects", "overdose_risk"},
        "chronic": {"carcinogenicity", "reproductive_toxicity", "organ_damage"},
        "special": {"drug_interactions", "contraindications", "genetic_factors"},
    }


@pytest.fixture
def predictor(
    tmp_path: Path,
    test_toxicity_mechanisms: Dict[str, Set[str]],
    test_organ_toxicity: Dict[str, Set[str]],
    test_safety_concerns: Dict[str, Set[str]],
) -> ToxicityPredictorEnhanced:
    """Create enhanced toxicity predictor."""
    return ToxicityPredictorEnhanced(
        model_dir=str(tmp_path / "models"),
        cache_dir=str(tmp_path / "cache"),
        log_level=logging.DEBUG,
        toxicity_mechanisms=test_toxicity_mechanisms,
        organ_toxicity=test_organ_toxicity,
        safety_concerns=test_safety_concerns,
    )


def test_predictor_initialization(
    predictor: ToxicityPredictorEnhanced,
    test_toxicity_mechanisms: Dict[str, Set[str]],
    test_organ_toxicity: Dict[str, Set[str]],
    test_safety_concerns: Dict[str, Set[str]],
):
    """Test predictor initialization."""
    assert predictor.ml_pipeline is not None
    assert predictor.feature_types == ["fingerprints", "descriptors", "enhanced"]
    assert predictor.toxicity_mechanisms == test_toxicity_mechanisms
    assert predictor.organ_toxicity == test_organ_toxicity
    assert predictor.safety_concerns == test_safety_concerns

    # Check ensemble initialization
    assert predictor.mechanism_ensemble is not None
    for category, mechanisms in test_toxicity_mechanisms.items():
        assert category in predictor.mechanism_ensembles
        for mechanism in mechanisms:
            assert mechanism in predictor.mechanism_ensembles[category]

    for organ, toxicities in test_organ_toxicity.items():
        assert organ in predictor.organ_ensembles
        for toxicity in toxicities:
            assert toxicity in predictor.organ_ensembles[organ]

    for category, concerns in test_safety_concerns.items():
        assert category in predictor.safety_ensembles
        for concern in concerns:
            assert concern in predictor.safety_ensembles[category]


def test_basic_prediction(
    predictor: ToxicityPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test basic toxicity prediction."""
    # Train predictor with both boolean and float labels
    predictor.train(
        test_compounds[:3],  # Known toxic compounds
        [1.0, 0.95, 0.9],
        test_compounds[3:8],  # Organ-specific toxicity compounds
        [0.7, 0.6, 0.8, 0.7, 0.8],
        test_compounds[8:],  # Non-toxic compounds
        [0.4, 0.3, 0.2],
    )

    # Test predictions
    scores = []
    confidences = []
    for compound in test_compounds:
        # Test float score prediction
        result = predictor.predict(compound)
        assert isinstance(result.value, ToxicityClass)
        assert isinstance(result.confidence, float)
        assert 0 <= result.confidence <= 1
        scores.append(float(result.value.value))
        confidences.append(result.confidence)

    # Validate score distribution
    scores = np.array(scores)
    confidences = np.array(confidences)
    assert np.all(scores[:3] > scores[8:])  # Known toxic compounds score higher
    assert np.all(scores[3:8] > scores[8:])  # Organ-toxic compounds score higher than non-toxic
    assert np.mean(confidences) > 0.5  # Average confidence should be reasonable


def test_mechanism_prediction(
    predictor: ToxicityPredictorEnhanced,
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
        result = predictor.predict(compound)
        mechanisms = result.supporting_data["mechanisms"]
        assert isinstance(mechanisms, Dict)
        for category, category_mechanisms in mechanisms.items():
            assert category in predictor.toxicity_mechanisms
            for mechanism, data in category_mechanisms.items():
                assert mechanism in predictor.toxicity_mechanisms[category]
                assert "score" in data
                assert "confidence" in data
                assert isinstance(data["score"], float)
                assert isinstance(data["confidence"], float)
                assert 0 <= data["score"] <= 1
                assert 0 <= data["confidence"] <= 1
                all_scores.append(data["score"])

    # Analyze mechanism score distributions
    score_matrix = np.array(all_scores)
    assert len(score_matrix) > 0
    assert np.all(score_matrix >= 0)
    assert np.all(score_matrix <= 1)


def test_organ_toxicity_prediction(
    predictor: ToxicityPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test organ toxicity prediction."""
    # Train predictor
    predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test organ toxicity predictions
    for compound in test_compounds:
        result = predictor.predict(compound)
        organs = result.supporting_data["organs"]
        assert isinstance(organs, Dict)
        for organ, toxicities in organs.items():
            assert organ in predictor.organ_toxicity
            for toxicity, data in toxicities.items():
                assert toxicity in predictor.organ_toxicity[organ]
                assert "severity" in data
                assert "confidence" in data
                assert isinstance(data["severity"], float)
                assert isinstance(data["confidence"], float)
                assert 0 <= data["severity"] <= 1
                assert 0 <= data["confidence"] <= 1


def test_safety_concern_prediction(
    predictor: ToxicityPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test safety concern prediction."""
    # Train predictor
    predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test safety concern predictions
    for compound in test_compounds:
        result = predictor.predict(compound)
        concerns = result.supporting_data["concerns"]
        assert isinstance(concerns, Dict)
        for category, category_concerns in concerns.items():
            assert category in predictor.safety_concerns
            for concern, data in category_concerns.items():
                assert concern in predictor.safety_concerns[category]
                assert "risk" in data
                assert "confidence" in data
                assert isinstance(data["risk"], float)
                assert isinstance(data["confidence"], float)
                assert 0 <= data["risk"] <= 1
                assert 0 <= data["confidence"] <= 1


def test_model_persistence(
    predictor: ToxicityPredictorEnhanced,
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
    new_predictor = ToxicityPredictorEnhanced(
        model_dir=str(save_dir),
        cache_dir=str(tmp_path / "new_cache"),
    )

    # Compare predictions
    for compound in test_compounds:
        result1 = predictor.predict(compound)
        result2 = new_predictor.predict(compound)

        assert result1.value == result2.value
        assert abs(result1.confidence - result2.confidence) < 1e-6

        # Compare mechanisms
        mech1 = result1.supporting_data["mechanisms"]
        mech2 = result2.supporting_data["mechanisms"]
        assert set(mech1.keys()) == set(mech2.keys())
        for category in mech1:
            assert set(mech1[category].keys()) == set(mech2[category].keys())


def test_error_handling(predictor: ToxicityPredictorEnhanced):
    """Test error handling."""
    # Test with invalid compound
    invalid_compound = CompoundData(
        name="Invalid",
        smiles="INVALID",
        cas_number="000-00-0",
    )
    result = predictor.predict(invalid_compound)

    assert result.value == ToxicityClass.UNKNOWN
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
    predictor: ToxicityPredictorEnhanced,
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
    assert "mechanism_metrics" in metrics
    assert "organ_metrics" in metrics
    assert "safety_metrics" in metrics
    assert all(isinstance(v, Dict) for v in metrics.values())
    assert all(all(isinstance(v, float) for v in m.values()) for m in metrics.values())
    assert all(all(v >= 0 for v in m.values()) for m in metrics.values())


def test_prediction_history(
    predictor: ToxicityPredictorEnhanced,
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


def test_bbb_integration(
    predictor: ToxicityPredictorEnhanced,
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
        for category, mechanisms in result.supporting_data["mechanisms"].items():
            for mechanism, data in mechanisms.items():
                # CNS mechanisms should be scaled by BBB permeability
                if category in ["neurological", "cognitive"]:
                    assert data["score"] <= data["bbb_factor"]
                assert 0 <= data["confidence"] <= 1

        # Check BBB impact on organ toxicity
        for organ, toxicities in result.supporting_data["organs"].items():
            for toxicity, data in toxicities.items():
                # Brain toxicity should be scaled by BBB permeability
                if organ == "brain":
                    assert data["severity"] <= data["bbb_factor"]
                assert 0 <= data["confidence"] <= 1


def test_feature_importance(
    predictor: ToxicityPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test feature importance analysis."""
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

    # Check structure
    assert isinstance(importances, Dict)
    assert len(importances) > 0

    # Check mechanism importances
    if "mechanism" in importances:
        for model_name, importance in importances["mechanism"].items():
            assert isinstance(importance, np.ndarray)
            assert len(importance) > 0
            assert np.all(importance >= 0)
            assert np.isclose(np.sum(importance), 1.0, rtol=1e-5)

    # Check ensemble importances
    for key, data in importances.items():
        if key != "mechanism":
            assert isinstance(data, Dict)
            assert "importance" in data
            assert isinstance(data["importance"], np.ndarray)
            assert len(data["importance"]) > 0
            assert np.all(data["importance"] >= 0)
            assert np.isclose(np.sum(data["importance"]), 1.0, rtol=1e-5)


def test_model_calibration(
    predictor: ToxicityPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test model calibration."""
    # Train predictor
    predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Create calibration data
    mechanism_data = {
        "cellular": {
            "oxidative_stress": [0.8, 0.7, 0.6, 0.4, 0.3],
            "mitochondrial_toxicity": [0.9, 0.8, 0.5, 0.3, 0.2],
        }
    }
    organ_data = {
        "liver": {
            "hepatocellular_damage": [0.9, 0.8, 0.6, 0.3, 0.2],
            "cholestasis": [0.7, 0.6, 0.5, 0.4, 0.3],
        }
    }
    safety_data = {
        "acute": {
            "ld50": [0.9, 0.8, 0.7, 0.3, 0.2],
            "overdose_risk": [0.8, 0.7, 0.6, 0.4, 0.3],
        }
    }

    # Calibrate models
    metrics = predictor.calibrate_models(
        test_compounds[:5],
        [ToxicityClass.HIGH] * 3 + [ToxicityClass.MEDIUM] * 2,
        mechanism_data,
        organ_data,
        safety_data,
    )

    # Check calibration metrics
    assert "mechanism_brier_score" in metrics
    assert metrics["mechanism_brier_score"] >= 0
    assert metrics["mechanism_brier_score"] <= 1

    # Verify calibrated predictions
    for compound in test_compounds:
        result = predictor.predict(compound)
        assert 0 <= result.confidence <= 1

        # Check mechanism confidences
        for category, mechanisms in result.supporting_data["mechanisms"].items():
            for mechanism, data in mechanisms.items():
                assert 0 <= data["confidence"] <= 1

        # Check organ toxicity confidences
        for organ, toxicities in result.supporting_data["organs"].items():
            for toxicity, data in toxicities.items():
                assert 0 <= data["confidence"] <= 1

        # Check safety concern confidences
        for category, concerns in result.supporting_data["concerns"].items():
            for concern, data in concerns.items():
                assert 0 <= data["confidence"] <= 1
