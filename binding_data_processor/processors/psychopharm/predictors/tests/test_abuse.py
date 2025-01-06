"""Tests for abuse potential prediction.

This module provides comprehensive tests for:
1. Basic abuse potential prediction
2. Mechanism prediction
3. Enhanced prediction with evidence
4. Model persistence
5. Error handling
6. Metrics calculation
7. Risk analysis
8. Regulatory analysis
"""

import logging
import pytest
from pathlib import Path
from typing import List, Dict

import numpy as np
from rdkit import Chem

from .....models.core import CompoundData
from ..abuse import AbusePotentialPredictor
from ..abuse_enhanced import AbusePotentialPredictorEnhanced


@pytest.fixture
def test_compounds() -> List[CompoundData]:
    """Create test compounds with known properties."""
    # Known high abuse potential
    morphine = CompoundData(
        name="Morphine",
        smiles="CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
        cas_number="57-27-2",
    )
    cocaine = CompoundData(
        name="Cocaine",
        smiles="CN1[C@H]2CC[C@H]1[C@@H](C(=O)OC)[C@H](OC(=O)c1ccccc1)C2",
        cas_number="50-36-2",
    )
    methamphetamine = CompoundData(
        name="Methamphetamine",
        smiles="CN[C@H](C)Cc1ccccc1",
        cas_number="537-46-2",
    )

    # Compounds with moderate abuse potential
    ketamine = CompoundData(
        name="Ketamine",
        smiles="CN1C(=O)CCC(N(C)C)(c2ccccc2)O",
        cas_number="6740-88-1",
    )
    tramadol = CompoundData(
        name="Tramadol",
        smiles="CC(=O)CC(c1ccccc1)N(C)C",
        cas_number="27203-92-5",
    )

    # Compounds with low/no abuse potential
    ibuprofen = CompoundData(
        name="Ibuprofen",
        smiles="CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O",
        cas_number="15687-27-1",
    )
    acetaminophen = CompoundData(
        name="Acetaminophen",
        smiles="CC(=O)Nc1ccc(O)cc1",
        cas_number="103-90-2",
    )

    # Calculate molecular properties for validation
    for compound in [morphine, cocaine, methamphetamine, ketamine, tramadol, ibuprofen, acetaminophen]:
        mol = Chem.MolFromSmiles(compound.smiles)
        if mol:
            compound.molecular_weight = Chem.Descriptors.ExactMolWt(mol)
            compound.logp = Chem.Crippen.MolLogP(mol)
            compound.hbd = Chem.rdMolDescriptors.CalcNumHBD(mol)
            compound.hba = Chem.rdMolDescriptors.CalcNumHBA(mol)
            compound.tpsa = Chem.Descriptors.TPSA(mol)
            compound.rotatable_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)

    return [morphine, cocaine, methamphetamine, ketamine, tramadol, ibuprofen, acetaminophen]


@pytest.fixture
def abuse_predictor(tmp_path: Path) -> AbusePotentialPredictor:
    """Create abuse potential predictor."""
    return AbusePotentialPredictor(
        model_dir=str(tmp_path / "models"),
        cache_dir=str(tmp_path / "cache"),
        log_level=logging.DEBUG,
    )


@pytest.fixture
def enhanced_predictor(tmp_path: Path) -> AbusePotentialPredictorEnhanced:
    """Create enhanced abuse potential predictor."""
    return AbusePotentialPredictorEnhanced(
        model_dir=str(tmp_path / "models"),
        cache_dir=str(tmp_path / "cache"),
        log_level=logging.DEBUG,
    )


def test_predictor_initialization(abuse_predictor: AbusePotentialPredictor):
    """Test predictor initialization."""
    assert abuse_predictor.ml_pipeline is not None
    assert abuse_predictor.feature_types == ["morgan", "maccs", "rdkit"]
    assert abuse_predictor.mechanism_types == [
        "DOPAMINE_RELEASE",
        "SEROTONIN_RELEASE",
        "NOREPINEPHRINE_RELEASE",
        "OPIOID_AGONISM",
        "GABA_MODULATION",
        "NMDA_ANTAGONISM",
        "CANNABINOID_AGONISM",
        "REWARD_PATHWAY",
        "TOLERANCE_DEVELOPMENT",
        "WITHDRAWAL_SEVERITY",
        "REINFORCEMENT",
    ]


def test_basic_prediction(
    abuse_predictor: AbusePotentialPredictor,
    test_compounds: List[CompoundData],
):
    """Test basic abuse potential prediction."""
    # Train predictor with both boolean and float labels
    abuse_predictor.train(
        test_compounds[:3],  # High abuse potential
        [1.0, 0.95, 0.9],
        test_compounds[3:5],  # Moderate abuse potential
        [0.6, 0.5],
        test_compounds[5:],  # Low/no abuse potential
        [0.1, 0.1],
    )

    # Test predictions
    scores = []
    confidences = []
    for compound in test_compounds:
        # Test float score prediction
        score, confidence = abuse_predictor.predict(compound)
        assert isinstance(score, float)
        assert isinstance(confidence, float)
        assert 0 <= score <= 1
        assert 0 <= confidence <= 1
        scores.append(score)
        confidences.append(confidence)

        # Test binary prediction
        is_high_risk = abuse_predictor.predict_binary(compound, threshold=0.5)
        assert isinstance(is_high_risk, bool)

    # Validate score distribution
    scores = np.array(scores)
    confidences = np.array(confidences)
    assert np.all(scores[:3] > scores[3:])  # High risk compounds score higher
    assert np.all(scores[3:5] > scores[5:])  # Moderate risk compounds score higher than low risk
    assert np.mean(confidences) > 0.5  # Average confidence should be reasonable


def test_mechanism_prediction(
    abuse_predictor: AbusePotentialPredictor,
    test_compounds: List[CompoundData],
):
    """Test mechanism prediction."""
    # Train predictor
    abuse_predictor.train(
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
        mechanisms = abuse_predictor.predict_mechanisms(compound)
        assert isinstance(mechanisms, Dict)
        assert all(m in abuse_predictor.mechanism_types for m in mechanisms)
        assert all(isinstance(v, float) for v in mechanisms.values())
        assert all(0 <= v <= 1 for v in mechanisms.values())
        all_scores.append(list(mechanisms.values()))

        # Test mechanism thresholds
        primary = abuse_predictor.get_primary_mechanisms(compound, threshold=0.5)
        assert isinstance(primary, List)
        assert all(m in abuse_predictor.mechanism_types for m in primary)
        assert len(primary) <= len(mechanisms)

    # Analyze mechanism score distributions
    score_matrix = np.array(all_scores)
    assert score_matrix.shape == (len(test_compounds), len(abuse_predictor.mechanism_types))
    assert np.all(np.mean(score_matrix[:3], axis=1) > np.mean(score_matrix[5:], axis=1))


def test_enhanced_prediction(
    enhanced_predictor: AbusePotentialPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test enhanced abuse potential prediction."""
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
        assert "abuse_potential" in prediction
        assert "confidence" in prediction
        assert "mechanisms" in prediction
        assert "literature_evidence" in prediction
        assert "regulatory_status" in prediction
        assert "risk_factors" in prediction
        assert "community_reports" in prediction
        assert isinstance(prediction["abuse_potential"], float)
        assert isinstance(prediction["confidence"], float)
        assert isinstance(prediction["mechanisms"], Dict)
        assert isinstance(prediction["literature_evidence"], List)
        assert isinstance(prediction["regulatory_status"], Dict)
        assert isinstance(prediction["risk_factors"], Dict)
        assert isinstance(prediction["community_reports"], List)

        # Check literature evidence structure
        for evidence in prediction["literature_evidence"]:
            assert "source" in evidence
            assert "title" in evidence
            assert "year" in evidence
            assert "relevance" in evidence
            assert "findings" in evidence

        # Check regulatory status structure
        assert "scheduling" in prediction["regulatory_status"]
        assert "controlled_status" in prediction["regulatory_status"]
        assert "restrictions" in prediction["regulatory_status"]
        assert "monitoring_requirements" in prediction["regulatory_status"]

        # Check risk factors structure
        assert "pharmacological" in prediction["risk_factors"]
        assert "behavioral" in prediction["risk_factors"]
        assert "social" in prediction["risk_factors"]
        assert "medical" in prediction["risk_factors"]


def test_model_persistence(
    abuse_predictor: AbusePotentialPredictor,
    test_compounds: List[CompoundData],
    tmp_path: Path,
):
    """Test model persistence."""
    # Train predictor
    abuse_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Save models
    save_dir = tmp_path / "saved_models"
    abuse_predictor.save_models(save_dir)

    # Create new predictor and load models
    new_predictor = AbusePotentialPredictor(
        model_dir=str(save_dir),
        cache_dir=str(tmp_path / "new_cache"),
    )

    # Compare predictions
    for compound in test_compounds:
        score1, conf1 = abuse_predictor.predict(compound)
        score2, conf2 = new_predictor.predict(compound)
        assert abs(score1 - score2) < 1e-6
        assert abs(conf1 - conf2) < 1e-6

        mech1 = abuse_predictor.predict_mechanisms(compound)
        mech2 = new_predictor.predict_mechanisms(compound)
        assert set(mech1.keys()) == set(mech2.keys())
        assert all(abs(mech1[k] - mech2[k]) < 1e-6 for k in mech1)


def test_error_handling(abuse_predictor: AbusePotentialPredictor):
    """Test error handling."""
    # Test with invalid compound
    invalid_compound = CompoundData(
        name="Invalid",
        smiles="INVALID",
        cas_number="000-00-0",
    )

    with pytest.raises(ValueError):
        abuse_predictor.predict(invalid_compound)

    with pytest.raises(ValueError):
        abuse_predictor.predict_mechanisms(invalid_compound)

    # Test with untrained model
    valid_compound = CompoundData(
        name="Valid",
        smiles="CCO",
        cas_number="64-17-5",
    )

    with pytest.raises(RuntimeError):
        abuse_predictor.predict(valid_compound)

    with pytest.raises(RuntimeError):
        abuse_predictor.predict_mechanisms(valid_compound)


def test_prediction_metrics(
    abuse_predictor: AbusePotentialPredictor,
    test_compounds: List[CompoundData],
):
    """Test prediction metrics."""
    # Train predictor
    metrics = abuse_predictor.train(
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
    binary_metrics = abuse_predictor.get_binary_metrics(threshold=0.5)
    assert "accuracy" in binary_metrics
    assert "precision" in binary_metrics
    assert "recall" in binary_metrics
    assert "f1" in binary_metrics
    assert all(isinstance(v, float) for v in binary_metrics.values())
    assert all(0 <= v <= 1 for v in binary_metrics.values())


def test_enhanced_metrics(
    enhanced_predictor: AbusePotentialPredictorEnhanced,
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
    assert "regulatory_metrics" in metrics
    assert "community_metrics" in metrics
    assert "combined_metrics" in metrics
    assert all(isinstance(v, Dict) for v in metrics.values())
    assert all(all(isinstance(v, float) for v in m.values()) for m in metrics.values())
    assert all(all(v >= 0 for v in m.values()) for m in metrics.values())


def test_risk_analysis(
    enhanced_predictor: AbusePotentialPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test risk analysis for abuse potential."""
    # Train predictor
    enhanced_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test risk analysis
    for compound in test_compounds:
        analysis = enhanced_predictor.analyze_risk(compound)
        assert isinstance(analysis, Dict)
        assert "abuse_score" in analysis
        assert "dependence_potential" in analysis
        assert "withdrawal_severity" in analysis
        assert "tolerance_profile" in analysis
        assert "contraindications" in analysis
        assert "drug_interactions" in analysis
        assert "overdose_risk" in analysis
        assert "social_factors" in analysis
        assert isinstance(analysis["abuse_score"], float)
        assert isinstance(analysis["dependence_potential"], float)
        assert isinstance(analysis["withdrawal_severity"], Dict)
        assert isinstance(analysis["tolerance_profile"], Dict)
        assert isinstance(analysis["contraindications"], List)
        assert isinstance(analysis["drug_interactions"], Dict)
        assert isinstance(analysis["overdose_risk"], Dict)
        assert isinstance(analysis["social_factors"], Dict)


def test_regulatory_analysis(
    enhanced_predictor: AbusePotentialPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test regulatory analysis for abuse potential."""
    # Train predictor
    enhanced_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1],
    )

    # Test regulatory analysis
    for compound in test_compounds:
        analysis = enhanced_predictor.analyze_regulatory_status(compound)
        assert isinstance(analysis, Dict)
        assert "scheduling_status" in analysis
        assert "control_measures" in analysis
        assert "prescribing_restrictions" in analysis
        assert "monitoring_requirements" in analysis
        assert "international_controls" in analysis
        assert "state_regulations" in analysis
        assert "reporting_requirements" in analysis
        assert "diversion_controls" in analysis
        assert isinstance(analysis["scheduling_status"], Dict)
        assert isinstance(analysis["control_measures"], List)
        assert isinstance(analysis["prescribing_restrictions"], Dict)
        assert isinstance(analysis["monitoring_requirements"], Dict)
        assert isinstance(analysis["international_controls"], Dict)
        assert isinstance(analysis["state_regulations"], Dict)
        assert isinstance(analysis["reporting_requirements"], List)
        assert isinstance(analysis["diversion_controls"], Dict)
