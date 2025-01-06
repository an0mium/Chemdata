"""Tests for toxicity prediction.

This module provides comprehensive tests for:
1. Basic toxicity prediction
2. Toxicity type prediction
3. Mechanism prediction
4. Enhanced prediction with evidence
5. Model persistence
6. Error handling
7. Metrics calculation
8. Risk analysis
9. Regulatory analysis
"""

import logging
import pytest
from pathlib import Path
from typing import List, Dict

import numpy as np
from rdkit import Chem

from .....models.core import CompoundData
from ..toxicity import ToxicityPredictor
from ..toxicity_enhanced import ToxicityPredictorEnhanced


@pytest.fixture
def test_compounds() -> List[CompoundData]:
    """Create test compounds with known toxicity profiles."""
    # Known high toxicity compounds
    strychnine = CompoundData(
        name="Strychnine",
        smiles="O=C1CC2OCC=C3CN4CCC56C(C(=O)CC5)C5(O)C(=O)OCC5C6C4C3C2C1",
        cas_number="57-24-9",
    )
    cyanide = CompoundData(
        name="Hydrogen Cyanide",
        smiles="C#N",
        cas_number="74-90-8",
    )
    tetrodotoxin = CompoundData(
        name="Tetrodotoxin",
        smiles="CC1=C(C(=O)N2CCCC2=O)NC(=O)NC1=O",
        cas_number="4368-28-9",
    )

    # Psychiatric drugs with toxicity concerns
    chlorpromazine = CompoundData(
        name="Chlorpromazine",
        smiles="CN(C)CCCN1C2=CC=CC=C2SC2=C1C=C(C=C2)Cl",
        cas_number="50-53-3",
    )
    amitriptyline = CompoundData(
        name="Amitriptyline",
        smiles="CN(C)CCC=C1C2=CC=CC=C2CCC2=CC=CC=C12",
        cas_number="50-48-6",
    )

    # Known safe compounds
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
    melatonin = CompoundData(
        name="Melatonin",
        smiles="CC(=O)NCCC1=CNc2c1cc(OC)cc2",
        cas_number="73-31-4",
    )

    # Calculate molecular properties for validation
    for compound in [strychnine, cyanide, tetrodotoxin, chlorpromazine, amitriptyline, glucose, glycine, melatonin]:
        mol = Chem.MolFromSmiles(compound.smiles)
        if mol:
            compound.molecular_weight = Chem.Descriptors.ExactMolWt(mol)
            compound.logp = Chem.Crippen.MolLogP(mol)
            compound.hbd = Chem.rdMolDescriptors.CalcNumHBD(mol)
            compound.hba = Chem.rdMolDescriptors.CalcNumHBA(mol)
            compound.tpsa = Chem.Descriptors.TPSA(mol)
            compound.rotatable_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)

    return [strychnine, cyanide, tetrodotoxin, chlorpromazine, amitriptyline, glucose, glycine, melatonin]


@pytest.fixture
def toxicity_predictor(tmp_path: Path) -> ToxicityPredictor:
    """Create toxicity predictor."""
    return ToxicityPredictor(
        model_dir=str(tmp_path / "models"),
        cache_dir=str(tmp_path / "cache"),
        log_level=logging.DEBUG,
    )


@pytest.fixture
def enhanced_predictor(tmp_path: Path) -> ToxicityPredictorEnhanced:
    """Create enhanced toxicity predictor."""
    return ToxicityPredictorEnhanced(
        model_dir=str(tmp_path / "models"),
        cache_dir=str(tmp_path / "cache"),
        log_level=logging.DEBUG,
    )


def test_predictor_initialization(toxicity_predictor: ToxicityPredictor):
    """Test predictor initialization."""
    assert toxicity_predictor.ml_pipeline is not None
    assert toxicity_predictor.feature_types == ["morgan", "maccs", "rdkit"]
    assert toxicity_predictor.toxicity_types == [
        "ACUTE_ORAL",
        "ACUTE_DERMAL",
        "ACUTE_INHALATION",
        "SKIN_IRRITATION",
        "EYE_IRRITATION",
        "SKIN_SENSITIZATION",
        "MUTAGENICITY",
        "CARCINOGENICITY",
        "REPRODUCTIVE_TOXICITY",
        "ORGAN_TOXICITY",
    ]
    assert toxicity_predictor.mechanism_types == [
        "CYTOTOXICITY",
        "HEPATOTOXICITY",
        "CARDIOTOXICITY",
        "NEUROTOXICITY",
        "NEPHROTOXICITY",
        "IMMUNOTOXICITY",
        "REPRODUCTIVE_TOXICITY",
        "GENOTOXICITY",
        "CARCINOGENICITY",
        "RESPIRATORY_TOXICITY",
        "ENDOCRINE_DISRUPTION",
    ]


def test_basic_prediction(
    toxicity_predictor: ToxicityPredictor,
    test_compounds: List[CompoundData],
):
    """Test basic toxicity prediction."""
    # Train predictor with both boolean and float labels
    toxicity_predictor.train(
        test_compounds[:3],  # High toxicity
        [1.0, 0.95, 0.9],
        test_compounds[3:5],  # Moderate toxicity
        [0.6, 0.5],
        test_compounds[5:],  # Low toxicity
        [0.1, 0.1, 0.1],
    )

    # Test predictions
    scores = []
    confidences = []
    for compound in test_compounds:
        # Test float score prediction
        score, confidence = toxicity_predictor.predict(compound)
        assert isinstance(score, float)
        assert isinstance(confidence, float)
        assert 0 <= score <= 1
        assert 0 <= confidence <= 1
        scores.append(score)
        confidences.append(confidence)

        # Test binary prediction
        is_toxic = toxicity_predictor.predict_binary(compound, threshold=0.5)
        assert isinstance(is_toxic, bool)

    # Validate score distribution
    scores = np.array(scores)
    confidences = np.array(confidences)
    assert np.all(scores[:3] > scores[3:])  # High toxicity compounds score higher
    assert np.all(scores[3:5] > scores[5:])  # Moderate toxicity compounds score higher than low toxicity
    assert np.mean(confidences) > 0.5  # Average confidence should be reasonable


def test_toxicity_type_prediction(
    toxicity_predictor: ToxicityPredictor,
    test_compounds: List[CompoundData],
):
    """Test toxicity type prediction."""
    # Train predictor
    toxicity_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1, 0.1],
    )

    # Test toxicity type predictions
    all_scores = []
    for compound in test_compounds:
        toxicities = toxicity_predictor.predict_toxicity_types(compound)
        assert isinstance(toxicities, Dict)
        assert all(t in toxicity_predictor.toxicity_types for t in toxicities)
        assert all(isinstance(v, float) for v in toxicities.values())
        assert all(0 <= v <= 1 for v in toxicities.values())
        all_scores.append(list(toxicities.values()))

        # Test toxicity type thresholds
        primary = toxicity_predictor.get_primary_toxicity_types(compound, threshold=0.5)
        assert isinstance(primary, List)
        assert all(t in toxicity_predictor.toxicity_types for t in primary)
        assert len(primary) <= len(toxicities)

    # Analyze toxicity type score distributions
    score_matrix = np.array(all_scores)
    assert score_matrix.shape == (len(test_compounds), len(toxicity_predictor.toxicity_types))
    assert np.all(np.mean(score_matrix[:3], axis=1) > np.mean(score_matrix[5:], axis=1))


def test_mechanism_prediction(
    toxicity_predictor: ToxicityPredictor,
    test_compounds: List[CompoundData],
):
    """Test mechanism prediction."""
    # Train predictor
    toxicity_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1, 0.1],
    )

    # Test mechanism predictions
    all_scores = []
    for compound in test_compounds:
        mechanisms = toxicity_predictor.predict_mechanisms(compound)
        assert isinstance(mechanisms, Dict)
        assert all(m in toxicity_predictor.mechanism_types for m in mechanisms)
        assert all(isinstance(v, float) for v in mechanisms.values())
        assert all(0 <= v <= 1 for v in mechanisms.values())
        all_scores.append(list(mechanisms.values()))

        # Test mechanism thresholds
        primary = toxicity_predictor.get_primary_mechanisms(compound, threshold=0.5)
        assert isinstance(primary, List)
        assert all(m in toxicity_predictor.mechanism_types for m in primary)
        assert len(primary) <= len(mechanisms)

    # Analyze mechanism score distributions
    score_matrix = np.array(all_scores)
    assert score_matrix.shape == (len(test_compounds), len(toxicity_predictor.mechanism_types))
    assert np.all(np.mean(score_matrix[:3], axis=1) > np.mean(score_matrix[5:], axis=1))


def test_enhanced_prediction(
    enhanced_predictor: ToxicityPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test enhanced toxicity prediction."""
    # Train predictor
    enhanced_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1, 0.1],
    )

    # Test predictions with evidence
    for compound in test_compounds:
        prediction = enhanced_predictor.predict_with_evidence(compound)
        assert isinstance(prediction, Dict)
        assert "toxicity_score" in prediction
        assert "confidence" in prediction
        assert "toxicity_types" in prediction
        assert "mechanisms" in prediction
        assert "literature_evidence" in prediction
        assert "regulatory_status" in prediction
        assert "risk_factors" in prediction
        assert "safety_data" in prediction
        assert isinstance(prediction["toxicity_score"], float)
        assert isinstance(prediction["confidence"], float)
        assert isinstance(prediction["toxicity_types"], Dict)
        assert isinstance(prediction["mechanisms"], Dict)
        assert isinstance(prediction["literature_evidence"], List)
        assert isinstance(prediction["regulatory_status"], Dict)
        assert isinstance(prediction["risk_factors"], Dict)
        assert isinstance(prediction["safety_data"], Dict)

        # Check literature evidence structure
        for evidence in prediction["literature_evidence"]:
            assert "source" in evidence
            assert "title" in evidence
            assert "year" in evidence
            assert "relevance" in evidence
            assert "findings" in evidence

        # Check regulatory status structure
        assert "warnings" in prediction["regulatory_status"]
        assert "contraindications" in prediction["regulatory_status"]
        assert "monitoring" in prediction["regulatory_status"]
        assert "reporting" in prediction["regulatory_status"]

        # Check risk factors structure
        assert "dose_dependent" in prediction["risk_factors"]
        assert "cumulative" in prediction["risk_factors"]
        assert "interactions" in prediction["risk_factors"]
        assert "populations" in prediction["risk_factors"]

        # Check safety data structure
        assert "ld50" in prediction["safety_data"]
        assert "noael" in prediction["safety_data"]
        assert "acute_toxicity" in prediction["safety_data"]
        assert "chronic_toxicity" in prediction["safety_data"]


def test_model_persistence(
    toxicity_predictor: ToxicityPredictor,
    test_compounds: List[CompoundData],
    tmp_path: Path,
):
    """Test model persistence."""
    # Train predictor
    toxicity_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1, 0.1],
    )

    # Save models
    save_dir = tmp_path / "saved_models"
    toxicity_predictor.save_models(save_dir)

    # Create new predictor and load models
    new_predictor = ToxicityPredictor(
        model_dir=str(save_dir),
        cache_dir=str(tmp_path / "new_cache"),
    )

    # Compare predictions
    for compound in test_compounds:
        score1, conf1 = toxicity_predictor.predict(compound)
        score2, conf2 = new_predictor.predict(compound)
        assert abs(score1 - score2) < 1e-6
        assert abs(conf1 - conf2) < 1e-6

        tox1 = toxicity_predictor.predict_toxicity_types(compound)
        tox2 = new_predictor.predict_toxicity_types(compound)
        assert set(tox1.keys()) == set(tox2.keys())
        assert all(abs(tox1[k] - tox2[k]) < 1e-6 for k in tox1)

        mech1 = toxicity_predictor.predict_mechanisms(compound)
        mech2 = new_predictor.predict_mechanisms(compound)
        assert set(mech1.keys()) == set(mech2.keys())
        assert all(abs(mech1[k] - mech2[k]) < 1e-6 for k in mech1)


def test_error_handling(toxicity_predictor: ToxicityPredictor):
    """Test error handling."""
    # Test with invalid compound
    invalid_compound = CompoundData(
        name="Invalid",
        smiles="INVALID",
        cas_number="000-00-0",
    )

    with pytest.raises(ValueError):
        toxicity_predictor.predict(invalid_compound)

    with pytest.raises(ValueError):
        toxicity_predictor.predict_toxicity_types(invalid_compound)

    with pytest.raises(ValueError):
        toxicity_predictor.predict_mechanisms(invalid_compound)

    # Test with untrained model
    valid_compound = CompoundData(
        name="Valid",
        smiles="CCO",
        cas_number="64-17-5",
    )

    with pytest.raises(RuntimeError):
        toxicity_predictor.predict(valid_compound)

    with pytest.raises(RuntimeError):
        toxicity_predictor.predict_toxicity_types(valid_compound)

    with pytest.raises(RuntimeError):
        toxicity_predictor.predict_mechanisms(valid_compound)


def test_prediction_metrics(
    toxicity_predictor: ToxicityPredictor,
    test_compounds: List[CompoundData],
):
    """Test prediction metrics."""
    # Train predictor
    metrics = toxicity_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1, 0.1],
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
    binary_metrics = toxicity_predictor.get_binary_metrics(threshold=0.5)
    assert "accuracy" in binary_metrics
    assert "precision" in binary_metrics
    assert "recall" in binary_metrics
    assert "f1" in binary_metrics
    assert all(isinstance(v, float) for v in binary_metrics.values())
    assert all(0 <= v <= 1 for v in binary_metrics.values())


def test_enhanced_metrics(
    enhanced_predictor: ToxicityPredictorEnhanced,
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
        [0.1, 0.1, 0.1],
    )

    # Check metrics
    assert isinstance(metrics, Dict)
    assert "ml_metrics" in metrics
    assert "literature_metrics" in metrics
    assert "regulatory_metrics" in metrics
    assert "safety_metrics" in metrics
    assert "combined_metrics" in metrics
    assert all(isinstance(v, Dict) for v in metrics.values())
    assert all(all(isinstance(v, float) for v in m.values()) for m in metrics.values())
    assert all(all(v >= 0 for v in m.values()) for m in metrics.values())


def test_risk_analysis(
    enhanced_predictor: ToxicityPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test risk analysis for toxicity."""
    # Train predictor
    enhanced_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1, 0.1],
    )

    # Test risk analysis
    for compound in test_compounds:
        analysis = enhanced_predictor.analyze_risk(compound)
        assert isinstance(analysis, Dict)
        assert "toxicity_score" in analysis
        assert "acute_risk" in analysis
        assert "chronic_risk" in analysis
        assert "organ_systems" in analysis
        assert "contraindications" in analysis
        assert "drug_interactions" in analysis
        assert "exposure_limits" in analysis
        assert "safety_margins" in analysis
        assert isinstance(analysis["toxicity_score"], float)
        assert isinstance(analysis["acute_risk"], float)
        assert isinstance(analysis["chronic_risk"], Dict)
        assert isinstance(analysis["organ_systems"], Dict)
        assert isinstance(analysis["contraindications"], List)
        assert isinstance(analysis["drug_interactions"], Dict)
        assert isinstance(analysis["exposure_limits"], Dict)
        assert isinstance(analysis["safety_margins"], Dict)


def test_regulatory_analysis(
    enhanced_predictor: ToxicityPredictorEnhanced,
    test_compounds: List[CompoundData],
):
    """Test regulatory analysis for toxicity."""
    # Train predictor
    enhanced_predictor.train(
        test_compounds[:3],
        [1.0, 0.95, 0.9],
        test_compounds[3:5],
        [0.6, 0.5],
        test_compounds[5:],
        [0.1, 0.1, 0.1],
    )

    # Test regulatory analysis
    for compound in test_compounds:
        analysis = enhanced_predictor.analyze_regulatory_status(compound)
        assert isinstance(analysis, Dict)
        assert "safety_warnings" in analysis
        assert "contraindications" in analysis
        assert "monitoring_requirements" in analysis
        assert "reporting_requirements" in analysis
        assert "exposure_limits" in analysis
        assert "safety_guidelines" in analysis
        assert "risk_mitigation" in analysis
        assert "regulatory_history" in analysis
        assert isinstance(analysis["safety_warnings"], List)
        assert isinstance(analysis["contraindications"], Dict)
        assert isinstance(analysis["monitoring_requirements"], Dict)
        assert isinstance(analysis["reporting_requirements"], List)
        assert isinstance(analysis["exposure_limits"], Dict)
        assert isinstance(analysis["safety_guidelines"], Dict)
        assert isinstance(analysis["risk_mitigation"], Dict)
        assert isinstance(analysis["regulatory_history"], List)
