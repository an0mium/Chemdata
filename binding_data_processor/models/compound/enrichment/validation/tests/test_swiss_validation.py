"""Tests for Swiss data validation functions."""

import pytest

from ..data import validate_swiss_data, ValidationError


@pytest.fixture
def valid_target_prediction():
    """Valid target prediction fixture."""
    return {
        "target": "5-HT2A",
        "probability": 0.95,
        "common_name": "Serotonin 2A receptor",
        "chembl_id": "CHEMBL214",
        "target_class": "G protein-coupled receptor",
        "gene_name": "HTR2A",
        "organism": "Homo sapiens",
        "confidence_score": 0.85,
    }


@pytest.fixture
def valid_adme_properties():
    """Valid ADME properties fixture."""
    return {
        "mw": 194.19,
        "logp": -0.07,
        "hbd": 0,
        "hba": 6,
        "tpsa": 58.44,
        "rotatable_bonds": 0,
        "aromatic_rings": 2,
        "solubility": -2.15,
        "permeability": -4.52,
        "drug_score": 0.85,
    }


@pytest.fixture
def valid_swiss_data(valid_target_prediction, valid_adme_properties):
    """Valid Swiss data fixture."""
    return {
        "target_predictions": [valid_target_prediction],
        "adme_properties": valid_adme_properties,
        "metadata": {
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "inchi": "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
            "timestamp": "2023-01-01T00:00:00Z",
            "version": "3.0.0",
        },
    }


def test_validate_swiss_data_success(valid_swiss_data):
    """Test successful validation of Swiss data."""
    # Should not raise any errors
    validate_swiss_data(valid_swiss_data)


def test_validate_swiss_data_missing_required_fields():
    """Test validation with missing required fields."""
    # Missing target predictions
    data = {
        "adme_properties": {
            "mw": 194.19,
            "logp": -0.07,
            "hbd": 0,
            "hba": 6,
            "tpsa": 58.44,
        },
        "metadata": {
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "timestamp": "2023-01-01T00:00:00Z",
        },
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_swiss_data(data)
    assert "Missing required field 'target_predictions'" in str(exc_info.value)


def test_validate_target_prediction_invalid(valid_adme_properties):
    """Test validation with invalid target prediction."""
    data = {
        "target_predictions": [
            {
                "target": "5-HT2A",
                "probability": "not a number",  # Should be number
                "common_name": "Serotonin 2A receptor",
                "chembl_id": "CHEMBL214",
            }
        ],
        "adme_properties": valid_adme_properties,
        "metadata": {
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "timestamp": "2023-01-01T00:00:00Z",
        },
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_swiss_data(data)
    assert "Invalid probability value" in str(exc_info.value)


def test_validate_adme_properties_invalid(valid_target_prediction):
    """Test validation with invalid ADME properties."""
    data = {
        "target_predictions": [valid_target_prediction],
        "adme_properties": {
            "mw": "not a number",  # Should be number
            "logp": -0.07,
            "hbd": 0,
            "hba": 6,
            "tpsa": 58.44,
        },
        "metadata": {
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "timestamp": "2023-01-01T00:00:00Z",
        },
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_swiss_data(data)
    assert "Invalid molecular weight value" in str(exc_info.value)


def test_validate_metadata_invalid(valid_target_prediction, valid_adme_properties):
    """Test validation with invalid metadata."""
    data = {
        "target_predictions": [valid_target_prediction],
        "adme_properties": valid_adme_properties,
        "metadata": {
            "smiles": 123,  # Should be string
            "timestamp": "not a date",  # Should be ISO date
        },
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_swiss_data(data)
    assert "Invalid SMILES format" in str(exc_info.value)


def test_validate_probability_range(valid_adme_properties):
    """Test validation of probability range."""
    data = {
        "target_predictions": [
            {
                "target": "5-HT2A",
                "probability": 1.5,  # Should be between 0 and 1
                "common_name": "Serotonin 2A receptor",
                "chembl_id": "CHEMBL214",
            }
        ],
        "adme_properties": valid_adme_properties,
        "metadata": {
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "timestamp": "2023-01-01T00:00:00Z",
        },
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_swiss_data(data)
    assert "Probability must be between 0 and 1" in str(exc_info.value)


def test_validate_confidence_range(valid_adme_properties):
    """Test validation of confidence score range."""
    data = {
        "target_predictions": [
            {
                "target": "5-HT2A",
                "probability": 0.95,
                "common_name": "Serotonin 2A receptor",
                "chembl_id": "CHEMBL214",
                "confidence_score": 1.5,  # Should be between 0 and 1
            }
        ],
        "adme_properties": valid_adme_properties,
        "metadata": {
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "timestamp": "2023-01-01T00:00:00Z",
        },
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_swiss_data(data)
    assert "Confidence score must be between 0 and 1" in str(exc_info.value)


def test_validate_chembl_id_format(valid_adme_properties):
    """Test validation of ChEMBL ID format."""
    data = {
        "target_predictions": [
            {
                "target": "5-HT2A",
                "probability": 0.95,
                "common_name": "Serotonin 2A receptor",
                "chembl_id": "invalid_id",  # Should match CHEMBL pattern
            }
        ],
        "adme_properties": valid_adme_properties,
        "metadata": {
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "timestamp": "2023-01-01T00:00:00Z",
        },
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_swiss_data(data)
    assert "Invalid ChEMBL ID format" in str(exc_info.value)


def test_validate_numeric_ranges(valid_target_prediction):
    """Test validation of numeric property ranges."""
    data = {
        "target_predictions": [valid_target_prediction],
        "adme_properties": {
            "mw": -1.0,  # Should be positive
            "logp": -0.07,
            "hbd": -1,  # Should be non-negative integer
            "hba": 6,
            "tpsa": -58.44,  # Should be non-negative
        },
        "metadata": {
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "timestamp": "2023-01-01T00:00:00Z",
        },
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_swiss_data(data)
    assert "must be positive" in str(exc_info.value)


def test_validate_integer_properties(valid_target_prediction):
    """Test validation of integer properties."""
    data = {
        "target_predictions": [valid_target_prediction],
        "adme_properties": {
            "mw": 194.19,
            "logp": -0.07,
            "hbd": 0.5,  # Should be integer
            "hba": 6.5,  # Should be integer
            "tpsa": 58.44,
            "rotatable_bonds": 1.5,  # Should be integer
            "aromatic_rings": 2.5,  # Should be integer
        },
        "metadata": {
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "timestamp": "2023-01-01T00:00:00Z",
        },
    }
    
    with pytest.raises(ValidationError) as exc_info:
        validate_swiss_data(data)
    assert "must be an integer" in str(exc_info.value)
