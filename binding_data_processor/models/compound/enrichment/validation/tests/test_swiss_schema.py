"""Tests for Swiss data JSON schema."""

import pytest
from jsonschema import validate, ValidationError as JsonSchemaError

from ..schema import SWISS_DATA_SCHEMA


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


def test_swiss_data_schema_valid(valid_swiss_data):
    """Test schema validation with valid data."""
    # Should not raise any errors
    validate(instance=valid_swiss_data, schema=SWISS_DATA_SCHEMA)


def test_swiss_data_schema_missing_required():
    """Test schema validation with missing required fields."""
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
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=SWISS_DATA_SCHEMA)


def test_swiss_data_schema_invalid_types():
    """Test schema validation with invalid types."""
    data = {
        "target_predictions": [
            {
                "target": 123,  # Should be string
                "probability": "0.95",  # Should be number
                "common_name": ["not a string"],  # Should be string
                "chembl_id": 123,  # Should be string
            }
        ],
        "adme_properties": {
            "mw": "194.19",  # Should be number
            "logp": "-0.07",  # Should be number
            "hbd": "0",  # Should be integer
            "hba": "6",  # Should be integer
            "tpsa": "58.44",  # Should be number
        },
        "metadata": {
            "smiles": 123,  # Should be string
            "timestamp": 123,  # Should be string
        },
    }
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=SWISS_DATA_SCHEMA)


def test_swiss_data_schema_invalid_arrays():
    """Test schema validation with invalid arrays."""
    data = {
        "target_predictions": "not an array",  # Should be array
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
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=SWISS_DATA_SCHEMA)


def test_swiss_data_schema_invalid_objects():
    """Test schema validation with invalid objects."""
    data = {
        "target_predictions": [],
        "adme_properties": "not an object",  # Should be object
        "metadata": "not an object",  # Should be object
    }
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=SWISS_DATA_SCHEMA)


def test_swiss_data_schema_invalid_ranges():
    """Test schema validation with invalid numeric ranges."""
    data = {
        "target_predictions": [
            {
                "target": "5-HT2A",
                "probability": 1.5,  # Should be between 0 and 1
                "common_name": "Serotonin 2A receptor",
                "chembl_id": "CHEMBL214",
                "confidence_score": -0.5,  # Should be between 0 and 1
            }
        ],
        "adme_properties": {
            "mw": -194.19,  # Should be positive
            "logp": -0.07,
            "hbd": -1,  # Should be non-negative
            "hba": -6,  # Should be non-negative
            "tpsa": -58.44,  # Should be non-negative
            "drug_score": 1.5,  # Should be between 0 and 1
        },
        "metadata": {
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "timestamp": "2023-01-01T00:00:00Z",
        },
    }
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=SWISS_DATA_SCHEMA)


def test_swiss_data_schema_invalid_integers():
    """Test schema validation with invalid integer values."""
    data = {
        "target_predictions": [
            {
                "target": "5-HT2A",
                "probability": 0.95,
                "common_name": "Serotonin 2A receptor",
                "chembl_id": "CHEMBL214",
            }
        ],
        "adme_properties": {
            "mw": 194.19,
            "logp": -0.07,
            "hbd": 1.5,  # Should be integer
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
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=SWISS_DATA_SCHEMA)


def test_swiss_data_schema_invalid_chembl():
    """Test schema validation with invalid ChEMBL IDs."""
    data = {
        "target_predictions": [
            {
                "target": "5-HT2A",
                "probability": 0.95,
                "common_name": "Serotonin 2A receptor",
                "chembl_id": "invalid_id",  # Should match CHEMBL pattern
            }
        ],
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
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=SWISS_DATA_SCHEMA)


def test_swiss_data_schema_invalid_timestamp():
    """Test schema validation with invalid timestamp format."""
    data = {
        "target_predictions": [
            {
                "target": "5-HT2A",
                "probability": 0.95,
                "common_name": "Serotonin 2A receptor",
                "chembl_id": "CHEMBL214",
            }
        ],
        "adme_properties": {
            "mw": 194.19,
            "logp": -0.07,
            "hbd": 0,
            "hba": 6,
            "tpsa": 58.44,
        },
        "metadata": {
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "timestamp": "not a timestamp",  # Should be ISO format
        },
    }
    
    with pytest.raises(JsonSchemaError):
        validate(instance=data, schema=SWISS_DATA_SCHEMA)
