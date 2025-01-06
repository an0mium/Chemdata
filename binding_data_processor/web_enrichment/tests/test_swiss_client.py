"""Tests for Swiss data enrichment client.

This module tests the SwissClient which integrates with SwissSimilarity,
SwissProt, SwissTargetPrediction, and SwissADME APIs.
"""

import pytest
from unittest.mock import Mock

from ..swiss import SwissClient
from ..base import WebClientError


@pytest.fixture
def mock_http_client():
    """Mock HTTP client fixture."""
    return Mock()


@pytest.fixture
def mock_logger():
    """Mock logger fixture."""
    return Mock()


@pytest.fixture
def valid_smiles():
    """Valid SMILES string fixture."""
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine


@pytest.fixture
def valid_accession():
    """Valid UniProt accession fixture."""
    return "P30542"  # Adenosine receptor A1


@pytest.fixture
def mock_swisssimilarity_response():
    """Mock SwissSimilarity API response fixture."""
    return {
        "data": {
            "query": {
                "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                "name": "Caffeine",
            },
            "similar_compounds": [
                {
                    "smiles": "CN1C=NC2=C1C(=O)NC(=O)N2C",
                    "name": "Theophylline",
                    "similarity": 0.95,
                    "tanimoto": 0.92,
                    "properties": {
                        "molecular_weight": 180.16,
                        "log_p": -0.02,
                    },
                },
                {
                    "smiles": "CN1C=NC2=C1C(=O)N(C)C(=O)N2",
                    "name": "Theobromine",
                    "similarity": 0.92,
                    "tanimoto": 0.89,
                    "properties": {
                        "molecular_weight": 180.16,
                        "log_p": -0.78,
                    },
                },
            ],
            "statistics": {
                "total_compounds": 2,
                "average_similarity": 0.935,
                "average_tanimoto": 0.905,
            },
        },
        "metadata": {
            "query_time": 0.15,
            "api_version": "2023.1",
        },
    }


@pytest.fixture
def mock_swissprot_response():
    """Mock SwissProt API response fixture."""
    return {
        "data": {
            "entry": {
                "accession": "P30542",
                "name": "AA1R_HUMAN",
                "protein": {
                    "recommendedName": {
                        "fullName": "Adenosine receptor A1",
                    },
                },
                "gene": [
                    {
                        "name": {
                            "primary": "ADORA1",
                        },
                    },
                ],
                "organism": {
                    "scientificName": "Homo sapiens",
                    "taxonId": 9606,
                },
                "proteinExistence": "Evidence at protein level",
                "keywords": [
                    {
                        "id": "KW-0025",
                        "category": "Molecular function",
                        "name": "G-protein coupled receptor",
                    },
                ],
                "features": [
                    {
                        "type": "BINDING",
                        "description": "Caffeine",
                        "begin": 123,
                        "end": 456,
                    },
                ],
                "references": [
                    {
                        "citation": {
                            "type": "journal article",
                            "title": "Example study of caffeine binding",
                            "authors": ["Smith J.", "Jones B."],
                            "publication": {
                                "journalName": "J Mol Biol",
                                "year": 2023,
                            },
                        },
                    },
                ],
            },
            "crossReferences": [
                {
                    "database": "PDB",
                    "id": "1ABC",
                },
                {
                    "database": "ChEMBL",
                    "id": "CHEMBL2034",
                },
            ],
        },
        "metadata": {
            "query_time": 0.08,
            "api_version": "2023.1",
        },
    }


@pytest.fixture
def mock_swissadme_response():
    """Mock SwissADME API response fixture."""
    return {
        "data": {
            "compound": {
                "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                "properties": {
                    "molecular": {
                        "formula": "C8H10N4O2",
                        "weight": 194.19,
                        "heavy_atoms": 14,
                        "arom_atoms": 9,
                        "rotatable_bonds": 0,
                        "hbd": 0,
                        "hba": 6,
                        "tpsa": 58.44,
                    },
                    "lipophilicity": {
                        "log_p": -0.07,
                        "consensus_log_p": -0.12,
                    },
                    "druglikeness": {
                        "lipinski": "Yes",
                        "ghose": "Yes",
                        "veber": "Yes",
                        "egan": "Yes",
                        "muegge": "Yes",
                    },
                    "medicinal_chemistry": {
                        "pains": 0,
                        "brenk": 0,
                        "synthetic_accessibility": 1.25,
                    },
                    "pharmacokinetics": {
                        "gi_absorption": "High",
                        "bbb_permeant": "Yes",
                        "pgp_substrate": "No",
                        "cyp_inhibition": {
                            "cyp1a2": "Yes",
                            "cyp2c19": "No",
                            "cyp2c9": "No",
                            "cyp2d6": "No",
                            "cyp3a4": "No",
                        },
                    },
                },
                "alerts": [
                    {
                        "type": "toxicity",
                        "name": "xanthine",
                        "description": "Potential CNS stimulant",
                    },
                ],
            },
            "predictions": [
                {
                    "property": "oral_bioavailability",
                    "value": 0.95,
                    "confidence": 0.85,
                },
            ],
        },
        "metadata": {
            "query_time": 0.12,
            "api_version": "2023.1",
        },
    }


@pytest.fixture
def mock_swisstargetprediction_response():
    """Mock SwissTargetPrediction API response fixture."""
    return {
        "data": {
            "query": {
                "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            },
            "predictions": [
                {
                    "target": {
                        "uniprot_id": "P30542",
                        "name": "Adenosine receptor A1",
                        "gene": "ADORA1",
                        "family": "7TM1",
                        "organism": "Human",
                    },
                    "probability": 0.95,
                    "known_active": True,
                    "details": {
                        "2d_similarity": 0.92,
                        "3d_similarity": 0.88,
                        "known_actives": 145,
                    },
                },
                {
                    "target": {
                        "uniprot_id": "P29274",
                        "name": "Adenosine receptor A2a",
                        "gene": "ADORA2A",
                        "family": "7TM1",
                        "organism": "Human",
                    },
                    "probability": 0.92,
                    "known_active": True,
                    "details": {
                        "2d_similarity": 0.90,
                        "3d_similarity": 0.85,
                        "known_actives": 156,
                    },
                },
            ],
        },
        "metadata": {
            "query_time": 1.25,
            "model_version": "2023.1",
        },
    }


@pytest.fixture
def swiss_client(mock_http_client, mock_logger):
    """Swiss client fixture with mocked dependencies."""
    return SwissClient(
        http_client=mock_http_client,
        logger=mock_logger,
    )


def test_init(mock_http_client, mock_logger):
    """Test client initialization."""
    client = SwissClient(
        http_client=mock_http_client,
        logger=mock_logger,
    )

    assert client.name == "swiss"
    assert client.http == mock_http_client
    assert client.logger == mock_logger
    assert client.swissprot_url == "https://rest.uniprot.org/uniprotkb"
    assert client.swissadme_url == "https://www.swissadme.ch/api"
    assert client.swisssimilarity_url == "https://www.swisssimilarity.ch/api"
    assert client.swisstargetprediction_url == "https://www.swisstargetprediction.ch/api"


def test_init_defaults():
    """Test initialization with default values."""
    client = SwissClient()

    assert client.name == "swiss"
    assert client.http is not None
    assert client.logger is not None


def test_get_swissprot_data_success(
    swiss_client,
    mock_http_client,
    mock_swissprot_response,
    valid_accession,
):
    """Test successful SwissProt data retrieval."""
    mock_http_client.get.return_value = mock_swissprot_response

    response = swiss_client.get_swissprot_data(
        accession=valid_accession,
        name="AA1R_HUMAN",
    )

    assert response == mock_swissprot_response
    mock_http_client.get.assert_called_with(
        url=f"{swiss_client.swissprot_url}/{valid_accession}",
        params={
            "name": "AA1R_HUMAN",
            "format": "json",
            "include_features": True,
        },
        headers={"Accept": "application/json"},
    )


def test_get_swissprot_data_error(swiss_client, mock_http_client):
    """Test SwissProt data error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "Entry not found",
        status_code=404,
    )

    with pytest.raises(WebClientError) as exc_info:
        swiss_client.get_swissprot_data(
            accession="invalid",
            name="invalid",
        )

    assert exc_info.value.status_code == 404
    assert "Entry not found" in str(exc_info.value)


def test_get_swissadme_data_success(
    swiss_client,
    mock_http_client,
    mock_swissadme_response,
    valid_smiles,
):
    """Test successful SwissADME data retrieval."""
    mock_http_client.get.return_value = mock_swissadme_response

    response = swiss_client.get_swissadme_data(smiles=valid_smiles)

    assert response == mock_swissadme_response
    mock_http_client.get.assert_called_with(
        url=f"{swiss_client.swissadme_url}/predict",
        params={
            "smiles": valid_smiles,
            "format": "json",
        },
        headers={"Accept": "application/json"},
    )


def test_get_swissadme_data_error(swiss_client, mock_http_client):
    """Test SwissADME data error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "Invalid SMILES",
        status_code=400,
    )

    with pytest.raises(WebClientError) as exc_info:
        swiss_client.get_swissadme_data(smiles="invalid")

    assert exc_info.value.status_code == 400
    assert "Invalid SMILES" in str(exc_info.value)


def test_get_swisssimilarity_data_success(
    swiss_client,
    mock_http_client,
    mock_swisssimilarity_response,
    valid_smiles,
):
    """Test successful SwissSimilarity data retrieval."""
    mock_http_client.get.return_value = mock_swisssimilarity_response

    response = swiss_client.get_swisssimilarity_data(
        smiles=valid_smiles,
        name="Caffeine",
        threshold=0.8,
        limit=10,
    )

    assert response == mock_swisssimilarity_response
    mock_http_client.get.assert_called_with(
        url=f"{swiss_client.swisssimilarity_url}/similar",
        params={
            "smiles": valid_smiles,
            "name": "Caffeine",
            "threshold": 0.8,
            "limit": 10,
            "format": "json",
        },
        headers={"Accept": "application/json"},
    )


def test_get_swisssimilarity_data_error(swiss_client, mock_http_client):
    """Test SwissSimilarity data error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "Invalid SMILES",
        status_code=400,
    )

    with pytest.raises(WebClientError) as exc_info:
        swiss_client.get_swisssimilarity_data(
            smiles="invalid",
            name="invalid",
        )

    assert exc_info.value.status_code == 400
    assert "Invalid SMILES" in str(exc_info.value)


def test_get_swisstargetprediction_data_success(
    swiss_client,
    mock_http_client,
    mock_swisstargetprediction_response,
    valid_smiles,
):
    """Test successful SwissTargetPrediction data retrieval."""
    mock_http_client.get.return_value = mock_swisstargetprediction_response

    response = swiss_client.get_swisstargetprediction_data(smiles=valid_smiles)

    assert response == mock_swisstargetprediction_response
    mock_http_client.get.assert_called_with(
        url=f"{swiss_client.swisstargetprediction_url}/predict",
        params={
            "smiles": valid_smiles,
            "format": "json",
        },
        headers={"Accept": "application/json"},
    )


def test_get_swisstargetprediction_data_error(swiss_client, mock_http_client):
    """Test SwissTargetPrediction data error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "Invalid SMILES",
        status_code=400,
    )

    with pytest.raises(WebClientError) as exc_info:
        swiss_client.get_swisstargetprediction_data(smiles="invalid")

    assert exc_info.value.status_code == 400
    assert "Invalid SMILES" in str(exc_info.value)


def test_get_all_data_success(
    swiss_client,
    mock_http_client,
    mock_swissprot_response,
    mock_swissadme_response,
    mock_swisssimilarity_response,
    mock_swisstargetprediction_response,
    valid_smiles,
    valid_accession,
):
    """Test successful retrieval of all Swiss data."""
    mock_http_client.get.side_effect = [
        mock_swissprot_response,
        mock_swissadme_response,
        mock_swisssimilarity_response,
        mock_swisstargetprediction_response,
    ]

    response = swiss_client.get_all_data(
        smiles=valid_smiles,
        name="Caffeine",
        accession=valid_accession,
        threshold=0.8,
        limit=10,
    )

    assert response == {
        "swissprot": mock_swissprot_response,
        "swissadme": mock_swissadme_response,
        "swisssimilarity": mock_swisssimilarity_response,
        "swisstargetprediction": mock_swisstargetprediction_response,
    }


def test_get_all_data_partial_success(
    swiss_client,
    mock_http_client,
    mock_swissprot_response,
    mock_swissadme_response,
    valid_smiles,
    valid_accession,
):
    """Test partial success in retrieving Swiss data."""
    mock_http_client.get.side_effect = [
        mock_swissprot_response,
        mock_swissadme_response,
        WebClientError("API error", status_code=500),
        WebClientError("API error", status_code=500),
    ]

    response = swiss_client.get_all_data(
        smiles=valid_smiles,
        name="Caffeine",
        accession=valid_accession,
    )

    assert response == {
        "swissprot": mock_swissprot_response,
        "swissadme": mock_swissadme_response,
        "swisssimilarity": None,
        "swisstargetprediction": None,
    }


# Validation Tests
def test_validate_accession_success(swiss_client, valid_accession):
    """Test successful accession validation."""
    # Should not raise any errors
    swiss_client._validate_accession(valid_accession)


def test_validate_accession_empty(swiss_client):
    """Test accession validation with empty string."""
    with pytest.raises(ValueError) as exc_info:
        swiss_client._validate_accession("")

    assert "Accession cannot be empty" in str(exc_info.value)


def test_validate_accession_invalid_type(swiss_client):
    """Test accession validation with invalid type."""
    with pytest.raises(TypeError) as exc_info:
        swiss_client._validate_accession(123)

    assert "Accession must be a string" in str(exc_info.value)


def test_validate_accession_invalid_format(swiss_client):
    """Test accession validation with invalid format."""
    with pytest.raises(ValueError) as exc_info:
        swiss_client._validate_accession("invalid")

    assert "Invalid UniProt accession format" in str(exc_info.value)


def test_validate_smiles_success(swiss_client, valid_smiles):
    """Test successful SMILES validation."""
    # Should not raise any errors
    swiss_client._validate_smiles(valid_smiles)


def test_validate_smiles_empty(swiss_client):
    """Test SMILES validation with empty string."""
    with pytest.raises(ValueError) as exc_info:
        swiss_client._validate_smiles("")

    assert "SMILES cannot be empty" in str(exc_info.value)


def test_validate_smiles_invalid_type(swiss_client):
    """Test SMILES validation with invalid type."""
    with pytest.raises(TypeError) as exc_info:
        swiss_client._validate_smiles(123)

    assert "SMILES must be a string" in str(exc_info.value)


def test_validate_name_success(swiss_client):
    """Test successful name validation."""
    valid_name = "Caffeine"

    # Should not raise any errors
    swiss_client._validate_name(valid_name)


def test_validate_name_empty(swiss_client):
    """Test name validation with empty string."""
    with pytest.raises(ValueError) as exc_info:
        swiss_client._validate_name("")

    assert "Name cannot be empty" in str(exc_info.value)


def test_validate_name_invalid_type(swiss_client):
    """Test name validation with invalid type."""
    with pytest.raises(TypeError) as exc_info:
        swiss_client._validate_name(123)

    assert "Name must be a string" in str(exc_info.value)


# Response Validation Tests
def test_validate_swissprot_response(swiss_client, mock_swissprot_response):
    """Test SwissProt response validation."""
    # Should not raise any errors
    swiss_client._validate_swissprot_response(mock_swissprot_response)

    # Missing required fields
    invalid_response = {"data": {"entry": {}}}
    with pytest.raises(WebClientError) as exc_info:
        swiss_client._validate_swissprot_response(invalid_response)
    assert "Missing required fields" in str(exc_info.value)

    # Invalid types
    invalid_response = {
        "data": {
            "entry": {
                "accession": 30542,  # Should be string
                "proteinExistence": 123,  # Should be string
            },
        },
    }
    with pytest.raises(WebClientError) as exc_info:
        swiss_client._validate_swissprot_response(invalid_response)
    assert "Invalid field types" in str(exc_info.value)


def test_validate_swissadme_response(swiss_client, mock_swissadme_response):
    """Test SwissADME response validation."""
    # Should not raise any errors
    swiss_client._validate_swissadme_response(mock_swissadme_response)

    # Missing required fields
    invalid_response = {"data": {"compound": {}}}
    with pytest.raises(WebClientError) as exc_info:
        swiss_client._validate_swissadme_response(invalid_response)
    assert "Missing required fields" in str(exc_info.value)

    # Invalid types
    invalid_response = {
        "data": {
            "compound": {
                "smiles": 123,  # Should be string
                "properties": "not an object",  # Should be object
            },
        },
    }
    with pytest.raises(WebClientError) as exc_info:
        swiss_client._validate_swissadme_response(invalid_response)
    assert "Invalid field types" in str(exc_info.value)


def test_validate_swisssimilarity_response(
    swiss_client,
    mock_swisssimilarity_response,
):
    """Test SwissSimilarity response validation."""
    # Should not raise any errors
    swiss_client._validate_swisssimilarity_response(mock_swisssimilarity_response)

    # Missing required fields
    invalid_response = {"data": {"query": {}}}
    with pytest.raises(WebClientError) as exc_info:
        swiss_client._validate_swisssimilarity_response(invalid_response)
    assert "Missing required fields" in str(exc_info.value)

    # Invalid types
    invalid_response = {
        "data": {
            "query": {
                "smiles": 123,  # Should be string
                "similar_compounds": "not a list",  # Should be list
            },
        },
    }
    with pytest.raises(WebClientError) as exc_info:
        swiss_client._validate_swisssimilarity_response(invalid_response)
    assert "Invalid field types" in str(exc_info.value)


def test_validate_swisstargetprediction_response(
    swiss_client,
    mock_swisstargetprediction_response,
):
    """Test SwissTargetPrediction response validation."""
    # Should not raise any errors
    swiss_client._validate_swisstargetprediction_response(mock_swisstargetprediction_response)

    # Missing required fields
    invalid_response = {"data": {"query": {}}}
    with pytest.raises(WebClientError) as exc_info:
        swiss_client._validate_swisstargetprediction_response(invalid_response)
    assert "Missing required fields" in str(exc_info.value)

    # Invalid types
    invalid_response = {
        "data": {
            "query": {
                "smiles": 123,  # Should be string
                "predictions": "not a list",  # Should be list
            },
        },
    }
    with pytest.raises(WebClientError) as exc_info:
        swiss_client._validate_swisstargetprediction_response(invalid_response)
    assert "Invalid field types" in str(exc_info.value)
