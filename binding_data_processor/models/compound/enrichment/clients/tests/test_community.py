"""Tests for community data enrichment client.

This module tests the CommunityClient which integrates with both scientific
(PubChem, ChEMBL) and community (PsychonautWiki, Erowid, TripSit) data sources,
including ML models for text analysis and entity recognition.
"""

import pytest
from unittest.mock import Mock, patch

from ..community import CommunityClient
from ..base import WebClientError


# Base Fixtures
@pytest.fixture
def mock_http_client():
    """Mock HTTP client fixture."""
    return Mock()


@pytest.fixture
def mock_logger():
    """Mock logger fixture."""
    return Mock()


# ML Model Fixtures
@pytest.fixture
def mock_text_classifier():
    """Mock text classifier fixture."""
    mock = Mock()
    mock.return_value = [{
        "label": "positive",
        "score": 0.95
    }]
    return mock


@pytest.fixture
def mock_ner_model():
    """Mock NER model fixture."""
    mock = Mock()
    mock.return_value = [{
        "word": "caffeine",
        "entity": "COMPOUND",
        "score": 0.9
    }]
    return mock


# Scientific Data Source Fixtures
@pytest.fixture
def mock_pubchem_response():
    """Mock PubChem API response fixture."""
    return {
        "data": {
            "compound": {
                "cid": "2519",
                "iupac_name": "1,3,7-trimethylpurine-2,6-dione",
                "molecular_formula": "C8H10N4O2",
                "molecular_weight": 194.19,
                "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                "inchi": "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
                "properties": {
                    "xlogp": -0.07,
                    "hbd": 0,
                    "hba": 6,
                    "rotatable_bonds": 0,
                    "tpsa": 58.44,
                    "complexity": 234,
                },
                "bioactivity": {
                    "assay_count": 543,
                    "active_assay_count": 234,
                    "activity_types": [
                        "Ki",
                        "IC50",
                        "EC50",
                    ],
                    "target_count": 156,
                },
            },
            "references": [
                {
                    "pubmed_id": "12345678",
                    "doi": "10.1234/example.123",
                    "title": "Example study of caffeine",
                    "journal": "Journal of Example Studies",
                    "year": 2023,
                },
            ],
        },
        "metadata": {
            "query_time": 0.08,
            "api_version": "2023.1",
        },
    }


@pytest.fixture
def mock_chembl_response():
    """Mock ChEMBL API response fixture."""
    return {
        "data": {
            "compound": {
                "chembl_id": "CHEMBL113",
                "pref_name": "CAFFEINE",
                "molecule_type": "Small molecule",
                "max_phase": 4,
                "molecule_properties": {
                    "alogp": -0.07,
                    "hbd": 0,
                    "hba": 6,
                    "psa": 58.44,
                    "rtb": 0,
                    "ro3_pass": "Y",
                    "num_ro5_violations": 0,
                },
                "molecule_hierarchy": {
                    "molecule_chembl_id": "CHEMBL113",
                    "parent_chembl_id": "CHEMBL113",
                },
                "cross_references": [
                    {
                        "xref_id": "2519",
                        "xref_name": "PubChem",
                    },
                    {
                        "xref_id": "P30542",
                        "xref_name": "UniProt",
                    },
                ],
            },
            "activities": [
                {
                    "activity_id": 1234567,
                    "assay_id": "CHEMBL1234567",
                    "target_chembl_id": "CHEMBL2034",
                    "standard_type": "Ki",
                    "standard_value": 12000.0,
                    "standard_units": "nM",
                    "pchembl_value": 4.92,
                },
            ],
            "targets": [
                {
                    "target_chembl_id": "CHEMBL2034",
                    "pref_name": "Adenosine A1 receptor",
                    "target_type": "SINGLE PROTEIN",
                    "organism": "Homo sapiens",
                    "tax_id": 9606,
                },
            ],
        },
        "metadata": {
            "query_time": 0.12,
            "api_version": "2023.1",
        },
    }


# Community Data Source Fixtures
@pytest.fixture
def mock_psychonautwiki_response():
    """Mock PsychonautWiki API response fixture."""
    return {
        "name": "Caffeine",
        "url": "https://psychonautwiki.org/wiki/Caffeine",
        "effects": [
            "Stimulation",
            "Anxiety",
            "Increased heart rate",
        ],
        "dosage": {
            "threshold": "20 mg",
            "light": "50-100 mg",
            "common": "100-200 mg",
            "strong": "200-500 mg",
            "heavy": "500+ mg",
        },
        "duration": {
            "onset": "15-45 minutes",
            "duration": "2-5 hours",
            "after_effects": "2-4 hours",
        },
        "risks": [
            "Anxiety",
            "Insomnia",
            "Dependence",
        ],
        "interactions": [
            {
                "substance": "Alcohol",
                "note": "May reduce perception of alcohol intoxication",
            },
        ],
    }


@pytest.fixture
def mock_erowid_response():
    """Mock Erowid API response fixture."""
    return {
        "name": "Caffeine",
        "url": "https://erowid.org/chemicals/caffeine/",
        "experiences": [
            {
                "title": "First Time Experience",
                "date": "2023-01-01",
                "rating": 4,
                "body": "Sample experience report",
                "author": "Anonymous",
                "tags": ["stimulant", "positive"],
            }
        ],
        "health": [
            "May cause anxiety in sensitive individuals",
            "Can lead to physical dependence",
        ],
        "effects": {
            "positive": [
                "Increased alertness",
                "Enhanced focus",
            ],
            "neutral": [
                "Mild stimulation",
            ],
            "negative": [
                "Anxiety",
                "Insomnia",
            ],
        },
        "dosage": {
            "notes": "Dosage information from user reports",
            "ranges": {
                "light": "50-100mg",
                "common": "100-200mg",
                "strong": "200-500mg",
            },
        },
    }


@pytest.fixture
def mock_tripsit_response():
    """Mock TripSit API response fixture."""
    return {
        "name": "Caffeine",
        "url": "https://drugs.tripsit.me/caffeine",
        "properties": {
            "dose_light": "50-100mg",
            "dose_common": "100-200mg",
            "dose_strong": "200-500mg",
            "duration_total": "2-5 hours",
            "roas": ["oral", "sublingual"],
        },
        "interactions": [
            {
                "substance": "Alcohol",
                "status": "Caution",
                "note": "May reduce perception of alcohol intoxication",
            },
        ],
        "categories": [
            "stimulant",
            "common",
            "legal",
        ],
        "summary": "Common stimulant found in coffee and tea",
        "effects": {
            "physical": [
                "Increased heart rate",
                "Increased blood pressure",
            ],
            "cognitive": [
                "Alertness",
                "Focus",
            ],
        },
    }


@pytest.fixture
def community_client(
    mock_http_client,
    mock_logger,
    mock_text_classifier,
    mock_ner_model,
    tmp_path,
):
    """Community client fixture with mocked dependencies."""
    with patch("transformers.pipeline") as mock_pipeline:
        mock_pipeline.return_value = mock_text_classifier
        
        client = CommunityClient(
            http_client=mock_http_client,
            logger=mock_logger,
            model_dir=tmp_path,
            cache_dir=tmp_path,
        )
        
        # Override NER model
        client.ner_model = mock_ner_model
        
        return client


# Initialization Tests
def test_init(mock_http_client, mock_logger):
    """Test client initialization."""
    client = CommunityClient(
        http_client=mock_http_client,
        logger=mock_logger,
    )
    
    assert client.name == "community"
    assert client.http == mock_http_client
    assert client.logger == mock_logger
    assert client.pubchem_url == "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    assert client.chembl_url == "https://www.ebi.ac.uk/chembl/api/data"
    assert client.psychonautwiki_url == "https://api.psychonautwiki.org/v1"
    assert client.erowid_url == "https://api.erowid.org/v1"
    assert client.tripsit_url == "https://api.tripsit.me/v1"


def test_init_defaults():
    """Test initialization with default values."""
    client = CommunityClient()
    
    assert client.name == "community"
    assert client.http is not None
    assert client.logger is not None
    assert client.text_classifier is None
    assert client.ner_model is None


# Scientific Data Source Tests
def test_get_pubchem_data_success(
    community_client,
    mock_http_client,
    mock_pubchem_response,
):
    """Test successful PubChem data retrieval."""
    mock_http_client.get.return_value = mock_pubchem_response
    
    response = community_client.get_pubchem_data(
        compound_id="2519",
        name="Caffeine",
    )
    
    assert response == mock_pubchem_response
    mock_http_client.get.assert_called_with(
        url=f"{community_client.pubchem_url}/compound/2519",
        params={
            "name": "Caffeine",
            "operation": "record",
            "format": "json",
        },
        headers={"Accept": "application/json"},
    )


def test_get_pubchem_data_error(community_client, mock_http_client):
    """Test PubChem data error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "Compound not found",
        status_code=404,
    )
    
    with pytest.raises(WebClientError) as exc_info:
        community_client.get_pubchem_data(
            compound_id="invalid",
            name="invalid",
        )
    
    assert exc_info.value.status_code == 404
    assert "Compound not found" in str(exc_info.value)


def test_get_chembl_data_success(
    community_client,
    mock_http_client,
    mock_chembl_response,
):
    """Test successful ChEMBL data retrieval."""
    mock_http_client.get.return_value = mock_chembl_response
    
    response = community_client.get_chembl_data(
        compound_id="CHEMBL113",
        name="Caffeine",
    )
    
    assert response == mock_chembl_response
    mock_http_client.get.assert_called_with(
        url=f"{community_client.chembl_url}/molecule/CHEMBL113",
        params={
            "name": "Caffeine",
            "format": "json",
            "include_references": True,
        },
        headers={"Accept": "application/json"},
    )


def test_get_chembl_data_error(community_client, mock_http_client):
    """Test ChEMBL data error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "Compound not found",
        status_code=404,
    )
    
    with pytest.raises(WebClientError) as exc_info:
        community_client.get_chembl_data(
            compound_id="invalid",
            name="invalid",
        )
    
    assert exc_info.value.status_code == 404
    assert "Compound not found" in str(exc_info.value)


# Community Data Source Tests
def test_get_psychonautwiki_data_success(
    community_client,
    mock_http_client,
    mock_psychonautwiki_response,
):
    """Test successful PsychonautWiki data retrieval."""
    mock_http_client.get.return_value = mock_psychonautwiki_response
    
    response = community_client.get_psychonautwiki_data(name="Caffeine")
    
    assert response == mock_psychonautwiki_response
    mock_http_client.get.assert_called_with(
        url=f"{community_client.psychonautwiki_url}/substances/Caffeine",
        headers={"Accept": "application/json"},
    )


def test_get_psychonautwiki_data_error(community_client, mock_http_client):
    """Test PsychonautWiki data error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "Substance not found",
        status_code=404,
    )
    
    with pytest.raises(WebClientError) as exc_info:
        community_client.get_psychonautwiki_data(name="invalid")
    
    assert exc_info.value.status_code == 404
    assert "Substance not found" in str(exc_info.value)


def test_get_erowid_data_success(
    community_client,
    mock_http_client,
    mock_erowid_response,
):
    """Test successful Erowid data retrieval."""
    mock_http_client.get.return_value = mock_erowid_response
    
    response = community_client.get_erowid_data(name="Caffeine")
    
    assert response == mock_erowid_response
    mock_http_client.get.assert_called_with(
        url=f"{community_client.erowid_url}/chemicals/Caffeine",
        headers={"Accept": "application/json"},
    )


def test_get_erowid_data_error(community_client, mock_http_client):
    """Test Erowid data error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "Chemical not found",
        status_code=404,
    )
    
    with pytest.raises(WebClientError) as exc_info:
        community_client.get_erowid_data(name="invalid")
    
    assert exc_info.value.status_code == 404
    assert "Chemical not found" in str(exc_info.value)


def test_get_tripsit_data_success(
    community_client,
    mock_http_client,
    mock_tripsit_response,
):
    """Test successful TripSit data retrieval."""
    mock_http_client.get.return_value = mock_tripsit_response
    
    response = community_client.get_tripsit_data(name="Caffeine")
    
    assert response == mock_tripsit_response
    mock_http_client.get.assert_called_with(
        url=f"{community_client.tripsit_url}/drugs/Caffeine",
        headers={"Accept": "application/json"},
    )


def test_get_tripsit_data_error(community_client, mock_http_client):
    """Test TripSit data error handling."""
    mock_http_client.get.side_effect = WebClientError(
        "Drug not found",
        status_code=404,
    )
    
    with pytest.raises(WebClientError) as exc_info:
        community_client.get_tripsit_data(name="invalid")
    
    assert exc_info.value.status_code == 404
    assert "Drug not found" in str(exc_info.value)


# Combined Data Tests
def test_get_all_data_success(
    community_client,
    mock_http_client,
    mock_pubchem_response,
    mock_chembl_response,
    mock_psychonautwiki_response,
    mock_erowid_response,
    mock_tripsit_response,
):
    """Test successful retrieval of all community data."""
    mock_http_client.get.side_effect = [
        mock_pubchem_response,
        mock_chembl_response,
        mock_psychonautwiki_response,
        mock_erowid_response,
        mock_tripsit_response,
    ]
    
    response = community_client.get_all_data(
        name="Caffeine",
        compound_id="2519",
        chembl_id="CHEMBL113",
    )
    
    assert response == {
        "pubchem": mock_pubchem_response,
        "chembl": mock_chembl_response,
        "psychonautwiki": mock_psychonautwiki_response,
        "erowid": mock_erowid_response,
        "tripsit": mock_tripsit_response,
    }


def test_get_all_data_partial_success(
    community_client,
    mock_http_client,
    mock_pubchem_response,
    mock_chembl_response,
    mock_psychonautwiki_response,
):
    """Test partial success in retrieving community data."""
    mock_http_client.get.side_effect = [
        mock_pubchem_response,
        mock_chembl_response,
        mock_psychonautwiki_response,
        WebClientError("API error", status_code=500),
        WebClientError("API error", status_code=500),
    ]
    
    response = community_client.get_all_data(
        name="Caffeine",
        compound_id="2519",
        chembl_id="CHEMBL113",
    )
    
    assert response == {
        "pubchem": mock_pubchem_response,
        "chembl": mock_chembl_response,
        "psychonautwiki": mock_psychonautwiki_response,
        "erowid": None,
        "tripsit": None,
    }


# ML Model Tests
def test_analyze_text_success(community_client):
    """Test successful text analysis."""
    text = "Caffeine is a widely used stimulant"
    
    classification = community_client.analyze_text(text)
    
    assert classification["label"] == "positive"
    assert classification["score"] == 0.95


def test_analyze_text_no_classifier(community_client):
    """Test text analysis with no classifier."""
    community_client.text_classifier = None
    text = "Test text"
    
    classification = community_client.analyze_text(text)
    
    assert classification is None


def test_extract_entities_success(community_client):
    """Test successful entity extraction."""
    text = "Caffeine is a stimulant"
    
    entities = community_client.extract_entities(text)
    
    assert len(entities) == 1
    assert entities[0]["word"] == "caffeine"
    assert entities[0]["entity"] == "COMPOUND"
    assert entities[0]["score"] == 0.9


def test_extract_entities_no_model(community_client):
    """Test entity extraction with no model."""
    community_client.ner_model = None
    text = "Test text"
    
    entities = community_client.extract_entities(text)
    
    assert entities == []


# Validation Tests
def test_validate_name_success(community_client):
    """Test successful name validation."""
    valid_name = "Caffeine"
    
    # Should not raise any errors
    community_client._validate_name(valid_name)


def test_validate_name_empty(community_client):
    """Test name validation with empty string."""
    with pytest.raises(ValueError) as exc_info:
        community_client._validate_name("")
    
    assert "Name cannot be empty" in str(exc_info.value)


def test_validate_name_invalid_type(community_client):
    """Test name validation with invalid type."""
    with pytest.raises(TypeError) as exc_info:
        community_client._validate_name(123)
    
    assert "Name must be a string" in str(exc_info.value)


def test_validate_compound_id_success(community_client):
    """Test successful compound ID validation."""
    valid_id = "2519"
    
    # Should not raise any errors
    community_client._validate_compound_id(valid_id)


def test_validate_compound_id_empty(community_client):
    """Test compound ID validation with empty string."""
    with pytest.raises(ValueError) as exc_info:
        community_client._validate_compound_id("")
    
    assert "Compound ID cannot be empty" in str(exc_info.value)


def test_validate_compound_id_invalid_type(community_client):
    """Test compound ID validation with invalid type."""
    with pytest.raises(TypeError) as exc_info:
        community_client._validate_compound_id(123)
    
    assert "Compound ID must be a string" in str(exc_info.value)


def test_validate_chembl_id_success(community_client):
    """Test successful ChEMBL ID validation."""
    valid_id = "CHEMBL113"
    
    # Should not raise any errors
    community_client._validate_chembl_id(valid_id)


def test_validate_chembl_id_empty(community_client):
    """Test ChEMBL ID validation with empty string."""
    with pytest.raises(ValueError) as exc_info:
        community_client._validate_chembl_id("")
    
    assert "ChEMBL ID cannot be empty" in str(exc_info.value)


def test_validate_chembl_id_invalid_type(community_client):
    """Test ChEMBL ID validation with invalid type."""
    with pytest.raises(TypeError) as exc_info:
        community_client._validate_chembl_id(123)
    
    assert "ChEMBL ID must be a string" in str(exc_info.value)


def test_validate_chembl_id_invalid_format(community_client):
    """Test ChEMBL ID validation with invalid format."""
    with pytest.raises(ValueError) as exc_info:
        community_client._validate_chembl_id("invalid")
    
    assert "Invalid ChEMBL ID format" in str(exc_info.value)


# Response Validation Tests
def test_validate_pubchem_response(community_client, mock_pubchem_response):
    """Test PubChem response validation."""
    # Should not raise any errors
    community_client._validate_pubchem_response(mock_pubchem_response)
    
    # Missing required fields
    invalid_response = {"data": {"compound": {}}}
    with pytest.raises(WebClientError) as exc_info:
        community_client._validate_pubchem_response(invalid_response)
    assert "Missing required fields" in str(exc_info.value)
    
    # Invalid types
    invalid_response = {
        "data": {
            "compound": {
                "cid": 2519,  # Should be string
                "molecular_weight": "194.19",  # Should be number
            },
        },
    }
    with pytest.raises(WebClientError) as exc_info:
        community_client._validate_pubchem_response(invalid_response)
    assert "Invalid field types" in str(exc_info.value)


def test_validate_chembl_response(community_client, mock_chembl_response):
    """Test ChEMBL response validation."""
    # Should not raise any errors
    community_client._validate_chembl_response(mock_chembl_response)
    
    # Missing required fields
    invalid_response = {"data": {"compound": {}}}
    with pytest.raises(WebClientError) as exc_info:
        community_client._validate_chembl_response(invalid_response)
    assert "Missing required fields" in str(exc_info.value)
    
    # Invalid types
    invalid_response = {
        "data": {
            "compound": {
                "chembl_id": 113,  # Should be string
                "max_phase": "4",  # Should be number
            },
        },
    }
    with pytest.raises(WebClientError) as exc_info:
        community_client._validate_chembl_response(invalid_response)
    assert "Invalid field types" in str(exc_info.value)


def test_validate_psychonautwiki_response(
    community_client,
    mock_psychonautwiki_response,
):
    """Test PsychonautWiki response validation."""
    # Should not raise any errors
    community_client._validate_psychonautwiki_response(
        mock_psychonautwiki_response
    )
    
    # Missing required fields
    invalid_response = {"name": "Caffeine"}
    with pytest.raises(WebClientError) as exc_info:
        community_client._validate_psychonautwiki_response(invalid_response)
    assert "Missing required fields" in str(exc_info.value)
    
    # Invalid types
    invalid_response = {
        "name": "Caffeine",
        "effects": "not a list",  # Should be list
        "dosage": "not an object",  # Should be object
    }
    with pytest.raises(WebClientError) as exc_info:
        community_client._validate_psychonautwiki_response(invalid_response)
    assert "Invalid field types" in str(exc_info.value)


def test_validate_erowid_response(community_client, mock_erowid_response):
    """Test Erowid response validation."""
    # Should not raise any errors
    community_client._validate_erowid_response(mock_erowid_response)
    
    # Missing required fields
    invalid_response = {"name": "Caffeine"}
    with pytest.raises(WebClientError) as exc_info:
        community_client._validate_erowid_response(invalid_response)
    assert "Missing required fields" in str(exc_info.value)
    
    # Invalid types
    invalid_response = {
        "name": "Caffeine",
        "experiences": "not a list",  # Should be list
        "effects": "not an object",  # Should be object
    }
    with pytest.raises(WebClientError) as exc_info:
        community_client._validate_erowid_response(invalid_response)
    assert "Invalid field types" in str(exc_info.value)


def test_validate_tripsit_response(community_client, mock_tripsit_response):
    """Test TripSit response validation."""
    # Should not raise any errors
    community_client._validate_tripsit_response(mock_tripsit_response)
    
    # Missing required fields
    invalid_response = {"name": "Caffeine"}
    with pytest.raises(WebClientError) as exc_info:
        community_client._validate_tripsit_response(invalid_response)
    assert "Missing required fields" in str(exc_info.value)
    
    # Invalid types
    invalid_response = {
        "name": "Caffeine",
        "properties": "not an object",  # Should be object
        "categories": "not a list",  # Should be list
    }
    with pytest.raises(WebClientError) as exc_info:
        community_client._validate_tripsit_response(invalid_response)
    assert "Invalid field types" in str(exc_info.value)
