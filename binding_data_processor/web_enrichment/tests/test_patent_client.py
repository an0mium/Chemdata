"""Tests for patent client."""

import pytest
from unittest.mock import Mock, patch
from datetime import datetime

from ..clients.patent import PatentClient
from ...models.compound import Compound


@pytest.fixture
def mock_http():
    """Mock HTTP client."""
    return Mock()


@pytest.fixture
def client(mock_http):
    """Patent client with mocked HTTP."""
    client = PatentClient(
        google_patents_key="test_key",
        lens_api_key="test_key",
    )
    client.http = mock_http
    return client


def test_search_patents(client, mock_http):
    """Test patent search."""
    # Mock Google Patents response
    mock_http.get.return_value.json.return_value = {
        "patents": [
            {
                "publication_number": "US123456",
                "title": "Test Patent",
                "abstract": "Test abstract",
                "filing_date": "2020-01-01",
                "publication_date": "2020-06-01",
                "assignee": "Test Company",
                "inventors": ["Test Inventor"],
                "description": "Test description",
                "claims": "Test claims",
            }
        ]
    }

    # Mock Lens.org response
    mock_http.post.return_value.json.return_value = {
        "data": [
            {
                "lens_id": "987654",
                "title": "Another Patent",
                "abstract": "Another abstract",
                "date_filed": "2020-02-01",
                "date_published": "2020-07-01",
                "assignee": "Another Company",
                "inventors": ["Another Inventor"],
                "description": "Another description",
                "claims": "Another claims",
            }
        ]
    }

    # Test search
    results = client.search_patents(
        query="test compound",
        chemical_structure="C1=CC=CC=C1",
    )

    assert len(results) == 2
    assert results[0]["patent_number"] == "US123456"
    assert results[1]["patent_number"] == "987654"

    # Verify API calls
    mock_http.get.assert_called_once()
    mock_http.post.assert_called_once()


def test_process_compounds(client):
    """Test compound processing."""
    # Create test compound
    compound = Compound(
        name="Test Compound",
        cas_number="123-45-6",
        smiles="C1=CC=CC=C1",
    )

    # Mock search results
    with patch.object(
        client,
        "search_patents",
        return_value=[
            {
                "patent_number": "US123456",
                "title": "Test Patent",
                "publication_date": "2020-01-01",
                "assignee": "Test Company",
            },
            {
                "patent_number": "US789012",
                "title": "Another Patent",
                "publication_date": "2020-02-01",
                "assignee": "Another Company",
            },
        ],
    ):
        client.process_compounds([compound])

    # Verify compound was updated
    assert hasattr(compound, "patents")
    assert len(compound.patents) == 2
    assert compound.patent_count == 2
    assert len(compound.patent_dates) == 2
    assert len(compound.patent_assignees) == 2
    assert "enrichment_metadata" in compound.__dict__
    assert "patent_search" in compound.enrichment_metadata


def test_get_patent_data(client, mock_http):
    """Test getting single patent data."""
    # Mock Google Patents response
    mock_http.get.return_value.json.return_value = {
        "publication_number": "US123456",
        "title": "Test Patent",
        "abstract": "Test abstract",
        "filing_date": "2020-01-01",
        "publication_date": "2020-06-01",
        "assignee": "Test Company",
        "inventors": ["Test Inventor"],
        "description": "Test description",
        "claims": "Test claims",
    }

    # Test getting patent
    patent = client.get_patent_data("US123456")

    assert patent is not None
    assert patent["patent_number"] == "US123456"
    assert "chemical_info" in patent
    assert "analysis" in patent

    # Verify API call
    mock_http.get.assert_called_once()


def test_error_handling(client, mock_http):
    """Test error handling."""
    # Mock API error
    mock_http.get.side_effect = Exception("API Error")

    # Test search error handling
    results = client.search_patents("test")
    assert len(results) == 0

    # Test get_patent_data error handling
    patent = client.get_patent_data("US123456")
    assert patent is None


def test_deduplication(client):
    """Test patent deduplication."""
    # Mock search results with duplicate
    with patch.object(
        client,
        "_search_google_patents",
        return_value=[
            {
                "patent_number": "US123456",
                "title": "Test Patent",
            },
            {
                "patent_number": "US123456",  # Duplicate
                "title": "Test Patent",
            },
        ],
    ), patch.object(
        client,
        "_search_lens",
        return_value=[
            {
                "patent_number": "US789012",
                "title": "Another Patent",
            },
            {
                "patent_number": "US123456",  # Duplicate
                "title": "Test Patent",
            },
        ],
    ):
        results = client.search_patents("test")

    # Verify duplicates were removed
    assert len(results) == 2
    patent_numbers = {r["patent_number"] for r in results}
    assert patent_numbers == {"US123456", "US789012"}


def test_metadata_timestamps(client):
    """Test metadata timestamps."""
    # Create test compound
    compound = Compound(
        name="Test Compound",
        cas_number="123-45-6",
    )

    # Mock search results
    with patch.object(
        client,
        "search_patents",
        return_value=[{"patent_number": "US123456"}],
    ), patch(
        "datetime.datetime",
        Mock(now=Mock(return_value=datetime(2023, 1, 1, 12, 0))),
    ):
        client.process_compounds([compound])

    # Verify timestamp
    assert compound.enrichment_metadata["patent_search"]["timestamp"] == "2023-01-01T12:00:00"


def test_empty_results(client):
    """Test handling of empty results."""
    # Mock empty responses
    mock_http.get.return_value.json.return_value = {"patents": []}
    mock_http.post.return_value.json.return_value = {"data": []}

    # Test search
    results = client.search_patents("test")
    assert len(results) == 0

    # Test compound processing
    compound = Compound(name="Test")
    client.process_compounds([compound])
    assert hasattr(compound, "patents")
    assert len(compound.patents) == 0
    assert compound.patent_count == 0
