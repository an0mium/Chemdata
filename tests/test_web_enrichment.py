"""Tests for web data enrichment functionality."""

import json
import tempfile
from pathlib import Path
from unittest import mock

import pytest
import requests

from binding_data_processor.config import Config
from web_enrichment.data_sources.community import CommunityDataSource
from web_enrichment.data_sources.pubchem import PubChemDataSource
from web_enrichment.data_sources.social import SocialDataSource
from web_enrichment.data_sources.swiss import SwissDataSource


@pytest.fixture
def test_config():
    """Create test configuration."""
    with tempfile.TemporaryDirectory() as temp_dir:
        config = Config()
        config.cache_dir = Path(temp_dir) / "cache"
        config.cache_dir.mkdir(parents=True)
        yield config


@pytest.fixture
def mock_pubchem_response():
    """Create mock PubChem API response."""
    return {
        "PC_Compounds": [
            {
                "id": {"id": {"cid": 123456}},
                "props": [
                    {
                        "urn": {"label": "IUPAC Name"},
                        "value": {"sval": "Test Compound"},
                    },
                    {
                        "urn": {"label": "Molecular Weight"},
                        "value": {"fval": 395.47},
                    },
                    {
                        "urn": {"label": "XLogP3"},
                        "value": {"fval": 3.2},
                    },
                ],
            }
        ]
    }


@pytest.fixture
def mock_swiss_response():
    """Create mock Swiss* API response."""
    return {
        "target_prediction": [
            {
                "target": "5-HT2A",
                "probability": 0.85,
                "common_name": "Serotonin 2A receptor",
            },
            {
                "target": "DRD2",
                "probability": 0.75,
                "common_name": "Dopamine D2 receptor",
            },
        ],
        "adme_prediction": {
            "absorption": 0.8,
            "distribution": 0.7,
            "metabolism": 0.6,
            "excretion": 0.5,
        },
    }


@pytest.fixture
def mock_community_data():
    """Create mock community data."""
    return {
        "reports": [
            {
                "source": "erowid",
                "url": "https://erowid.org/experiences/123",
                "title": "Test Report 1",
                "text": "Test experience report content 1",
                "date": "2023-01-01",
            },
            {
                "source": "psychonautwiki",
                "url": "https://psychonautwiki.org/wiki/Test",
                "title": "Test Report 2",
                "text": "Test experience report content 2",
                "date": "2023-02-01",
            },
        ],
        "safety": {
            "warnings": ["May cause drowsiness"],
            "interactions": ["MAOIs", "SSRIs"],
            "contraindications": ["Heart conditions"],
        },
    }


@pytest.fixture
def mock_social_data():
    """Create mock social media data."""
    return {
        "reddit": [
            {
                "subreddit": "DrugNerds",
                "title": "Test Post 1",
                "text": "Test content 1",
                "url": "https://reddit.com/r/DrugNerds/123",
                "date": "2023-01-01",
            },
        ],
        "twitter": [
            {
                "text": "Test tweet about compound",
                "url": "https://twitter.com/user/123",
                "date": "2023-02-01",
            },
        ],
    }


def test_pubchem_enrichment(test_config, mock_pubchem_response):
    """Test PubChem data enrichment."""
    source = PubChemDataSource(config=test_config)

    # Mock PubChem API response
    with mock.patch("requests.get") as mock_get:
        mock_get.return_value.json.return_value = mock_pubchem_response
        mock_get.return_value.status_code = 200

        # Test successful enrichment
        data = source.get_compound_data(
            "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
        )

        assert data["name"] == "Test Compound"
        assert data["molecular_weight"] == 395.47
        assert data["logp"] == 3.2

        # Test rate limiting
        assert mock_get.call_count == 1  # Should use cache for subsequent calls
        data = source.get_compound_data(
            "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
        )
        assert mock_get.call_count == 1

        # Test error handling
        mock_get.return_value.status_code = 404
        with pytest.raises(requests.exceptions.RequestException):
            source.get_compound_data("invalid_smiles")


def test_swiss_enrichment(test_config, mock_swiss_response):
    """Test Swiss* data enrichment."""
    source = SwissDataSource(config=test_config)

    # Mock Swiss* API response
    with mock.patch("requests.post") as mock_post:
        mock_post.return_value.json.return_value = mock_swiss_response
        mock_post.return_value.status_code = 200

        # Test successful enrichment
        data = source.get_compound_data(
            "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
        )

        assert len(data["target_prediction"]) == 2
        assert data["target_prediction"][0]["target"] == "5-HT2A"
        assert data["target_prediction"][0]["probability"] == 0.85
        assert "adme_prediction" in data

        # Test rate limiting
        assert mock_post.call_count == 1  # Should use cache for subsequent calls
        data = source.get_compound_data(
            "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
        )
        assert mock_post.call_count == 1

        # Test error handling
        mock_post.return_value.status_code = 500
        with pytest.raises(requests.exceptions.RequestException):
            source.get_compound_data("invalid_smiles")


def test_community_enrichment(test_config, mock_community_data):
    """Test community data enrichment."""
    source = CommunityDataSource(config=test_config)

    # Mock web scraping
    with mock.patch(
        "web_enrichment.data_sources.community.CommunityDataSource._scrape_data"
    ) as mock_scrape:
        mock_scrape.return_value = mock_community_data

        # Test successful enrichment
        data = source.get_compound_data(
            "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
        )

        assert len(data["reports"]) == 2
        assert data["reports"][0]["source"] == "erowid"
        assert "safety" in data
        assert len(data["safety"]["warnings"]) == 1

        # Test caching
        assert mock_scrape.call_count == 1  # Should use cache for subsequent calls
        data = source.get_compound_data(
            "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
        )
        assert mock_scrape.call_count == 1

        # Test error handling
        mock_scrape.side_effect = Exception("Scraping error")
        with pytest.raises(Exception):
            source.get_compound_data("invalid_smiles")


def test_social_enrichment(test_config, mock_social_data):
    """Test social media data enrichment."""
    source = SocialDataSource(
        config=test_config,
        reddit_client_id="test_id",
        reddit_client_secret="test_secret",
        twitter_api_key="test_key",
        twitter_api_secret="test_secret",
    )

    # Mock API responses
    with mock.patch(
        "web_enrichment.data_sources.social.SocialDataSource._get_reddit_data"
    ) as mock_reddit:
        mock_reddit.return_value = mock_social_data["reddit"]

        with mock.patch(
            "web_enrichment.data_sources.social.SocialDataSource._get_twitter_data"
        ) as mock_twitter:
            mock_twitter.return_value = mock_social_data["twitter"]

            # Test successful enrichment
            data = source.get_compound_data(
                "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
            )

            assert len(data["reddit"]) == 1
            assert data["reddit"][0]["subreddit"] == "DrugNerds"
            assert len(data["twitter"]) == 1

            # Test caching
            assert mock_reddit.call_count == 1  # Should use cache
            assert mock_twitter.call_count == 1
            data = source.get_compound_data(
                "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
            )
            assert mock_reddit.call_count == 1
            assert mock_twitter.call_count == 1

            # Test error handling
            mock_reddit.side_effect = Exception("Reddit API error")
            mock_twitter.side_effect = Exception("Twitter API error")
            with pytest.raises(Exception):
                source.get_compound_data("invalid_smiles")


def test_data_validation(test_config):
    """Test data validation."""
    source = PubChemDataSource(config=test_config)

    # Test invalid SMILES
    with pytest.raises(ValueError):
        source.get_compound_data("invalid_smiles")

    # Test missing required fields
    with mock.patch("requests.get") as mock_get:
        mock_get.return_value.json.return_value = {"invalid": "response"}
        mock_get.return_value.status_code = 200

        with pytest.raises(ValueError):
            source.get_compound_data(
                "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
            )


def test_cache_management(test_config, mock_pubchem_response):
    """Test cache management."""
    source = PubChemDataSource(config=test_config)

    # Mock API response
    with mock.patch("requests.get") as mock_get:
        mock_get.return_value.json.return_value = mock_pubchem_response
        mock_get.return_value.status_code = 200

        # Initial request should hit API
        data1 = source.get_compound_data(
            "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
        )
        assert mock_get.call_count == 1

        # Subsequent request should use cache
        data2 = source.get_compound_data(
            "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
        )
        assert mock_get.call_count == 1
        assert data1 == data2

        # Check cache file
        cache_file = test_config.cache_dir / "pubchem_cache.json"
        assert cache_file.exists()
        with open(cache_file) as f:
            cache = json.load(f)
            assert len(cache) == 1


def test_rate_limiting(test_config, mock_pubchem_response):
    """Test rate limiting."""
    source = PubChemDataSource(config=test_config)

    # Mock API response
    with mock.patch("requests.get") as mock_get:
        mock_get.return_value.json.return_value = mock_pubchem_response
        mock_get.return_value.status_code = 200

        # Make multiple requests
        for _ in range(5):
            source.get_compound_data(
                "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
            )

        # Should only make one actual API call due to caching
        assert mock_get.call_count == 1

        # Test rate limit error
        mock_get.return_value.status_code = 429
        with pytest.raises(requests.exceptions.RequestException):
            source.get_compound_data("different_smiles")


def test_error_recovery(test_config, mock_pubchem_response):
    """Test error recovery."""
    source = PubChemDataSource(config=test_config)

    # Mock API response
    with mock.patch("requests.get") as mock_get:
        # First call fails
        mock_get.return_value.status_code = 500
        with pytest.raises(requests.exceptions.RequestException):
            source.get_compound_data(
                "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
            )

        # Second call succeeds
        mock_get.return_value.json.return_value = mock_pubchem_response
        mock_get.return_value.status_code = 200
        data = source.get_compound_data(
            "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC"
        )

        assert data["name"] == "Test Compound"
