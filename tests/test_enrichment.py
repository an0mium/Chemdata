"""Tests for compound enrichment script."""

import os
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest
import pandas as pd

from binding_data_processor.data_sources.bindingdb import BindingDBSource
from binding_data_processor.web_enrichment.manager import (
    WebEnrichmentManager,
    EnrichmentConfig,
)


@pytest.fixture
def test_data():
    """Create test data file."""
    data = pd.DataFrame({
        "Ligand Name": ["Caffeine", "Amphetamine"],
        "Ligand SMILES": [
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "CC(N)CC1=CC=CC=C1",
        ],
        "Ligand InChI": ["", ""],
        "Ligand InChI Key": ["", ""],
        "Target Name": ["A2A", "DAT"],
        "Target Source Organism": ["Human", "Human"],
        "Ki (nM)": [0.8, 0.05],
        "IC50 (nM)": [None, None],
        "Kd (nM)": [None, None],
        "EC50 (nM)": [None, None],
        "DOI": ["10.1021/test1", "10.1021/test2"],
        "PubMed ID": [12345678, 23456789],
    })

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        data.to_csv(f, sep="\t", index=False)
        return Path(f.name)


@pytest.fixture
def mock_clients():
    """Create mock web clients."""
    with patch("binding_data_processor.web_enrichment.manager.HTTPClient") as mock_http, \
         patch("binding_data_processor.web_enrichment.manager.SwissClient") as mock_swiss, \
         patch("binding_data_processor.web_enrichment.manager.CommunityClient") as mock_community, \
         patch("binding_data_processor.web_enrichment.manager.SocialClient") as mock_social:
        
        yield {
            "http": mock_http,
            "swiss": mock_swiss,
            "community": mock_community,
            "social": mock_social,
        }


def test_enrichment_basic(test_data, mock_clients, tmp_path):
    """Test basic enrichment functionality."""
    # Set up test environment
    os.environ["REDDIT_CLIENT_ID"] = "test_id"
    os.environ["REDDIT_CLIENT_SECRET"] = "test_secret"
    os.environ["TWITTER_BEARER_TOKEN"] = "test_token"

    # Create output file
    output_file = tmp_path / "enriched.tsv"

    # Create config
    config = EnrichmentConfig(
        model_dir=tmp_path / "models",
        cache_dir=tmp_path / "cache",
        n_workers=1,
        batch_size=10,
    )

    # Load compounds
    bindingdb = BindingDBSource()
    compounds = bindingdb.load_compounds(test_data)

    # Create manager
    manager = WebEnrichmentManager(config)

    # Mock client methods
    manager.swiss_client.process_compounds = Mock()
    manager.community_client.process_compounds = Mock()
    manager.social_client.process_compounds = Mock()

    # Enrich compounds
    manager.enrich_compounds(compounds)

    # Check client calls
    assert manager.swiss_client.process_compounds.called
    assert manager.community_client.process_compounds.called
    assert manager.social_client.process_compounds.called

    # Check compound data
    for compound in compounds:
        assert isinstance(compound.swiss_data, dict)
        assert isinstance(compound.community_data, dict)
        assert isinstance(compound.social_data, dict)
        assert isinstance(compound.enrichment_metadata, dict)
        assert "timestamp" in compound.enrichment_metadata
        assert "sources" in compound.enrichment_metadata

    # Save enriched compounds
    df = pd.DataFrame([{
        "name": c.name,
        "smiles": c.smiles,
        "swiss_data": c.swiss_data,
        "community_data": c.community_data,
        "social_data": c.social_data,
    } for c in compounds])
    df.to_csv(output_file, sep="\t", index=False)

    # Verify output file
    assert output_file.exists()
    df_loaded = pd.read_csv(output_file, sep="\t")
    assert len(df_loaded) == len(compounds)
    assert all(c.name in df_loaded["name"].values for c in compounds)


def test_enrichment_skip_predictions(test_data, mock_clients, tmp_path):
    """Test enrichment with predictions skipped."""
    # Set up test environment
    os.environ["REDDIT_CLIENT_ID"] = "test_id"
    os.environ["REDDIT_CLIENT_SECRET"] = "test_secret"
    os.environ["TWITTER_BEARER_TOKEN"] = "test_token"

    # Create config
    config = EnrichmentConfig(
        model_dir=tmp_path / "models",
        cache_dir=tmp_path / "cache",
        skip_predictions=True,
    )

    # Load compounds
    bindingdb = BindingDBSource()
    compounds = bindingdb.load_compounds(test_data)

    # Create manager
    manager = WebEnrichmentManager(config)

    # Mock client methods
    manager.swiss_client.process_compounds = Mock()
    manager.community_client.process_compounds = Mock()
    manager.social_client.process_compounds = Mock()

    # Enrich compounds
    manager.enrich_compounds(compounds)

    # Check client calls
    assert not manager.swiss_client.process_compounds.called
    assert manager.community_client.process_compounds.called
    assert manager.social_client.process_compounds.called


def test_enrichment_skip_web_data(test_data, mock_clients, tmp_path):
    """Test enrichment with web data skipped."""
    # Set up test environment
    os.environ["REDDIT_CLIENT_ID"] = "test_id"
    os.environ["REDDIT_CLIENT_SECRET"] = "test_secret"
    os.environ["TWITTER_BEARER_TOKEN"] = "test_token"

    # Create config
    config = EnrichmentConfig(
        model_dir=tmp_path / "models",
        cache_dir=tmp_path / "cache",
        skip_web_data=True,
    )

    # Load compounds
    bindingdb = BindingDBSource()
    compounds = bindingdb.load_compounds(test_data)

    # Create manager
    manager = WebEnrichmentManager(config)

    # Mock client methods
    manager.swiss_client.process_compounds = Mock()
    manager.community_client.process_compounds = Mock()
    manager.social_client.process_compounds = Mock()

    # Enrich compounds
    manager.enrich_compounds(compounds)

    # Check client calls
    assert manager.swiss_client.process_compounds.called
    assert not manager.community_client.process_compounds.called
    assert not manager.social_client.process_compounds.called


def test_enrichment_error_handling(test_data, mock_clients, tmp_path):
    """Test enrichment error handling."""
    # Set up test environment
    os.environ["REDDIT_CLIENT_ID"] = "test_id"
    os.environ["REDDIT_CLIENT_SECRET"] = "test_secret"
    os.environ["TWITTER_BEARER_TOKEN"] = "test_token"

    # Create config
    config = EnrichmentConfig(
        model_dir=tmp_path / "models",
        cache_dir=tmp_path / "cache",
    )

    # Load compounds
    bindingdb = BindingDBSource()
    compounds = bindingdb.load_compounds(test_data)

    # Create manager
    manager = WebEnrichmentManager(config)

    # Mock client methods to raise errors
    manager.swiss_client.process_compounds = Mock(side_effect=Exception("Swiss error"))
    manager.community_client.process_compounds = Mock(side_effect=Exception("Community error"))
    manager.social_client.process_compounds = Mock(side_effect=Exception("Social error"))

    # Enrich compounds (should not raise errors)
    manager.enrich_compounds(compounds)

    # Check compound data still initialized
    for compound in compounds:
        assert isinstance(compound.swiss_data, dict)
        assert isinstance(compound.community_data, dict)
        assert isinstance(compound.social_data, dict)
        assert isinstance(compound.enrichment_metadata, dict)


def test_enrichment_resume(test_data, mock_clients, tmp_path):
    """Test enrichment resume functionality."""
    # Set up test environment
    os.environ["REDDIT_CLIENT_ID"] = "test_id"
    os.environ["REDDIT_CLIENT_SECRET"] = "test_secret"
    os.environ["TWITTER_BEARER_TOKEN"] = "test_token"

    # Create checkpoint file
    checkpoint_dir = tmp_path / "checkpoints"
    checkpoint_dir.mkdir()
    checkpoint_file = checkpoint_dir / "enrichment.json"
    checkpoint_file.write_text('{"completed": ["Caffeine"]}')

    # Create config
    config = EnrichmentConfig(
        model_dir=tmp_path / "models",
        cache_dir=tmp_path / "cache",
        checkpoint_dir=checkpoint_dir,
        resume=True,
    )

    # Load compounds
    bindingdb = BindingDBSource()
    compounds = bindingdb.load_compounds(test_data)

    # Create manager
    manager = WebEnrichmentManager(config)

    # Mock client methods
    manager.swiss_client.process_compounds = Mock()
    manager.community_client.process_compounds = Mock()
    manager.social_client.process_compounds = Mock()

    # Enrich compounds
    manager.enrich_compounds(compounds)

    # Check only Amphetamine was processed
    processed_compounds = []
    for call in manager.swiss_client.process_compounds.call_args_list:
        processed_compounds.extend(call[0][0])
    
    assert len(processed_compounds) == 1
    assert processed_compounds[0].name == "Amphetamine"
