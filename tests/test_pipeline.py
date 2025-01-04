"""Tests for the data processing pipeline."""

import os
import tempfile
from pathlib import Path
from unittest import mock

import pandas as pd
import pytest

from binding_data_processor.config import Config
from binding_data_processor.pipeline import PipelineManager


@pytest.fixture
def test_config():
    """Create test configuration."""
    with tempfile.TemporaryDirectory() as temp_dir:
        config = Config()
        config.data_dir = Path(temp_dir) / "data"
        config.model_dir = Path(temp_dir) / "models"
        config.cache_dir = Path(temp_dir) / "cache"
        config.log_dir = Path(temp_dir) / "logs"

        # Create directories
        for directory in [
            config.data_dir,
            config.model_dir,
            config.cache_dir,
            config.log_dir,
        ]:
            directory.mkdir(parents=True)

        yield config


@pytest.fixture
def test_data():
    """Create test BindingDB data."""
    data = pd.DataFrame(
        {
            "Ligand SMILES": [
                "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC",
                "CC1=CC=C(C=C1)NC(=O)CCCN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)F",
            ],
            "Target Name": [
                "5-HT2A receptor",
                "Dopamine D2 receptor",
            ],
            "Activity Type": [
                "Ki",
                "IC50",
            ],
            "Activity Value": [
                "10.5",
                "25.3",
            ],
            "Activity Units": [
                "nM",
                "nM",
            ],
            "Activity Relation": [
                "=",
                "=",
            ],
            "Reference": [
                "10.1021/jm123456",
                "10.1021/jm234567",
            ],
        }
    )
    return data


@pytest.fixture
def mock_web_data():
    """Create mock web data."""
    return {
        "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC": {
            "pubchem": {
                "name": "Test Compound 1",
                "molecular_weight": 395.47,
                "logp": 3.2,
            },
            "community": {
                "reports": [
                    {
                        "source": "erowid",
                        "text": "Test report 1",
                    }
                ],
            },
        },
        "CC1=CC=C(C=C1)NC(=O)CCCN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)F": {
            "pubchem": {
                "name": "Test Compound 2",
                "molecular_weight": 401.45,
                "logp": 2.8,
            },
            "community": {
                "reports": [
                    {
                        "source": "psychonautwiki",
                        "text": "Test report 2",
                    }
                ],
            },
        },
    }


@pytest.fixture
def mock_predictions():
    """Create mock ML predictions."""
    return {
        "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC": {
            "toxicity": {
                "hepatotoxicity": {
                    "probability": 0.2,
                    "confidence": 0.8,
                },
            },
            "abuse": {
                "potential": {
                    "probability": 0.3,
                    "confidence": 0.7,
                },
            },
        },
        "CC1=CC=C(C=C1)NC(=O)CCCN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)F": {
            "toxicity": {
                "hepatotoxicity": {
                    "probability": 0.4,
                    "confidence": 0.6,
                },
            },
            "abuse": {
                "potential": {
                    "probability": 0.5,
                    "confidence": 0.7,
                },
            },
        },
    }


def test_pipeline_initialization(test_config):
    """Test pipeline initialization."""
    pipeline = PipelineManager(
        data_dir=str(test_config.data_dir),
        model_dir=str(test_config.model_dir),
        n_workers=1,
        batch_size=10,
    )

    assert pipeline.data_dir == test_config.data_dir
    assert pipeline.model_dir == test_config.model_dir
    assert pipeline.n_workers == 1
    assert pipeline.batch_size == 10


def test_pipeline_processing(test_config, test_data, mock_web_data, mock_predictions):
    """Test pipeline data processing."""
    # Save test data
    test_data_file = test_config.data_dir / "test_bindingdb.tsv"
    test_data.to_csv(test_data_file, sep="\t", index=False)

    # Initialize pipeline
    pipeline = PipelineManager(
        data_dir=str(test_config.data_dir),
        model_dir=str(test_config.model_dir),
        n_workers=1,
        batch_size=10,
    )

    # Mock web data enrichment
    with mock.patch(
        "binding_data_processor.pipeline.WebEnrichmentProcessor.enrich_compound",
        side_effect=lambda smiles: mock_web_data.get(smiles, {}),
    ):
        # Mock ML predictions
        with mock.patch(
            "binding_data_processor.pipeline.MLPredictor.predict_compound",
            side_effect=lambda smiles: mock_predictions.get(smiles, {}),
        ):
            # Run pipeline
            stats = pipeline.run_pipeline(
                bindingdb_file=str(test_data_file),
                target_patterns={
                    "serotonin": r"5-HT\d*[A-Z]?",
                    "dopamine": r"D\d+",
                },
                skip_predictions=False,
                skip_web_data=False,
                output_dir=str(test_config.data_dir / "output"),
            )

    # Check statistics
    assert stats["total_compounds"] == 2
    assert stats["bindingdb_compounds"] == 2
    assert stats["web_compounds"] == 2
    assert stats["with_predictions"] == 2
    assert stats["with_web_data"] == 2
    assert not stats["errors"]

    # Check output files
    output_dir = test_config.data_dir / "output"
    assert (output_dir / "compounds.tsv").exists()
    assert (output_dir / "activities.tsv").exists()
    assert (output_dir / "predictions.json").exists()
    assert (output_dir / "web_data.json").exists()

    # Check compound data
    compounds = pd.read_csv(output_dir / "compounds.tsv", sep="\t")
    assert len(compounds) == 2
    assert "molecular_weight" in compounds.columns
    assert "logp" in compounds.columns

    # Check activity data
    activities = pd.read_csv(output_dir / "activities.tsv", sep="\t")
    assert len(activities) == 2
    assert "target_type" in activities.columns
    assert "activity_type" in activities.columns
    assert "activity_value" in activities.columns


def test_pipeline_error_handling(test_config):
    """Test pipeline error handling."""
    pipeline = PipelineManager(
        data_dir=str(test_config.data_dir),
        model_dir=str(test_config.model_dir),
        n_workers=1,
        batch_size=10,
    )

    # Test missing input file
    with pytest.raises(FileNotFoundError):
        pipeline.run_pipeline(
            bindingdb_file="nonexistent.tsv",
            output_dir=str(test_config.data_dir / "output"),
        )

    # Test invalid target patterns
    with pytest.raises(ValueError):
        pipeline.run_pipeline(
            bindingdb_file=str(test_config.data_dir / "test.tsv"),
            target_patterns={"invalid": "["},  # Invalid regex
            output_dir=str(test_config.data_dir / "output"),
        )

    # Test invalid output directory
    with pytest.raises(NotADirectoryError):
        pipeline.run_pipeline(
            bindingdb_file=str(test_config.data_dir / "test.tsv"),
            output_dir="/nonexistent/directory",
        )
