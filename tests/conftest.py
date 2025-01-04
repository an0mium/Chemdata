"""Test configuration and shared fixtures."""

import json
import os
import tempfile
from pathlib import Path
from typing import Dict, Any

import numpy as np
import pandas as pd
import pytest
import torch
from rdkit import Chem

from binding_data_processor.config import Config


@pytest.fixture(scope="session")
def test_dir():
    """Create temporary test directory."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield Path(temp_dir)


@pytest.fixture(scope="session")
def test_config(test_dir):
    """Create test configuration."""
    config = Config()
    config.data_dir = test_dir / "data"
    config.model_dir = test_dir / "models"
    config.cache_dir = test_dir / "cache"
    config.log_dir = test_dir / "logs"

    # Create directories
    for directory in [
        config.data_dir,
        config.model_dir,
        config.cache_dir,
        config.log_dir,
    ]:
        directory.mkdir(parents=True)

    return config


@pytest.fixture(scope="session")
def test_smiles():
    """Create test SMILES strings."""
    return [
        "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC",
        "CC1=CC=C(C=C1)NC(=O)CCCN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)F",
    ]


@pytest.fixture(scope="session")
def test_molecules(test_smiles):
    """Create test RDKit molecules."""
    return [Chem.MolFromSmiles(smiles) for smiles in test_smiles]


@pytest.fixture(scope="session")
def test_fingerprints(test_molecules):
    """Create test Morgan fingerprints."""
    return [
        np.array(
            list(
                Chem.AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048).ToBitString()
            ),
            dtype=np.float32,
        )
        for mol in test_molecules
    ]


@pytest.fixture(scope="session")
def test_data():
    """Create test compound data."""
    return pd.DataFrame(
        {
            "name": ["Test Compound 1", "Test Compound 2"],
            "smiles": [
                "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC",
                "CC1=CC=C(C=C1)NC(=O)CCCN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)F",
            ],
            "molecular_weight": [395.47, 401.45],
            "logp": [3.2, 2.8],
            "target_type": ["serotonin", "dopamine"],
            "activity_type": ["Ki", "IC50"],
            "activity_value": [10.5, 25.3],
            "activity_unit": ["nM", "nM"],
            "reference": ["10.1021/jm123456", "10.1021/jm234567"],
        }
    )


@pytest.fixture(scope="session")
def test_predictions():
    """Create test prediction data."""
    return {
        "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC": {
            "toxicity": {
                "hepatotoxicity": {"probability": 0.2, "confidence": 0.8},
                "cardiotoxicity": {"probability": 0.1, "confidence": 0.9},
            },
            "abuse": {
                "potential": {"probability": 0.3, "confidence": 0.7},
                "dependence": {"probability": 0.2, "confidence": 0.8},
            },
            "activity": {
                "5-HT2A": {"probability": 0.9, "confidence": 0.8},
                "D2": {"probability": 0.3, "confidence": 0.7},
            },
        },
        "CC1=CC=C(C=C1)NC(=O)CCCN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)F": {
            "toxicity": {
                "hepatotoxicity": {"probability": 0.4, "confidence": 0.6},
                "cardiotoxicity": {"probability": 0.3, "confidence": 0.7},
            },
            "abuse": {
                "potential": {"probability": 0.5, "confidence": 0.7},
                "dependence": {"probability": 0.4, "confidence": 0.6},
            },
            "activity": {
                "5-HT2A": {"probability": 0.4, "confidence": 0.7},
                "D2": {"probability": 0.8, "confidence": 0.8},
            },
        },
    }


@pytest.fixture(scope="session")
def test_web_data():
    """Create test web data."""
    return {
        "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC": {
            "community": {
                "reports": [
                    {
                        "source": "erowid",
                        "url": "https://erowid.org/experiences/123",
                        "title": "Test Report 1",
                        "text": "Test experience report content 1",
                        "date": "2023-01-01",
                    }
                ],
                "safety": {
                    "warnings": ["May cause drowsiness"],
                    "interactions": ["MAOIs", "SSRIs"],
                    "contraindications": ["Heart conditions"],
                },
            },
            "social": {
                "reddit": [
                    {
                        "subreddit": "DrugNerds",
                        "title": "Test Post 1",
                        "text": "Test content 1",
                        "url": "https://reddit.com/r/DrugNerds/123",
                        "date": "2023-01-01",
                    }
                ],
                "twitter": [
                    {
                        "text": "Test tweet about compound",
                        "url": "https://twitter.com/user/123",
                        "date": "2023-02-01",
                    }
                ],
            },
        },
        "CC1=CC=C(C=C1)NC(=O)CCCN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)F": {
            "community": {
                "reports": [
                    {
                        "source": "psychonautwiki",
                        "url": "https://psychonautwiki.org/wiki/Test",
                        "title": "Test Report 2",
                        "text": "Test experience report content 2",
                        "date": "2023-02-01",
                    }
                ],
                "safety": {
                    "warnings": ["May cause dizziness"],
                    "interactions": ["Alcohol", "CNS depressants"],
                    "contraindications": ["Liver disease"],
                },
            },
            "social": {
                "reddit": [
                    {
                        "subreddit": "Nootropics",
                        "title": "Test Post 2",
                        "text": "Test content 2",
                        "url": "https://reddit.com/r/Nootropics/456",
                        "date": "2023-02-01",
                    }
                ],
                "twitter": [
                    {
                        "text": "Another test tweet",
                        "url": "https://twitter.com/user/456",
                        "date": "2023-03-01",
                    }
                ],
            },
        },
    }


@pytest.fixture(scope="session")
def mock_model_state():
    """Create mock ML model state."""
    return {
        "input_size": 2048,
        "hidden_sizes": [1024, 512, 256],
        "output_size": 128,
        "dropout": 0.2,
        "state_dict": {
            "fc1.weight": torch.randn(1024, 2048),
            "fc1.bias": torch.randn(1024),
            "fc2.weight": torch.randn(512, 1024),
            "fc2.bias": torch.randn(512),
            "fc3.weight": torch.randn(256, 512),
            "fc3.bias": torch.randn(256),
            "output.weight": torch.randn(128, 256),
            "output.bias": torch.randn(128),
        },
    }


@pytest.fixture(scope="session")
def mock_gnn_state():
    """Create mock GNN model state."""
    return {
        "node_features": 64,
        "edge_features": 32,
        "hidden_size": 128,
        "num_layers": 4,
        "dropout": 0.2,
        "state_dict": {
            "node_embedding.weight": torch.randn(64, 32),
            "edge_embedding.weight": torch.randn(32, 16),
            "gnn_layers.0.weight": torch.randn(128, 64),
            "gnn_layers.0.bias": torch.randn(128),
            "output.weight": torch.randn(128, 128),
            "output.bias": torch.randn(128),
        },
    }


def save_test_data(config: Config, data: Dict[str, Any]) -> None:
    """Save test data to files.

    Args:
        config: Test configuration
        data: Dictionary containing test data
    """
    # Save compound data
    if "compounds" in data:
        data["compounds"].to_csv(
            config.data_dir / "compounds.tsv", sep="\t", index=False
        )

    # Save predictions
    if "predictions" in data:
        with open(config.data_dir / "predictions.json", "w") as f:
            json.dump(data["predictions"], f, indent=2)

    # Save web data
    if "web_data" in data:
        with open(config.data_dir / "web_data.json", "w") as f:
            json.dump(data["web_data"], f, indent=2)

    # Save model states
    if "model_states" in data:
        for name, state in data["model_states"].items():
            torch.save(state, config.model_dir / f"{name}.pt")


def clean_test_data(config: Config) -> None:
    """Clean up test data files.

    Args:
        config: Test configuration
    """
    # Remove data files
    for file in config.data_dir.glob("*"):
        file.unlink()

    # Remove model files
    for file in config.model_dir.glob("*"):
        file.unlink()

    # Remove cache files
    for file in config.cache_dir.glob("*"):
        file.unlink()

    # Remove log files
    for file in config.log_dir.glob("*"):
        file.unlink()
