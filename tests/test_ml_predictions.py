"""Tests for machine learning prediction functionality."""

import json
import tempfile
from pathlib import Path
from unittest import mock

import numpy as np
import pytest
import torch
from rdkit import Chem

from binding_data_processor.config import Config
from binding_data_processor.processors.structure.ml.base import BaseMLModel
from binding_data_processor.processors.structure.ml.ensemble import EnsembleModel
from binding_data_processor.processors.structure.ml.feature_processing import (
    FeatureProcessor,
)
from binding_data_processor.processors.structure.ml.models.gnn import GNNModel
from binding_data_processor.processors.structure.ml.predictors.activity import (
    ActivityPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.affinity import (
    AffinityPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.toxicity import (
    ToxicityPredictor,
)


@pytest.fixture
def test_config():
    """Create test configuration."""
    with tempfile.TemporaryDirectory() as temp_dir:
        config = Config()
        config.model_dir = Path(temp_dir) / "models"
        config.cache_dir = Path(temp_dir) / "cache"

        # Create directories
        for directory in [config.model_dir, config.cache_dir]:
            directory.mkdir(parents=True)

        yield config


@pytest.fixture
def test_molecule():
    """Create test molecule."""
    return Chem.MolFromSmiles("CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC(=O)NC3=CC=C(C=C3)OC")


@pytest.fixture
def mock_model_state():
    """Create mock model state."""
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


@pytest.fixture
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


def test_feature_processing(test_config, test_molecule):
    """Test feature processing."""
    processor = FeatureProcessor(config=test_config)

    # Test Morgan fingerprint generation
    fingerprint = processor.generate_fingerprint(test_molecule)
    assert isinstance(fingerprint, np.ndarray)
    assert fingerprint.shape == (2048,)
    assert fingerprint.dtype == np.float32

    # Test graph features
    node_features, edge_features, edge_index = processor.generate_graph_features(
        test_molecule
    )
    assert isinstance(node_features, torch.Tensor)
    assert isinstance(edge_features, torch.Tensor)
    assert isinstance(edge_index, torch.Tensor)
    assert node_features.dim() == 2
    assert edge_features.dim() == 2
    assert edge_index.dim() == 2

    # Test caching
    cached_fingerprint = processor.generate_fingerprint(test_molecule)
    assert np.array_equal(fingerprint, cached_fingerprint)


def test_base_model(test_config, mock_model_state):
    """Test base ML model."""
    model = BaseMLModel(
        input_size=mock_model_state["input_size"],
        hidden_sizes=mock_model_state["hidden_sizes"],
        output_size=mock_model_state["output_size"],
        dropout=mock_model_state["dropout"],
    )

    # Save and load model
    model_path = test_config.model_dir / "test_model.pt"
    torch.save(mock_model_state, model_path)

    model.load_state_dict(torch.load(model_path)["state_dict"])
    model.eval()

    # Test forward pass
    x = torch.randn(1, mock_model_state["input_size"])
    output = model(x)
    assert output.shape == (1, mock_model_state["output_size"])


def test_gnn_model(test_config, mock_gnn_state):
    """Test GNN model."""
    model = GNNModel(
        node_features=mock_gnn_state["node_features"],
        edge_features=mock_gnn_state["edge_features"],
        hidden_size=mock_gnn_state["hidden_size"],
        num_layers=mock_gnn_state["num_layers"],
        dropout=mock_gnn_state["dropout"],
    )

    # Save and load model
    model_path = test_config.model_dir / "test_gnn.pt"
    torch.save(mock_gnn_state, model_path)

    model.load_state_dict(torch.load(model_path)["state_dict"])
    model.eval()

    # Test forward pass
    node_features = torch.randn(10, mock_gnn_state["node_features"])
    edge_features = torch.randn(15, mock_gnn_state["edge_features"])
    edge_index = torch.randint(0, 10, (2, 15))

    output = model(node_features, edge_features, edge_index)
    assert output.shape == (1, mock_gnn_state["hidden_size"])


def test_ensemble_model(test_config, mock_model_state):
    """Test ensemble model."""
    # Create ensemble of 3 models
    models = []
    for i in range(3):
        model = BaseMLModel(
            input_size=mock_model_state["input_size"],
            hidden_sizes=mock_model_state["hidden_sizes"],
            output_size=mock_model_state["output_size"],
            dropout=mock_model_state["dropout"],
        )
        model_path = test_config.model_dir / f"model_{i}.pt"
        torch.save(mock_model_state, model_path)
        model.load_state_dict(torch.load(model_path)["state_dict"])
        models.append(model)

    ensemble = EnsembleModel(models)
    ensemble.eval()

    # Test forward pass
    x = torch.randn(1, mock_model_state["input_size"])
    output, confidence = ensemble(x)
    assert output.shape == (1, mock_model_state["output_size"])
    assert confidence.shape == (1, 1)
    assert 0 <= confidence.item() <= 1


def test_activity_predictor(test_config, test_molecule, mock_model_state):
    """Test activity prediction."""
    predictor = ActivityPredictor(config=test_config)

    # Mock model loading
    with mock.patch("torch.load") as mock_load:
        mock_load.return_value = mock_model_state
        predictor.load_models()

    # Test prediction
    predictions = predictor.predict_compound(test_molecule)
    assert "target_predictions" in predictions
    assert len(predictions["target_predictions"]) > 0
    for pred in predictions["target_predictions"]:
        assert "target" in pred
        assert "probability" in pred
        assert "confidence" in pred
        assert 0 <= pred["probability"] <= 1
        assert 0 <= pred["confidence"] <= 1


def test_affinity_predictor(test_config, test_molecule, mock_gnn_state):
    """Test affinity prediction."""
    predictor = AffinityPredictor(config=test_config)

    # Mock model loading
    with mock.patch("torch.load") as mock_load:
        mock_load.return_value = mock_gnn_state
        predictor.load_models()

    # Test prediction
    predictions = predictor.predict_compound(test_molecule)
    assert "binding_predictions" in predictions
    assert len(predictions["binding_predictions"]) > 0
    for pred in predictions["binding_predictions"]:
        assert "target" in pred
        assert "affinity" in pred
        assert "confidence" in pred
        assert pred["affinity"] > 0
        assert 0 <= pred["confidence"] <= 1


def test_toxicity_predictor(test_config, test_molecule, mock_model_state):
    """Test toxicity prediction."""
    predictor = ToxicityPredictor(config=test_config)

    # Mock model loading
    with mock.patch("torch.load") as mock_load:
        mock_load.return_value = mock_model_state
        predictor.load_models()

    # Test prediction
    predictions = predictor.predict_compound(test_molecule)
    assert "toxicity_predictions" in predictions
    assert len(predictions["toxicity_predictions"]) > 0
    for pred in predictions["toxicity_predictions"]:
        assert "type" in pred
        assert "probability" in pred
        assert "confidence" in pred
        assert 0 <= pred["probability"] <= 1
        assert 0 <= pred["confidence"] <= 1


def test_model_caching(test_config, test_molecule, mock_model_state):
    """Test model prediction caching."""
    predictor = ActivityPredictor(config=test_config)

    # Mock model loading
    with mock.patch("torch.load") as mock_load:
        mock_load.return_value = mock_model_state
        predictor.load_models()

    # Initial prediction
    pred1 = predictor.predict_compound(test_molecule)

    # Cache should be used for subsequent predictions
    pred2 = predictor.predict_compound(test_molecule)
    assert pred1 == pred2

    # Check cache file
    cache_file = test_config.cache_dir / "activity_predictions.json"
    assert cache_file.exists()
    with open(cache_file) as f:
        cache = json.load(f)
        assert len(cache) == 1


def test_confidence_scoring(test_config, test_molecule, mock_model_state):
    """Test confidence score calculation."""
    predictor = ActivityPredictor(config=test_config)

    # Mock model loading with different predictions
    def mock_predict(x):
        return torch.tensor(
            [[0.8, 0.1, 0.1], [0.7, 0.2, 0.1], [0.9, 0.05, 0.05]]
        )  # High agreement

    with mock.patch("torch.load") as mock_load:
        mock_load.return_value = mock_model_state
        predictor.load_models()

        # Mock ensemble predictions
        predictor.model.models = [mock.MagicMock() for _ in range(3)]
        for model in predictor.model.models:
            model.forward.side_effect = mock_predict

        # Test prediction with high confidence
        predictions = predictor.predict_compound(test_molecule)
        for pred in predictions["target_predictions"]:
            assert pred["confidence"] > 0.8  # High confidence due to agreement

        # Mock divergent predictions
        def mock_predict_divergent(x):
            return torch.tensor(
                [[0.4, 0.3, 0.3], [0.3, 0.4, 0.3], [0.3, 0.3, 0.4]]
            )  # Low agreement

        for model in predictor.model.models:
            model.forward.side_effect = mock_predict_divergent

        # Test prediction with low confidence
        predictions = predictor.predict_compound(test_molecule)
        for pred in predictions["target_predictions"]:
            assert pred["confidence"] < 0.5  # Low confidence due to disagreement


def test_error_handling(test_config, test_molecule):
    """Test error handling in predictions."""
    predictor = ActivityPredictor(config=test_config)

    # Test missing model file
    with pytest.raises(FileNotFoundError):
        predictor.load_models()

    # Test invalid molecule
    with pytest.raises(ValueError):
        predictor.predict_compound(None)

    # Test model loading error
    with mock.patch("torch.load", side_effect=Exception("Model loading error")):
        with pytest.raises(Exception):
            predictor.load_models()

    # Test prediction error
    predictor.model = mock.MagicMock()
    predictor.model.forward.side_effect = Exception("Prediction error")
    with pytest.raises(Exception):
        predictor.predict_compound(test_molecule)
