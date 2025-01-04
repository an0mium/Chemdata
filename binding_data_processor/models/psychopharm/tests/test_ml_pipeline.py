"""Tests for machine learning pipeline functionality."""

import pytest
import numpy as np
import pandas as pd
from unittest.mock import patch, Mock
from pathlib import Path

from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..ml import (
    MLPipeline,
    ModelConfig,
    PipelineConfig,
    PipelineStage,
    PipelineResult,
    DataLoader,
    DataPreprocessor,
    FeatureEngineer,
    ModelTrainer,
    ModelEvaluator,
    PredictionGenerator,
    CheckpointManager,
    ValidationConfig,
    EnsembleConfig,
)


@pytest.fixture
def test_compounds():
    """Create test compounds fixture."""
    compounds = []
    
    # Caffeine
    caffeine = PsychoactiveCompound(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    caffeine.psychoactive_class = PsychoactiveClass.STIMULANT
    caffeine.add_receptor_binding(
        "A2A",
        affinity=0.8,
        confidence=0.95,
        activity="antagonist"
    )
    caffeine.add_receptor_binding(
        "A1",
        affinity=0.7,
        confidence=0.9,
        activity="antagonist"
    )
    caffeine.effect_profile = {
        "stimulation": (0.8, 0.9),
        "focus": (0.7, 0.8),
        "wakefulness": (0.9, 0.95),
    }
    caffeine.safety_alerts = {
        "anxiety": RiskLevel.MODERATE,
        "insomnia": RiskLevel.HIGH,
        "tachycardia": RiskLevel.LOW,
    }
    compounds.append(caffeine)
    
    # Amphetamine
    amphetamine = PsychoactiveCompound(
        name="Amphetamine",
        smiles="CC(N)CC1=CC=CC=C1",
        cas_number="300-62-9",
    )
    amphetamine.psychoactive_class = PsychoactiveClass.STIMULANT
    amphetamine.add_receptor_binding(
        "DAT",
        affinity=0.05,
        confidence=0.95,
        activity="inhibitor"
    )
    amphetamine.add_receptor_binding(
        "NET",
        affinity=0.1,
        confidence=0.9,
        activity="inhibitor"
    )
    amphetamine.effect_profile = {
        "stimulation": (0.9, 0.95),
        "euphoria": (0.8, 0.9),
        "focus": (0.85, 0.9),
    }
    amphetamine.safety_alerts = {
        "addiction": RiskLevel.HIGH,
        "cardiovascular": RiskLevel.HIGH,
        "neurotoxicity": RiskLevel.MODERATE,
    }
    compounds.append(amphetamine)
    
    return compounds


@pytest.fixture
def pipeline_config():
    """Create test pipeline configuration fixture."""
    return PipelineConfig(
        stages=[
            PipelineStage.PREPROCESS,
            PipelineStage.FEATURE_ENGINEERING,
            PipelineStage.MODEL_TRAINING,
            PipelineStage.MODEL_EVALUATION,
            PipelineStage.PREDICTION,
        ],
        input_format="tsv",
        output_format="tsv",
        model_config=ModelConfig(
            model_type="ensemble",
            architectures=["gnn", "transformer"],
            feature_types=["structural", "physicochemical"],
            training={"epochs": 100, "batch_size": 32},
            validation={"split": 0.2, "k_folds": 5}
        ),
        ensemble_config=EnsembleConfig(
            num_models=5,
            voting_method="soft",
            diversity_metric="correlation",
        ),
        validation_config=ValidationConfig(
            test_size=0.2,
            cv_folds=5,
            metrics=["accuracy", "f1", "auroc"],
        ),
        parallel=True,
        n_jobs=2,
        checkpoint_interval=10,
        progress_callback=lambda x: None,
        training_dir=Path("training"),
        model_dir=Path("models"),
        cache_dir=Path("cache"),
    )


@pytest.fixture
def test_data(tmp_path, test_compounds):
    """Create test input data fixture."""
    # Create test TSV file
    input_file = tmp_path / "input.tsv"
    data = pd.DataFrame({
        "name": [c.name for c in test_compounds],
        "smiles": [c.smiles for c in test_compounds],
        "cas_number": [c.cas_number for c in test_compounds]
    })
    data.to_csv(input_file, sep="\t", index=False)
    return input_file


@pytest.fixture
def ml_pipeline(pipeline_config):
    """Create test ML pipeline fixture."""
    return MLPipeline(config=pipeline_config)


class TestPipelineBasics:
    """Tests for basic pipeline functionality."""

    def test_initialization(self, ml_pipeline, pipeline_config):
        """Test initialization of MLPipeline."""
        assert isinstance(ml_pipeline.loader, DataLoader)
        assert isinstance(ml_pipeline.preprocessor, DataPreprocessor)
        assert isinstance(ml_pipeline.feature_engineer, FeatureEngineer)
        assert isinstance(ml_pipeline.trainer, ModelTrainer)
        assert isinstance(ml_pipeline.evaluator, ModelEvaluator)
        assert isinstance(ml_pipeline.predictor, PredictionGenerator)
        assert isinstance(ml_pipeline.checkpoint_manager, CheckpointManager)
        assert ml_pipeline.config == pipeline_config
        assert ml_pipeline.stats == {}

    def test_load_data(self, ml_pipeline, test_data):
        """Test data loading."""
        # Load data
        result = ml_pipeline.load_data(test_data)
        
        # Check result
        assert isinstance(result, PipelineResult)
        assert result.success
        assert len(result.compounds) == 2
        assert all(isinstance(c, PsychoactiveCompound) for c in result.compounds)
        assert "data_loading" in ml_pipeline.stats

    def test_preprocess_data(self, ml_pipeline, test_compounds):
        """Test data preprocessing."""
        # Preprocess data
        result = ml_pipeline.preprocess_data(test_compounds)
        
        # Check result
        assert isinstance(result, PipelineResult)
        assert result.success
        assert len(result.compounds) == len(test_compounds)
        assert "preprocessing" in ml_pipeline.stats
        assert "cleaned_compounds" in result.data
        assert "standardized_compounds" in result.data


class TestFeatureEngineering:
    """Tests for feature engineering."""

    @patch('rdkit.Chem.AllChem.GetMorganFingerprintAsBitVect')
    def test_feature_engineering(self, mock_fingerprint, ml_pipeline, test_compounds):
        """Test feature engineering."""
        # Mock fingerprint generation
        mock_fingerprint.return_value = np.zeros(2048, dtype=np.int8)
        
        # Engineer features
        result = ml_pipeline.engineer_features(test_compounds)
        
        # Check result
        assert isinstance(result, PipelineResult)
        assert result.success
        assert isinstance(result.features, np.ndarray)
        assert result.features.shape[0] == len(test_compounds)
        assert not np.any(np.isnan(result.features))
        assert "feature_engineering" in ml_pipeline.stats
        assert "structural_features" in result.data
        assert "physicochemical_features" in result.data

    def test_molecular_features(self, test_compounds, ml_pipeline):
        """Test molecular feature extraction."""
        # Extract features
        features = ml_pipeline.feature_engineer.extract_molecular_features(
            compounds=test_compounds
        )
        
        # Check features
        assert isinstance(features, np.ndarray)
        assert features.shape[0] == len(test_compounds)
        assert features.shape[1] > 0  # Multiple features per compound
        assert not np.any(np.isnan(features))  # No missing values

    def test_graph_features(self, test_compounds, ml_pipeline):
        """Test molecular graph feature extraction."""
        # Extract features
        features = ml_pipeline.feature_engineer.extract_graph_features(
            compounds=test_compounds
        )
        
        # Check features
        assert isinstance(features, list)
        assert len(features) == len(test_compounds)
        assert all(isinstance(f, dict) for f in features)
        assert all("nodes" in f and "edges" in f for f in features)

    def test_pharmacophore_features(self, test_compounds, ml_pipeline):
        """Test pharmacophore feature extraction."""
        # Extract features
        features = ml_pipeline.feature_engineer.extract_pharmacophore_features(
            compounds=test_compounds
        )
        
        # Check features
        assert isinstance(features, np.ndarray)
        assert features.shape[0] == len(test_compounds)
        assert features.shape[1] > 0
        assert all(f in ["aromatic", "hbd", "hba"] for f in features.dtype.names)


class TestModelTraining:
    """Tests for model training."""

    def test_train_binding_model(self, test_compounds, ml_pipeline):
        """Test binding affinity model training."""
        # Mock training data
        X = np.random.randn(100, 256)  # Features
        y = np.random.randn(100)  # Binding affinities
        
        with patch.object(ml_pipeline.trainer, "_prepare_data") as mock_prep:
            mock_prep.return_value = (X, y)
            
            # Train model
            result = ml_pipeline.trainer.train_binding_model(
                compounds=test_compounds,
                target="A2A"
            )
            
            # Check result
            assert result.success
            assert result.model is not None
            assert result.metrics["train_loss"] < result.metrics["initial_loss"]
            assert 0 <= result.metrics["val_accuracy"] <= 1

    def test_train_activity_model(self, test_compounds, ml_pipeline):
        """Test activity prediction model training."""
        # Mock training data
        X = np.random.randn(100, 256)  # Features
        y = np.random.randint(0, 2, 100)  # Binary activity labels
        
        with patch.object(ml_pipeline.trainer, "_prepare_data") as mock_prep:
            mock_prep.return_value = (X, y)
            
            # Train model
            result = ml_pipeline.trainer.train_activity_model(
                compounds=test_compounds,
                activity_type="stimulant"
            )
            
            # Check result
            assert result.success
            assert result.model is not None
            assert result.metrics["train_loss"] < result.metrics["initial_loss"]
            assert 0 <= result.metrics["val_accuracy"] <= 1

    def test_train_safety_model(self, test_compounds, ml_pipeline):
        """Test safety prediction model training."""
        # Mock training data
        X = np.random.randn(100, 256)  # Features
        y = np.random.randint(0, 3, 100)  # Risk level labels
        
        with patch.object(ml_pipeline.trainer, "_prepare_data") as mock_prep:
            mock_prep.return_value = (X, y)
            
            # Train model
            result = ml_pipeline.trainer.train_safety_model(
                compounds=test_compounds,
                risk_type="addiction"
            )
            
            # Check result
            assert result.success
            assert result.model is not None
            assert result.metrics["train_loss"] < result.metrics["initial_loss"]
            assert 0 <= result.metrics["val_accuracy"] <= 1


class TestModelEvaluation:
    """Tests for model evaluation."""

    def test_cross_validation(self, test_compounds, ml_pipeline):
        """Test cross-validation evaluation."""
        # Mock data and model
        X = np.random.randn(100, 256)
        y = np.random.randn(100)
        mock_model = Mock()
        mock_model.predict.return_value = np.random.randn(20)  # Predictions for each fold
        
        with patch.object(ml_pipeline.evaluator, "_get_data") as mock_get:
            mock_get.return_value = (X, y)
            
            # Evaluate model
            result = ml_pipeline.evaluator.cross_validate(
                model=mock_model,
                compounds=test_compounds,
                target="binding"
            )
            
            # Check result
            assert result.success
            assert len(result.fold_metrics) == ml_pipeline.config.validation_config.cv_folds
            assert all(0 <= m["accuracy"] <= 1 for m in result.fold_metrics)
            assert "mean_accuracy" in result.aggregate_metrics
            assert "std_accuracy" in result.aggregate_metrics

    def test_holdout_validation(self, test_compounds, ml_pipeline):
        """Test holdout set validation."""
        # Mock data and model
        X = np.random.randn(100, 256)
        y = np.random.randn(100)
        mock_model = Mock()
        mock_model.predict.return_value = np.random.randn(20)
        
        with patch.object(ml_pipeline.evaluator, "_get_data") as mock_get:
            mock_get.return_value = (X, y)
            
            # Evaluate model
            result = ml_pipeline.evaluator.validate_holdout(
                model=mock_model,
                compounds=test_compounds,
                target="activity"
            )
            
            # Check result
            assert result.success
            assert "test_accuracy" in result.metrics
            assert "test_f1" in result.metrics
            assert "test_auroc" in result.metrics
            assert all(0 <= v <= 1 for v in result.metrics.values())


class TestModelEnsemble:
    """Tests for model ensemble."""

    def test_ensemble_training(self, test_compounds, ml_pipeline):
        """Test ensemble model training."""
        # Mock base models
        mock_models = [Mock() for _ in range(ml_pipeline.config.ensemble_config.num_models)]
        for model in mock_models:
            model.predict.return_value = np.random.randn(len(test_compounds))
        
        with patch.object(ml_pipeline.trainer, "_train_base_model") as mock_train:
            mock_train.side_effect = mock_models
            
            # Train ensemble
            result = ml_pipeline.trainer.train_ensemble(
                compounds=test_compounds,
                target="binding"
            )
            
            # Check result
            assert result.success
            assert len(result.models) == ml_pipeline.config.ensemble_config.num_models
            assert result.diversity_score > 0
            assert all(m is not None for m in result.models)

    def test_ensemble_prediction(self, test_compounds, ml_pipeline):
        """Test ensemble prediction."""
        # Mock ensemble models
        mock_models = [Mock() for _ in range(ml_pipeline.config.ensemble_config.num_models)]
        for model in mock_models:
            model.predict.return_value = np.random.randn(len(test_compounds))
        
        with patch.object(ml_pipeline.predictor, "models", mock_models):
            # Get ensemble predictions
            result = ml_pipeline.predictor.predict_ensemble(
                compounds=test_compounds
            )
            
            # Check result
            assert result.success
            assert len(result.predictions) == len(test_compounds)
            assert len(result.uncertainties) == len(test_compounds)
            assert all(0 <= u <= 1 for u in result.uncertainties)


class TestUncertaintyEstimation:
    """Tests for uncertainty estimation."""

    def test_ensemble_uncertainty(self, test_compounds, ml_pipeline):
        """Test ensemble-based uncertainty estimation."""
        # Mock ensemble predictions
        mock_predictions = [
            np.random.randn(len(test_compounds))
            for _ in range(ml_pipeline.config.ensemble_config.num_models)
        ]
        
        # Calculate uncertainty
        uncertainties = ml_pipeline.predictor.estimate_ensemble_uncertainty(
            predictions=mock_predictions
        )
        
        # Check uncertainties
        assert len(uncertainties) == len(test_compounds)
        assert all(0 <= u <= 1 for u in uncertainties)
        assert any(u > 0 for u in uncertainties)  # Some uncertainty present

    def test_dropout_uncertainty(self, test_compounds, ml_pipeline):
        """Test dropout-based uncertainty estimation."""
        # Mock model with dropout
        mock_model = Mock()
        mock_model.predict.side_effect = [
            np.random.randn(len(test_compounds))
            for _ in range(10)  # 10 forward passes
        ]
        
        # Calculate uncertainty
        uncertainties = ml_pipeline.predictor.estimate_dropout_uncertainty(
            model=mock_model,
            compounds=test_compounds,
            num_passes=10
        )
        
        # Check uncertainties
        assert len(uncertainties) == len(test_compounds)
        assert all(0 <= u <= 1 for u in uncertainties)
        assert any(u > 0 for u in uncertainties)  # Some uncertainty present


class TestPipelineIntegration:
    """Tests for full pipeline integration."""

    def test_run_pipeline(self, ml_pipeline, test_data, tmp_path):
        """Test full pipeline execution."""
        # Configure output
        output_file = tmp_path / "output.tsv"
        
        # Run pipeline
        result = ml_pipeline.run(
            input_file=test_data,
            output_file=output_file
        )
        
        # Check result
        assert isinstance(result, PipelineResult)
        assert result.success
        assert output_file.exists()
        assert "pipeline_time" in ml_pipeline.stats
        completed_stages = result.data["completed_stages"]
        assert all(stage in completed_stages for stage in ml_pipeline.config.stages)

    def test_checkpointing(self, ml_pipeline, test_data, tmp_path):
        """Test pipeline checkpointing."""
        # Configure checkpointing
        checkpoint_dir = tmp_path / "checkpoints"
        checkpoint_dir.mkdir()
        ml_pipeline.config.checkpoint_dir = checkpoint_dir
        ml_pipeline.config.checkpoint_interval = 1
        
        # Run pipeline
        result = ml_pipeline.run(input_file=test_data)
        
        # Check checkpoints
        assert result.success
        assert checkpoint_dir.exists()
        checkpoints = list(checkpoint_dir.glob("*.pkl"))
        assert len(checkpoints) > 0
        
        # Test checkpoint recovery
        recovered_pipeline = MLPipeline(config=ml_pipeline.config)
        recovery_result = recovered_pipeline.recover_from_checkpoint(
            checkpoint_file=checkpoints[-1]
        )
        assert recovery_result.success
        assert len(recovery_result.compounds) == 2

    def test_error_handling(self, ml_pipeline, test_data):
        """Test pipeline error handling."""
        # Configure invalid settings
        ml_pipeline.config.model_config.architectures = ["invalid_model"]
        
        # Run pipeline
        result = ml_pipeline.run(input_file=test_data)
        
        # Check error handling
        assert not result.success
        assert "invalid model architecture" in str(result.error).lower()
        assert "failed_stages" in ml_pipeline.stats

    def test_progress_monitoring(self, ml_pipeline, test_data):
        """Test progress monitoring."""
        # Configure progress monitoring
        progress = []
        ml_pipeline.config.progress_callback = lambda x: progress.append(x)
        
        # Run pipeline
        result = ml_pipeline.run(input_file=test_data)
        
        # Check progress monitoring
        assert result.success
        assert len(progress) > 0
        assert progress[-1] == 100
        assert "stage_progress" in ml_pipeline.stats


if __name__ == "__main__":
    pytest.main([__file__])
