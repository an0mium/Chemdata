"""Tests for machine learning model functionality."""

import pytest
import numpy as np
from unittest.mock import patch, Mock

from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..ml import (
    MLModel,
    BindingPredictor,
    ActivityPredictor,
    SafetyPredictor,
    ModelConfig,
    ModelResult,
    ModelEnsemble,
    FeatureExtractor,
    ModelTrainer,
    ModelValidator,
    ModelEvaluator,
    UncertaintyEstimator,
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
    caffeine.effect_profile = {
        "stimulation": (0.8, 0.9),
        "focus": (0.7, 0.8),
    }
    caffeine.safety_alerts = {
        "anxiety": RiskLevel.MODERATE,
        "insomnia": RiskLevel.HIGH,
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
    amphetamine.effect_profile = {
        "stimulation": (0.9, 0.95),
        "euphoria": (0.8, 0.9),
    }
    amphetamine.safety_alerts = {
        "addiction": RiskLevel.HIGH,
        "cardiovascular": RiskLevel.HIGH,
    }
    compounds.append(amphetamine)
    
    return compounds


@pytest.fixture
def model_config():
    """Create test model configuration fixture."""
    return ModelConfig(
        model_type="ensemble",
        architectures=["gnn", "transformer"],
        feature_types=["structural", "physicochemical"],
        training={"epochs": 100, "batch_size": 32},
        validation={"split": 0.2, "k_folds": 5}
    )


@pytest.fixture
def feature_extractor():
    """Create test feature extractor fixture."""
    return FeatureExtractor()


@pytest.fixture
def model_trainer():
    """Create test model trainer fixture."""
    return ModelTrainer()


@pytest.fixture
def model_validator():
    """Create test model validator fixture."""
    return ModelValidator()


@pytest.fixture
def model_evaluator():
    """Create test model evaluator fixture."""
    return ModelEvaluator()


@pytest.fixture
def model_ensemble():
    """Create test model ensemble fixture."""
    return ModelEnsemble(
        models=[
            MLModel(name="model1"),
            MLModel(name="model2"),
            MLModel(name="model3"),
        ]
    )


@pytest.fixture
def uncertainty_estimator():
    """Create test uncertainty estimator fixture."""
    return UncertaintyEstimator()


class TestFeatureExtractor:
    """Tests for FeatureExtractor class."""

    def test_initialization(self, feature_extractor):
        """Test initialization of FeatureExtractor."""
        assert feature_extractor.fingerprint_size == 2048
        assert feature_extractor.include_3d is False
        assert feature_extractor.stats == {}

    def test_structural_features(self, feature_extractor, test_compounds):
        """Test structural feature extraction."""
        features = feature_extractor.extract_structural_features(test_compounds[0])
        assert isinstance(features, np.ndarray)
        assert features.shape[1] >= 100  # Minimum feature dimension
        assert not np.any(np.isnan(features))  # No NaN values

    def test_physicochemical_features(self, feature_extractor, test_compounds):
        """Test physicochemical feature extraction."""
        features = feature_extractor.extract_physicochemical_features(test_compounds[0])
        assert isinstance(features, np.ndarray)
        assert features.shape[1] >= 50  # Minimum feature dimension
        assert not np.any(np.isnan(features))  # No NaN values

    def test_combined_features(self, feature_extractor, test_compounds):
        """Test combined feature extraction."""
        features = feature_extractor.extract_features(
            test_compounds[0],
            feature_types=["structural", "physicochemical"]
        )
        assert isinstance(features, np.ndarray)
        assert features.shape[1] >= 150  # Combined feature dimensions
        assert not np.any(np.isnan(features))  # No NaN values

    def test_3d_features(self, feature_extractor, test_compounds):
        """Test 3D feature extraction."""
        # Enable 3D features
        feature_extractor.include_3d = True
        
        # Extract features
        features = feature_extractor.extract_features(test_compounds)
        
        # Check 3D features
        assert features.shape[1] > feature_extractor.fingerprint_size
        assert "3d_features" in feature_extractor.stats

    def test_feature_importance(self, feature_extractor, test_compounds):
        """Test feature importance analysis."""
        # Extract features and analyze importance
        features = feature_extractor.extract_features(test_compounds)
        importance = feature_extractor.analyze_feature_importance(features)
        
        # Check importance scores
        assert isinstance(importance, np.ndarray)
        assert len(importance) == features.shape[1]
        assert np.all(importance >= 0)
        assert np.isclose(np.sum(importance), 1.0)


class TestModelTrainer:
    """Tests for ModelTrainer class."""

    def test_initialization(self, model_trainer):
        """Test initialization of ModelTrainer."""
        assert model_trainer.stats == {}

    @patch.object(MLModel, "train")
    def test_train_model(self, mock_train, model_trainer, test_compounds, model_config):
        """Test model training."""
        # Mock training
        mock_train.return_value = Mock(score=0.85)
        
        # Train model
        result = model_trainer.train(test_compounds, config=model_config)
        
        # Check result
        assert result.success
        assert result.metrics["score"] == 0.85
        assert "training_time" in result.stats

    def test_train_gnn(self, model_trainer, test_compounds, model_config):
        """Test GNN model training."""
        # Prepare data
        X = np.random.randn(len(test_compounds), 100)  # Mock features
        y = np.array([0.8, 0.05])  # Mock binding affinities
        
        # Train model
        model = model_trainer.train_gnn(X=X, y=y, config=model_config)
        
        # Check model
        assert model is not None
        assert hasattr(model, "predict")
        
        # Test prediction
        pred = model.predict(X)
        assert pred.shape == (len(test_compounds),)
        assert np.all((pred >= 0) & (pred <= 1))  # Valid probabilities

    def test_train_transformer(self, model_trainer, test_compounds, model_config):
        """Test transformer model training."""
        # Prepare data
        X = np.random.randn(len(test_compounds), 100)  # Mock features
        y = np.array([0.8, 0.05])  # Mock binding affinities
        
        # Train model
        model = model_trainer.train_transformer(X=X, y=y, config=model_config)
        
        # Check model
        assert model is not None
        assert hasattr(model, "predict")
        
        # Test prediction
        pred = model.predict(X)
        assert pred.shape == (len(test_compounds),)
        assert np.all((pred >= 0) & (pred <= 1))  # Valid probabilities

    def test_train_ensemble(self, model_trainer, test_compounds, model_config):
        """Test ensemble model training."""
        # Prepare data
        X = np.random.randn(len(test_compounds), 100)  # Mock features
        y = np.array([0.8, 0.05])  # Mock binding affinities
        
        # Train ensemble
        ensemble = model_trainer.train_ensemble(
            X=X,
            y=y,
            config=model_config
        )
        
        # Check ensemble
        assert isinstance(ensemble, ModelEnsemble)
        assert len(ensemble.models) == len(model_config.architectures)
        
        # Test prediction
        pred = ensemble.predict(X)
        assert pred.shape == (len(test_compounds),)
        assert np.all((pred >= 0) & (pred <= 1))  # Valid probabilities
        
        # Check uncertainty estimation
        uncertainty = ensemble.estimate_uncertainty(X)
        assert uncertainty.shape == (len(test_compounds),)
        assert np.all(uncertainty >= 0)  # Valid uncertainty values


class TestModelValidator:
    """Tests for ModelValidator class."""

    def test_initialization(self, model_validator):
        """Test initialization of ModelValidator."""
        assert model_validator.stats == {}

    def test_cross_validation(self, model_validator, test_compounds, model_config):
        """Test cross-validation."""
        # Prepare data
        X = np.random.randn(len(test_compounds), 100)  # Mock features
        y = np.array([0.8, 0.05])  # Mock binding affinities
        
        # Run cross-validation
        cv_results = model_validator.cross_validate(X=X, y=y, config=model_config)
        
        # Check results
        assert "mean_score" in cv_results
        assert "std_score" in cv_results
        assert "fold_scores" in cv_results
        assert len(cv_results["fold_scores"]) == model_config.validation["k_folds"]

    def test_validation_metrics(self, model_validator, test_compounds):
        """Test validation metrics."""
        # Prepare data
        y_true = np.array([0.8, 0.05])
        y_pred = np.array([0.75, 0.1])
        
        # Calculate metrics
        metrics = model_validator.calculate_metrics(y_true=y_true, y_pred=y_pred)
        
        # Check metrics
        assert "mse" in metrics
        assert "mae" in metrics
        assert "r2" in metrics
        assert all(v >= 0 for v in metrics.values())  # Valid metric values


class TestBindingPredictor:
    """Tests for BindingPredictor class."""

    def test_initialization(self, binding_predictor):
        """Test initialization of BindingPredictor."""
        assert isinstance(binding_predictor.model, MLModel)
        assert binding_predictor.receptors == ["A2A", "DAT", "NET"]
        assert binding_predictor.stats == {}

    def test_predict_binding(self, test_compounds, model_config):
        """Test binding affinity prediction."""
        predictor = BindingPredictor(config=model_config)
        
        # Make prediction
        result = predictor.predict_binding(test_compounds[0])
        
        # Check result
        assert isinstance(result, ModelResult)
        assert result.success
        assert "binding_affinity" in result.predictions
        assert "confidence" in result.predictions
        assert 0 <= result.predictions["binding_affinity"] <= 1
        assert 0 <= result.predictions["confidence"] <= 1

    def test_predict_selectivity(self, test_compounds, model_config):
        """Test binding selectivity prediction."""
        predictor = BindingPredictor(config=model_config)
        
        # Make prediction
        result = predictor.predict_selectivity(test_compounds[0])
        
        # Check result
        assert isinstance(result, ModelResult)
        assert result.success
        assert "selectivity_score" in result.predictions
        assert "target_preferences" in result.predictions
        assert 0 <= result.predictions["selectivity_score"] <= 1

    def test_uncertainty_estimation(self, binding_predictor, test_compounds):
        """Test prediction uncertainty estimation."""
        # Make predictions with uncertainty
        predictions = binding_predictor.predict_binding(
            test_compounds,
            estimate_uncertainty=True
        )
        
        # Check uncertainty estimates
        for compound_id, receptor_preds in predictions.items():
            for receptor, pred in receptor_preds.items():
                assert isinstance(pred, tuple)  # (mean, std)
                assert 0 <= pred[0] <= 1  # mean
                assert pred[1] >= 0  # std


class TestActivityPredictor:
    """Tests for ActivityPredictor class."""

    def test_initialization(self, activity_predictor):
        """Test initialization of ActivityPredictor."""
        assert isinstance(activity_predictor.model, MLModel)
        assert activity_predictor.effects == ["stimulation", "euphoria", "focus"]
        assert activity_predictor.stats == {}

    def test_predict_activity(self, test_compounds, model_config):
        """Test activity prediction."""
        predictor = ActivityPredictor(config=model_config)
        
        # Make prediction
        result = predictor.predict_activity(test_compounds[0])
        
        # Check result
        assert isinstance(result, ModelResult)
        assert result.success
        assert "activity_score" in result.predictions
        assert "mechanism" in result.predictions
        assert 0 <= result.predictions["activity_score"] <= 1

    def test_predict_effects(self, test_compounds, model_config):
        """Test effect prediction."""
        predictor = ActivityPredictor(config=model_config)
        
        # Make prediction
        result = predictor.predict_effects(test_compounds[0])
        
        # Check result
        assert isinstance(result, ModelResult)
        assert result.success
        assert "effects" in result.predictions
        assert all(0 <= v <= 1 for v in result.predictions["effects"].values())


class TestSafetyPredictor:
    """Tests for SafetyPredictor class."""

    def test_initialization(self, safety_predictor):
        """Test initialization of SafetyPredictor."""
        assert isinstance(safety_predictor.model, MLModel)
        assert safety_predictor.risks == ["addiction", "cardiovascular", "anxiety"]
        assert safety_predictor.stats == {}

    def test_predict_toxicity(self, test_compounds, model_config):
        """Test toxicity prediction."""
        predictor = SafetyPredictor(config=model_config)
        
        # Make prediction
        result = predictor.predict_toxicity(test_compounds[0])
        
        # Check result
        assert isinstance(result, ModelResult)
        assert result.success
        assert "toxicity_score" in result.predictions
        assert "risk_factors" in result.predictions
        assert 0 <= result.predictions["toxicity_score"] <= 1

    def test_predict_abuse_potential(self, test_compounds, model_config):
        """Test abuse potential prediction."""
        predictor = SafetyPredictor(config=model_config)
        
        # Make prediction
        result = predictor.predict_abuse_potential(test_compounds[0])
        
        # Check result
        assert isinstance(result, ModelResult)
        assert result.success
        assert "abuse_score" in result.predictions
        assert "risk_factors" in result.predictions
        assert 0 <= result.predictions["abuse_score"] <= 1


class TestModelEvaluator:
    """Tests for ModelEvaluator class."""

    def test_initialization(self, model_evaluator):
        """Test initialization of ModelEvaluator."""
        assert model_evaluator.stats == {}

    def test_evaluate_model(self, model_evaluator, model_ensemble, test_compounds):
        """Test model evaluation."""
        # Evaluate model
        result = model_evaluator.evaluate_model(model_ensemble, test_compounds)
        
        # Check result
        assert result.success
        assert "accuracy" in result.metrics
        assert "precision" in result.metrics
        assert "recall" in result.metrics
        assert "f1_score" in result.metrics
        assert "evaluation_time" in result.stats

    def test_cross_validate(self, model_evaluator, model_ensemble, test_compounds):
        """Test cross validation."""
        # Cross validate
        result = model_evaluator.cross_validate(model_ensemble, test_compounds)
        
        # Check result
        assert result.success
        assert len(result.fold_scores) == 5  # Default 5-fold CV
        assert "mean_score" in result.metrics
        assert "std_score" in result.metrics
        assert "cross_validation_time" in result.stats


class TestModelEnsemble:
    """Tests for ModelEnsemble class."""

    def test_initialization(self, model_ensemble):
        """Test initialization of ModelEnsemble."""
        assert len(model_ensemble.models) == 3
        assert all(isinstance(m, MLModel) for m in model_ensemble.models)
        assert model_ensemble.stats == {}

    def test_ensemble_prediction(self, model_ensemble, test_compounds):
        """Test ensemble prediction."""
        # Make ensemble predictions
        predictions = model_ensemble.predict(test_compounds)
        
        # Check predictions
        assert isinstance(predictions, dict)
        assert len(predictions) == len(test_compounds)
        for compound_id, pred in predictions.items():
            assert isinstance(pred, tuple)  # (mean, std)
            assert 0 <= pred[0] <= 1  # mean
            assert pred[1] >= 0  # std

    def test_model_weighting(self, model_ensemble, test_compounds):
        """Test model weighting in ensemble."""
        # Set model weights
        weights = [0.5, 0.3, 0.2]
        model_ensemble.set_weights(weights)
        
        # Make weighted predictions
        predictions = model_ensemble.predict(test_compounds)
        
        # Check weighted predictions
        assert model_ensemble.weights == weights
        assert np.isclose(sum(weights), 1.0)
        assert "weighted_predictions" in model_ensemble.stats
        assert len(predictions) == len(test_compounds)


class TestUncertaintyEstimator:
    """Tests for UncertaintyEstimator class."""

    def test_initialization(self, uncertainty_estimator):
        """Test initialization of UncertaintyEstimator."""
        assert uncertainty_estimator.stats == {}

    def test_estimate_uncertainty(self, uncertainty_estimator, model_ensemble, test_compounds):
        """Test uncertainty estimation."""
        # Estimate uncertainty
        result = uncertainty_estimator.estimate_uncertainty(model_ensemble, test_compounds[0])
        
        # Check result
        assert result.success
        assert "aleatoric" in result.uncertainty
        assert "epistemic" in result.uncertainty
        assert 0 <= result.uncertainty["aleatoric"] <= 1
        assert 0 <= result.uncertainty["epistemic"] <= 1
        assert "estimation_time" in result.stats

    def test_calibrate_uncertainty(self, uncertainty_estimator, model_ensemble, test_compounds):
        """Test uncertainty calibration."""
        # Calibrate uncertainty
        result = uncertainty_estimator.calibrate_uncertainty(model_ensemble, test_compounds)
        
        # Check result
        assert result.success
        assert "calibration_score" in result.metrics
        assert 0 <= result.metrics["calibration_score"] <= 1
        assert "calibration_time" in result.stats


if __name__ == "__main__":
    pytest.main([__file__])
