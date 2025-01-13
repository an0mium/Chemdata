"""Machine learning manager for pipeline.

This module provides the MLManager class that handles:
1. Model loading and management
2. Feature extraction
3. Prediction generation
4. Model retraining
5. Prediction caching

The manager supports multiple model types:
- Receptor binding prediction
- Psychoactive effects prediction
- Nootropic activity prediction
- Abuse potential prediction
- Toxicity prediction
- Blood-brain barrier penetration
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, field
import json
import numpy as np
import torch
from transformers import pipeline
from sentence_transformers import SentenceTransformer

from ..models import CompoundData
from ..processors.psychopharm.predictors import (
    ReceptorPredictor,
    PsychoactivePredictor,
    NootropicPredictor,
    AbusePredictor,
    ToxicityPredictor,
    BBBPredictor,
)
from ..processors.structure.ml.models.ensemble import EnsembleModel


@dataclass
class ModelConfig:
    """ML model configuration."""

    # Model paths
    receptor_model: str = "models/receptor"
    psychoactive_model: str = "models/psychoactive"
    nootropic_model: str = "models/nootropic"
    abuse_model: str = "models/abuse"
    toxicity_model: str = "models/toxicity"
    bbb_model: str = "models/bbb"

    # Pipeline models
    text_classifier: str = "models/pipeline/compound_classifier"
    activity_classifier: str = "models/pipeline/activity_classifier"
    safety_classifier: str = "models/pipeline/safety_classifier"
    mechanism_classifier: str = "models/pipeline/mechanism_classifier"

    # Transformer models
    similarity_model: str = "sentence-transformers/all-MiniLM-L6-v2"

    # Prediction settings
    confidence_threshold: float = 0.5
    batch_size: int = 32
    use_cache: bool = True
    device: str = "cuda" if torch.cuda.is_available() else "cpu"


@dataclass
class PredictionStats:
    """ML prediction statistics."""

    # Prediction counts
    total_predictions: int = 0
    successful_predictions: int = 0
    failed_predictions: int = 0
    cached_predictions: int = 0

    # Model stats
    model_stats: Dict[str, Dict[str, float]] = field(default_factory=dict)

    # Error tracking
    errors: List[Dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary format."""
        return {
            "predictions": {
                "total": self.total_predictions,
                "successful": self.successful_predictions,
                "failed": self.failed_predictions,
                "cached": self.cached_predictions,
                "success_rate": self._get_success_rate(),
            },
            "models": self.model_stats,
            "errors": self.errors,
        }

    def _get_success_rate(self) -> Optional[float]:
        """Get prediction success rate."""
        if not self.total_predictions:
            return None
        return self.successful_predictions / self.total_predictions


class MLManager:
    """Manager for ML models and predictions."""

    def __init__(
        self,
        model_dir: Union[str, Path],
        cache_dir: Optional[Union[str, Path]] = None,
        config: Optional[ModelConfig] = None,
    ):
        """Initialize ML manager.

        Args:
            model_dir: Directory containing ML models
            cache_dir: Optional directory for caching
            config: Optional model configuration
        """
        self.logger = logging.getLogger(self.__class__.__name__)
        self.model_dir = Path(model_dir)
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.config = config or ModelConfig()

        # Initialize models
        self._init_models()

        # Initialize cache
        self._prediction_cache: Dict[str, Dict] = {}

        # Initialize stats
        self.stats = PredictionStats()

    def _init_models(self) -> None:
        """Initialize ML models."""
        try:
            # Base predictors
            self.predictors = {
                "receptor": ReceptorPredictor(
                    model_dir=self.model_dir / self.config.receptor_model,
                    device=self.config.device,
                ),
                "psychoactive": PsychoactivePredictor(
                    model_dir=self.model_dir / self.config.psychoactive_model,
                    device=self.config.device,
                ),
                "nootropic": NootropicPredictor(
                    model_dir=self.model_dir / self.config.nootropic_model,
                    device=self.config.device,
                ),
                "abuse": AbusePredictor(
                    model_dir=self.model_dir / self.config.abuse_model,
                    device=self.config.device,
                ),
                "toxicity": ToxicityPredictor(
                    model_dir=self.model_dir / self.config.toxicity_model,
                    device=self.config.device,
                ),
                "bbb": BBBPredictor(
                    model_dir=self.model_dir / self.config.bbb_model,
                    device=self.config.device,
                ),
            }

            # Pipeline models
            self.pipelines = {
                "text": pipeline(
                    "text-classification",
                    model=self.model_dir / self.config.text_classifier,
                    device=self.config.device,
                    from_pt=True,
                ),
                "activity": pipeline(
                    "text-classification",
                    model=self.model_dir / self.config.activity_classifier,
                    device=self.config.device,
                    from_pt=True,
                ),
                "safety": pipeline(
                    "text-classification",
                    model=self.model_dir / self.config.safety_classifier,
                    device=self.config.device,
                    from_pt=True,
                ),
                "mechanism": pipeline(
                    "text-classification",
                    model=self.model_dir / self.config.mechanism_classifier,
                    device=self.config.device,
                    from_pt=True,
                ),
            }

            # Transformer models
            self.similarity_model = SentenceTransformer(self.config.similarity_model).to(self.config.device)

            # Ensemble models
            self.ensembles = {
                "activity": EnsembleModel(
                    [
                        self.predictors["psychoactive"],
                        self.predictors["nootropic"],
                    ],
                    weights=[0.5, 0.5],
                ),
                "safety": EnsembleModel(
                    [
                        self.predictors["toxicity"],
                        self.predictors["abuse"],
                    ],
                    weights=[0.5, 0.5],
                ),
            }

            self.logger.info("Successfully initialized ML models")

        except Exception as e:
            self.logger.error(f"Failed to initialize ML models: {str(e)}")
            raise

    def predict_compound(
        self,
        compound: CompoundData,
        use_cache: Optional[bool] = None,
    ) -> CompoundData:
        """Generate predictions for a compound.

        Args:
            compound: CompoundData instance to generate predictions for
            use_cache: Whether to use prediction cache

        Returns:
            CompoundData with predictions
        """
        try:
            self.stats.total_predictions += 1

            # Check cache
            if use_cache if use_cache is not None else self.config.use_cache:
                cached = self._get_cached_predictions(compound)
                if cached:
                    self.stats.cached_predictions += 1
                    return self._apply_cached_predictions(compound, cached)

            # Run base predictions
            for name, predictor in self.predictors.items():
                try:
                    result = predictor.predict(compound)
                    compound.set_prediction(name, result)
                    self._update_model_stats(name, result)
                except Exception as e:
                    self.logger.error(f"Failed to run {name} prediction: {str(e)}")
                    self.stats.errors.append(
                        {
                            "type": f"{name}_prediction_error",
                            "compound": compound.name,
                            "error": str(e),
                        }
                    )

            # Run ensemble predictions
            for name, ensemble in self.ensembles.items():
                try:
                    result = ensemble.predict(compound)
                    compound.set_prediction(f"{name}_ensemble", result)
                    self._update_model_stats(f"{name}_ensemble", result)
                except Exception as e:
                    self.logger.error(f"Failed to run {name} ensemble: {str(e)}")
                    self.stats.errors.append(
                        {
                            "type": f"{name}_ensemble_error",
                            "compound": compound.name,
                            "error": str(e),
                        }
                    )

            # Cache predictions
            if self.config.use_cache:
                self._cache_predictions(compound)

            self.stats.successful_predictions += 1
            return compound

        except Exception as e:
            self.logger.error(f"Failed to generate predictions for {compound.name}: {str(e)}")
            self.stats.failed_predictions += 1
            self.stats.errors.append(
                {
                    "type": "prediction_error",
                    "compound": compound.name,
                    "error": str(e),
                }
            )
            raise

    def _get_cached_predictions(
        self,
        compound: CompoundData,
    ) -> Optional[Dict]:
        """Get cached predictions for compound."""
        return self._prediction_cache.get(compound.smiles)

    def _apply_cached_predictions(
        self,
        compound: CompoundData,
        cached: Dict,
    ) -> CompoundData:
        """Apply cached predictions to compound."""
        for name, prediction in cached.items():
            compound.set_prediction(name, prediction)
        return compound

    def _cache_predictions(
        self,
        compound: CompoundData,
    ) -> None:
        """Cache predictions for compound."""
        self._prediction_cache[compound.smiles] = compound.get_predictions()

    def _update_model_stats(
        self,
        model: str,
        result: Dict,
    ) -> None:
        """Update model statistics."""
        if model not in self.stats.model_stats:
            self.stats.model_stats[model] = {
                "total": 0,
                "confidence_sum": 0,
                "high_confidence": 0,
            }

        stats = self.stats.model_stats[model]
        stats["total"] += 1
        stats["confidence_sum"] += result.get("confidence", 0)
        if result.get("confidence", 0) > self.config.confidence_threshold:
            stats["high_confidence"] += 1

    def get_model_versions(self) -> Dict[str, str]:
        """Get model version information."""
        versions = {}

        # Base predictors
        for name, predictor in self.predictors.items():
            versions[name] = predictor.version

        # Pipeline models
        for name, pipeline in self.pipelines.items():
            versions[f"{name}_pipeline"] = pipeline.model.config.version

        # Transformer models
        versions["similarity"] = self.similarity_model.get_version()

        return versions

    def clear_cache(self) -> None:
        """Clear prediction cache."""
        self._prediction_cache.clear()
        self.stats.cached_predictions = 0
