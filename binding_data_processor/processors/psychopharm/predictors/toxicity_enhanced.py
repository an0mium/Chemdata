"""Enhanced toxicity predictor using ensemble models.

This module enhances the ToxicityPredictorBase with:
1. Additional training methods
2. Enhanced model calibration
3. Improved prediction confidence
4. Better model versioning
5. Enhanced feature extraction
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Any, TypeAlias

from ....models.core import CompoundData
from ....models.psychopharm import ToxicityClass
from .types import PredictionResult
from .base import ToxicityPredictorBase

ImportanceDict: TypeAlias = Dict[str, Dict[str, float]]


class ToxicityPredictorEnhanced(ToxicityPredictorBase):
    """Enhanced toxicity predictor with additional training capabilities."""

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        toxicity_mechanisms: Optional[Dict[str, Set[str]]] = None,
        organ_toxicity: Optional[Dict[str, Set[str]]] = None,
        safety_concerns: Optional[Dict[str, Set[str]]] = None,
        bbb_model_dir: Optional[str] = None,
    ):
        """Initialize enhanced toxicity predictor.

        Args:
            model_dir: Optional directory containing trained models
            cache_dir: Optional directory for caching
            log_level: Logging level
            toxicity_mechanisms: Optional custom toxicity mechanisms
            organ_toxicity: Optional custom organ toxicity categories
            safety_concerns: Optional custom safety concerns
            bbb_model_dir: Optional directory containing BBB models
        """
        super().__init__(
            model_dir=model_dir,
            cache_dir=cache_dir,
            log_level=log_level,
            toxicity_mechanisms=toxicity_mechanisms,
            organ_toxicity=organ_toxicity,
            safety_concerns=safety_concerns,
            bbb_model_dir=bbb_model_dir,
        )

        self.logger.info("ToxicityPredictorEnhanced initialized successfully")

    def retrain(
        self,
        compounds: List[CompoundData],
        mechanism_labels: List[ToxicityClass],
        mechanism_data: Dict[str, Dict[str, List[float]]],
        organ_data: Dict[str, Dict[str, List[float]]],
        safety_data: Dict[str, Dict[str, List[float]]],
        **kwargs,
    ) -> Dict[str, float]:
        """Retrain ensemble models with enhanced metrics.

        Args:
            compounds: List of compounds to train on
            mechanism_labels: Toxicity mechanism labels
            mechanism_data: Mechanism scores for training
            organ_data: Organ toxicity scores for training
            safety_data: Safety concern scores for training
            **kwargs: Additional training parameters

        Returns:
            Dictionary of training metrics
        """
        self.logger.info(f"Retraining ensemble models with {len(compounds)} compounds")

        try:
            metrics = {}

            # Train mechanism ensemble
            self.logger.debug("Training mechanism ensemble")
            mechanism_metrics = self.mechanism_ensemble.fit(
                "mechanism",
                compounds,
                mechanism_labels,
            )
            # Update metrics with mechanism results
            mechanism_dict = {f"mechanism_{k}": v for k, v in mechanism_metrics.items()}
            metrics.update(mechanism_dict)

            # Train mechanism, organ, and safety ensembles
            mechanism_metrics = self._train_mechanism_ensembles(compounds, mechanism_data)
            organ_metrics = self._train_organ_ensembles(compounds, organ_data)
            safety_metrics = self._train_safety_ensembles(compounds, safety_data)
            metrics.update(mechanism_metrics)
            metrics.update(organ_metrics)
            metrics.update(safety_metrics)

            # Save updated models
            self.save_models()

            self.logger.info("Ensemble model retraining completed successfully")
            return metrics

        except Exception as e:
            self.logger.error(f"Error retraining ensemble models: {str(e)}")
            return {"error": str(e)}

    def _train_safety_ensembles(
        self,
        compounds: List[CompoundData],
        safety_data: Dict[str, Dict[str, List[float]]],
    ) -> Dict[str, float]:
        """Train safety concern ensembles with enhanced metrics.

        Args:
            compounds: List of compounds to train on
            safety_data: Safety concern scores for training

        Returns:
            Dictionary of training metrics
        """
        metrics = {}
        for category, concerns in safety_data.items():
            for concern, risks in concerns.items():
                if concern in self.safety_concerns.get(category, []):
                    self.logger.debug(f"Training safety ensemble for {concern}")
                    concern_metrics = self.safety_ensembles[category][concern].fit(
                        concern,
                        compounds,
                        risks,
                    )
                    # Create metrics dict with prefixed keys
                    concern_metrics_dict = {}
                    for k, v in concern_metrics.items():
                        key = f"concern_{concern}_{k}"
                        concern_metrics_dict[key] = v
                    metrics.update(concern_metrics_dict)
        return metrics
