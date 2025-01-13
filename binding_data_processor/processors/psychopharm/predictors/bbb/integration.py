"""BBB permeability predictor with integrated predictions.

This module provides:
1. BBBPredictorBase - Core BBB permeability prediction functionality
2. BBBIntegrationMixin - Mixin for integrating BBB prediction into other predictors
3. BBBPredictorEnhanced - Enhanced BBB predictor with integrated predictions from:
   - Base BBB permeability prediction
   - Transporter analysis
   - Abuse potential prediction 
   - Toxicity prediction
   - Receptor binding prediction
   - Psychoactive effects prediction
   - Nootropic activity prediction
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Any, Tuple
import pandas as pd
import numpy as np
from sklearn.calibration import CalibratedClassifierCV
from sklearn.isotonic import IsotonicRegression

from .....models.core import CompoundData
from .....models.psychopharm import NootropicMechanism, BBBPermeability
from ..types import PredictionResult
from ..abuse import AbusePotentialPredictor
from ..receptors import ReceptorPredictor
from ..psychoactive import PsychoactivePredictor
from ..nootropic.base import NootropicPredictor
from .predictors import BBBPredictor
from .....models.compound.ml.ensemble import ClassifierEnsemble, RegressorEnsemble
from ..base import EnsemblePredictor


class BBBPredictorBase(EnsemblePredictor):
    """Base class for BBB permeability prediction."""

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        web_enrichment: bool = True,
    ):
        """Initialize BBB predictor.

        Args:
            model_dir: Optional directory containing trained models
            cache_dir: Optional directory for caching
            log_level: Logging level
            web_enrichment: Whether to use web enrichment
        """
        super().__init__(
            model_dir=model_dir,
            cache_dir=cache_dir,
            feature_types=["fingerprints", "descriptors", "enhanced"],
            log_level=log_level,
        )

        self.web_enrichment = web_enrichment
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)

        # Add file handler if model_dir provided
        if model_dir:
            log_path = Path(model_dir) / "bbb_predictor.log"
            fh = logging.FileHandler(log_path)
            fh.setLevel(log_level)
            formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            fh.setFormatter(formatter)
            self.logger.addHandler(fh)

        # Initialize ensemble models
        self.ensemble = ClassifierEnsemble(
            model_dir=self.model_dir,
            cache_dir=self.cache_dir,
            feature_types=self.feature_types,
            log_level=log_level,
        )

    def predict(self, compound: CompoundData) -> PredictionResult:
        """Predict BBB permeability for a compound.

        Args:
            compound: Compound to predict BBB permeability for

        Returns:
            PredictionResult containing BBB prediction
        """
        # Get base prediction from ensemble
        prediction, confidence = self.ensemble.predict("bbb", compound)

        # Get predictions from individual models for uncertainty
        predictions = []
        probabilities = []
        for model in self.ensemble.models["bbb"].estimators_:
            pred = model.predict([compound])[0]
            prob = model.predict_proba([compound])[0]
            predictions.append(pred)
            probabilities.append(prob)

        # Calculate uncertainty metrics
        predictions = np.array(predictions)
        probabilities = np.array(probabilities)

        # Model disagreement
        disagreement = 1 - (np.sum(predictions == prediction) / len(predictions))

        # Probability variance
        prob_variance = np.var(probabilities, axis=0)

        # Combine metrics for final confidence
        adjusted_confidence = confidence * (1 - disagreement) * (1 - np.mean(prob_variance))

        return PredictionResult(
            value=prediction,
            confidence=adjusted_confidence,
            supporting_data={
                "model_disagreement": disagreement,
                "probability_variance": float(np.mean(prob_variance)),
            },
        )

    def _load_models(self) -> Dict:
        """Load ensemble models."""
        self.logger.info("Loading BBB ensemble models")

        try:
            self.ensemble.load_models()
            self.logger.info("BBB ensemble models loaded successfully")

        except Exception as e:
            self.logger.error(f"Error loading BBB ensemble models: {str(e)}")
            self.logger.info("Initializing new BBB ensemble models")
            self._initialize_models()

        return {}  # Models stored in ensemble

    def _initialize_models(self) -> Dict:
        """Initialize ensemble models."""
        self.logger.info("Initializing BBB ensemble models")

        self.ensemble.create_ensemble(
            "bbb",
            ensemble_type="voting",
        )

        return {}  # Models stored in ensemble

    def save_models(self, save_dir: Optional[str] = None) -> None:
        """Save ensemble models."""
        save_dir = Path(save_dir) if save_dir else self.model_dir
        if not save_dir:
            return

        save_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Saving BBB models to {save_dir}")

        try:
            self.ensemble.save_models(save_dir)
            self.logger.info("BBB ensemble models saved successfully")

        except Exception as e:
            self.logger.error(f"Error saving BBB ensemble models: {str(e)}")


class BBBIntegrationMixin:
    """Mixin for integrating BBB prediction into other predictors."""

    def __init__(
        self,
        bbb_model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize BBB integration.

        Args:
            bbb_model_dir: Optional directory containing BBB models
            cache_dir: Optional directory for caching
            log_level: Logging level
        """
        self.bbb_predictor = BBBPredictorBase(
            model_dir=bbb_model_dir,
            cache_dir=cache_dir,
            log_level=log_level,
        )

    def predict_bbb(self, compound: CompoundData) -> PredictionResult:
        """Predict BBB permeability for a compound.

        Args:
            compound: Compound to predict BBB permeability for

        Returns:
            PredictionResult containing BBB prediction
        """
        return self.bbb_predictor.predict(compound)


class BBBPredictorEnhanced(BBBPredictor):
    """BBB permeability predictor with integrated predictions."""

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        transporters: Optional[Dict[str, Set[str]]] = None,
        receptor_transporters: Optional[Dict[str, Set[str]]] = None,
        predictors: Optional[Dict[str, Any]] = None,
    ):
        """Initialize enhanced BBB predictor.

        Args:
            model_dir: Optional directory containing trained models
            cache_dir: Optional directory for caching
            log_level: Logging level
            transporters: Optional custom transporter definitions
            receptor_transporters: Optional custom receptor transporter definitions
            predictors: Optional custom predictor instances
        """
        super().__init__(
            model_dir=model_dir,
            cache_dir=cache_dir,
            log_level=log_level,
            transporters=transporters,
            receptor_transporters=receptor_transporters,
        )

        # Initialize other predictors
        self.predictors = predictors or {
            "abuse": AbusePotentialPredictor(
                model_dir=str(Path(model_dir) / "abuse") if model_dir else None,
                cache_dir=cache_dir,
                log_level=log_level,
            ),
            "receptors": ReceptorPredictor(
                model_dir=str(Path(model_dir) / "receptors") if model_dir else None,
                cache_dir=cache_dir,
                log_level=log_level,
            ),
            "psychoactive": PsychoactivePredictor(
                model_dir=str(Path(model_dir) / "psychoactive") if model_dir else None,
                cache_dir=cache_dir,
                log_level=log_level,
            ),
            "nootropic": NootropicPredictor(
                model_dir=str(Path(model_dir) / "nootropic") if model_dir else None,
                cache_dir=cache_dir,
                log_level=log_level,
            ),
        }

        # Add integrated prediction columns to history
        self.prediction_history = pd.concat(
            [
                self.prediction_history,
                pd.DataFrame(
                    columns=[
                        "abuse_potential",
                        "abuse_confidence",
                        "receptor_binding",
                        "receptor_confidence",
                        "psychoactive_class",
                        "psychoactive_confidence",
                        "nootropic_mechanisms",
                        "nootropic_confidence",
                    ]
                ),
            ],
            axis=1,
        )

        self.logger.info("BBBPredictorEnhanced initialized successfully")

    def predict(self, compound: CompoundData) -> PredictionResult:
        """Generate comprehensive BBB permeability predictions.

        Args:
            compound: Compound to predict BBB permeability for

        Returns:
            PredictionResult containing:
            - BBB permeability class
            - Confidence score
            - Supporting data including:
                - Permeability score
                - Transporter interactions
                - Receptor-mediated transport
                - Duration metrics
                - Abuse potential
                - Receptor binding
                - Psychoactive effects
                - Nootropic activity
        """
        self.logger.debug(f"Generating predictions for {compound.name}")
        try:
            # Get base predictions
            result = super().predict(compound)

            # Get integrated predictions
            integrated_predictions = self._get_integrated_predictions(compound)

            # Update supporting data
            result.supporting_data.update(integrated_predictions)

            # Update prediction history
            self._update_integrated_prediction_history(integrated_predictions, compound)

            return result

        except Exception as e:
            self.logger.error("Error predicting BBB properties", f"Compound {compound.name}: {str(e)}", exc_info=True)
            return PredictionResult(
                value=BBBPermeability.UNKNOWN,
                confidence=0.0,
                supporting_data={"error": str(e)},
            )

    def _get_integrated_predictions(self, compound: CompoundData) -> Dict[str, Dict[str, Any]]:
        """Get predictions from integrated predictors."""
        predictions = {}

        # Get abuse potential predictions
        abuse_result = self.predictors["abuse"].predict(compound)
        predictions["abuse"] = {
            "potential": abuse_result.value,
            "confidence": abuse_result.confidence,
            "supporting_data": abuse_result.supporting_data,
        }

        # Get receptor predictions
        receptor_result = self.predictors["receptors"].predict(compound)
        predictions["receptors"] = {
            "binding": receptor_result.value,
            "confidence": receptor_result.confidence,
            "supporting_data": receptor_result.supporting_data,
        }

        # Get psychoactive predictions
        psychoactive_result = self.predictors["psychoactive"].predict(compound)
        predictions["psychoactive"] = {
            "class": psychoactive_result.value,
            "confidence": psychoactive_result.confidence,
            "supporting_data": psychoactive_result.supporting_data,
        }

        # Get nootropic predictions
        nootropic_result = self.predictors["nootropic"].predict(compound)
        predictions["nootropic"] = {
            "mechanisms": nootropic_result.value,
            "confidence": nootropic_result.confidence,
            "supporting_data": nootropic_result.supporting_data,
        }

        return predictions

    def _update_integrated_prediction_history(
        self,
        predictions: Dict[str, Dict[str, Any]],
        compound: CompoundData,
    ) -> None:
        """Update prediction history with integrated predictions."""
        # Get latest prediction index
        latest_idx = self.prediction_history.index[-1]

        # Update abuse predictions
        self.prediction_history.at[latest_idx, "abuse_potential"] = predictions["abuse"]["potential"].value
        self.prediction_history.at[latest_idx, "abuse_confidence"] = predictions["abuse"]["confidence"]

        # Update receptor predictions
        self.prediction_history.at[latest_idx, "receptor_binding"] = str(predictions["receptors"]["binding"])
        self.prediction_history.at[latest_idx, "receptor_confidence"] = predictions["receptors"]["confidence"]

        # Update psychoactive predictions
        self.prediction_history.at[latest_idx, "psychoactive_class"] = predictions["psychoactive"]["class"].value
        self.prediction_history.at[latest_idx, "psychoactive_confidence"] = predictions["psychoactive"]["confidence"]

        # Update nootropic predictions
        self.prediction_history.at[latest_idx, "nootropic_mechanisms"] = ",".join(m.value for m in predictions["nootropic"]["mechanisms"])
        self.prediction_history.at[latest_idx, "nootropic_confidence"] = predictions["nootropic"]["confidence"]

    def retrain(
        self,
        compounds: List[CompoundData],
        labels: List[BBBPermeability],
        scores: Optional[List[float]] = None,
        transporter_data: Optional[Dict[str, List[bool]]] = None,
        receptor_data: Optional[Dict[str, List[bool]]] = None,
        duration_labels: Optional[List[str]] = None,
        abuse_labels: Optional[List[Any]] = None,
        receptor_binding_labels: Optional[List[Any]] = None,
        psychoactive_labels: Optional[List[Any]] = None,
        nootropic_labels: Optional[List[Any]] = None,
        **kwargs,
    ) -> Dict[str, float]:
        """Retrain models with new data.

        Args:
            compounds: List of compounds to train on
            labels: BBB permeability class labels
            scores: Optional permeability scores for regression
            transporter_data: Optional transporter substrate labels
            receptor_data: Optional receptor transport labels
            duration_labels: Optional duration class labels
            abuse_labels: Optional abuse potential labels
            receptor_binding_labels: Optional receptor binding labels
            psychoactive_labels: Optional psychoactive class labels
            nootropic_labels: Optional nootropic mechanism labels
            **kwargs: Additional training parameters

        Returns:
            Dictionary of training metrics
        """
        # Train base models
        metrics = super().retrain(
            compounds=compounds,
            labels=labels,
            scores=scores,
            transporter_data=transporter_data,
            receptor_data=receptor_data,
            duration_labels=duration_labels,
            **kwargs,
        )

        # Train abuse potential model
        if abuse_labels is not None:
            abuse_metrics = self.predictors["abuse"].retrain(
                compounds=compounds,
                labels=abuse_labels,
            )
            metrics.update({f"abuse_{k}": v for k, v in abuse_metrics.items()})

        # Train receptor binding model
        if receptor_binding_labels is not None:
            receptor_metrics = self.predictors["receptors"].retrain(
                compounds=compounds,
                labels=receptor_binding_labels,
            )
            metrics.update({f"receptor_{k}": v for k, v in receptor_metrics.items()})

        # Train psychoactive model
        if psychoactive_labels is not None:
            psychoactive_metrics = self.predictors["psychoactive"].retrain(
                compounds=compounds,
                labels=psychoactive_labels,
            )
            metrics.update({f"psychoactive_{k}": v for k, v in psychoactive_metrics.items()})

        # Train nootropic model
        if nootropic_labels is not None:
            nootropic_metrics = self.predictors["nootropic"].retrain(
                compounds=compounds,
                labels=nootropic_labels,
            )
            metrics.update({f"nootropic_{k}": v for k, v in nootropic_metrics.items()})

        self.logger.info("Model retraining completed successfully")
        return metrics
