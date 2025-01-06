"""Enhanced toxicity predictor using ensemble models.

This module enhances the ToxicityPredictor with:
1. Ensemble model integration
2. BBB permeability integration
3. Improved prediction confidence
4. Better model versioning
5. Enhanced feature extraction
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Any, TypeAlias

import numpy as np
from sklearn.calibration import CalibratedClassifierCV
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import brier_score_loss

from ....models.core import CompoundData
from ....models.psychopharm import ToxicityClass
from ....models.compound.ml.ensemble import ClassifierEnsemble, RegressorEnsemble
from ..base import PredictionResult
from .toxicity import ToxicityPredictor
from ..bbb.integration import BBBPredictorEnhanced

ImportanceDict: TypeAlias = Dict[str, Dict[str, float]]


class ToxicityPredictorEnhanced(ToxicityPredictor):
    """Enhanced toxicity predictor using ensemble models and BBB integration."""

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
        )

        # Initialize BBB predictor with web enrichment
        self.bbb_predictor = BBBPredictorEnhanced(
            model_dir=bbb_model_dir,
            cache_dir=cache_dir,
            log_level=log_level,
            web_enrichment=True,
        )

        # Initialize ensemble models
        self.mechanism_ensemble = ClassifierEnsemble(
            model_dir=self.model_dir,
            cache_dir=self.cache_dir,
            feature_types=self.feature_types,
            log_level=log_level,
        )

        # Initialize mechanism ensembles
        self.mechanism_ensembles = {}
        for category, mechanisms in self.toxicity_mechanisms.items():
            category_ensembles = {}
            for mechanism in mechanisms:
                category_ensembles[mechanism] = RegressorEnsemble(
                    model_dir=self.model_dir,
                    cache_dir=self.cache_dir,
                    feature_types=self.feature_types,
                    log_level=log_level,
                )
            self.mechanism_ensembles[category] = category_ensembles

        # Initialize organ toxicity ensembles
        self.organ_ensembles = {}
        for organ, toxicities in self.organ_toxicity.items():
            organ_ensembles = {}
            for toxicity in toxicities:
                organ_ensembles[toxicity] = RegressorEnsemble(
                    model_dir=self.model_dir,
                    cache_dir=self.cache_dir,
                    feature_types=self.feature_types,
                    log_level=log_level,
                )
            self.organ_ensembles[organ] = organ_ensembles

        # Initialize safety concern ensembles
        self.safety_ensembles = {}
        for category, concerns in self.safety_concerns.items():
            category_ensembles = {}
            for concern in concerns:
                category_ensembles[concern] = RegressorEnsemble(
                    model_dir=self.model_dir,
                    cache_dir=self.cache_dir,
                    feature_types=self.feature_types,
                    log_level=log_level,
                )
            self.safety_ensembles[category] = category_ensembles

        self.logger.info("ToxicityPredictorEnhanced initialized successfully")

    def _load_models(self) -> Dict:
        """Load ensemble models."""
        self.logger.info("Loading ensemble models")

        try:
            # Load BBB predictor models
            self.bbb_predictor._load_models()

            # Load mechanism ensemble
            self.mechanism_ensemble.load_models()

            # Load mechanism ensembles
            for category, ensembles in self.mechanism_ensembles.items():
                for mechanism, ensemble in ensembles.items():
                    ensemble.load_models()

            # Load organ toxicity ensembles
            for organ, ensembles in self.organ_ensembles.items():
                for toxicity, ensemble in ensembles.items():
                    ensemble.load_models()

            # Load safety concern ensembles
            for category, ensembles in self.safety_ensembles.items():
                for concern, ensemble in ensembles.items():
                    ensemble.load_models()

            self.logger.info("Ensemble models loaded successfully")

        except Exception as e:
            self.logger.error(f"Error loading ensemble models: {str(e)}")
            self.logger.info("Initializing new ensemble models")
            self._initialize_models()

    def _initialize_models(self) -> Dict:
        """Initialize ensemble models."""
        self.logger.info("Initializing ensemble models")

        # Initialize mechanism ensemble
        self.mechanism_ensemble.create_ensemble(
            "mechanism",
            ensemble_type="voting",
        )

        # Initialize mechanism ensembles
        for category, mechanisms in self.toxicity_mechanisms.items():
            for mechanism in mechanisms:
                self.mechanism_ensembles[category][mechanism].create_ensemble(
                    mechanism,
                    ensemble_type="stacking",
                )

        # Initialize organ toxicity ensembles
        for organ, toxicities in self.organ_toxicity.items():
            for toxicity in toxicities:
                self.organ_ensembles[organ][toxicity].create_ensemble(
                    toxicity,
                    ensemble_type="stacking",
                )

        # Initialize safety concern ensembles
        for category, concerns in self.safety_concerns.items():
            for concern in concerns:
                self.safety_ensembles[category][concern].create_ensemble(
                    concern,
                    ensemble_type="stacking",
                )

        return {}  # Models stored in ensembles

    def save_models(self, save_dir: Optional[str] = None) -> None:
        """Save ensemble models."""
        save_dir = Path(save_dir) if save_dir else self.model_dir
        if not save_dir:
            return

        save_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Saving models to {save_dir}")

        try:
            # Save BBB predictor models
            self.bbb_predictor.save_models(save_dir / "bbb")

            # Save mechanism ensemble
            self.mechanism_ensemble.save_models(save_dir)

            # Save mechanism ensembles
            for category, ensembles in self.mechanism_ensembles.items():
                category_dir = save_dir / category
                for mechanism, ensemble in ensembles.items():
                    ensemble.save_models(category_dir)

            # Save organ toxicity ensembles
            for organ, ensembles in self.organ_ensembles.items():
                organ_dir = save_dir / organ
                for toxicity, ensemble in ensembles.items():
                    ensemble.save_models(organ_dir)

            # Save safety concern ensembles
            for category, ensembles in self.safety_ensembles.items():
                category_dir = save_dir / category
                for concern, ensemble in ensembles.items():
                    ensemble.save_models(category_dir)

            # Save prediction history
            history_path = save_dir / "prediction_history.csv"
            self.prediction_history.to_csv(history_path, index=False)

            self.logger.info("Ensemble models saved successfully")

        except Exception as e:
            self.logger.error(f"Error saving ensemble models: {str(e)}")

    def _predict_mechanism(self, compound: CompoundData) -> Tuple[str, float]:
        """Predict mechanism using ensemble with uncertainty estimation.

        Uses:
        1. Ensemble disagreement for uncertainty
        2. Calibrated probabilities for confidence
        3. Feature importance analysis
        """
        # Get raw predictions and probabilities from ensemble
        mechanism, raw_confidence = self.mechanism_ensemble.predict("mechanism", compound)

        # Get predictions from individual models for uncertainty
        predictions = []
        probabilities = []
        importances = {}

        for name, model in self.mechanism_ensemble.models["mechanism"].estimators_:
            pred = model.predict([compound])[0]
            prob = model.predict_proba([compound])[0]
            predictions.append(pred)
            probabilities.append(prob)

            # Get feature importances if available
            if hasattr(model, "feature_importances_"):
                importances[name] = model.feature_importances_

        # Calculate uncertainty metrics
        predictions = np.array(predictions)
        probabilities = np.array(probabilities)

        # Model disagreement
        disagreement = 1 - (np.sum(predictions == mechanism) / len(predictions))

        # Probability variance
        prob_variance = np.var(probabilities, axis=0)

        # Combine metrics for final confidence
        confidence = raw_confidence * (1 - disagreement) * (1 - np.mean(prob_variance))

        # Store feature importance info
        if importances:
            self.mechanism_importances = importances

        return mechanism, confidence

    def calibrate_models(
        self,
        compounds: List[CompoundData],
        mechanism_labels: List[str],
        mechanism_data: Dict[str, Dict[str, List[float]]],
        organ_data: Dict[str, Dict[str, List[float]]],
        safety_data: Dict[str, Dict[str, List[float]]],
        cv: int = 5,
    ) -> Dict[str, float]:
        """Calibrate model probabilities using held-out data.

        Args:
            compounds: Compounds for calibration
            mechanism_labels: True mechanism labels
            mechanism_data: Mechanism scores for calibration
            organ_data: Organ toxicity scores for calibration
            safety_data: Safety concern scores for calibration
            cv: Number of cross-validation folds

        Returns:
            Dictionary of calibration metrics
        """
        metrics = {}

        # Calibrate mechanism ensemble
        base_model = self.mechanism_ensemble.models["mechanism"]
        calibrated = CalibratedClassifierCV(base_model, cv=cv, method="sigmoid")
        calibrated.fit(compounds, mechanism_labels)
        self.mechanism_ensemble.models["mechanism"] = calibrated

        # Calculate calibration metrics
        mechanism_probs = calibrated.predict_proba(compounds)
        metrics["mechanism_brier_score"] = brier_score_loss(mechanism_labels, mechanism_probs[:, 1])

        # Calibrate mechanism ensembles
        for category, mechanisms in self.mechanism_ensembles.items():
            for mechanism, ensemble in mechanisms.items():
                if mechanism in mechanism_data.get(category, {}):
                    true_scores = mechanism_data[category][mechanism]
                    calibrated = self._calibrate_regressor(ensemble, compounds, true_scores, cv=cv)
                    self.mechanism_ensembles[category][mechanism] = calibrated

        # Calibrate organ toxicity ensembles
        for organ, toxicities in self.organ_ensembles.items():
            for toxicity, ensemble in toxicities.items():
                if toxicity in organ_data.get(organ, {}):
                    true_scores = organ_data[organ][toxicity]
                    calibrated = self._calibrate_regressor(ensemble, compounds, true_scores, cv=cv)
                    self.organ_ensembles[organ][toxicity] = calibrated

        # Calibrate safety concern ensembles
        for category, concerns in self.safety_ensembles.items():
            for concern, ensemble in concerns.items():
                if concern in safety_data.get(category, {}):
                    true_scores = safety_data[category][concern]
                    calibrated = self._calibrate_regressor(ensemble, compounds, true_scores, cv=cv)
                    self.safety_ensembles[category][concern] = calibrated

        return metrics

    def _calibrate_regressor(
        self,
        ensemble: Any,
        compounds: List[CompoundData],
        true_scores: List[float],
        cv: int = 5,
    ) -> Any:
        """Calibrate regression ensemble predictions."""
        # Get uncalibrated predictions
        pred_scores = []
        for compound in compounds:
            score, _ = ensemble.predict(compound)
            pred_scores.append(score)

        # Fit isotonic regression calibration
        calibrator = IsotonicRegression(out_of_bounds="clip")
        calibrator.fit(pred_scores, true_scores)

        # Create calibrated ensemble
        calibrated_ensemble = ensemble
        calibrated_ensemble.calibrator = calibrator
        return calibrated_ensemble

    def _predict_mechanisms(
        self,
        compound: CompoundData,
        bbb_result: PredictionResult,
    ) -> Dict[str, Dict[str, Dict[str, float]]]:
        """Predict toxicity mechanisms using ensembles and BBB data."""
        mechanism_predictions = {}
        for category, ensembles in self.mechanism_ensembles.items():
            category_predictions = {}
            for mechanism, ensemble in ensembles.items():
                # Get base prediction
                score, conf = ensemble.predict(mechanism, compound)

                # Adjust based on BBB permeability for CNS mechanisms
                if conf > 0.5:  # Base confidence threshold
                    if category in ["neurological", "cognitive"]:
                        # Scale effect by BBB permeability
                        adjusted_score = score * bbb_result.confidence
                        # Combine confidences
                        adjusted_conf = conf * bbb_result.confidence
                    else:
                        # Non-CNS mechanisms less dependent on BBB
                        adjusted_score = score
                        adjusted_conf = conf

                    if adjusted_conf > 0.3:  # Adjusted confidence threshold
                        category_predictions[mechanism] = {
                            "score": adjusted_score,
                            "confidence": adjusted_conf,
                            "bbb_factor": bbb_result.confidence,
                        }
            if category_predictions:
                mechanism_predictions[category] = category_predictions
        return mechanism_predictions

    def _predict_organ_effects(
        self,
        compound: CompoundData,
        bbb_result: PredictionResult,
    ) -> Dict[str, Dict[str, Dict[str, float]]]:
        """Predict organ-specific toxicity using ensembles and BBB data."""
        organ_predictions = {}
        for organ, ensembles in self.organ_ensembles.items():
            organ_predictions_inner = {}
            for toxicity, ensemble in ensembles.items():
                # Get base prediction
                severity, conf = ensemble.predict(toxicity, compound)

                # Adjust based on BBB permeability for CNS effects
                if conf > 0.5:  # Base confidence threshold
                    if organ == "brain":
                        # Scale severity by BBB permeability
                        adjusted_severity = severity * bbb_result.confidence
                        # Combine confidences
                        adjusted_conf = conf * bbb_result.confidence
                    else:
                        # Non-CNS effects less dependent on BBB
                        adjusted_severity = severity
                        adjusted_conf = conf

                    if adjusted_conf > 0.3:  # Adjusted confidence threshold
                        organ_predictions_inner[toxicity] = {
                            "severity": adjusted_severity,
                            "confidence": adjusted_conf,
                            "bbb_factor": bbb_result.confidence,
                        }
            if organ_predictions_inner:
                organ_predictions[organ] = organ_predictions_inner
        return organ_predictions

    def _predict_safety_concerns(
        self,
        compound: CompoundData,
        bbb_result: PredictionResult,
    ) -> Dict[str, Dict[str, Dict[str, float]]]:
        """Predict safety concerns using ensembles and BBB data."""
        concern_predictions = {}
        for category, ensembles in self.safety_ensembles.items():
            category_predictions = {}
            for concern, ensemble in ensembles.items():
                # Get base prediction
                risk, conf = ensemble.predict(concern, compound)

                # Adjust based on BBB permeability for CNS-related concerns
                if conf > 0.5:  # Base confidence threshold
                    if category in ["neurological", "cognitive", "behavioral"]:
                        # Scale risk by BBB permeability
                        adjusted_risk = risk * bbb_result.confidence
                        # Combine confidences
                        adjusted_conf = conf * bbb_result.confidence
                    else:
                        # Non-CNS concerns less dependent on BBB
                        adjusted_risk = risk
                        adjusted_conf = conf

                    if adjusted_conf > 0.3:  # Adjusted confidence threshold
                        category_predictions[concern] = {
                            "risk": adjusted_risk,
                            "confidence": adjusted_conf,
                            "bbb_factor": bbb_result.confidence,
                        }
            if category_predictions:
                concern_predictions[category] = category_predictions
        return concern_predictions

    def predict(self, compound: CompoundData) -> PredictionResult:
        """Generate toxicity predictions using ensembles and BBB data."""
        self.logger.debug(f"Generating ensemble predictions for {compound.name}")

        try:
            # Get BBB prediction first
            bbb_result = self.bbb_predictor.predict(compound)
            value = bbb_result.value.value
            conf = bbb_result.confidence
            bbb_msg = f"BBB prediction: {value} " f"(confidence: {conf:.3f})"
            self.logger.debug(bbb_msg)

            # Predict mechanism
            mechanism, mechanism_confidence = self._predict_mechanism(compound)

            # Predict toxicity mechanisms with BBB integration
            mechanism_predictions = self._predict_mechanisms(compound, bbb_result)

            # Update prediction history for mechanisms
            for category, mechanisms in mechanism_predictions.items():
                for mechanism, pred in mechanisms.items():
                    self._update_prediction_history(
                        compound=compound,
                        mechanism=mechanism,
                        mechanism_confidence=mechanism_confidence,
                        mechanism_category=category,
                        mechanism_score=pred["score"],
                        mechanism_confidence=pred["confidence"],
                        bbb_permeability=bbb_result.value.value,
                        bbb_confidence=bbb_result.confidence,
                    )

            # Predict organ toxicity with BBB integration
            organ_predictions = self._predict_organ_effects(compound, bbb_result)

            # Update prediction history for organ toxicity
            for organ, toxicities in organ_predictions.items():
                for toxicity, pred in toxicities.items():
                    self._update_prediction_history(
                        compound=compound,
                        mechanism=mechanism,
                        mechanism_confidence=mechanism_confidence,
                        organ=organ,
                        toxicity_type=toxicity,
                        severity=pred["severity"],
                        organ_confidence=pred["confidence"],
                        bbb_permeability=bbb_result.value.value,
                        bbb_confidence=bbb_result.confidence,
                    )

            # Predict safety concerns with BBB integration
            concern_predictions = self._predict_safety_concerns(compound, bbb_result)

            # Update prediction history for safety concerns
            for category, concerns in concern_predictions.items():
                for concern, pred in concerns.items():
                    self._update_prediction_history(
                        compound=compound,
                        mechanism=mechanism,
                        mechanism_confidence=mechanism_confidence,
                        concern_category=category,
                        concern=concern,
                        risk_level=pred["risk"],
                        risk_confidence=pred["confidence"],
                        bbb_permeability=bbb_result.value.value,
                        bbb_confidence=bbb_result.confidence,
                    )

            # Format result
            result = PredictionResult(
                value=ToxicityClass(mechanism),
                confidence=mechanism_confidence,
                supporting_data={
                    "bbb_prediction": {
                        "value": bbb_result.value.value,
                        "confidence": bbb_result.confidence,
                        **bbb_result.supporting_data,
                    },
                    "mechanisms": mechanism_predictions,
                    "organs": organ_predictions,
                    "concerns": concern_predictions,
                    **self._format_supporting_data(
                        {},  # Features handled by ensembles
                        [
                            (mechanism, pred["confidence"])
                            for preds in mechanism_predictions.values()
                            for mechanism, pred in preds.items()
                        ],
                    ),
                },
            )

            self.logger.info(
                f"Generated ensemble predictions for {compound.name}: "
                f"{result.value.value} (confidence: {result.confidence:.3f})"
            )

            return result

        except Exception as e:
            self.logger.error(
                "Error generating ensemble predictions",
                f"Compound {compound.name}: {str(e)}",
                exc_info=True,
            )
            return PredictionResult(
                value=ToxicityClass.UNKNOWN,
                confidence=0.0,
                supporting_data={"error": str(e)},
            )

    def _get_mechanism_importances(self) -> ImportanceDict:
        """Get feature importances from mechanism models."""
        importances = {}
        if hasattr(self, "mechanism_importances"):
            importances["mechanism"] = self.mechanism_importances
        return importances

    def _get_ensemble_importances(
        self,
        ensembles: Dict[str, Dict[str, Any]],
        prefix: str,
    ) -> ImportanceDict:
        """Get feature importances from ensemble models.

        Args:
            ensembles: Dictionary of ensemble models
            prefix: Prefix for importance keys

        Returns:
            Dictionary of feature importance scores
        """
        importances = {}
        for category, category_ensembles in ensembles.items():
            for name, ensemble in category_ensembles.items():
                if hasattr(ensemble, "feature_importances_"):
                    key = f"{prefix}_{category}_{name}"
                    importances[key] = {"importance": ensemble.feature_importances_}
        return importances

    def get_feature_importance(self, mechanism: Optional[str] = None) -> ImportanceDict:
        """Get feature importance scores for models.

        Args:
            mechanism: Optional specific mechanism to get importances for

        Returns:
            Dictionary of feature importance scores by model/mechanism
        """
        importances = {}

        # Get importances from different model types
        importances.update(self._get_mechanism_importances())
        importances.update(self._get_ensemble_importances(self.mechanism_ensembles, "mechanism"))
        importances.update(self._get_ensemble_importances(self.organ_ensembles, "organ"))
        importances.update(self._get_ensemble_importances(self.safety_ensembles, "concern"))

        return importances

    def _train_mechanism_ensembles(
        self,
        compounds: List[CompoundData],
        mechanism_data: Dict[str, Dict[str, List[float]]],
    ) -> Dict[str, float]:
        """Train toxicity mechanism ensembles."""
        metrics = {}
        for category, mechanisms in mechanism_data.items():
            for mechanism, scores in mechanisms.items():
                if mechanism in self.toxicity_mechanisms.get(category, []):
                    self.logger.debug(f"Training mechanism ensemble for {mechanism}")
                    mechanism_metrics = self.mechanism_ensembles[category][mechanism].fit(
                        mechanism,
                        compounds,
                        scores,
                    )
                    # Create metrics dict with prefixed keys
                    mechanism_metrics_dict = {}
                    for k, v in mechanism_metrics.items():
                        key = f"mechanism_{mechanism}_{k}"
                        mechanism_metrics_dict[key] = v
                    metrics.update(mechanism_metrics_dict)
        return metrics

    def _train_organ_ensembles(
        self,
        compounds: List[CompoundData],
        organ_data: Dict[str, Dict[str, List[float]]],
    ) -> Dict[str, float]:
        """Train organ toxicity ensembles."""
        metrics = {}
        for organ, toxicities in organ_data.items():
            for toxicity, severities in toxicities.items():
                if toxicity in self.organ_toxicity.get(organ, []):
                    self.logger.debug(f"Training organ ensemble for {toxicity}")
                    organ_metrics = self.organ_ensembles[organ][toxicity].fit(
                        toxicity,
                        compounds,
                        severities,
                    )
                    # Create metrics dict with prefixed keys
                    organ_metrics_dict = {}
                    for k, v in organ_metrics.items():
                        key = f"organ_{toxicity}_{k}"
                        organ_metrics_dict[key] = v
                    metrics.update(organ_metrics_dict)
        return metrics

    def _train_safety_ensembles(
        self,
        compounds: List[CompoundData],
        safety_data: Dict[str, Dict[str, List[float]]],
    ) -> Dict[str, float]:
        """Train safety concern ensembles."""
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

    def retrain(
        self,
        compounds: List[CompoundData],
        mechanism_labels: List[ToxicityClass],
        mechanism_data: Dict[str, Dict[str, List[float]]],
        organ_data: Dict[str, Dict[str, List[float]]],
        safety_data: Dict[str, Dict[str, List[float]]],
        **kwargs,
    ) -> Dict[str, float]:
        """Retrain ensemble models."""
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
