"""Enhanced nootropic effects predictor using ensemble models.

This module enhances the NootropicPredictor with:
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
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import brier_score_loss

from ....models.core import CompoundData
from ....models.psychopharm import NootropicMechanism
from ....models.compound.ml.ensemble import ClassifierEnsemble, RegressorEnsemble
from ..base import PredictionResult
from .nootropic import NootropicPredictor
from ..bbb.integration import BBBPredictorEnhanced

ImportanceDict: TypeAlias = Dict[str, Dict[str, float]]


class NootropicPredictorEnhanced(NootropicPredictor):
    """Enhanced nootropic predictor using ensemble models and BBB integration."""

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        cognitive_domains: Optional[Dict[str, Set[str]]] = None,
        side_effects: Optional[Dict[str, Set[str]]] = None,
        bbb_model_dir: Optional[str] = None,
    ):
        """Initialize enhanced nootropic predictor.

        Args:
            model_dir: Optional directory containing trained models
            cache_dir: Optional directory for caching
            log_level: Logging level
            cognitive_domains: Optional custom cognitive domains
            side_effects: Optional custom side effects
            bbb_model_dir: Optional directory containing BBB models
        """
        super().__init__(
            model_dir=model_dir,
            cache_dir=cache_dir,
            log_level=log_level,
            cognitive_domains=cognitive_domains,
            side_effects=side_effects,
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

        self.effect_ensembles = {}
        for domain, effects in self.cognitive_domains.items():
            domain_ensembles = {}
            for effect in effects:
                domain_ensembles[effect] = RegressorEnsemble(
                    model_dir=self.model_dir,
                    cache_dir=self.cache_dir,
                    feature_types=self.feature_types,
                    log_level=log_level,
                )
            self.effect_ensembles[domain] = domain_ensembles

        self.side_effect_ensembles = {}
        for category, effects in self.side_effects.items():
            category_ensembles = {}
            for effect in effects:
                category_ensembles[effect] = RegressorEnsemble(
                    model_dir=self.model_dir,
                    cache_dir=self.cache_dir,
                    feature_types=self.feature_types,
                    log_level=log_level,
                )
            self.side_effect_ensembles[category] = category_ensembles

        self.logger.info("NootropicPredictorEnhanced initialized successfully")

    def _load_models(self) -> Dict:
        """Load ensemble models."""
        self.logger.info("Loading ensemble models")

        try:
            # Load BBB predictor models
            self.bbb_predictor._load_models()

            # Load mechanism ensemble
            self.mechanism_ensemble.load_models()

            # Load effect ensembles
            for domain, ensembles in self.effect_ensembles.items():
                for effect, ensemble in ensembles.items():
                    ensemble.load_models()

            # Load side effect ensembles
            for category, ensembles in self.side_effect_ensembles.items():
                for effect, ensemble in ensembles.items():
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

        # Initialize effect ensembles
        for domain, effects in self.cognitive_domains.items():
            for effect in effects:
                self.effect_ensembles[domain][effect].create_ensemble(
                    effect,
                    ensemble_type="stacking",
                )

        # Initialize side effect ensembles
        for category, effects in self.side_effects.items():
            for effect in effects:
                self.side_effect_ensembles[category][effect].create_ensemble(
                    effect,
                    ensemble_type="stacking",
                )

        return {}  # Models stored in ensembles

    def save_models(self, save_dir: Optional[str] = None) -> None:
        """Save ensemble models."""
        save_dir = Path(save_dir) if save_dir else self.model_dir
        if not save_dir:
            save_dir = self.DEFAULT_MODEL_DIR

        save_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Saving ensemble models to {save_dir}")

        try:
            # Save BBB predictor models
            self.bbb_predictor.save_models(save_dir / "bbb")

            # Save mechanism ensemble
            self.mechanism_ensemble.save_models(save_dir)

            # Save effect ensembles
            for domain, ensembles in self.effect_ensembles.items():
                domain_dir = save_dir / domain
                for effect, ensemble in ensembles.items():
                    ensemble.save_models(domain_dir)

            # Save side effect ensembles
            for category, ensembles in self.side_effect_ensembles.items():
                category_dir = save_dir / category
                for effect, ensemble in ensembles.items():
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
        effect_data: Dict[str, Dict[str, List[float]]],
        side_effect_data: Dict[str, Dict[str, List[float]]],
        cv: int = 5,
    ) -> Dict[str, float]:
        """Calibrate model probabilities using held-out data.

        Args:
            compounds: Compounds for calibration
            mechanism_labels: True mechanism labels
            effect_data: Effect scores for calibration
            side_effect_data: Side effect scores for calibration
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

        # Calibrate effect ensembles
        for domain, effects in self.effect_ensembles.items():
            for effect, ensemble in effects.items():
                if effect in effect_data.get(domain, {}):
                    true_scores = effect_data[domain][effect]
                    calibrated = self._calibrate_regressor(ensemble, compounds, true_scores, cv=cv)
                    self.effect_ensembles[domain][effect] = calibrated

        # Calibrate side effect ensembles
        for category, effects in self.side_effect_ensembles.items():
            for effect, ensemble in effects.items():
                if effect in side_effect_data.get(category, {}):
                    true_scores = side_effect_data[category][effect]
                    calibrated = self._calibrate_regressor(ensemble, compounds, true_scores, cv=cv)
                    self.side_effect_ensembles[category][effect] = calibrated

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

    def _predict_effects(
        self,
        compound: CompoundData,
        bbb_result: PredictionResult,
    ) -> Dict[str, Dict[str, Dict[str, float]]]:
        """Predict cognitive effects using ensembles and BBB data."""
        effect_predictions = {}
        for domain, ensembles in self.effect_ensembles.items():
            domain_predictions = {}
            for effect, ensemble in ensembles.items():
                # Get base prediction
                score, conf = ensemble.predict(effect, compound)

                # Adjust based on BBB permeability
                if conf > 0.5:  # Base confidence threshold
                    # Scale effect by BBB permeability
                    adjusted_score = score * bbb_result.confidence
                    # Combine confidences
                    adjusted_conf = conf * bbb_result.confidence

                    if adjusted_conf > 0.3:  # Adjusted confidence threshold
                        domain_predictions[effect] = {
                            "score": adjusted_score,
                            "confidence": adjusted_conf,
                            "bbb_factor": bbb_result.confidence,
                        }
            if domain_predictions:
                effect_predictions[domain] = domain_predictions
        return effect_predictions

    def get_feature_importance(self, mechanism: Optional[str] = None) -> ImportanceDict:
        """Get feature importance scores for models.

        Args:
            mechanism: Optional specific mechanism to get importances for

        Returns:
            Dictionary of feature importance scores by model/mechanism
        """
        importances = {}

        # Get mechanism model importances
        if hasattr(self, "mechanism_importances"):
            importances["mechanism"] = self.mechanism_importances

        # Get effect model importances
        for domain, ensembles in self.effect_ensembles.items():
            for effect, ensemble in ensembles.items():
                if hasattr(ensemble, "feature_importances_"):
                    key = f"effect_{domain}_{effect}"
                    importances[key] = {"importance": ensemble.feature_importances_}

        # Get side effect model importances
        for category, ensembles in self.side_effect_ensembles.items():
            for effect, ensemble in ensembles.items():
                if hasattr(ensemble, "feature_importances_"):
                    key = f"side_effect_{category}_{effect}"
                    importances[key] = {"importance": ensemble.feature_importances_}

        return importances

    def _predict_side_effects(
        self,
        compound: CompoundData,
        bbb_result: PredictionResult,
    ) -> Dict[str, Dict[str, Dict[str, float]]]:
        """Predict side effects using ensembles and BBB data."""
        side_effect_predictions = {}
        for category, ensembles in self.side_effect_ensembles.items():
            category_predictions = {}
            for effect, ensemble in ensembles.items():
                # Get base prediction
                risk, conf = ensemble.predict(effect, compound)

                # Adjust based on BBB permeability
                if conf > 0.5:  # Base confidence threshold
                    # Scale risk by BBB permeability for CNS effects
                    if category in ["cognitive", "neurological"]:
                        adjusted_risk = risk * bbb_result.confidence
                        adjusted_conf = conf * bbb_result.confidence
                    else:
                        # Non-CNS effects less dependent on BBB
                        adjusted_risk = risk
                        adjusted_conf = conf

                    if adjusted_conf > 0.3:  # Adjusted confidence threshold
                        category_predictions[effect] = {
                            "risk": adjusted_risk,
                            "confidence": adjusted_conf,
                            "bbb_factor": bbb_result.confidence,
                        }
            if category_predictions:
                side_effect_predictions[category] = category_predictions
        return side_effect_predictions

    def predict(self, compound: CompoundData) -> PredictionResult:
        """Generate nootropic predictions using ensembles and BBB data."""
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

            # Predict cognitive effects with BBB integration
            effect_predictions = self._predict_effects(compound, bbb_result)

            # Update prediction history for effects
            for domain, effects in effect_predictions.items():
                for effect, pred in effects.items():
                    self._update_prediction_history(
                        compound=compound,
                        mechanism=mechanism,
                        mechanism_confidence=mechanism_confidence,
                        domain=domain,
                        effect=effect,
                        effect_score=pred["score"],
                        effect_confidence=pred["confidence"],
                        bbb_permeability=bbb_result.value.value,
                        bbb_confidence=bbb_result.confidence,
                    )

            # Predict side effects with BBB integration
            side_effect_predictions = self._predict_side_effects(compound, bbb_result)

            # Update prediction history for side effects
            for category, effects in side_effect_predictions.items():
                for effect, pred in effects.items():
                    self._update_prediction_history(
                        compound=compound,
                        mechanism=mechanism,
                        mechanism_confidence=mechanism_confidence,
                        side_effect_category=category,
                        side_effect=effect,
                        risk_score=pred["risk"],
                        risk_confidence=pred["confidence"],
                        bbb_permeability=bbb_result.value.value,
                        bbb_confidence=bbb_result.confidence,
                    )

            # Format result
            result = PredictionResult(
                value=NootropicMechanism(mechanism),
                confidence=mechanism_confidence,
                supporting_data={
                    "bbb_prediction": {
                        "value": bbb_result.value.value,
                        "confidence": bbb_result.confidence,
                        **bbb_result.supporting_data,
                    },
                    "effects": effect_predictions,
                    "side_effects": side_effect_predictions,
                    **self._format_supporting_data(
                        {},  # Features handled by ensembles
                        [
                            (effect, pred["confidence"])
                            for preds in effect_predictions.values()
                            for effect, pred in preds.items()
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
                value=NootropicMechanism.UNKNOWN,
                confidence=0.0,
                supporting_data={"error": str(e)},
            )

    def _update_prediction_history(
        self,
        compound: CompoundData,
        mechanism: str,
        mechanism_confidence: float,
        domain: Optional[str] = None,
        effect: Optional[str] = None,
        effect_score: Optional[float] = None,
        effect_confidence: Optional[float] = None,
        side_effect_category: Optional[str] = None,
        side_effect: Optional[str] = None,
        risk_score: Optional[float] = None,
        risk_confidence: Optional[float] = None,
        bbb_permeability: Optional[str] = None,
        bbb_confidence: Optional[float] = None,
    ) -> None:
        """Update prediction history with new prediction."""
        self.prediction_history = pd.concat(
            [
                self.prediction_history,
                pd.DataFrame(
                    [
                        {
                            "compound_name": compound.name,
                            "mechanism": mechanism,
                            "mechanism_confidence": mechanism_confidence,
                            "domain": domain,
                            "effect": effect,
                            "effect_score": effect_score,
                            "effect_confidence": effect_confidence,
                            "side_effect_category": side_effect_category,
                            "side_effect": side_effect,
                            "risk_score": risk_score,
                            "risk_confidence": risk_confidence,
                            "bbb_permeability": bbb_permeability,
                            "bbb_confidence": bbb_confidence,
                            "timestamp": pd.Timestamp.now(),
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

    def _train_effect_ensembles(
        self,
        compounds: List[CompoundData],
        effect_data: Dict[str, Dict[str, List[float]]],
    ) -> Dict[str, float]:
        """Train cognitive effect ensembles."""
        metrics = {}
        for domain, effects in effect_data.items():
            for effect, scores in effects.items():
                if effect in self.cognitive_domains.get(domain, []):
                    self.logger.debug(f"Training effect ensemble for {effect}")
                    effect_metrics = self.effect_ensembles[domain][effect].fit(
                        effect,
                        compounds,
                        scores,
                    )
                    # Create metrics dict with prefixed keys
                    effect_metrics_dict = {}
                    for k, v in effect_metrics.items():
                        key = f"effect_{effect}_{k}"
                        effect_metrics_dict[key] = v
                    metrics.update(effect_metrics_dict)
        return metrics

    def _train_side_effect_ensembles(
        self,
        compounds: List[CompoundData],
        side_effect_data: Dict[str, Dict[str, List[float]]],
    ) -> Dict[str, float]:
        """Train side effect ensembles."""
        metrics = {}
        for category, effects in side_effect_data.items():
            for effect, risks in effects.items():
                if effect in self.side_effects.get(category, []):
                    self.logger.debug(f"Training side effect ensemble for {effect}")
                    effect_metrics = self.side_effect_ensembles[category][effect].fit(
                        effect,
                        compounds,
                        risks,
                    )
                    # Create metrics dict with prefixed keys
                    side_effect_metrics_dict = {}
                    for k, v in effect_metrics.items():
                        key = f"side_effect_{effect}_{k}"
                        side_effect_metrics_dict[key] = v
                    metrics.update(side_effect_metrics_dict)
        return metrics

    def retrain(
        self,
        compounds: List[CompoundData],
        mechanism_labels: List[NootropicMechanism],
        effect_data: Dict[str, Dict[str, List[float]]],
        side_effect_data: Dict[str, Dict[str, List[float]]],
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

            # Train effect and side effect ensembles
            # Train and update metrics for effects and side effects
            effect_metrics = self._train_effect_ensembles(compounds, effect_data)
            side_effect_metrics = self._train_side_effect_ensembles(compounds, side_effect_data)
            metrics.update(effect_metrics)
            metrics.update(side_effect_metrics)

            # Save updated models
            self.save_models()

            self.logger.info("Ensemble model retraining completed successfully")
            return metrics

        except Exception as e:
            self.logger.error(f"Error retraining ensemble models: {str(e)}")
            return {"error": str(e)}
