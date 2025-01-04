"""Toxicity predictor.

This module provides the ToxicityPredictor class that:
1. Predicts toxicological risks and hazards
2. Predicts organ-specific toxicity
3. Predicts toxicity mechanisms
4. Uses ensemble of ML models for robust predictions
5. Provides confidence scores and supporting data
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler

from ....models.core import CompoundData
from ....models.psychopharm import ToxicityClass, ToxicityMechanism, OrganToxicity
from ..base import PredictionResult
from .base import PredictorBase


class ToxicityPredictor(PredictorBase):
    """Predict toxicological risks and safety concerns."""

    # Toxicity mechanisms
    TOXICITY_MECHANISMS = {
        "cellular": {
            "oxidative_stress", "mitochondrial_toxicity", "dna_damage",
            "protein_adducts", "lipid_peroxidation", "membrane_disruption",
            "enzyme_inhibition", "receptor_overstimulation",
        },
        "molecular": {
            "reactive_metabolites", "free_radical_formation", "protein_crosslinking",
            "covalent_binding", "ion_channel_blockade", "transporter_inhibition",
        },
        "systemic": {
            "immune_activation", "inflammation", "apoptosis_induction",
            "necrosis", "fibrosis", "organ_failure", "metabolic_disruption",
        },
    }

    # Organ toxicity types
    ORGAN_TOXICITY = {
        "liver": {
            "hepatocellular_damage", "cholestasis", "steatosis",
            "fibrosis", "enzyme_elevation", "metabolic_dysfunction",
        },
        "kidney": {
            "tubular_damage", "glomerular_damage", "crystal_formation",
            "filtration_impairment", "electrolyte_imbalance",
        },
        "heart": {
            "arrhythmia", "contractility_changes", "qt_prolongation",
            "conduction_abnormalities", "structural_changes",
        },
        "brain": {
            "neurotoxicity", "cognitive_impairment", "seizure_risk",
            "behavioral_changes", "neurotransmitter_imbalance",
        },
        "blood": {
            "bone_marrow_suppression", "hemolysis", "coagulation_changes",
            "platelet_dysfunction", "immune_suppression",
        },
    }

    # Safety concerns
    SAFETY_CONCERNS = {
        "acute": {
            "ld50", "acute_toxicity", "immediate_effects",
            "overdose_risk", "emergency_concerns",
        },
        "chronic": {
            "carcinogenicity", "mutagenicity", "reproductive_toxicity",
            "developmental_toxicity", "organ_damage",
        },
        "special": {
            "drug_interactions", "contraindications", "vulnerable_populations",
            "genetic_factors", "environmental_risks",
        },
    }

    # Default model paths
    DEFAULT_MODEL_DIR = Path("models/toxicity")
    MODEL_FILENAMES = {
        "class": "class_predictor.pkl",
        "mechanism": "mechanism_predictor.pkl",
        "organ": "organ_predictor.pkl",
    }

    def __init__(
        self,
        model_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        toxicity_mechanisms: Optional[Dict[str, Set[str]]] = None,
        organ_toxicity: Optional[Dict[str, Set[str]]] = None,
        safety_concerns: Optional[Dict[str, Set[str]]] = None,
    ):
        """Initialize toxicity predictor."""
        # Setup logging
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        
        # Add file handler if model_dir provided
        if model_dir:
            log_path = Path(model_dir) / "toxicity_predictor.log"
            fh = logging.FileHandler(log_path)
            fh.setLevel(log_level)
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            fh.setFormatter(formatter)
            self.logger.addHandler(fh)

        # Use custom categories if provided
        self.toxicity_mechanisms = toxicity_mechanisms or self.TOXICITY_MECHANISMS
        self.organ_toxicity = organ_toxicity or self.ORGAN_TOXICITY
        self.safety_concerns = safety_concerns or self.SAFETY_CONCERNS
        
        # Initialize base class
        super().__init__(
            model_dir=model_dir,
            cache_dir=cache_dir,
            feature_types=["fingerprints", "descriptors", "enhanced"],
        )

        # Initialize feature scalers
        self.scalers = self._initialize_scalers()

        # Initialize prediction history
        self.prediction_history = pd.DataFrame(
            columns=[
                'compound_name',
                'toxicity_class',
                'class_confidence',
                'mechanism_category',
                'mechanism',
                'mechanism_score',
                'mechanism_confidence',
                'organ',
                'toxicity_type',
                'severity',
                'organ_confidence',
                'concern_category',
                'concern',
                'risk_level',
                'risk_confidence',
                'timestamp',
            ]
        )

        self.logger.info("ToxicityPredictor initialized successfully")

    def _initialize_scalers(self) -> Dict[str, StandardScaler]:
        """Initialize feature scalers for each feature type."""
        self.logger.debug("Initializing feature scalers")
        return {
            feature_type: StandardScaler()
            for feature_type in self.feature_types
        }

    def _load_models(self) -> Dict:
        """Load toxicity prediction models from disk."""
        models = {}
        model_dir = Path(self.model_dir) if self.model_dir else self.DEFAULT_MODEL_DIR

        if model_dir.exists():
            self.logger.info(f"Loading models from {model_dir}")
            try:
                # Load toxicity class model
                class_path = model_dir / "class_predictor.pkl"
                if class_path.exists():
                    self.logger.debug("Loading toxicity class model")
                    models["class"] = np.load(class_path, allow_pickle=True)
                else:
                    self.logger.warning(
                        "Class predictor not found, initializing new model"
                    )
                    models["class"] = self._initialize_class_model()

                # Load mechanism prediction models
                mechanism_models = {}
                for category, mechanisms in self.toxicity_mechanisms.items():
                    category_models = {}
                    for mechanism in mechanisms:
                        model_path = model_dir / f"mechanism_{mechanism}.pkl"
                        if model_path.exists():
                            self.logger.debug(f"Loading mechanism model for {mechanism}")
                            category_models[mechanism] = np.load(
                                model_path, allow_pickle=True
                            )
                        else:
                            self.logger.warning(
                                f"Model not found for {mechanism}, initializing new model"
                            )
                            category_models[mechanism] = self._initialize_mechanism_model()
                    mechanism_models[category] = category_models
                models["mechanisms"] = mechanism_models

                # Load organ toxicity models
                organ_models = {}
                for organ, toxicities in self.organ_toxicity.items():
                    organ_type_models = {}
                    for toxicity in toxicities:
                        model_path = model_dir / f"organ_{organ}_{toxicity}.pkl"
                        if model_path.exists():
                            self.logger.debug(f"Loading organ model for {organ}_{toxicity}")
                            organ_type_models[toxicity] = np.load(
                                model_path, allow_pickle=True
                            )
                        else:
                            self.logger.warning(
                                f"Model not found for {organ}_{toxicity}, initializing new model"
                            )
                            organ_type_models[toxicity] = self._initialize_organ_model()
                    organ_models[organ] = organ_type_models
                models["organs"] = organ_models

                # Load safety concern models
                concern_models = {}
                for category, concerns in self.safety_concerns.items():
                    category_models = {}
                    for concern in concerns:
                        model_path = model_dir / f"concern_{concern}.pkl"
                        if model_path.exists():
                            self.logger.debug(f"Loading concern model for {concern}")
                            category_models[concern] = np.load(
                                model_path, allow_pickle=True
                            )
                        else:
                            self.logger.warning(
                                f"Model not found for {concern}, initializing new model"
                            )
                            category_models[concern] = self._initialize_concern_model()
                    concern_models[category] = category_models
                models["concerns"] = concern_models

            except Exception as e:
                self.logger.error(f"Error loading models: {str(e)}")
                self.logger.info("Initializing new models")
                models = self._initialize_models()
        else:
            self.logger.info(f"Model directory not found: {model_dir}")
            self.logger.info("Initializing new models")
            models = self._initialize_models()

        return models

    def _initialize_class_model(self) -> RandomForestClassifier:
        """Initialize new toxicity class prediction model."""
        return RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            random_state=42,
            n_jobs=-1,
            verbose=1,
        )

    def _initialize_mechanism_model(self) -> GradientBoostingRegressor:
        """Initialize new mechanism prediction model."""
        return GradientBoostingRegressor(
            n_estimators=100,
            max_depth=5,
            random_state=42,
            verbose=1,
        )

    def _initialize_organ_model(self) -> GradientBoostingRegressor:
        """Initialize new organ toxicity prediction model."""
        return GradientBoostingRegressor(
            n_estimators=100,
            max_depth=5,
            random_state=42,
            verbose=1,
        )

    def _initialize_concern_model(self) -> GradientBoostingRegressor:
        """Initialize new safety concern prediction model."""
        return GradientBoostingRegressor(
            n_estimators=100,
            max_depth=5,
            random_state=42,
            verbose=1,
        )

    def _initialize_models(self) -> Dict:
        """Initialize all toxicity prediction models."""
        self.logger.info("Initializing new models")
        
        models = {
            "class": self._initialize_class_model(),
            "mechanisms": {},
            "organs": {},
            "concerns": {},
        }

        # Initialize mechanism models
        for category, mechanisms in self.toxicity_mechanisms.items():
            category_models = {
                mechanism: self._initialize_mechanism_model()
                for mechanism in mechanisms
            }
            models["mechanisms"][category] = category_models

        # Initialize organ toxicity models
        for organ, toxicities in self.organ_toxicity.items():
            organ_type_models = {
                toxicity: self._initialize_organ_model()
                for toxicity in toxicities
            }
            models["organs"][organ] = organ_type_models

        # Initialize safety concern models
        for category, concerns in self.safety_concerns.items():
            category_models = {
                concern: self._initialize_concern_model()
                for concern in concerns
            }
            models["concerns"][category] = category_models

        return models

    def save_models(self, save_dir: Optional[str] = None) -> None:
        """Save models to disk with versioning."""
        save_dir = Path(save_dir) if save_dir else self.model_dir
        if not save_dir:
            save_dir = self.DEFAULT_MODEL_DIR

        save_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Saving models to {save_dir}")

        # Save models with versioning
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        version_dir = save_dir / f"version_{timestamp}"
        version_dir.mkdir(exist_ok=True)

        try:
            # Save toxicity class model
            class_path = save_dir / "class_predictor.pkl"
            np.save(class_path, self.models["class"])
            version_path = version_dir / "class_predictor.pkl"
            np.save(version_path, self.models["class"])
            self.logger.debug("Saved toxicity class model")

            # Save mechanism prediction models
            for category, mechanisms in self.models["mechanisms"].items():
                for mechanism, model in mechanisms.items():
                    # Save current version
                    model_path = save_dir / f"mechanism_{mechanism}.pkl"
                    np.save(model_path, model)
                    
                    # Save versioned copy
                    version_path = version_dir / f"mechanism_{mechanism}.pkl"
                    np.save(version_path, model)
                    
                    self.logger.debug(f"Saved mechanism model for {mechanism}")

            # Save organ toxicity models
            for organ, toxicities in self.models["organs"].items():
                for toxicity, model in toxicities.items():
                    # Save current version
                    model_path = save_dir / f"organ_{organ}_{toxicity}.pkl"
                    np.save(model_path, model)
                    
                    # Save versioned copy
                    version_path = version_dir / f"organ_{organ}_{toxicity}.pkl"
                    np.save(version_path, model)
                    
                    self.logger.debug(f"Saved organ model for {organ}_{toxicity}")

            # Save safety concern models
            for category, concerns in self.models["concerns"].items():
                for concern, model in concerns.items():
                    # Save current version
                    model_path = save_dir / f"concern_{concern}.pkl"
                    np.save(model_path, model)
                    
                    # Save versioned copy
                    version_path = version_dir / f"concern_{concern}.pkl"
                    np.save(version_path, model)
                    
                    self.logger.debug(f"Saved concern model for {concern}")

            # Save scalers
            scaler_path = save_dir / "scalers.pkl"
            np.save(scaler_path, self.scalers)
            self.logger.debug("Saved feature scalers")

            # Save prediction history
            history_path = save_dir / "prediction_history.csv"
            self.prediction_history.to_csv(history_path, index=False)
            self.logger.debug("Saved prediction history")

        except Exception as e:
            self.logger.error(f"Error saving models: {str(e)}")

    def predict(self, compound: CompoundData) -> PredictionResult:
        """Generate toxicity predictions for a compound."""
        self.logger.debug(f"Generating predictions for {compound.name}")
        try:
            # Extract features
            features = {}
            for feature_type in self.feature_types:
                features[feature_type] = self._extract_features(
                    compound, feature_type
                )
                self.logger.debug(
                    f"Extracted {feature_type} features: shape={features[feature_type].shape}"
                )

            # Predict toxicity class
            tox_class, class_confidence = self._predict_class(
                self.models["class"],
                features,
            )
            self.logger.debug(
                f"Toxicity class: {tox_class} (confidence: {class_confidence:.3f})"
            )

            # Predict toxicity mechanisms
            mechanism_predictions = {}
            for category, mechanisms in self.toxicity_mechanisms.items():
                category_predictions = {}
                for mechanism in mechanisms:
                    # Predict mechanism score
                    score, conf = self._predict_mechanism(
                        self.models["mechanisms"][category][mechanism],
                        features,
                        mechanism,
                    )
                    
                    if conf > 0.5:  # Confidence threshold
                        category_predictions[mechanism] = {
                            "score": score,
                            "confidence": conf,
                        }

                        # Update prediction history
                        self.prediction_history = pd.concat([
                            self.prediction_history,
                            pd.DataFrame([{
                                'compound_name': compound.name,
                                'toxicity_class': tox_class,
                                'class_confidence': class_confidence,
                                'mechanism_category': category,
                                'mechanism': mechanism,
                                'mechanism_score': score,
                                'mechanism_confidence': conf,
                                'organ': None,
                                'toxicity_type': None,
                                'severity': None,
                                'organ_confidence': None,
                                'concern_category': None,
                                'concern': None,
                                'risk_level': None,
                                'risk_confidence': None,
                                'timestamp': pd.Timestamp.now(),
                            }])
                        ], ignore_index=True)

                if category_predictions:
                    mechanism_predictions[category] = category_predictions

            # Predict organ toxicity
            organ_predictions = {}
            for organ, toxicities in self.organ_toxicity.items():
                organ_type_predictions = {}
                for toxicity in toxicities:
                    # Predict toxicity severity
                    severity, conf = self._predict_organ_toxicity(
                        self.models["organs"][organ][toxicity],
                        features,
                        organ,
                        toxicity,
                    )
                    
                    if conf > 0.5:  # Confidence threshold
                        organ_type_predictions[toxicity] = {
                            "severity": severity,
                            "confidence": conf,
                        }

                        # Update prediction history
                        self.prediction_history = pd.concat([
                            self.prediction_history,
                            pd.DataFrame([{
                                'compound_name': compound.name,
                                'toxicity_class': tox_class,
                                'class_confidence': class_confidence,
                                'mechanism_category': None,
                                'mechanism': None,
                                'mechanism_score': None,
                                'mechanism_confidence': None,
                                'organ': organ,
                                'toxicity_type': toxicity,
                                'severity': severity,
                                'organ_confidence': conf,
                                'concern_category': None,
                                'concern': None,
                                'risk_level': None,
                                'risk_confidence': None,
                                'timestamp': pd.Timestamp.now(),
                            }])
                        ], ignore_index=True)

                if organ_type_predictions:
                    organ_predictions[organ] = organ_type_predictions

            # Predict safety concerns
            concern_predictions = {}
            for category, concerns in self.safety_concerns.items():
                category_predictions = {}
                for concern in concerns:
                    # Predict risk level
                    risk, conf = self._predict_safety_concern(
                        self.models["concerns"][category][concern],
                        features,
                        concern,
                    )
                    
                    if conf > 0.5:  # Confidence threshold
                        category_predictions[concern] = {
                            "risk": risk,
                            "confidence": conf,
                        }

                        # Update prediction history
                        self.prediction_history = pd.concat([
                            self.prediction_history,
                            pd.DataFrame([{
                                'compound_name': compound.name,
                                'toxicity_class': tox_class,
                                'class_confidence': class_confidence,
                                'mechanism_category': None,
                                'mechanism': None,
                                'mechanism_score': None,
                                'mechanism_confidence': None,
                                'organ': None,
                                'toxicity_type': None,
                                'severity': None,
                                'organ_confidence': None,
                                'concern_category': category,
                                'concern': concern,
                                'risk_level': risk,
                                'risk_confidence': conf,
                                'timestamp': pd.Timestamp.now(),
                            }])
                        ], ignore_index=True)

                if category_predictions:
                    concern_predictions[category] = category_predictions

            # Format result
            result = PredictionResult(
                value=ToxicityClass(tox_class),
                confidence=class_confidence,
                supporting_data={
                    "mechanisms": mechanism_predictions,
                    "organs": organ_predictions,
                    "concerns": concern_predictions,
                    **self._format_supporting_data(
                        features,
                        [(mechanism, pred["confidence"])
                         for preds in mechanism_predictions.values()
                         for mechanism, pred in preds.items()],
                    ),
                },
            )
            
            self.logger.info(
                f"Generated predictions for {compound.name}: "
                f"{result.value.value} (confidence: {result.confidence:.3f})"
            )
            
            return result

        except Exception as e:
            self.logger.error(
                "Error predicting toxicity",
                f"Compound {compound.name}: {str(e)}",
                exc_info=True
            )
            return PredictionResult(
                value=ToxicityClass.UNKNOWN,
                confidence=0.0,
                supporting_data={"error": str(e)},
            )

    def _predict_class(
        self,
        model: RandomForestClassifier,
        features: Dict[str, np.ndarray],
    ) -> Tuple[str, float]:
        """Predict toxicity class."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers["class"].transform(X)

        # Get class probabilities
        probs = model.predict_proba(X_scaled)[0]
        pred_idx = np.argmax(probs)
        
        # Map to ToxicityClass
        toxicity_class = ToxicityClass(
            model.classes_[pred_idx]
        ).value
        
        return toxicity_class, float(probs[pred_idx])

    def _predict_mechanism(
        self,
        model: GradientBoostingRegressor,
        features: Dict[str, np.ndarray],
        mechanism: str,
    ) -> Tuple[float, float]:
        """Predict toxicity mechanism score."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers[f"mechanism_{mechanism}"].transform(X)

        # Get prediction and confidence
        score = model.predict(X_scaled)[0]
        confidence = 1.0 - np.std(
            [est.predict(X_scaled)[0] for est in model.estimators_]
        )

        return float(score), float(confidence)

    def _predict_organ_toxicity(
        self,
        model: GradientBoostingRegressor,
        features: Dict[str, np.ndarray],
        organ: str,
        toxicity: str,
    ) -> Tuple[float, float]:
        """Predict organ toxicity severity."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers[f"organ_{organ}_{toxicity}"].transform(X)

        # Get prediction and confidence
        severity = model.predict(X_scaled)[0]
        confidence = 1.0 - np.std(
            [est.predict(X_scaled)[0] for est in model.estimators_]
        )

        return float(severity), float(confidence)

    def _predict_safety_concern(
        self,
        model: GradientBoostingRegressor,
        features: Dict[str, np.ndarray],
        concern: str,
    ) -> Tuple[float, float]:
        """Predict safety concern risk level."""
        # Combine and scale features
        X = np.hstack([features[ft] for ft in self.feature_types])
        X_scaled = self.scalers[f"concern_{concern}"].transform(X)

        # Get prediction and confidence
        risk = model.predict(X_scaled)[0]
        confidence = 1.0 - np.std(
            [est.predict(X_scaled)[0] for est in model.estimators_]
        )

        return float(risk), float(confidence)

    def get_prediction_statistics(self) -> pd.DataFrame:
        """Get statistics about predictions made so far."""
        stats = pd.DataFrame()
        
        # Toxicity class distribution
        stats['class_dist'] = (
            self.prediction_history['toxicity_class'].value_counts(normalize=True)
        )
        
        # Average class confidence
        stats['class_confidence'] = (
            self.prediction_history.groupby('toxicity_class')['class_confidence'].mean()
        )
        
        # Mechanism score distribution
        stats['mechanism_score'] = (
            self.prediction_history.groupby(['mechanism_category', 'mechanism'])['mechanism_score'].mean()
        )
        
        # Mechanism confidence distribution
        stats['mechanism_confidence'] = (
            self.prediction_history.groupby(['mechanism_category', 'mechanism'])['mechanism_confidence'].mean()
        )
        
        # Organ toxicity severity distribution
        stats['organ_severity'] = (
            self.prediction_history.groupby(['organ', 'toxicity_type'])['severity'].mean()
        )
        
        # Organ toxicity confidence distribution
        stats['organ_confidence'] = (
            self.prediction_history.groupby(['organ', 'toxicity_type'])['organ_confidence'].mean()
        )
        
        # Safety concern risk distribution
        stats['risk_level'] = (
            self.prediction_history.groupby(['concern_category', 'concern'])['risk_level'].mean()
        )
        
        # Safety concern confidence distribution
        stats['risk_confidence'] = (
            self.prediction_history.groupby(['concern_category', 'concern'])['risk_confidence'].mean()
        )
        
        return stats

    def retrain(
        self,
        compounds: List[CompoundData],
        class_labels: List[ToxicityClass],
        mechanism_data: Dict[str, Dict[str, List[float]]],
        organ_data: Dict[str, Dict[str, List[float]]],
        concern_data: Dict[str, Dict[str, List[float]]],
        **kwargs,
    ) -> Dict[str, float]:
        """Retrain models with new data."""
        self.logger.info(
            f"Retraining models with {len(compounds)} compounds"
        )
        
        # Extract features
        X = []
        for compound in compounds:
            features = []
            for feature_type in self.feature_types:
                feat = self._extract_features(compound, feature_type)
                features.append(feat)
            X.append(np.hstack(features))
        X = np.vstack(X)

        metrics = {}

        # Train toxicity class model
        self.logger.debug("Training toxicity class model")
        self.models["class"].fit(X, class_labels)
        score = self.models["class"].score(X, class_labels)
        metrics["class_score"] = score
        self.logger.debug(f"Toxicity class prediction score: {score:.3f}")

        # Train mechanism prediction models
        for category, mechanisms in mechanism_data.items():
            for mechanism, scores in mechanisms.items():
                if mechanism in self.toxicity_mechanisms.get(category, []):
                    self.logger.debug(f"Training mechanism model for {mechanism}")
                    model = self.models["mechanisms"][category][mechanism]
                    model.fit(X, scores)
                    score = model.score(X, scores)
                    metrics[f"mechanism_{mechanism}_score"] = score
                    self.logger.debug(
                        f"Mechanism model for {mechanism} training score: {score:.3f}"
                    )

        # Train organ toxicity models
        for organ, toxicities in organ_data.items():
            for toxicity, severities in toxicities.items():
                if toxicity in self.organ_toxicity.get(organ, []):
                    self.logger.debug(f"Training organ model for {organ}_{toxicity}")
                    model = self.models["organs"][organ][toxicity]
                    model.fit(X, severities)
                    score = model.score(X, severities)
                    metrics[f"organ_{organ}_{toxicity}_score"] = score
                    self.logger.debug(
                        f"Organ model for {organ}_{toxicity} training score: {score:.3f}"
                    )

        # Train safety concern models
        for category, concerns in concern_data.items():
            for concern, risks in concerns.items():
                if concern in self.safety_concerns.get(category, []):
                    self.logger.debug(f"Training concern model for {concern}")
                    model = self.models["concerns"][category][concern]
                    model.fit(X, risks)
                    score = model.score(X, risks)
                    metrics[f"concern_{concern}_score"] = score
                    self.logger.debug(
                        f"Concern model for {concern} training score: {score:.3f}"
                    )

        # Save updated models
        self.save_models()
        
        self.logger.info("Model retraining completed successfully")
        return metrics
