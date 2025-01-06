"""Base toxicity predictor.

This module provides the core functionality for predicting toxicological properties
of chemical compounds. It handles:
1. Model loading and management
2. Feature extraction and scaling
3. Basic toxicity predictions
4. Prediction history tracking
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Set

import pandas as pd

from .....models.core import CompoundData
from .....models.psychopharm import ToxicityClass
from ...base import PredictorBase
from . import model_loading, prediction, features


class ToxicityPredictorBase(PredictorBase):
    """Base class for toxicity prediction."""

    # Default categories
    TOXICITY_MECHANISMS = {
        "cellular": {
            "oxidative_stress",
            "mitochondrial_toxicity",
            "dna_damage",
            "protein_adducts",
            "lipid_peroxidation",
            "membrane_disruption",
            "enzyme_inhibition",
            "receptor_overstimulation",
        },
        "molecular": {
            "reactive_metabolites",
            "free_radical_formation",
            "protein_crosslinking",
            "covalent_binding",
            "ion_channel_blockade",
            "transporter_inhibition",
        },
        "systemic": {
            "immune_activation",
            "inflammation",
            "apoptosis_induction",
            "necrosis",
            "fibrosis",
            "organ_failure",
            "metabolic_disruption",
        },
    }

    ORGAN_TOXICITY = {
        "liver": {
            "hepatocellular_damage",
            "cholestasis",
            "steatosis",
            "fibrosis",
            "enzyme_elevation",
            "metabolic_dysfunction",
        },
        "kidney": {
            "tubular_damage",
            "glomerular_damage",
            "crystal_formation",
            "filtration_impairment",
            "electrolyte_imbalance",
        },
        "heart": {
            "arrhythmia",
            "contractility_changes",
            "qt_prolongation",
            "conduction_abnormalities",
            "structural_changes",
        },
        "brain": {
            "neurotoxicity",
            "cognitive_impairment",
            "seizure_risk",
            "behavioral_changes",
            "neurotransmitter_imbalance",
        },
        "blood": {
            "bone_marrow_suppression",
            "hemolysis",
            "coagulation_changes",
            "platelet_dysfunction",
            "immune_suppression",
        },
    }

    SAFETY_CONCERNS = {
        "acute": {
            "ld50",
            "acute_toxicity",
            "immediate_effects",
            "overdose_risk",
            "emergency_concerns",
        },
        "chronic": {
            "carcinogenicity",
            "mutagenicity",
            "reproductive_toxicity",
            "developmental_toxicity",
            "organ_damage",
        },
        "special": {
            "drug_interactions",
            "contraindications",
            "vulnerable_populations",
            "genetic_factors",
            "environmental_risks",
        },
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
        """Initialize toxicity predictor.

        Args:
            model_dir: Directory containing model files
            cache_dir: Directory for caching predictions
            log_level: Logging level
            toxicity_mechanisms: Custom toxicity mechanism categories
            organ_toxicity: Custom organ toxicity categories
            safety_concerns: Custom safety concern categories
        """
        # Setup logging
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)

        # Add file handler if model_dir provided
        if model_dir:
            log_path = Path(model_dir) / "toxicity_predictor.log"
            fh = logging.FileHandler(log_path)
            fh.setLevel(log_level)
            formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
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

        # Initialize models and scalers
        self.models = model_loading.load_models(
            model_dir=self.model_dir,
            toxicity_mechanisms=self.toxicity_mechanisms,
            organ_toxicity=self.organ_toxicity,
            safety_concerns=self.safety_concerns,
            logger=self.logger,
        )
        self.scalers = model_loading.initialize_scalers(
            feature_types=self.feature_types,
            toxicity_mechanisms=self.toxicity_mechanisms,
            organ_toxicity=self.organ_toxicity,
            safety_concerns=self.safety_concerns,
            logger=self.logger,
        )

        # Initialize prediction history
        self.prediction_history = pd.DataFrame(
            columns=[
                "compound_name",
                "toxicity_class",
                "class_confidence",
                "mechanism_category",
                "mechanism",
                "mechanism_score",
                "mechanism_confidence",
                "organ",
                "toxicity_type",
                "severity",
                "organ_confidence",
                "concern_category",
                "concern",
                "risk_level",
                "risk_confidence",
                "timestamp",
            ]
        )

        self.logger.info("ToxicityPredictorBase initialized successfully")

    def predict(self, compound: CompoundData) -> Dict:
        """Generate toxicity predictions for a compound.

        Args:
            compound: Compound to predict toxicity for

        Returns:
            Dictionary containing:
            - toxicity_class: Overall toxicity classification
            - class_confidence: Confidence in classification
            - mechanisms: Predicted toxicity mechanisms
            - organs: Predicted organ-specific toxicity
            - concerns: Predicted safety concerns
        """
        self.logger.debug(f"Generating predictions for {compound.name}")

        try:
            # Extract features
            compound_features = {}
            for feature_type in self.feature_types:
                compound_features[feature_type] = features.extract_features(compound, feature_type)
                self.logger.debug(
                    f"Extracted {feature_type} features: " f"shape={compound_features[feature_type].shape}"
                )

            # Predict toxicity class
            tox_class, class_confidence = prediction.predict_class(
                self.models["class"],
                compound_features,
                self.feature_types,
                self.scalers,
                self.logger,
            )
            self.logger.debug(f"Toxicity class: {tox_class} (confidence: {class_confidence:.3f})")

            # Predict mechanisms
            mechanism_predictions = prediction.predict_all_mechanisms(
                compound,
                self.models,
                self.toxicity_mechanisms,
                self.feature_types,
                self.scalers,
                tox_class,
                class_confidence,
                self.logger,
            )

            # Predict organ toxicity
            organ_predictions = prediction.predict_all_organ_effects(
                compound,
                self.models,
                self.organ_toxicity,
                self.feature_types,
                self.scalers,
                tox_class,
                class_confidence,
                self.logger,
            )

            # Predict safety concerns
            concern_predictions = prediction.predict_all_safety_concerns(
                compound,
                self.models,
                self.safety_concerns,
                self.feature_types,
                self.scalers,
                tox_class,
                class_confidence,
                self.logger,
            )

            # Update prediction history
            self._update_history(
                compound,
                tox_class,
                class_confidence,
                mechanism_predictions,
                organ_predictions,
                concern_predictions,
            )

            return {
                "toxicity_class": ToxicityClass(tox_class),
                "class_confidence": class_confidence,
                "mechanisms": mechanism_predictions,
                "organs": organ_predictions,
                "concerns": concern_predictions,
            }

        except Exception as e:
            self.logger.error(
                f"Error predicting toxicity for {compound.name}: {str(e)}",
                exc_info=True,
            )
            return {
                "toxicity_class": ToxicityClass.UNKNOWN,
                "class_confidence": 0.0,
                "error": str(e),
            }

    def _update_history(
        self,
        compound: CompoundData,
        tox_class: str,
        class_confidence: float,
        mechanism_predictions: Dict,
        organ_predictions: Dict,
        concern_predictions: Dict,
    ) -> None:
        """Update prediction history with new predictions."""
        # Add mechanism predictions
        for category, mechanisms in mechanism_predictions.items():
            for mechanism, data in mechanisms.items():
                self.prediction_history = pd.concat(
                    [
                        self.prediction_history,
                        pd.DataFrame(
                            [
                                {
                                    "compound_name": compound.name,
                                    "toxicity_class": tox_class,
                                    "class_confidence": class_confidence,
                                    "mechanism_category": category,
                                    "mechanism": mechanism,
                                    "mechanism_score": data["score"],
                                    "mechanism_confidence": data["confidence"],
                                    "organ": None,
                                    "toxicity_type": None,
                                    "severity": None,
                                    "organ_confidence": None,
                                    "concern_category": None,
                                    "concern": None,
                                    "risk_level": None,
                                    "risk_confidence": None,
                                    "timestamp": pd.Timestamp.now(),
                                }
                            ]
                        ),
                    ],
                    ignore_index=True,
                )

        # Add organ toxicity predictions
        for organ, toxicities in organ_predictions.items():
            for toxicity, data in toxicities.items():
                self.prediction_history = pd.concat(
                    [
                        self.prediction_history,
                        pd.DataFrame(
                            [
                                {
                                    "compound_name": compound.name,
                                    "toxicity_class": tox_class,
                                    "class_confidence": class_confidence,
                                    "mechanism_category": None,
                                    "mechanism": None,
                                    "mechanism_score": None,
                                    "mechanism_confidence": None,
                                    "organ": organ,
                                    "toxicity_type": toxicity,
                                    "severity": data["severity"],
                                    "organ_confidence": data["confidence"],
                                    "concern_category": None,
                                    "concern": None,
                                    "risk_level": None,
                                    "risk_confidence": None,
                                    "timestamp": pd.Timestamp.now(),
                                }
                            ]
                        ),
                    ],
                    ignore_index=True,
                )

        # Add safety concern predictions
        for category, concerns in concern_predictions.items():
            for concern, data in concerns.items():
                self.prediction_history = pd.concat(
                    [
                        self.prediction_history,
                        pd.DataFrame(
                            [
                                {
                                    "compound_name": compound.name,
                                    "toxicity_class": tox_class,
                                    "class_confidence": class_confidence,
                                    "mechanism_category": None,
                                    "mechanism": None,
                                    "mechanism_score": None,
                                    "mechanism_confidence": None,
                                    "organ": None,
                                    "toxicity_type": None,
                                    "severity": None,
                                    "organ_confidence": None,
                                    "concern_category": category,
                                    "concern": concern,
                                    "risk_level": data["risk"],
                                    "risk_confidence": data["confidence"],
                                    "timestamp": pd.Timestamp.now(),
                                }
                            ]
                        ),
                    ],
                    ignore_index=True,
                )

    def get_prediction_statistics(self) -> pd.DataFrame:
        """Get statistics about predictions made so far."""
        stats = pd.DataFrame()

        # Toxicity class distribution
        stats["class_dist"] = self.prediction_history["toxicity_class"].value_counts(normalize=True)

        # Average class confidence
        stats["class_confidence"] = self.prediction_history.groupby("toxicity_class")["class_confidence"].mean()

        # Mechanism score distribution
        stats["mechanism_score"] = self.prediction_history.groupby(["mechanism_category", "mechanism"])[
            "mechanism_score"
        ].mean()

        # Mechanism confidence distribution
        stats["mechanism_confidence"] = self.prediction_history.groupby(["mechanism_category", "mechanism"])[
            "mechanism_confidence"
        ].mean()

        # Organ toxicity severity distribution
        stats["organ_severity"] = self.prediction_history.groupby(["organ", "toxicity_type"])["severity"].mean()

        # Organ toxicity confidence distribution
        stats["organ_confidence"] = self.prediction_history.groupby(["organ", "toxicity_type"])[
            "organ_confidence"
        ].mean()

        # Safety concern risk distribution
        stats["risk_level"] = self.prediction_history.groupby(["concern_category", "concern"])["risk_level"].mean()

        # Safety concern confidence distribution
        stats["risk_confidence"] = self.prediction_history.groupby(["concern_category", "concern"])[
            "risk_confidence"
        ].mean()

        return stats
