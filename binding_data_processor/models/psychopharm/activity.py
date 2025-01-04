"""Psychoactive and nootropic activity analysis functionality."""

from dataclasses import field
from typing import Dict, List, Set

import numpy as np

from .base import (
    PsychoactiveClass,
    NootropicMechanism,
    EffectScore,
    RiskLevel,
)


class ActivityProfileMixin:
    """Mixin class providing psychoactive and nootropic activity analysis."""

    # Psychoactive classification
    psychoactive_class: PsychoactiveClass = PsychoactiveClass.UNKNOWN
    secondary_classes: Set[PsychoactiveClass] = field(default_factory=set)
    effect_profile: Dict[str, EffectScore] = field(default_factory=dict)
    
    # Nootropic properties
    nootropic_mechanisms: Set[NootropicMechanism] = field(default_factory=set)
    cognitive_effects: Dict[str, EffectScore] = field(default_factory=dict)
    side_effects: Dict[str, EffectScore] = field(default_factory=dict)

    # Effect mappings
    effect_mechanisms: Dict[str, List[str]] = {
        "psychedelic": [
            "5-HT2A_agonist",
            "5-HT2C_agonist",
        ],
        "empathogen": [
            "SERT_inhibitor",
            "5-HT1A_agonist",
            "5-HT2A_partial",
        ],
        "stimulant": [
            "DAT_inhibitor",
            "NET_inhibitor",
            "D1_agonist",
            "D2_agonist",
        ],
        "dissociative": [
            "NMDA_antagonist",
            "kappa_agonist",
        ],
        "anxiolytic": [
            "GABA-A_modulator",
            "5-HT1A_agonist",
        ],
        "antipsychotic": [
            "D2_antagonist",
            "5-HT2A_antagonist",
        ],
        "cognitive_enhancement": [
            "acetylcholine_increase",
            "glutamate_modulation",
            "dopamine_modulation",
        ],
        "mood_enhancement": [
            "serotonin_increase",
            "dopamine_increase",
            "endorphin_increase",
        ],
    }

    # Nootropic mechanism mappings
    nootropic_pathways: Dict[NootropicMechanism, List[str]] = {
        NootropicMechanism.CHOLINERGIC: [
            "nAChR_agonist",
            "mAChR_agonist",
            "AChE_inhibitor",
        ],
        NootropicMechanism.GLUTAMATERGIC: [
            "AMPA_modulator",
            "NMDA_modulator",
            "mGluR_modulator",
        ],
        NootropicMechanism.DOPAMINERGIC: [
            "D1_agonist",
            "DAT_inhibitor",
            "COMT_inhibitor",
        ],
        NootropicMechanism.SEROTONERGIC: [
            "5-HT1A_agonist",
            "5-HT6_antagonist",
            "5-HT7_agonist",
        ],
        NootropicMechanism.GABA: [
            "GABA-A_modulator",
            "GABA-B_modulator",
        ],
        NootropicMechanism.AMPAKINE: [
            "AMPA_positive_modulator",
        ],
        NootropicMechanism.BDNF: [
            "BDNF_increase",
            "TrkB_activation",
        ],
        NootropicMechanism.NGF: [
            "NGF_increase",
            "TrkA_activation",
        ],
    }

    def get_activity_dict(self) -> Dict:
        """Get dictionary of activity data."""
        return {
            "psychoactive": self._format_psychoactive_data(),
            "nootropic": self._format_nootropic_data(),
            "effects": self._analyze_effects(),
            "mechanisms": self._analyze_mechanisms(),
            "predictions": self._predict_activities(),
        }

    def _format_psychoactive_data(self) -> Dict:
        """Format psychoactive classification data."""
        return {
            "primary_class": self.psychoactive_class.value,
            "secondary_classes": [cls.value for cls in self.secondary_classes],
            "effect_profile": {
                effect: {"magnitude": mag, "confidence": conf}
                for effect, (mag, conf) in self.effect_profile.items()
            },
        }

    def _format_nootropic_data(self) -> Dict:
        """Format nootropic property data."""
        return {
            "mechanisms": [mech.value for mech in self.nootropic_mechanisms],
            "cognitive_effects": {
                effect: {"magnitude": mag, "confidence": conf}
                for effect, (mag, conf) in self.cognitive_effects.items()
            },
            "side_effects": {
                effect: {"magnitude": mag, "confidence": conf}
                for effect, (mag, conf) in self.side_effects.items()
            },
        }

    def _analyze_effects(self) -> Dict:
        """Analyze reported and predicted effects."""
        effects = {}
        
        # Analyze each effect type
        for effect, mechanisms in self.effect_mechanisms.items():
            # Check binding data for relevant mechanisms
            mechanism_scores = []
            for mechanism in mechanisms:
                receptor, activity = mechanism.split("_")
                if receptor in self.receptor_profiles:
                    affinity, confidence, actual_activity = self.receptor_profiles[receptor]
                    if actual_activity.lower() == activity.lower():
                        score = confidence * (1 / (1 + affinity))
                        mechanism_scores.append(score)
            
            if mechanism_scores:
                # Calculate effect probability
                probability = np.mean(mechanism_scores)
                confidence = len(mechanism_scores) / len(mechanisms)
                
                effects[effect] = {
                    "probability": probability,
                    "confidence": confidence,
                    "mechanisms": [m for m in mechanisms if any(
                        m.startswith(r) for r in self.receptor_profiles
                    )],
                }
                
        return effects

    def _analyze_mechanisms(self) -> Dict:
        """Analyze active mechanisms."""
        mechanisms = {}
        
        # Analyze each nootropic pathway
        for mechanism, pathways in self.nootropic_pathways.items():
            # Check binding data for relevant pathways
            pathway_scores = []
            active_pathways = []
            
            for pathway in pathways:
                receptor, activity = pathway.split("_")
                if receptor in self.receptor_profiles:
                    affinity, confidence, actual_activity = self.receptor_profiles[receptor]
                    if actual_activity.lower() == activity.lower():
                        score = confidence * (1 / (1 + affinity))
                        pathway_scores.append(score)
                        active_pathways.append(pathway)
            
            if pathway_scores:
                # Calculate mechanism strength
                strength = np.mean(pathway_scores)
                confidence = len(pathway_scores) / len(pathways)
                
                mechanisms[mechanism.value] = {
                    "strength": strength,
                    "confidence": confidence,
                    "active_pathways": active_pathways,
                }
                
        return mechanisms

    def _predict_activities(self) -> Dict:
        """Predict likely activities and effects."""
        predictions = {
            "psychoactive": self._predict_psychoactive_class(),
            "nootropic": self._predict_nootropic_effects(),
            "risks": self._predict_activity_risks(),
        }
        
        # Add confidence scores
        predictions["confidence"] = {
            "psychoactive": self._calculate_prediction_confidence(
                predictions["psychoactive"]
            ),
            "nootropic": self._calculate_prediction_confidence(
                predictions["nootropic"]
            ),
            "risks": self._calculate_prediction_confidence(
                predictions["risks"]
            ),
        }
        
        return predictions

    def _predict_psychoactive_class(self) -> Dict:
        """Predict primary psychoactive classification."""
        scores = {}
        
        # Score each class based on receptor binding
        for effect, mechanisms in self.effect_mechanisms.items():
            effect_score = 0
            effect_confidence = 0
            
            for mechanism in mechanisms:
                receptor, activity = mechanism.split("_")
                if receptor in self.receptor_profiles:
                    affinity, confidence, actual_activity = self.receptor_profiles[receptor]
                    if actual_activity.lower() == activity.lower():
                        effect_score += 1 / (1 + affinity)
                        effect_confidence += confidence
            
            if effect_score > 0:
                scores[effect] = {
                    "score": effect_score / len(mechanisms),
                    "confidence": effect_confidence / len(mechanisms),
                }
        
        # Determine primary and secondary classes
        if scores:
            sorted_effects = sorted(
                scores.items(),
                key=lambda x: (x[1]["score"], x[1]["confidence"]),
                reverse=True
            )
            
            return {
                "primary": sorted_effects[0][0],
                "secondary": [e[0] for e in sorted_effects[1:3]],
                "scores": scores,
            }
            
        return {}

    def _predict_nootropic_effects(self) -> Dict:
        """Predict nootropic effects and mechanisms."""
        predictions = {}
        
        # Analyze each nootropic mechanism
        for mechanism, pathways in self.nootropic_pathways.items():
            mechanism_score = 0
            mechanism_confidence = 0
            active_pathways = []
            
            for pathway in pathways:
                receptor, activity = pathway.split("_")
                if receptor in self.receptor_profiles:
                    affinity, confidence, actual_activity = self.receptor_profiles[receptor]
                    if actual_activity.lower() == activity.lower():
                        score = 1 / (1 + affinity)
                        mechanism_score += score
                        mechanism_confidence += confidence
                        active_pathways.append({
                            "pathway": pathway,
                            "score": score,
                            "confidence": confidence,
                        })
            
            if mechanism_score > 0:
                predictions[mechanism.value] = {
                    "score": mechanism_score / len(pathways),
                    "confidence": mechanism_confidence / len(pathways),
                    "active_pathways": active_pathways,
                }
                
        return predictions

    def _predict_activity_risks(self) -> Dict:
        """Predict risks associated with activities."""
        risks = {}
        
        # Check for risky activity patterns
        if self.psychoactive_class != PsychoactiveClass.UNKNOWN:
            # Class-specific risks
            class_risks = {
                PsychoactiveClass.PSYCHEDELIC: {
                    "psychosis": RiskLevel.HIGH,
                    "hppd": RiskLevel.MODERATE,
                },
                PsychoactiveClass.STIMULANT: {
                    "addiction": RiskLevel.HIGH,
                    "cardiovascular": RiskLevel.HIGH,
                },
                PsychoactiveClass.DEPRESSANT: {
                    "dependence": RiskLevel.HIGH,
                    "respiratory": RiskLevel.HIGH,
                },
                PsychoactiveClass.DISSOCIATIVE: {
                    "cognitive": RiskLevel.HIGH,
                    "bladder": RiskLevel.MODERATE,
                },
            }
            
            if self.psychoactive_class in class_risks:
                risks["class_risks"] = {
                    risk: level.value
                    for risk, level in class_risks[self.psychoactive_class].items()
                }
        
        # Mechanism-based risks
        mechanism_risks = []
        for mechanism in self.nootropic_mechanisms:
            if mechanism in [
                NootropicMechanism.GLUTAMATERGIC,
                NootropicMechanism.DOPAMINERGIC,
            ]:
                mechanism_risks.append({
                    "mechanism": mechanism.value,
                    "risk": "excitotoxicity",
                    "level": RiskLevel.MODERATE.value,
                })
                
        if mechanism_risks:
            risks["mechanism_risks"] = mechanism_risks
            
        return risks

    def _calculate_prediction_confidence(self, predictions: Dict) -> float:
        """Calculate overall confidence for predictions."""
        if not predictions:
            return 0.0
            
        confidences = []
        
        def extract_confidence(data: Dict) -> None:
            if isinstance(data, dict):
                if "confidence" in data:
                    confidences.append(data["confidence"])
                for value in data.values():
                    extract_confidence(value)
                    
        extract_confidence(predictions)
        return np.mean(confidences) if confidences else 0.0

    def _merge_effect_data(
        self,
        target_dict: Dict[str, EffectScore],
        source_dict: Dict[str, EffectScore]
    ) -> None:
        """Merge effect data with confidence weighting."""
        for effect, (magnitude, confidence) in source_dict.items():
            if effect not in target_dict:
                target_dict[effect] = (magnitude, confidence)
            else:
                old_mag, old_conf = target_dict[effect]
                if confidence > old_conf:
                    target_dict[effect] = (magnitude, confidence)

    def _merge_psychoactive_data(self, other: "ActivityProfileMixin") -> None:
        """Merge psychoactive classification data."""
        if other.psychoactive_class != PsychoactiveClass.UNKNOWN:
            if self.psychoactive_class == PsychoactiveClass.UNKNOWN:
                self.psychoactive_class = other.psychoactive_class
            else:
                self.secondary_classes.add(self.psychoactive_class)
                self.psychoactive_class = other.psychoactive_class
                
        self.secondary_classes.update(other.secondary_classes)

    def merge_activity_data(self, other: "ActivityProfileMixin") -> None:
        """Merge activity data from another instance."""
        # Merge classification data
        self._merge_psychoactive_data(other)
        
        # Merge all effect data using helper method
        for target, source in [
            (self.effect_profile, other.effect_profile),
            (self.cognitive_effects, other.cognitive_effects),
            (self.side_effects, other.side_effects),
        ]:
            self._merge_effect_data(target, source)
        
        # Merge mechanisms
        self.nootropic_mechanisms.update(other.nootropic_mechanisms)
