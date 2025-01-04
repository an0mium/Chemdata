"""Psychopharmacological property functionality.

This module provides the PsychopharmMixin class that implements methods for
handling psychopharmacological properties including:
1. Receptor binding profiles
2. Blood-brain barrier penetration
3. Psychoactive effects classification
4. Nootropic activity metrics
5. Abuse potential assessment
6. Risk analysis
7. Community data integration
"""

from dataclasses import field
from typing import Dict, List, Optional, Set, Tuple
from enum import Enum
import numpy as np


class PsychoactiveClass(Enum):
    """Classification of psychoactive effects."""
    
    PSYCHEDELIC = "psychedelic"
    EMPATHOGEN = "empathogen"
    STIMULANT = "stimulant"
    DEPRESSANT = "depressant"
    DISSOCIATIVE = "dissociative"
    DELIRIANT = "deliriant"
    NOOTROPIC = "nootropic"
    ANXIOLYTIC = "anxiolytic"
    ANTIPSYCHOTIC = "antipsychotic"
    ANTIDEPRESSANT = "antidepressant"
    MOOD_STABILIZER = "mood_stabilizer"
    UNKNOWN = "unknown"


class NootropicMechanism(Enum):
    """Mechanisms of nootropic activity."""
    
    CHOLINERGIC = "cholinergic"
    GLUTAMATERGIC = "glutamatergic"
    DOPAMINERGIC = "dopaminergic"
    SEROTONERGIC = "serotonergic"
    GABA = "gaba_modulation"
    AMPAKINE = "ampakine"
    BDNF = "bdnf_modulation"
    NGF = "ngf_modulation"
    NEUROPLASTICITY = "neuroplasticity"
    ANTI_INFLAMMATORY = "anti_inflammatory"
    ANTIOXIDANT = "antioxidant"
    UNKNOWN = "unknown"


class BBBPermeability(Enum):
    """Blood-brain barrier permeability classification."""
    
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    NEGLIGIBLE = "negligible"
    UNKNOWN = "unknown"


class RiskLevel(Enum):
    """Risk level classification."""
    
    SEVERE = "severe"
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    MINIMAL = "minimal"
    UNKNOWN = "unknown"


class PsychopharmMixin:
    """Mixin class providing psychopharmacological property methods."""

    # Receptor binding profiles
    receptor_profiles: Dict[str, Dict] = field(default_factory=dict)
    
    # Blood-brain barrier properties
    bbb_permeability: BBBPermeability = BBBPermeability.UNKNOWN
    bbb_score: float = 0.0
    p_glycoprotein_substrate: bool = False
    
    # Psychoactive classification
    psychoactive_class: PsychoactiveClass = PsychoactiveClass.UNKNOWN
    secondary_classes: Set[PsychoactiveClass] = field(default_factory=set)
    effect_profile: Dict[str, float] = field(default_factory=dict)
    
    # Nootropic properties
    nootropic_mechanisms: Set[NootropicMechanism] = field(default_factory=set)
    cognitive_effects: Dict[str, Dict] = field(default_factory=dict)
    side_effects: Dict[str, Dict] = field(default_factory=dict)
    
    # Duration metrics
    onset_time: Optional[str] = None
    duration: Optional[str] = None
    half_life: Optional[str] = None
    
    # Tolerance and withdrawal
    tolerance_profile: Dict[str, str] = field(default_factory=dict)
    withdrawal_profile: Dict[str, str] = field(default_factory=dict)
    cross_tolerance: Set[str] = field(default_factory=set)

    # Risk assessment
    risk_factors: Dict[str, RiskLevel] = field(default_factory=dict)
    contraindications: List[Dict] = field(default_factory=list)
    interaction_risks: Dict[str, Dict] = field(default_factory=dict)
    overdose_risk: RiskLevel = RiskLevel.UNKNOWN
    addiction_risk: RiskLevel = RiskLevel.UNKNOWN

    # Community data
    experience_reports: List[Dict] = field(default_factory=list)
    reported_effects: Dict[str, int] = field(default_factory=dict)
    reported_interactions: Dict[str, List[Dict]] = field(default_factory=dict)
    reported_risks: Dict[str, List[Dict]] = field(default_factory=dict)

    def get_psychopharm_dict(self) -> Dict:
        """Get dictionary of psychopharmacological properties."""
        return {
            "receptor_profiles": self._format_receptor_profiles(),
            "bbb_properties": self._format_bbb_properties(),
            "psychoactive_properties": self._format_psychoactive_properties(),
            "nootropic_properties": self._format_nootropic_properties(),
            "duration_metrics": self._format_duration_metrics(),
            "tolerance_data": self._format_tolerance_data(),
            "risk_assessment": self._format_risk_assessment(),
            "community_data": self._format_community_data(),
        }

    def _format_receptor_profiles(self) -> Dict:
        """Format receptor binding profiles."""
        formatted = {}
        for receptor, data in self.receptor_profiles.items():
            formatted[receptor] = {
                "affinity": data.get("affinity", 0),
                "activity": data.get("activity", "unknown"),
                "confidence": data.get("confidence", 0),
                "source": data.get("source", "predicted"),
                "notes": data.get("notes", ""),
                "binding_pattern": self._analyze_binding_pattern(receptor),
            }
        return formatted

    def _format_bbb_properties(self) -> Dict:
        """Format blood-brain barrier properties."""
        return {
            "permeability_class": self.bbb_permeability.value,
            "permeability_score": self.bbb_score,
            "p_glycoprotein_substrate": self.p_glycoprotein_substrate,
            "predicted_concentration": self._predict_brain_concentration(),
        }

    def _format_psychoactive_properties(self) -> Dict:
        """Format psychoactive classification data."""
        return {
            "primary_class": self.psychoactive_class.value,
            "secondary_classes": [cls.value for cls in self.secondary_classes],
            "effect_profile": self.effect_profile,
            "predicted_effects": self._predict_psychoactive_effects(),
        }

    def _format_nootropic_properties(self) -> Dict:
        """Format nootropic property data."""
        return {
            "mechanisms": [mech.value for mech in self.nootropic_mechanisms],
            "cognitive_effects": self.cognitive_effects,
            "side_effects": self.side_effects,
            "efficacy_analysis": self._analyze_nootropic_efficacy(),
        }

    def _format_duration_metrics(self) -> Dict:
        """Format duration metrics."""
        return {
            "onset": self.onset_time,
            "duration": self.duration,
            "half_life": self.half_life,
            "predicted_timeline": self._predict_effect_timeline(),
        }

    def _format_tolerance_data(self) -> Dict:
        """Format tolerance and withdrawal data."""
        return {
            "tolerance_profile": self.tolerance_profile,
            "withdrawal_profile": self.withdrawal_profile,
            "cross_tolerance": list(self.cross_tolerance),
            "risk_analysis": self._analyze_tolerance_risks(),
        }

    def _format_risk_assessment(self) -> Dict:
        """Format risk assessment data."""
        return {
            "risk_factors": {k: v.value for k, v in self.risk_factors.items()},
            "contraindications": self.contraindications,
            "interaction_risks": self.interaction_risks,
            "overdose_risk": self.overdose_risk.value,
            "addiction_risk": self.addiction_risk.value,
            "risk_analysis": self._analyze_overall_risks(),
        }

    def _format_community_data(self) -> Dict:
        """Format community data."""
        return {
            "experience_reports": self._summarize_experience_reports(),
            "reported_effects": self.reported_effects,
            "reported_interactions": self.reported_interactions,
            "reported_risks": self.reported_risks,
            "community_analysis": self._analyze_community_data(),
        }

    def merge_psychopharm_data(self, other: "PsychopharmMixin") -> None:
        """Merge psychopharmacological data from another instance."""
        self._merge_receptor_profiles(other)
        self._merge_bbb_properties(other)
        self._merge_psychoactive_properties(other)
        self._merge_nootropic_properties(other)
        self._merge_duration_metrics(other)
        self._merge_tolerance_data(other)
        self._merge_risk_data(other)
        self._merge_community_data(other)

    def _merge_receptor_profiles(self, other: "PsychopharmMixin") -> None:
        """Merge receptor binding profiles."""
        for receptor, data in other.receptor_profiles.items():
            if receptor not in self.receptor_profiles:
                self.receptor_profiles[receptor] = data
            else:
                # Keep data with higher confidence
                old_conf = self.receptor_profiles[receptor].get("confidence", 0)
                if data.get("confidence", 0) > old_conf:
                    self.receptor_profiles[receptor] = data

    def _merge_bbb_properties(self, other: "PsychopharmMixin") -> None:
        """Merge blood-brain barrier properties."""
        if other.bbb_score > self.bbb_score:
            self.bbb_permeability = other.bbb_permeability
            self.bbb_score = other.bbb_score
            self.p_glycoprotein_substrate = other.p_glycoprotein_substrate

    def _merge_psychoactive_properties(self, other: "PsychopharmMixin") -> None:
        """Merge psychoactive classification data."""
        if other.psychoactive_class != PsychoactiveClass.UNKNOWN:
            if self.psychoactive_class == PsychoactiveClass.UNKNOWN:
                self.psychoactive_class = other.psychoactive_class
            else:
                self.secondary_classes.add(self.psychoactive_class)
                self.psychoactive_class = other.psychoactive_class
        self.secondary_classes.update(other.secondary_classes)
        self.effect_profile.update(other.effect_profile)

    def _merge_nootropic_properties(self, other: "PsychopharmMixin") -> None:
        """Merge nootropic properties."""
        self.nootropic_mechanisms.update(other.nootropic_mechanisms)
        self.cognitive_effects.update(other.cognitive_effects)
        self.side_effects.update(other.side_effects)

    def _merge_duration_metrics(self, other: "PsychopharmMixin") -> None:
        """Merge duration metrics."""
        if not self.onset_time:
            self.onset_time = other.onset_time
        if not self.duration:
            self.duration = other.duration
        if not self.half_life:
            self.half_life = other.half_life

    def _merge_tolerance_data(self, other: "PsychopharmMixin") -> None:
        """Merge tolerance and withdrawal data."""
        self.tolerance_profile.update(other.tolerance_profile)
        self.withdrawal_profile.update(other.withdrawal_profile)
        self.cross_tolerance.update(other.cross_tolerance)

    def _merge_risk_data(self, other: "PsychopharmMixin") -> None:
        """Merge risk assessment data."""
        self.risk_factors.update(other.risk_factors)
        self.contraindications.extend(
            c for c in other.contraindications
            if c not in self.contraindications
        )
        self.interaction_risks.update(other.interaction_risks)
        if other.overdose_risk.value > self.overdose_risk.value:
            self.overdose_risk = other.overdose_risk
        if other.addiction_risk.value > self.addiction_risk.value:
            self.addiction_risk = other.addiction_risk

    def _merge_community_data(self, other: "PsychopharmMixin") -> None:
        """Merge community data."""
        self.experience_reports.extend(
            report for report in other.experience_reports
            if report not in self.experience_reports
        )
        for effect, count in other.reported_effects.items():
            self.reported_effects[effect] = (
                self.reported_effects.get(effect, 0) + count
            )
        self.reported_interactions.update(other.reported_interactions)
        self.reported_risks.update(other.reported_risks)

    def _analyze_binding_pattern(self, receptor: str) -> Dict:
        """Analyze binding pattern for receptor."""
        data = self.receptor_profiles[receptor]
        return {
            "binding_type": self._predict_binding_type(data),
            "selectivity": self._calculate_selectivity(receptor),
            "predicted_effects": self._predict_receptor_effects(receptor),
        }

    def _predict_binding_type(self, data: Dict) -> str:
        """Predict binding type from receptor data."""
        affinity = data.get("affinity", 0)
        if affinity <= 1:  # Sub-nanomolar
            return "very_strong"
        elif affinity <= 10:  # Low nanomolar
            return "strong"
        elif affinity <= 100:  # High nanomolar
            return "moderate"
        elif affinity <= 1000:  # Micromolar
            return "weak"
        return "negligible"

    def _calculate_selectivity(self, receptor: str) -> float:
        """Calculate receptor selectivity ratio."""
        target_affinity = self.receptor_profiles[receptor].get("affinity", 0)
        if not target_affinity:
            return 0.0
            
        other_affinities = [
            data.get("affinity", 0)
            for r, data in self.receptor_profiles.items()
            if r != receptor and data.get("affinity", 0) > 0
        ]
        
        if not other_affinities:
            return 1.0
            
        return target_affinity / np.mean(other_affinities)

    def _predict_receptor_effects(self, receptor: str) -> List[str]:
        """Predict effects based on receptor binding."""
        data = self.receptor_profiles[receptor]
        effects = []
        
        # Example predictions based on receptor and activity
        if "5-HT2A" in receptor:
            if data.get("activity") == "agonist":
                effects.extend(["psychedelic", "cognitive_enhancement"])
            elif data.get("activity") == "antagonist":
                effects.extend(["antipsychotic", "neuroprotective"])
                
        elif "D2" in receptor:
            if data.get("activity") == "agonist":
                effects.extend(["stimulant", "euphoric"])
            elif data.get("activity") == "antagonist":
                effects.extend(["antipsychotic", "sedating"])
                
        return effects

    def _predict_brain_concentration(self) -> Dict:
        """Predict brain concentration based on properties."""
        if self.bbb_permeability == BBBPermeability.HIGH:
            factor = 0.8
        elif self.bbb_permeability == BBBPermeability.MODERATE:
            factor = 0.5
        elif self.bbb_permeability == BBBPermeability.LOW:
            factor = 0.2
        else:
            factor = 0.0
            
        if self.p_glycoprotein_substrate:
            factor *= 0.5
            
        return {
            "relative_concentration": factor,
            "confidence": self.bbb_score,
        }

    def _predict_psychoactive_effects(self) -> Dict:
        """Predict psychoactive effects based on receptor profile."""
        effects = {}
        
        # Analyze receptor contributions
        for receptor, data in self.receptor_profiles.items():
            predicted = self._predict_receptor_effects(receptor)
            for effect in predicted:
                effects[effect] = effects.get(effect, 0) + (
                    data.get("confidence", 0) * 
                    (1 / (1 + data.get("affinity", 1000)))
                )
                
        # Normalize scores
        if effects:
            max_score = max(effects.values())
            effects = {
                k: v/max_score
                for k, v in effects.items()
            }
            
        return effects

    def _analyze_nootropic_efficacy(self) -> Dict:
        """Analyze nootropic efficacy based on mechanisms."""
        efficacy = {}
        
        for mechanism in self.nootropic_mechanisms:
            # Calculate efficacy score based on receptor profiles
            score = self._calculate_mechanism_efficacy(mechanism)
            
            # Get supporting evidence from cognitive effects
            evidence = [
                effect
                for effect in self.cognitive_effects
                if effect in self._get_mechanism_effects(mechanism)
            ]
            
            efficacy[mechanism.value] = {
                "score": score,
                "supporting_evidence": evidence,
                "confidence": len(evidence) / 10,  # Scale 0-1
            }
            
        return efficacy

    def _calculate_mechanism_efficacy(self, mechanism: NootropicMechanism) -> float:
        """Calculate efficacy score for nootropic mechanism."""
        relevant_receptors = self._get_mechanism_receptors(mechanism)
        
        scores = []
        for receptor in relevant_receptors:
            if receptor in self.receptor_profiles:
                data = self.receptor_profiles[receptor]
                score = (
                    data.get("confidence", 0) * 
                    (1 / (1 + data.get("affinity", 1000)))
                )
                scores.append(score)
                
        return np.mean(scores) if scores else 0.0

    def _get_mechanism_receptors(self, mechanism: NootropicMechanism) -> List[str]:
        """Get relevant receptors for mechanism."""
        # Example mapping
        mapping = {
            NootropicMechanism.CHOLINERGIC: ["nAChR", "mAChR"],
            NootropicMechanism.GLUTAMATERGIC: ["NMDA", "AMPA", "mGluR"],
            NootropicMechanism.DOPAMINERGIC: ["D1", "D2", "D3", "D4", "D5"],
            NootropicMechanism.SEROTONERGIC: [
                "5-HT1A", "5-HT2A", "5-HT2C", "5-HT6", "5-HT7"
            ],
        }
        return mapping.get(mechanism, [])

    def _get_mechanism_effects(self, mechanism: NootropicMechanism) -> List[str]:
        """Get cognitive effects associated with mechanism."""
        # Example mapping
        mapping = {
            NootropicMechanism.CHOLINERGIC: [
                "memory_enhancement",
                "attention_enhancement",
            ],
            NootropicMechanism.GLUTAMATERGIC: [
                "learning_enhancement",
                "neuroplasticity",
            ],
            NootropicMechanism.DOPAMINERGIC: [
                "motivation_enhancement",
                "focus_enhancement",
            ],
            NootropicMechanism.SEROTONERGIC: [
                "mood_enhancement",
                "anxiety_reduction",
            ],
        }
        return mapping.get(mechanism, [])

    def _predict_effect_timeline(self) -> Dict:
        """Predict detailed effect timeline."""
        if not self.onset_time or not self.duration or not self.half_life:
            return {}
            
        try:
            # Parse duration values
            onset = self._parse_duration(self.onset_time)
            duration = self._parse_duration(self.duration)
            half_life = self._parse_duration(self.half_life)
            
            # Generate timeline
            timeline = {
                "onset": {
                    "start": 0,
                    "peak": onset,
                },
                "peak": {
                    "start": onset,
                    "end": onset + duration,
                },
                "offset": {
                    "start": onset + duration,
                    "end": onset + duration + (half_life * 5),
                },
            }
            
            return timeline
            
        except ValueError:
            return {}

    def _parse_duration(self, duration_str: str) -> float:
        """Parse duration string to hours."""
        # Example: "2-4 hours" -> 3.0
        try:
            parts = duration_str.split()
            if len(parts) != 2:
                return 0.0
                
            time_range = parts[0].split("-")
            unit = parts[1].lower()
            
            if len(time_range) == 2:
                time = (float(time_range[0]) + float(time_range[1])) / 2
            else:
                time = float(time_range[0])
                
            # Convert to hours
            if unit.startswith("minute"):
                return time / 60
            elif unit.startswith("hour"):
                return time
            elif unit.startswith("day"):
                return time * 24
                
            return 0.0
            
        except (ValueError, IndexError):
            return 0.0

    def _analyze_tolerance_risks(self) -> Dict:
        """Analyze tolerance development risks."""
        risk_factors = []
        
        # Check receptor profiles
        for receptor, data in self.receptor_profiles.items():
            if self._has_tolerance_risk(receptor, data):
                risk_factors.append({
                    "receptor": receptor,
                    "mechanism": "receptor_downregulation",
                    "confidence": data.get("confidence", 0),
                })
                
        # Check reported tolerance
        tolerance_speed = self._analyze_tolerance_speed()
        
        return {
            "risk_factors": risk_factors,
            "tolerance_speed": tolerance_speed,
            "cross_tolerance_risks": self._analyze_cross_tolerance(),
            "withdrawal_risk": self._analyze_withdrawal_risk(),
        }

    def _has_tolerance_risk(self, receptor: str, data: Dict) -> bool:
        """Check if receptor binding suggests tolerance risk."""
        high_risk_receptors = {
            "5-HT2A", "D2", "mu-opioid", "GABA-A", "CB1"
        }
        
        if any(r in receptor for r in high_risk_receptors):
            affinity = data.get("affinity", 0)
            activity = data.get("activity", "")
            
            # Strong agonists have higher tolerance risk
            return (
                affinity < 100 and  # Strong binding
                activity.lower() == "agonist"
            )
            
        return False

    def _analyze_tolerance_speed(self) -> str:
        """Analyze speed of tolerance development."""
        if not self.tolerance_profile:
            return "unknown"
            
        rapid_indicators = [
            "rapid", "fast", "quick", "immediate"
        ]
        
        profile_text = " ".join(self.tolerance_profile.values()).lower()
        
        if any(i in profile_text for i in rapid_indicators):
            return "rapid"
        return "gradual"

    def _analyze_cross_tolerance(self) -> List[Dict]:
        """Analyze cross-tolerance risks."""
        risks = []
        
        for compound in self.cross_tolerance:
            # Example analysis
            risks.append({
                "compound": compound,
                "mechanism": "shared_receptor_adaptation",
                "risk_level": "high",
            })
            
        return risks

    def _analyze_withdrawal_risk(self) -> Dict:
        """Analyze withdrawal risks."""
        risk_level = RiskLevel.UNKNOWN
        risk_factors = []
        
        # Analyze receptor profiles
        for receptor, data in self.receptor_profiles.items():
            if self._has_withdrawal_risk(receptor, data):
                risk_factors.append({
                    "receptor": receptor,
                    "mechanism": "receptor_adaptation",
                    "severity": "high",
                })
                
        # Consider tolerance profile
        if self.tolerance_profile:
            profile_risk = self._assess_tolerance_withdrawal_risk()
            risk_factors.extend(profile_risk)
            
        # Set overall risk level
        if risk_factors:
            risk_count = len(risk_factors)
            if risk_count > 3:
                risk_level = RiskLevel.SEVERE
            elif risk_count > 1:
                risk_level = RiskLevel.HIGH
            else:
                risk_level = RiskLevel.MODERATE
                
        return {
            "risk_level": risk_level.value,
            "risk_factors": risk_factors,
            "withdrawal_profile": self.withdrawal_profile,
        }

    def _has_withdrawal_risk(self, receptor: str, data: Dict) -> bool:
        """Check if receptor binding suggests withdrawal risk."""
        high_risk_receptors = {
            "GABA-A", "mu-opioid", "D2", "5-HT2A"
        }
        
        if any(r in receptor for r in high_risk_receptors):
            affinity = data.get("affinity", 0)
            activity = data.get("activity", "")
            
            # Strong agonists of certain receptors have higher withdrawal risk
            return (
                affinity < 100 and  # Strong binding
                activity.lower() == "agonist"
            )
            
        return False

    def _assess_tolerance_withdrawal_risk(self) -> List[Dict]:
        """Assess withdrawal risk based on tolerance profile."""
        risk_factors = []
        
        for aspect, description in self.tolerance_profile.items():
            if "rapid" in description.lower():
                risk_factors.append({
                    "factor": f"rapid_{aspect}_tolerance",
                    "mechanism": "rapid_adaptation",
                    "severity": "high",
                })
            elif "significant" in description.lower():
                risk_factors.append({
                    "factor": f"significant_{aspect}_tolerance",
                    "mechanism": "substantial_adaptation",
                    "severity": "moderate",
                })
                
        return risk_factors

    def _analyze_overall_risks(self) -> Dict:
        """Analyze overall risk profile."""
        return {
            "acute_risks": self._analyze_acute_risks(),
            "chronic_risks": self._analyze_chronic_risks(),
            "interaction_risks": self._analyze_interaction_risks(),
            "population_risks": self._analyze_population_risks(),
        }

    def _analyze_acute_risks(self) -> List[Dict]:
        """Analyze acute risk factors."""
        risks = []
        
        # Check receptor profiles for dangerous combinations
        for receptor, data in self.receptor_profiles.items():
            if self._has_acute_risk(receptor, data):
                risks.append({
                    "type": "receptor_mediated",
                    "receptor": receptor,
                    "mechanism": self._get_risk_mechanism(receptor),
                    "severity": RiskLevel.HIGH.value,
                })
                
        # Add overdose risks
        if self.overdose_risk != RiskLevel.UNKNOWN:
            risks.append({
                "type": "overdose",
                "mechanism": "dose_dependent_toxicity",
                "severity": self.overdose_risk.value,
            })
            
        return risks

    def _has_acute_risk(self, receptor: str, data: Dict) -> bool:
        """Check if receptor binding presents acute risks."""
        high_risk_receptors = {
            "5-HT2B": "cardiotoxicity",
            "mu-opioid": "respiratory_depression",
            "NMDA": "excitotoxicity",
        }
        
        if any(r in receptor for r in high_risk_receptors):
            affinity = data.get("affinity", 0)
            activity = data.get("activity", "")
            
            # Strong agonists/antagonists of certain receptors are risky
            return (
                affinity < 100 and  # Strong binding
                activity.lower() in ["agonist", "antagonist"]
            )
            
        return False

    def _get_risk_mechanism(self, receptor: str) -> str:
        """Get risk mechanism for receptor."""
        mechanisms = {
            "5-HT2B": "cardiac_valve_proliferation",
            "mu-opioid": "respiratory_depression",
            "NMDA": "excitotoxicity",
            "5-HT2A": "serotonin_syndrome",
            "D2": "dopamine_dysregulation",
        }
        
        for r, mechanism in mechanisms.items():
            if r in receptor:
                return mechanism
                
        return "unknown"

    def _analyze_chronic_risks(self) -> List[Dict]:
        """Analyze chronic risk factors."""
        risks = []
        
        # Analyze tolerance risks
        tolerance_risks = self._analyze_tolerance_risks()
        if tolerance_risks.get("risk_factors"):
            risks.extend([
                {
                    "type": "tolerance",
                    "mechanism": factor["mechanism"],
                    "severity": RiskLevel.MODERATE.value,
                }
                for factor in tolerance_risks["risk_factors"]
            ])
            
        # Analyze withdrawal risks
        withdrawal_analysis = self._analyze_withdrawal_risk()
        if withdrawal_analysis["risk_level"] != RiskLevel.UNKNOWN.value:
            risks.append({
                "type": "withdrawal",
                "mechanism": "physiological_dependence",
                "severity": withdrawal_analysis["risk_level"],
            })
            
        # Check addiction risk
        if self.addiction_risk != RiskLevel.UNKNOWN:
            risks.append({
                "type": "addiction",
                "mechanism": "reward_pathway_adaptation",
                "severity": self.addiction_risk.value,
            })
            
        return risks

    def _analyze_interaction_risks(self) -> List[Dict]:
        """Analyze drug interaction risks."""
        risks = []
        
        # Analyze reported interactions
        for drug, interactions in self.reported_interactions.items():
            for interaction in interactions:
                risks.append({
                    "interacting_drug": drug,
                    "mechanism": interaction.get("mechanism", "unknown"),
                    "severity": interaction.get("severity", RiskLevel.UNKNOWN.value),
                    "evidence": interaction.get("evidence", []),
                })
                
        # Predict additional interactions based on receptor profile
        predicted = self._predict_interactions()
        risks.extend(predicted)
        
        return risks

    def _predict_interactions(self) -> List[Dict]:
        """Predict potential drug interactions."""
        predictions = []
        
        # Example interaction predictions based on receptor profile
        for receptor, data in self.receptor_profiles.items():
            # Serotonergic interactions
            if "5-HT" in receptor:
                predictions.append({
                    "interacting_drug": "SSRIs",
                    "mechanism": "serotonin_syndrome",
                    "severity": RiskLevel.SEVERE.value,
                    "confidence": data.get("confidence", 0),
                })
                predictions.append({
                    "interacting_drug": "MAOIs",
                    "mechanism": "serotonin_syndrome",
                    "severity": RiskLevel.SEVERE.value,
                    "confidence": data.get("confidence", 0),
                })
                
            # Dopaminergic interactions    
            elif "D2" in receptor:
                predictions.append({
                    "interacting_drug": "Antipsychotics",
                    "mechanism": "dopamine_modulation",
                    "severity": RiskLevel.HIGH.value,
                    "confidence": data.get("confidence", 0),
                })
                predictions.append({
                    "interacting_drug": "Stimulants",
                    "mechanism": "dopamine_potentiation",
                    "severity": RiskLevel.HIGH.value,
                    "confidence": data.get("confidence", 0),
                })
                
            # GABAergic interactions
            elif "GABA" in receptor:
                predictions.append({
                    "interacting_drug": "Benzodiazepines",
                    "mechanism": "gaba_potentiation",
                    "severity": RiskLevel.HIGH.value,
                    "confidence": data.get("confidence", 0),
                })
                predictions.append({
                    "interacting_drug": "Alcohol",
                    "mechanism": "respiratory_depression",
                    "severity": RiskLevel.SEVERE.value,
                    "confidence": data.get("confidence", 0),
                })
                
            # NMDA interactions
            elif "NMDA" in receptor:
                predictions.append({
                    "interacting_drug": "Dissociatives",
                    "mechanism": "nmda_antagonism",
                    "severity": RiskLevel.HIGH.value,
                    "confidence": data.get("confidence", 0),
                })
                predictions.append({
                    "interacting_drug": "Alcohol",
                    "mechanism": "cognitive_impairment",
                    "severity": RiskLevel.HIGH.value,
                    "confidence": data.get("confidence", 0),
                })
                
            # Opioid interactions
            elif "mu-opioid" in receptor:
                predictions.append({
                    "interacting_drug": "Opioids",
                    "mechanism": "respiratory_depression",
                    "severity": RiskLevel.SEVERE.value,
                    "confidence": data.get("confidence", 0),
                })
                predictions.append({
                    "interacting_drug": "Benzodiazepines",
                    "mechanism": "respiratory_depression",
                    "severity": RiskLevel.SEVERE.value,
                    "confidence": data.get("confidence", 0),
                })
                
        return predictions

    def _analyze_population_risks(self) -> List[Dict]:
        """Analyze population-specific risk factors."""
        risks = []
        
        # Psychiatric risks
        if any("5-HT2A" in r for r in self.receptor_profiles):
            risks.append({
                "population": "psychiatric_history",
                "condition": "psychosis",
                "mechanism": "serotonin_modulation",
                "severity": RiskLevel.SEVERE.value,
            })
            risks.append({
                "population": "bipolar_disorder",
                "condition": "mania",
                "mechanism": "serotonin_modulation",
                "severity": RiskLevel.HIGH.value,
            })
            
        # Cardiovascular risks    
        if any("D2" in r for r in self.receptor_profiles):
            risks.append({
                "population": "cardiovascular_disease",
                "condition": "arrhythmia",
                "mechanism": "autonomic_effects",
                "severity": RiskLevel.HIGH.value,
            })
            risks.append({
                "population": "hypertension",
                "condition": "blood_pressure",
                "mechanism": "sympathetic_activation",
                "severity": RiskLevel.HIGH.value,
            })
            
        # Respiratory risks
        if any("GABA" in r for r in self.receptor_profiles):
            risks.append({
                "population": "respiratory_conditions",
                "condition": "sleep_apnea",
                "mechanism": "respiratory_depression",
                "severity": RiskLevel.HIGH.value,
            })
            risks.append({
                "population": "asthma",
                "condition": "breathing_difficulty",
                "mechanism": "respiratory_depression",
                "severity": RiskLevel.HIGH.value,
            })
            
        # Neurological risks
        if any("NMDA" in r for r in self.receptor_profiles):
            risks.append({
                "population": "epilepsy",
                "condition": "seizures",
                "mechanism": "glutamate_modulation",
                "severity": RiskLevel.HIGH.value,
            })
            risks.append({
                "population": "traumatic_brain_injury",
                "condition": "cognitive_impairment",
                "mechanism": "glutamate_dysregulation",
                "severity": RiskLevel.HIGH.value,
            })
            
        # Liver risks
        if any("CYP" in r for r in self.receptor_profiles):
            risks.append({
                "population": "liver_disease",
                "condition": "metabolism_impairment",
                "mechanism": "enzyme_inhibition",
                "severity": RiskLevel.HIGH.value,
            })
            risks.append({
                "population": "hepatitis",
                "condition": "liver_toxicity",
                "mechanism": "metabolic_stress",
                "severity": RiskLevel.HIGH.value,
            })
            
        return risks

    def _summarize_experience_reports(self) -> Dict:
        """Summarize and analyze experience reports."""
        if not self.experience_reports:
            return {}
            
        summary = {
            "total_reports": len(self.experience_reports),
            "effect_frequencies": {},
            "dose_ranges": {},
            "common_combinations": {},
            "reported_issues": {},
        }
        
        for report in self.experience_reports:
            # Analyze effects
            for effect in report.get("effects", []):
                summary["effect_frequencies"][effect] = (
                    summary["effect_frequencies"].get(effect, 0) + 1
                )
                
            # Analyze doses
            dose = report.get("dose", {})
            if dose:
                route = dose.get("route", "unknown")
                amount = dose.get("amount")
                if amount:
                    if route not in summary["dose_ranges"]:
                        summary["dose_ranges"][route] = {
                            "min": amount,
                            "max": amount,
                            "counts": {},
                        }
                    else:
                        summary["dose_ranges"][route]["min"] = min(
                            summary["dose_ranges"][route]["min"],
                            amount
                        )
                        summary["dose_ranges"][route]["max"] = max(
                            summary["dose_ranges"][route]["max"],
                            amount
                        )
                    summary["dose_ranges"][route]["counts"][amount] = (
                        summary["dose_ranges"][route]["counts"].get(amount, 0) + 1
                    )
                    
            # Analyze combinations
            for combo in report.get("combinations", []):
                summary["common_combinations"][combo] = (
                    summary["common_combinations"].get(combo, 0) + 1
                )
                
            # Analyze issues
            for issue in report.get("issues", []):
                summary["reported_issues"][issue] = (
                    summary["reported_issues"].get(issue, 0) + 1
                )
                
        return summary

    def _analyze_community_data(self) -> Dict:
        """Analyze community-reported data."""
        return {
            "effect_analysis": self._analyze_reported_effects(),
            "safety_analysis": self._analyze_reported_safety(),
            "usage_patterns": self._analyze_usage_patterns(),
            "risk_factors": self._analyze_reported_risks(),
        }

    def _analyze_reported_effects(self) -> Dict:
        """Analyze reported effects."""
        if not self.reported_effects:
            return {}
            
        total_reports = sum(self.reported_effects.values())
        
        return {
            "total_reports": total_reports,
            "effect_frequencies": {
                effect: {
                    "count": count,
                    "frequency": count / total_reports,
                }
                for effect, count in self.reported_effects.items()
            },
            "primary_effects": [
                effect for effect, count in self.reported_effects.items()
                if count / total_reports > 0.5
            ],
            "rare_effects": [
                effect for effect, count in self.reported_effects.items()
                if count / total_reports < 0.1
            ],
        }

    def _analyze_reported_safety(self) -> Dict:
        """Analyze reported safety information."""
        if not self.reported_risks:
            return {}
            
        severity_counts = {
            "low": 0,
            "moderate": 0,
            "high": 0,
            "severe": 0,
        }
        
        risk_categories = {}
        
        for category, risks in self.reported_risks.items():
            risk_categories[category] = {
                "count": len(risks),
                "severity_distribution": {},
                "common_factors": [],
            }
            
            # Analyze severity distribution
            for risk in risks:
                severity = risk.get("severity", "unknown").lower()
                if severity in severity_counts:
                    severity_counts[severity] += 1
                    risk_categories[category]["severity_distribution"][severity] = (
                        risk_categories[category]["severity_distribution"].get(severity, 0) + 1
                    )
                    
            # Identify common factors
            factors = [
                risk.get("factor")
                for risk in risks
                if risk.get("factor")
            ]
            if factors:
                from collections import Counter
                common = Counter(factors).most_common(3)
                risk_categories[category]["common_factors"] = [
                    {"factor": f, "count": c}
                    for f, c in common
                ]
                
        return {
            "severity_distribution": severity_counts,
            "risk_categories": risk_categories,
            "high_priority_risks": [
                risk
                for risks in self.reported_risks.values()
                for risk in risks
                if risk.get("severity", "").lower() in ["high", "severe"]
            ],
        }

    def _analyze_usage_patterns(self) -> Dict:
        """Analyze reported usage patterns."""
        if not self.experience_reports:
            return {}
            
        patterns = {
            "routes": {},
            "doses": {},
            "frequencies": {},
            "contexts": {},
        }
        
        for report in self.experience_reports:
            # Analyze administration routes
            route = report.get("route")
            if route:
                patterns["routes"][route] = patterns["routes"].get(route, 0) + 1
                
            # Analyze dosing
            dose = report.get("dose")
            if dose:
                amount = dose.get("amount")
                if amount:
                    if route not in patterns["doses"]:
                        patterns["doses"][route] = []
                    patterns["doses"][route].append(amount)
                    
            # Analyze frequency
            frequency = report.get("frequency")
            if frequency:
                patterns["frequencies"][frequency] = (
                    patterns["frequencies"].get(frequency, 0) + 1
                )
                
            # Analyze context
            context = report.get("context")
            if context:
                patterns["contexts"][context] = (
                    patterns["contexts"].get(context, 0) + 1
                )
                
        # Calculate statistics for doses
        for route in patterns["doses"]:
            doses = patterns["doses"][route]
            patterns["doses"][route] = {
                "min": min(doses),
                "max": max(doses),
                "mean": sum(doses) / len(doses),
                "count": len(doses),
            }
            
        return patterns

    def _analyze_reported_risks(self) -> Dict:
        """Analyze reported risk factors."""
        if not self.reported_risks:
            return {}
            
        analysis = {
            "total_risks": sum(len(risks) for risks in self.reported_risks.values()),
            "risk_categories": {},
            "severity_distribution": {},
            "common_factors": [],
            "high_priority": [],
        }
        
        all_factors = []
        
        for category, risks in self.reported_risks.items():
            category_analysis = {
                "count": len(risks),
                "severity_distribution": {},
                "factors": [],
            }
            
            for risk in risks:
                # Analyze severity
                severity = risk.get("severity", "unknown")
                category_analysis["severity_distribution"][severity] = (
                    category_analysis["severity_distribution"].get(severity, 0) + 1
                )
                analysis["severity_distribution"][severity] = (
                    analysis["severity_distribution"].get(severity, 0) + 1
                )
                
                # Track factors
                factor = risk.get("factor")
                if factor:
                    category_analysis["factors"].append(factor)
                    all_factors.append(factor)
                    
                # Track high priority risks
                if severity.lower() in ["high", "severe"]:
                    analysis["high_priority"].append(risk)
                    
            analysis["risk_categories"][category] = category_analysis
            
        # Identify common factors
        if all_factors:
            from collections import Counter
            analysis["common_factors"] = [
                {"factor": f, "count": c}
                for f, c in Counter(all_factors).most_common(5)
            ]
            
        return analysis
