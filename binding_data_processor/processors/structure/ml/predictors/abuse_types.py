"""Abuse potential categories and risk factors.

This module defines:
1. Abuse potential categories and mechanisms
2. Addiction risk factors
3. Withdrawal effects and severity levels
4. Receptor systems involved in abuse
5. Risk level definitions
"""

from typing import Dict, List

# Abuse potential categories with detailed mechanisms
ABUSE_CATEGORIES = {
    "reward_potential": {
        "mechanisms": [
            "dopamine_release",
            "serotonin_release",
            "opioid_activation",
            "gaba_modulation",
            "glutamate_modulation",
        ],
        "risk_levels": ["low", "moderate", "high", "very_high"],
        "description": "Potential for producing rewarding/euphoric effects",
    },
    "addiction_risk": {
        "mechanisms": [
            "reward_sensitization",
            "craving_induction",
            "compulsive_use",
            "loss_of_control",
        ],
        "risk_levels": ["low", "moderate", "high", "very_high"],
        "description": "Risk of developing addiction/dependence",
    },
    "dependence_liability": {
        "mechanisms": [
            "physical_dependence",
            "psychological_dependence",
            "withdrawal_syndrome",
            "tolerance_development",
        ],
        "risk_levels": ["low", "moderate", "high", "very_high"],
        "description": "Likelihood of developing physical/psychological dependence",
    },
    "withdrawal_risk": {
        "mechanisms": [
            "receptor_adaptation",
            "neurotransmitter_depletion",
            "homeostatic_disruption",
            "stress_response",
        ],
        "risk_levels": ["mild", "moderate", "severe", "life_threatening"],
        "description": "Severity of withdrawal symptoms upon discontinuation",
    },
    "tolerance_risk": {
        "mechanisms": [
            "receptor_downregulation",
            "metabolic_adaptation",
            "synaptic_plasticity",
            "behavioral_adaptation",
        ],
        "risk_levels": ["low", "moderate", "high", "very_high"],
        "description": "Risk of developing tolerance to effects",
    },
}

# Detailed abuse effects
ABUSE_EFFECTS = {
    "euphoria": ["mood_elevation", "pleasure", "wellbeing"],
    "stimulation": ["arousal", "energy", "motivation"],
    "sedation": ["relaxation", "anxiety_reduction", "sleep"],
    "dissociation": ["detachment", "depersonalization", "derealization"],
    "hallucination": ["visual", "auditory", "tactile"],
    "cognitive": ["focus", "memory", "executive_function"],
    "emotional": ["empathy", "sociability", "mood"],
    "physical": ["analgesia", "motor_effects", "autonomic"],
}

# Receptor systems involved in abuse
RECEPTOR_SYSTEMS = {
    "dopamine": ["d1", "d2", "d3", "d4", "d5"],
    "serotonin": ["5ht1a", "5ht2a", "5ht2b", "5ht2c", "5ht3"],
    "opioid": ["mu", "kappa", "delta", "nociceptin"],
    "gaba": ["gaba_a", "gaba_b"],
    "glutamate": ["nmda", "ampa", "kainate", "mglur"],
    "cannabinoid": ["cb1", "cb2"],
    "cholinergic": ["alpha7", "alpha4beta2", "muscarinic"],
    "adrenergic": ["alpha1", "alpha2", "beta"],
    "sigma": ["sigma1", "sigma2"],
}

# Risk level descriptions
RISK_LEVELS = {
    "very_low": "Minimal risk of abuse or dependence",
    "low": "Some risk but generally manageable",
    "moderate": "Significant risk requiring monitoring",
    "high": "High risk requiring careful management",
    "very_high": "Extreme risk with serious concerns",
    "mild": "Mild effects/symptoms",
    "severe": "Severe effects requiring intervention",
    "life_threatening": "Critical risk requiring emergency care",
}

# Withdrawal severity levels
WITHDRAWAL_SEVERITY = {
    "mild": {
        "duration": "days",
        "symptoms": ["anxiety", "irritability", "sleep_changes"],
        "management": "outpatient monitoring",
    },
    "moderate": {
        "duration": "weeks",
        "symptoms": ["depression", "cravings", "physical_symptoms"],
        "management": "close monitoring",
    },
    "severe": {
        "duration": "weeks-months",
        "symptoms": ["severe_symptoms", "medical_complications"],
        "management": "medical supervision",
    },
    "life_threatening": {
        "duration": "variable",
        "symptoms": ["seizures", "delirium", "cardiovascular"],
        "management": "emergency care",
    },
}


def get_abuse_category_info(category: str) -> Dict:
    """Get detailed information for an abuse category."""
    return ABUSE_CATEGORIES.get(category, {})


def get_receptor_info(receptor: str) -> Dict:
    """Get receptor system information."""
    for system, receptors in RECEPTOR_SYSTEMS.items():
        if receptor in receptors:
            return {
                "system": system,
                "subtype": receptor,
                "family": receptors,
            }
    return {}


def get_risk_level_info(level: str) -> str:
    """Get risk level description."""
    return RISK_LEVELS.get(level, "Unknown risk level")


def get_withdrawal_info(severity: str) -> Dict:
    """Get withdrawal severity information."""
    return WITHDRAWAL_SEVERITY.get(severity, {})
