"""Toxicity type definitions and constants.

This module defines:
1. General toxicity categories
2. Organ-specific toxicity
3. Toxicity mechanisms
4. Safety concerns
5. Adverse effects
6. Drug interactions
"""

from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass
class ToxicityType:
    """Toxicity type definition."""

    id: str
    name: str
    description: str
    severity_levels: List[str]
    mechanisms: Optional[List[str]] = None
    related_effects: Optional[List[str]] = None


# General toxicity categories
TOXICITY_TYPES = {
    "acute_toxicity": ToxicityType(
        id="acute_toxicity",
        name="Acute Toxicity",
        description="Adverse effects occurring shortly after exposure",
        severity_levels=["minimal", "mild", "moderate", "severe", "lethal"],
        mechanisms=[
            "direct cellular damage",
            "organ system failure",
            "metabolic disruption",
        ],
    ),
    "chronic_toxicity": ToxicityType(
        id="chronic_toxicity",
        name="Chronic Toxicity",
        description="Long-term adverse effects from repeated exposure",
        severity_levels=["minimal", "mild", "moderate", "severe", "life_threatening"],
        mechanisms=[
            "cumulative damage",
            "chronic inflammation",
            "tissue degeneration",
        ],
    ),
}

# Organ-specific toxicity
ORGAN_TOXICITY = {
    "cardiotoxicity": ToxicityType(
        id="cardiotoxicity",
        name="Cardiotoxicity",
        description="Toxic effects on the heart and cardiovascular system",
        severity_levels=["minimal", "mild", "moderate", "severe", "critical"],
        mechanisms=[
            "ion channel disruption",
            "oxidative stress",
            "mitochondrial dysfunction",
        ],
        related_effects=[
            "arrhythmia",
            "heart failure",
            "hypertension",
        ],
    ),
    "hepatotoxicity": ToxicityType(
        id="hepatotoxicity",
        name="Hepatotoxicity",
        description="Liver damage and dysfunction",
        severity_levels=["minimal", "mild", "moderate", "severe", "liver_failure"],
        mechanisms=[
            "metabolic stress",
            "oxidative damage",
            "immune-mediated injury",
        ],
        related_effects=[
            "elevated enzymes",
            "cholestasis",
            "fibrosis",
        ],
    ),
    "neurotoxicity": ToxicityType(
        id="neurotoxicity",
        name="Neurotoxicity",
        description="Damage to nervous system structure or function",
        severity_levels=["minimal", "mild", "moderate", "severe", "critical"],
        mechanisms=[
            "neurotransmitter disruption",
            "axonal degeneration",
            "synaptic dysfunction",
        ],
        related_effects=[
            "cognitive impairment",
            "motor dysfunction",
            "seizures",
        ],
    ),
    "nephrotoxicity": ToxicityType(
        id="nephrotoxicity",
        name="Nephrotoxicity",
        description="Kidney damage and dysfunction",
        severity_levels=["minimal", "mild", "moderate", "severe", "renal_failure"],
        mechanisms=[
            "tubular injury",
            "glomerular damage",
            "vascular injury",
        ],
        related_effects=[
            "proteinuria",
            "electrolyte imbalance",
            "reduced filtration",
        ],
    ),
}

# Specific toxicity mechanisms
TOXICITY_MECHANISMS = {
    "oxidative_stress": ToxicityType(
        id="oxidative_stress",
        name="Oxidative Stress",
        description="Generation of harmful reactive oxygen species",
        severity_levels=["low", "moderate", "high", "severe"],
        mechanisms=[
            "free radical generation",
            "lipid peroxidation",
            "protein oxidation",
        ],
    ),
    "mitochondrial_toxicity": ToxicityType(
        id="mitochondrial_toxicity",
        name="Mitochondrial Toxicity",
        description="Disruption of cellular energy production",
        severity_levels=["minimal", "mild", "moderate", "severe"],
        mechanisms=[
            "electron transport disruption",
            "ATP depletion",
            "membrane damage",
        ],
    ),
    "immunotoxicity": ToxicityType(
        id="immunotoxicity",
        name="Immunotoxicity",
        description="Adverse effects on immune system function",
        severity_levels=["minimal", "mild", "moderate", "severe"],
        mechanisms=[
            "immune suppression",
            "hypersensitivity",
            "autoimmunity",
        ],
    ),
}

# Drug interaction risks
INTERACTION_RISKS = {
    "cyp_inhibition": ToxicityType(
        id="cyp_inhibition",
        name="CYP450 Enzyme Inhibition",
        description="Interference with drug metabolism",
        severity_levels=["low", "moderate", "high", "severe"],
        mechanisms=[
            "competitive inhibition",
            "mechanism-based inhibition",
            "mixed inhibition",
        ],
    ),
    "transporter_interaction": ToxicityType(
        id="transporter_interaction",
        name="Drug Transporter Interaction",
        description="Altered drug absorption or distribution",
        severity_levels=["low", "moderate", "high", "severe"],
        mechanisms=[
            "p-glycoprotein inhibition",
            "organic anion transport",
            "active efflux",
        ],
    ),
}

# Safety concerns
SAFETY_CONCERNS = {
    "reproductive_toxicity": ToxicityType(
        id="reproductive_toxicity",
        name="Reproductive Toxicity",
        description="Effects on reproductive function and development",
        severity_levels=["minimal", "mild", "moderate", "severe"],
        mechanisms=[
            "endocrine disruption",
            "developmental toxicity",
            "genetic damage",
        ],
    ),
    "genotoxicity": ToxicityType(
        id="genotoxicity",
        name="Genotoxicity",
        description="DNA damage and mutation risk",
        severity_levels=["low", "moderate", "high", "severe"],
        mechanisms=[
            "DNA strand breaks",
            "chromosomal damage",
            "point mutations",
        ],
    ),
    "carcinogenicity": ToxicityType(
        id="carcinogenicity",
        name="Carcinogenicity",
        description="Cancer-causing potential",
        severity_levels=["minimal", "low", "moderate", "high"],
        mechanisms=[
            "mutagenesis",
            "epigenetic changes",
            "cell proliferation",
        ],
    ),
}

# Combined toxicity categories
ALL_TOXICITY_TYPES = {
    **TOXICITY_TYPES,
    **ORGAN_TOXICITY,
    **TOXICITY_MECHANISMS,
    **INTERACTION_RISKS,
    **SAFETY_CONCERNS,
}

# Severity level definitions
SEVERITY_LEVELS = {
    "minimal": "Minimal risk or effects, generally reversible",
    "mild": "Mild effects, fully reversible, no lasting damage",
    "moderate": "Moderate effects, may require intervention",
    "severe": "Severe effects, significant intervention needed",
    "critical": "Life-threatening effects requiring immediate action",
    "lethal": "Fatal or potentially fatal effects",
}

# Risk level definitions
RISK_LEVELS = {
    "very_low": "Very low probability of adverse effects",
    "low": "Low probability of adverse effects",
    "moderate": "Moderate probability of adverse effects",
    "high": "High probability of adverse effects",
    "very_high": "Very high probability of adverse effects",
}
