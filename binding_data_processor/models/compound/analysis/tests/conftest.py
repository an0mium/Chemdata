"""Shared test fixtures for compound analysis tests."""

import pytest


@pytest.fixture
def lsd_data():
    """Create test data for LSD."""
    return {
        "smiles": "CCN(CC)C(=O)C1CN(C)C2CC3=CNC4=CC=CC(=C34)C2=C1",
        "name": "LSD",
        "cas_number": "50-37-3",
        "targets": [
            {
                "common_name": "5-HT2A",
                "affinity_value": 1.2,
                "affinity_type": "Ki",
                "affinity_unit": "nM",
                "activity_type": "agonist",
                "confidence": 0.95,
                "assay_details": {
                    "family": "serotonin",
                    "effects": ["psychedelic", "mood"],
                    "mechanisms": ["receptor activation"],
                },
            },
            {
                "common_name": "5-HT1A",
                "affinity_value": 8.5,
                "affinity_type": "Ki",
                "affinity_unit": "nM",
                "activity_type": "partial agonist",
                "confidence": 0.9,
                "assay_details": {
                    "family": "serotonin",
                    "effects": ["mood", "anxiety"],
                    "mechanisms": ["receptor modulation"],
                },
            },
        ],
        "primary_activity": 1.2,  # nM at 5-HT2A
        "mechanism_of_action": "serotonin 5-HT2A receptor agonist",
        "effect_profile": {
            "psychedelic": (0.9, 0.95),
            "mood": (0.8, 0.9),
            "cognition": (0.6, 0.7),
            "perception": (0.85, 0.9),
            "anxiety": (-0.3, 0.6),
        },
        "safety_profile": {
            "known_risks": [
                {
                    "description": "Psychological distress",
                    "severity": "moderate",
                    "confidence": 0.9,
                    "mechanism": "5-HT2A activation",
                    "conditions": ["predisposed individuals"],
                },
            ],
            "known_interactions": [
                {
                    "description": "Lithium",
                    "severity": "severe",
                    "confidence": 0.95,
                    "mechanism": "serotonergic effects",
                    "recommendations": ["contraindicated"],
                },
            ],
        },
        "experimental_data": {
            "bioavailability": 0.7,
            "half_life": 3.6,  # hours
            "volume_distribution": 0.8,  # L/kg
        },
    }


@pytest.fixture
def reference_compounds():
    """Create test data for reference compounds."""
    return [
        {
            "name": "Psilocin",
            "smiles": "CN(C)CCc1c[nH]c2cccc(O)c12",
            "primary_activity": 10.0,  # nM at 5-HT2A
            "targets": [
                {
                    "common_name": "5-HT2A",
                    "affinity_value": 10.0,
                    "affinity_type": "Ki",
                    "affinity_unit": "nM",
                    "activity_type": "agonist",
                    "confidence": 0.9,
                },
            ],
        },
        {
            "name": "DMT",
            "smiles": "CN(C)CCc1c[nH]c2ccccc12",
            "primary_activity": 100.0,  # nM at 5-HT2A
            "targets": [
                {
                    "common_name": "5-HT2A",
                    "affinity_value": 100.0,
                    "affinity_type": "Ki",
                    "affinity_unit": "nM",
                    "activity_type": "agonist",
                    "confidence": 0.85,
                },
            ],
        },
        {
            "name": "5-MeO-DMT",
            "smiles": "CN(C)CCc1c[nH]c2ccc(OC)cc12",
            "primary_activity": 15.0,  # nM at 5-HT2A
            "targets": [
                {
                    "common_name": "5-HT2A",
                    "affinity_value": 15.0,
                    "affinity_type": "Ki",
                    "affinity_unit": "nM",
                    "activity_type": "agonist",
                    "confidence": 0.9,
                },
            ],
        },
    ]


@pytest.fixture
def drug_like_properties():
    """Create test data for drug-like compound properties."""
    return {
        "molecular_weight": 320.0,
        "logp": 2.8,
        "hbd": 2,
        "hba": 5,
        "tpsa": 90.0,
        "rotatable_bonds": 6,
        "charge": 0,
        "stereocenter_count": 2,
        "ring_count": 3,
        "atom_count": 45,
    }


@pytest.fixture
def non_drug_like_properties():
    """Create test data for non-drug-like compound properties."""
    return {
        "molecular_weight": 650.0,
        "logp": 6.5,
        "hbd": 8,
        "hba": 12,
        "tpsa": 180.0,
        "rotatable_bonds": 15,
        "charge": 2,
        "stereocenter_count": 0,
        "ring_count": 8,
        "atom_count": 85,
    }


@pytest.fixture
def experimental_data():
    """Create test experimental data."""
    return {
        "bioavailability": 0.7,
        "half_life": 3.6,  # hours
        "volume_distribution": 0.8,  # L/kg
        "clearance": 1.2,  # L/h/kg
        "protein_binding": 0.85,  # fraction bound
        "metabolism": {
            "cyp_inhibition": {
                "CYP2D6": 0.3,  # IC50 μM
                "CYP3A4": 1.5,
            },
            "cyp_induction": {
                "CYP1A2": False,
                "CYP3A4": False,
            },
        },
    }
