"""Tests for SAR analysis functionality."""

import pytest
from rdkit import Chem

from ..sar_analysis import SARAnalysisMixin


class TestCompound(SARAnalysisMixin):
    """Test class implementing SARAnalysisMixin."""

    def __init__(self, smiles=None, reference_compounds=None, primary_activity=None):
        self._sar_analysis = {}
        self.smiles = smiles
        self.reference_compounds = reference_compounds or []
        self.primary_activity = primary_activity


@pytest.fixture
def test_compounds():
    """Create test compound data."""
    return {
        # Main test compound (LSD)
        "main": {
            "smiles": "CCN(CC)C(=O)C1CN(C)C2CC3=CNC4=CC=CC(=C34)C2=C1",
            "primary_activity": 1.2,  # nM at 5-HT2A
        },
        # Reference compounds
        "references": [
            {
                # Psilocin (similar core, different substituents)
                "name": "Psilocin",
                "smiles": "CN(C)CCc1c[nH]c2cccc(O)c12",
                "primary_activity": 10.0,  # nM at 5-HT2A
            },
            {
                # DMT (similar core, minimal substituents)
                "name": "DMT",
                "smiles": "CN(C)CCc1c[nH]c2ccccc12",
                "primary_activity": 100.0,  # nM at 5-HT2A
            },
            {
                # 5-MeO-DMT (similar core, different position)
                "name": "5-MeO-DMT",
                "smiles": "CN(C)CCc1c[nH]c2ccc(OC)cc12",
                "primary_activity": 15.0,  # nM at 5-HT2A
            },
        ],
    }


def test_pharmacophore_analysis(test_compounds):
    """Test pharmacophore pattern analysis."""
    compound = TestCompound(smiles=test_compounds["main"]["smiles"])
    pharmacophores = compound._analyze_pharmacophores()
    
    # Check basic pharmacophore features
    features = {p["type"] for p in pharmacophores}
    assert "aromatic" in features
    assert "hbd" in features
    assert "hba" in features

    # Check specific features
    aromatic = next(p for p in pharmacophores if p["type"] == "aromatic")
    assert aromatic["count"] >= 2  # Indole + phenyl
    assert len(aromatic["positions"]) >= 2
    assert 0.8 <= aromatic["confidence"] <= 1.0

    # Check confidence calculation
    for feature in pharmacophores:
        assert 0 <= feature["confidence"] <= 1.0


def test_similarity_analysis(test_compounds):
    """Test structural similarity analysis."""
    compound = TestCompound(
        smiles=test_compounds["main"]["smiles"],
        reference_compounds=test_compounds["references"],
        primary_activity=test_compounds["main"]["primary_activity"]
    )
    similarity = compound._analyze_similarity()
    
    # Check overall similarity metrics
    assert "most_similar" in similarity
    assert "all_similarities" in similarity
    assert "average_similarity" in similarity

    # Check individual similarities
    similarities = similarity["all_similarities"]
    assert len(similarities) == len(test_compounds["references"])
    
    # Psilocin should be quite similar
    psilocin = next(s for s in similarities if s["compound"] == "Psilocin")
    assert psilocin["score"] >= 0.6
    assert "indole" in psilocin["shared_features"]
    assert psilocin["activity_ratio"] == pytest.approx(1.2 / 10.0)


def test_activity_cliffs(test_compounds):
    """Test activity cliff detection."""
    compound = TestCompound(
        smiles=test_compounds["main"]["smiles"],
        reference_compounds=test_compounds["references"],
        primary_activity=test_compounds["main"]["primary_activity"]
    )
    cliffs = compound._analyze_activity_cliffs()
    
    # Should find cliff with DMT (similar structure, much lower activity)
    dmt_cliff = next(
        c for c in cliffs 
        if c["compound"] == "DMT"
    )
    assert dmt_cliff["similarity"] >= 0.6
    assert dmt_cliff["activity_ratio"] >= 10  # Order of magnitude difference
    assert "structural_differences" in dmt_cliff
    assert dmt_cliff["significance"] > 0


def test_sar_patterns(test_compounds):
    """Test SAR pattern analysis."""
    compound = TestCompound(
        smiles=test_compounds["main"]["smiles"],
        reference_compounds=test_compounds["references"],
        primary_activity=test_compounds["main"]["primary_activity"]
    )
    patterns = compound._analyze_sar_patterns()
    
    # Check activity trends
    trends = patterns["activity_trends"]
    assert len(trends) > 0
    for trend in trends:
        assert "feature" in trend
        assert "avg_activity" in trend
        assert "activity_range" in trend
        assert "correlation" in trend

    # Check feature importance
    assert "feature_importance" in patterns
    assert "activity_switches" in patterns
    assert "optimal_features" in patterns


def test_substructure_analysis(test_compounds):
    """Test substructure analysis."""
    compound = TestCompound(smiles=test_compounds["main"]["smiles"])
    substructures = compound._analyze_substructures()
    
    # Check ring systems
    rings = substructures["rings"]
    assert len(rings) >= 2  # Should find indole and other rings
    for ring in rings:
        assert "size" in ring
        assert "aromatic" in ring
        assert "heteroatoms" in ring
        assert "substitution_count" in ring

    # Check functional groups
    groups = substructures["functional_groups"]
    assert "amine" in groups  # Should find tertiary amine
    assert "amide" in groups  # Should find amide

    # Check scaffolds
    scaffolds = substructures["scaffolds"]
    assert len(scaffolds) > 0
    for scaffold in scaffolds:
        assert "type" in scaffold
        assert "smiles" in scaffold
        assert "complexity" in scaffold
        assert "ring_count" in scaffold


def test_no_reference_compounds():
    """Test analysis with no reference compounds."""
    compound = TestCompound(
        smiles="CCN(CC)C(=O)C1CN(C)C2CC3=CNC4=CC=CC(=C34)C2=C1"
    )
    
    # Should still analyze pharmacophores
    pharmacophores = compound._analyze_pharmacophores()
    assert len(pharmacophores) > 0

    # Should handle missing references gracefully
    similarity = compound._analyze_similarity()
    assert not similarity["all_similarities"]
    assert similarity["average_similarity"] == 0

    # Should have empty activity cliffs
    cliffs = compound._analyze_activity_cliffs()
    assert not cliffs

    # Should still analyze substructures
    substructures = compound._analyze_substructures()
    assert substructures["rings"]
    assert substructures["functional_groups"]


def test_invalid_smiles():
    """Test analysis with invalid SMILES."""
    compound = TestCompound(smiles="invalid_smiles")
    
    # Should handle invalid structure gracefully
    pharmacophores = compound._analyze_pharmacophores()
    assert not pharmacophores

    similarity = compound._analyze_similarity()
    assert similarity["average_similarity"] == 0

    substructures = compound._analyze_substructures()
    assert not substructures["rings"]
    assert not substructures["functional_groups"]
    assert not substructures["scaffolds"]


def test_chain_analysis(test_compounds):
    """Test chain system analysis."""
    # Add compound with significant chain
    compound = TestCompound(
        smiles="CCCCCN1CCN(CCCN2C=C(C)C(=O)NC2=O)CC1"  # Buspirone
    )
    substructures = compound._analyze_substructures()
    
    # Check chain systems
    chains = substructures["chains"]
    assert len(chains) > 0
    
    # Check longest chain
    longest = max(chains, key=lambda x: x["length"])
    assert longest["length"] >= 5
    assert "C" in longest["composition"]
    assert longest["branching"] > 0


def test_functional_group_counting():
    """Test functional group detection and counting."""
    # Compound with multiple functional groups
    compound = TestCompound(
        smiles="CC(=O)NC1=CC=C(OC(=O)N(CC)CC)C=C1"
    )
    substructures = compound._analyze_substructures()
    groups = substructures["functional_groups"]
    
    # Check specific groups
    assert groups["amide"] >= 1
    assert groups["ether"] >= 1
    assert groups["ester"] >= 1
