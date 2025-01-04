"""Structure-Activity Relationship (SAR) analysis functionality for compound data.

This module provides the SARAnalysisMixin class that adds SAR analysis capabilities:
- Pharmacophore analysis
- Similarity searching
- Activity cliff detection
- SAR pattern analysis
- Substructure analysis
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Pharm2D import Generate

from ...compound import CompoundData


@dataclass
class SARAnalysisMixin:
    """Mixin class adding SAR analysis capabilities."""

    _sar_analysis: Dict = field(default_factory=dict)

    def analyze_sar(self) -> Dict:
        """Analyze structure-activity relationships."""
        analysis = {
            "pharmacophores": self._analyze_pharmacophores(),
            "similarity": self._analyze_similarity(),
            "activity_cliffs": self._analyze_activity_cliffs(),
            "sar_patterns": self._analyze_sar_patterns(),
            "substructures": self._analyze_substructures(),
        }
        self._sar_analysis = analysis
        return analysis

    def _analyze_pharmacophores(self) -> List[Dict]:
        """Analyze pharmacophore patterns."""
        pharmacophores = []

        # Get molecule object
        mol = Chem.MolFromSmiles(self.smiles)
        if not mol:
            return pharmacophores

        # Generate 2D pharmacophore fingerprint
        fp = Generate.Gen2DFingerprint(mol)
        if not fp:
            return pharmacophores

        # Analyze pharmacophore features
        features = {
            "aromatic": self._find_aromatic_centers(mol),
            "hbd": self._find_hb_donors(mol),
            "hba": self._find_hb_acceptors(mol),
            "charged": self._find_charged_groups(mol),
            "hydrophobic": self._find_hydrophobic_regions(mol),
        }

        # Generate pharmacophore hypotheses
        for feature_type, positions in features.items():
            if positions:
                pharmacophores.append({
                    "type": feature_type,
                    "count": len(positions),
                    "positions": positions,
                    "confidence": self._calculate_feature_confidence(feature_type),
                })

        return pharmacophores

    def _find_aromatic_centers(self, mol: Chem.Mol) -> List[int]:
        """Find aromatic ring centers."""
        aromatic_atoms = []
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic():
                aromatic_atoms.append(atom.GetIdx())
        return aromatic_atoms

    def _find_hb_donors(self, mol: Chem.Mol) -> List[int]:
        """Find hydrogen bond donors."""
        donors = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in ["N", "O"] and atom.GetTotalNumHs() > 0:
                donors.append(atom.GetIdx())
        return donors

    def _find_hb_acceptors(self, mol: Chem.Mol) -> List[int]:
        """Find hydrogen bond acceptors."""
        acceptors = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in ["N", "O", "F"]:
                acceptors.append(atom.GetIdx())
        return acceptors

    def _find_charged_groups(self, mol: Chem.Mol) -> List[int]:
        """Find charged groups."""
        charged = []
        for atom in mol.GetAtoms():
            if atom.GetFormalCharge() != 0:
                charged.append(atom.GetIdx())
        return charged

    def _find_hydrophobic_regions(self, mol: Chem.Mol) -> List[int]:
        """Find hydrophobic regions."""
        hydrophobic = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "C" and not atom.GetIsAromatic():
                hydrophobic.append(atom.GetIdx())
        return hydrophobic

    def _calculate_feature_confidence(self, feature_type: str) -> float:
        """Calculate confidence in pharmacophore feature."""
        # Base confidence
        confidence = 0.8

        # Adjust based on experimental validation
        if hasattr(self, "experimental_data"):
            confidence *= 1.2

        # Adjust based on feature type
        type_adjustments = {
            "aromatic": 1.1,  # More reliable
            "hbd": 1.0,
            "hba": 1.0,
            "charged": 1.1,  # More reliable
            "hydrophobic": 0.9,  # Less reliable
        }
        confidence *= type_adjustments.get(feature_type, 1.0)

        return min(confidence, 1.0)

    def _analyze_similarity(self) -> Dict:
        """Analyze structural similarity to known actives."""
        if not hasattr(self, "reference_compounds"):
            return {}

        similarities = []
        for ref in self.reference_compounds:
            similarity = {
                "compound": ref.name,
                "score": self._calculate_similarity(ref),
                "shared_features": self._find_shared_features(ref),
                "activity_ratio": self._calculate_activity_ratio(ref),
            }
            similarities.append(similarity)

        return {
            "most_similar": max(similarities, key=lambda x: x["score"]),
            "all_similarities": similarities,
            "average_similarity": sum(s["score"] for s in similarities) / len(similarities),
        }

    def _calculate_similarity(self, other: 'CompoundData') -> float:
        """Calculate structural similarity between compounds."""
        mol1 = Chem.MolFromSmiles(self.smiles)
        mol2 = Chem.MolFromSmiles(other.smiles)
        if not mol1 or not mol2:
            return 0.0

        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
        return DataStructs.TanimotoSimilarity(fp1, fp2)

    def _find_shared_features(self, other: 'CompoundData') -> List[str]:
        """Find pharmacophore features shared with another compound."""
        shared = []
        mol1 = Chem.MolFromSmiles(self.smiles)
        mol2 = Chem.MolFromSmiles(other.smiles)
        if not mol1 or not mol2:
            return shared

        # Compare features
        features1 = {
            "aromatic": set(self._find_aromatic_centers(mol1)),
            "hbd": set(self._find_hb_donors(mol1)),
            "hba": set(self._find_hb_acceptors(mol1)),
            "charged": set(self._find_charged_groups(mol1)),
            "hydrophobic": set(self._find_hydrophobic_regions(mol1)),
        }

        features2 = {
            "aromatic": set(self._find_aromatic_centers(mol2)),
            "hbd": set(self._find_hb_donors(mol2)),
            "hba": set(self._find_hb_acceptors(mol2)),
            "charged": set(self._find_charged_groups(mol2)),
            "hydrophobic": set(self._find_hydrophobic_regions(mol2)),
        }

        for feature_type in features1:
            if features1[feature_type] & features2[feature_type]:
                shared.append(feature_type)

        return shared

    def _calculate_activity_ratio(self, other: 'CompoundData') -> Optional[float]:
        """Calculate activity ratio between compounds."""
        if not hasattr(self, "primary_activity") or not hasattr(other, "primary_activity"):
            return None

        if self.primary_activity <= 0 or other.primary_activity <= 0:
            return None

        return self.primary_activity / other.primary_activity

    def _analyze_activity_cliffs(self) -> List[Dict]:
        """Analyze activity cliffs with similar compounds."""
        if not hasattr(self, "reference_compounds"):
            return []

        cliffs = []
        for ref in self.reference_compounds:
            similarity = self._calculate_similarity(ref)
            if similarity >= 0.7:  # High structural similarity
                activity_ratio = self._calculate_activity_ratio(ref)
                if activity_ratio and abs(activity_ratio) >= 10:  # Significant activity difference
                    cliffs.append({
                        "compound": ref.name,
                        "similarity": similarity,
                        "activity_ratio": activity_ratio,
                        "structural_differences": self._find_structural_differences(ref),
                        "significance": self._calculate_cliff_significance(similarity, activity_ratio),
                    })

        return sorted(cliffs, key=lambda x: x["significance"], reverse=True)

    def _find_structural_differences(self, other: 'CompoundData') -> List[str]:
        """Find key structural differences with another compound."""
        differences = []
        mol1 = Chem.MolFromSmiles(self.smiles)
        mol2 = Chem.MolFromSmiles(other.smiles)
        if not mol1 or not mol2:
            return differences

        # Compare basic properties
        if abs(mol1.GetNumAtoms() - mol2.GetNumAtoms()) > 2:
            differences.append("size")
        if abs(mol1.GetNumRotatableBonds() - mol2.GetNumRotatableBonds()) > 2:
            differences.append("flexibility")
        if abs(mol1.GetNumHBA() - mol2.GetNumHBA()) > 1:
            differences.append("h-bond acceptors")
        if abs(mol1.GetNumHBD() - mol2.GetNumHBD()) > 1:
            differences.append("h-bond donors")

        return differences

    def _calculate_cliff_significance(self, similarity: float, activity_ratio: float) -> float:
        """Calculate significance of activity cliff."""
        return similarity * abs(activity_ratio) / 10

    def _analyze_sar_patterns(self) -> Dict:
        """Analyze SAR patterns across similar compounds."""
        if not hasattr(self, "reference_compounds"):
            return {}

        patterns = {
            "activity_trends": self._analyze_activity_trends(),
            "feature_importance": self._analyze_feature_importance(),
            "activity_switches": self._find_activity_switches(),
            "optimal_features": self._identify_optimal_features(),
        }

        return patterns

    def _analyze_activity_trends(self) -> List[Dict]:
        """Analyze trends in activity vs structural features."""
        trends = []
        if not hasattr(self, "reference_compounds"):
            return trends

        # Group compounds by feature presence
        feature_groups = {}
        for ref in self.reference_compounds:
            features = self._get_compound_features(ref)
            for feature in features:
                if feature not in feature_groups:
                    feature_groups[feature] = []
                feature_groups[feature].append(ref)

        # Analyze activity distribution for each feature
        for feature, compounds in feature_groups.items():
            if len(compounds) < 2:
                continue

            activities = [c.primary_activity for c in compounds if hasattr(c, "primary_activity")]
            if not activities:
                continue

            trends.append({
                "feature": feature,
                "avg_activity": sum(activities) / len(activities),
                "activity_range": (min(activities), max(activities)),
                "compound_count": len(compounds),
                "correlation": self._calculate_feature_correlation(feature, compounds),
            })

        return sorted(trends, key=lambda x: abs(x["correlation"]), reverse=True)

    def _get_compound_features(self, compound: 'CompoundData') -> Set[str]:
        """Get set of structural features for a compound."""
        features = set()
        mol = Chem.MolFromSmiles(compound.smiles)
        if not mol:
            return features

        # Add basic features
        if self._find_aromatic_centers(mol):
            features.add("aromatic")
        if self._find_hb_donors(mol):
            features.add("hbd")
        if self._find_hb_acceptors(mol):
            features.add("hba")
        if self._find_charged_groups(mol):
            features.add("charged")
        if self._find_hydrophobic_regions(mol):
            features.add("hydrophobic")

        return features

    def _calculate_feature_correlation(self, feature: str, compounds: List['CompoundData']) -> float:
        """Calculate correlation between feature and activity."""
        activities = []
        has_feature = []
        for compound in compounds:
            if not hasattr(compound, "primary_activity"):
                continue
            activities.append(compound.primary_activity)
            features = self._get_compound_features(compound)
            has_feature.append(1.0 if feature in features else 0.0)

        if len(activities) < 2:
            return 0.0

        return self._calculate_correlation(has_feature, activities)

    def _calculate_correlation(self, x: List[float], y: List[float]) -> float:
        """Calculate Pearson correlation coefficient."""
        if len(x) != len(y) or len(x) < 2:
            return 0.0

        n = len(x)
        sum_x = sum(x)
        sum_y = sum(y)
        sum_xy = sum(i * j for i, j in zip(x, y))
        sum_x2 = sum(i * i for i in x)
        sum_y2 = sum(i * i for i in y)

        numerator = n * sum_xy - sum_x * sum_y
        denominator = ((n * sum_x2 - sum_x * sum_x) * (n * sum_y2 - sum_y * sum_y)) ** 0.5

        return numerator / denominator if denominator != 0 else 0.0

    def _analyze_substructures(self) -> Dict:
        """Analyze important substructures."""
        mol = Chem.MolFromSmiles(self.smiles)
        if not mol:
            return {}

        return {
            "rings": self._analyze_ring_systems(mol),
            "chains": self._analyze_chain_systems(mol),
            "functional_groups": self._analyze_functional_groups(mol),
            "scaffolds": self._analyze_scaffolds(mol),
        }

    def _analyze_ring_systems(self, mol: Chem.Mol) -> List[Dict]:
        """Analyze ring systems."""
        rings = []
        ring_info = mol.GetRingInfo()
        
        for ring_atoms in ring_info.AtomRings():
            ring = {
                "size": len(ring_atoms),
                "aromatic": all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring_atoms),
                "heteroatoms": [
                    mol.GetAtomWithIdx(i).GetSymbol()
                    for i in ring_atoms
                    if mol.GetAtomWithIdx(i).GetSymbol() != "C"
                ],
                "substitution_count": sum(
                    len(set(mol.GetAtomWithIdx(i).GetNeighbors()) - set(ring_atoms))
                    for i in ring_atoms
                ),
            }
            rings.append(ring)

        return rings

    def _analyze_chain_systems(self, mol: Chem.Mol) -> List[Dict]:
        """Analyze chain systems."""
        chains = []
        for atom in mol.GetAtoms():
            if atom.IsInRing():
                continue
            
            chain = self._trace_chain(mol, atom)
            if chain:
                chains.append(chain)

        return chains

    def _trace_chain(self, mol: Chem.Mol, start_atom: Chem.Atom) -> Optional[Dict]:
        """Trace a chain starting from an atom."""
        visited = set()
        current = start_atom
        chain_atoms = []

        while current and not current.IsInRing():
            visited.add(current.GetIdx())
            chain_atoms.append(current.GetIdx())
            
            # Find next non-ring neighbor
            next_atom = None
            for neighbor in current.GetNeighbors():
                if neighbor.GetIdx() not in visited and not neighbor.IsInRing():
                    next_atom = neighbor
                    break
            current = next_atom

        if len(chain_atoms) < 2:
            return None

        return {
            "length": len(chain_atoms),
            "composition": [
                mol.GetAtomWithIdx(i).GetSymbol()
                for i in chain_atoms
            ],
            "branching": sum(
                len(mol.GetAtomWithIdx(i).GetNeighbors()) > 2
                for i in chain_atoms
            ),
        }

    def _analyze_functional_groups(self, mol: Chem.Mol) -> Dict:
        """Analyze functional groups."""
        groups = {
            "alcohol": self._count_pattern(mol, "OH"),
            "amine": self._count_pattern(mol, "N"),
            "carbonyl": self._count_pattern(mol, "C(=O)"),
            "carboxyl": self._count_pattern(mol, "C(=O)O"),
            "ether": self._count_pattern(mol, "COC"),
            "ester": self._count_pattern(mol, "C(=O)OC"),
            "amide": self._count_pattern(mol, "C(=O)N"),
        }
        return {k: v for k, v in groups.items() if v > 0}

    def _count_pattern(self, mol: Chem.Mol, pattern: str) -> int:
        """Count occurrences of a SMARTS pattern."""
        pattern_mol = Chem.MolFromSmarts(pattern)
        if not pattern_mol:
            return 0
        return len(mol.GetSubstructMatches(pattern_mol))

    def _analyze_scaffolds(self, mol: Chem.Mol) -> List[Dict]:
        """Analyze molecular scaffolds."""
        from rdkit.Chem.Scaffolds import MurckoScaffold

        scaffolds = []
        
        # Get Murcko scaffold
        scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold_mol:
            scaffolds.append({
                "type": "murcko",
                "smiles": Chem.MolToSmiles(scaffold_mol),
                "complexity": scaffold_mol.GetNumAtoms(),
                "ring_count": scaffold_mol.GetRingInfo().NumRings(),
            })

        return scaffolds
