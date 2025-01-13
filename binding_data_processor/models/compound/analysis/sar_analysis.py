"""Structure-Activity Relationship (SAR) analysis functionality for compound data.

This module provides the SARAnalyzer class that adds SAR analysis capabilities:
- Pharmacophore analysis
- Similarity searching
- Activity cliff detection
- SAR pattern analysis
- Substructure analysis
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, TYPE_CHECKING

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdRGroupDecomposition
from rdkit.Chem.Pharm2D import Generate

if TYPE_CHECKING:
    from ..base.core import CompoundData

logger = logging.getLogger(__name__)


def default_sar_analysis() -> Dict:
    """Default empty SAR analysis dictionary."""
    return {}


@dataclass
class SARAnalyzer:
    """Class providing SAR analysis capabilities."""

    _sar_analysis: Dict = field(default_factory=default_sar_analysis)

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
                pharmacophores.append(
                    {
                        "type": feature_type,
                        "count": len(positions),
                        "positions": positions,
                        "confidence": self._calculate_feature_confidence(feature_type),
                    }
                )

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

    def _calculate_similarity(self, other: "CompoundData") -> float:
        """Calculate structural similarity between compounds."""
        mol1 = Chem.MolFromSmiles(self.smiles)
        mol2 = Chem.MolFromSmiles(other.smiles)
        if not mol1 or not mol2:
            return 0.0

        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
        return DataStructs.TanimotoSimilarity(fp1, fp2)

    def _find_shared_features(self, other: "CompoundData") -> List[str]:
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

    def _calculate_activity_ratio(self, other: "CompoundData") -> Optional[float]:
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
                    cliffs.append(
                        {
                            "compound": ref.name,
                            "similarity": similarity,
                            "activity_ratio": activity_ratio,
                            "structural_differences": self._find_structural_differences(ref),
                            "significance": self._calculate_cliff_significance(similarity, activity_ratio),
                        }
                    )

        return sorted(cliffs, key=lambda x: x["significance"], reverse=True)

    def _find_structural_differences(self, other: "CompoundData") -> List[str]:
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

            trends.append(
                {
                    "feature": feature,
                    "avg_activity": sum(activities) / len(activities),
                    "activity_range": (min(activities), max(activities)),
                    "compound_count": len(compounds),
                    "correlation": self._calculate_feature_correlation(feature, compounds),
                }
            )

        return sorted(trends, key=lambda x: abs(x["correlation"]), reverse=True)

    def _analyze_feature_importance(self) -> List[Dict]:
        """Analyze importance of structural features for activity."""
        if not hasattr(self, "reference_compounds"):
            return []

        features = {}
        for ref in self.reference_compounds:
            if not hasattr(ref, "primary_activity"):
                continue

            compound_features = self._get_compound_features(ref)
            for feature in compound_features:
                if feature not in features:
                    features[feature] = {"compounds": [], "activities": []}
                features[feature]["compounds"].append(ref)
                features[feature]["activities"].append(ref.primary_activity)

        importance = []
        for feature, data in features.items():
            if len(data["activities"]) < 2:
                continue

            avg_activity = sum(data["activities"]) / len(data["activities"])
            importance.append(
                {
                    "feature": feature,
                    "avg_activity": avg_activity,
                    "compound_count": len(data["compounds"]),
                    "activity_range": (min(data["activities"]), max(data["activities"])),
                    "correlation": self._calculate_feature_correlation(feature, data["compounds"]),
                }
            )

        return sorted(importance, key=lambda x: abs(x["correlation"]), reverse=True)

    def _find_activity_switches(self) -> List[Dict]:
        """Find structural changes that cause significant activity changes."""
        if not hasattr(self, "reference_compounds"):
            return []

        switches = []
        for i, ref1 in enumerate(self.reference_compounds):
            for ref2 in self.reference_compounds[i + 1 :]:
                similarity = self._calculate_similarity(ref1)
                if similarity < 0.7:  # Only consider similar compounds
                    continue

                activity_ratio = self._calculate_activity_ratio(ref2)
                if not activity_ratio or abs(activity_ratio) < 5:  # Significant activity difference
                    continue

                differences = self._find_structural_differences(ref2)
                if differences:
                    switches.append(
                        {
                            "compounds": [ref1.name, ref2.name],
                            "similarity": similarity,
                            "activity_ratio": activity_ratio,
                            "structural_changes": differences,
                            "significance": abs(activity_ratio) * similarity,
                        }
                    )

        return sorted(switches, key=lambda x: x["significance"], reverse=True)

    def _identify_optimal_features(self) -> Dict:
        """Identify structural features associated with optimal activity."""
        if not hasattr(self, "reference_compounds"):
            return {}

        # Group compounds by activity level
        high_activity = []
        low_activity = []
        threshold = self._calculate_activity_threshold()

        for ref in self.reference_compounds:
            if not hasattr(ref, "primary_activity"):
                continue

            if ref.primary_activity > threshold:
                high_activity.append(ref)
            else:
                low_activity.append(ref)

        # Analyze features in high vs low activity compounds
        high_features = self._analyze_group_features(high_activity)
        low_features = self._analyze_group_features(low_activity)

        # Identify distinguishing features
        optimal_features = []
        for feature, high_stats in high_features.items():
            if feature not in low_features:
                continue

            low_stats = low_features[feature]
            if high_stats["avg_activity"] > 2 * low_stats["avg_activity"]:
                optimal_features.append(
                    {
                        "feature": feature,
                        "high_activity_stats": high_stats,
                        "low_activity_stats": low_stats,
                        "enrichment": high_stats["avg_activity"] / low_stats["avg_activity"],
                    }
                )

        return {
            "optimal_features": sorted(optimal_features, key=lambda x: x["enrichment"], reverse=True),
            "activity_threshold": threshold,
            "high_activity_count": len(high_activity),
            "low_activity_count": len(low_activity),
        }

    def _calculate_activity_threshold(self) -> float:
        """Calculate activity threshold for optimal feature analysis."""
        activities = []
        for ref in self.reference_compounds:
            if hasattr(ref, "primary_activity"):
                activities.append(ref.primary_activity)

        if not activities:
            return 0.0

        return sum(activities) / len(activities)  # Use mean as threshold

    def _analyze_group_features(self, compounds: List["CompoundData"]) -> Dict:
        """Analyze features for a group of compounds."""
        features = {}
        for compound in compounds:
            compound_features = self._get_compound_features(compound)
            for feature in compound_features:
                if feature not in features:
                    features[feature] = {
                        "count": 0,
                        "total_activity": 0.0,
                        "compounds": [],
                    }
                features[feature]["count"] += 1
                features[feature]["total_activity"] += compound.primary_activity
                features[feature]["compounds"].append(compound)

        # Calculate statistics
        for feature in features:
            count = features[feature]["count"]
            features[feature]["avg_activity"] = features[feature]["total_activity"] / count
            features[feature]["frequency"] = count / len(compounds)

        return features

    def _get_compound_features(self, compound: "CompoundData") -> Set[str]:
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

    def _calculate_feature_correlation(self, feature: str, compounds: List["CompoundData"]) -> float:
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
                "heteroatoms": [mol.GetAtomWithIdx(i).GetSymbol() for i in ring_atoms if mol.GetAtomWithIdx(i).GetSymbol() != "C"],
                "substitution_count": sum(len(set(mol.GetAtomWithIdx(i).GetNeighbors()) - set(ring_atoms)) for i in ring_atoms),
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
            "composition": [mol.GetAtomWithIdx(i).GetSymbol() for i in chain_atoms],
            "branching": sum(len(mol.GetAtomWithIdx(i).GetNeighbors()) > 2 for i in chain_atoms),
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
        """Analyze molecular scaffolds and R-groups."""
        from rdkit.Chem.Scaffolds import MurckoScaffold

        scaffolds = []

        # Get Murcko scaffold
        scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold_mol:
            scaffold_info = {
                "type": "murcko",
                "smiles": Chem.MolToSmiles(scaffold_mol),
                "complexity": scaffold_mol.GetNumAtoms(),
                "ring_count": scaffold_mol.GetRingInfo().NumRings(),
            }

            # Analyze R-groups
            try:
                rgroup_params = rdRGroupDecomposition.RGroupDecompositionParameters()
                rgroup_params.removeHydrogensPostMatch = True
                rgroup_params.onlyMatchAtRGroups = False

                decomp = rdRGroupDecomposition.RGroupDecomposition(scaffold_mol, rgroup_params)
                decomp.Add(mol)
                if decomp.Process():
                    rgroups = decomp.GetRGroupsAsColumns()
                    if rgroups:
                        scaffold_info["rgroups"] = {}
                        for rgroup_label, rgroup_mols in rgroups.items():
                            if rgroup_mols and rgroup_mols[0] is not None:
                                scaffold_info["rgroups"][rgroup_label] = {
                                    "smiles": Chem.MolToSmiles(rgroup_mols[0]),
                                    "size": rgroup_mols[0].GetNumAtoms(),
                                    "complexity": self._calculate_rgroup_complexity(rgroup_mols[0]),
                                }
            except Exception as e:
                logger.debug(f"Error in R-group decomposition: {str(e)}")

            scaffolds.append(scaffold_info)

        return scaffolds

    def _calculate_rgroup_complexity(self, rgroup_mol: Chem.Mol) -> float:
        """Calculate complexity score for an R-group."""
        try:
            # Basic complexity factors
            num_atoms = rgroup_mol.GetNumAtoms()
            num_bonds = rgroup_mol.GetNumBonds()
            num_rings = rgroup_mol.GetRingInfo().NumRings()

            # Normalize and combine scores
            atom_score = min(1.0, num_atoms / 10)  # Smaller scale for R-groups
            bond_score = min(1.0, num_bonds / 12)
            ring_score = min(1.0, num_rings / 2)

            # Weighted average
            weights = [1.0, 1.0, 1.5]
            total_weight = sum(weights)
            score = sum([atom_score * weights[0], bond_score * weights[1], ring_score * weights[2]])
            return float(score / total_weight)
        except Exception:
            return 0.0
