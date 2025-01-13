"""Base classes for pharmacophore detection and generation."""

from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

from .features import PharmacophoreFeature, FEATURE_PATTERNS


@dataclass
class PharmacophorePoint:
    """A single point in a pharmacophore model."""

    feature_type: str
    x: float
    y: float
    z: float
    radius: float = 1.0
    enabled: bool = True
    weight: float = 1.0
    vector: Optional[Tuple[float, float, float]] = None


class PharmacophoreGenerator:
    """Base class for generating pharmacophore models from molecules."""

    def __init__(self):
        """Initialize the pharmacophore generator."""
        self.features = FEATURE_PATTERNS

    def generate(self, mol: Chem.Mol) -> List[PharmacophoreFeature]:
        """Generate pharmacophore features for a molecule.

        Args:
            mol: RDKit molecule to analyze

        Returns:
            List of pharmacophore features found in the molecule
        """
        if not mol:
            return []

        # Generate 3D coordinates if not present
        if not mol.GetNumConformers():
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)

        features = []
        for feature_type, patterns in self.features.items():
            for pattern in patterns:
                # Find matches for the SMARTS pattern
                pattern_mol = Chem.MolFromSmarts(pattern)
                if not pattern_mol:
                    continue

                matches = mol.GetSubstructMatches(pattern_mol)
                for match in matches:
                    # Calculate centroid of matched atoms
                    if not match:
                        continue

                    positions = []
                    for atom_idx in match:
                        pos = mol.GetConformer().GetAtomPosition(atom_idx)
                        positions.append((pos.x, pos.y, pos.z))

                    if not positions:
                        continue

                    # Calculate average position
                    x = sum(p[0] for p in positions) / len(positions)
                    y = sum(p[1] for p in positions) / len(positions)
                    z = sum(p[2] for p in positions) / len(positions)

                    feature = PharmacophoreFeature(feature_type=feature_type, atoms=match, position=(x, y, z))
                    features.append(feature)

        return features

    def get_feature_types(self) -> Set[str]:
        """Get the set of all available feature types.

        Returns:
            Set of feature type strings
        """
        return set(self.features.keys())

    def get_patterns(self, feature_type: str) -> List[str]:
        """Get SMARTS patterns for a specific feature type.

        Args:
            feature_type: Type of pharmacophore feature

        Returns:
            List of SMARTS patterns for that feature
        """
        return self.features.get(feature_type, [])
