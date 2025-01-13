"""Molecular descriptor calculation."""

import logging
from typing import Dict, List, Optional, Set, Union

import numpy as np
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    rdMolDescriptors,
    Scaffolds,
)

from .base import MLProcessor


class DescriptorGenerator(MLProcessor):
    """Generate molecular descriptors."""

    # Descriptor definitions with metadata
    DESCRIPTOR_TYPES = {
        "constitutional": {
            "description": "Basic molecular properties",
            "functions": [
                (Descriptors.ExactMolWt, "molecular_weight"),
                (rdMolDescriptors.CalcNumAtoms, "num_atoms"),
                (lambda mol: mol.GetNumBonds(), "num_bonds"),
                (Descriptors.NumRotatableBonds, "num_rotatable_bonds"),
                (Descriptors.NumHAcceptors, "num_h_acceptors"),
                (Descriptors.NumHDonors, "num_h_donors"),
                (Descriptors.RingCount, "ring_count"),
                (Descriptors.NumAromaticRings, "num_aromatic_rings"),
                (Descriptors.NumAliphaticRings, "num_aliphatic_rings"),
            ],
        },
        "topological": {
            "description": "2D molecular shape and connectivity",
            "functions": [
                (Descriptors.BertzCT, "complexity"),
                (Descriptors.HallKierAlpha, "hall_kier_alpha"),
                (Descriptors.Kappa1, "kappa1"),
                (Descriptors.Kappa2, "kappa2"),
                (Descriptors.Kappa3, "kappa3"),
                (Descriptors.Chi0v, "chi0v"),
                (Descriptors.Chi1v, "chi1v"),
                (Descriptors.Chi2v, "chi2v"),
                (Descriptors.Chi3v, "chi3v"),
                (Descriptors.Chi4v, "chi4v"),
            ],
        },
        "electronic": {
            "description": "Electronic properties",
            "functions": [
                (Descriptors.TPSA, "tpsa"),
                (Descriptors.MolLogP, "logp"),
                (Descriptors.MolMR, "molar_refractivity"),
                (Descriptors.MaxPartialCharge, "max_partial_charge"),
                (Descriptors.MinPartialCharge, "min_partial_charge"),
                (Descriptors.MaxAbsPartialCharge, "max_abs_partial_charge"),
                (Descriptors.MinAbsPartialCharge, "min_abs_partial_charge"),
            ],
        },
        "geometric": {
            "description": "3D molecular shape",
            "functions": [
                (AllChem.ComputeMolVolume, "molecular_volume"),
                (rdMolDescriptors.CalcAsphericity, "asphericity"),
                (rdMolDescriptors.CalcEccentricity, "eccentricity"),
                (rdMolDescriptors.CalcInertialShapeFactor, "inertial_shape_factor"),
                (rdMolDescriptors.CalcRadiusOfGyration, "radius_of_gyration"),
                (rdMolDescriptors.CalcSpherocityIndex, "spherocity"),
            ],
            "requires_3d": True,
        },
        "fragment": {
            "description": "Molecular fragments and scaffolds",
            "functions": [
                (rdMolDescriptors.CalcNumSpiroAtoms, "num_spiro_atoms"),
                (rdMolDescriptors.CalcNumBridgeheadAtoms, "num_bridgehead_atoms"),
                (Scaffolds.MurckoScaffold.GetScaffoldForMol, "scaffold"),
            ],
        },
    }

    def __init__(
        self,
        descriptor_types: Optional[List[str]] = None,
        include_3d: bool = False,
        normalize: bool = True,
    ):
        """Initialize descriptor generator.

        Args:
            descriptor_types: Types of descriptors to generate
            include_3d: Whether to include 3D descriptors
            normalize: Whether to normalize descriptor values
        """
        super().__init__()
        self.descriptor_types = descriptor_types or [
            "constitutional",
            "topological",
            "electronic",
        ]
        if include_3d:
            self.descriptor_types.append("geometric")
        self.descriptor_types.append("fragment")

        self.include_3d = include_3d
        self.normalize = normalize

        # Initialize descriptor functions
        self._setup_descriptors()

    def _setup_descriptors(self):
        """Set up descriptor calculation functions."""
        self.descriptors = {}
        for desc_type in self.descriptor_types:
            if desc_type in self.DESCRIPTOR_TYPES:
                if not self.include_3d and self.DESCRIPTOR_TYPES[desc_type].get("requires_3d", False):
                    continue
                self.descriptors[desc_type] = self.DESCRIPTOR_TYPES[desc_type]["functions"]

    def generate(
        self,
        mols: List[Chem.Mol],
        descriptor_types: Optional[List[str]] = None,
    ) -> Dict[str, np.ndarray]:
        """Generate descriptors for molecules.

        Args:
            mols: List of RDKit molecules
            descriptor_types: Types of descriptors to generate

        Returns:
            Dictionary mapping descriptor types to feature arrays
        """
        if descriptor_types is None:
            descriptor_types = self.descriptor_types

        descriptors = {}
        try:
            for desc_type in descriptor_types:
                if desc_type in self.descriptors:
                    desc_values = []
                    for mol in mols:
                        try:
                            values = []
                            for func, _ in self.descriptors[desc_type]:
                                try:
                                    value = func(mol)
                                    if isinstance(value, (int, float)):
                                        values.append(float(value))
                                    else:
                                        values.append(0.0)
                                except Exception as e:
                                    self.logger.debug(f"Error calculating descriptor {func.__name__}: {str(e)}")
                                    values.append(0.0)
                            desc_values.append(values)
                        except Exception as e:
                            self.logger.debug(f"Error processing molecule for {desc_type}: {str(e)}")
                            desc_values.append([0.0] * len(self.descriptors[desc_type]))
                    descriptors[desc_type] = np.array(desc_values)

            return descriptors

        except Exception as e:
            self.logger.error(f"Error generating descriptors: {str(e)}")
            return {}

    def get_info(self) -> Dict[str, Dict]:
        """Get information about available descriptors."""
        info = {}
        for desc_type in self.descriptor_types:
            if desc_type in self.DESCRIPTOR_TYPES:
                info[desc_type] = {
                    "description": self.DESCRIPTOR_TYPES[desc_type]["description"],
                    "names": [name for _, name in self.DESCRIPTOR_TYPES[desc_type]["functions"]],
                    "requires_3d": self.DESCRIPTOR_TYPES[desc_type].get("requires_3d", False),
                }
        return info

    def get_descriptor_names(self, desc_type: str) -> List[str]:
        """Get names of descriptors for a given type."""
        if desc_type in self.DESCRIPTOR_TYPES:
            return [name for _, name in self.DESCRIPTOR_TYPES[desc_type]["functions"]]
        return []
