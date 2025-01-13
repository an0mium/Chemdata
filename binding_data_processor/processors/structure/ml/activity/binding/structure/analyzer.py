"""Enhanced structure analyzer with integrated analysis capabilities.

This module provides a comprehensive structure analysis framework that:
1. Coordinates specialized analysis components
2. Integrates with binding site detection
3. Integrates with pharmacophore analysis
4. Integrates with AlphaFold predictions
5. Provides ML-based scoring and analysis
"""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union, Any

import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.DSSP import dssp_dict_from_pdb_file, residue_max_acc

from .geometry import GeometryCalculator
from .conservation import ConservationAnalyzer
from .dynamics import DynamicsAnalyzer
from .surface import SurfaceAnalyzer
from .pockets import PocketDetector
from ....pharmacophore.conformers import ConformerGenerator
from ..alphafold import AlphaFoldIntegrator
from ..constants import BINDING_SITE_PARAMS

logger = logging.getLogger(__name__)


class StructureAnalyzer:
    """Enhanced structure analyzer coordinating all analysis components."""

    def __init__(
        self,
        use_alphafold: bool = True,
        use_ml: bool = True,
        use_conservation: bool = True,
        use_dynamics: bool = True,
        use_shape_analysis: bool = True,
        use_pharmacophore: bool = True,
        **kwargs,
    ):
        """Initialize structure analyzer with all components.

        Args:
            use_alphafold: Whether to use AlphaFold integration
            use_ml: Whether to use machine learning models
            use_conservation: Whether to include conservation analysis
            use_dynamics: Whether to include dynamics analysis
            use_shape_analysis: Whether to include shape analysis
            use_pharmacophore: Whether to include pharmacophore detection
            **kwargs: Additional configuration parameters
        """
        self.version = "0.4.0"
        self.params = {**BINDING_SITE_PARAMS, **kwargs}

        # Initialize core analyzers
        self.geometry = GeometryCalculator()
        self.surface = SurfaceAnalyzer()
        self.pockets = PocketDetector(
            use_ml=use_ml,
            use_shape_analysis=use_shape_analysis,
            use_pharmacophore=use_pharmacophore,
        )

        # Optional analyzers based on features
        self.conservation = ConservationAnalyzer() if use_conservation else None
        self.dynamics = DynamicsAnalyzer() if use_dynamics else None

        # Optional integrations
        self.pharmacophore = ConformerGenerator() if use_pharmacophore else None
        self.alphafold = AlphaFoldIntegrator() if use_alphafold else None

        # Store settings
        self.settings = {
            "use_alphafold": use_alphafold,
            "use_ml": use_ml,
            "use_conservation": use_conservation,
            "use_dynamics": use_dynamics,
            "use_shape_analysis": use_shape_analysis,
            "use_pharmacophore": use_pharmacophore,
        }

    def get_structure_properties(
        self,
        structure: Structure,
        include_dssp: bool = True,
        include_conservation: bool = True,
        include_dynamics: bool = True,
        include_pharmacophore: bool = True,
    ) -> Dict[str, Any]:
        """Get comprehensive structure properties.

        Args:
            structure: BioPython Structure object
            include_dssp: Whether to include secondary structure
            include_conservation: Whether to include conservation analysis
            include_dynamics: Whether to include dynamics analysis
            include_pharmacophore: Whether to include pharmacophore features

        Returns:
            Dictionary of structure properties
        """
        try:
            # Get basic geometric properties
            properties = self.geometry.get_properties(structure)

            # Add surface properties
            surface_props = self.surface.analyze_surface(structure)
            properties.update(surface_props)

            # Secondary structure if requested
            if include_dssp:
                ss_props = self._get_secondary_structure(structure)
                properties["secondary_structure"] = ss_props

            # Optional analyses based on settings
            if include_conservation and self.conservation:
                conservation = self.conservation.analyze_conservation(structure)
                properties["conservation"] = conservation

            if include_dynamics and self.dynamics:
                dynamics = self.dynamics.analyze_dynamics(structure)
                properties["dynamics"] = dynamics

            if include_pharmacophore and self.pharmacophore:
                pharm_features = self.pharmacophore.get_structure_features(structure)
                properties["pharmacophore"] = pharm_features

            # Enhanced properties
            properties.update(
                {
                    "residues": self._get_residue_properties(structure),
                    "chains": self._get_chain_properties(structure),
                    "packing_density": self.geometry.calculate_packing_density(structure),
                    "electrostatics": self.geometry.calculate_electrostatics(structure),
                    "shape_descriptors": self.geometry.calculate_shape_descriptors(structure),
                    "residue_composition": self._get_residue_composition(structure),
                    "surface_patches": {
                        "hydrophobic": self._get_hydrophobic_patches(structure),
                        "charged": self._get_charged_patches(structure),
                    },
                }
            )

            return properties

        except Exception as e:
            logger.error(f"Error getting structure properties: {str(e)}")
            return {}

    def find_binding_sites(
        self,
        structure: Structure,
        min_volume: float = 100.0,
        max_sites: int = 5,
        include_scores: bool = True,
        include_pharmacophore: bool = True,
    ) -> List[Dict[str, Any]]:
        """Find potential binding sites in the structure.

        Args:
            structure: BioPython Structure object
            min_volume: Minimum pocket volume in Å³
            max_sites: Maximum number of sites to return
            include_scores: Whether to include detailed scoring
            include_pharmacophore: Whether to include pharmacophore analysis

        Returns:
            List of dictionaries containing binding site properties
        """
        try:
            # Get structure properties first
            properties = self.get_structure_properties(structure)

            # Find potential binding sites
            sites = self.pockets.find_pockets(
                structure,
                min_volume=min_volume,
                properties=properties,
            )

            # Score and rank sites
            scored_sites = self.pockets.score_pockets(
                sites,
                structure,
                properties,
                include_scores=include_scores,
            )

            # Get top N sites
            top_sites = sorted(
                scored_sites,
                key=lambda x: x["score"],
                reverse=True,
            )[:max_sites]

            # Add detailed analysis for top sites
            for site in top_sites:
                site.update(
                    self._analyze_binding_site(
                        structure,
                        site,
                        properties,
                        include_pharmacophore=include_pharmacophore,
                    )
                )

            return top_sites

        except Exception as e:
            logger.error(f"Error finding binding sites: {str(e)}")
            return []

    def predict_structure(
        self,
        sequence: str,
        **kwargs,
    ) -> Optional[Structure]:
        """Predict protein structure using AlphaFold.

        Args:
            sequence: Protein sequence
            **kwargs: Additional arguments passed to AlphaFold

        Returns:
            Predicted structure or None if prediction fails
        """
        try:
            if not self.alphafold:
                logger.warning("AlphaFold integration not enabled")
                return None

            return self.alphafold.predict_structure(sequence, **kwargs)

        except Exception as e:
            logger.error(f"Error predicting structure: {str(e)}")
            return None

    def _analyze_binding_site(
        self,
        structure: Structure,
        site: Dict[str, Any],
        properties: Dict[str, Any],
        include_pharmacophore: bool = True,
    ) -> Dict[str, Any]:
        """Perform detailed analysis of a binding site.

        Args:
            structure: BioPython Structure object
            site: Dictionary containing binding site info
            properties: Pre-calculated structure properties
            include_pharmacophore: Whether to include pharmacophore analysis

        Returns:
            Dictionary of binding site properties
        """
        try:
            analysis = {}

            # Get geometric properties
            geometry = self.geometry.analyze_site(structure, site["residues"])
            analysis["geometry"] = geometry

            # Get surface properties
            surface = self.surface.analyze_site(structure, site["residues"])
            analysis["surface"] = surface

            # Get conservation if available
            if self.conservation and "conservation" in properties:
                conservation = self.conservation.analyze_site(
                    site["residues"],
                    properties["conservation"],
                )
                analysis["conservation"] = conservation

            # Get dynamics if available
            if self.dynamics and "dynamics" in properties:
                dynamics = self.dynamics.analyze_site(
                    site["residues"],
                    properties["dynamics"],
                )
                analysis["dynamics"] = dynamics

            # Enhanced analysis
            analysis.update(
                {
                    "volume": self.geometry.calculate_volume(site["residues"]),
                    "surface_area": self.surface.calculate_surface_area(site["residues"]),
                    "depth": self.surface.calculate_depth(site["residues"]),
                    "hydrophobicity": self.geometry.calculate_hydrophobicity(site["residues"]),
                    "charge": self.geometry.calculate_charge(site["residues"]),
                    "shape": self.geometry.calculate_shape_descriptors(site["residues"]),
                    "electrostatics": self.geometry.calculate_electrostatics(site["residues"]),
                    "residue_composition": self._get_residue_composition(site["residues"]),
                    "secondary_structure": self._analyze_secondary_structure(site["residues"]),
                }
            )

            # Optional pharmacophore analysis
            if include_pharmacophore and self.pharmacophore:
                analysis["pharmacophore"] = self.pharmacophore.get_site_features(
                    structure,
                    site["residues"],
                )

            return analysis

        except Exception as e:
            logger.error(f"Error analyzing binding site: {str(e)}")
            return {}

    def _get_residue_properties(self, structure: Structure) -> Dict[str, Any]:
        """Get residue-level properties.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of residue properties
        """
        try:
            properties = {}

            for residue in structure.get_residues():
                res_id = residue.get_id()[1]
                properties[res_id] = {
                    "name": residue.get_resname(),
                    "chain": residue.get_parent().get_id(),
                    "center": self.geometry.get_center(residue),
                    "num_atoms": len(list(residue.get_atoms())),
                    "surface_area": self.surface.get_residue_area(residue),
                    "is_surface": self.surface.is_surface_residue(residue),
                }

            return properties

        except Exception as e:
            logger.error(f"Error getting residue properties: {str(e)}")
            return {}

    def _get_chain_properties(self, structure: Structure) -> Dict[str, Any]:
        """Get chain-level properties.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of chain properties
        """
        try:
            properties = {}

            for chain in structure.get_chains():
                chain_id = chain.get_id()
                residues = list(chain.get_residues())

                properties[chain_id] = {
                    "num_residues": len(residues),
                    "sequence": self._get_chain_sequence(chain),
                    "center": self.geometry.get_chain_center(chain),
                    "secondary_structure": self._get_chain_ss(chain),
                }

            return properties

        except Exception as e:
            logger.error(f"Error getting chain properties: {str(e)}")
            return {}

    def _get_chain_sequence(self, chain: Chain) -> str:
        """Get amino acid sequence.

        Args:
            chain: BioPython Chain object

        Returns:
            Amino acid sequence string
        """
        try:
            three_to_one = {
                "ALA": "A",
                "CYS": "C",
                "ASP": "D",
                "GLU": "E",
                "PHE": "F",
                "GLY": "G",
                "HIS": "H",
                "ILE": "I",
                "LYS": "K",
                "LEU": "L",
                "MET": "M",
                "ASN": "N",
                "PRO": "P",
                "GLN": "Q",
                "ARG": "R",
                "SER": "S",
                "THR": "T",
                "VAL": "V",
                "TRP": "W",
                "TYR": "Y",
            }

            sequence = ""
            for residue in chain:
                if residue.get_id()[0] == " ":  # Standard amino acid
                    sequence += three_to_one.get(residue.get_resname(), "X")

            return sequence

        except Exception as e:
            logger.error(f"Error getting chain sequence: {str(e)}")
            return ""

    def _get_secondary_structure(self, structure: Structure) -> Dict[str, Any]:
        """Get secondary structure information.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of secondary structure properties
        """
        try:
            # Get DSSP data
            dssp_dict = dssp_dict_from_pdb_file(structure)[0]

            # Calculate statistics
            ss_stats = {"H": 0, "E": 0, "C": 0}  # Helix, Sheet, Coil
            accessibility = {}
            phi_psi = {}

            for residue in structure.get_residues():
                res_id = residue.get_id()
                chain_id = residue.get_parent().get_id()
                key = (chain_id, res_id)

                if key in dssp_dict:
                    dssp_data = dssp_dict[key]

                    # Secondary structure
                    ss = dssp_data[2]
                    if ss in ["H", "G", "I"]:  # Various helix types
                        ss_stats["H"] += 1
                    elif ss in ["E", "B"]:  # Various sheet types
                        ss_stats["E"] += 1
                    else:
                        ss_stats["C"] += 1

                    # Accessibility
                    acc = float(dssp_data[3])
                    max_acc = residue_max_acc[residue.get_resname()]
                    rel_acc = acc / max_acc if max_acc > 0 else 0
                    accessibility[res_id[1]] = rel_acc

                    # Phi/Psi angles
                    phi_psi[res_id[1]] = (float(dssp_data[4]), float(dssp_data[5]))

            # Calculate percentages
            total = sum(ss_stats.values())
            if total > 0:
                ss_stats = {k: v / total for k, v in ss_stats.items()}

            return {
                "secondary_structure": ss_stats,
                "accessibility": accessibility,
                "phi_psi": phi_psi,
            }

        except Exception as e:
            logger.error(f"Error getting secondary structure: {str(e)}")
            return {}

    def _get_residue_composition(self, residues: List[Residue]) -> Dict[str, float]:
        """Get amino acid composition.

        Args:
            residues: List of residues

        Returns:
            Dictionary mapping residue types to percentages
        """
        try:
            composition = {}
            total = len(residues)

            for res in residues:
                res_name = res.get_resname()
                composition[res_name] = composition.get(res_name, 0) + 1

            # Convert to percentages
            if total > 0:
                composition = {k: (v / total) * 100 for k, v in composition.items()}

            return composition

        except Exception as e:
            logger.error(f"Error getting residue composition: {str(e)}")
            return {}

    def _analyze_secondary_structure(self, residues: List[Residue]) -> Dict[str, float]:
        """Analyze secondary structure composition.

        Args:
            residues: List of residues

        Returns:
            Dictionary of secondary structure statistics
        """
        try:
            ss_counts = {"H": 0, "E": 0, "C": 0}  # Helix, Sheet, Coil
            total = len(residues)

            for res in residues:
                ss = self.geometry.get_secondary_structure(res)
                ss_counts[ss] += 1

            # Convert to percentages
            if total > 0:
                ss_counts = {k: (v / total) * 100 for k, v in ss_counts.items()}

            return ss_counts

        except Exception as e:
            logger.error(f"Error analyzing secondary structure: {str(e)}")
            return {"H": 0.0, "E": 0.0, "C": 0.0}

    def _get_hydrophobic_patches(self, structure: Structure) -> List[Dict[str, Any]]:
        """Find hydrophobic surface patches.

        Args:
            structure: BioPython Structure object

        Returns:
            List of hydrophobic patch properties
        """
        try:
            return self.surface.find_hydrophobic_patches(structure)

        except Exception as e:
            logger.error(f"Error finding hydrophobic patches: {str(e)}")
            return []

    def _get_charged_patches(self, structure: Structure) -> List[Dict[str, Any]]:
        """Find charged surface patches.

        Args:
            structure: BioPython Structure object

        Returns:
            List of charged patch properties
        """
        try:
            return self.surface.find_charged_patches(structure)

        except Exception as e:
            logger.error(f"Error finding charged patches: {str(e)}")
            return []

    def _get_chain_ss(self, chain: Chain) -> Dict[str, float]:
        """Get chain secondary structure composition.

        Args:
            chain: BioPython Chain object

        Returns:
            Dictionary of secondary structure percentages
        """
        try:
            residues = list(chain.get_residues())
            return self._analyze_secondary_structure(residues)

        except Exception as e:
            logger.error(f"Error getting chain secondary structure: {str(e)}")
            return {"H": 0.0, "E": 0.0, "C": 0.0}
