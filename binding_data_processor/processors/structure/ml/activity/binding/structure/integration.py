"""Integration module for coordinating structure analysis components.

This module provides:
1. Coordination between analysis components
2. Integration with external tools and services
3. Data sharing and transformation between components
4. Error handling and validation
"""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union, Any

import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

from .analyzer import StructureAnalyzer
from .geometry import GeometryCalculator
from .conservation import ConservationAnalyzer
from .dynamics import DynamicsAnalyzer
from .surface import SurfaceAnalyzer
from .pockets import PocketDetector
from ..pharmacophore import PharmacophoreGenerator
from ..alphafold import AlphaFoldIntegrator
from ..constants import BINDING_SITE_PARAMS

logger = logging.getLogger(__name__)


class StructureIntegrator:
    """Coordinates interactions between structure analysis components."""

    def __init__(
        self,
        analyzer: Optional[StructureAnalyzer] = None,
        geometry: Optional[GeometryCalculator] = None,
        conservation: Optional[ConservationAnalyzer] = None,
        dynamics: Optional[DynamicsAnalyzer] = None,
        surface: Optional[SurfaceAnalyzer] = None,
        pockets: Optional[PocketDetector] = None,
        pharmacophore: Optional[PharmacophoreGenerator] = None,
        alphafold: Optional[AlphaFoldIntegrator] = None,
    ):
        """Initialize structure integrator.

        Args:
            analyzer: Structure analyzer instance
            geometry: Geometry calculator instance
            conservation: Conservation analyzer instance
            dynamics: Dynamics analyzer instance
            surface: Surface analyzer instance
            pockets: Pocket detector instance
            pharmacophore: Pharmacophore generator instance
            alphafold: AlphaFold integrator instance
        """
        # Core components
        self.analyzer = analyzer or StructureAnalyzer()
        self.geometry = geometry or GeometryCalculator()
        self.surface = surface or SurfaceAnalyzer()
        self.pockets = pockets or PocketDetector()

        # Optional components
        self.conservation = conservation
        self.dynamics = dynamics
        self.pharmacophore = pharmacophore
        self.alphafold = alphafold

        # Cache for intermediate results
        self._cache = {}

    def analyze_structure(
        self,
        structure: Structure,
        include_conservation: bool = True,
        include_dynamics: bool = True,
        include_pharmacophore: bool = True,
        cache_results: bool = True,
    ) -> Dict[str, Any]:
        """Coordinate comprehensive structure analysis.

        Args:
            structure: BioPython Structure object
            include_conservation: Whether to include conservation analysis
            include_dynamics: Whether to include dynamics analysis
            include_pharmacophore: Whether to include pharmacophore features
            cache_results: Whether to cache intermediate results

        Returns:
            Dictionary of structure properties
        """
        try:
            # Clear cache if not reusing
            if not cache_results:
                self._cache.clear()

            # Basic analysis
            properties = self._analyze_basic_properties(structure)

            # Optional analyses
            if include_conservation and self.conservation:
                properties.update(self._analyze_conservation(structure))

            if include_dynamics and self.dynamics:
                properties.update(self._analyze_dynamics(structure))

            if include_pharmacophore and self.pharmacophore:
                properties.update(self._analyze_pharmacophore(structure))

            # Cache results if requested
            if cache_results:
                self._cache["structure_properties"] = properties

            return properties

        except Exception as e:
            logger.error(f"Error in structure analysis: {str(e)}")
            return {}

    def find_binding_sites(
        self,
        structure: Structure,
        properties: Optional[Dict[str, Any]] = None,
        min_volume: float = 100.0,
        max_sites: int = 5,
        include_scores: bool = True,
        include_pharmacophore: bool = True,
        reuse_cache: bool = True,
    ) -> List[Dict[str, Any]]:
        """Coordinate binding site detection and analysis.

        Args:
            structure: BioPython Structure object
            properties: Pre-calculated structure properties
            min_volume: Minimum pocket volume in Å³
            max_sites: Maximum number of sites to return
            include_scores: Whether to include detailed scoring
            include_pharmacophore: Whether to include pharmacophore analysis
            reuse_cache: Whether to reuse cached results

        Returns:
            List of dictionaries containing binding site properties
        """
        try:
            # Get or calculate properties
            if properties is None:
                if reuse_cache and "structure_properties" in self._cache:
                    properties = self._cache["structure_properties"]
                else:
                    properties = self.analyze_structure(structure)

            # Find and analyze sites
            sites = self._find_and_analyze_sites(
                structure,
                properties,
                min_volume=min_volume,
                max_sites=max_sites,
                include_scores=include_scores,
                include_pharmacophore=include_pharmacophore,
            )

            return sites

        except Exception as e:
            logger.error(f"Error finding binding sites: {str(e)}")
            return []

    def predict_structure(
        self,
        sequence: str,
        **kwargs,
    ) -> Optional[Structure]:
        """Coordinate structure prediction with AlphaFold.

        Args:
            sequence: Protein sequence
            **kwargs: Additional arguments passed to AlphaFold

        Returns:
            Predicted structure or None if prediction fails
        """
        try:
            if not self.alphafold:
                logger.warning("AlphaFold integration not available")
                return None

            structure = self.alphafold.predict_structure(sequence, **kwargs)
            if structure is not None:
                # Pre-calculate common properties
                properties = self.analyze_structure(structure)
                self._cache["structure_properties"] = properties

            return structure

        except Exception as e:
            logger.error(f"Error predicting structure: {str(e)}")
            return None

    def _analyze_basic_properties(self, structure: Structure) -> Dict[str, Any]:
        """Coordinate basic property analysis.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of basic properties
        """
        try:
            properties = {}

            # Geometric properties
            geom_props = self.geometry.get_properties(structure)
            properties["geometry"] = geom_props

            # Surface properties
            surface_props = self.surface.analyze_surface(structure)
            properties["surface"] = surface_props

            # Residue properties
            res_props = self._analyze_residues(structure)
            properties["residues"] = res_props

            # Chain properties
            chain_props = self._analyze_chains(structure)
            properties["chains"] = chain_props

            return properties

        except Exception as e:
            logger.error(f"Error analyzing basic properties: {str(e)}")
            return {}

    def _analyze_conservation(self, structure: Structure) -> Dict[str, Any]:
        """Coordinate conservation analysis.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of conservation properties
        """
        try:
            if not self.conservation:
                return {}

            conservation = self.conservation.analyze_conservation(structure)
            return {"conservation": conservation}

        except Exception as e:
            logger.error(f"Error analyzing conservation: {str(e)}")
            return {}

    def _analyze_dynamics(self, structure: Structure) -> Dict[str, Any]:
        """Coordinate dynamics analysis.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of dynamics properties
        """
        try:
            if not self.dynamics:
                return {}

            dynamics = self.dynamics.analyze_dynamics(structure)
            return {"dynamics": dynamics}

        except Exception as e:
            logger.error(f"Error analyzing dynamics: {str(e)}")
            return {}

    def _analyze_pharmacophore(self, structure: Structure) -> Dict[str, Any]:
        """Coordinate pharmacophore analysis.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of pharmacophore properties
        """
        try:
            if not self.pharmacophore:
                return {}

            features = self.pharmacophore.get_structure_features(structure)
            return {"pharmacophore": features}

        except Exception as e:
            logger.error(f"Error analyzing pharmacophore: {str(e)}")
            return {}

    def _analyze_residues(self, structure: Structure) -> Dict[int, Dict[str, Any]]:
        """Coordinate residue-level analysis.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary mapping residue numbers to properties
        """
        try:
            properties = {}

            for residue in structure.get_residues():
                res_id = residue.get_id()[1]
                properties[res_id] = {
                    "name": residue.get_resname(),
                    "chain": residue.get_parent().get_id(),
                    "center": self.geometry.get_center(residue),
                    "surface_area": self.surface.get_residue_area(residue),
                    "is_surface": self.surface.is_surface_residue(residue),
                }

                # Add conservation if available
                if self.conservation:
                    cons_score = self.conservation.get_residue_conservation(residue)
                    properties[res_id]["conservation"] = cons_score

                # Add dynamics if available
                if self.dynamics:
                    dyn_props = self.dynamics.get_residue_dynamics(residue)
                    properties[res_id]["dynamics"] = dyn_props

            return properties

        except Exception as e:
            logger.error(f"Error analyzing residues: {str(e)}")
            return {}

    def _analyze_chains(self, structure: Structure) -> Dict[str, Dict[str, Any]]:
        """Coordinate chain-level analysis.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary mapping chain IDs to properties
        """
        try:
            properties = {}

            for chain in structure.get_chains():
                chain_id = chain.get_id()
                properties[chain_id] = {
                    "sequence": self._get_sequence(chain),
                    "num_residues": len(list(chain.get_residues())),
                    "center": self.geometry.get_chain_center(chain),
                }

                # Add conservation if available
                if self.conservation:
                    cons_profile = self.conservation.get_chain_profile(chain)
                    properties[chain_id]["conservation_profile"] = cons_profile

                # Add dynamics if available
                if self.dynamics:
                    dyn_props = self.dynamics.get_chain_dynamics(chain)
                    properties[chain_id]["dynamics"] = dyn_props

            return properties

        except Exception as e:
            logger.error(f"Error analyzing chains: {str(e)}")
            return {}

    def _find_and_analyze_sites(
        self,
        structure: Structure,
        properties: Dict[str, Any],
        min_volume: float = 100.0,
        max_sites: int = 5,
        include_scores: bool = True,
        include_pharmacophore: bool = True,
    ) -> List[Dict[str, Any]]:
        """Coordinate binding site detection and analysis.

        Args:
            structure: BioPython Structure object
            properties: Structure properties
            min_volume: Minimum pocket volume
            max_sites: Maximum number of sites
            include_scores: Whether to include scores
            include_pharmacophore: Whether to include pharmacophore analysis

        Returns:
            List of binding site properties
        """
        try:
            # Find potential sites
            sites = self.pockets.find_pockets(
                structure,
                min_volume=min_volume,
                properties=properties,
            )

            # Score and rank sites
            scored_sites = self._score_sites(
                sites,
                structure,
                properties,
                include_scores=include_scores,
            )

            # Get top sites
            top_sites = sorted(
                scored_sites,
                key=lambda x: x["score"],
                reverse=True,
            )[:max_sites]

            # Analyze top sites
            analyzed_sites = []
            for site in top_sites:
                analysis = self._analyze_site(
                    structure,
                    site,
                    properties,
                    include_pharmacophore=include_pharmacophore,
                )
                analyzed_sites.append({**site, **analysis})

            return analyzed_sites

        except Exception as e:
            logger.error(f"Error finding and analyzing sites: {str(e)}")
            return []

    def _score_sites(
        self,
        sites: List[Dict[str, Any]],
        structure: Structure,
        properties: Dict[str, Any],
        include_scores: bool = True,
    ) -> List[Dict[str, Any]]:
        """Coordinate binding site scoring.

        Args:
            sites: List of binding sites
            structure: BioPython Structure object
            properties: Structure properties
            include_scores: Whether to include detailed scores

        Returns:
            List of scored binding sites
        """
        try:
            scored_sites = []
            for site in sites:
                scores = {}

                # Geometric score
                geom_score = self.geometry.score_site(site, structure)
                scores["geometry"] = float(geom_score)

                # Conservation score if available
                if self.conservation:
                    cons_score = self.conservation.score_site(
                        site["residues"],
                        properties.get("conservation", {}),
                    )
                    scores["conservation"] = float(cons_score)

                # Dynamics score if available
                if self.dynamics:
                    dyn_score = self.dynamics.score_site(
                        site["residues"],
                        properties.get("dynamics", {}),
                    )
                    scores["dynamics"] = float(dyn_score)

                # Pharmacophore score if available
                if self.pharmacophore:
                    pharm_score = self.pharmacophore.score_site(
                        structure,
                        site["residues"],
                    )
                    scores["pharmacophore"] = float(pharm_score)

                # Calculate total score
                weights = {
                    "geometry": 0.4,
                    "conservation": 0.2,
                    "dynamics": 0.2,
                    "pharmacophore": 0.2,
                }
                total = 0.0
                weight_sum = 0.0

                for key, weight in weights.items():
                    if key in scores:
                        total += weight * scores[key]
                        weight_sum += weight

                if weight_sum > 0:
                    total /= weight_sum

                # Add scores to site
                site_copy = site.copy()
                site_copy["score"] = float(total)
                if include_scores:
                    site_copy["score_components"] = scores

                scored_sites.append(site_copy)

            return scored_sites

        except Exception as e:
            logger.error(f"Error scoring sites: {str(e)}")
            return []

    def _analyze_site(
        self,
        structure: Structure,
        site: Dict[str, Any],
        properties: Dict[str, Any],
        include_pharmacophore: bool = True,
    ) -> Dict[str, Any]:
        """Coordinate detailed binding site analysis.

        Args:
            structure: BioPython Structure object
            site: Binding site information
            properties: Structure properties
            include_pharmacophore: Whether to include pharmacophore analysis

        Returns:
            Dictionary of binding site properties
        """
        try:
            analysis = {}

            # Basic properties
            analysis.update(
                {
                    "volume": self.geometry.calculate_volume(site["residues"]),
                    "surface_area": self.surface.calculate_surface_area(site["residues"]),
                    "depth": self.surface.calculate_depth(site["residues"]),
                }
            )

            # Residue properties
            res_props = self._analyze_site_residues(
                structure,
                site["residues"],
                properties,
            )
            analysis["residue_properties"] = res_props

            # Pharmacophore analysis
            if include_pharmacophore and self.pharmacophore:
                pharm_props = self.pharmacophore.analyze_site(
                    structure,
                    site["residues"],
                )
                analysis["pharmacophore"] = pharm_props

            return analysis

        except Exception as e:
            logger.error(f"Error analyzing site: {str(e)}")
            return {}

    def _analyze_site_residues(
        self,
        structure: Structure,
        residue_ids: List[int],
        properties: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Analyze binding site residues.

        Args:
            structure: BioPython Structure object
            residue_ids: List of residue IDs
            properties: Structure properties

        Returns:
            Dictionary of residue properties
        """
        try:
            analysis = {}

            # Get residues
            residues = []
            for res_id in residue_ids:
                for residue in structure.get_residues():
                    if residue.get_id()[1] == res_id:
                        residues.append(residue)
                        break

            # Analyze composition
            composition = {}
            for residue in residues:
                res_name = residue.get_resname()
                composition[res_name] = composition.get(res_name, 0) + 1

            # Convert to percentages
            total = len(residues)
            if total > 0:
                composition = {k: (v / total) * 100 for k, v in composition.items()}

            analysis["composition"] = composition

            # Add conservation if available
            if self.conservation:
                cons_stats = self._get_conservation_stats(residue_ids, properties)
                analysis["conservation"] = cons_stats

            # Add dynamics if available
            if self.dynamics:
                dyn_stats = self._get_dynamics_stats(residue_ids, properties)
                analysis["dynamics"] = dyn_stats

            return analysis

        except Exception as e:
            logger.error(f"Error analyzing site residues: {str(e)}")
            return {}

    def _get_sequence(self, chain: Chain) -> str:
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
            logger.error(f"Error getting sequence: {str(e)}")
            return ""

    def _get_conservation_stats(
        self,
        residue_ids: List[int],
        properties: Dict[str, Any],
    ) -> Dict[str, float]:
        """Get conservation statistics.

        Args:
            residue_ids: List of residue IDs
            properties: Structure properties

        Returns:
            Dictionary of conservation statistics
        """
        try:
            if not self.conservation:
                return {}

            conservation = properties.get("conservation", {})
            scores = [conservation.get(res_id, 0.0) for res_id in residue_ids]

            if not scores:
                return {}

            return {
                "mean": float(np.mean(scores)),
                "std": float(np.std(scores)),
                "min": float(np.min(scores)),
                "max": float(np.max(scores)),
            }

        except Exception as e:
            logger.error(f"Error getting conservation stats: {str(e)}")
            return {}

    def _get_dynamics_stats(
        self,
        residue_ids: List[int],
        properties: Dict[str, Any],
    ) -> Dict[str, float]:
        """Get dynamics statistics.

        Args:
            residue_ids: List of residue IDs
            properties: Structure properties

        Returns:
            Dictionary of dynamics statistics
        """
        try:
            if not self.dynamics:
                return {}

            dynamics = properties.get("dynamics", {})
            flexibility = [dynamics.get(res_id, {}).get("flexibility", 0.0) for res_id in residue_ids]

            if not flexibility:
                return {}

            return {
                "mean_flexibility": float(np.mean(flexibility)),
                "std_flexibility": float(np.std(flexibility)),
                "min_flexibility": float(np.min(flexibility)),
                "max_flexibility": float(np.max(flexibility)),
            }

        except Exception as e:
            logger.error(f"Error getting dynamics stats: {str(e)}")
            return {}
