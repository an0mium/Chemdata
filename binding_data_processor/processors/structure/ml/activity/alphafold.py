"""AlphaFold integration for protein structure prediction and analysis."""

import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
import json
import requests
import numpy as np
from Bio import PDB, SeqIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial import ConvexHull

logger = logging.getLogger(__name__)


@dataclass
class AlphaFoldConfig:
    """Configuration for AlphaFold integration."""

    # Model paths and parameters
    model_preset: str = "monomer"  # monomer, monomer_casp14, monomer_ptm, multimer
    model_directory: Optional[str] = None
    max_template_date: str = "2022-01-01"
    db_preset: str = "full_dbs"

    # Prediction parameters
    num_ensemble: int = 1
    num_recycle: int = 3
    recycle_early_stop_tolerance: float = 0.5
    num_seeds: int = 1

    # Resource constraints
    max_memory_gb: float = 14.0
    max_gpu_memory_gb: float = 13.0

    # Output configuration
    output_dir: Optional[str] = None
    save_all: bool = False
    save_intermediates: bool = False

    # API configuration
    use_api: bool = True
    api_base_url: str = "https://alphafold.ebi.ac.uk/api"
    api_timeout: int = 300


class AlphaFoldResult:
    """Container for AlphaFold prediction results."""

    def __init__(
        self,
        structure: Structure,
        confidence: Dict[str, float],
        plddt_scores: Dict[str, float],
        pae_matrix: Optional[np.ndarray] = None,
        intermediates: Optional[Dict[str, Any]] = None,
    ):
        """Initialize AlphaFold results.

        Args:
            structure: Predicted protein structure
            confidence: Overall confidence metrics
            plddt_scores: Per-residue confidence scores
            pae_matrix: Predicted aligned error matrix
            intermediates: Intermediate prediction results
        """
        self.structure = structure
        self.confidence = confidence
        self.plddt_scores = plddt_scores
        self.pae_matrix = pae_matrix
        self.intermediates = intermediates or {}


class AlphaFoldPredictor:
    """Interface for AlphaFold protein structure prediction and analysis."""

    # AlphaFold confidence score thresholds
    CONFIDENCE_THRESHOLDS = {"high": 90.0, "medium": 70.0, "low": 50.0}

    def __init__(
        self,
        config: Optional[Union[Dict, AlphaFoldConfig]] = None,
        cache_dir: Optional[Union[str, Path]] = None,
    ):
        """Initialize AlphaFold predictor.

        Args:
            config: Predictor configuration
            cache_dir: Directory to cache predictions
        """
        if isinstance(config, dict):
            config = AlphaFoldConfig(**config)
        elif config is None:
            config = AlphaFoldConfig()

        self.config = config
        self.cache_dir = Path(cache_dir) if cache_dir else Path.home() / ".alphafold_cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        if self.config.output_dir:
            Path(self.config.output_dir).mkdir(parents=True, exist_ok=True)

        self.logger = logging.getLogger(self.__class__.__name__)
        self.pdb_parser = PDB.PDBParser(QUIET=True)
        self.structure_cache = {}

    async def predict(
        self,
        sequence: str,
        use_templates: bool = True,
        force_download: bool = False,
        **kwargs,
    ) -> Optional[AlphaFoldResult]:
        """Predict protein structure using AlphaFold.

        Args:
            sequence: Amino acid sequence
            use_templates: Whether to use templates
            force_download: Whether to force download even if cached
            **kwargs: Additional arguments passed to AlphaFold

        Returns:
            Prediction results or None if prediction fails
        """
        try:
            if self.config.use_api:
                return await self._predict_with_api(sequence, force_download)
            else:
                return await self._predict_with_local_model(sequence, use_templates, **kwargs)

        except Exception as e:
            self.logger.error(f"Error in AlphaFold prediction: {str(e)}")
            return None

    async def _predict_with_api(
        self,
        sequence: str,
        force_download: bool = False,
    ) -> Optional[AlphaFoldResult]:
        """Predict structure using AlphaFold API."""
        try:
            # Search AlphaFold database
            search_params = {"query": sequence}
            response = requests.get(
                f"{self.config.api_base_url}/prediction/search",
                params=search_params,
                timeout=self.config.api_timeout,
            )
            response.raise_for_status()
            results = response.json()

            if not results:
                self.logger.warning(f"No AlphaFold structure found for sequence")
                return None

            # Get best model by confidence score
            best_model = max(results, key=lambda x: x["confidence_score"])
            if best_model["confidence_score"] < self.CONFIDENCE_THRESHOLDS["medium"]:
                self.logger.warning(f"Best model confidence ({best_model['confidence_score']:.1f}) " f"below threshold ({self.CONFIDENCE_THRESHOLDS['medium']})")
                return None

            # Download structure
            download_params = {"id": best_model["id"]}
            response = requests.get(
                f"{self.config.api_base_url}/prediction/download",
                params=download_params,
                timeout=self.config.api_timeout,
            )
            response.raise_for_status()

            # Parse structure and create result
            structure = self.pdb_parser.get_structure("prediction", response.content)

            # Extract confidence scores
            confidence = {
                "global": best_model["confidence_score"],
                "local": self._extract_plddt_scores(structure),
            }

            result = AlphaFoldResult(
                structure=structure,
                confidence=confidence,
                plddt_scores=confidence["local"],
                pae_matrix=best_model.get("pae_matrix"),
            )

            return result

        except Exception as e:
            self.logger.error(f"Error in API prediction: {str(e)}")
            return None

    async def _predict_with_local_model(
        self,
        sequence: str,
        use_templates: bool = True,
        **kwargs,
    ) -> Optional[AlphaFoldResult]:
        """Predict structure using local AlphaFold installation."""
        try:
            # TODO: Implement local model prediction
            self.logger.warning("Local model prediction not yet implemented")
            return None

        except Exception as e:
            self.logger.error(f"Error in local model prediction: {str(e)}")
            return None

    async def analyze_structure(
        self,
        structure: Structure,
        include_binding_sites: bool = True,
        **kwargs,
    ) -> Dict[str, Any]:
        """Analyze protein structure comprehensively.

        Args:
            structure: BioPython Structure object
            include_binding_sites: Whether to analyze binding sites
            **kwargs: Additional analysis parameters

        Returns:
            Dictionary containing analysis results
        """
        try:
            analysis = {
                "basic_properties": self._analyze_basic_properties(structure),
                "secondary_structure": self._analyze_secondary_structure(structure),
                "surface_properties": self._analyze_surface_properties(structure),
                "stability": self._analyze_stability(structure),
            }

            if include_binding_sites:
                analysis["binding_sites"] = self._analyze_binding_sites(structure)

            return analysis

        except Exception as e:
            self.logger.error(f"Error analyzing structure: {str(e)}")
            return {}

    def _analyze_basic_properties(self, structure: Structure) -> Dict[str, Any]:
        """Analyze basic structural properties."""
        try:
            properties = {
                "num_residues": self._count_residues(structure),
                "num_atoms": self._count_atoms(structure),
                "molecular_weight": self._calculate_molecular_weight(structure),
                "radius_of_gyration": self._calculate_radius_of_gyration(structure),
                "residue_composition": self._get_residue_composition(structure),
            }
            return properties
        except Exception as e:
            self.logger.error(f"Error analyzing basic properties: {str(e)}")
            return {}

    def _analyze_binding_sites(
        self,
        structure: Structure,
        radius: float = 10.0,
    ) -> List[Dict[str, Any]]:
        """Analyze potential binding sites."""
        try:
            binding_sites = []

            # Get surface pockets
            pockets = self._find_surface_pockets(structure)

            for pocket in pockets:
                site = {
                    "residues": pocket["residues"],
                    "volume": self._calculate_volume(pocket["residues"]),
                    "surface_area": self._calculate_surface_area(pocket["residues"]),
                    "hydrophobicity": self._calculate_hydrophobicity(pocket["residues"]),
                    "charge": self._calculate_charge(pocket["residues"]),
                    "conservation": self._analyze_conservation(pocket["residues"]),
                    "flexibility": self._analyze_flexibility(pocket["residues"]),
                }
                binding_sites.append(site)

            return binding_sites

        except Exception as e:
            self.logger.error(f"Error analyzing binding sites: {str(e)}")
            return []

    def _analyze_surface_properties(self, structure: Structure) -> Dict[str, Any]:
        """Analyze protein surface properties."""
        try:
            properties = {
                "total_area": self._calculate_total_surface_area(structure),
                "hydrophobic_patches": self._find_hydrophobic_patches(structure),
                "electrostatic_potential": self._calculate_electrostatic_potential(structure),
                "surface_residues": self._identify_surface_residues(structure),
            }
            return properties
        except Exception as e:
            self.logger.error(f"Error analyzing surface properties: {str(e)}")
            return {}

    def _analyze_stability(self, structure: Structure) -> Dict[str, Any]:
        """Analyze protein stability indicators."""
        try:
            stability = {
                "hydrogen_bonds": self._analyze_hydrogen_bonds(structure),
                "salt_bridges": self._analyze_salt_bridges(structure),
                "disulfide_bonds": self._find_disulfide_bonds(structure),
                "buried_surface_area": self._calculate_buried_surface_area(structure),
                "packing_density": self._calculate_packing_density(structure),
            }
            return stability
        except Exception as e:
            self.logger.error(f"Error analyzing stability: {str(e)}")
            return {}

    def _extract_plddt_scores(self, structure: Structure) -> Dict[str, float]:
        """Extract per-residue pLDDT scores from structure."""
        scores = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    ca_atom = residue.get("CA")
                    if ca_atom:
                        scores[f"{chain.id}_{residue.id[1]}"] = ca_atom.get_bfactor()
        return scores

    # Helper methods from AlphaFoldIntegration
    def _calculate_volume(self, residues: List[PDB.Residue.Residue]) -> float:
        """Calculate volume of residue selection."""
        coords = []
        for res in residues:
            for atom in res:
                coords.append(atom.get_coord())
        coords = np.array(coords)
        hull = ConvexHull(coords)
        return hull.volume

    def _calculate_surface_area(self, residues: List[PDB.Residue.Residue]) -> float:
        """Calculate solvent accessible surface area."""
        coords = []
        radii = []
        for res in residues:
            for atom in res:
                coords.append(atom.get_coord())
                radii.append(self._get_atom_radius(atom))
        from Bio.PDB.SASA import ShrakeRupley

        sr = ShrakeRupley()
        return sr.compute(coords, radii)

    def _calculate_hydrophobicity(self, residues: List[PDB.Residue.Residue]) -> float:
        """Calculate average hydrophobicity of residues."""
        # Kyte-Doolittle hydrophobicity scale
        hydrophobicity_scale = {
            "ILE": 4.5,
            "VAL": 4.2,
            "LEU": 3.8,
            "PHE": 2.8,
            "CYS": 2.5,
            "MET": 1.9,
            "ALA": 1.8,
            "GLY": -0.4,
            "THR": -0.7,
            "SER": -0.8,
            "TRP": -0.9,
            "TYR": -1.3,
            "PRO": -1.6,
            "HIS": -3.2,
            "GLU": -3.5,
            "GLN": -3.5,
            "ASP": -3.5,
            "ASN": -3.5,
            "LYS": -3.9,
            "ARG": -4.5,
        }
        scores = [hydrophobicity_scale.get(res.resname, 0) for res in residues]
        return np.mean(scores)

    def _calculate_charge(self, residues: List[PDB.Residue.Residue]) -> float:
        """Calculate net charge of residues."""
        charge_scale = {"ARG": 1, "LYS": 1, "ASP": -1, "GLU": -1, "HIS": 0.5}
        charges = [charge_scale.get(res.resname, 0) for res in residues]
        return sum(charges)

    def _analyze_conservation(self, residues: List[PDB.Residue.Residue]) -> Dict[int, float]:
        """Analyze residue conservation scores."""
        # This would typically use pre-computed conservation scores
        # from multiple sequence alignments
        return {res.id[1]: 0.5 for res in residues}  # Placeholder

    def _analyze_flexibility(self, residues: List[PDB.Residue.Residue]) -> Dict[int, float]:
        """Analyze residue flexibility from B-factors."""
        b_factors = {}
        for res in residues:
            ca_atom = res.get("CA")
            if ca_atom:
                b_factors[res.id[1]] = ca_atom.get_bfactor()
        return b_factors

    def _get_atom_radius(self, atom: PDB.Atom.Atom) -> float:
        """Get van der Waals radius for atom."""
        # Standard vdW radii in Angstroms
        radii = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "P": 1.8, "H": 1.2, "F": 1.47, "Cl": 1.75, "Br": 1.85, "I": 1.98}
        return radii.get(atom.element, 1.5)  # Default radius
