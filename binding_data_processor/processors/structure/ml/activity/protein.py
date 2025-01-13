"""Enhanced protein structure prediction and analysis with AlphaFold integration."""

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
import json
import requests
import numpy as np
import torch
from Bio import PDB, SeqIO
from Bio.PDB import *
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeHelpers, rdMolTransforms
from scipy.spatial import ConvexHull

from ..base import MLProcessor
from ...base import BaseStructureProcessor
from .alphafold import AlphaFoldPredictor, AlphaFoldConfig, AlphaFoldResult

logger = logging.getLogger(__name__)


@dataclass
class ProteinConfig:
    """Configuration for protein structure prediction and analysis."""

    # Model configuration
    model_path: Optional[str] = None
    device: str = "cuda" if torch.cuda.is_available() else "cpu"
    cache_dir: Optional[str] = None

    # AlphaFold integration
    use_alphafold: bool = True
    alphafold_config: Optional[AlphaFoldConfig] = None
    use_alphafold_api: bool = True
    alphafold_api_base_url: str = "https://alphafold.ebi.ac.uk/api"
    alphafold_api_timeout: int = 300

    # Analysis configuration
    analyze_pockets: bool = True
    analyze_dynamics: bool = True
    analyze_conservation: bool = True
    analyze_interfaces: bool = True
    analyze_quality: bool = True

    # Feature configuration
    use_pssm: bool = True
    use_dssp: bool = True
    use_surface_features: bool = True

    # Performance
    batch_size: int = 1
    num_workers: int = 4

    # Confidence thresholds
    confidence_thresholds: Dict[str, float] = field(default_factory=lambda: {"high": 90.0, "medium": 70.0, "low": 50.0})


class ProteinStructureAnalyzer:
    """Analyzes protein structures and their properties."""

    def __init__(self, config: Optional[Union[Dict, ProteinConfig]] = None):
        """Initialize protein structure analyzer.

        Args:
            config: Configuration for analysis
        """
        if isinstance(config, dict):
            config = ProteinConfig(**config)
        elif config is None:
            config = ProteinConfig()

        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)
        self.pdb_parser = PDB.PDBParser(QUIET=True)

    def get_structure_properties(self, structure: Structure) -> Dict[str, Any]:
        """Get comprehensive properties of protein structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of structure properties
        """
        try:
            properties = {
                "basic": {
                    "num_residues": self._count_residues(structure),
                    "num_atoms": self._count_atoms(structure),
                    "radius_of_gyration": self._calc_radius_of_gyration(structure),
                    "residue_composition": self._get_residue_composition(structure),
                },
                "surface": self._analyze_surface(structure),
                "secondary_structure": self._analyze_secondary_structure(structure),
                "pockets": self._analyze_binding_pockets(structure) if self.config.analyze_pockets else [],
                "interfaces": self._analyze_interfaces(structure) if self.config.analyze_interfaces else [],
                "dynamics": self._analyze_dynamics(structure) if self.config.analyze_dynamics else {},
                "conservation": self._analyze_conservation(structure) if self.config.analyze_conservation else {},
                "quality": self._calculate_quality_metrics(structure) if self.config.analyze_quality else {},
            }
            return properties

        except Exception as e:
            self.logger.error(f"Error getting structure properties: {str(e)}")
            return {}

    def get_molecular_surface(self, structure: Structure) -> Dict[str, Any]:
        """Get molecular surface representation.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary with surface points and properties
        """
        try:
            # Get atomic coordinates
            coords = []
            radii = []
            for atom in structure.get_atoms():
                coords.append(atom.get_coord())
                radii.append(self._get_atom_radius(atom))

            coords = np.array(coords)
            radii = np.array(radii)

            # Calculate surface points
            from Bio.PDB.SASA import ShrakeRupley

            sr = ShrakeRupley()
            surface_points = sr.get_surface_points(coords, radii)

            return {"points": surface_points, "area": sr.compute(coords, radii), "volume": self._calculate_volume(coords)}

        except Exception as e:
            self.logger.error(f"Error getting molecular surface: {str(e)}")
            return {}

    def find_surface_cavities(self, surface: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Find cavities in molecular surface.

        Args:
            surface: Surface dictionary from get_molecular_surface()

        Returns:
            List of cavity properties
        """
        try:
            cavities = []
            points = surface["points"]

            # Use alpha shape to find cavities
            from scipy.spatial import Delaunay

            tri = Delaunay(points)

            # Find cavity points
            for simplex in tri.simplices:
                center = points[simplex].mean(axis=0)
                radius = np.linalg.norm(points[simplex[0]] - center)

                if radius > 2.0:  # Minimum cavity size
                    cavity = {
                        "center": center,
                        "radius": radius,
                        "volume": (4 / 3) * np.pi * radius**3,
                        "surface_points": points[simplex],
                        "depth": self._calculate_cavity_depth(center, points),
                    }
                    cavities.append(cavity)

            return cavities

        except Exception as e:
            self.logger.error(f"Error finding surface cavities: {str(e)}")
            return []

    def _get_atom_radius(self, atom: PDB.Atom.Atom) -> float:
        """Get van der Waals radius for atom."""
        radii = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "P": 1.8, "H": 1.2, "F": 1.47, "Cl": 1.75, "Br": 1.85, "I": 1.98}
        return radii.get(atom.element, 1.5)

    def _calculate_volume(self, coords: np.ndarray) -> float:
        """Calculate volume of point cloud."""
        try:
            hull = ConvexHull(coords)
            return hull.volume
        except Exception as e:
            self.logger.error(f"Error calculating volume: {str(e)}")
            return 0.0

    def _calculate_cavity_depth(self, center: np.ndarray, surface_points: np.ndarray) -> float:
        """Calculate depth of cavity from surface."""
        try:
            distances = np.linalg.norm(surface_points - center, axis=1)
            return float(np.min(distances))
        except Exception as e:
            self.logger.error(f"Error calculating cavity depth: {str(e)}")
            return 0.0

    def _count_residues(self, structure: Structure) -> int:
        """Count residues in structure."""
        try:
            count = 0
            for model in structure:
                for chain in model:
                    count += len(list(chain.get_residues()))
            return count
        except Exception as e:
            self.logger.error(f"Error counting residues: {str(e)}")
            return 0

    def _count_atoms(self, structure: Structure) -> int:
        """Count atoms in structure."""
        try:
            return len(list(structure.get_atoms()))
        except Exception as e:
            self.logger.error(f"Error counting atoms: {str(e)}")
            return 0

    def _calc_radius_of_gyration(self, structure: Structure) -> float:
        """Calculate radius of gyration."""
        try:
            coords = []
            masses = []
            for atom in structure.get_atoms():
                coords.append(atom.get_coord())
                masses.append(atom.mass if hasattr(atom, "mass") else 1.0)

            coords = np.array(coords)
            masses = np.array(masses)
            center = np.average(coords, weights=masses, axis=0)

            r2 = np.sum(masses * np.sum((coords - center) ** 2, axis=1))
            total_mass = np.sum(masses)

            return np.sqrt(r2 / total_mass)

        except Exception as e:
            self.logger.error(f"Error calculating radius of gyration: {str(e)}")
            return 0.0

    def _get_residue_composition(self, structure: Structure) -> Dict[str, int]:
        """Get amino acid composition."""
        try:
            composition = {}
            for residue in structure.get_residues():
                resname = residue.get_resname()
                composition[resname] = composition.get(resname, 0) + 1
            return composition
        except Exception as e:
            self.logger.error(f"Error getting residue composition: {str(e)}")
            return {}

    def _analyze_surface(self, structure: Structure) -> Dict[str, Any]:
        """Analyze surface properties."""
        try:
            surface = self.get_molecular_surface(structure)
            return {"area": surface["area"], "volume": surface["volume"], "properties": self._analyze_surface_properties(structure)}
        except Exception as e:
            self.logger.error(f"Error analyzing surface: {str(e)}")
            return {}

    def _analyze_surface_properties(self, structure: Structure) -> Dict[str, float]:
        """Analyze chemical properties of surface."""
        try:
            properties = {
                "hydrophobicity": self._calculate_surface_hydrophobicity(structure),
                "charge": self._calculate_surface_charge(structure),
                "polarity": self._calculate_surface_polarity(structure),
            }
            return properties
        except Exception as e:
            self.logger.error(f"Error analyzing surface properties: {str(e)}")
            return {}

    def _analyze_secondary_structure(self, structure: Structure) -> Dict[str, float]:
        """Analyze secondary structure composition."""
        try:
            dssp = dssp_dict_from_pdb_file(structure.id)[0]
            ss_counts = {"H": 0, "B": 0, "E": 0, "G": 0, "I": 0, "T": 0, "S": 0}

            for residue in dssp:
                ss = dssp[residue][2]
                if ss in ss_counts:
                    ss_counts[ss] += 1

            total = sum(ss_counts.values())
            if total > 0:
                return {k: v / total for k, v in ss_counts.items()}
            return ss_counts

        except Exception as e:
            self.logger.error(f"Error analyzing secondary structure: {str(e)}")
            return {}

    def _analyze_binding_pockets(self, structure: Structure) -> List[Dict[str, Any]]:
        """Analyze potential binding pockets."""
        try:
            surface = self.get_molecular_surface(structure)
            cavities = self.find_surface_cavities(surface)

            pockets = []
            for cavity in cavities:
                if cavity["volume"] > 100:  # Minimum pocket volume
                    pocket = {
                        "center": cavity["center"],
                        "volume": cavity["volume"],
                        "depth": cavity["depth"],
                        "residues": self._get_pocket_residues(structure, cavity["center"]),
                        "properties": self._analyze_pocket_properties(structure, cavity),
                    }
                    pockets.append(pocket)

            return pockets

        except Exception as e:
            self.logger.error(f"Error analyzing binding pockets: {str(e)}")
            return []

    def _analyze_interfaces(self, structure: Structure) -> List[Dict[str, Any]]:
        """Analyze protein-protein interfaces."""
        try:
            interfaces = []
            for model in structure:
                chains = list(model.get_chains())
                for i in range(len(chains)):
                    for j in range(i + 1, len(chains)):
                        interface = self._analyze_chain_interface(chains[i], chains[j])
                        if interface:
                            interfaces.append(interface)
            return interfaces

        except Exception as e:
            self.logger.error(f"Error analyzing interfaces: {str(e)}")
            return []

    def _analyze_dynamics(self, structure: Structure) -> Dict[str, Any]:
        """Analyze protein dynamics."""
        try:
            dynamics = {"b_factors": self._analyze_b_factors(structure), "flexibility": self._analyze_flexibility(structure), "domains": self._analyze_domains(structure)}
            return dynamics

        except Exception as e:
            self.logger.error(f"Error analyzing dynamics: {str(e)}")
            return {}

    def _analyze_conservation(self, structure: Structure) -> Dict[str, float]:
        """Analyze sequence conservation."""
        try:
            # Placeholder - would need sequence alignment data
            return {}
        except Exception as e:
            self.logger.error(f"Error analyzing conservation: {str(e)}")
            return {}

    def _calculate_quality_metrics(self, structure: Structure) -> Dict[str, float]:
        """Calculate structure quality metrics."""
        try:
            metrics = {
                "clashes": self._count_clashes(structure),
                "rama_outliers": self._count_ramachandran_outliers(structure),
                "rotamer_outliers": self._count_rotamer_outliers(structure),
            }
            return metrics

        except Exception as e:
            self.logger.error(f"Error calculating quality metrics: {str(e)}")
            return {}

    def _get_pocket_residues(self, structure: Structure, center: np.ndarray, cutoff: float = 8.0) -> List[int]:
        """Get residues within cutoff distance of pocket center."""
        try:
            pocket_residues = []
            for residue in structure.get_residues():
                ca = residue["CA"]
                if ca:
                    dist = np.linalg.norm(ca.get_coord() - center)
                    if dist <= cutoff:
                        pocket_residues.append(residue.get_id()[1])
            return pocket_residues

        except Exception as e:
            self.logger.error(f"Error getting pocket residues: {str(e)}")
            return []

    def _analyze_pocket_properties(self, structure: Structure, cavity: Dict[str, Any]) -> Dict[str, Any]:
        """Analyze chemical properties of binding pocket."""
        try:
            residues = self._get_pocket_residues(structure, cavity["center"])
            return {
                "hydrophobicity": self._calculate_pocket_hydrophobicity(structure, residues),
                "charge": self._calculate_pocket_charge(structure, residues),
                "volume": cavity["volume"],
                "depth": cavity["depth"],
            }
        except Exception as e:
            self.logger.error(f"Error analyzing pocket properties: {str(e)}")
            return {}

    def _analyze_chain_interface(self, chain1: Chain, chain2: Chain) -> Optional[Dict[str, Any]]:
        """Analyze interface between two chains."""
        try:
            contacts = []
            for res1 in chain1:
                for res2 in chain2:
                    if self._residues_in_contact(res1, res2):
                        contacts.append((res1.get_id()[1], res2.get_id()[1]))

            if contacts:
                return {"chains": (chain1.get_id(), chain2.get_id()), "contacts": contacts, "area": self._calculate_interface_area(chain1, chain2)}
            return None

        except Exception as e:
            self.logger.error(f"Error analyzing chain interface: {str(e)}")
            return None

    def _residues_in_contact(self, res1: PDB.Residue, res2: PDB.Residue, cutoff: float = 5.0) -> bool:
        """Check if two residues are in contact."""
        try:
            for atom1 in res1:
                for atom2 in res2:
                    dist = np.linalg.norm(atom1.get_coord() - atom2.get_coord())
                    if dist <= cutoff:
                        return True
            return False

        except Exception as e:
            self.logger.error(f"Error checking residue contact: {str(e)}")
            return False

    def _calculate_interface_area(self, chain1: Chain, chain2: Chain) -> float:
        """Calculate interface surface area between chains."""
        try:
            # Placeholder - would need more sophisticated SASA calculation
            return 0.0
        except Exception as e:
            self.logger.error(f"Error calculating interface area: {str(e)}")
            return 0.0

    def _calculate_surface_hydrophobicity(self, structure: Structure) -> float:
        """Calculate surface hydrophobicity."""
        try:
            # Placeholder - would need surface residue identification
            return 0.0
        except Exception as e:
            self.logger.error(f"Error calculating surface hydrophobicity: {str(e)}")
            return 0.0

    def _calculate_surface_charge(self, structure: Structure) -> float:
        """Calculate surface charge."""
        try:
            # Placeholder - would need surface residue identification
            return 0.0
        except Exception as e:
            self.logger.error(f"Error calculating surface charge: {str(e)}")
            return 0.0

    def _calculate_surface_polarity(self, structure: Structure) -> float:
        """Calculate surface polarity."""
        try:
            # Placeholder - would need surface residue identification
            return 0.0
        except Exception as e:
            self.logger.error(f"Error calculating surface polarity: {str(e)}")
            return 0.0

    def _calculate_pocket_hydrophobicity(self, structure: Structure, residues: List[int]) -> float:
        """Calculate pocket hydrophobicity."""
        try:
            # Placeholder - would need hydrophobicity scale
            return 0.0
        except Exception as e:
            self.logger.error(f"Error calculating pocket hydrophobicity: {str(e)}")
            return 0.0

    def _calculate_pocket_charge(self, structure: Structure, residues: List[int]) -> float:
        """Calculate pocket charge."""
        try:
            # Placeholder - would need charge scale
            return 0.0
        except Exception as e:
            self.logger.error(f"Error calculating pocket charge: {str(e)}")
            return 0.0

    def _analyze_b_factors(self, structure: Structure) -> Dict[str, float]:
        """Analyze B-factors."""
        try:
            # Placeholder - would need B-factor analysis
            return {}
        except Exception as e:
            self.logger.error(f"Error analyzing B-factors: {str(e)}")
            return {}

    def _analyze_flexibility(self, structure: Structure) -> Dict[str, float]:
        """Analyze flexibility."""
        try:
            # Placeholder - would need flexibility analysis
            return {}
        except Exception as e:
            self.logger.error(f"Error analyzing flexibility: {str(e)}")
            return {}

    def _analyze_domains(self, structure: Structure) -> Dict[str, Any]:
        """Analyze protein domains."""
        try:
            # Placeholder - would need domain analysis
            return {}
        except Exception as e:
            self.logger.error(f"Error analyzing domains: {str(e)}")
            return {}

    def _count_clashes(self, structure: Structure) -> int:
        """Count atomic clashes."""
        try:
            # Placeholder - would need clash detection
            return 0
        except Exception as e:
            self.logger.error(f"Error counting clashes: {str(e)}")
            return 0

    def _count_ramachandran_outliers(self, structure: Structure) -> int:
        """Count Ramachandran plot outliers."""
        try:
            # Placeholder - would need Ramachandran analysis
            return 0
        except Exception as e:
            self.logger.error(f"Error counting Ramachandran outliers: {str(e)}")
            return 0

    def _count_rotamer_outliers(self, structure: Structure) -> int:
        """Count rotamer outliers."""
        try:
            # Placeholder - would need rotamer analysis
            return 0
        except Exception as e:
            self.logger.error(f"Error counting rotamer outliers: {str(e)}")
            return 0


class ProteinPredictor(MLProcessor, BaseStructureProcessor):
    """Advanced protein structure prediction and analysis with AlphaFold integration."""

    def __init__(
        self,
        config: Optional[Union[Dict, ProteinConfig]] = None,
        **kwargs,
    ):
        """Initialize protein predictor.

        Args:
            config: Predictor configuration
            **kwargs: Additional arguments passed to parent classes
        """
        if isinstance(config, dict):
            config = ProteinConfig(**config)
        elif config is None:
            config = ProteinConfig()

        super().__init__(
            model_path=config.model_path,
            device=config.device,
            **kwargs,
        )
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)

        if self.config.cache_dir:
            Path(self.config.cache_dir).mkdir(parents=True, exist_ok=True)

        if self.config.use_alphafold and not self.config.use_alphafold_api:
            self.alphafold = AlphaFoldPredictor(
                config=self.config.alphafold_config,
                cache_dir=self.config.cache_dir,
            )

        self.pdb_parser = PDB.PDBParser(QUIET=True)
        self.structure_builder = PDB.StructureBuilder()
        self.analyzer = ProteinStructureAnalyzer(config=config)

    async def predict_structure(
        self,
        sequence: str,
        use_cache: bool = True,
        use_alphafold: Optional[bool] = None,
        confidence_threshold: Optional[float] = None,
        **kwargs,
    ) -> Optional[Structure]:
        """Predict protein structure from sequence.

        Args:
            sequence: Amino acid sequence
            use_cache: Whether to use cached predictions
            use_alphafold: Whether to use AlphaFold (overrides config)
            confidence_threshold: Minimum confidence score threshold
            **kwargs: Additional arguments

        Returns:
            Predicted structure or None if prediction fails
        """
        try:
            if use_cache:
                cached = self._load_cached_prediction(sequence)
                if cached:
                    return cached

            if use_alphafold is None:
                use_alphafold = self.config.use_alphafold

            if not confidence_threshold:
                confidence_threshold = self.config.confidence_thresholds["medium"]

            if use_alphafold:
                if self.config.use_alphafold_api:
                    return await self._predict_with_alphafold_api(sequence, confidence_threshold)
                else:
                    result = await self._predict_with_alphafold(sequence)
                    if result:
                        return result
            else:
                return await self._predict_with_local_model(sequence, **kwargs)

            return None

        except Exception as e:
            self.logger.error(f"Error predicting structure: {str(e)}")
            return None

    async def _predict_with_alphafold_api(
        self,
        sequence: str,
        confidence_threshold: float,
    ) -> Optional[Structure]:
        """Predict structure using AlphaFold API."""
        try:
            # Search AlphaFold database
            search_params = {"query": sequence}
            response = requests.get(
                f"{self.config.alphafold_api_base_url}/prediction/search",
                params=search_params,
                timeout=self.config.alphafold_api_timeout,
            )
            response.raise_for_status()
            results = response.json()

            if not results:
                self.logger.warning(f"No AlphaFold structure found for sequence")
                return None

            # Get best model by confidence score
            best_model = max(results, key=lambda x: x["confidence_score"])
            if best_model["confidence_score"] < confidence_threshold:
                self.logger.warning(f"Best model confidence ({best_model['confidence_score']:.1f}) " f"below threshold ({confidence_threshold})")
                return None

            # Download structure
            download_params = {"id": best_model["id"]}
            response = requests.get(
                f"{self.config.alphafold_api_base_url}/prediction/download",
                params=download_params,
                timeout=self.config.alphafold_api_timeout,
            )
            response.raise_for_status()

            # Parse structure
            structure = self.pdb_parser.get_structure("prediction", response.content)

            # Cache if enabled
            if self.config.cache_dir:
                cache_file = self.config.cache_dir / f"{sequence_hash}.pdb"
                with open(cache_file, "wb") as f:
                    f.write(response.content)

            return structure

        except Exception as e:
            self.logger.error(f"Error in AlphaFold API prediction: {str(e)}")
            return None

    async def _predict_with_alphafold(self, sequence: str) -> Optional[Structure]:
        """Predict structure using local AlphaFold."""
        try:
            if not hasattr(self, "alphafold"):
                self.logger.warning("Local AlphaFold not initialized")
                return None

            result = await self.alphafold.predict(sequence)
            if result and isinstance(result, AlphaFoldResult):
                return result.structure
            return None

        except Exception as e:
            self.logger.error(f"Error in local AlphaFold prediction: {str(e)}")
            return None

    async def _predict_with_local_model(self, sequence: str, **kwargs) -> Optional[Structure]:
        """Predict structure using local model."""
        try:
            # TODO: Implement local model prediction
            self.logger.warning("Local model prediction not yet implemented")
            return None

        except Exception as e:
            self.logger.error(f"Error in local model prediction: {str(e)}")
            return None

    def _load_cached_prediction(self, sequence: str) -> Optional[Structure]:
        """Load cached structure prediction."""
        if not self.config.cache_dir:
            return None

        try:
            # TODO: Implement cache loading
            return None
        except Exception as e:
            self.logger.error(f"Error loading cached prediction: {str(e)}")
            return None

    def analyze_structure(self, structure: Structure, **kwargs) -> Dict[str, Any]:
        """Analyze protein structure using the analyzer."""
        return self.analyzer.get_structure_properties(structure)
