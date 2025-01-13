"""Advanced binding site prediction and analysis using ML and structural approaches."""

import logging
from typing import Dict, List, Optional, Tuple, Union, Any
from pathlib import Path
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms, rdShapeHelpers
import torch
import torch.nn.functional as F
from scipy.spatial import ConvexHull, Voronoi
from scipy.spatial.distance import cdist

import requests
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from Bio.PDB.Polypeptide import three_to_one, is_aa
from Bio.PDB.vectors import calc_angle, calc_dihedral
from Bio.PDB.SASA import calculate_sasa
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom


class BindingSitePredictor:
    """Advanced binding site prediction and analysis with ML integration."""

    # Comprehensive residue property scales
    HYDROPHOBICITY = {
        "ALA": 1.8,
        "ARG": -4.5,
        "ASN": -3.5,
        "ASP": -3.5,
        "CYS": 2.5,
        "GLN": -3.5,
        "GLU": -3.5,
        "GLY": -0.4,
        "HIS": -3.2,
        "ILE": 4.5,
        "LEU": 3.8,
        "LYS": -3.9,
        "MET": 1.9,
        "PHE": 2.8,
        "PRO": -1.6,
        "SER": -0.8,
        "THR": -0.7,
        "TRP": -0.9,
        "TYR": -1.3,
        "VAL": 4.2,
    }

    VOLUME = {
        "ALA": 88.6,
        "ARG": 173.4,
        "ASN": 114.1,
        "ASP": 111.1,
        "CYS": 108.5,
        "GLN": 143.8,
        "GLU": 138.4,
        "GLY": 60.1,
        "HIS": 153.2,
        "ILE": 166.7,
        "LEU": 166.7,
        "LYS": 168.6,
        "MET": 162.9,
        "PHE": 189.9,
        "PRO": 112.7,
        "SER": 89.0,
        "THR": 116.1,
        "TRP": 227.8,
        "TYR": 193.6,
        "VAL": 140.0,
    }

    CHARGE = {"ARG": 1, "LYS": 1, "ASP": -1, "GLU": -1, "HIS": 0.1}  # pKa dependent

    def __init__(
        self,
        model_dir: Optional[str] = None,
        use_ml: bool = True,
        device: Optional[str] = None,
    ):
        """Initialize binding site predictor with ML capabilities.

        Args:
            model_dir: Directory containing trained models
            use_ml: Whether to use ML-based prediction
            device: Device to run ML models on
        """
        self.logger = logging.getLogger(__name__)
        self.model_dir = model_dir
        self.use_ml = use_ml
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")

        # Detection parameters
        self.pocket_params = {
            "min_volume": 100.0,  # Minimum pocket volume in Å³
            "max_exposure": 0.5,  # Maximum solvent exposure
            "min_depth": 4.0,  # Minimum pocket depth in Å
            "conservation_cutoff": 0.7,  # Minimum conservation score
            "interaction_distance": 4.5,  # Maximum distance for interactions
            "probe_radius": 1.4,  # Solvent probe radius in Å
            "grid_spacing": 0.5,  # Grid spacing for energy calculations
            "energy_cutoff": -2.0,  # Energy cutoff for favorable regions
        }

        if use_ml and model_dir:
            self._load_models()

    def predict_binding_sites(
        self,
        structure: Structure,
        ligand: Optional[Union[str, Chem.Mol]] = None,
        include_ml_scores: bool = True,
        min_score: float = 0.5,
        max_pockets: int = 5,
    ) -> List[Dict[str, Any]]:
        """Predict potential binding sites using multiple approaches.

        Args:
            structure: BioPython Structure object
            ligand: Optional ligand as SMILES or RDKit mol
            include_ml_scores: Whether to include ML model predictions
            min_score: Minimum confidence score threshold
            max_pockets: Maximum number of pockets to return

        Returns:
            List of dictionaries containing binding site predictions and properties
        """
        try:
            # Get structure properties
            coords, properties = self._get_structure_properties(structure)

            # Geometric pocket detection
            geometric_pockets = self._detect_geometric_pockets(
                coords,
                min_volume=self.pocket_params["min_volume"],
                probe_radius=self.pocket_params["probe_radius"],
            )

            # Energy-based detection
            energy_pockets = self._detect_energy_pockets(
                coords,
                grid_spacing=self.pocket_params["grid_spacing"],
                energy_cutoff=self.pocket_params["energy_cutoff"],
            )

            # Conservation analysis
            conservation = self._analyze_conservation(structure)

            # Combine predictions
            combined_pockets = self._combine_predictions(
                geometric_pockets,
                energy_pockets,
                conservation,
                properties,
            )

            # ML scoring if requested
            if include_ml_scores and self.use_ml and hasattr(self, "ml_model"):
                combined_pockets = self._add_ml_scores(combined_pockets, structure)

            # Filter and rank pockets
            scored_pockets = sorted(
                [p for p in combined_pockets if p["score"] >= min_score],
                key=lambda x: x["score"],
                reverse=True,
            )[:max_pockets]

            # Add detailed analysis
            for pocket in scored_pockets:
                pocket["properties"] = self._analyze_pocket_properties(
                    structure,
                    pocket["residues"],
                    ligand=ligand if isinstance(ligand, Chem.Mol) else None,
                )

            return scored_pockets

        except Exception as e:
            self.logger.error(f"Error predicting binding sites: {str(e)}")
            return []

    def analyze_pocket(
        self,
        structure: Structure,
        pocket_residues: List[int],
        include_conservation: bool = True,
        include_dynamics: bool = False,
    ) -> Dict[str, Any]:
        """Comprehensive analysis of binding pocket properties.

        Args:
            structure: BioPython Structure object
            pocket_residues: List of residue numbers in pocket
            include_conservation: Whether to include conservation analysis
            include_dynamics: Whether to include dynamic properties

        Returns:
            Dictionary containing detailed pocket properties
        """
        try:
            # Get pocket residues
            residues = [structure[0]["A"][res_num] for res_num in pocket_residues]

            # Basic properties
            properties = {
                "volume": self._calculate_volume(residues),
                "surface_area": self._calculate_surface_area(residues),
                "depth": self._calculate_pocket_depth(residues),
                "hydrophobicity": self._calculate_hydrophobicity(residues),
                "charge": self._calculate_charge(residues),
                "residue_composition": self._get_residue_composition(residues),
                "secondary_structure": self._analyze_secondary_structure(residues),
                "shape_descriptors": self._calculate_shape_descriptors(residues),
                "electrostatics": self._calculate_electrostatics(residues),
            }

            # Conservation analysis
            if include_conservation:
                properties["conservation"] = self._calculate_conservation(residues)

            # Dynamic properties
            if include_dynamics:
                properties["dynamics"] = self._analyze_pocket_dynamics(
                    structure,
                    residues,
                )

            return properties

        except Exception as e:
            self.logger.error(f"Error analyzing pocket: {str(e)}")
            return {}

    def _get_structure_properties(
        self,
        structure: Structure,
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        """Extract key properties from protein structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Tuple of (coordinates array, property dictionary)
        """
        coords = []
        properties = {
            "residues": [],
            "hydrophobicity": [],
            "charge": [],
            "conservation": [],
        }

        for residue in structure.get_residues():
            res_name = residue.get_resname()
            properties["residues"].append(residue)
            properties["hydrophobicity"].append(self.HYDROPHOBICITY.get(res_name, 0.0))
            properties["charge"].append(self.CHARGE.get(res_name, 0.0))

            for atom in residue:
                coords.append(atom.get_coord())

        properties["coords"] = np.array(coords)
        properties["hydrophobicity"] = np.array(properties["hydrophobicity"])
        properties["charge"] = np.array(properties["charge"])

        return np.array(coords), properties

    def _load_models(self) -> None:
        """Load ML models for pocket prediction and scoring."""
        try:
            model_path = Path(self.model_dir) / "pocket_scorer.pt"
            if model_path.exists():
                self.ml_model = torch.load(model_path, map_location=self.device)
                self.ml_model.eval()
            else:
                self.logger.warning(f"Model not found at {model_path}")
                self.ml_model = None

        except Exception as e:
            self.logger.warning(f"Failed to load ML models: {e}")
            self.ml_model = None

    def _analyze_pocket_dynamics(
        self,
        structure: Structure,
        residues: List[Residue],
    ) -> Dict[str, Any]:
        """Analyze dynamic properties of binding pocket.

        Args:
            structure: BioPython Structure object
            residues: List of pocket residues

        Returns:
            Dictionary of dynamic properties
        """
        try:
            # Calculate B-factors
            b_factors = []
            for res in residues:
                for atom in res:
                    b_factors.append(atom.get_bfactor())

            # Analyze flexibility
            dynamics = {
                "average_bfactor": float(np.mean(b_factors)),
                "bfactor_std": float(np.std(b_factors)),
                "relative_flexibility": float(np.mean(b_factors) / np.mean([atom.get_bfactor() for atom in structure.get_atoms()])),
            }

            return dynamics

        except Exception as e:
            self.logger.error(f"Error analyzing dynamics: {str(e)}")
            return {}

    def _detect_geometric_pockets(
        self,
        coords: np.ndarray,
        min_volume: float = 100.0,
        probe_radius: float = 1.4,
    ) -> List[Dict[str, Any]]:
        """Detect pockets using geometric criteria and alpha shapes.

        Args:
            coords: Atomic coordinates
            min_volume: Minimum pocket volume in Å³
            probe_radius: Probe radius for surface calculation

        Returns:
            List of detected pockets with properties
        """
        try:
            # Calculate alpha shape
            alpha_shape = self._calculate_alpha_shape(coords, probe_radius)

            # Get Voronoi diagram
            vor = Voronoi(coords)

            # Find cavities using alpha shape and Voronoi vertices
            cavities = []
            for v in vor.vertices:
                if self._is_cavity(v, alpha_shape, coords, probe_radius):
                    cavity = {
                        "center": v,
                        "volume": self._estimate_cavity_volume(v, vor, coords),
                        "residues": self._get_cavity_residues(v, coords),
                        "depth": self._calculate_cavity_depth(v, coords),
                        "exposure": self._calculate_exposure(v, coords),
                    }
                    if cavity["volume"] >= min_volume:
                        cavities.append(cavity)

            return cavities

        except Exception as e:
            self.logger.error(f"Error in geometric detection: {str(e)}")
            return []

    def _detect_energy_pockets(
        self,
        coords: np.ndarray,
        grid_spacing: float = 0.5,
        energy_cutoff: float = -2.0,
    ) -> List[Dict[str, Any]]:
        """Detect pockets using energy-based approaches.

        Args:
            coords: Atomic coordinates
            grid_spacing: Grid spacing for energy calculations
            energy_cutoff: Energy cutoff for favorable regions

        Returns:
            List of detected pockets with properties
        """
        try:
            # Calculate energy grid
            grid_coords, energy_grid = self._calculate_energy_grid(
                coords,
                spacing=grid_spacing,
            )

            # Find favorable regions
            regions = []
            for i, point in enumerate(grid_coords):
                if energy_grid[i] < energy_cutoff:
                    region = {
                        "center": point,
                        "energy": float(energy_grid[i]),
                        "residues": self._get_region_residues(point, coords),
                        "interactions": self._analyze_interactions(point, coords),
                    }
                    regions.append(region)

            # Cluster nearby regions
            clustered = self._cluster_regions(regions)
            return clustered

        except Exception as e:
            self.logger.error(f"Error in energy detection: {str(e)}")
            return []

    def _calculate_energy_grid(
        self,
        coords: np.ndarray,
        spacing: float = 0.5,
        padding: float = 5.0,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate energy grid around protein.

        Args:
            coords: Atomic coordinates
            spacing: Grid spacing
            padding: Padding around protein

        Returns:
            Tuple of (grid coordinates, energy values)
        """
        try:
            # Define grid boundaries
            min_coords = coords.min(axis=0) - padding
            max_coords = coords.max(axis=0) + padding

            # Create grid points
            x = np.arange(min_coords[0], max_coords[0], spacing)
            y = np.arange(min_coords[1], max_coords[1], spacing)
            z = np.arange(min_coords[2], max_coords[2], spacing)

            grid_coords = np.array([(x, y, z) for x in x for y in y for z in z])

            # Calculate energy at each point
            energies = np.zeros(len(grid_coords))
            for i, point in enumerate(grid_coords):
                energies[i] = self._calculate_point_energy(point, coords)

            return grid_coords, energies

        except Exception as e:
            self.logger.error(f"Error calculating energy grid: {str(e)}")
            return np.array([]), np.array([])

    def _calculate_point_energy(
        self,
        point: np.ndarray,
        coords: np.ndarray,
        vdw_cutoff: float = 8.0,
    ) -> float:
        """Calculate energy at a point considering protein atoms.

        Args:
            point: 3D coordinates of point
            coords: Atomic coordinates
            vdw_cutoff: Cutoff for van der Waals interactions

        Returns:
            Energy value at point
        """
        try:
            # Calculate distances to all atoms
            distances = np.linalg.norm(coords - point, axis=1)

            # Van der Waals energy
            vdw_energy = np.sum(4.0 * ((3.5 / distances[distances < vdw_cutoff]) ** 12 - (3.5 / distances[distances < vdw_cutoff]) ** 6))

            # Electrostatic energy could be added here

            return float(vdw_energy)

        except Exception as e:
            self.logger.error(f"Error calculating point energy: {str(e)}")
            return 0.0

    def _cluster_regions(
        self,
        regions: List[Dict[str, Any]],
        distance_cutoff: float = 2.0,
    ) -> List[Dict[str, Any]]:
        """Cluster nearby regions into pockets.

        Args:
            regions: List of favorable regions
            distance_cutoff: Maximum distance for clustering

        Returns:
            List of clustered pockets
        """
        try:
            if not regions:
                return []

            # Get region centers
            centers = np.array([r["center"] for r in regions])

            # Calculate distance matrix
            distances = cdist(centers, centers)

            # Cluster regions
            clusters = []
            used = set()

            for i in range(len(regions)):
                if i in used:
                    continue

                # Find nearby regions
                nearby = np.where(distances[i] < distance_cutoff)[0]
                cluster_regions = [regions[j] for j in nearby]

                # Create cluster
                cluster = {
                    "center": np.mean([r["center"] for r in cluster_regions], axis=0),
                    "energy": np.mean([r["energy"] for r in cluster_regions]),
                    "residues": list(set().union(*[r["residues"] for r in cluster_regions])),
                    "size": len(cluster_regions),
                }

                clusters.append(cluster)
                used.update(nearby)

            return clusters

        except Exception as e:
            self.logger.error(f"Error clustering regions: {str(e)}")
            return []

    def _calculate_alpha_shape(
        self,
        coords: np.ndarray,
        probe_radius: float = 1.4,
    ) -> Any:
        """Calculate alpha shape of protein surface.

        Args:
            coords: Atomic coordinates
            probe_radius: Probe radius for surface calculation

        Returns:
            Alpha shape object
        """
        try:
            from alphashape import alphashape

            return alphashape(coords, alpha=1.0 / probe_radius)
        except Exception as e:
            self.logger.error(f"Error calculating alpha shape: {str(e)}")
            return None

    def _is_cavity(
        self,
        point: np.ndarray,
        alpha_shape: Any,
        coords: np.ndarray,
        probe_radius: float,
    ) -> bool:
        """Check if a point represents a cavity.

        Args:
            point: Point to check
            alpha_shape: Alpha shape of protein
            coords: Atomic coordinates
            probe_radius: Probe radius

        Returns:
            True if point is in a cavity
        """
        try:
            # Check if point is inside alpha shape
            if not alpha_shape.contains(point):
                return False

            # Check distances to atoms
            distances = np.linalg.norm(coords - point, axis=1)
            return np.all(distances > probe_radius)

        except Exception as e:
            self.logger.error(f"Error checking cavity: {str(e)}")
            return False

    def _estimate_cavity_volume(
        self,
        center: np.ndarray,
        vor: Voronoi,
        coords: np.ndarray,
        max_radius: float = 10.0,
    ) -> float:
        """Estimate cavity volume using Voronoi cells.

        Args:
            center: Cavity center point
            vor: Voronoi diagram
            coords: Atomic coordinates
            max_radius: Maximum radius to consider

        Returns:
            Estimated cavity volume in Å³
        """
        try:
            # Find nearby Voronoi vertices
            distances = np.linalg.norm(vor.vertices - center, axis=1)
            nearby = vor.vertices[distances < max_radius]

            if len(nearby) < 4:
                return 0.0

            # Calculate convex hull volume
            hull = ConvexHull(nearby)
            return hull.volume

        except Exception as e:
            self.logger.error(f"Error estimating cavity volume: {str(e)}")
            return 0.0

    def _get_cavity_residues(
        self,
        center: np.ndarray,
        coords: np.ndarray,
        cutoff: float = 8.0,
    ) -> List[int]:
        """Get residues forming a cavity.

        Args:
            center: Cavity center point
            coords: Atomic coordinates
            cutoff: Distance cutoff for residue inclusion

        Returns:
            List of residue indices
        """
        try:
            # Calculate distances to cavity center
            distances = np.linalg.norm(coords - center, axis=1)

            # Get indices of nearby atoms
            nearby = np.where(distances < cutoff)[0]

            # Convert to residue indices
            residue_indices = list(set([i // 10 for i in nearby]))  # Approximate residue grouping
            return residue_indices

        except Exception as e:
            self.logger.error(f"Error getting cavity residues: {str(e)}")
            return []

    def _calculate_cavity_depth(
        self,
        center: np.ndarray,
        coords: np.ndarray,
    ) -> float:
        """Calculate depth of cavity from protein surface.

        Args:
            center: Cavity center point
            coords: Atomic coordinates

        Returns:
            Cavity depth in Å
        """
        try:
            # Get surface atoms
            surface_atoms = self._get_surface_atoms(coords)

            if len(surface_atoms) == 0:
                return 0.0

            # Calculate minimum distance to surface
            distances = np.linalg.norm(surface_atoms - center, axis=1)
            return float(np.min(distances))

        except Exception as e:
            self.logger.error(f"Error calculating cavity depth: {str(e)}")
            return 0.0

    def _calculate_exposure(
        self,
        center: np.ndarray,
        coords: np.ndarray,
        num_directions: int = 100,
    ) -> float:
        """Calculate solvent exposure of cavity.

        Args:
            center: Cavity center point
            coords: Atomic coordinates
            num_directions: Number of directions to check

        Returns:
            Exposure score (0-1)
        """
        try:
            # Generate random directions
            directions = np.random.randn(num_directions, 3)
            directions /= np.linalg.norm(directions, axis=1)[:, np.newaxis]

            # Ray casting
            exposed = 0
            for direction in directions:
                ray = center + direction * np.arange(0, 20, 0.5)[:, np.newaxis]
                distances = cdist(ray, coords)
                if np.all(distances.min(axis=1) > 2.0):  # No collisions
                    exposed += 1

            return exposed / num_directions

        except Exception as e:
            self.logger.error(f"Error calculating exposure: {str(e)}")
            return 1.0

    def _get_surface_atoms(
        self,
        coords: np.ndarray,
        probe_radius: float = 1.4,
    ) -> np.ndarray:
        """Identify surface atoms using rolling ball algorithm.

        Args:
            coords: Atomic coordinates
            probe_radius: Probe radius for surface calculation

        Returns:
            Array of surface atom coordinates
        """
        try:
            # Calculate distances between all atoms
            distances = cdist(coords, coords)

            # Find atoms with exposed surface area
            surface_atoms = []
            for i, atom_coords in enumerate(coords):
                # Get neighboring atoms
                neighbors = coords[distances[i] < (probe_radius * 3)]

                if len(neighbors) < 12:  # Exposed atom
                    surface_atoms.append(atom_coords)

            return np.array(surface_atoms)

        except Exception as e:
            self.logger.error(f"Error identifying surface atoms: {str(e)}")
            return np.array([])

    def _get_region_residues(
        self,
        point: np.ndarray,
        coords: np.ndarray,
        cutoff: float = 8.0,
    ) -> List[int]:
        """Get residues around a point.

        Args:
            point: Point to analyze
            coords: Atomic coordinates
            cutoff: Distance cutoff

        Returns:
            List of residue indices
        """
        try:
            # Calculate distances to point
            distances = np.linalg.norm(coords - point, axis=1)

            # Get nearby atoms
            nearby = np.where(distances < cutoff)[0]

            # Convert to residue indices (approximate)
            residue_indices = list(set([i // 10 for i in nearby]))
            return residue_indices

        except Exception as e:
            self.logger.error(f"Error getting region residues: {str(e)}")
            return []

    def _analyze_conservation(
        self,
        structure: Structure,
        window_size: int = 5,
    ) -> Dict[str, float]:
        """Analyze sequence conservation using sliding window.

        Args:
            structure: Protein structure
            window_size: Window size for conservation calculation

        Returns:
            Dictionary mapping residue IDs to conservation scores
        """
        try:
            conservation = {}

            # Get sequence
            sequence = ""
            residue_ids = []
            for residue in structure.get_residues():
                sequence += three_to_one(residue.get_resname())
                residue_ids.append(residue.get_id()[1])

            # Calculate conservation scores
            for i in range(len(sequence)):
                # Get window
                start = max(0, i - window_size // 2)
                end = min(len(sequence), i + window_size // 2 + 1)
                window = sequence[start:end]

                # Calculate conservation score
                score = self._calculate_window_conservation(window)
                conservation[residue_ids[i]] = score

            return conservation

        except Exception as e:
            self.logger.error(f"Error analyzing conservation: {str(e)}")
            return {}

    def _calculate_window_conservation(
        self,
        window: str,
        blosum_cutoff: int = 62,
    ) -> float:
        """Calculate conservation score for sequence window.

        Args:
            window: Sequence window
            blosum_cutoff: BLOSUM matrix cutoff

        Returns:
            Conservation score (0-1)
        """
        try:
            # Use BLOSUM matrix for scoring
            from Bio.Align import substitution_matrices

            matrix = substitution_matrices.load("BLOSUM" + str(blosum_cutoff))

            # Calculate average substitution score
            scores = []
            for i in range(len(window)):
                for j in range(i + 1, len(window)):
                    score = matrix[window[i]][window[j]]
                    scores.append(score)

            if not scores:
                return 0.0

            # Normalize score
            avg_score = np.mean(scores)
            max_score = max(matrix[aa][aa] for aa in matrix)
            return float(avg_score / max_score)

        except Exception as e:
            self.logger.error(f"Error calculating conservation: {str(e)}")
            return 0.0

    def _combine_predictions(
        self,
        geometric_pockets: List[Dict[str, Any]],
        energy_pockets: List[Dict[str, Any]],
        conservation: Dict[str, float],
        properties: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Combine predictions from different methods.

        Args:
            geometric_pockets: Geometric pocket predictions
            energy_pockets: Energy-based pocket predictions
            conservation: Conservation scores
            properties: Structure properties

        Returns:
            List of combined pocket predictions
        """
        try:
            combined = []

            # Match overlapping pockets
            for g_pocket in geometric_pockets:
                g_center = g_pocket["center"]

                # Find overlapping energy pockets
                for e_pocket in energy_pockets:
                    e_center = e_pocket["center"]

                    # Check overlap
                    distance = np.linalg.norm(g_center - e_center)
                    if distance < 5.0:  # Overlapping pockets
                        # Combine scores
                        geometric_score = self._score_geometric_pocket(g_pocket)
                        energy_score = self._score_energy_pocket(e_pocket)
                        conservation_score = self._score_conservation(g_pocket["residues"], conservation)

                        # Calculate combined score
                        score = 0.4 * geometric_score + 0.4 * energy_score + 0.2 * conservation_score

                        # Create combined pocket
                        pocket = {
                            "score": score,
                            "center": (g_center + e_center) / 2,
                            "residues": list(set(g_pocket["residues"]) | set(e_pocket["residues"])),
                            "volume": g_pocket["volume"],
                            "energy": e_pocket["energy"],
                            "geometric_features": g_pocket,
                            "energy_features": e_pocket,
                        }
                        combined.append(pocket)

            return combined

        except Exception as e:
            self.logger.error(f"Error combining predictions: {str(e)}")
            return []

    def _score_geometric_pocket(self, pocket: Dict[str, Any]) -> float:
        """Score geometric pocket properties.

        Args:
            pocket: Geometric pocket data

        Returns:
            Score between 0 and 1
        """
        try:
            # Score different properties
            volume_score = min(1.0, pocket["volume"] / 1000)  # Volume up to 1000 Å³
            depth_score = min(1.0, pocket["depth"] / 10)  # Depth up to 10 Å
            exposure_score = 1.0 - pocket["exposure"]  # Lower exposure is better

            # Combine scores
            return (volume_score + depth_score + exposure_score) / 3

        except Exception as e:
            self.logger.error(f"Error scoring geometric pocket: {str(e)}")
            return 0.0

    def _score_energy_pocket(self, pocket: Dict[str, Any]) -> float:
        """Score energy-based pocket properties.

        Args:
            pocket: Energy pocket data

        Returns:
            Score between 0 and 1
        """
        try:
            # Score energy
            energy_score = 1.0 - min(1.0, abs(pocket["energy"]) / 10)

            # Score interactions
            interactions = pocket["interactions"]
            interaction_score = min(1.0, (interactions["hydrophobic"] + interactions["hbond"] * 2 + interactions["ionic"] * 2) / 20)

            # Combine scores
            return (energy_score + interaction_score) / 2

        except Exception as e:
            self.logger.error(f"Error scoring energy pocket: {str(e)}")
            return 0.0

    def _score_conservation(
        self,
        residues: List[int],
        conservation: Dict[str, float],
    ) -> float:
        """Score pocket conservation.

        Args:
            residues: List of residue indices
            conservation: Conservation scores

        Returns:
            Score between 0 and 1
        """
        try:
            if not residues or not conservation:
                return 0.0

            # Get conservation scores for pocket residues
            scores = [conservation.get(res, 0.0) for res in residues]

            # Calculate average score
            return float(np.mean(scores))

        except Exception as e:
            self.logger.error(f"Error scoring conservation: {str(e)}")
            return 0.0

    def _get_residue_center(self, residue: Residue) -> np.ndarray:
        """Calculate geometric center of residue.

        Args:
            residue: BioPython residue object

        Returns:
            Array of center coordinates
        """
        try:
            coords = []
            for atom in residue:
                coords.append(atom.get_coord())
            return np.mean(coords, axis=0)
        except Exception as e:
            self.logger.error(f"Error calculating residue center: {str(e)}")
            return np.zeros(3)

    def _get_ring_parameters(self, residue: Residue) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """Get aromatic ring center and normal vector.

        Args:
            residue: BioPython residue object

        Returns:
            Tuple of (ring center, normal vector) or (None, None)
        """
        try:
            if residue.get_resname() not in ["PHE", "TYR", "TRP", "HIS"]:
                return None, None

            # Get ring atoms
            ring_atoms = []
            if residue.get_resname() in ["PHE", "TYR"]:
                ring_atoms = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
            elif residue.get_resname() == "TRP":
                ring_atoms = ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"]
            elif residue.get_resname() == "HIS":
                ring_atoms = ["CG", "ND1", "CD2", "CE1", "NE2"]

            # Get coordinates
            coords = []
            for atom_name in ring_atoms:
                if atom_name in residue:
                    coords.append(residue[atom_name].get_coord())

            if len(coords) < 3:
                return None, None

            coords = np.array(coords)
            center = np.mean(coords, axis=0)

            # Calculate normal vector using cross products
            v1 = coords[1] - coords[0]
            v2 = coords[2] - coords[0]
            normal = np.cross(v1, v2)
            normal = normal / np.linalg.norm(normal)

            return center, normal

        except Exception as e:
            self.logger.error(f"Error getting ring parameters: {str(e)}")
            return None, None

    def _get_hbond_features(self, residue: Residue) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """Get H-bond donor and acceptor features.

        Args:
            residue: BioPython residue object

        Returns:
            Tuple of (donor features, acceptor features)
        """
        try:
            donors = []
            acceptors = []
            center = self._get_residue_center(residue)

            # Add backbone features
            if "N" in residue and "H" in residue:
                donors.append(
                    {
                        "center": residue["N"].get_coord(),
                        "vector": residue["H"].get_coord() - residue["N"].get_coord(),
                        "radius": 1.0,
                        "residue": residue.get_id()[1],
                    }
                )

            if "O" in residue:
                acceptors.append(
                    {
                        "center": residue["O"].get_coord(),
                        "vector": center - residue["O"].get_coord(),
                        "radius": 1.0,
                        "residue": residue.get_id()[1],
                    }
                )

            # Add sidechain features
            res_name = residue.get_resname()
            if res_name in ["SER", "THR", "TYR"]:
                if "OH" in residue:
                    donors.append(
                        {
                            "center": residue["OH"].get_coord(),
                            "vector": center - residue["OH"].get_coord(),
                            "radius": 1.0,
                            "residue": residue.get_id()[1],
                        }
                    )
                    acceptors.append(
                        {
                            "center": residue["OH"].get_coord(),
                            "vector": center - residue["OH"].get_coord(),
                            "radius": 1.0,
                            "residue": residue.get_id()[1],
                        }
                    )

            elif res_name in ["ASN", "GLN"]:
                if "ND2" in residue or "NE2" in residue:
                    atom_name = "ND2" if "ND2" in residue else "NE2"
                    donors.append(
                        {
                            "center": residue[atom_name].get_coord(),
                            "vector": center - residue[atom_name].get_coord(),
                            "radius": 1.0,
                            "residue": residue.get_id()[1],
                        }
                    )

            elif res_name in ["LYS", "ARG"]:
                if "NZ" in residue:
                    donors.append(
                        {
                            "center": residue["NZ"].get_coord(),
                            "vector": center - residue["NZ"].get_coord(),
                            "radius": 1.0,
                            "residue": residue.get_id()[1],
                        }
                    )

            return donors, acceptors

        except Exception as e:
            self.logger.error(f"Error getting H-bond features: {str(e)}")
            return [], []

    def _get_secondary_structure(self, residue: Residue) -> str:
        """Get secondary structure assignment for residue.

        Args:
            residue: BioPython residue object

        Returns:
            Secondary structure type (H, E, C, etc.)
        """
        try:
            # Calculate phi/psi angles
            phi, psi = self._calculate_backbone_angles(residue)

            if phi is None or psi is None:
                return "C"  # Coil

            # Ramachandran plot regions
            if -140 < phi < -60 and -70 < psi < -15:
                return "H"  # Alpha helix
            elif -150 < phi < -50 and 100 < psi < 180:
                return "E"  # Beta sheet
            else:
                return "C"  # Coil

        except Exception as e:
            self.logger.error(f"Error getting secondary structure: {str(e)}")
            return "C"

    def _calculate_backbone_angles(self, residue: Residue) -> Tuple[Optional[float], Optional[float]]:
        """Calculate backbone phi/psi angles.

        Args:
            residue: BioPython residue object

        Returns:
            Tuple of (phi angle, psi angle) in degrees
        """
        try:
            # Get required atoms
            if not all(atom in residue for atom in ["N", "CA", "C"]):
                return None, None

            # Get previous and next residues
            prev_res = residue.get_previous_residue()
            next_res = residue.get_next_residue()

            # Calculate phi angle (requires previous residue)
            phi = None
            if prev_res and "C" in prev_res:
                phi = calc_dihedral(
                    prev_res["C"].get_vector(),
                    residue["N"].get_vector(),
                    residue["CA"].get_vector(),
                    residue["C"].get_vector(),
                )

            # Calculate psi angle (requires next residue)
            psi = None
            if next_res and "N" in next_res:
                psi = calc_dihedral(
                    residue["N"].get_vector(),
                    residue["CA"].get_vector(),
                    residue["C"].get_vector(),
                    next_res["N"].get_vector(),
                )

            return phi, psi

        except Exception as e:
            self.logger.error(f"Error calculating backbone angles: {str(e)}")
            return None, None

    def _analyze_pocket_properties(
        self,
        structure: Structure,
        residues: List[int],
        ligand: Optional[Chem.Mol] = None,
    ) -> Dict[str, Any]:
        """Comprehensive analysis of binding pocket properties.

        Args:
            structure: BioPython Structure object
            residues: List of residue numbers
            ligand: Optional ligand molecule

        Returns:
            Dictionary of pocket properties
        """
        try:
            # Get pocket residues
            pocket_residues = [structure[0]["A"][res_num] for res_num in residues]

            # Basic properties
            properties = {
                "volume": self._calculate_volume(pocket_residues),
                "surface_area": self._calculate_surface_area(pocket_residues),
                "depth": self._calculate_pocket_depth(pocket_residues),
                "hydrophobicity": self._calculate_hydrophobicity(pocket_residues),
                "charge": self._calculate_charge(pocket_residues),
                "residue_composition": self._get_residue_composition(pocket_residues),
                "secondary_structure": self._analyze_secondary_structure(pocket_residues),
                "shape_descriptors": self._calculate_shape_descriptors(pocket_residues),
                "electrostatics": self._calculate_electrostatics(pocket_residues),
                "conservation": self._calculate_conservation(pocket_residues),
                "dynamics": self._analyze_pocket_dynamics(structure, pocket_residues),
                "alphafold_confidence": self._get_alphafold_confidence(structure, residues),
                "structural_features": self._analyze_structural_features(pocket_residues),
                "pharmacophore": self._generate_pharmacophore(pocket_residues),
                "druggability": self._assess_druggability(pocket_residues),
            }

            # Add ligand-based analysis if ligand provided
            if ligand is not None:
                properties.update(
                    {
                        "shape_complementarity": self._calculate_shape_complementarity(structure, residues, ligand),
                        "interaction_energy": self._calculate_interaction_energy(pocket_residues, ligand),
                        "binding_mode": self._predict_binding_mode(pocket_residues, ligand),
                    }
                )

            return properties

        except Exception as e:
            self.logger.error(f"Error analyzing pocket properties: {str(e)}")
            return {}

    def _get_alphafold_confidence(
        self,
        structure: Structure,
        residues: List[int],
        confidence_cutoff: float = 70.0,
    ) -> Dict[str, Any]:
        """Get AlphaFold confidence scores for pocket residues.

        Args:
            structure: BioPython Structure object
            residues: List of residue numbers
            confidence_cutoff: Minimum confidence score

        Returns:
            Dictionary of confidence metrics
        """
        try:
            # Get UniProt ID from structure
            uniprot_id = structure.header.get("idcode", "").split("_")[0]

            if not uniprot_id:
                return {}

            # Query AlphaFold DB
            url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
            response = requests.get(url)

            if not response.ok:
                return {}

            data = response.json()

            # Extract confidence scores
            confidence = {
                "plddt": [],  # Per-residue confidence
                "pae": [],  # Predicted aligned error
                "ptm": None,  # Predicted TM-score
            }

            for res_num in residues:
                if str(res_num) in data["plddt"]:
                    confidence["plddt"].append(data["plddt"][str(res_num)])

            if confidence["plddt"]:
                confidence["average_plddt"] = float(np.mean(confidence["plddt"]))
                confidence["min_plddt"] = float(np.min(confidence["plddt"]))
                confidence["reliable"] = confidence["min_plddt"] >= confidence_cutoff

            if "predicted_tm_score" in data:
                confidence["ptm"] = float(data["predicted_tm_score"])

            return confidence

        except Exception as e:
            self.logger.error(f"Error getting AlphaFold confidence: {str(e)}")
            return {}

    def _analyze_structural_features(
        self,
        residues: List[Residue],
    ) -> Dict[str, Any]:
        """Analyze detailed structural features of pocket residues.

        Args:
            residues: List of pocket residues

        Returns:
            Dictionary of structural features
        """
        try:
            features = {
                "secondary_structure": {},
                "accessibility": [],
                "phi_psi": [],
                "hbonds": [],
                "salt_bridges": [],
                "disulfides": [],
            }

            # Calculate secondary structure
            for res in residues:
                ss = self._get_secondary_structure(res)
                features["secondary_structure"][res.get_id()[1]] = ss

                # Calculate accessibility
                sasa = calculate_sasa(res)
                features["accessibility"].append(sasa)

                # Get backbone angles
                phi, psi = self._calculate_backbone_angles(res)
                if phi is not None and psi is not None:
                    features["phi_psi"].append((phi, psi))

            # Find interactions
            features["hbonds"] = self._find_hbonds(residues)
            features["salt_bridges"] = self._find_salt_bridges(residues)
            features["disulfides"] = self._find_disulfides(residues)

            # Calculate statistics
            features["average_sasa"] = float(np.mean(features["accessibility"]))
            features["buried_residues"] = sum(1 for sasa in features["accessibility"] if sasa < 20)
            features["exposed_residues"] = sum(1 for sasa in features["accessibility"] if sasa >= 20)

            return features

        except Exception as e:
            self.logger.error(f"Error analyzing structural features: {str(e)}")
            return {}

    def _generate_pharmacophore(
        self,
        residues: List[Residue],
        include_vectors: bool = True,
    ) -> Dict[str, Any]:
        """Generate pharmacophore model from pocket residues.

        Args:
            residues: List of pocket residues
            include_vectors: Whether to include directional features

        Returns:
            Dictionary containing pharmacophore features
        """
        try:
            features = {
                "donors": [],
                "acceptors": [],
                "aromatic": [],
                "hydrophobic": [],
                "positive": [],
                "negative": [],
                "excluded_volumes": [],
            }

            for res in residues:
                # Get residue center
                center = self._get_residue_center(res)

                # Analyze residue type
                if res.get_resname() in ["ARG", "LYS", "HIS"]:
                    features["positive"].append(
                        {
                            "center": center,
                            "radius": 2.0,
                            "residue": res.get_id()[1],
                        }
                    )

                elif res.get_resname() in ["ASP", "GLU"]:
                    features["negative"].append(
                        {
                            "center": center,
                            "radius": 2.0,
                            "residue": res.get_id()[1],
                        }
                    )

                elif res.get_resname() in ["PHE", "TYR", "TRP", "HIS"]:
                    ring_center, normal = self._get_ring_parameters(res)
                    if ring_center is not None:
                        features["aromatic"].append(
                            {
                                "center": ring_center,
                                "normal": normal,
                                "radius": 2.5,
                                "residue": res.get_id()[1],
                            }
                        )

                # Add H-bond donors/acceptors
                donors, acceptors = self._get_hbond_features(res)
                features["donors"].extend(donors)
                features["acceptors"].extend(acceptors)

                # Add hydrophobic features
                if res.get_resname() in ["ALA", "VAL", "LEU", "ILE", "MET", "PRO"]:
                    features["hydrophobic"].append(
                        {
                            "center": center,
                            "radius": 2.0,
                            "residue": res.get_id()[1],
                        }
                    )

            return features

        except Exception as e:
            self.logger.error(f"Error generating pharmacophore: {str(e)}")
            return {}

    def _assess_druggability(
        self,
        residues: List[Residue],
    ) -> Dict[str, float]:
        """Assess pocket druggability using multiple metrics.

        Args:
            residues: List of pocket residues

        Returns:
            Dictionary of druggability scores
        """
        try:
            scores = {}

            # Volume-based score
            volume = self._calculate_volume(residues)
            scores["volume_score"] = min(1.0, volume / 1000)  # Normalize to 0-1

            # Hydrophobicity score
            hydrophobicity = self._calculate_hydrophobicity(residues)
            scores["hydrophobicity_score"] = (hydrophobicity + 4.5) / 9.0  # Normalize

            # Composition score
            composition = self._get_residue_composition(residues)
            aromatic = sum(composition.get(aa, 0) for aa in ["PHE", "TYR", "TRP"])
            charged = sum(composition.get(aa, 0) for aa in ["ARG", "LYS", "ASP", "GLU"])
            scores["composition_score"] = min(1.0, (aromatic * 0.4 + charged * 0.3) / len(residues))

            # Accessibility score
            accessibility = [calculate_sasa(res) for res in residues]
            scores["accessibility_score"] = 1.0 - min(1.0, np.mean(accessibility) / 100)

            # Shape score
            shape = self._calculate_shape_descriptors(residues)
            scores["shape_score"] = shape.get("sphericity", 0.0)

            # Combined score (weighted average)
            weights = {
                "volume_score": 0.25,
                "hydrophobicity_score": 0.2,
                "composition_score": 0.2,
                "accessibility_score": 0.15,
                "shape_score": 0.2,
            }

            scores["total_score"] = sum(score * weights[name] for name, score in scores.items() if name in weights)

            return scores

        except Exception as e:
            self.logger.error(f"Error assessing druggability: {str(e)}")
            return {}

    def _calculate_volume(self, residues: List[Residue]) -> float:
        """Calculate volume of residue selection.

        Args:
            residues: List of residues

        Returns:
            Volume in Å³
        """
        try:
            # Get all atom coordinates
            coords = []
            for res in residues:
                for atom in res:
                    coords.append(atom.get_coord())
            coords = np.array(coords)

            if len(coords) < 4:
                return 0.0

            # Calculate convex hull volume
            hull = ConvexHull(coords)
            return hull.volume

        except Exception as e:
            self.logger.error(f"Error calculating volume: {str(e)}")
            return 0.0

    def _calculate_surface_area(self, residues: List[Residue]) -> float:
        """Calculate solvent accessible surface area.

        Args:
            residues: List of residues

        Returns:
            Surface area in Å²
        """
        try:
            # Get coordinates and radii for all atoms
            coords = []
            radii = []
            for res in residues:
                for atom in res:
                    coords.append(atom.get_coord())
                    radii.append(get_atom_radius(atom))

            if not coords:
                return 0.0

            # Calculate SASA using utils.calculate_sasa
            from ..protein.utils import calculate_sasa

            return calculate_sasa(np.array(coords), np.array(radii))

        except Exception as e:
            self.logger.error(f"Error calculating surface area: {str(e)}")
            return 0.0

    def _calculate_hydrophobicity(self, residues: List[Residue]) -> float:
        """Calculate average hydrophobicity.

        Args:
            residues: List of residues

        Returns:
            Average hydrophobicity score
        """
        try:
            scores = []
            for res in residues:
                score = self.HYDROPHOBICITY.get(res.get_resname(), 0.0)
                scores.append(score)
            return float(np.mean(scores)) if scores else 0.0

        except Exception as e:
            self.logger.error(f"Error calculating hydrophobicity: {str(e)}")
            return 0.0

    def _calculate_charge(self, residues: List[Residue]) -> float:
        """Calculate net charge.

        Args:
            residues: List of residues

        Returns:
            Net charge
        """
        try:
            total_charge = 0.0
            for res in residues:
                charge = self.CHARGE.get(res.get_resname(), 0.0)
                total_charge += charge
            return total_charge

        except Exception as e:
            self.logger.error(f"Error calculating charge: {str(e)}")
            return 0.0

    def _get_residue_composition(self, residues: List[Residue]) -> Dict[str, int]:
        """Get amino acid composition.

        Args:
            residues: List of residues

        Returns:
            Dictionary mapping residue types to counts
        """
        try:
            composition = {}
            for res in residues:
                res_name = res.get_resname()
                composition[res_name] = composition.get(res_name, 0) + 1
            return composition

        except Exception as e:
            self.logger.error(f"Error getting residue composition: {str(e)}")
            return {}

    def _analyze_secondary_structure(self, residues: List[Residue]) -> Dict[str, float]:
        """Analyze secondary structure composition.

        Args:
            residues: List of residues

        Returns:
            Dictionary of secondary structure percentages
        """
        try:
            ss_counts = {"H": 0, "E": 0, "C": 0}  # Helix, Sheet, Coil
            total = 0

            for res in residues:
                ss = self._get_secondary_structure(res)
                ss_counts[ss] += 1
                total += 1

            if total == 0:
                return {"H": 0.0, "E": 0.0, "C": 0.0}

            return {k: v / total for k, v in ss_counts.items()}

        except Exception as e:
            self.logger.error(f"Error analyzing secondary structure: {str(e)}")
            return {"H": 0.0, "E": 0.0, "C": 0.0}

    def _calculate_shape_descriptors(self, residues: List[Residue]) -> Dict[str, float]:
        """Calculate shape descriptors for pocket.

        Args:
            residues: List of residues

        Returns:
            Dictionary of shape descriptors
        """
        try:
            # Get coordinates
            coords = []
            for res in residues:
                for atom in res:
                    coords.append(atom.get_coord())
            coords = np.array(coords)

            if len(coords) < 4:
                return {
                    "sphericity": 0.0,
                    "asphericity": 0.0,
                    "eccentricity": 0.0,
                }

            # Calculate principal components
            centered = coords - np.mean(coords, axis=0)
            cov = np.cov(centered.T)
            eigenvals = np.linalg.eigvals(cov)
            eigenvals.sort()

            # Shape descriptors
            a, b, c = eigenvals
            descriptors = {
                "sphericity": c / a,  # 1 for perfect sphere
                "asphericity": (a - (b + c) / 2) / a,  # 0 for perfect sphere
                "eccentricity": np.sqrt(1 - c / a),  # 0 for perfect sphere
            }

            return descriptors

        except Exception as e:
            self.logger.error(f"Error calculating shape descriptors: {str(e)}")
            return {"sphericity": 0.0, "asphericity": 0.0, "eccentricity": 0.0}

    def _calculate_electrostatics(self, residues: List[Residue]) -> Dict[str, Any]:
        """Calculate electrostatic properties.

        Args:
            residues: List of residues

        Returns:
            Dictionary of electrostatic properties
        """
        try:
            # Count charged residues
            pos_charged = sum(1 for res in residues if res.get_resname() in ["ARG", "LYS", "HIS"])
            neg_charged = sum(1 for res in residues if res.get_resname() in ["ASP", "GLU"])

            # Calculate charge density
            volume = self._calculate_volume(residues)
            charge_density = (pos_charged - neg_charged) / volume if volume > 0 else 0.0

            return {
                "positive_charges": pos_charged,
                "negative_charges": neg_charged,
                "net_charge": pos_charged - neg_charged,
                "charge_density": charge_density,
            }

        except Exception as e:
            self.logger.error(f"Error calculating electrostatics: {str(e)}")
            return {
                "positive_charges": 0,
                "negative_charges": 0,
                "net_charge": 0,
                "charge_density": 0.0,
            }

    def _find_hbonds(self, residues: List[Residue]) -> List[Dict[str, Any]]:
        """Find hydrogen bonds between residues.

        Args:
            residues: List of residues

        Returns:
            List of hydrogen bond descriptions
        """
        try:
            hbonds = []
            for i, res1 in enumerate(residues):
                for res2 in residues[i + 1 :]:
                    # Get donors and acceptors
                    donors1, acceptors1 = self._get_hbond_features(res1)
                    donors2, acceptors2 = self._get_hbond_features(res2)

                    # Check all donor-acceptor pairs
                    for donor in donors1:
                        for acceptor in acceptors2:
                            if self._check_hbond(donor, acceptor):
                                hbonds.append(
                                    {
                                        "donor_res": res1.get_id()[1],
                                        "acceptor_res": res2.get_id()[1],
                                        "distance": np.linalg.norm(donor["center"] - acceptor["center"]),
                                    }
                                )

                    for donor in donors2:
                        for acceptor in acceptors1:
                            if self._check_hbond(donor, acceptor):
                                hbonds.append(
                                    {
                                        "donor_res": res2.get_id()[1],
                                        "acceptor_res": res1.get_id()[1],
                                        "distance": np.linalg.norm(donor["center"] - acceptor["center"]),
                                    }
                                )

            return hbonds

        except Exception as e:
            self.logger.error(f"Error finding hydrogen bonds: {str(e)}")
            return []

    def _check_hbond(
        self,
        donor: Dict[str, Any],
        acceptor: Dict[str, Any],
        distance_cutoff: float = 3.5,
        angle_cutoff: float = 30.0,
    ) -> bool:
        """Check if donor-acceptor pair forms hydrogen bond.

        Args:
            donor: Donor feature
            acceptor: Acceptor feature
            distance_cutoff: Maximum H-bond distance in Å
            angle_cutoff: Maximum H-bond angle deviation in degrees

        Returns:
            True if H-bond criteria are met
        """
        try:
            # Check distance
            distance = np.linalg.norm(donor["center"] - acceptor["center"])
            if distance > distance_cutoff:
                return False

            # Check angle if vectors available
            if donor.get("vector") is not None and acceptor.get("vector") is not None:
                donor_vec = donor["vector"]
                acceptor_vec = acceptor["vector"]
                angle = np.degrees(np.arccos(np.dot(donor_vec, acceptor_vec)))
                if abs(180 - angle) > angle_cutoff:
                    return False

            return True

        except Exception as e:
            self.logger.error(f"Error checking hydrogen bond: {str(e)}")
            return False

    def _find_salt_bridges(self, residues: List[Residue]) -> List[Dict[str, Any]]:
        """Find salt bridges between residues.

        Args:
            residues: List of residues

        Returns:
            List of salt bridge descriptions
        """
        try:
            salt_bridges = []
            pos_residues = [res for res in residues if res.get_resname() in ["ARG", "LYS"]]
            neg_residues = [res for res in residues if res.get_resname() in ["ASP", "GLU"]]

            for pos_res in pos_residues:
                pos_center = self._get_residue_center(pos_res)
                for neg_res in neg_residues:
                    neg_center = self._get_residue_center(neg_res)
                    distance = np.linalg.norm(pos_center - neg_center)
                    if distance < 4.0:  # Salt bridge cutoff
                        salt_bridges.append(
                            {
                                "positive_res": pos_res.get_id()[1],
                                "negative_res": neg_res.get_id()[1],
                                "distance": distance,
                            }
                        )

            return salt_bridges

        except Exception as e:
            self.logger.error(f"Error finding salt bridges: {str(e)}")
            return []

    def _find_disulfides(self, residues: List[Residue]) -> List[Dict[str, Any]]:
        """Find disulfide bonds between cysteines.

        Args:
            residues: List of residues

        Returns:
            List of disulfide bond descriptions
        """
        try:
            disulfides = []
            cysteines = [res for res in residues if res.get_resname() == "CYS"]

            for i, cys1 in enumerate(cysteines):
                if "SG" not in cys1:
                    continue
                sg1 = cys1["SG"]
                for cys2 in cysteines[i + 1 :]:
                    if "SG" not in cys2:
                        continue
                    sg2 = cys2["SG"]
                    distance = sg1 - sg2
                    if distance < 2.2:  # Typical S-S bond length
                        disulfides.append(
                            {
                                "residue1": cys1.get_id()[1],
                                "residue2": cys2.get_id()[1],
                                "distance": distance,
                            }
                        )

            return disulfides

        except Exception as e:
            self.logger.error(f"Error finding disulfide bonds: {str(e)}")
            return []

    def _calculate_shape_complementarity(
        self,
        structure: Structure,
        residues: List[int],
        ligand: Chem.Mol,
    ) -> float:
        """Calculate shape complementarity between pocket and ligand.

        Args:
            structure: BioPython Structure object
            residues: List of residue numbers
            ligand: RDKit molecule

        Returns:
            Shape complementarity score (0-1)
        """
        try:
            # Get pocket surface
            pocket_residues = [structure[0]["A"][res_num] for res_num in residues]
            pocket_coords = []
            for res in pocket_residues:
                for atom in res:
                    pocket_coords.append(atom.get_coord())
            pocket_coords = np.array(pocket_coords)

            # Get ligand surface
            conf = ligand.GetConformer()
            ligand_coords = []
            for i in range(ligand.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                ligand_coords.append([pos.x, pos.y, pos.z])
            ligand_coords = np.array(ligand_coords)

            # Calculate shape overlap using RDKit
            score = rdShapeHelpers.ShapeTanimotoDist(
                pocket_coords,
                ligand_coords,
            )

            return float(score)

        except Exception as e:
            self.logger.error(f"Error calculating shape complementarity: {str(e)}")
            return 0.0

    def _calculate_interaction_energy(
        self,
        residues: List[Residue],
        ligand: Chem.Mol,
    ) -> Dict[str, float]:
        """Calculate interaction energy between pocket and ligand.

        Args:
            residues: List of residues
            ligand: RDKit molecule

        Returns:
            Dictionary of interaction energies
        """
        try:
            # Get pocket atoms
            pocket_coords = []
            pocket_charges = []
            for res in residues:
                for atom in res:
                    pocket_coords.append(atom.get_coord())
                    # Approximate charges based on atom type
                    if atom.get_name().startswith("N"):
                        pocket_charges.append(0.5)
                    elif atom.get_name().startswith("O"):
                        pocket_charges.append(-0.5)
                    else:
                        pocket_charges.append(0.0)

            pocket_coords = np.array(pocket_coords)
            pocket_charges = np.array(pocket_charges)

            # Get ligand atoms
            conf = ligand.GetConformer()
            ligand_coords = []
            ligand_charges = []
            for i in range(ligand.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                ligand_coords.append([pos.x, pos.y, pos.z])
                # Get formal charges
                atom = ligand.GetAtomWithIdx(i)
                ligand_charges.append(float(atom.GetFormalCharge()))

            ligand_coords = np.array(ligand_coords)
            ligand_charges = np.array(ligand_charges)

            # Calculate energies
            vdw_energy = self._calculate_vdw_energy(pocket_coords, ligand_coords)
            elec_energy = self._calculate_electrostatic_energy(
                pocket_coords,
                pocket_charges,
                ligand_coords,
                ligand_charges,
            )

            return {
                "vdw_energy": vdw_energy,
                "electrostatic_energy": elec_energy,
                "total_energy": vdw_energy + elec_energy,
            }

        except Exception as e:
            self.logger.error(f"Error calculating interaction energy: {str(e)}")
            return {
                "vdw_energy": 0.0,
                "electrostatic_energy": 0.0,
                "total_energy": 0.0,
            }

    def _calculate_vdw_energy(
        self,
        coords1: np.ndarray,
        coords2: np.ndarray,
        epsilon: float = 0.1,
        sigma: float = 3.5,
    ) -> float:
        """Calculate van der Waals energy between two sets of coordinates.

        Args:
            coords1: First set of coordinates
            coords2: Second set of coordinates
            epsilon: Energy well depth
            sigma: Distance at zero energy

        Returns:
            Van der Waals energy
        """
        try:
            # Calculate all pairwise distances
            distances = cdist(coords1, coords2)

            # Calculate Lennard-Jones potential
            sr6 = (sigma / distances) ** 6
            energy = 4 * epsilon * (sr6**2 - sr6)

            return float(np.sum(energy))

        except Exception as e:
            self.logger.error(f"Error calculating van der Waals energy: {str(e)}")
            return 0.0

    def _calculate_electrostatic_energy(
        self,
        coords1: np.ndarray,
        charges1: np.ndarray,
        coords2: np.ndarray,
        charges2: np.ndarray,
        dielectric: float = 80.0,
    ) -> float:
        """Calculate electrostatic energy between charged particles.

        Args:
            coords1: First set of coordinates
            charges1: First set of charges
            coords2: Second set of coordinates
            charges2: Second set of charges
            dielectric: Dielectric constant

        Returns:
            Electrostatic energy
        """
        try:
            # Calculate all pairwise distances
            distances = cdist(coords1, coords2)

            # Calculate Coulomb energy
            charge_products = charges1[:, np.newaxis] * charges2
            energy = 332.0 * charge_products / (dielectric * distances)  # 332 kcal/mol conversion

            return float(np.sum(energy))

        except Exception as e:
            self.logger.error(f"Error calculating electrostatic energy: {str(e)}")
            return 0.0

    def _predict_binding_mode(
        self,
        residues: List[Residue],
        ligand: Chem.Mol,
    ) -> Dict[str, Any]:
        """Predict binding mode of ligand in pocket.

        Args:
            residues: List of pocket residues
            ligand: RDKit molecule

        Returns:
            Dictionary describing predicted binding mode
        """
        try:
            # Get pocket pharmacophore
            pocket_features = self._generate_pharmacophore(residues)

            # Get ligand pharmacophore
            ligand_features = self._get_ligand_pharmacophore(ligand)

            # Match pharmacophores
            matches = self._match_pharmacophores(pocket_features, ligand_features)

            # Score matches
            scored_matches = []
            for match in matches:
                score = self._score_pharmacophore_match(match)
                scored_matches.append((score, match))

            # Get best match
            if scored_matches:
                best_score, best_match = max(scored_matches, key=lambda x: x[0])
                return {
                    "score": best_score,
                    "interactions": best_match,
                    "binding_energy": self._estimate_binding_energy(best_match),
                }
            else:
                return {
                    "score": 0.0,
                    "interactions": [],
                    "binding_energy": 0.0,
                }

        except Exception as e:
            self.logger.error(f"Error predicting binding mode: {str(e)}")
            return {
                "score": 0.0,
                "interactions": [],
                "binding_energy": 0.0,
            }
