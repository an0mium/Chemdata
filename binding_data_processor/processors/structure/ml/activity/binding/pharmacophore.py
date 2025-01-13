"""Pharmacophore analysis for binding sites."""

from typing import Dict, List, Optional, Set, Tuple, Union, Any
import logging
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from sklearn.cluster import DBSCAN
from scipy.spatial.distance import cdist

from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeHelpers, ChemicalFeatures
from Bio.PDB import Structure, Residue

logger = logging.getLogger(__name__)


@dataclass
class PharmacophoreFeature:
    """Representation of a pharmacophore feature."""

    feature_type: str
    position: np.ndarray
    radius: float = 1.0
    weight: float = 1.0
    vector: Optional[np.ndarray] = None
    excluded: bool = False
    residue_id: Optional[str] = None
    score: float = 1.0


class PharmacophoreGenerator:
    """Generate and analyze pharmacophore features for binding prediction."""

    def __init__(
        self,
        feature_radius: float = 1.0,
        clustering_eps: float = 2.0,
        min_cluster_size: int = 2,
        consensus_threshold: float = 0.5,
        cache_dir: Optional[Path] = None,
    ):
        """Initialize pharmacophore generator.

        Args:
            feature_radius: Default radius for pharmacophore features
            clustering_eps: DBSCAN clustering distance threshold
            min_cluster_size: Minimum points for DBSCAN cluster
            consensus_threshold: Fraction of structures needed for consensus
            cache_dir: Directory to cache pharmacophore models
        """
        self.feature_radius = feature_radius
        self.clustering_eps = clustering_eps
        self.min_cluster_size = min_cluster_size
        self.consensus_threshold = consensus_threshold
        self.cache_dir = cache_dir

        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Initialize feature factory
        self.factory = ChemicalFeatures.BuildFeatureFactory()

        # Feature definitions
        self.feature_types = {
            "hydrophobic": self._is_hydrophobic,
            "hbond_donor": self._is_hbond_donor,
            "hbond_acceptor": self._is_hbond_acceptor,
            "aromatic": self._is_aromatic,
            "charged_pos": self._is_positively_charged,
            "charged_neg": self._is_negatively_charged,
            "metal": self._is_metal_binding,
        }

        # QSAR model parameters
        self.qsar_descriptors = {
            "hydrophobic_score": self._calc_hydrophobic_score,
            "hbond_score": self._calc_hbond_score,
            "shape_score": self._calc_shape_score,
            "electrostatic_score": self._calc_electrostatic_score,
        }

    def get_features(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Get pharmacophore features for a molecule.

        Args:
            mol: Input molecule

        Returns:
            Dictionary of pharmacophore features
        """
        try:
            # Generate 3D conformer if needed
            if not mol.GetNumConformers():
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol)

            # Get chemical features
            features = {}
            feats = self.factory.GetFeaturesForMol(mol)

            # Count feature types
            feature_counts = {}
            for feat in feats:
                feat_type = feat.GetFamily()
                feature_counts[feat_type] = feature_counts.get(feat_type, 0) + 1

            features["counts"] = feature_counts

            # Get feature positions
            feature_positions = {}
            for feat in feats:
                feat_type = feat.GetFamily()
                if feat_type not in feature_positions:
                    feature_positions[feat_type] = []
                feature_positions[feat_type].append(feat.GetPos())

            features["positions"] = feature_positions

            # Calculate feature distances
            features["distances"] = self._calculate_feature_distances(feats)

            # Calculate feature angles
            features["angles"] = self._calculate_feature_angles(feats)

            return features

        except Exception as e:
            logger.error(f"Error getting pharmacophore features: {str(e)}")
            return {}

    def generate_pharmacophore(
        self,
        structure: Structure,
        binding_site_residues: Optional[List[Residue]] = None,
        include_vectors: bool = True,
        include_excluded_volumes: bool = True,
        min_feature_score: float = 0.5,
    ) -> List[PharmacophoreFeature]:
        """Generate pharmacophore model from binding site.

        Args:
            structure: Protein structure
            binding_site_residues: List of binding site residues, or all if None
            include_vectors: Whether to include directional vectors
            include_excluded_volumes: Whether to add excluded volume spheres
            min_feature_score: Minimum score threshold for features

        Returns:
            List of pharmacophore features
        """
        features = []
        residues = binding_site_residues or []

        # Generate basic features
        for residue in residues:
            coords = self._get_centroid(residue)
            if coords is None:
                continue

            for feat_type, check_func in self.feature_types.items():
                if check_func(residue):
                    score = self._score_feature(residue, feat_type)
                    if score >= min_feature_score:
                        vector = self._get_feature_vector(residue, feat_type) if include_vectors else None
                        features.append(
                            PharmacophoreFeature(feature_type=feat_type, position=coords, radius=self.feature_radius, vector=vector, residue_id=str(residue.id[1]), score=score)
                        )

        # Add excluded volumes
        if include_excluded_volumes:
            excluded = self._generate_excluded_volumes(structure, residues)
            features.extend(excluded)

        # Cluster similar features
        features = self._cluster_features(features)

        return features

    def generate_consensus_pharmacophore(
        self,
        structures: List[Structure],
        binding_sites: List[List[Residue]],
        min_fraction: Optional[float] = None,
    ) -> List[PharmacophoreFeature]:
        """Generate consensus pharmacophore from multiple structures.

        Args:
            structures: List of protein structures
            binding_sites: List of binding site residues for each structure
            min_fraction: Minimum fraction of structures needed (overrides default)

        Returns:
            List of consensus pharmacophore features
        """
        min_fraction = min_fraction or self.consensus_threshold
        all_features = []

        # Generate features for each structure
        for struct, site in zip(structures, binding_sites):
            features = self.generate_pharmacophore(struct, site)
            all_features.append(features)

        # Find consensus features
        consensus = []
        ref_features = all_features[0]

        for ref_feat in ref_features:
            matching = [1]  # Count reference structure

            # Look for matching features in other structures
            for other_feats in all_features[1:]:
                for feat in other_feats:
                    if self._features_match(ref_feat, feat):
                        matching.append(1)
                        break

            # Add to consensus if enough structures have this feature
            if len(matching) / len(structures) >= min_fraction:
                consensus.append(ref_feat)

        return consensus

    def get_pocket_features(
        self,
        residues: List[int],
        properties: Dict[str, Any],
    ) -> List[float]:
        """Get pharmacophore features for binding pocket residues.

        Args:
            residues: List of residue numbers
            properties: Structure properties

        Returns:
            List of pharmacophore feature scores
        """
        try:
            features = []

            # Get residue properties
            for res in residues:
                if str(res) in properties.get("residue_properties", {}):
                    props = properties["residue_properties"][str(res)]

                    # Add pharmacophore-relevant properties
                    features.extend(
                        [
                            float(props.get("hydrophobicity", 0.0)),
                            float(props.get("aromaticity", 0.0)),
                            float(props.get("hbond_donor", 0.0)),
                            float(props.get("hbond_acceptor", 0.0)),
                            float(props.get("charge", 0.0)),
                        ]
                    )

            if not features:
                return [0.0] * 5  # Default feature vector

            return features

        except Exception as e:
            logger.error(f"Error getting pocket features: {str(e)}")
            return [0.0] * 5

    def score_features(self, features: List[float]) -> float:
        """Score pharmacophore features.

        Args:
            features: List of pharmacophore feature scores

        Returns:
            Combined feature score (0-1)
        """
        try:
            if not features:
                return 0.0

            # Normalize features
            features = np.array(features)
            features = np.clip(features, 0, 1)

            # Weight and combine features
            weights = [0.3, 0.2, 0.2, 0.2, 0.1]  # Importance weights
            score = np.sum(features * weights) / sum(weights)

            return float(score)

        except Exception as e:
            logger.error(f"Error scoring features: {str(e)}")
            return 0.0

    def score_pharmacophore_match(
        self,
        pharmacophore: List[PharmacophoreFeature],
        molecule: Chem.Mol,
        conformer_id: int = -1,
    ) -> float:
        """Score how well a molecule matches a pharmacophore.

        Args:
            pharmacophore: List of pharmacophore features
            molecule: Molecule to score
            conformer_id: Which conformer to use

        Returns:
            Match score between 0 and 1
        """
        if not pharmacophore:
            return 0.0

        # Get molecule features
        mol_features = self._get_molecule_features(molecule, conformer_id)

        # Score each pharmacophore feature
        scores = []
        for pharm_feat in pharmacophore:
            feat_scores = []

            for mol_feat in mol_features:
                if mol_feat.feature_type == pharm_feat.feature_type:
                    # Distance score
                    dist = np.linalg.norm(mol_feat.position - pharm_feat.position)
                    dist_score = np.exp(-dist / pharm_feat.radius)

                    # Vector score if applicable
                    vec_score = 1.0
                    if pharm_feat.vector is not None and mol_feat.vector is not None:
                        cos_angle = np.dot(pharm_feat.vector, mol_feat.vector)
                        vec_score = (cos_angle + 1) / 2

                    feat_scores.append(dist_score * vec_score)

            # Take best matching molecular feature
            if feat_scores:
                scores.append(max(feat_scores) * pharm_feat.weight)
            else:
                scores.append(0.0)

        # Combine scores
        if not scores:
            return 0.0
        return np.mean(scores)

    def _calculate_feature_distances(
        self,
        features: List[Any],
    ) -> Dict[str, List[float]]:
        """Calculate distances between pharmacophore features.

        Args:
            features: List of chemical features

        Returns:
            Dictionary of feature pair distances
        """
        distances = {}
        for i, feat1 in enumerate(features):
            for j, feat2 in enumerate(features[i + 1 :], i + 1):
                key = f"{feat1.GetFamily()}-{feat2.GetFamily()}"
                if key not in distances:
                    distances[key] = []
                pos1 = feat1.GetPos()
                pos2 = feat2.GetPos()
                dist = np.linalg.norm(np.array(pos1) - np.array(pos2))
                distances[key].append(float(dist))
        return distances

    def _calculate_feature_angles(
        self,
        features: List[Any],
    ) -> Dict[str, List[float]]:
        """Calculate angles between pharmacophore feature triplets.

        Args:
            features: List of chemical features

        Returns:
            Dictionary of feature triplet angles
        """
        angles = {}
        for i, feat1 in enumerate(features):
            for j, feat2 in enumerate(features[i + 1 :], i + 1):
                for k, feat3 in enumerate(features[j + 1 :], j + 1):
                    key = f"{feat1.GetFamily()}-{feat2.GetFamily()}-{feat3.GetFamily()}"
                    if key not in angles:
                        angles[key] = []

                    # Calculate angle between vectors
                    pos1 = np.array(feat1.GetPos())
                    pos2 = np.array(feat2.GetPos())
                    pos3 = np.array(feat3.GetPos())

                    v1 = pos2 - pos1
                    v2 = pos3 - pos2

                    # Normalize vectors
                    v1_norm = np.linalg.norm(v1)
                    v2_norm = np.linalg.norm(v2)
                    if v1_norm > 0 and v2_norm > 0:
                        v1 = v1 / v1_norm
                        v2 = v2 / v2_norm
                        angle = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
                        angles[key].append(float(np.degrees(angle)))

        return angles

    def _cluster_features(
        self,
        features: List[PharmacophoreFeature],
    ) -> List[PharmacophoreFeature]:
        """Cluster similar pharmacophore features."""
        if not features:
            return []

        # Group by feature type
        by_type = {}
        for feat in features:
            if feat.feature_type not in by_type:
                by_type[feat.feature_type] = []
            by_type[feat.feature_type].append(feat)

        # Cluster each type
        clustered = []
        for feat_type, type_features in by_type.items():
            if len(type_features) < self.min_cluster_size:
                clustered.extend(type_features)
                continue

            # Get coordinates for clustering
            coords = np.array([f.position for f in type_features])

            # Run DBSCAN
            clustering = DBSCAN(eps=self.clustering_eps, min_samples=self.min_cluster_size).fit(coords)

            # Create consensus features for each cluster
            for cluster_id in set(clustering.labels_):
                if cluster_id == -1:  # Noise points
                    continue

                # Get features in this cluster
                cluster_feats = [f for i, f in enumerate(type_features) if clustering.labels_[i] == cluster_id]

                # Create consensus feature
                avg_pos = np.mean([f.position for f in cluster_feats], axis=0)
                avg_score = np.mean([f.score for f in cluster_feats])

                clustered.append(PharmacophoreFeature(feature_type=feat_type, position=avg_pos, score=avg_score, weight=len(cluster_feats) / len(type_features)))

        return clustered

    def _generate_excluded_volumes(
        self,
        structure: Structure,
        binding_site_residues: List[Residue],
    ) -> List[PharmacophoreFeature]:
        """Generate excluded volume spheres."""
        excluded = []

        # Get binding site atoms
        binding_atoms = set()
        for res in binding_site_residues:
            binding_atoms.update(res.get_atoms())

        # Find cavities
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue not in binding_site_residues:
                        # Check if atoms are near binding site
                        for atom in residue.get_atoms():
                            near_binding = False
                            for bind_atom in binding_atoms:
                                dist = np.linalg.norm(atom.get_coord() - bind_atom.get_coord())
                                if dist < 4.0:  # Cutoff distance
                                    near_binding = True
                                    break

                            if near_binding:
                                excluded.append(PharmacophoreFeature(feature_type="excluded", position=atom.get_coord(), radius=2.0, excluded=True))

        return excluded

    def _get_molecule_features(
        self,
        molecule: Chem.Mol,
        conformer_id: int = -1,
    ) -> List[PharmacophoreFeature]:
        """Extract pharmacophore features from molecule."""
        features = []
        conf = molecule.GetConformer(conformer_id)

        # Add donors
        for idx in AllChem.FindAllHBondDonorAtoms(molecule):
            pos = conf.GetAtomPosition(idx)
            features.append(PharmacophoreFeature(feature_type="hbond_donor", position=np.array([pos.x, pos.y, pos.z])))

        # Add acceptors
        for idx in AllChem.FindAllHBondAcceptorAtoms(molecule):
            pos = conf.GetAtomPosition(idx)
            features.append(PharmacophoreFeature(feature_type="hbond_acceptor", position=np.array([pos.x, pos.y, pos.z])))

        # Add aromatic
        for idx in self._find_aromatic_atoms(molecule):
            pos = conf.GetAtomPosition(idx)
            features.append(PharmacophoreFeature(feature_type="aromatic", position=np.array([pos.x, pos.y, pos.z])))

        return features

    # Feature type checks
    def _is_hydrophobic(self, residue: Residue) -> bool:
        """Check if residue is hydrophobic."""
        hydrophobic = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
        return residue.get_resname() in hydrophobic

    def _is_hbond_donor(self, residue: Residue) -> bool:
        """Check if residue can be H-bond donor."""
        donors = {"LYS", "ARG", "HIS", "SER", "THR", "ASN", "GLN", "TYR", "TRP"}
        return residue.get_resname() in donors

    def _is_hbond_acceptor(self, residue: Residue) -> bool:
        """Check if residue can be H-bond acceptor."""
        acceptors = {"ASP", "GLU", "ASN", "GLN", "HIS", "SER", "THR", "TYR"}
        return residue.get_resname() in acceptors

    def _is_aromatic(self, residue: Residue) -> bool:
        """Check if residue has aromatic ring."""
        aromatic = {"PHE", "TYR", "TRP", "HIS"}
        return residue.get_resname() in aromatic

    def _is_positively_charged(self, residue: Residue) -> bool:
        """Check if residue is positively charged."""
        positive = {"LYS", "ARG", "HIS"}
        return residue.get_resname() in positive

    def _is_negatively_charged(self, residue: Residue) -> bool:
        """Check if residue is negatively charged."""
        negative = {"ASP", "GLU"}
        return residue.get_resname() in negative

    def _is_metal_binding(self, residue: Residue) -> bool:
        """Check if residue can bind metals."""
        metal_binding = {"HIS", "CYS", "ASP", "GLU"}
        return residue.get_resname() in metal_binding

    # QSAR descriptors
    def _calc_hydrophobic_score(self, residue: Residue) -> float:
        """Calculate hydrophobicity score."""
        # Normalized hydrophobicity values
        hydrophobicity = {
            "ILE": 1.0,
            "VAL": 0.9,
            "LEU": 0.9,
            "PHE": 0.8,
            "MET": 0.7,
            "ALA": 0.6,
            "TRP": 0.5,
            "GLY": 0.4,
            "PRO": 0.3,
            "THR": 0.2,
            "SER": 0.1,
            "TYR": 0.1,
            "HIS": 0.1,
            "GLN": 0.0,
            "ASN": 0.0,
            "GLU": 0.0,
            "LYS": 0.0,
            "ASP": 0.0,
            "ARG": 0.0,
        }
        return hydrophobicity.get(residue.get_resname(), 0.0)

    def _calc_hbond_score(self, residue: Residue) -> float:
        """Calculate H-bond potential score."""
        hbond_potential = {
            "SER": 1.0,
            "THR": 1.0,
            "ASN": 1.0,
            "GLN": 1.0,
            "TYR": 0.8,
            "TRP": 0.7,
            "HIS": 0.7,
            "LYS": 0.6,
            "ARG": 0.6,
            "ASP": 0.5,
            "GLU": 0.5,
        }
        return hbond_potential.get(residue.get_resname(), 0.0)

    def _calc_shape_score(self, residue: Residue) -> float:
        """Calculate shape importance score."""
        shape_importance = {
            "TRP": 1.0,
            "PHE": 0.9,
            "TYR": 0.9,
            "HIS": 0.8,
            "ARG": 0.7,
            "LYS": 0.6,
            "MET": 0.5,
            "LEU": 0.4,
            "ILE": 0.4,
            "GLN": 0.3,
            "GLU": 0.3,
            "ASN": 0.2,
            "ASP": 0.2,
            "PRO": 0.2,
            "VAL": 0.1,
            "THR": 0.1,
            "ALA": 0.0,
            "GLY": 0.0,
        }
        return shape_importance.get(residue.get_resname(), 0.0)

    def _calc_electrostatic_score(self, residue: Residue) -> float:
        """Calculate electrostatic importance score."""
        electrostatic = {
            "ARG": 1.0,
            "LYS": 1.0,
            "ASP": 1.0,
            "GLU": 1.0,
            "HIS": 0.7,
            "TYR": 0.3,
            "THR": 0.2,
            "SER": 0.2,
            "ASN": 0.2,
            "GLN": 0.2,
        }
        return electrostatic.get(residue.get_resname(), 0.0)

    def _get_centroid(self, residue: Residue) -> Optional[np.ndarray]:
        """Get centroid coordinates of residue."""
        coords = []
        for atom in residue:
            coords.append(atom.get_coord())
        if not coords:
            return None
        return np.mean(coords, axis=0)

    def _get_feature_vector(
        self,
        residue: Residue,
        feature_type: str,
    ) -> Optional[np.ndarray]:
        """Get directional vector for a feature."""
        if feature_type == "hbond_donor":
            # Point from heavy atom to H
            return self._get_hbond_vector(residue, donor=True)
        elif feature_type == "hbond_acceptor":
            # Point from heavy atom with lone pair
            return self._get_hbond_vector(residue, donor=False)
        elif feature_type == "aromatic":
            # Normal to ring plane
            return self._get_aromatic_vector(residue)
        return None

    def _get_hbond_vector(
        self,
        residue: Residue,
        donor: bool,
    ) -> Optional[np.ndarray]:
        """Get H-bond directional vector."""
        # This would use atom positions to determine vector
        # Placeholder implementation
        return np.array([0, 0, 1])

    def _get_aromatic_vector(
        self,
        residue: Residue,
    ) -> Optional[np.ndarray]:
        """Get aromatic ring normal vector."""
        # This would calculate ring plane normal
        # Placeholder implementation
        return np.array([0, 0, 1])

    def _find_aromatic_atoms(self, mol: Chem.Mol) -> Set[int]:
        """Find aromatic atoms in molecule."""
        aromatic = set()
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic():
                aromatic.add(atom.GetIdx())
        return aromatic

    def _features_match(
        self,
        feat1: PharmacophoreFeature,
        feat2: PharmacophoreFeature,
    ) -> bool:
        """Check if two features match."""
        if feat1.feature_type != feat2.feature_type:
            return False

        dist = np.linalg.norm(feat1.position - feat2.position)
        if dist > (feat1.radius + feat2.radius):
            return False

        if feat1.vector is not None and feat2.vector is not None:
            cos_angle = np.dot(feat1.vector, feat2.vector)
            if cos_angle < 0.7:  # About 45 degrees
                return False

        return True
