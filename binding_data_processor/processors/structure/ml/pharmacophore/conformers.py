"""Conformer generation and optimization for pharmacophore analysis.

This module provides functionality for:
1. 3D conformer generation using RDKit
2. Energy minimization and geometry optimization 
3. Conformer sampling, clustering and analysis
4. Conformer validation and filtering
5. Torsion analysis and fingerprinting
6. Ensemble management and optimization
"""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints

logger = logging.getLogger(__name__)


class ConformerGenerator:
    """Generate and optimize 3D conformers for molecules."""

    def __init__(
        self,
        num_confs: int = 50,
        max_attempts: int = 200,
        prune_rms_thresh: float = 0.5,
        energy_window: float = 10.0,
        random_seed: Optional[int] = None,
    ):
        """Initialize conformer generator.

        Args:
            num_confs: Number of conformers to generate
            max_attempts: Maximum number of attempts to generate conformers
            prune_rms_thresh: RMSD threshold for pruning similar conformers
            energy_window: Energy window (kcal/mol) for keeping conformers
            random_seed: Random seed for reproducibility
        """
        self.num_confs = num_confs
        self.max_attempts = max_attempts
        self.prune_rms_thresh = prune_rms_thresh
        self.energy_window = energy_window
        self.random_seed = random_seed

    def generate(
        self,
        mol: Union[str, Chem.Mol],
        optimize: bool = True,
        return_energies: bool = False,
    ) -> Union[Chem.Mol, Tuple[Chem.Mol, np.ndarray]]:
        """Generate 3D conformers for a molecule.

        Args:
            mol: Input molecule (SMILES string or RDKit mol)
            optimize: Whether to optimize conformer geometries
            return_energies: Whether to return conformer energies

        Returns:
            Molecule with conformers and optionally their energies
        """
        # Convert SMILES to mol if needed
        if isinstance(mol, str):
            mol = Chem.MolFromSmiles(mol)
            if mol is None:
                raise ValueError("Invalid SMILES string")

        # Prepare molecule
        mol = Chem.AddHs(mol)

        # Generate conformers
        conf_ids = []
        for attempt in range(self.max_attempts):
            try:
                # Use random seed if provided
                seed = self.random_seed + attempt if self.random_seed is not None else None

                # Generate conformers using ETKDGv3
                params = AllChem.ETKDGv3()
                if seed is not None:
                    params.randomSeed = seed
                params.pruneRmsThresh = self.prune_rms_thresh
                params.enforceChirality = True
                params.useExpTorsionAnglePrefs = True
                params.useBasicKnowledge = True
                params.ETversion = 2
                params.maxIterations = 1000

                conf_ids = AllChem.EmbedMultipleConfs(
                    mol,
                    numConfs=self.num_confs,
                    params=params,
                    clearConfs=True,
                    numThreads=0,  # Use all available cores
                )

                if len(conf_ids) > 0:
                    break

            except Exception as e:
                logger.warning(f"Conformer generation attempt {attempt} failed: {str(e)}")

        if len(conf_ids) == 0:
            raise RuntimeError(f"Failed to generate conformers after {self.max_attempts} attempts")

        # Optimize geometries if requested
        if optimize:
            self._optimize_conformers(mol)

        # Calculate energies and filter
        energies = self._get_conformer_energies(mol)
        mol = self._filter_conformers(mol, energies)

        if return_energies:
            return mol, np.array(energies)
        return mol

    def _optimize_conformers(self, mol: Chem.Mol) -> None:
        """Optimize conformers using MMFF94s force field.

        Args:
            mol: Molecule with conformers
        """
        for conf_id in range(mol.GetNumConformers()):
            try:
                # Optimize using MMFF94s
                AllChem.MMFFOptimizeMolecule(mol, confId=conf_id, maxIters=500)
            except Exception as e:
                logger.warning(f"Failed to optimize conformer {conf_id}: {str(e)}")

    def _get_conformer_energies(self, mol: Chem.Mol) -> List[float]:
        """Get MMFF94s energies for conformers.

        Args:
            mol: Molecule with conformers

        Returns:
            List of conformer energies
        """
        energies = []
        for conf_id in range(mol.GetNumConformers()):
            try:
                # Get MMFF94s energy
                props = AllChem.MMFFGetMoleculeProperties(mol)
                ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
                if ff is None:
                    raise ValueError("Failed to get force field")
                energy = ff.CalcEnergy()
                energies.append(energy)
            except Exception as e:
                logger.warning(f"Failed to get energy for conformer {conf_id}: {str(e)}")
                energies.append(float("inf"))
        return energies

    def _filter_conformers(
        self,
        mol: Chem.Mol,
        energies: List[float],
    ) -> Chem.Mol:
        """Filter conformers by energy and RMSD.

        Args:
            mol: Molecule with conformers
            energies: Conformer energies

        Returns:
            Filtered molecule
        """
        if not energies:
            return mol

        # Sort conformers by energy
        conf_energies = list(zip(range(len(energies)), energies))
        conf_energies.sort(key=lambda x: x[1])

        # Keep conformers within energy window
        min_energy = conf_energies[0][1]
        kept_conf_ids = []
        for conf_id, energy in conf_energies:
            if energy - min_energy <= self.energy_window:
                kept_conf_ids.append(conf_id)

        # Create new molecule with filtered conformers
        new_mol = Chem.Mol(mol)
        new_mol.RemoveAllConformers()
        for conf_id in kept_conf_ids:
            conf = mol.GetConformer(conf_id)
            new_mol.AddConformer(conf, assignId=True)

        return new_mol

    def cluster_conformers(
        self,
        mol: Chem.Mol,
        n_clusters: int = 5,
        cutoff: float = 2.0,
    ) -> Tuple[List[int], List[List[int]]]:
        """Cluster conformers by RMSD.

        Args:
            mol: Molecule with conformers
            n_clusters: Number of clusters
            cutoff: RMSD cutoff for clustering

        Returns:
            Tuple of (representative conformer IDs, cluster assignments)
        """
        if mol.GetNumConformers() == 0:
            return [], []

        # Calculate pairwise RMSD matrix
        rmsd_matrix = self.get_conformer_rmsd(mol)

        # Cluster using average linkage
        from sklearn.cluster import AgglomerativeClustering

        clustering = AgglomerativeClustering(
            n_clusters=min(n_clusters, mol.GetNumConformers()),
            affinity="precomputed",
            linkage="complete",
        )
        labels = clustering.fit_predict(rmsd_matrix)

        # Get cluster assignments
        clusters = [[] for _ in range(max(labels) + 1)]
        for i, label in enumerate(labels):
            clusters[label].append(i)

        # Get representative conformers (closest to cluster center)
        representatives = []
        for cluster in clusters:
            if not cluster:
                continue
            # Get mean RMSD to all other conformers in cluster
            mean_rmsds = []
            for conf_id in cluster:
                cluster_rmsds = [rmsd_matrix[conf_id, j] for j in cluster if j != conf_id]
                mean_rmsds.append(np.mean(cluster_rmsds) if cluster_rmsds else 0)
            # Select conformer with minimum mean RMSD
            rep_idx = np.argmin(mean_rmsds)
            representatives.append(cluster[rep_idx])

        return representatives, clusters

    def get_conformer_rmsd(self, mol: Chem.Mol) -> np.ndarray:
        """Calculate RMSD matrix between all conformer pairs.

        Args:
            mol: Molecule with conformers

        Returns:
            NxN RMSD matrix
        """
        n_confs = mol.GetNumConformers()
        rmsd_matrix = np.zeros((n_confs, n_confs))

        for i in range(n_confs):
            for j in range(i + 1, n_confs):
                try:
                    rmsd = AllChem.GetBestRMS(mol, mol, i, j)
                    rmsd_matrix[i, j] = rmsd
                    rmsd_matrix[j, i] = rmsd
                except Exception as e:
                    logger.warning(f"Failed to calculate RMSD between conformers {i} and {j}: {str(e)}")
                    rmsd_matrix[i, j] = float("inf")
                    rmsd_matrix[j, i] = float("inf")

        return rmsd_matrix

    def get_conformer_torsions(self, mol: Chem.Mol) -> List[np.ndarray]:
        """Get torsion fingerprints for conformers.

        Args:
            mol: Molecule with conformers

        Returns:
            List of torsion fingerprint arrays
        """
        torsions = []
        for conf_id in range(mol.GetNumConformers()):
            try:
                # Get torsion fingerprint
                torsion_fp = TorsionFingerprints.GetTFDMatrix(
                    [mol],
                    useWeights=True,
                    maxDev="equal",
                    confId1=[conf_id],
                )[0]
                torsions.append(torsion_fp)
            except Exception as e:
                logger.warning(f"Failed to get torsions for conformer {conf_id}: {str(e)}")
                torsions.append(np.array([]))
        return torsions

    def align_conformers(self, mol: Chem.Mol, ref_conf_id: int = 0) -> None:
        """Align all conformers to a reference conformer.

        Args:
            mol: Molecule with conformers
            ref_conf_id: Reference conformer ID
        """
        try:
            AllChem.AlignMolConformers(mol, confId=ref_conf_id)
        except Exception as e:
            logger.warning(f"Failed to align conformers: {str(e)}")

    def get_conformer_coordinates(self, mol: Chem.Mol, conf_id: int = 0) -> np.ndarray:
        """Get atomic coordinates for a conformer.

        Args:
            mol: Molecule with conformers
            conf_id: Conformer ID

        Returns:
            Nx3 array of atomic coordinates
        """
        try:
            conf = mol.GetConformer(conf_id)
            coords = []
            for i in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                coords.append([pos.x, pos.y, pos.z])
            return np.array(coords)
        except Exception as e:
            logger.error(f"Failed to get conformer coordinates: {str(e)}")
            return np.array([])

    def get_conformer_features(
        self,
        mol: Chem.Mol,
        conf_id: int = 0,
    ) -> Dict[str, np.ndarray]:
        """Get geometric features for a conformer.

        Args:
            mol: Molecule with conformers
            conf_id: Conformer ID

        Returns:
            Dictionary of feature arrays including:
            - center_of_mass: Center of mass coordinates
            - principal_axes: Principal axes of the molecule
            - moments_of_inertia: Principal moments of inertia
            - radius_of_gyration: Radius of gyration
        """
        try:
            coords = self.get_conformer_coordinates(mol, conf_id)
            if len(coords) == 0:
                return {}

            features = {}

            # Center of mass
            features["center_of_mass"] = np.mean(coords, axis=0)

            # Principal axes and moments of inertia
            coords_centered = coords - features["center_of_mass"]
            cov = np.dot(coords_centered.T, coords_centered)
            eigenvals, eigenvecs = np.linalg.eigh(cov)
            features["principal_axes"] = eigenvecs.T
            features["moments_of_inertia"] = eigenvals

            # Radius of gyration
            features["radius_of_gyration"] = np.sqrt(np.mean(np.sum(coords_centered**2, axis=1)))

            return features

        except Exception as e:
            logger.error(f"Failed to calculate conformer features: {str(e)}")
            return {}
