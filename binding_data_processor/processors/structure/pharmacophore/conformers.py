"""Conformer generation and optimization for pharmacophore analysis."""

from typing import Dict, List, Optional, Tuple, Union
import logging
from dataclasses import dataclass

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit.Chem.rdMolTransforms import GetDihedralDeg, SetDihedralDeg
from rdkit.Chem.MolStandardize import rdMolStandardize


@dataclass
class ConformerStats:
    """Statistics for generated conformers."""

    num_conformers: int
    energy_min: float
    energy_max: float
    energy_mean: float
    energy_std: float
    rmsd_min: float
    rmsd_max: float
    rmsd_mean: float
    rmsd_std: float


class ConformerGenerator:
    """Enhanced class for generating and optimizing 3D conformers."""

    def __init__(
        self,
        max_conformers: int = 50,
        energy_window: float = 10.0,
        rmsd_threshold: float = 0.5,
        force_field: str = "MMFF94s",
        random_seed: int = 42,
    ):
        """Initialize the conformer generator with enhanced options.

        Args:
            max_conformers: Maximum number of conformers to generate
            energy_window: Energy window in kcal/mol for pruning conformers
            rmsd_threshold: RMSD threshold for clustering conformers
            force_field: Force field to use ("MMFF94s", "MMFF94", or "UFF")
            random_seed: Random seed for reproducibility
        """
        self.max_conformers = max_conformers
        self.energy_window = energy_window
        self.rmsd_threshold = rmsd_threshold
        self.force_field = force_field
        self.random_seed = random_seed
        self.logger = logging.getLogger(__name__)

        # Enhanced ETKDGv3 parameters
        self.etkdg_params = {
            "useExpTorsionAnglePrefs": True,  # Use experimental torsion angle preferences
            "useBasicKnowledge": True,  # Use basic chemical knowledge
            "enforceChirality": True,  # Enforce input molecule chirality
            "useSmallRingTorsions": True,  # Enhanced small ring conformations
            "useMacrocycleTorsions": True,  # Enhanced macrocycle handling
            "ETversion": 3,  # Use ETKDGv3
            "pruneRmsThresh": rmsd_threshold,  # RMSD threshold for pruning
            "numThreads": 0,  # Use all available CPUs
        }

        # Force field parameters
        self.ff_params = {
            "maxIters": 500,  # Maximum optimization iterations
            "nonBondedThresh": 100.0,  # Non-bonded interactions threshold
            "vdwThresh": 10.0,  # Van der Waals interactions threshold
            "ignoreInterfragInteractions": False,  # Consider fragment interactions
        }

    def prepare_molecule(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Prepare molecule for conformer generation.

        Args:
            mol: Input molecule

        Returns:
            Prepared molecule or None if preparation fails
        """
        if mol is None:
            return None

        try:
            # Make a copy
            mol = Chem.Mol(mol)

            # Clean up structure
            mol = rdMolStandardize.Cleanup(mol)
            mol = rdMolStandardize.Normalize(mol)

            # Add hydrogens
            mol = Chem.AddHs(mol)

            return mol

        except Exception as e:
            self.logger.error(f"Error preparing molecule: {str(e)}")
            return None

    def generate(
        self,
        mol: Chem.Mol,
        num_conformers: Optional[int] = None,
        optimize: bool = True,
        cluster: bool = True,
    ) -> Optional[Chem.Mol]:
        """Generate 3D conformers with enhanced options.

        Args:
            mol: Input molecule
            num_conformers: Number of conformers to generate (default: self.max_conformers)
            optimize: Whether to optimize conformers
            cluster: Whether to cluster similar conformers

        Returns:
            Molecule with generated conformers
        """
        mol = self.prepare_molecule(mol)
        if mol is None:
            return None

        try:
            # Set up ETKDG parameters
            params = AllChem.ETKDGv3()
            for key, value in self.etkdg_params.items():
                if hasattr(params, key):
                    setattr(params, key, value)
            params.randomSeed = self.random_seed

            # Generate conformers
            n = num_conformers or self.max_conformers
            cids = AllChem.EmbedMultipleConfs(
                mol,
                numConfs=n,
                params=params,
                randomSeed=self.random_seed,
                numThreads=0,
            )

            if len(cids) == 0:
                self.logger.warning("Failed to generate conformers with ETKDG, trying distance geometry")
                cids = AllChem.EmbedMultipleConfs(
                    mol,
                    numConfs=n,
                    randomSeed=self.random_seed,
                    numThreads=0,
                )

            if len(cids) == 0:
                self.logger.error("Failed to generate any conformers")
                return None

            # Optimize conformers
            if optimize:
                self.optimize(mol)

            # Cluster conformers
            if cluster:
                mol = self.cluster_conformers(mol)

            return mol

        except Exception as e:
            self.logger.error(f"Error generating conformers: {str(e)}")
            return None

    def optimize(self, mol: Chem.Mol, max_iters: Optional[int] = None) -> bool:
        """Optimize conformers using selected force field.

        Args:
            mol: Molecule with conformers
            max_iters: Maximum optimization iterations (overrides default)

        Returns:
            Whether optimization succeeded
        """
        if not mol or not mol.GetNumConformers():
            return False

        try:
            # Set up force field parameters
            ff_params = self.ff_params.copy()
            if max_iters is not None:
                ff_params["maxIters"] = max_iters

            # Optimize using selected force field
            if self.force_field == "MMFF94s":
                results = AllChem.MMFFOptimizeMoleculeConfs(
                    mol,
                    mmffVariant="MMFF94s",
                    **ff_params,
                )
            elif self.force_field == "MMFF94":
                results = AllChem.MMFFOptimizeMoleculeConfs(
                    mol,
                    **ff_params,
                )
            else:  # UFF
                results = AllChem.UFFOptimizeMoleculeConfs(
                    mol,
                    **ff_params,
                )

            # Check results
            return all(res[0] == 0 for res in results)

        except Exception as e:
            self.logger.error(f"Error optimizing conformers: {str(e)}")
            return False

    def cluster_conformers(self, mol: Chem.Mol) -> Chem.Mol:
        """Cluster conformers based on RMSD and energy.

        Args:
            mol: Molecule with conformers

        Returns:
            Molecule with clustered conformers
        """
        if not mol or mol.GetNumConformers() <= 1:
            return mol

        try:
            # Calculate energies
            energies = self.get_conformer_energies(mol)
            if not energies:
                return mol

            # Sort by energy
            energy_sorted = sorted(enumerate(energies), key=lambda x: x[1])
            min_energy = energy_sorted[0][1]

            # Keep conformers within energy window
            keep_conf_ids = []
            for conf_id, energy in energy_sorted:
                if energy - min_energy <= self.energy_window:
                    keep_conf_ids.append(conf_id)

            # Calculate RMSD matrix for remaining conformers
            rmsd_clusters = []
            for i, conf_id in enumerate(keep_conf_ids):
                if not any(conf_id in cluster for cluster in rmsd_clusters):
                    # Start new cluster
                    cluster = {conf_id}
                    for j in range(i + 1, len(keep_conf_ids)):
                        other_id = keep_conf_ids[j]
                        if other_id not in cluster:
                            rmsd = AllChem.GetConformerRMS(mol, conf_id, other_id)
                            if rmsd < self.rmsd_threshold:
                                cluster.add(other_id)
                    rmsd_clusters.append(cluster)

            # Keep lowest energy conformer from each cluster
            final_conf_ids = []
            for cluster in rmsd_clusters:
                cluster_energies = [(cid, energies[cid]) for cid in cluster]
                best_conf_id = min(cluster_energies, key=lambda x: x[1])[0]
                final_conf_ids.append(best_conf_id)

            # Create new molecule with kept conformers
            new_mol = Chem.Mol(mol)
            new_mol.RemoveAllConformers()
            for conf_id in sorted(final_conf_ids):
                conf = mol.GetConformer(conf_id)
                new_mol.AddConformer(conf, assignId=True)

            return new_mol

        except Exception as e:
            self.logger.error(f"Error clustering conformers: {str(e)}")
            return mol

    def get_conformer_energies(self, mol: Chem.Mol) -> List[float]:
        """Get energies for all conformers using selected force field.

        Args:
            mol: Molecule with conformers

        Returns:
            List of energies in kcal/mol
        """
        if not mol or not mol.GetNumConformers():
            return []

        try:
            energies = []
            for conf_id in range(mol.GetNumConformers()):
                if self.force_field.startswith("MMFF94"):
                    # Try MMFF94s first
                    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
                    if mp is None:
                        # Fall back to MMFF94
                        mp = AllChem.MMFFGetMoleculeProperties(mol)
                    if mp is not None:
                        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=conf_id)
                        if ff is not None:
                            energy = ff.CalcEnergy()
                            energies.append(energy)
                            continue
                # Fall back to UFF
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                if ff is not None:
                    energy = ff.CalcEnergy()
                    energies.append(energy)

            return energies

        except Exception as e:
            self.logger.error(f"Error calculating energies: {str(e)}")
            return []

    def get_conformer_stats(self, mol: Chem.Mol) -> Optional[ConformerStats]:
        """Calculate statistics for conformer ensemble.

        Args:
            mol: Molecule with conformers

        Returns:
            ConformerStats object or None if calculation fails
        """
        if not mol or mol.GetNumConformers() <= 1:
            return None

        try:
            # Get energies
            energies = self.get_conformer_energies(mol)
            if not energies:
                return None

            # Calculate energy statistics
            energy_array = np.array(energies)
            energy_stats = {
                "min": float(np.min(energy_array)),
                "max": float(np.max(energy_array)),
                "mean": float(np.mean(energy_array)),
                "std": float(np.std(energy_array)),
            }

            # Calculate RMSD matrix
            n_conf = mol.GetNumConformers()
            rmsd_values = []
            for i in range(n_conf):
                for j in range(i + 1, n_conf):
                    rmsd = AllChem.GetConformerRMS(mol, i, j)
                    rmsd_values.append(rmsd)

            # Calculate RMSD statistics
            rmsd_array = np.array(rmsd_values)
            rmsd_stats = {
                "min": float(np.min(rmsd_array)),
                "max": float(np.max(rmsd_array)),
                "mean": float(np.mean(rmsd_array)),
                "std": float(np.std(rmsd_array)),
            }

            return ConformerStats(
                num_conformers=n_conf,
                energy_min=energy_stats["min"],
                energy_max=energy_stats["max"],
                energy_mean=energy_stats["mean"],
                energy_std=energy_stats["std"],
                rmsd_min=rmsd_stats["min"],
                rmsd_max=rmsd_stats["max"],
                rmsd_mean=rmsd_stats["mean"],
                rmsd_std=rmsd_stats["std"],
            )

        except Exception as e:
            self.logger.error(f"Error calculating conformer statistics: {str(e)}")
            return None

    def get_torsion_angles(self, mol: Chem.Mol, conf_id: int = 0) -> List[Tuple[Tuple[int, int, int, int], float]]:
        """Get rotatable bond torsion angles.

        Args:
            mol: Molecule with conformers
            conf_id: Conformer ID to analyze

        Returns:
            List of ((atom1, atom2, atom3, atom4), angle) tuples
        """
        if not mol or not mol.GetNumConformers():
            return []

        try:
            # Get rotatable bonds
            rot_bonds = TorsionFingerprints.FindRotatableBonds(mol)

            # Get torsion angles
            torsions = []
            for bond in rot_bonds:
                # Get atoms defining torsion
                a1, a2 = bond
                # Get neighboring atoms
                for n1 in mol.GetAtomWithIdx(a1).GetNeighbors():
                    n1_idx = n1.GetIdx()
                    if n1_idx != a2:
                        for n2 in mol.GetAtomWithIdx(a2).GetNeighbors():
                            n2_idx = n2.GetIdx()
                            if n2_idx != a1:
                                # Get torsion angle
                                angle = GetDihedralDeg(mol.GetConformer(conf_id), n1_idx, a1, a2, n2_idx)
                                torsions.append(((n1_idx, a1, a2, n2_idx), angle))
                                break
                        break

            return torsions

        except Exception as e:
            self.logger.error(f"Error calculating torsion angles: {str(e)}")
            return []

    def set_torsion_angles(self, mol: Chem.Mol, torsion_angles: List[Tuple[Tuple[int, int, int, int], float]], conf_id: int = 0) -> bool:
        """Set torsion angles for rotatable bonds.

        Args:
            mol: Molecule with conformers
            torsion_angles: List of ((atom1, atom2, atom3, atom4), angle) tuples
            conf_id: Conformer ID to modify

        Returns:
            Whether setting angles succeeded
        """
        if not mol or not mol.GetNumConformers():
            return False

        try:
            conf = mol.GetConformer(conf_id)
            for atoms, angle in torsion_angles:
                SetDihedralDeg(conf, *atoms, angle)
            return True

        except Exception as e:
            self.logger.error(f"Error setting torsion angles: {str(e)}")
            return False
