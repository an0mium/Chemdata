"""Base functionality for chemical structure property calculations."""

from dataclasses import dataclass
from typing import Dict, Any, Optional, List, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdDepictor

from ....logger import LogManager


@dataclass
class PropertyCalculationConfig:
    """Configuration for property calculations."""
    
    generate_3d: bool = True           # Generate 3D conformers
    optimize_3d: bool = True           # Optimize 3D geometry
    num_conformers: int = 10           # Number of conformers to generate
    random_seed: int = 42              # Random seed for reproducibility
    energy_threshold: float = 10.0     # Energy threshold for conformer pruning (kcal/mol)
    rmsd_threshold: float = 0.5        # RMSD threshold for conformer clustering (Å)


class PropertyCalculator:
    """Base class for chemical property calculations."""
    
    def __init__(self, config: Optional[PropertyCalculationConfig] = None):
        """
        Initialize property calculator.
        
        Args:
            config: Optional configuration, uses defaults if not provided
        """
        self.config = config or PropertyCalculationConfig()
        self.logger = LogManager().get_logger("property_calculator")

    def prepare_molecule(self, smiles: str) -> Optional[Chem.Mol]:
        """
        Prepare molecule for property calculations.
        
        Args:
            smiles: SMILES string of molecule
            
        Returns:
            Prepared RDKit molecule or None if preparation fails
        """
        try:
            # Create molecule from SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                self.logger.warning(f"Failed to create molecule from SMILES: {smiles}")
                return None
                
            # Generate 2D coordinates if needed
            if not mol.GetNumConformers():
                rdDepictor.Compute2DCoords(mol)
                
            # Generate 3D conformer if requested
            if self.config.generate_3d:
                mol = self._generate_3d_conformer(mol)
                if mol is None:
                    return None
                    
            return mol
            
        except Exception as e:
            self.logger.error(f"Error preparing molecule: {str(e)}")
            return None

    def _generate_3d_conformer(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """
        Generate 3D conformer for molecule.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Molecule with 3D conformer or None if generation fails
        """
        try:
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Generate conformers
            params = AllChem.ETKDGv3()
            params.randomSeed = self.config.random_seed
            params.numThreads = 0  # Use all available threads
            
            # Embed multiple conformers
            n_confs = AllChem.EmbedMultipleConfs(
                mol,
                numConfs=self.config.num_conformers,
                params=params,
                pruneRmsThresh=self.config.rmsd_threshold
            )
            
            if n_confs == -1:
                self.logger.warning("Failed to generate any conformers")
                return None
                
            # Optimize conformers if requested
            if self.config.optimize_3d:
                self._optimize_conformers(mol)
                
            # Select lowest energy conformer
            energies = self._calculate_conformer_energies(mol)
            if not energies:
                return None
                
            # Keep only lowest energy conformer
            lowest_conf = np.argmin(energies)
            conf_ids = list(range(mol.GetNumConformers()))
            conf_ids.remove(lowest_conf)
            for conf_id in conf_ids:
                mol.RemoveConformer(conf_id)
                
            return mol
            
        except Exception as e:
            self.logger.error(f"Error generating 3D conformer: {str(e)}")
            return None

    def _optimize_conformers(self, mol: Chem.Mol) -> None:
        """
        Optimize conformer geometries using MMFF94s force field.
        
        Args:
            mol: RDKit molecule with conformers
        """
        for conf_id in range(mol.GetNumConformers()):
            try:
                # Try MMFF94s first
                if AllChem.MMFFHasAllMoleculeParams(mol):
                    AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
                # Fall back to UFF if MMFF94s fails
                else:
                    AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
            except Exception as e:
                self.logger.warning(
                    f"Failed to optimize conformer {conf_id}: {str(e)}"
                )

    def _calculate_conformer_energies(self, mol: Chem.Mol) -> List[float]:
        """
        Calculate energies for all conformers.
        
        Args:
            mol: RDKit molecule with conformers
            
        Returns:
            List of conformer energies
        """
        energies = []
        for conf_id in range(mol.GetNumConformers()):
            try:
                # Try MMFF94s first
                if AllChem.MMFFHasAllMoleculeParams(mol):
                    ff = AllChem.MMFFGetMoleculeForceField(
                        mol,
                        AllChem.MMFFGetMoleculeProperties(mol),
                        confId=conf_id
                    )
                # Fall back to UFF if MMFF94s fails
                else:
                    ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                    
                if ff:
                    energy = ff.CalcEnergy()
                    energies.append(energy)
                    
            except Exception as e:
                self.logger.warning(
                    f"Failed to calculate energy for conformer {conf_id}: {str(e)}"
                )
                
        return energies

    def calculate_basic_properties(self, mol: Chem.Mol) -> Dict[str, Any]:
        """
        Calculate basic molecular properties.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Dictionary of calculated properties
        """
        try:
            props = {
                'molecular_weight': Descriptors.ExactMolWt(mol),
                'heavy_atoms': mol.GetNumHeavyAtoms(),
                'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
                'rings': Descriptors.RingCount(mol),
                'aromatic_rings': Descriptors.NumAromaticRings(mol),
                'hbd': Descriptors.NumHDonors(mol),
                'hba': Descriptors.NumHAcceptors(mol),
                'tpsa': Descriptors.TPSA(mol),
                'logp': Descriptors.MolLogP(mol),
                'charge': Chem.GetFormalCharge(mol)
            }
            
            # Add 3D properties if available
            if mol.GetNumConformers() > 0:
                props.update(self._calculate_3d_properties(mol))
                
            return props
            
        except Exception as e:
            self.logger.error(f"Error calculating basic properties: {str(e)}")
            return {}

    def _calculate_3d_properties(self, mol: Chem.Mol) -> Dict[str, Any]:
        """
        Calculate 3D molecular properties.
        
        Args:
            mol: RDKit molecule with 3D conformer
            
        Returns:
            Dictionary of 3D properties
        """
        props = {}
        
        # Calculate surface area and volume
        self._add_surface_properties(mol, props)
        
        # Calculate shape properties
        self._add_shape_properties(mol, props)
        
        return props

    def _add_surface_properties(self, mol: Chem.Mol, props: Dict[str, Any]) -> None:
        """Add surface area and volume to properties."""
        try:
            if hasattr(AllChem, 'ComputeMolSurf'):
                props['surface_area'] = AllChem.ComputeMolSurf(mol)
            if hasattr(Descriptors, 'ComputeMolVolume'):
                props['volume'] = Descriptors.ComputeMolVolume(mol)
        except Exception as e:
            self.logger.error(f"Error calculating surface properties: {str(e)}")

    def _add_shape_properties(self, mol: Chem.Mol, props: Dict[str, Any]) -> None:
        """Add shape-related properties."""
        try:
            # Principal moments
            moments = self._calculate_principal_moments(mol)
            if moments:
                props['principal_moments'] = moments
                
            # Radius of gyration
            rg = self._calculate_radius_of_gyration(mol)
            if rg:
                props['radius_of_gyration'] = rg
                
        except Exception as e:
            self.logger.error(f"Error calculating shape properties: {str(e)}")

    def _calculate_principal_moments(
        self,
        mol: Chem.Mol
    ) -> Optional[Tuple[float, float, float]]:
        """
        Calculate principal moments of inertia.
        
        Args:
            mol: RDKit molecule with 3D conformer
            
        Returns:
            Tuple of principal moments or None if calculation fails
        """
        try:
            conf = mol.GetConformer()
            
            # Get atomic masses and coordinates
            masses = []
            coords = []
            for i, atom in enumerate(mol.GetAtoms()):
                masses.append(atom.GetMass())
                pos = conf.GetAtomPosition(i)
                coords.append([pos.x, pos.y, pos.z])
                
            # Convert to numpy arrays
            masses = np.array(masses)
            coords = np.array(coords)
            
            # Calculate center of mass
            com = np.average(coords, weights=masses, axis=0)
            
            # Translate to center of mass
            coords -= com
            
            # Calculate inertia tensor
            inertia_tensor = np.zeros((3, 3))
            for m, r in zip(masses, coords):
                x, y, z = r
                inertia_tensor[0, 0] += m * (y * y + z * z)
                inertia_tensor[1, 1] += m * (x * x + z * z)
                inertia_tensor[2, 2] += m * (x * x + y * y)
                inertia_tensor[0, 1] -= m * x * y
                inertia_tensor[0, 2] -= m * x * z
                inertia_tensor[1, 2] -= m * y * z
            inertia_tensor[1, 0] = inertia_tensor[0, 1]
            inertia_tensor[2, 0] = inertia_tensor[0, 2]
            inertia_tensor[2, 1] = inertia_tensor[1, 2]
            
            # Get eigenvalues (principal moments)
            moments = np.sort(np.linalg.eigvals(inertia_tensor))
            return tuple(moments)
            
        except Exception as e:
            self.logger.error(f"Error calculating principal moments: {str(e)}")
            return None

    def _calculate_radius_of_gyration(self, mol: Chem.Mol) -> Optional[float]:
        """
        Calculate radius of gyration.
        
        Args:
            mol: RDKit molecule with 3D conformer
            
        Returns:
            Radius of gyration in Angstroms or None if calculation fails
        """
        try:
            conf = mol.GetConformer()
            
            # Get atomic masses and coordinates
            masses = []
            coords = []
            for i, atom in enumerate(mol.GetAtoms()):
                masses.append(atom.GetMass())
                pos = conf.GetAtomPosition(i)
                coords.append([pos.x, pos.y, pos.z])
                
            # Convert to numpy arrays
            masses = np.array(masses)
            coords = np.array(coords)
            
            # Calculate center of mass
            com = np.average(coords, weights=masses, axis=0)
            
            # Calculate radius of gyration
            rg2 = np.sum(masses * np.sum((coords - com)**2, axis=1)) / np.sum(masses)
            return np.sqrt(rg2)
            
        except Exception as e:
            self.logger.error(f"Error calculating radius of gyration: {str(e)}")
            return None
