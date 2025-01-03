"""Chemical property calculations using RDKit."""

from typing import Dict, Any, Optional, Tuple, List
import numpy as np

from rdkit import Chem
from rdkit.Chem import (
    Descriptors, AllChem, Crippen, rdMolDescriptors, 
    rdDepictor, rdMolTransforms, rdFMCS, rdMolAlign
)
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem.rdDepictor import Compute2DCoords

from logger import LogManager


class ChemicalProperties:
    """Handles chemical property calculations using RDKit."""
    
    def __init__(self):
        """Initialize chemical properties calculator."""
        self.logger = LogManager().get_logger("chemical_properties")

    def calculate_properties(self, smiles: str, generate_3d: bool = True) -> Dict[str, Any]:
        """
        Calculate molecular properties from SMILES using RDKit.
        
        Args:
            smiles: SMILES string of the compound
            generate_3d: Whether to generate and optimize 3D conformation
            
        Returns:
            Dictionary of calculated properties
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Invalid SMILES string: {smiles}")
                
            # Generate 3D conformation if requested
            if generate_3d:
                mol = Chem.AddHs(mol)  # Add hydrogens
                AllChem.EmbedMolecule(mol, randomSeed=42)  # Generate 3D coords
                AllChem.MMFFOptimizeMolecule(mol)  # Optimize geometry
                
            # Basic properties
            props = {
                'molecular_weight': Descriptors.ExactMolWt(mol),
                'logp': Crippen.MolLogP(mol),
                'hbd': rdMolDescriptors.CalcNumHBD(mol),
                'hba': rdMolDescriptors.CalcNumHBA(mol),
                'tpsa': Descriptors.TPSA(mol),
                'rotatable_bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
                'rings': rdMolDescriptors.CalcNumRings(mol),
                'aromatic_rings': rdMolDescriptors.CalcNumAromaticRings(mol),
                'heavy_atoms': mol.GetNumHeavyAtoms(),
                'fraction_sp3': rdMolDescriptors.CalcFractionCSP3(mol),
                'complexity': Descriptors.BertzCT(mol)
            }
            
            # 3D properties if available
            if generate_3d:
                props.update({
                    'surface_area': AllChem.ComputeMolSurf(mol),
                    'volume': AllChem.ComputeMolVolume(mol),
                    'principal_moments': self._calculate_principal_moments(mol),
                    'radius_of_gyration': self._calculate_radius_of_gyration(mol),
                    'shape_factors': self._calculate_shape_factors(mol),
                    'conformer_energies': self._calculate_conformer_energies(mol)
                })
            
            # Generate InChI and InChIKey
            inchi = Chem.MolToInchi(mol)
            if inchi:
                props['inchi'] = inchi
                props['inchi_key'] = Chem.MolToInchiKey(mol)
                
            # Additional descriptors
            props.update({
                'qed': Descriptors.qed(mol),  # Drug-likeness score
                'sas': Descriptors.sas(mol),  # Synthetic accessibility score
                'charge': Chem.GetFormalCharge(mol),
                'stereocenters': rdMolDescriptors.CalcNumAtomStereoCenters(mol),
                'unspecified_stereocenters': (
                    rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)
                ),
                'topological_polar_surface_area': Descriptors.TPSA(mol),
                'van_der_waals_volume': Descriptors.ComputeMolVolume(mol)
            })
            
            # Clean up
            if generate_3d:
                mol = Chem.RemoveHs(mol)
                
            return props
            
        except Exception as e:
            self.logger.error(f"Error calculating properties: {str(e)}")
            return {}

    def _calculate_principal_moments(self, mol: Chem.Mol) -> Tuple[float, float, float]:
        """
        Calculate principal moments of inertia.
        
        Args:
            mol: RDKit molecule with 3D coordinates
            
        Returns:
            Tuple of principal moments (I1, I2, I3)
        """
        try:
            # Get atomic masses and coordinates
            masses = []
            coords = []
            conf = mol.GetConformer()
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
            I = np.zeros((3, 3))
            for m, r in zip(masses, coords):
                x, y, z = r
                I[0, 0] += m * (y * y + z * z)
                I[1, 1] += m * (x * x + z * z)
                I[2, 2] += m * (x * x + y * y)
                I[0, 1] -= m * x * y
                I[0, 2] -= m * x * z
                I[1, 2] -= m * y * z
            I[1, 0] = I[0, 1]
            I[2, 0] = I[0, 2]
            I[2, 1] = I[1, 2]
            
            # Get eigenvalues (principal moments)
            moments = np.sort(np.linalg.eigvals(I))
            return tuple(moments)
            
        except Exception as e:
            self.logger.error(f"Error calculating principal moments: {str(e)}")
            return (0.0, 0.0, 0.0)

    def _calculate_radius_of_gyration(self, mol: Chem.Mol) -> float:
        """
        Calculate radius of gyration.
        
        Args:
            mol: RDKit molecule with 3D coordinates
            
        Returns:
            Radius of gyration in Angstroms
        """
        try:
            # Get atomic masses and coordinates
            masses = []
            coords = []
            conf = mol.GetConformer()
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
            return 0.0

    def _calculate_shape_factors(self, mol: Chem.Mol) -> Dict[str, float]:
        """
        Calculate molecular shape factors.
        
        Args:
            mol: RDKit molecule with 3D coordinates
            
        Returns:
            Dictionary of shape factors
        """
        try:
            moments = self._calculate_principal_moments(mol)
            if not any(moments):
                return {}
                
            # Sort moments in descending order
            I1, I2, I3 = sorted(moments, reverse=True)
            
            # Calculate shape factors
            return {
                'asphericity': I1 - 0.5 * (I2 + I3),
                'acylindricity': I2 - I3,
                'relative_shape_anisotropy': (
                    (I1 - I2)**2 + (I2 - I3)**2 + (I1 - I3)**2
                ) / (2 * (I1 + I2 + I3)**2)
            }
            
        except Exception as e:
            self.logger.error(f"Error calculating shape factors: {str(e)}")
            return {}

    def standardize_structure(self, smiles: str) -> Optional[str]:
        """
        Standardize chemical structure.
        
        Args:
            smiles: Input SMILES string
            
        Returns:
            Standardized SMILES string or None if failed
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
                
            # Remove hydrogens
            mol = Chem.RemoveHs(mol)
            
            # Kekulize
            Chem.Kekulize(mol)
            
            # Generate 2D coordinates with improved layout
            Compute2DCoords(mol)
            rdDepictor.GenerateDepictionMatching2DStructure(mol)
            
            # Return canonical SMILES
            return Chem.MolToSmiles(
                mol,
                isomericSmiles=True,
                canonical=True,
                kekuleSmiles=True
            )
            
        except Exception as e:
            self.logger.error(f"Error standardizing structure: {str(e)}")
            return None

    def validate_structure(self, smiles: str) -> Tuple[bool, str]:
        """
        Validate chemical structure.
        
        Args:
            smiles: SMILES string to validate
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False, "Invalid SMILES string"
                
            # Check for valence errors
            try:
                Chem.SanitizeMol(mol)
            except Exception as e:
                return False, f"Structure sanitization failed: {str(e)}"
                
            # Check for disconnected fragments
            if len(Chem.GetMolFrags(mol)) > 1:
                return False, "Structure contains disconnected fragments"
                
            # Check for unusual valences
            for atom in mol.GetAtoms():
                if atom.GetImplicitValence() == -1:
                    return False, f"Unusual valence for atom {atom.GetSymbol()}"
                    
            return True, "Valid structure"
            
        except Exception as e:
            return False, f"Validation error: {str(e)}"

    def generate_conformers(self, mol: Chem.Mol, n_confs: int = 10) -> List[float]:
        """
        Generate multiple conformers and return their energies.
        
        Args:
            mol: RDKit molecule
            n_confs: Number of conformers to generate
            
        Returns:
            List of conformer energies
        """
        try:
            # Add hydrogens and generate conformers
            mol = Chem.AddHs(mol)
            AllChem.EmbedMultipleConfs(
                mol,
                numConfs=n_confs,
                randomSeed=42,
                pruneRmsThresh=0.5  # Remove similar conformers
            )
            
            # Optimize all conformers
            energies = []
            for conf_id in range(mol.GetNumConformers()):
                # Optimize geometry
                AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
                
                # Calculate energy
                props = AllChem.MMFFGetMoleculeProperties(mol)
                energy = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
                if energy:
                    energies.append(energy.CalcEnergy())
            
            return energies
            
        except Exception as e:
            self.logger.error(f"Error generating conformers: {str(e)}")
            return []

    def _calculate_conformer_energies(self, mol: Chem.Mol) -> Dict[str, float]:
        """
        Calculate conformer energies and statistics.
        
        Args:
            mol: RDKit molecule with 3D coordinates
            
        Returns:
            Dictionary of energy statistics
        """
        try:
            energies = self.generate_conformers(mol)
            if not energies:
                return {}
                
            return {
                'min_energy': min(energies),
                'max_energy': max(energies),
                'mean_energy': np.mean(energies),
                'energy_range': max(energies) - min(energies),
                'energy_std': np.std(energies)
            }
            
        except Exception as e:
            self.logger.error(f"Error calculating conformer energies: {str(e)}")
            return {}

    def calculate_similarity(self, mol1: Chem.Mol, mol2: Chem.Mol) -> Dict[str, float]:
        """
        Calculate various similarity metrics between two molecules.
        
        Args:
            mol1: First RDKit molecule
            mol2: Second RDKit molecule
            
        Returns:
            Dictionary of similarity metrics
        """
        try:
            # Find maximum common substructure
            mcs = rdFMCS.FindMCS([mol1, mol2])
            mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
            
            # Get matches for MCS in both molecules
            matches1 = mol1.GetSubstructMatches(mcs_mol)
            matches2 = mol2.GetSubstructMatches(mcs_mol)
            
            # Calculate RMSD after alignment using MCS
            if matches1 and matches2:
                # Transform molecules to align MCS
                conf1 = mol1.GetConformer()
                conf2 = mol2.GetConformer()
                match1 = matches1[0]
                match2 = matches2[0]
                
                # Use rdMolTransforms to align molecules
                rdMolTransforms.AlignMolConformers(mol1, mol2, atomIds1=match1, atomIds2=match2)
                rmsd = rdMolTransforms.GetBestRMS(mol1, mol2, match1, match2)
            else:
                rmsd = float('inf')
            
            # Calculate various similarity metrics
            similarities = {
                'mcs_size': mcs.numAtoms,
                'mcs_fraction': mcs.numAtoms / min(mol1.GetNumAtoms(), mol2.GetNumAtoms()),
                'shape_similarity': rdMolAlign.AlignMol(mol1, mol2),
                'property_similarity': self._calculate_property_similarity(mol1, mol2),
                'rmsd_after_alignment': rmsd,
                'num_mcs_matches': len(matches1) * len(matches2)
            }
            
            return similarities
            
        except Exception as e:
            self.logger.error(f"Error calculating similarity: {str(e)}")
            return {}

    def _calculate_property_similarity(self, mol1: Chem.Mol, mol2: Chem.Mol) -> float:
        """
        Calculate similarity based on molecular properties.
        
        Args:
            mol1: First RDKit molecule
            mol2: Second RDKit molecule
            
        Returns:
            Property-based similarity score
        """
        try:
            # Calculate properties for both molecules
            props1 = {
                'mw': Descriptors.ExactMolWt(mol1),
                'logp': Crippen.MolLogP(mol1),
                'tpsa': Descriptors.TPSA(mol1),
                'hbd': rdMolDescriptors.CalcNumHBD(mol1),
                'hba': rdMolDescriptors.CalcNumHBA(mol1)
            }
            
            props2 = {
                'mw': Descriptors.ExactMolWt(mol2),
                'logp': Crippen.MolLogP(mol2),
                'tpsa': Descriptors.TPSA(mol2),
                'hbd': rdMolDescriptors.CalcNumHBD(mol2),
                'hba': rdMolDescriptors.CalcNumHBA(mol2)
            }
            
            # Calculate normalized differences
            diffs = []
            weights = {'mw': 0.2, 'logp': 0.2, 'tpsa': 0.2, 'hbd': 0.2, 'hba': 0.2}
            
            for prop in props1:
                max_val = max(abs(props1[prop]), abs(props2[prop]))
                if max_val > 0:
                    diff = abs(props1[prop] - props2[prop]) / max_val
                    diffs.append(diff * weights[prop])
            
            # Return similarity score (1 - average difference)
            return 1 - sum(diffs)
            
        except Exception as e:
            self.logger.error(f"Error calculating property similarity: {str(e)}")
            return 0.0
