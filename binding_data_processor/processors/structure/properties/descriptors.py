"""Molecular descriptor calculations using RDKit."""

from typing import Dict, Any

from rdkit import Chem
from rdkit.Chem import (
    Descriptors, rdMolDescriptors, Crippen, GraphDescriptors, QED
)

from ....logger import LogManager


class DescriptorCalculator:
    """Calculate molecular descriptors using RDKit."""
    
    def __init__(self):
        """Initialize descriptor calculator."""
        self.logger = LogManager().get_logger("descriptor_calculator")

    def calculate_all_descriptors(self, mol: Chem.Mol) -> Dict[str, Any]:
        """
        Calculate all available descriptors.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Dictionary of descriptor values
        """
        descriptors = {}
        descriptors.update(self.calculate_constitutional_descriptors(mol))
        descriptors.update(self.calculate_topological_descriptors(mol))
        descriptors.update(self.calculate_physicochemical_descriptors(mol))
        descriptors.update(self.calculate_druglikeness_descriptors(mol))
        return descriptors

    def calculate_constitutional_descriptors(self, mol: Chem.Mol) -> Dict[str, Any]:
        """
        Calculate constitutional descriptors.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Dictionary of descriptor values
        """
        try:
            return {
                'num_atoms': mol.GetNumAtoms(),
                'num_heavy_atoms': mol.GetNumHeavyAtoms(),
                'num_bonds': mol.GetNumBonds(),
                'num_rotatable_bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
                'num_rings': rdMolDescriptors.CalcNumRings(mol),
                'num_aromatic_rings': rdMolDescriptors.CalcNumAromaticRings(mol),
                'num_aliphatic_rings': rdMolDescriptors.CalcNumAliphaticRings(mol),
                'num_saturated_rings': rdMolDescriptors.CalcNumSaturatedRings(mol),
                'num_heterocycles': rdMolDescriptors.CalcNumHeterocycles(mol),
                'num_stereocenters': rdMolDescriptors.CalcNumAtomStereoCenters(mol),
                'fraction_sp3': rdMolDescriptors.CalcFractionCSP3(mol),
                'num_hbd': rdMolDescriptors.CalcNumHBD(mol),
                'num_hba': rdMolDescriptors.CalcNumHBA(mol)
            }
        except Exception as e:
            self.logger.error(f"Error calculating constitutional descriptors: {str(e)}")
            return {}

    def calculate_topological_descriptors(self, mol: Chem.Mol) -> Dict[str, Any]:
        """
        Calculate topological descriptors.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Dictionary of descriptor values
        """
        try:
            return {
                'bertz_ct': GraphDescriptors.BertzCT(mol),
                'chi0v': GraphDescriptors.Chi0v(mol),
                'chi1v': GraphDescriptors.Chi1v(mol),
                'chi2v': GraphDescriptors.Chi2v(mol),
                'chi3v': GraphDescriptors.Chi3v(mol),
                'chi4v': GraphDescriptors.Chi4v(mol),
                'hall_kier_alpha': GraphDescriptors.HallKierAlpha(mol),
                'kappa1': GraphDescriptors.Kappa1(mol),
                'kappa2': GraphDescriptors.Kappa2(mol),
                'kappa3': GraphDescriptors.Kappa3(mol)
            }
        except Exception as e:
            self.logger.error(f"Error calculating topological descriptors: {str(e)}")
            return {}

    def calculate_physicochemical_descriptors(self, mol: Chem.Mol) -> Dict[str, Any]:
        """
        Calculate physicochemical descriptors.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Dictionary of descriptor values
        """
        try:
            return {
                'molecular_weight': Descriptors.ExactMolWt(mol),
                'heavy_atom_weight': Descriptors.HeavyAtomMolWt(mol),
                'logp': Crippen.MolLogP(mol),
                'mr': Crippen.MolMR(mol),
                'tpsa': Descriptors.TPSA(mol),
                'formal_charge': Chem.GetFormalCharge(mol),
                'num_valence_electrons': Descriptors.NumValenceElectrons(mol),
                'polar_surface_area': Descriptors.TPSA(mol),
                'van_der_waals_volume': Descriptors.ComputeMolVolume(mol)
            }
        except Exception as e:
            self.logger.error(f"Error calculating physicochemical descriptors: {str(e)}")
            return {}

    def calculate_druglikeness_descriptors(self, mol: Chem.Mol) -> Dict[str, Any]:
        """
        Calculate drug-likeness descriptors.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Dictionary of descriptor values
        """
        try:
            descriptors = {
                'qed': QED.default(mol),
                'molecular_weight_violations': self._check_mw_violations(mol),
                'logp_violations': self._check_logp_violations(mol),
                'hbd_violations': self._check_hbd_violations(mol),
                'hba_violations': self._check_hba_violations(mol),
                'rotatable_bonds_violations': self._check_rotatable_violations(mol),
                'tpsa_violations': self._check_tpsa_violations(mol)
            }
            
            # Add Lipinski descriptors
            lipinski = self._calculate_lipinski(mol)
            descriptors.update(lipinski)
            
            # Add Veber descriptors
            veber = self._calculate_veber(mol)
            descriptors.update(veber)
            
            return descriptors
            
        except Exception as e:
            self.logger.error(f"Error calculating druglikeness descriptors: {str(e)}")
            return {}

    def _calculate_lipinski(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Calculate Lipinski's Rule of Five descriptors."""
        try:
            mw = Descriptors.ExactMolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            
            violations = 0
            if mw > 500:
                violations += 1
            if logp > 5:
                violations += 1
            if hbd > 5:
                violations += 1
            if hba > 10:
                violations += 1
            
            return {
                'lipinski_violations': violations,
                'lipinski_pass': violations <= 1
            }
            
        except Exception as e:
            self.logger.error(f"Error calculating Lipinski descriptors: {str(e)}")
            return {}

    def _calculate_veber(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Calculate Veber's druglikeness descriptors."""
        try:
            rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
            tpsa = Descriptors.TPSA(mol)
            
            violations = 0
            if rotatable_bonds > 10:
                violations += 1
            if tpsa > 140:
                violations += 1
            
            return {
                'veber_violations': violations,
                'veber_pass': violations == 0
            }
            
        except Exception as e:
            self.logger.error(f"Error calculating Veber descriptors: {str(e)}")
            return {}

    def _check_mw_violations(self, mol: Chem.Mol) -> bool:
        """Check molecular weight violations."""
        try:
            return Descriptors.ExactMolWt(mol) > 500
        except Exception as e:
            self.logger.error(f"Error checking MW violations: {str(e)}")
            return False

    def _check_logp_violations(self, mol: Chem.Mol) -> bool:
        """Check logP violations."""
        try:
            return Crippen.MolLogP(mol) > 5
        except Exception as e:
            self.logger.error(f"Error checking logP violations: {str(e)}")
            return False

    def _check_hbd_violations(self, mol: Chem.Mol) -> bool:
        """Check H-bond donor violations."""
        try:
            return rdMolDescriptors.CalcNumHBD(mol) > 5
        except Exception as e:
            self.logger.error(f"Error checking HBD violations: {str(e)}")
            return False

    def _check_hba_violations(self, mol: Chem.Mol) -> bool:
        """Check H-bond acceptor violations."""
        try:
            return rdMolDescriptors.CalcNumHBA(mol) > 10
        except Exception as e:
            self.logger.error(f"Error checking HBA violations: {str(e)}")
            return False

    def _check_rotatable_violations(self, mol: Chem.Mol) -> bool:
        """Check rotatable bonds violations."""
        try:
            return rdMolDescriptors.CalcNumRotatableBonds(mol) > 10
        except Exception as e:
            self.logger.error(f"Error checking rotatable bond violations: {str(e)}")
            return False

    def _check_tpsa_violations(self, mol: Chem.Mol) -> bool:
        """Check TPSA violations."""
        try:
            return Descriptors.TPSA(mol) > 140
        except Exception as e:
            self.logger.error(f"Error checking TPSA violations: {str(e)}")
            return False
