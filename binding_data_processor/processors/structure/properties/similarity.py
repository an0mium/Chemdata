"""Structure similarity calculations using RDKit."""

from typing import Dict, Any, List, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import (
    AllChem, DataStructs, rdFMCS, rdMolAlign, rdMolDescriptors,
    rdShapeHelpers
)

from ....logger import LogManager


class SimilarityCalculator:
    """Calculate molecular similarity using various methods."""
    
    def __init__(self):
        """Initialize similarity calculator."""
        self.logger = LogManager().get_logger("similarity_calculator")

    def calculate_all_similarities(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol
    ) -> Dict[str, Any]:
        """
        Calculate all similarity metrics between two molecules.
        
        Args:
            mol1: First RDKit molecule
            mol2: Second RDKit molecule
            
        Returns:
            Dictionary of similarity metrics
        """
        similarities = {}
        
        # 2D similarities
        similarities.update(self.calculate_2d_similarities(mol1, mol2))
        
        # 3D similarities (if conformers available)
        if (mol1.GetNumConformers() > 0 and mol2.GetNumConformers() > 0):
            similarities.update(self.calculate_3d_similarities(mol1, mol2))
            
        # Substructure similarity
        similarities.update(self.calculate_substructure_similarity(mol1, mol2))
        
        return similarities

    def calculate_2d_similarities(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol
    ) -> Dict[str, float]:
        """
        Calculate 2D fingerprint-based similarities.
        
        Args:
            mol1: First RDKit molecule
            mol2: Second RDKit molecule
            
        Returns:
            Dictionary of similarity scores
        """
        try:
            # Morgan fingerprints (ECFP4)
            fp1_morgan = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
            fp2_morgan = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
            morgan_sim = DataStructs.TanimotoSimilarity(fp1_morgan, fp2_morgan)
            
            # MACCS keys
            fp1_maccs = rdMolDescriptors.GetMACCSKeysFingerprint(mol1)
            fp2_maccs = rdMolDescriptors.GetMACCSKeysFingerprint(mol2)
            maccs_sim = DataStructs.TanimotoSimilarity(fp1_maccs, fp2_maccs)
            
            # Topological fingerprints
            fp1_topo = Chem.RDKFingerprint(mol1)
            fp2_topo = Chem.RDKFingerprint(mol2)
            topo_sim = DataStructs.TanimotoSimilarity(fp1_topo, fp2_topo)
            
            # Atom pairs
            fp1_pairs = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol1)
            fp2_pairs = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol2)
            pairs_sim = DataStructs.TanimotoSimilarity(fp1_pairs, fp2_pairs)
            
            # Torsions
            fp1_torsion = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol1)
            fp2_torsion = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol2)
            torsion_sim = DataStructs.TanimotoSimilarity(fp1_torsion, fp2_torsion)
            
            return {
                'morgan_similarity': morgan_sim,
                'maccs_similarity': maccs_sim,
                'topological_similarity': topo_sim,
                'atom_pairs_similarity': pairs_sim,
                'torsion_similarity': torsion_sim,
                'mean_2d_similarity': np.mean([
                    morgan_sim, maccs_sim, topo_sim, pairs_sim, torsion_sim
                ])
            }
            
        except Exception as e:
            self.logger.error(f"Error calculating 2D similarities: {str(e)}")
            return {}

    def calculate_3d_similarities(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol
    ) -> Dict[str, float]:
        """
        Calculate 3D shape-based similarities.
        
        Args:
            mol1: First RDKit molecule with 3D conformer
            mol2: Second RDKit molecule with 3D conformer
            
        Returns:
            Dictionary of similarity scores
        """
        try:
            # Align molecules
            rmsd = rdMolAlign.AlignMol(mol1, mol2)
            
            # Shape similarity
            protrude = rdShapeHelpers.ShapeProtrudeDist(
                mol1, mol2,
                allowReordering=True
            )
            volume = rdShapeHelpers.ShapeTanimotoDist(
                mol1, mol2,
                allowReordering=True
            )
            
            return {
                'rmsd': rmsd,
                'shape_protrude': protrude,
                'shape_volume': volume,
                'mean_3d_similarity': np.mean([
                    1.0 / (1.0 + rmsd),
                    1.0 - protrude,
                    1.0 - volume
                ])
            }
            
        except Exception as e:
            self.logger.error(f"Error calculating 3D similarities: {str(e)}")
            return {}

    def calculate_substructure_similarity(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol
    ) -> Dict[str, Any]:
        """
        Calculate similarity based on maximum common substructure.
        
        Args:
            mol1: First RDKit molecule
            mol2: Second RDKit molecule
            
        Returns:
            Dictionary with MCS information
        """
        try:
            # Find maximum common substructure
            mcs = rdFMCS.FindMCS(
                [mol1, mol2],
                bondCompare=rdFMCS.BondCompare.CompareOrder,
                ringMatchesRingOnly=True,
                completeRingsOnly=True
            )
            
            if mcs.numAtoms == 0:
                return {}
                
            # Create MCS molecule
            mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
            if mcs_mol is None:
                return {}
                
            # Get matches
            matches1 = mol1.GetSubstructMatches(mcs_mol)
            matches2 = mol2.GetSubstructMatches(mcs_mol)
            
            # Calculate similarity metrics
            mcs_size = mcs.numAtoms
            similarity = mcs_size / min(mol1.GetNumAtoms(), mol2.GetNumAtoms())
            
            return {
                'mcs_size': mcs_size,
                'mcs_similarity': similarity,
                'mcs_smarts': mcs.smartsString,
                'num_matches': len(matches1) * len(matches2)
            }
            
        except Exception as e:
            self.logger.error(f"Error calculating substructure similarity: {str(e)}")
            return {}

    def find_similar_compounds(
        self,
        query_mol: Chem.Mol,
        library_mols: List[Chem.Mol],
        similarity_threshold: float = 0.7
    ) -> List[Tuple[int, float]]:
        """
        Find similar compounds in a library.
        
        Args:
            query_mol: Query molecule
            library_mols: List of molecules to search
            similarity_threshold: Minimum similarity score (0-1)
            
        Returns:
            List of (compound_index, similarity_score) tuples
        """
        try:
            # Generate query fingerprints
            query_morgan = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2)
            query_maccs = rdMolDescriptors.GetMACCSKeysFingerprint(query_mol)
            
            similar_compounds = []
            for i, mol in enumerate(library_mols):
                try:
                    # Calculate fingerprint similarities
                    mol_morgan = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
                    mol_maccs = rdMolDescriptors.GetMACCSKeysFingerprint(mol)
                    
                    morgan_sim = DataStructs.TanimotoSimilarity(
                        query_morgan, mol_morgan
                    )
                    maccs_sim = DataStructs.TanimotoSimilarity(
                        query_maccs, mol_maccs
                    )
                    
                    # Use mean similarity
                    mean_sim = (morgan_sim + maccs_sim) / 2
                    
                    if mean_sim >= similarity_threshold:
                        similar_compounds.append((i, mean_sim))
                        
                except Exception as e:
                    self.logger.warning(
                        f"Error processing compound {i}: {str(e)}"
                    )
                    continue
                    
            # Sort by similarity score
            similar_compounds.sort(key=lambda x: x[1], reverse=True)
            return similar_compounds
            
        except Exception as e:
            self.logger.error(f"Error finding similar compounds: {str(e)}")
            return []

    def find_substructure_matches(
        self,
        query_mol: Chem.Mol,
        library_mols: List[Chem.Mol]
    ) -> List[int]:
        """
        Find compounds containing the query as a substructure.
        
        Args:
            query_mol: Query substructure
            library_mols: List of molecules to search
            
        Returns:
            List of matching compound indices
        """
        try:
            matches = []
            for i, mol in enumerate(library_mols):
                try:
                    if mol.HasSubstructMatch(query_mol):
                        matches.append(i)
                except Exception as e:
                    self.logger.warning(
                        f"Error processing compound {i}: {str(e)}"
                    )
                    continue
            return matches
            
        except Exception as e:
            self.logger.error(f"Error finding substructure matches: {str(e)}")
            return []
