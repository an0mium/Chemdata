"""Structure standardization and validation using RDKit."""

from typing import Tuple, Optional

from rdkit import Chem
from rdkit.Chem import rdDepictor

from ....logger import LogManager


class StructureStandardizer:
    """Handle chemical structure standardization and validation."""

    def __init__(self):
        """Initialize structure standardizer."""
        self.logger = LogManager().get_logger("structure_standardizer")

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
            rdDepictor.Compute2DCoords(mol)
            rdDepictor.GenerateDepictionMatching2DStructure(mol)

            # Return canonical SMILES
            return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True, kekuleSmiles=True)

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

    def generate_inchi(self, mol: Chem.Mol) -> Optional[Tuple[str, str]]:
        """
        Generate InChI and InChIKey.

        Args:
            mol: RDKit molecule

        Returns:
            Tuple of (InChI, InChIKey) or None if failed
        """
        try:
            inchi = Chem.MolToInchi(mol)
            if not inchi:
                return None

            inchi_key = Chem.MolToInchiKey(mol)
            if not inchi_key:
                return None

            return inchi, inchi_key

        except Exception as e:
            self.logger.error(f"Error generating InChI: {str(e)}")
            return None

    def standardize_stereo(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """
        Standardize stereochemistry.

        Args:
            mol: RDKit molecule

        Returns:
            Standardized molecule or None if failed
        """
        try:
            # Find potential stereocenters
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

            # Count unspecified stereocenters
            unspec = Chem.CalcNumUnspecifiedAtomStereoCenters(mol)
            if unspec > 0:
                self.logger.warning(f"Molecule has {unspec} unspecified stereocenters")

            return mol

        except Exception as e:
            self.logger.error(f"Error standardizing stereochemistry: {str(e)}")
            return None

    def standardize_tautomer(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """
        Standardize tautomeric form.

        Args:
            mol: RDKit molecule

        Returns:
            Standardized molecule or None if failed
        """
        try:
            # Enumerate tautomers and pick canonical one
            from rdkit.Chem.MolStandardize import rdMolStandardize

            enumerator = rdMolStandardize.TautomerEnumerator()
            canon_taut = enumerator.Canonicalize(mol)

            return canon_taut

        except Exception as e:
            self.logger.error(f"Error standardizing tautomer: {str(e)}")
            return None

    def standardize_charges(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """
        Standardize formal charges.

        Args:
            mol: RDKit molecule

        Returns:
            Standardized molecule or None if failed
        """
        try:
            # Reionize to ensure standard charge states
            from rdkit.Chem.MolStandardize import rdMolStandardize

            reionizer = rdMolStandardize.Reionizer()
            reionized_mol = reionizer.reionize(mol)

            return reionized_mol

        except Exception as e:
            self.logger.error(f"Error standardizing charges: {str(e)}")
            return None

    def standardize_all(self, smiles: str) -> Optional[str]:
        """
        Apply all standardization steps.

        Args:
            smiles: Input SMILES string

        Returns:
            Fully standardized SMILES string or None if failed
        """
        try:
            # Create molecule
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            # Add hydrogens
            mol = Chem.AddHs(mol)

            # Apply standardizations
            mol = self.standardize_charges(mol)
            if mol is None:
                return None

            mol = self.standardize_tautomer(mol)
            if mol is None:
                return None

            mol = self.standardize_stereo(mol)
            if mol is None:
                return None

            # Remove hydrogens
            mol = Chem.RemoveHs(mol)

            # Generate 2D coordinates
            rdDepictor.Compute2DCoords(mol)
            rdDepictor.GenerateDepictionMatching2DStructure(mol)

            # Return canonical SMILES
            return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True, kekuleSmiles=True)

        except Exception as e:
            self.logger.error(f"Error in complete standardization: {str(e)}")
            return None
