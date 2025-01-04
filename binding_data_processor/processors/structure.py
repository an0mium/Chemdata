"""Chemical structure processing and analysis.

This module provides comprehensive functionality for:
1. Structure validation and standardization
2. 3D conformation generation and optimization 
3. Molecular descriptor calculation
4. Structural feature and pharmacophore detection
5. Chemical similarity calculations
6. Structure depiction generation
"""

import logging
import io
from typing import Dict, List, Optional, Set, Tuple, Union
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Draw,
    Fragments,
    rdDecomposition,
    rdMolDescriptors,
    rdMolStandardize,
)

logger = logging.getLogger(__name__)


class StructureProcessor:
    """Handles chemical structure processing and analysis."""

    # Combined SMARTS patterns for structural features
    STRUCTURE_PATTERNS = {
        # Basic functional groups
        "basic_amine": "[NX3;H2,H1;!$(NC=O)]",
        "acidic_oh": "[OH;$(O[#6]);!$(OC=O)]",
        "phenol": "[OH]c1[c,n][c,n][c,n][c,n][c,n]1",
        "carboxylic_acid": "[CX3](=O)[OX2H1]",
        "sulfonamide": "[SX4](=[OX1])(=[OX1])[NX3H2]",
        "amide": "[NX3][CX3](=[OX1])[#6]",
        "ester": "[#6][CX3](=O)[OX2H0][#6]",
        "phosphate": "[PX4](=[OX1])([OX2H1])([OX2H1])[OX2H1]",
        "sulfate": "[SX4](=[OX1])(=[OX1])([OX2H1])[OX2H1]",
        # Ring systems
        "aromatic_ring": "[a;r5,r6]1[a]:[a]:[a]:[a]:[a]1",
        "aliphatic_ring": "[A;R]1[A;R][A;R][A;R][A;R][A;R]1",
        "heterocycle": "[a;!c]",
        # Other features
        "halogen": "[F,Cl,Br,I]",
        "alkene": "[CX3]=[CX3]",
        "alkyne": "[CX2]#[CX2]",
        "ether": "[OX2]([#6])[#6]",
        "ketone": "[#6][CX3](=O)[#6]",
        "alcohol": "[OX2H][CX4]",
        "guanidine": "[N;H2]C(=N)N",
        "alkyl_chain": "[CR0;H2,H3][CR0;H2][CR0;H2,H3]",
        "michael_acceptor": "[CX3]=[CX3][$([CX3]=[OX1]),$([CX3](=O)[OX2])]",
        "acyl_halide": "[CX3](=[OX1])[F,Cl,Br,I]",
        "hydrazine": "[NX3][NX3]",
        "azide": "[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]",
        "thiol": "[SX2H]",
        "thioether": "[#16X2][#6]",
        "urea": "[NX3][CX3](=[OX1])[NX3]",
        "carbonate": "[OX2][CX3](=[OX1])[OX2]",
    }

    # Pharmacophore feature definitions
    PHARMACOPHORE_FEATURES = {
        "hydrophobic": [
            "[C;D3,D4;!$(C=*);!$(C[O,N,F,Cl,Br,I,S,P])]",
            "[C;r3,r4,r5,r6;!$(C=*);!$(C[O,N,F,Cl,Br,I,S,P])]",
        ],
        "aromatic": [
            "[a;r5,r6]1[a;r5,r6][a;r5,r6][a;r5,r6][a;r5,r6][a;r5,r6]1",
            "[a;r5]1[a;r5][a;r5][a;r5][a;r5]1",
        ],
        "hbond_donor": [
            "[N,O;H1,H2]-[C,c,S]",
            "[n;H1]",
        ],
        "hbond_acceptor": [
            "[N,O,F;!H1;v2]",
            "[n;H0]",
            "[O;H0;v2]",
        ],
        "positive": [
            "[N;H0;+1;D3]",
            "[N;H0;+1;D4]",
            "[n;H0;+1]",
        ],
        "negative": [
            "[C,S,P](=[O])[O-]",
            "[O;H1;-1]",
        ],
    }

    def __init__(self):
        """Initialize structure processor."""
        self.logger = logging.getLogger(__name__)

        # Pre-compile SMARTS patterns
        self.structure_patterns = {
            name: Chem.MolFromSmarts(smarts)
            for name, smarts in self.STRUCTURE_PATTERNS.items()
        }

        self.pharmacophore_patterns = {
            feature: [Chem.MolFromSmarts(smarts) for smarts in patterns]
            for feature, patterns in self.PHARMACOPHORE_FEATURES.items()
        }

        # Initialize standardizer
        self.standardizer = rdMolStandardize.CleanupParameters()
        self.standardizer.normalization = True
        self.standardizer.reionization = True
        self.standardizer.tautomerization = True

    def validate_structure(self, smiles: str) -> Tuple[bool, Optional[str]]:
        """
        Validate chemical structure.

        Args:
            smiles: SMILES string

        Returns:
            Tuple of (is_valid, error_message)
        """
        if not smiles:
            return False, "Empty SMILES string"

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False, "Invalid SMILES syntax"

            # Check for valid atoms
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    return False, f"Invalid atom at index {atom.GetIdx()}"

            # Check for valid bonds
            for bond in mol.GetBonds():
                if bond.GetBondType() == Chem.BondType.UNSPECIFIED:
                    return False, f"Invalid bond at index {bond.GetIdx()}"

            # Check for valid valences
            try:
                Chem.SanitizeMol(mol)
            except ValueError as e:
                return False, f"Structure sanitization failed: {str(e)}"

            return True, None

        except Exception as e:
            self.logger.error(f"Structure validation error: {str(e)}")
            return False, str(e)

    def standardize_structure(self, mol: Chem.Mol) -> Optional[str]:
        """
        Standardize chemical structure.

        Args:
            mol: RDKit molecule

        Returns:
            Standardized SMILES or None if invalid
        """
        try:
            # Remove hydrogens
            mol = Chem.RemoveHs(mol)

            # Standardize using RDKit's MolStandardize
            mol = rdMolStandardize.Cleanup(mol, self.standardizer)
            mol = rdMolStandardize.FragmentParent(mol)

            # Neutralize charges where possible
            for atom in mol.GetAtoms():
                charge = atom.GetFormalCharge()
                if charge != 0:
                    # Add/remove H to neutralize if possible
                    if charge > 0 and atom.GetNumImplicitHs() > 0:
                        atom.SetFormalCharge(charge - 1)
                        atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                    elif charge < 0:
                        atom.SetFormalCharge(charge + 1)
                        atom.SetNumExplicitHs(max(0, atom.GetNumExplicitHs() - 1))

            # Generate canonical SMILES
            return Chem.MolToSmiles(mol, canonical=True)

        except Exception as e:
            self.logger.error(f"Structure standardization error: {str(e)}")
            return None

    def generate_3d_conformation(
        self, mol: Chem.Mol, num_confs: int = 1
    ) -> Tuple[bool, Optional[List[float]]]:
        """
        Generate and optimize 3D conformations.

        Args:
            mol: RDKit molecule
            num_confs: Number of conformers to generate

        Returns:
            Tuple of (success, energies)
        """
        try:
            # Add hydrogens
            mol = Chem.AddHs(mol)

            # Generate conformers
            conf_ids = AllChem.EmbedMultipleConfs(
                mol,
                numConfs=num_confs,
                randomSeed=42,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True,
            )

            if not conf_ids:
                return False, None

            # Optimize conformers
            energies = []
            for conf_id in conf_ids:
                # Try MMFF94s first, fall back to UFF
                props = AllChem.MMFFGetMoleculeProperties(mol)
                if props is None:
                    success = AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
                    if success == 0:
                        energy = AllChem.UFFGetMoleculeForceField(
                            mol, confId=conf_id
                        ).CalcEnergy()
                else:
                    success = AllChem.MMFFOptimizeMolecule(
                        mol, confId=conf_id, mmffVariant="MMFF94s"
                    )
                    if success == 0:
                        energy = AllChem.MMFFGetMoleculeForceField(
                            mol, mmffVariant="MMFF94s", confId=conf_id
                        ).CalcEnergy()

                if success == 0:
                    energies.append(energy)

            return bool(energies), energies if energies else None

        except Exception as e:
            self.logger.error(f"3D conformation generation error: {str(e)}")
            return False, None

    def calculate_descriptors(self, mol: Chem.Mol) -> Dict[str, float]:
        """
        Calculate molecular descriptors.

        Args:
            mol: RDKit molecule

        Returns:
            Dictionary of descriptor values
        """
        try:
            descriptors = {
                # Physical properties
                "molecular_weight": Descriptors.ExactMolWt(mol),
                "heavy_atom_count": Descriptors.HeavyAtomCount(mol),
                "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
                "ring_count": Descriptors.RingCount(mol),
                "aromatic_rings": Descriptors.NumAromaticRings(mol),
                # Topological properties
                "tpsa": Descriptors.TPSA(mol),
                "logp": Descriptors.MolLogP(mol),
                "molar_refractivity": Descriptors.MolMR(mol),
                "balaban_j": Descriptors.BalabanJ(mol),
                "bertz_ct": Descriptors.BertzCT(mol),
                # Shape descriptors
                "asphericity": rdMolDescriptors.CalcAsphericity(mol),
                "eccentricity": rdMolDescriptors.CalcEccentricity(mol),
                "inertial_shape_factor": rdMolDescriptors.CalcInertialShapeFactor(mol),
                # Electronic properties
                "max_partial_charge": Descriptors.MaxPartialCharge(mol),
                "min_partial_charge": Descriptors.MinPartialCharge(mol),
                "max_abs_partial_charge": Descriptors.MaxAbsPartialCharge(mol),
                # Surface properties
                "sasa": Descriptors.SASA(mol),
                "polar_surface_area": Descriptors.TPSA(mol),
                # Fragment counts
                "hydrogen_bond_donors": Descriptors.NumHDonors(mol),
                "hydrogen_bond_acceptors": Descriptors.NumHAcceptors(mol),
                "sp3_carbons": Descriptors.NumSaturatedCarbocycles(mol),
                "sp2_carbons": Descriptors.NumAromaticCarbocycles(mol),
                "aliphatic_rings": Descriptors.NumAliphaticRings(mol),
                # Drug-likeness
                "qed": Descriptors.qed(mol),
            }

            # Add structural feature counts
            for name, pattern in self.structure_patterns.items():
                if pattern is not None:
                    descriptors[f"count_{name}"] = len(mol.GetSubstructMatches(pattern))

            # Add fragment counts
            descriptors.update(
                {
                    "alkyl_halide": Fragments.fr_alkyl_halide(mol),
                    "allylic_oxid": Fragments.fr_allylic_oxid(mol),
                    "aryl_methyl": Fragments.fr_aryl_methyl(mol),
                    "benzene": Fragments.fr_benzene(mol),
                    "bicyclic": Fragments.fr_bicyclic(mol),
                    "diazo": Fragments.fr_diazo(mol),
                    "halogen": Fragments.fr_halogen(mol),
                    "hdrzine": Fragments.fr_hdrzine(mol),
                    "hdrzone": Fragments.fr_hdrzone(mol),
                    "imide": Fragments.fr_imide(mol),
                    "isocyan": Fragments.fr_isocyan(mol),
                    "isothiocyan": Fragments.fr_isothiocyan(mol),
                    "ketone": Fragments.fr_ketone(mol),
                    "lactam": Fragments.fr_lactam(mol),
                    "lactone": Fragments.fr_lactone(mol),
                    "methoxy": Fragments.fr_methoxy(mol),
                    "morpholine": Fragments.fr_morpholine(mol),
                    "nitrile": Fragments.fr_nitrile(mol),
                    "nitro": Fragments.fr_nitro(mol),
                    "oxazole": Fragments.fr_oxazole(mol),
                    "oxime": Fragments.fr_oxime(mol),
                    "para_hydroxylation": Fragments.fr_para_hydroxylation(mol),
                    "phenol": Fragments.fr_phenol(mol),
                    "phenol_noOrthoHbond": Fragments.fr_phenol_noOrthoHbond(mol),
                    "phos_acid": Fragments.fr_phos_acid(mol),
                    "phos_ester": Fragments.fr_phos_ester(mol),
                    "piperdine": Fragments.fr_piperdine(mol),
                    "piperzine": Fragments.fr_piperzine(mol),
                    "priamide": Fragments.fr_priamide(mol),
                    "prisulfonamd": Fragments.fr_prisulfonamd(mol),
                    "pyridine": Fragments.fr_pyridine(mol),
                    "quatN": Fragments.fr_quatN(mol),
                    "sulfide": Fragments.fr_sulfide(mol),
                    "sulfonamd": Fragments.fr_sulfonamd(mol),
                    "sulfone": Fragments.fr_sulfone(mol),
                    "term_acetylene": Fragments.fr_term_acetylene(mol),
                    "tetrazole": Fragments.fr_tetrazole(mol),
                    "thiazole": Fragments.fr_thiazole(mol),
                    "thiocyan": Fragments.fr_thiocyan(mol),
                    "thiophene": Fragments.fr_thiophene(mol),
                    "unbrch_alkane": Fragments.fr_unbrch_alkane(mol),
                }
            )

            return {k: float(v) for k, v in descriptors.items()}

        except Exception as e:
            self.logger.error(f"Descriptor calculation error: {str(e)}")
            return {}

    def detect_pharmacophores(self, mol: Chem.Mol) -> Dict[str, List[Tuple[int, ...]]]:
        """
        Detect pharmacophore features.

        Args:
            mol: RDKit molecule

        Returns:
            Dictionary mapping feature types to atom indices
        """
        try:
            features = {}

            for feature_type, patterns in self.pharmacophore_patterns.items():
                matches = set()
                for pattern in patterns:
                    if pattern is not None:
                        # Get all matches and their atom indices
                        for match in mol.GetSubstructMatches(pattern):
                            matches.add(tuple(sorted(match)))
                if matches:
                    features[feature_type] = list(matches)

            return features

        except Exception as e:
            self.logger.error(f"Pharmacophore detection error: {str(e)}")
            return {}

    def calculate_similarity(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol,
        method: str = "tanimoto",
        fp_type: str = "morgan",
    ) -> Optional[float]:
        """
        Calculate chemical similarity.

        Args:
            mol1: First molecule
            mol2: Second molecule
            method: Similarity method ("tanimoto", "dice", "cosine")
            fp_type: Fingerprint type ("morgan", "maccs", "topological")

        Returns:
            Similarity score between 0 and 1
        """
        try:
            # Generate fingerprints
            if fp_type == "morgan":
                fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
                fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
            elif fp_type == "maccs":
                fp1 = AllChem.GetMACCSKeysFingerprint(mol1)
                fp2 = AllChem.GetMACCSKeysFingerprint(mol2)
            elif fp_type == "topological":
                fp1 = Chem.RDKFingerprint(mol1)
                fp2 = Chem.RDKFingerprint(mol2)
            else:
                raise ValueError(f"Unknown fingerprint type: {fp_type}")

            # Calculate similarity
            if method == "tanimoto":
                return DataStructs.TanimotoSimilarity(fp1, fp2)
            elif method == "dice":
                return DataStructs.DiceSimilarity(fp1, fp2)
            elif method == "cosine":
                return DataStructs.CosineSimilarity(fp1, fp2)
            else:
                raise ValueError(f"Unknown similarity method: {method}")

        except Exception as e:
            self.logger.error(f"Similarity calculation error: {str(e)}")
            return None

    def find_mcs(self, mols: List[Chem.Mol], min_atoms: int = 5) -> Optional[Chem.Mol]:
        """
        Find maximum common substructure.

        Args:
            mols: List of molecules
            min_atoms: Minimum number of atoms in MCS

        Returns:
            MCS as RDKit mol or None
        """
        try:
            if len(mols) < 2:
                return None

            # Find MCS
            mcs = rdDecomposition.FindMCS(
                mols,
                atomCompare=rdDecomposition.AtomCompare.CompareElements,
                bondCompare=rdDecomposition.BondCompare.CompareOrder,
                matchValences=True,
                ringMatchesRingOnly=True,
                completeRingsOnly=True,
            )

            if mcs.numAtoms >= min_atoms:
                return Chem.MolFromSmarts(mcs.smartsString)
            return None

        except Exception as e:
            self.logger.error(f"MCS calculation error: {str(e)}")
            return None

    def generate_depiction(
        self,
        mol: Chem.Mol,
        highlight_atoms: Optional[List[int]] = None,
        highlight_bonds: Optional[List[int]] = None,
        legend: Optional[str] = None,
        size: Tuple[int, int] = (400, 400),
    ) -> Optional[bytes]:
        """
        Generate 2D depiction of molecule.

        Args:
            mol: RDKit molecule
            highlight_atoms: List of atom indices to highlight
            highlight_bonds: List of bond indices to highlight
            legend: Optional structure legend
            size: Image size as (width, height)

        Returns:
            PNG image data as bytes
        """
        try:
            # Set up drawing options
            drawer = Draw.rdDepictor.Compute2DCoords(mol)

            # Create drawing options
            opts = Draw.DrawingOptions()
            opts.addStereoAnnotation = True
            opts.addAtomIndices = True
            opts.legendFontSize = 16

            # Generate image
            if highlight_atoms or highlight_bonds:
                img = Draw.MolToImage(
                    mol,
                    size=size,
                    legend=legend,
                    highlightAtoms=highlight_atoms,
                    highlightBonds=highlight_bonds,
                    highlightColor=(0.2, 0.8, 0.2),
                )
            else:
                img = Draw.MolToImage(
                    mol,
                    size=size,
                    legend=legend,
                )

            # Convert to PNG bytes
            img_bytes = io.BytesIO()
            img.save(img_bytes, format="PNG")
            return img_bytes.getvalue()

        except Exception as e:
            self.logger.error(f"Depiction generation error: {str(e)}")
            return None
