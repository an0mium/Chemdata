"""2D molecular structure depiction and visualization.

This module provides functionality for:
1. 2D structure visualization
2. Substructure highlighting
3. Atom and bond annotations
4. Custom coloring and styling
5. Grid and template layouts
"""

import io
import logging
from typing import Dict, List, Optional, Tuple, Union
import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Draw,
    rdDepictor,
)
from rdkit.Chem.Draw import rdMolDraw2D

from .base_depiction import BaseStructureDepiction


class Structure2DDepiction(BaseStructureDepiction):
    """Handles 2D molecular structure depiction."""

    # Additional 2D-specific options
    DRAWING_OPTIONS_2D = {
        "kekulize": True,
        "wedge_bonds": True,
        "show_atom_numbers": False,
        "show_bond_numbers": False,
        "annotate_stereo": True,
        "use_svg": False,
    }

    def __init__(self, config: Optional[Dict] = None):
        """Initialize 2D structure depiction."""
        super().__init__(config)
        # Add 2D-specific options
        self.drawing_options.update(self.DRAWING_OPTIONS_2D)

    def prepare_mol_2d(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Prepare molecule for 2D depiction.

        Args:
            mol: Input molecule

        Returns:
            Prepared molecule
        """
        try:
            if not self.validate_mol(mol):
                return None

            # Add hydrogens if configured
            mol = self.add_hydrogens(mol)
            if mol is None:
                return None

            # Kekulize if configured
            if self.drawing_options["kekulize"]:
                Chem.Kekulize(mol)

            # Generate 2D coordinates if needed
            if not mol.GetNumConformers():
                rdDepictor.Compute2DCoords(mol)

            # Wedge bonds if configured
            if self.drawing_options["wedge_bonds"]:
                Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

            return mol

        except Exception as e:
            self.logger.error(f"Error preparing molecule: {str(e)}")
            return None

    def depict_2d(
        self,
        mol: Chem.Mol,
        highlight_atoms: Optional[List[int]] = None,
        highlight_bonds: Optional[List[int]] = None,
        atom_colors: Optional[Dict[int, Tuple[float, float, float]]] = None,
        bond_colors: Optional[Dict[int, Tuple[float, float, float]]] = None,
        atom_labels: Optional[Dict[int, str]] = None,
        return_pil: bool = False,
    ) -> Union[bytes, Image.Image]:
        """Generate 2D depiction of molecule.

        Args:
            mol: Molecule to depict
            highlight_atoms: Optional list of atom indices to highlight
            highlight_bonds: Optional list of bond indices to highlight
            atom_colors: Optional dict mapping atom indices to RGB colors
            bond_colors: Optional dict mapping bond indices to RGB colors
            atom_labels: Optional dict mapping atom indices to custom labels
            return_pil: Whether to return PIL Image

        Returns:
            Image data
        """
        try:
            mol = self.prepare_mol_2d(mol)
            if mol is None:
                return None

            # Create drawing object
            if self.drawing_options["use_svg"]:
                drawer = rdMolDraw2D.MolDraw2DSVG(
                    self.drawing_options["size"][0],
                    self.drawing_options["size"][1],
                )
            else:
                drawer = rdMolDraw2D.MolDraw2DCairo(
                    self.drawing_options["size"][0],
                    self.drawing_options["size"][1],
                )

            # Set drawing options
            opts = drawer.drawOptions()
            for k, v in self.get_drawing_options_dict().items():
                setattr(opts, k, v)

            # Additional 2D options
            opts.addAtomIndices = self.drawing_options["show_atom_numbers"]
            opts.addBondIndices = self.drawing_options["show_bond_numbers"]
            opts.includeAtomTags = bool(atom_labels)

            # Draw molecule
            drawer.DrawMolecule(
                mol,
                highlightAtoms=highlight_atoms,
                highlightBonds=highlight_bonds,
                highlightAtomColors=atom_colors,
                highlightBondColors=bond_colors,
                atomLabels=atom_labels,
            )
            drawer.FinishDrawing()

            # Get image data
            img_data = drawer.GetDrawingText()
            return self.create_image(img_data, return_pil)

        except Exception as e:
            self.logger.error(f"Error depicting 2D structure: {str(e)}")
            return None

    def depict_grid(
        self,
        mols: List[Chem.Mol],
        legends: Optional[List[str]] = None,
        mols_per_row: int = 3,
        sub_img_size: Optional[Tuple[int, int]] = None,
        return_pil: bool = False,
    ) -> Union[bytes, Image.Image]:
        """Generate grid depiction of multiple molecules.

        Args:
            mols: List of molecules to depict
            legends: Optional list of legends for each molecule
            mols_per_row: Number of molecules per row
            sub_img_size: Optional size for each sub-image
            return_pil: Whether to return PIL Image

        Returns:
            Grid image data
        """
        try:
            if not mols:
                return None

            # Prepare molecules
            prepared_mols = []
            for mol in mols:
                prep_mol = self.prepare_mol_2d(mol)
                if prep_mol is not None:
                    prepared_mols.append(prep_mol)

            if not prepared_mols:
                return None

            # Create grid image
            img = Draw.MolsToGridImage(
                prepared_mols,
                legends=legends,
                molsPerRow=mols_per_row,
                subImgSize=sub_img_size or self.drawing_options["size"],
                returnPNG=not return_pil,
            )

            return img

        except Exception as e:
            self.logger.error(f"Error depicting molecule grid: {str(e)}")
            return None

    def depict_with_highlights(
        self,
        mol: Chem.Mol,
        query: Chem.Mol,
        return_pil: bool = False,
    ) -> Union[bytes, Image.Image]:
        """Generate depiction with highlighted substructure match.

        Args:
            mol: Molecule to depict
            query: Query substructure to highlight
            return_pil: Whether to return PIL Image

        Returns:
            Image data with highlighted substructure
        """
        try:
            mol = self.prepare_mol_2d(mol)
            if mol is None:
                return None

            # Find matches
            matches = mol.GetSubstructMatches(query)
            if not matches:
                return self.depict_2d(mol, return_pil=return_pil)

            # Get atoms and bonds to highlight
            highlight_atoms = set()
            highlight_bonds = set()
            match_colors = {}

            for i, match in enumerate(matches):
                # Generate distinct color for each match
                hue = i / len(matches)
                color = self._hsl_to_rgb(hue, 1.0, 0.5)

                # Add atoms
                for atom_idx in match:
                    highlight_atoms.add(atom_idx)
                    match_colors[atom_idx] = color

                # Add bonds between matched atoms
                for bond in mol.GetBonds():
                    begin_idx = bond.GetBeginAtomIdx()
                    end_idx = bond.GetEndAtomIdx()
                    if begin_idx in match and end_idx in match:
                        highlight_bonds.add(bond.GetIdx())

            return self.depict_2d(
                mol,
                highlight_atoms=list(highlight_atoms),
                highlight_bonds=list(highlight_bonds),
                atom_colors=match_colors,
                return_pil=return_pil,
            )

        except Exception as e:
            self.logger.error(f"Error depicting with highlights: {str(e)}")
            return None

    def _hsl_to_rgb(self, h: float, s: float, l: float) -> Tuple[float, float, float]:
        """Convert HSL color to RGB.

        Args:
            h: Hue (0-1)
            s: Saturation (0-1)
            l: Lightness (0-1)

        Returns:
            RGB color tuple
        """

        def hue_to_rgb(p: float, q: float, t: float) -> float:
            if t < 0:
                t += 1
            if t > 1:
                t -= 1
            if t < 1 / 6:
                return p + (q - p) * 6 * t
            if t < 1 / 2:
                return q
            if t < 2 / 3:
                return p + (q - p) * (2 / 3 - t) * 6
            return p

        if s == 0:
            return (l, l, l)

        q = l * (1 + s) if l < 0.5 else l + s - l * s
        p = 2 * l - q

        r = hue_to_rgb(p, q, h + 1 / 3)
        g = hue_to_rgb(p, q, h)
        b = hue_to_rgb(p, q, h - 1 / 3)

        return (r, g, b)
