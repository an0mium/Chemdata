"""Molecular structure depiction and visualization.

This module provides comprehensive functionality for:
1. 2D and 3D structure visualization
2. Substructure highlighting and pattern matching
3. Grid layouts and reaction depiction
4. ML-enhanced visualization features
5. Image generation and export
"""

import io
import logging
from typing import Dict, List, Optional, Set, Tuple, Union
import numpy as np
from PIL import Image
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    AllChem,
    Draw,
    rdDepictor,
    rdMolDescriptors,
    rdMolTransforms,
)
from rdkit.Chem.Draw import rdMolDraw2D, IPythonConsole

# Optional ML imports
try:
    from sklearn.manifold import TSNE
    from sklearn.decomposition import PCA
    from sklearn.cluster import DBSCAN

    ML_AVAILABLE = True
except ImportError:
    ML_AVAILABLE = False


class StructureDepiction:
    """Handles molecular structure depiction and visualization."""

    # Drawing options
    DRAWING_OPTIONS = {
        "size": (400, 400),
        "legend_font_size": 12,
        "atom_label_font_size": 10,
        "bond_line_width": 2,
        "add_atom_indices": False,
        "add_bond_indices": False,
        "highlight_color": (0.7, 0.7, 1.0),
        "background_color": (1.0, 1.0, 1.0),
    }

    # Highlight colors
    HIGHLIGHT_COLORS = {
        "red": (1.0, 0.0, 0.0),
        "green": (0.0, 1.0, 0.0),
        "blue": (0.0, 0.0, 1.0),
        "yellow": (1.0, 1.0, 0.0),
        "cyan": (0.0, 1.0, 1.0),
        "magenta": (1.0, 0.0, 1.0),
        "orange": (1.0, 0.5, 0.0),
        "purple": (0.5, 0.0, 0.5),
    }

    def __init__(self, config: Optional[Dict] = None):
        """Initialize structure depiction."""
        self.logger = logging.getLogger(__name__)
        self.config = config or {}
        self.drawing_options = self.DRAWING_OPTIONS.copy()

        # Initialize ML components if available
        if ML_AVAILABLE:
            self.tsne = TSNE(n_components=2, random_state=42)
            self.pca = PCA(n_components=2)
            self.clustering = DBSCAN(eps=0.5, min_samples=5)

    def set_drawing_options(self, **kwargs):
        """Set drawing options."""
        self.drawing_options.update(kwargs)

    def generate_2d_coords(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Generate 2D coordinates for molecule."""
        try:
            if mol is None:
                return None

            mol = Chem.Mol(mol)
            if not mol.GetNumConformers():
                rdDepictor.Compute2DCoords(mol)
            return mol

        except Exception as e:
            self.logger.error(f"Error generating 2D coords: {str(e)}")
            return None

    def depict_molecule(
        self,
        mol: Chem.Mol,
        highlight_atoms: Optional[List[int]] = None,
        highlight_bonds: Optional[List[int]] = None,
        highlight_color: Optional[Tuple[float, float, float]] = None,
        legend: Optional[str] = None,
        return_pil: bool = False,
    ) -> Union[bytes, Image.Image]:
        """Generate 2D depiction of molecule."""
        try:
            if mol is None:
                return None

            # Generate 2D coordinates if needed
            mol = self.generate_2d_coords(mol)
            if mol is None:
                return None

            # Create drawing object
            drawer = rdMolDraw2D.MolDraw2DCairo(
                self.drawing_options["size"][0],
                self.drawing_options["size"][1],
            )

            # Set drawing options
            opts = drawer.drawOptions()
            opts.legendFontSize = self.drawing_options["legend_font_size"]
            opts.atomLabelFontSize = self.drawing_options["atom_label_font_size"]
            opts.bondLineWidth = self.drawing_options["bond_line_width"]
            opts.addAtomIndices = self.drawing_options["add_atom_indices"]
            opts.addBondIndices = self.drawing_options["add_bond_indices"]
            opts.clearBackground = False

            # Set highlighting
            highlight_atoms = highlight_atoms or []
            highlight_bonds = highlight_bonds or []
            highlight_color = highlight_color or self.drawing_options["highlight_color"]
            highlight_radii = {i: 0.5 for i in highlight_atoms}

            # Draw molecule
            drawer.DrawMolecule(
                mol,
                legend=legend or "",
                highlightAtoms=highlight_atoms,
                highlightBonds=highlight_bonds,
                highlightColor=highlight_color,
                highlightRadii=highlight_radii,
            )
            drawer.FinishDrawing()

            # Get image data
            if return_pil:
                img_data = drawer.GetDrawingText()
                return Image.open(io.BytesIO(img_data))
            else:
                return drawer.GetDrawingText()

        except Exception as e:
            self.logger.error(f"Error depicting molecule: {str(e)}")
            return None

    def depict_substructure_match(
        self,
        mol: Chem.Mol,
        pattern: Union[str, Chem.Mol],
        all_matches: bool = False,
        return_pil: bool = False,
    ) -> Union[bytes, Image.Image]:
        """Generate depiction with highlighted substructure match."""
        try:
            if mol is None:
                return None

            # Convert pattern to query molecule if needed
            if isinstance(pattern, str):
                query = Chem.MolFromSmarts(pattern)
                if query is None:
                    return None
            else:
                query = pattern

            # Find matches
            if all_matches:
                matches = mol.GetSubstructMatches(query)
            else:
                match = mol.GetSubstructMatch(query)
                matches = [match] if match else []

            if not matches:
                return self.depict_molecule(mol, return_pil=return_pil)

            # Highlight each match with a different color
            highlight_atoms = set()
            highlight_bonds = set()
            colors = {}

            for i, match in enumerate(matches):
                color = list(self.HIGHLIGHT_COLORS.values())[
                    i % len(self.HIGHLIGHT_COLORS)
                ]
                for atom_idx in match:
                    highlight_atoms.add(atom_idx)
                    colors[atom_idx] = color

                # Find bonds between matched atoms
                for bond in query.GetBonds():
                    aid1 = match[bond.GetBeginAtomIdx()]
                    aid2 = match[bond.GetEndAtomIdx()]
                    bond_idx = mol.GetBondBetweenAtoms(aid1, aid2).GetIdx()
                    highlight_bonds.add(bond_idx)

            return self.depict_molecule(
                mol,
                highlight_atoms=list(highlight_atoms),
                highlight_bonds=list(highlight_bonds),
                highlight_color=list(colors.values())[0],
                return_pil=return_pil,
            )

        except Exception as e:
            self.logger.error(f"Error depicting substructure match: {str(e)}")
            return None

    def depict_reaction(
        self,
        rxn: Chem.ChemicalReaction,
        highlight_by_reactant: bool = True,
        return_pil: bool = False,
    ) -> Union[bytes, Image.Image]:
        """Generate depiction of chemical reaction."""
        try:
            if rxn is None:
                return None

            # Create drawing object with larger size for reactions
            width = self.drawing_options["size"][0] * 2
            height = self.drawing_options["size"][1]
            drawer = rdMolDraw2D.MolDraw2DCairo(width, height)

            # Set drawing options
            opts = drawer.drawOptions()
            opts.legendFontSize = self.drawing_options["legend_font_size"]
            opts.atomLabelFontSize = self.drawing_options["atom_label_font_size"]
            opts.bondLineWidth = self.drawing_options["bond_line_width"]
            opts.addAtomIndices = self.drawing_options["add_atom_indices"]
            opts.addBondIndices = self.drawing_options["add_bond_indices"]
            opts.clearBackground = False

            # Draw reaction
            colors = None
            if highlight_by_reactant:
                colors = []
                for i in range(rxn.GetNumReactantTemplates()):
                    color = list(self.HIGHLIGHT_COLORS.values())[
                        i % len(self.HIGHLIGHT_COLORS)
                    ]
                    colors.append(color)

            drawer.DrawReaction(rxn, highlightByReactant=highlight_by_reactant)
            drawer.FinishDrawing()

            # Get image data
            if return_pil:
                img_data = drawer.GetDrawingText()
                return Image.open(io.BytesIO(img_data))
            else:
                return drawer.GetDrawingText()

        except Exception as e:
            self.logger.error(f"Error depicting reaction: {str(e)}")
            return None

    def depict_grid(
        self,
        mols: List[Chem.Mol],
        legends: Optional[List[str]] = None,
        mols_per_row: int = 3,
        sub_img_size: Tuple[int, int] = (200, 200),
        return_pil: bool = False,
    ) -> Union[bytes, Image.Image]:
        """Generate grid depiction of multiple molecules."""
        try:
            if not mols:
                return None

            # Filter out None molecules
            valid_mols = []
            valid_legends = []
            for i, mol in enumerate(mols):
                if mol is not None:
                    valid_mols.append(mol)
                    if legends and i < len(legends):
                        valid_legends.append(legends[i])
                    else:
                        valid_legends.append("")

            if not valid_mols:
                return None

            # Create grid image
            img = Draw.MolsToGridImage(
                valid_mols,
                legends=valid_legends,
                molsPerRow=mols_per_row,
                subImgSize=sub_img_size,
                returnPNG=not return_pil,
            )

            return img

        except Exception as e:
            self.logger.error(f"Error depicting molecule grid: {str(e)}")
            return None

    def depict_similarity_map(
        self,
        mols: List[Chem.Mol],
        method: str = "tsne",
        return_pil: bool = False,
    ) -> Union[bytes, Image.Image, Tuple[bytes, np.ndarray]]:
        """
        Generate 2D similarity map using ML dimensionality reduction.

        Args:
            mols: List of molecules
            method: Dimensionality reduction method ('tsne' or 'pca')
            return_pil: Return PIL Image

        Returns:
            Image data and optionally coordinates
        """
        if not ML_AVAILABLE:
            self.logger.warning("ML functionality not available")
            return None

        try:
            if not mols:
                return None

            # Generate fingerprints
            fps = []
            valid_mols = []
            for mol in mols:
                if mol is not None:
                    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2)
                    if fp is not None:
                        fps.append(fp)
                        valid_mols.append(mol)

            if not fps:
                return None

            # Convert fingerprints to numpy array
            fp_array = []
            for fp in fps:
                arr = np.zeros((1,))
                DataStructs.ConvertToNumpyArray(fp, arr)
                fp_array.append(arr)
            X = np.vstack(fp_array)

            # Reduce dimensionality
            if method == "tsne":
                coords = self.tsne.fit_transform(X)
            else:
                coords = self.pca.fit_transform(X)

            # Scale coordinates to image size
            x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
            y_min, y_max = coords[:, 1].min(), coords[:, 1].max()

            scaled_coords = np.zeros_like(coords)
            scaled_coords[:, 0] = (coords[:, 0] - x_min) / (x_max - x_min) * 0.8 + 0.1
            scaled_coords[:, 1] = (coords[:, 1] - y_min) / (y_max - y_min) * 0.8 + 0.1

            # Create image
            width, height = self.drawing_options["size"]
            drawer = rdMolDraw2D.MolDraw2DCairo(width, height)

            # Draw molecules at scaled coordinates
            for mol, (x, y) in zip(valid_mols, scaled_coords):
                drawer.SetOffset(int(x * width), int(y * height))
                drawer.DrawMolecule(mol)

            drawer.FinishDrawing()

            # Get image data
            if return_pil:
                img_data = drawer.GetDrawingText()
                return Image.open(io.BytesIO(img_data)), coords
            else:
                return drawer.GetDrawingText(), coords

        except Exception as e:
            self.logger.error(f"Error generating similarity map: {str(e)}")
            return None

    def save_image(
        self,
        image_data: Union[bytes, Image.Image],
        filename: str,
        img_format: str = "png",
    ) -> bool:
        """Save image data to file."""
        try:
            if image_data is None:
                return False

            if isinstance(image_data, Image.Image):
                image_data.save(filename, format=img_format.upper())
            else:
                with open(filename, "wb") as f:
                    f.write(image_data)

            return True

        except Exception as e:
            self.logger.error(f"Error saving image: {str(e)}")
            return False
