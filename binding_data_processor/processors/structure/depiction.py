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
from rdkit.Chem.AllChem import ChemicalReaction
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
        "white": (1.0, 1.0, 1.0),
        "black": (0.0, 0.0, 0.0),
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
        size: Optional[Tuple[int, int]] = None,
        highlight_atoms: Optional[List[int]] = None,
        highlight_bonds: Optional[List[int]] = None,
        highlight_colors: Optional[Dict[int, Tuple[float, float, float]]] = None,
        legend: Optional[str] = None,
        return_pil: bool = False,
    ) -> Union[str, bytes, Image.Image, None]:
        """Generate depiction of molecule.

        Args:
            mol: RDKit molecule
            size: Image size (width, height)
            highlight_atoms: List of atom indices to highlight
            highlight_bonds: List of bond indices to highlight
            highlight_colors: Dict mapping indices to RGB colors
            legend: Optional legend text
            return_pil: Return PIL Image instead of SVG/bytes

        Returns:
            SVG string, image bytes, PIL Image, or None if error
        """
        try:
            if mol is None:
                return None

            # Generate 2D coordinates if needed
            mol = self.generate_2d_coords(mol)
            if mol is None:
                return None

            # Set up size
            size = size or self.drawing_options["size"]

            # Create drawing object
            if return_pil:
                drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            else:
                drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])

            # Set drawing options
            opts = drawer.drawOptions()
            opts.legendFontSize = self.drawing_options["legend_font_size"]
            opts.atomLabelFontSize = self.drawing_options["atom_label_font_size"]
            opts.bondLineWidth = self.drawing_options["bond_line_width"]
            opts.addAtomIndices = self.drawing_options["add_atom_indices"]
            opts.addBondIndices = self.drawing_options["add_bond_indices"]
            opts.clearBackground = False

            # Set up highlighting
            highlight_atoms = highlight_atoms or []
            highlight_bonds = highlight_bonds or []
            colors = {}
            if highlight_colors:
                for idx, color in highlight_colors.items():
                    colors[idx] = color

            # Draw molecule
            drawer.DrawMolecule(
                mol,
                legend=legend or "",
                highlightAtoms=highlight_atoms,
                highlightBonds=highlight_bonds,
                highlightAtomColors=colors,
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
                color = list(self.HIGHLIGHT_COLORS.values())[i % len(self.HIGHLIGHT_COLORS)]
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
        rxn: ChemicalReaction,
        size: Optional[Tuple[int, int]] = None,
        highlight_by_reactant: bool = True,
        return_pil: bool = False,
    ) -> Union[str, bytes, Image.Image, None]:
        """Generate depiction of chemical reaction.

        Args:
            rxn: RDKit reaction
            size: Image size (width, height)
            highlight_by_reactant: Highlight atoms by reactant
            return_pil: Return PIL Image instead of SVG/bytes

        Returns:
            SVG string, image bytes, PIL Image, or None if error
        """
        try:
            if rxn is None:
                return None

            # Set up size - wider for reactions
            if size is None:
                width = self.drawing_options["size"][0] * 2
                height = self.drawing_options["size"][1]
                size = (width, height)

            # Create drawing object
            if return_pil:
                drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            else:
                drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])

            # Set drawing options
            opts = drawer.drawOptions()
            opts.legendFontSize = self.drawing_options["legend_font_size"]
            opts.atomLabelFontSize = self.drawing_options["atom_label_font_size"]
            opts.bondLineWidth = self.drawing_options["bond_line_width"]
            opts.addAtomIndices = self.drawing_options["add_atom_indices"]
            opts.addBondIndices = self.drawing_options["add_bond_indices"]
            opts.clearBackground = False

            # Generate 2D coordinates for all molecules
            for mol in rxn.GetReactants() + rxn.GetProducts():
                if not mol.GetNumConformers():
                    rdDepictor.Compute2DCoords(mol)

            # Set up highlighting colors
            colors = None
            if highlight_by_reactant:
                colors = []
                for i in range(rxn.GetNumReactantTemplates()):
                    color = list(self.HIGHLIGHT_COLORS.values())[i % len(self.HIGHLIGHT_COLORS)]
                    colors.append(color)

            # Draw reaction
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

    def depict_conformer(
        self,
        mol: Chem.Mol,
        conf_id: int = -1,
        size: Optional[Tuple[int, int]] = None,
        highlight_atoms: Optional[List[int]] = None,
        highlight_bonds: Optional[List[int]] = None,
        highlight_colors: Optional[Dict[int, Tuple[float, float, float]]] = None,
        legend: Optional[str] = None,
        return_pil: bool = False,
    ) -> Union[str, bytes, Image.Image, None]:
        """Generate depiction of 3D conformer.

        Args:
            mol: RDKit molecule
            conf_id: Conformer ID (-1 for current)
            size: Image size (width, height)
            highlight_atoms: List of atom indices to highlight
            highlight_bonds: List of bond indices to highlight
            highlight_colors: Dict mapping indices to RGB colors
            legend: Optional legend text
            return_pil: Return PIL Image instead of SVG/bytes

        Returns:
            SVG string, image bytes, PIL Image, or None if error
        """
        try:
            if mol is None:
                return None

            # Check conformer exists
            if not mol.GetNumConformers():
                return None

            # Set up size
            size = size or self.drawing_options["size"]

            # Create drawing object
            if return_pil:
                drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            else:
                drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])

            # Set drawing options
            opts = drawer.drawOptions()
            opts.legendFontSize = self.drawing_options["legend_font_size"]
            opts.atomLabelFontSize = self.drawing_options["atom_label_font_size"]
            opts.bondLineWidth = self.drawing_options["bond_line_width"]
            opts.addAtomIndices = self.drawing_options["add_atom_indices"]
            opts.addBondIndices = self.drawing_options["add_bond_indices"]
            opts.clearBackground = False

            # Set up highlighting
            highlight_atoms = highlight_atoms or []
            highlight_bonds = highlight_bonds or []
            colors = {}
            if highlight_colors:
                for idx, color in highlight_colors.items():
                    colors[idx] = color

            # Draw conformer
            drawer.DrawMolecule(
                mol,
                confId=conf_id,
                highlightAtoms=highlight_atoms,
                highlightBonds=highlight_bonds,
                highlightAtomColors=colors,
            )
            drawer.FinishDrawing()

            # Get image data
            if return_pil:
                img_data = drawer.GetDrawingText()
                return Image.open(io.BytesIO(img_data))
            else:
                return drawer.GetDrawingText()

        except Exception as e:
            self.logger.error(f"Error depicting conformer: {str(e)}")
            return None

    def depict_grid(
        self,
        mols: List[Chem.Mol],
        legends: Optional[List[str]] = None,
        mols_per_row: int = 3,
        sub_img_size: Optional[Tuple[int, int]] = None,
        return_pil: bool = False,
    ) -> Union[bytes, Image.Image, None]:
        """Generate grid depiction of multiple molecules.

        Args:
            mols: List of molecules
            legends: Optional list of legends
            mols_per_row: Number of molecules per row
            sub_img_size: Size of each molecule image
            return_pil: Return PIL Image instead of bytes

        Returns:
            Image bytes, PIL Image, or None if error
        """
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

            # Set up sub-image size
            sub_img_size = sub_img_size or (200, 200)

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
    ) -> Union[bytes, Image.Image, Tuple[bytes, np.ndarray], None]:
        """Generate 2D similarity map using ML dimensionality reduction.

        Args:
            mols: List of molecules
            method: Dimensionality reduction method ('tsne' or 'pca')
            return_pil: Return PIL Image instead of bytes

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
            if return_pil:
                drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
            else:
                drawer = rdMolDraw2D.MolDraw2DSVG(width, height)

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
        """Save image data to file.

        Args:
            image_data: Image data to save
            filename: Output filename
            img_format: Image format (e.g., 'png', 'svg')

        Returns:
            Success status
        """
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
