"""Visualization utilities for binding site analysis."""

import io
import logging
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from PIL import Image
from Bio.PDB import Structure, Residue
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


class SiteVisualizer:
    """Visualize binding sites and their properties."""

    def __init__(self):
        """Initialize site visualizer."""
        self.logger = logging.getLogger(self.__class__.__name__)

        # Color schemes for different visualization types
        self.color_schemes = {
            "conservation": {
                "low": (1.0, 0.0, 0.0),  # Red
                "medium": (1.0, 1.0, 0.0),  # Yellow
                "high": (0.0, 1.0, 0.0),  # Green
            },
            "pharmacophore": {
                "hydrophobic": (0.7, 0.7, 0.7),  # Gray
                "hbond_donor": (0.0, 0.0, 1.0),  # Blue
                "hbond_acceptor": (1.0, 0.0, 0.0),  # Red
                "aromatic": (1.0, 0.5, 0.0),  # Orange
                "charged_pos": (0.0, 0.0, 1.0),  # Blue
                "charged_neg": (1.0, 0.0, 0.0),  # Red
            },
            "score": {
                "low": (1.0, 0.0, 0.0),  # Red
                "medium": (1.0, 1.0, 0.0),  # Yellow
                "high": (0.0, 1.0, 0.0),  # Green
            },
        }

    def visualize_binding_site(
        self,
        structure: Structure.Structure,
        site_residues: List[Residue.Residue],
        ligand: Optional[Chem.Mol] = None,
        color_by: str = "score",
        scores: Optional[Dict[str, float]] = None,
        output_path: Optional[str] = None,
        width: int = 800,
        height: int = 800,
    ) -> Optional[Image.Image]:
        """Generate visualization of binding site.

        Args:
            structure: Full protein structure
            site_residues: List of binding site residues
            ligand: Optional ligand molecule
            color_by: How to color residues ('conservation', 'pharmacophore', 'score')
            scores: Optional dict of residue scores for coloring
            output_path: Path to save image, or None for in-memory
            width: Image width in pixels
            height: Image height in pixels

        Returns:
            PIL Image if output_path is None, else None
        """
        try:
            # This is a placeholder that would use PyMOL or similar
            # to generate actual 3D visualization
            self.logger.warning("3D visualization not implemented")
            return None

        except Exception as e:
            self.logger.error(f"Error visualizing binding site: {str(e)}")
            return None

    def draw_site(
        self,
        mol: Chem.Mol,
        site_atoms: List[int],
        highlight_colors: Optional[Dict[int, Tuple[float, float, float]]] = None,
        size: Tuple[int, int] = (400, 400),
        legend: Optional[str] = None,
    ) -> Optional[Image.Image]:
        """Draw molecule with highlighted binding site.

        Args:
            mol: RDKit molecule
            site_atoms: List of atom indices in binding site
            highlight_colors: Optional color mapping for atoms
            size: Image size (width, height)
            legend: Optional legend text

        Returns:
            PIL Image or None if error
        """
        try:
            if mol is None:
                return None

            # Generate 2D coordinates if needed
            if not mol.GetNumConformers():
                AllChem.Compute2DCoords(mol)

            # Set up drawing
            drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            opts = drawer.drawOptions()
            opts.legendFontSize = 12
            opts.atomLabelFontSize = 10
            opts.bondLineWidth = 2
            opts.addAtomIndices = False
            opts.addBondIndices = False
            opts.clearBackground = False

            # Set up highlighting
            highlight_bonds = []
            if highlight_colors is None:
                highlight_colors = {}
                for idx in site_atoms:
                    highlight_colors[idx] = (0.7, 0.7, 1.0)

            # Find bonds between highlighted atoms
            for bond in mol.GetBonds():
                if bond.GetBeginAtomIdx() in site_atoms and bond.GetEndAtomIdx() in site_atoms:
                    highlight_bonds.append(bond.GetIdx())

            # Draw molecule
            drawer.DrawMolecule(
                mol,
                legend=legend or "",
                highlightAtoms=site_atoms,
                highlightBonds=highlight_bonds,
                highlightAtomColors=highlight_colors,
            )
            drawer.FinishDrawing()

            # Convert to PIL Image
            img_data = drawer.GetDrawingText()
            return Image.open(io.BytesIO(img_data))

        except Exception as e:
            self.logger.error(f"Error drawing binding site: {str(e)}")
            return None

    def draw_interaction_map(
        self,
        mol: Chem.Mol,
        interactions: Dict[int, List[Dict]],
        size: Tuple[int, int] = (400, 400),
    ) -> Optional[Image.Image]:
        """Draw molecule with interaction map.

        Args:
            mol: RDKit molecule
            interactions: Dictionary mapping atom indices to interaction types
            size: Image size (width, height)

        Returns:
            PIL Image or None if error
        """
        try:
            if mol is None:
                return None

            # Generate 2D coordinates if needed
            if not mol.GetNumConformers():
                AllChem.Compute2DCoords(mol)

            # Set up drawing
            drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            opts = drawer.drawOptions()
            opts.legendFontSize = 12
            opts.atomLabelFontSize = 10
            opts.bondLineWidth = 2
            opts.addAtomIndices = False
            opts.addBondIndices = False
            opts.clearBackground = False

            # Set up highlighting colors
            highlight_colors = {}
            for atom_idx, atom_interactions in interactions.items():
                # Color based on interaction type
                if any(i["type"] == "hbond_donor" for i in atom_interactions):
                    highlight_colors[atom_idx] = self.color_schemes["pharmacophore"]["hbond_donor"]
                elif any(i["type"] == "hbond_acceptor" for i in atom_interactions):
                    highlight_colors[atom_idx] = self.color_schemes["pharmacophore"]["hbond_acceptor"]
                elif any(i["type"] == "aromatic" for i in atom_interactions):
                    highlight_colors[atom_idx] = self.color_schemes["pharmacophore"]["aromatic"]
                elif any(i["type"] == "hydrophobic" for i in atom_interactions):
                    highlight_colors[atom_idx] = self.color_schemes["pharmacophore"]["hydrophobic"]
                else:
                    highlight_colors[atom_idx] = (0.7, 0.7, 0.7)  # Gray

            # Draw molecule
            drawer.DrawMolecule(
                mol,
                highlightAtoms=list(interactions.keys()),
                highlightAtomColors=highlight_colors,
            )
            drawer.FinishDrawing()

            # Convert to PIL Image
            img_data = drawer.GetDrawingText()
            return Image.open(io.BytesIO(img_data))

        except Exception as e:
            self.logger.error(f"Error drawing interaction map: {str(e)}")
            return None

    def draw_pharmacophore(
        self,
        mol: Chem.Mol,
        features: List[Dict],
        size: Tuple[int, int] = (400, 400),
    ) -> Optional[Image.Image]:
        """Draw molecule with pharmacophore features.

        Args:
            mol: RDKit molecule
            features: List of pharmacophore feature dictionaries
            size: Image size (width, height)

        Returns:
            PIL Image or None if error
        """
        try:
            if mol is None:
                return None

            # Generate 2D coordinates if needed
            if not mol.GetNumConformers():
                AllChem.Compute2DCoords(mol)

            # Set up drawing
            drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            opts = drawer.drawOptions()
            opts.legendFontSize = 12
            opts.atomLabelFontSize = 10
            opts.bondLineWidth = 2
            opts.addAtomIndices = False
            opts.addBondIndices = False
            opts.clearBackground = False

            # Set up highlighting colors
            highlight_atoms = []
            highlight_colors = {}
            for feature in features:
                atoms = feature["atoms"]
                highlight_atoms.extend(atoms)
                color = self.color_schemes["pharmacophore"].get(feature["type"], (0.7, 0.7, 0.7))
                for atom_idx in atoms:
                    highlight_colors[atom_idx] = color

            # Draw molecule
            drawer.DrawMolecule(
                mol,
                highlightAtoms=highlight_atoms,
                highlightAtomColors=highlight_colors,
            )
            drawer.FinishDrawing()

            # Convert to PIL Image
            img_data = drawer.GetDrawingText()
            return Image.open(io.BytesIO(img_data))

        except Exception as e:
            self.logger.error(f"Error drawing pharmacophore: {str(e)}")
            return None

    def get_residue_color(
        self,
        residue: Residue.Residue,
        color_by: str,
        score: Optional[float] = None,
    ) -> Tuple[float, float, float]:
        """Get color for a residue based on properties/scores.

        Args:
            residue: Residue to color
            color_by: How to color ('conservation', 'pharmacophore', 'score')
            score: Optional score value for coloring

        Returns:
            RGB color tuple
        """
        if color_by == "score" and score is not None:
            if score < 0.4:
                return self.color_schemes["score"]["low"]
            elif score < 0.7:
                return self.color_schemes["score"]["medium"]
            else:
                return self.color_schemes["score"]["high"]

        elif color_by == "pharmacophore":
            resname = residue.get_resname()
            if resname in {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}:
                return self.color_schemes["pharmacophore"]["hydrophobic"]
            elif resname in {"LYS", "ARG", "HIS"}:
                return self.color_schemes["pharmacophore"]["charged_pos"]
            elif resname in {"ASP", "GLU"}:
                return self.color_schemes["pharmacophore"]["charged_neg"]
            elif resname in {"SER", "THR", "ASN", "GLN", "TYR", "TRP"}:
                return self.color_schemes["pharmacophore"]["hbond_donor"]
            elif resname in {"PHE", "TYR", "TRP", "HIS"}:
                return self.color_schemes["pharmacophore"]["aromatic"]

        # Default color
        return (0.5, 0.5, 0.5)  # Gray

    def create_ligand_interaction_diagram(
        self,
        ligand: Chem.Mol,
        site_residues: List[Residue.Residue],
        interactions: Dict[str, List[Tuple[int, str]]],
        output_path: Optional[str] = None,
        width: int = 800,
        height: int = 800,
    ) -> Optional[Image.Image]:
        """Create 2D diagram showing ligand-residue interactions.

        Args:
            ligand: Ligand molecule
            site_residues: List of binding site residues
            interactions: Dict mapping interaction types to lists of
                        (atom_idx, residue_name) tuples
            output_path: Path to save image, or None for in-memory
            width: Image width in pixels
            height: Image height in pixels

        Returns:
            PIL Image if output_path is None, else None
        """
        try:
            # This is a placeholder that would use RDKit's drawing
            # capabilities to create interaction diagrams
            self.logger.warning("Interaction diagram not implemented")
            return None

        except Exception as e:
            self.logger.error(f"Error creating interaction diagram: {str(e)}")
            return None

    def create_pharmacophore_diagram(
        self,
        site_residues: List[Residue.Residue],
        features: Dict[str, List[Tuple[float, float, float]]],
        output_path: Optional[str] = None,
        width: int = 800,
        height: int = 800,
    ) -> Optional[Image.Image]:
        """Create 3D pharmacophore feature diagram.

        Args:
            site_residues: List of binding site residues
            features: Dict mapping feature types to lists of 3D coordinates
            output_path: Path to save image, or None for in-memory
            width: Image width in pixels
            height: Image height in pixels

        Returns:
            PIL Image if output_path is None, else None
        """
        try:
            # This is a placeholder that would use PyMOL or similar
            # to create pharmacophore visualizations
            self.logger.warning("Pharmacophore diagram not implemented")
            return None

        except Exception as e:
            self.logger.error(f"Error creating pharmacophore diagram: {str(e)}")
            return None

    def save_image(
        self,
        image: Image.Image,
        path: str,
        format: str = "PNG",
    ) -> bool:
        """Save image to file.

        Args:
            image: PIL Image to save
            path: Output path
            format: Image format

        Returns:
            True if successful
        """
        try:
            if image is None:
                return False
            image.save(path, format=format)
            return True
        except Exception as e:
            self.logger.error(f"Error saving image: {str(e)}")
            return False
