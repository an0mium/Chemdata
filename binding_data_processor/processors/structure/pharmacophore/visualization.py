"""2D and 3D visualization methods for pharmacophore analysis with enhanced RDKit features."""

from typing import Dict, List, Optional, Tuple, Union
import json

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdDepictor, rdMolTransforms
from rdkit.Chem.Draw import rdMolDraw2D, IPythonConsole
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Chem.Draw.rdDepictor import Compute2DCoords, SetPreferCoordGen
from rdkit.Geometry import Point3D
from rdkit.Chem import rdShapeHelpers, rdMolDescriptors, rdMolTransforms

from .features import COLOR_SCHEMES, PharmacophoreFeature

# Enable modern 2D coordinate generation
SetPreferCoordGen(True)


class PharmacophoreVisualizer:
    """Enhanced class for visualizing pharmacophore features and alignments."""

    def __init__(self, color_scheme: str = "default"):
        """Initialize the pharmacophore visualizer.

        Args:
            color_scheme: Color scheme to use ("default" or "colorblind")
        """
        self.colors = COLOR_SCHEMES[color_scheme]
        self.drawing_options = self._get_enhanced_drawing_options()

    def _get_enhanced_drawing_options(self) -> DrawingOptions:
        """Get enhanced drawing options for better visualization."""
        opts = DrawingOptions()
        opts.atomLabelFontSize = 14
        opts.bondLineWidth = 2.0
        opts.dblBondOffset = 0.3
        opts.includeAtomNumbers = False
        opts.additionalAtomLabelPadding = 0.3
        opts.explicitMethyl = False
        opts.coordScale = 1.5
        opts.comicMode = False
        opts.addStereoAnnotation = True
        opts.addAtomIndices = False
        opts.legendFontSize = 16
        opts.multipleBondOffset = 0.3
        opts.padding = 0.05
        return opts

    def draw_2d_with_features(
        self,
        mol: Chem.Mol,
        features: List[PharmacophoreFeature],
        size: Tuple[int, int] = (400, 400),
        legend: bool = True,
        interactive: bool = True,
        include_atom_indices: bool = False,
        highlight_style: str = "circles",
        show_vectors: bool = True,
    ) -> str:
        """Draw molecule with highlighted pharmacophore features in 2D.

        Args:
            mol: Molecule to visualize
            features: List of pharmacophore features
            size: Image size in pixels
            legend: Whether to include feature legend
            interactive: Whether to generate interactive SVG with tooltips
            include_atom_indices: Whether to show atom indices
            highlight_style: Style for highlighting ("circles", "radial", or "filled")
            show_vectors: Whether to show feature vectors

        Returns:
            SVG string of the drawing
        """
        # Generate 2D coordinates if needed
        if not mol.GetNumConformers():
            rdDepictor.Compute2DCoords(mol)

        # Create drawing object
        drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
        for opt_name, opt_value in vars(self.drawing_options).items():
            if hasattr(drawer.drawOptions(), opt_name):
                setattr(drawer.drawOptions(), opt_name, opt_value)
        drawer.drawOptions().addAtomIndices = include_atom_indices

        # Prepare atom highlights
        atom_colors = {}
        atom_radii = {}
        legend_entries = set()
        atom_tooltips = {}
        vectors = []

        for feature in features:
            if not feature.enabled:
                continue

            color = self.colors.get(feature.feature_type, (0.5, 0.5, 0.5))
            for atom_idx in feature.atoms:
                atom_colors[atom_idx] = color
                atom_radii[atom_idx] = feature.radius * 0.5
                if interactive:
                    atom_tooltips[atom_idx] = f"{feature.feature_type} ({', '.join(map(str, feature.atoms))})"
                if legend:
                    legend_entries.add(feature.feature_type)

            # Add feature vector if available
            if show_vectors and feature.vector and feature.position:
                vectors.append((Point3D(*feature.position), Point3D(*(p + v * 2.0 for p, v in zip(feature.position, feature.vector))), color))

        # Draw molecule with highlights
        if highlight_style == "radial":
            drawer.drawOptions().fillHighlights = True
            for atom_idx, color in atom_colors.items():
                drawer.SetFillRadialGradient(atom_idx, color, (1, 1, 1), 0.3)
        elif highlight_style == "filled":
            drawer.drawOptions().fillHighlights = True

        drawer.DrawMolecule(mol, highlightAtoms=list(atom_colors.keys()), highlightAtomColors=atom_colors, highlightAtomRadii=atom_radii)

        # Draw feature vectors
        for start, end, color in vectors:
            drawer.SetColour(color)
            drawer.DrawArrow(start, end)

        # Add legend if requested
        if legend and legend_entries:
            drawer.DrawText("Pharmacophore Features:", (10, 20))
            y = 40
            for feature_type in sorted(legend_entries):
                color = self.colors.get(feature_type, (0.5, 0.5, 0.5))
                drawer.SetColour(color)
                drawer.DrawText(feature_type, (20, y))
                y += 20

        # Add tooltips if interactive
        if interactive and atom_tooltips:
            drawer.AddMoleculeMetadata(mol)
            for atom_idx, tooltip in atom_tooltips.items():
                drawer.AddMetadata("atom_tooltip", atom_idx, tooltip)

        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        # Add JavaScript for interactivity if requested
        if interactive:
            svg = self._add_svg_interactivity(svg)

        return svg

    def draw_2d_alignment(self, mol1: Chem.Mol, mol2: Chem.Mol, size: Tuple[int, int] = (800, 400), show_indices: bool = False, show_shape_overlap: bool = True) -> str:
        """Draw two molecules side by side for alignment comparison.

        Args:
            mol1: First molecule
            mol2: Second molecule
            size: Image size in pixels
            show_indices: Whether to show atom indices
            show_shape_overlap: Whether to show shape overlap score

        Returns:
            SVG string of the drawing
        """
        # Generate 2D coordinates if needed
        if not mol1.GetNumConformers():
            rdDepictor.Compute2DCoords(mol1)
        if not mol2.GetNumConformers():
            rdDepictor.Compute2DCoords(mol2)

        # Create panel with two molecules
        drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
        drawer.drawOptions().addAtomIndices = show_indices
        drawer.drawOptions().prepareMolsBeforeDrawing = False

        # Draw molecules side by side
        panel_width = size[0] // 2
        drawer.SetOffset(0, 0)
        drawer.SetScale(panel_width, size[1])
        drawer.DrawMolecule(mol1, legend="Molecule 1")

        drawer.SetOffset(panel_width, 0)
        drawer.SetScale(panel_width, size[1])
        drawer.DrawMolecule(mol2, legend="Molecule 2")

        # Add shape overlap if requested
        if show_shape_overlap and mol1.GetNumConformers() and mol2.GetNumConformers():
            shape_overlap = rdShapeHelpers.ShapeTanimotoDist(mol1, mol2)
            volume_overlap = rdShapeHelpers.ShapeProtrudeDist(mol1, mol2)

            drawer.DrawText(f"Shape Overlap: {shape_overlap:.2f}\nVolume Overlap: {volume_overlap:.2f}", (10, size[1] - 30))

        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    def draw_2d_pharmacophore_match(
        self,
        mol: Chem.Mol,
        match_atoms: List[Tuple[int, ...]],
        match_types: List[str],
        size: Tuple[int, int] = (400, 400),
        show_indices: bool = False,
        highlight_style: str = "circles",
    ) -> str:
        """Draw molecule with pharmacophore matches highlighted.

        Args:
            mol: Molecule to visualize
            match_atoms: List of atom index tuples for each match
            match_types: List of feature types for each match
            size: Image size in pixels
            show_indices: Whether to show atom indices
            highlight_style: Style for highlighting ("circles", "radial", or "filled")

        Returns:
            SVG string of the drawing
        """
        # Generate 2D coordinates if needed
        if not mol.GetNumConformers():
            rdDepictor.Compute2DCoords(mol)

        # Create drawing object
        drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
        drawer.drawOptions().addStereoAnnotation = True
        drawer.drawOptions().addAtomIndices = show_indices

        # Prepare atom highlights
        atom_colors = {}
        atom_radii = {}
        for atoms, feature_type in zip(match_atoms, match_types):
            color = self.colors.get(feature_type, (0.5, 0.5, 0.5))
            for atom_idx in atoms:
                atom_colors[atom_idx] = color
                atom_radii[atom_idx] = 0.5

        # Set highlight style
        if highlight_style == "radial":
            drawer.drawOptions().fillHighlights = True
            for atom_idx, color in atom_colors.items():
                drawer.SetFillRadialGradient(atom_idx, color, (1, 1, 1), 0.3)
        elif highlight_style == "filled":
            drawer.drawOptions().fillHighlights = True

        # Draw molecule with highlights
        drawer.DrawMolecule(mol, highlightAtoms=list(atom_colors.keys()), highlightAtomColors=atom_colors, highlightAtomRadii=atom_radii)

        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    def draw_3d_conformer(
        self,
        mol: Chem.Mol,
        features: Optional[List[PharmacophoreFeature]] = None,
        conf_id: int = -1,
        size: Tuple[int, int] = (400, 400),
        surface_opacity: float = 0.3,
        show_vectors: bool = True,
        show_surface: bool = True,
        highlight_style: str = "circles",
    ) -> str:
        """Draw molecule conformer in 3D with enhanced visualization.

        Args:
            mol: Molecule to visualize
            features: Optional list of pharmacophore features to highlight
            conf_id: Conformer ID to visualize
            size: Image size in pixels
            surface_opacity: Opacity of molecular surface (0-1)
            show_vectors: Whether to show feature vectors
            show_surface: Whether to show molecular surface
            highlight_style: Style for highlighting ("circles", "radial", or "filled")

        Returns:
            SVG string of the drawing
        """
        if not mol.GetNumConformers():
            return "No 3D conformer available"

        # Create drawing object
        drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
        drawer.drawOptions().addStereoAnnotation = True
        drawer.drawOptions().addAtomIndices = False
        drawer.drawOptions().clearBackground = False
        drawer.drawOptions().continuousHighlight = True

        # Calculate and add molecular surface if requested
        if show_surface:
            surface = rdShapeHelpers.GenerateShape(mol, confId=conf_id)
            drawer.AddMoleculeMetadata(surface, "surface")
            drawer.SetSurfaceOpacity(surface_opacity)

        # Prepare atom highlights and vectors
        atom_colors = {}
        atom_radii = {}
        vectors = []

        if features:
            for feature in features:
                if not feature.enabled:
                    continue

                color = self.colors.get(feature.feature_type, (0.5, 0.5, 0.5))
                for atom_idx in feature.atoms:
                    atom_colors[atom_idx] = color
                    atom_radii[atom_idx] = feature.radius

                # Add feature vector if available
                if show_vectors and feature.vector and feature.position:
                    vectors.append((Point3D(*feature.position), Point3D(*(p + v * 2.0 for p, v in zip(feature.position, feature.vector))), color))

        # Set highlight style
        if highlight_style == "radial":
            drawer.drawOptions().fillHighlights = True
            for atom_idx, color in atom_colors.items():
                drawer.SetFillRadialGradient(atom_idx, color, (1, 1, 1), 0.3)
        elif highlight_style == "filled":
            drawer.drawOptions().fillHighlights = True

        # Draw molecule with highlights
        drawer.DrawMolecule(mol, confId=conf_id, highlightAtoms=list(atom_colors.keys()), highlightAtomColors=atom_colors, highlightAtomRadii=atom_radii)

        # Draw feature vectors
        for start, end, color in vectors:
            drawer.SetColour(color)
            drawer.DrawArrow(start, end)

        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    def draw_3d_pharmacophore(self, features: List[PharmacophoreFeature], size: Tuple[int, int] = (400, 400), show_vectors: bool = True, highlight_style: str = "circles") -> str:
        """Draw pharmacophore features in 3D space.

        Args:
            features: List of pharmacophore features to visualize
            size: Image size in pixels
            show_vectors: Whether to show feature vectors
            highlight_style: Style for highlighting ("circles", "radial", or "filled")

        Returns:
            SVG string of the drawing
        """
        # Create empty molecule to hold pharmacophore points
        mol = Chem.MolFromSmiles("")
        conf = Chem.Conformer(len(features))
        mol.AddConformer(conf)

        # Add dummy atoms at feature positions
        atom_colors = {}
        atom_radii = {}
        vectors = []

        for i, feature in enumerate(features):
            if not feature.enabled or not feature.position:
                continue

            # Add dummy atom
            atom = Chem.Atom(0)  # Dummy atom
            mol.AddAtom(atom)

            # Set 3D position
            pos = conf.GetAtomPosition(i)
            pos.x = feature.position[0]
            pos.y = feature.position[1]
            pos.z = feature.position[2]

            # Set visualization properties
            color = self.colors.get(feature.feature_type, (0.5, 0.5, 0.5))
            atom_colors[i] = color
            atom_radii[i] = feature.radius

            # Add feature vector if available
            if show_vectors and feature.vector:
                vectors.append((Point3D(*feature.position), Point3D(*(p + v * 2.0 for p, v in zip(feature.position, feature.vector))), color))

        # Create drawing object
        drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
        drawer.drawOptions().addStereoAnnotation = False
        drawer.drawOptions().addAtomIndices = False
        drawer.drawOptions().clearBackground = False

        # Set highlight style
        if highlight_style == "radial":
            drawer.drawOptions().fillHighlights = True
            for atom_idx, color in atom_colors.items():
                drawer.SetFillRadialGradient(atom_idx, color, (1, 1, 1), 0.3)
        elif highlight_style == "filled":
            drawer.drawOptions().fillHighlights = True

        # Draw pharmacophore points
        drawer.DrawMolecule(mol, highlightAtoms=list(atom_colors.keys()), highlightAtomColors=atom_colors, highlightAtomRadii=atom_radii)

        # Draw feature vectors
        for start, end, color in vectors:
            drawer.SetColour(color)
            drawer.DrawArrow(start, end)

        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    def _add_svg_interactivity(self, svg: str) -> str:
        """Add JavaScript for interactive SVG features."""
        js = """
        <script type="text/javascript">
        function showTooltip(evt, text) {
            let tooltip = document.getElementById("tooltip");
            if (!tooltip) {
                tooltip = document.createElementNS("http://www.w3.org/2000/svg", "title");
                tooltip.setAttribute("id", "tooltip");
                document.querySelector("svg").appendChild(tooltip);
            }
            tooltip.textContent = text;
            evt.target.appendChild(tooltip);
        }

        function hideTooltip(evt) {
            let tooltip = document.getElementById("tooltip");
            if (tooltip) {
                tooltip.remove();
            }
        }

        function atomClick(evt, atomIdx) {
            console.log("Clicked atom:", atomIdx);
            // Add custom click handling here
        }
        </script>
        """
        return svg.replace("</svg>", f"{js}</svg>")
