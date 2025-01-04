"""Base molecular structure depiction functionality.

This module provides base functionality for:
1. Common drawing options and configuration
2. Image format handling and saving
3. Error handling and logging
4. Shared utility functions
"""

import io
import logging
from typing import Dict, List, Optional, Tuple, Union
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw


class BaseStructureDepiction:
    """Base class for molecular structure depiction."""

    # Common drawing options
    DEFAULT_DRAWING_OPTIONS = {
        "size": (400, 400),
        "background_color": (1.0, 1.0, 1.0),
        "atom_colors": None,  # Use RDKit defaults
        "bond_line_width": 2,
        "atom_label_font_size": 12,
        "add_hydrogens": True,
        "random_seed": 42,
    }

    # Feature color schemes
    FEATURE_COLORS = {
        "donor": (0.0, 1.0, 0.0),  # Green
        "acceptor": (1.0, 0.0, 0.0),  # Red
        "aromatic": (0.0, 0.0, 1.0),  # Blue
        "hydrophobic": (1.0, 1.0, 0.0),  # Yellow
        "positive": (1.0, 0.5, 0.0),  # Orange
        "negative": (0.5, 0.0, 1.0),  # Purple
        "other": (0.5, 0.5, 0.5),  # Gray
    }

    def __init__(self, config: Optional[Dict] = None):
        """Initialize base structure depiction.

        Args:
            config: Optional configuration dictionary
        """
        self.logger = logging.getLogger(__name__)
        self.config = config or {}
        self.drawing_options = self.DEFAULT_DRAWING_OPTIONS.copy()
        self.drawing_options.update(self.config.get("drawing_options", {}))

    def set_drawing_options(self, **kwargs) -> None:
        """Set drawing options.

        Args:
            **kwargs: Drawing option key-value pairs
        """
        self.drawing_options.update(kwargs)

    def get_drawing_size(self) -> Tuple[int, int]:
        """Get current drawing size.

        Returns:
            Tuple of (width, height)
        """
        return self.drawing_options["size"]

    def validate_mol(self, mol: Optional[Chem.Mol]) -> bool:
        """Validate RDKit molecule object.

        Args:
            mol: RDKit molecule to validate

        Returns:
            True if valid, False otherwise
        """
        if mol is None:
            self.logger.warning("Molecule is None")
            return False
        try:
            Chem.SanitizeMol(mol)
            return True
        except Exception as e:
            self.logger.warning(f"Invalid molecule: {str(e)}")
            return False

    def add_hydrogens(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Add hydrogens to molecule if configured.

        Args:
            mol: Input molecule

        Returns:
            Molecule with hydrogens added
        """
        try:
            if mol is None:
                return None
            if self.drawing_options["add_hydrogens"]:
                return Chem.AddHs(mol)
            return mol
        except Exception as e:
            self.logger.error(f"Error adding hydrogens: {str(e)}")
            return None

    def create_image(
        self, img_data: Union[bytes, Image.Image], return_pil: bool
    ) -> Union[bytes, Image.Image, None]:
        """Create image from raw data.

        Args:
            img_data: Raw image data
            return_pil: Whether to return PIL Image

        Returns:
            Image data in requested format
        """
        try:
            if img_data is None:
                return None

            if return_pil and not isinstance(img_data, Image.Image):
                return Image.open(io.BytesIO(img_data))
            elif not return_pil and isinstance(img_data, Image.Image):
                img_buffer = io.BytesIO()
                img_data.save(img_buffer, format="PNG")
                return img_buffer.getvalue()
            return img_data

        except Exception as e:
            self.logger.error(f"Error creating image: {str(e)}")
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
            img_format: Image format (default: png)

        Returns:
            True if successful, False otherwise
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

    def get_feature_color(self, feature_type: str) -> Tuple[float, float, float]:
        """Get color for feature type.

        Args:
            feature_type: Type of structural/pharmacophore feature

        Returns:
            RGB color tuple
        """
        return self.FEATURE_COLORS.get(feature_type, self.FEATURE_COLORS["other"])

    def create_legend(self, text: str, font_size: Optional[int] = None) -> str:
        """Create legend text.

        Args:
            text: Legend text
            font_size: Optional font size override

        Returns:
            Formatted legend text
        """
        if font_size is None:
            font_size = self.drawing_options["atom_label_font_size"]
        return f'<font size="{font_size}">{text}</font>'

    def get_drawing_options_dict(self) -> Dict:
        """Get drawing options as dictionary.

        Returns:
            Dictionary of current drawing options
        """
        return {
            "bondLineWidth": self.drawing_options["bond_line_width"],
            "atomLabelFontSize": self.drawing_options["atom_label_font_size"],
            "backgroundColor": self.drawing_options["background_color"],
            "atomColors": self.drawing_options["atom_colors"],
        }
