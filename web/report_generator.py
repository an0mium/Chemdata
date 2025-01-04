"""Report generation functionality for ChemData application."""

import base64
from datetime import datetime
from typing import Optional

from jinja2 import Environment, FileSystemLoader
from rdkit import Chem
from rdkit.Chem import Draw

from binding_data_processor.models.compound import CompoundData


class ReportGenerator:
    """Generator for detailed compound reports."""

    def __init__(self, template_dir: str = "templates"):
        """Initialize report generator.

        Args:
            template_dir: Directory containing report templates
        """
        self.env = Environment(
            loader=FileSystemLoader(template_dir),
            trim_blocks=True,
            lstrip_blocks=True,
        )
        self.template = self.env.get_template("report.html")

        # Add custom filters
        self.env.filters["zip"] = zip
        self.env.filters["slice"] = lambda value, length: value[:length]

    def generate_report(
        self,
        compound: CompoundData,
        data_version: Optional[str] = None,
        model_version: Optional[str] = None,
    ) -> str:
        """Generate detailed HTML report for compound.

        Args:
            compound: Compound data to generate report for
            data_version: Optional version of data sources
            model_version: Optional version of ML models

        Returns:
            HTML report as string
        """
        # Generate structure image
        structure_image = self._get_structure_image(compound.smiles)

        # Prepare template variables
        template_vars = {
            "compound": compound,
            "structure_image": structure_image,
            "generation_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "data_version": data_version or "Unknown",
            "model_version": model_version or "Unknown",
        }

        # Render template
        return self.template.render(**template_vars)

    def _get_structure_image(self, smiles: str) -> Optional[str]:
        """Convert SMILES to base64 encoded PNG image.

        Args:
            smiles: SMILES string to convert

        Returns:
            Base64 encoded PNG image or None if conversion fails
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Generate 2D depiction
                img = Draw.MolToImage(
                    mol,
                    size=(600, 600),  # Higher resolution for reports
                    imageType="PNG",
                    fitImage=True,
                )

                # Convert to base64
                import io

                img_buffer = io.BytesIO()
                img.save(img_buffer, format="PNG")
                img_str = base64.b64encode(img_buffer.getvalue()).decode()

                return img_str
        except Exception as e:
            print(f"Error generating structure image: {e}")
            return None
