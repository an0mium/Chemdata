"""Enhanced compound export component.

This module provides an enhanced web interface for exporting compound data with:
- Multiple export formats
- Column selection
- Data validation
- Export history
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
from datetime import datetime

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from flask import render_template, request, jsonify

from ...models.compound import Compound
from ..base import BaseComponent, ViewResult


class CompoundExportEnhanced(BaseComponent):
    """Enhanced compound export component."""

    def __init__(
        self,
        template_dir: Optional[Path] = None,
        static_dir: Optional[Path] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize export component.
        
        Args:
            template_dir: Optional template directory
            static_dir: Optional static files directory
            logger: Optional logger instance
        """
        super().__init__(template_dir, static_dir, logger)

        # Initialize tracking
        self.export_stats = {
            "total_exports": 0,
            "tsv_exports": 0,
            "sdf_exports": 0,
            "json_exports": 0,
            "export_history": [],
        }

        # Define available columns
        self.available_columns = {
            "basic": [
                "name",
                "smiles",
                "cas_number",
            ],
            "binding": [
                "binding_data",
                "binding_targets",
                "binding_affinities",
                "binding_confidences",
            ],
            "social": [
                "social_data",
                "reddit_posts",
                "twitter_mentions",
                "post_dates",
                "post_counts",
            ],
            "predictions": [
                "binding_predictions",
                "activity_predictions",
                "safety_predictions",
                "prediction_confidences",
            ],
        }

    def render_export(
        self,
        compounds: List[Compound],
        format: str = "tsv",
        columns: Optional[List[str]] = None,
        validate: bool = True,
    ) -> ViewResult:
        """Render export interface.
        
        Args:
            compounds: List of compounds to export
            format: Export format (tsv/sdf/json)
            columns: Optional list of columns to export
            validate: Whether to validate data before export
            
        Returns:
            ViewResult containing rendered HTML and metadata
        """
        try:
            # Validate data
            if validate:
                validation_errors = self._validate_data(compounds, columns)
                if validation_errors:
                    return ViewResult(
                        success=False,
                        error="Validation errors",
                        data={"validation_errors": validation_errors},
                    )

            # Export data
            if format == "tsv":
                export_data = self._export_tsv(compounds, columns)
                self.export_stats["tsv_exports"] += 1
            elif format == "sdf":
                export_data = self._export_sdf(compounds)
                self.export_stats["sdf_exports"] += 1
            elif format == "json":
                export_data = self._export_json(compounds, columns)
                self.export_stats["json_exports"] += 1
            else:
                raise ValueError(f"Unsupported format: {format}")

            # Update stats
            self.export_stats["total_exports"] += 1
            self.export_stats["export_history"].append({
                "timestamp": datetime.now().isoformat(),
                "format": format,
                "columns": columns,
                "compound_count": len(compounds),
            })

            # Render template
            html = render_template(
                "compound_export.html",
                compounds=compounds,
                format=format,
                columns=columns or [],
                available_columns=self.available_columns,
                stats=self.export_stats,
            )

            return ViewResult(
                success=True,
                data={
                    "html": html,
                    "export_data": export_data,
                    "format": format,
                    "columns": columns,
                    "export_stats": self.export_stats,
                },
            )

        except Exception as e:
            self.logger.error(f"Error rendering export: {str(e)}")
            return ViewResult(
                success=False,
                error=str(e),
            )

    def _validate_data(
        self,
        compounds: List[Compound],
        columns: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """Validate compound data before export.
        
        Args:
            compounds: List of compounds to validate
            columns: Optional list of columns to validate
            
        Returns:
            List of validation errors
        """
        errors = []

        for compound in compounds:
            # Validate basic fields
            if not compound.name:
                errors.append({
                    "compound": compound.cas_number,
                    "field": "name",
                    "error": "Missing name",
                })
            if not compound.smiles:
                errors.append({
                    "compound": compound.name,
                    "field": "smiles",
                    "error": "Missing SMILES",
                })
            if not compound.cas_number:
                errors.append({
                    "compound": compound.name,
                    "field": "cas_number",
                    "error": "Missing CAS number",
                })

            # Validate structure
            mol = Chem.MolFromSmiles(compound.smiles)
            if not mol:
                errors.append({
                    "compound": compound.name,
                    "field": "smiles",
                    "error": "Invalid SMILES",
                })

            # Validate requested columns
            if columns:
                for column in columns:
                    if column in self.available_columns["binding"]:
                        if not hasattr(compound, "binding_data"):
                            errors.append({
                                "compound": compound.name,
                                "field": column,
                                "error": "Missing binding data",
                            })
                    elif column in self.available_columns["social"]:
                        if not hasattr(compound, "social_data"):
                            errors.append({
                                "compound": compound.name,
                                "field": column,
                                "error": "Missing social data",
                            })
                    elif column in self.available_columns["predictions"]:
                        if not hasattr(compound, "predictions"):
                            errors.append({
                                "compound": compound.name,
                                "field": column,
                                "error": "Missing predictions",
                            })

        return errors

    def _export_tsv(
        self,
        compounds: List[Compound],
        columns: Optional[List[str]] = None,
    ) -> str:
        """Export compounds as TSV.
        
        Args:
            compounds: List of compounds to export
            columns: Optional list of columns to export
            
        Returns:
            TSV string
        """
        # Build data
        data = []
        for compound in compounds:
            row = {
                "name": compound.name,
                "smiles": compound.smiles,
                "cas_number": compound.cas_number,
            }

            # Add binding data
            if hasattr(compound, "binding_data"):
                row["binding_data"] = compound.binding_data
                row["binding_targets"] = [b["target"] for b in compound.binding_data]
                row["binding_affinities"] = [b["affinity"] for b in compound.binding_data]
                row["binding_confidences"] = [b.get("confidence", 1.0) for b in compound.binding_data]

            # Add social data
            if hasattr(compound, "social_data"):
                row["social_data"] = compound.social_data
                row["reddit_posts"] = len(compound.social_data.get("reddit", {}).get("posts", []))
                row["twitter_mentions"] = len(compound.social_data.get("twitter", {}).get("tweets", []))

            # Add predictions
            if hasattr(compound, "predictions"):
                row["binding_predictions"] = compound.predictions.get("binding", {})
                row["activity_predictions"] = compound.predictions.get("activity", {})
                row["safety_predictions"] = compound.predictions.get("safety", {})

            data.append(row)

        # Convert to DataFrame
        df = pd.DataFrame(data)

        # Filter columns
        if columns:
            df = df[columns]

        # Export
        return df.to_csv(sep="\t", index=False)

    def _export_sdf(
        self,
        compounds: List[Compound],
    ) -> str:
        """Export compounds as SDF.
        
        Args:
            compounds: List of compounds to export
            
        Returns:
            SDF string
        """
        sdf_writer = Chem.SDWriter("temp.sdf")

        for compound in compounds:
            mol = Chem.MolFromSmiles(compound.smiles)
            if not mol:
                continue

            # Add properties
            mol.SetProp("_Name", compound.name)
            mol.SetProp("CAS", compound.cas_number)

            # Add binding data
            if hasattr(compound, "binding_data"):
                mol.SetProp("BindingData", str(compound.binding_data))

            # Add social data
            if hasattr(compound, "social_data"):
                mol.SetProp("SocialData", str(compound.social_data))

            # Add predictions
            if hasattr(compound, "predictions"):
                mol.SetProp("Predictions", str(compound.predictions))

            sdf_writer.write(mol)

        sdf_writer.close()
        with open("temp.sdf", "r") as f:
            sdf_data = f.read()

        return sdf_data

    def _export_json(
        self,
        compounds: List[Compound],
        columns: Optional[List[str]] = None,
    ) -> str:
        """Export compounds as JSON.
        
        Args:
            compounds: List of compounds to export
            columns: Optional list of columns to export
            
        Returns:
            JSON string
        """
        # Build data
        data = []
        for compound in compounds:
            row = {
                "name": compound.name,
                "smiles": compound.smiles,
                "cas_number": compound.cas_number,
            }

            # Add binding data
            if hasattr(compound, "binding_data"):
                row["binding_data"] = compound.binding_data

            # Add social data
            if hasattr(compound, "social_data"):
                row["social_data"] = compound.social_data

            # Add predictions
            if hasattr(compound, "predictions"):
                row["predictions"] = compound.predictions

            # Filter columns
            if columns:
                row = {k: v for k, v in row.items() if k in columns}

            data.append(row)

        # Convert to JSON
        return pd.Series(data).to_json()

    def get_metrics(self) -> Dict[str, Any]:
        """Get component metrics."""
        return {
            "export_stats": self.export_stats,
        }
