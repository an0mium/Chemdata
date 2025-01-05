"""Enhanced compound detail view component.

This module provides an enhanced web interface for displaying detailed compound information with:
- Interactive structure visualization
- Binding affinity plots
- Social media analysis
- Export capabilities
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
from datetime import datetime

import pandas as pd
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import Draw
from flask import render_template, request, jsonify

from ...models.compound import Compound
from ..base import BaseComponent, ViewResult


class CompoundDetailEnhanced(BaseComponent):
    """Enhanced compound detail view component."""

    def __init__(
        self,
        template_dir: Optional[Path] = None,
        static_dir: Optional[Path] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize compound detail component.
        
        Args:
            template_dir: Optional template directory
            static_dir: Optional static files directory
            logger: Optional logger instance
        """
        super().__init__(template_dir, static_dir, logger)

        # Initialize tracking
        self.view_stats = {
            "total_views": 0,
            "structure_views": 0,
            "binding_views": 0,
            "social_views": 0,
            "exported_views": 0,
        }

    def render_detail(
        self,
        compound: Compound,
        include_structure: bool = True,
        include_binding: bool = True,
        include_social: bool = True,
    ) -> ViewResult:
        """Render compound detail view.
        
        Args:
            compound: Compound to display
            include_structure: Whether to include structure visualization
            include_binding: Whether to include binding data
            include_social: Whether to include social data
            
        Returns:
            ViewResult containing rendered HTML and metadata
        """
        try:
            # Generate visualizations
            visualizations = {}
            if include_structure:
                structure_svg = self._generate_structure(compound)
                visualizations["structure"] = structure_svg
                self.view_stats["structure_views"] += 1

            if include_binding and hasattr(compound, "binding_data"):
                binding_plot = self._generate_binding_plot(compound)
                visualizations["binding"] = binding_plot
                self.view_stats["binding_views"] += 1

            if include_social and hasattr(compound, "social_data"):
                social_plots = self._generate_social_plots(compound)
                visualizations["social"] = social_plots
                self.view_stats["social_views"] += 1

            # Update stats
            self.view_stats["total_views"] += 1

            # Render template
            html = render_template(
                "compound_detail.html",
                compound=compound,
                visualizations=visualizations,
                stats=self.view_stats,
            )

            return ViewResult(
                success=True,
                data={
                    "html": html,
                    "compound": compound,
                    "visualizations": visualizations,
                },
            )

        except Exception as e:
            self.logger.error(f"Error rendering compound detail: {str(e)}")
            return ViewResult(
                success=False,
                error=str(e),
            )

    def handle_export(
        self,
        compound: Compound,
        format: str = "json",
        include_structure: bool = True,
        include_binding: bool = True,
        include_social: bool = True,
    ) -> ViewResult:
        """Handle export request.
        
        Args:
            compound: Compound to export
            format: Export format (json/sdf)
            include_structure: Whether to include structure
            include_binding: Whether to include binding data
            include_social: Whether to include social data
            
        Returns:
            ViewResult containing export data
        """
        try:
            # Build export data
            data = {
                "name": compound.name,
                "smiles": compound.smiles,
                "cas_number": compound.cas_number,
            }

            if include_binding and hasattr(compound, "binding_data"):
                data["binding_data"] = compound.binding_data

            if include_social and hasattr(compound, "social_data"):
                data["social_data"] = compound.social_data

            # Export
            if format == "json":
                export_data = pd.Series(data).to_json()
            elif format == "sdf":
                mol = Chem.MolFromSmiles(compound.smiles)
                for key, value in data.items():
                    if key != "smiles":
                        mol.SetProp(key, str(value))
                export_data = Chem.SDWriter(mol)
            else:
                raise ValueError(f"Unsupported format: {format}")

            # Update stats
            self.view_stats["exported_views"] += 1

            return ViewResult(
                success=True,
                data={
                    "export_data": export_data,
                    "format": format,
                },
            )

        except Exception as e:
            self.logger.error(f"Error exporting compound: {str(e)}")
            return ViewResult(
                success=False,
                error=str(e),
            )

    def _generate_structure(
        self,
        compound: Compound,
    ) -> str:
        """Generate structure visualization.
        
        Args:
            compound: Compound to visualize
            
        Returns:
            SVG string of structure visualization
        """
        mol = Chem.MolFromSmiles(compound.smiles)
        return Draw.MolToImage(mol).tostring()

    def _generate_binding_plot(
        self,
        compound: Compound,
    ) -> str:
        """Generate binding affinity plot.
        
        Args:
            compound: Compound to visualize
            
        Returns:
            JSON string of binding plot
        """
        if not hasattr(compound, "binding_data"):
            return None

        # Extract binding data
        data = []
        for binding in compound.binding_data:
            data.append({
                "target": binding["target"],
                "affinity": float(binding["affinity"]),
                "confidence": binding.get("confidence", 1.0),
            })

        # Create plot
        df = pd.DataFrame(data)
        fig = go.Figure(data=[
            go.Bar(
                x=df["target"],
                y=df["affinity"],
                error_y=dict(
                    type="data",
                    array=1 - df["confidence"],
                    visible=True,
                ),
                name="Binding Affinity",
            )
        ])
        fig.update_layout(
            title="Binding Affinity by Target",
            yaxis_title="Affinity (nM)",
            showlegend=True,
        )
        return fig.to_json()

    def _generate_social_plots(
        self,
        compound: Compound,
    ) -> Dict[str, str]:
        """Generate social media analysis plots.
        
        Args:
            compound: Compound to visualize
            
        Returns:
            Dictionary of plot names to JSON strings
        """
        if not hasattr(compound, "social_data"):
            return None

        plots = {}

        # Reddit activity
        if "reddit" in compound.social_data:
            reddit_data = compound.social_data["reddit"]
            posts_by_date = pd.DataFrame(reddit_data["posts"])
            if not posts_by_date.empty:
                posts_by_date["date"] = pd.to_datetime(posts_by_date["created_utc"], unit="s")
                posts_by_date = posts_by_date.groupby("date").size()
                
                fig = go.Figure(data=[
                    go.Scatter(
                        x=posts_by_date.index,
                        y=posts_by_date.values,
                        mode="lines+markers",
                        name="Reddit Posts",
                    )
                ])
                fig.update_layout(
                    title="Reddit Activity Over Time",
                    yaxis_title="Number of Posts",
                    showlegend=True,
                )
                plots["reddit_activity"] = fig.to_json()

        # Twitter activity
        if "twitter" in compound.social_data:
            twitter_data = compound.social_data["twitter"]
            tweets_by_date = pd.DataFrame(twitter_data["tweets"])
            if not tweets_by_date.empty:
                tweets_by_date["date"] = pd.to_datetime(tweets_by_date["created_at"])
                tweets_by_date = tweets_by_date.groupby("date").size()
                
                fig = go.Figure(data=[
                    go.Scatter(
                        x=tweets_by_date.index,
                        y=tweets_by_date.values,
                        mode="lines+markers",
                        name="Twitter Mentions",
                    )
                ])
                fig.update_layout(
                    title="Twitter Activity Over Time",
                    yaxis_title="Number of Mentions",
                    showlegend=True,
                )
                plots["twitter_activity"] = fig.to_json()

        return plots

    def get_metrics(self) -> Dict[str, Any]:
        """Get component metrics."""
        return {
            "view_stats": self.view_stats,
        }
