"""Enhanced compound visualization component.

This module provides an enhanced web interface for visualizing compound data with:
- Interactive structure viewer
- Binding profile plots
- Activity plots
- Safety plots
- Social data visualizations
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
from datetime import datetime

import pandas as pd
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from flask import render_template, request, jsonify

from ...models.compound import Compound
from ..base import BaseComponent, ViewResult


class CompoundVisualizationEnhanced(BaseComponent):
    """Enhanced compound visualization component."""

    def __init__(
        self,
        template_dir: Optional[Path] = None,
        static_dir: Optional[Path] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize visualization component.
        
        Args:
            template_dir: Optional template directory
            static_dir: Optional static files directory
            logger: Optional logger instance
        """
        super().__init__(template_dir, static_dir, logger)

        # Initialize tracking
        self.visualization_stats = {
            "total_visualizations": 0,
            "structure_views": 0,
            "binding_plots": 0,
            "activity_plots": 0,
            "safety_plots": 0,
            "social_plots": 0,
            "visualization_history": [],
        }

    def render_visualization(
        self,
        compound: Compound,
        plot_types: Optional[List[str]] = None,
        interactive: bool = True,
    ) -> ViewResult:
        """Render visualization interface.
        
        Args:
            compound: Compound to visualize
            plot_types: Optional list of plot types to show
            interactive: Whether to make plots interactive
            
        Returns:
            ViewResult containing rendered HTML and metadata
        """
        try:
            # Generate visualizations
            visualizations = {}

            # Structure visualization
            structure_svg = self._render_structure(compound)
            visualizations["structure"] = structure_svg
            self.visualization_stats["structure_views"] += 1

            # Plot types
            if not plot_types:
                plot_types = ["binding", "activity", "safety", "social"]

            # Generate requested plots
            for plot_type in plot_types:
                if plot_type == "binding" and hasattr(compound, "binding_data"):
                    plot = self._create_binding_plot(compound, interactive)
                    visualizations["binding_plot"] = plot
                    self.visualization_stats["binding_plots"] += 1

                elif plot_type == "activity" and hasattr(compound, "predictions"):
                    plot = self._create_activity_plot(compound, interactive)
                    visualizations["activity_plot"] = plot
                    self.visualization_stats["activity_plots"] += 1

                elif plot_type == "safety" and hasattr(compound, "predictions"):
                    plot = self._create_safety_plot(compound, interactive)
                    visualizations["safety_plot"] = plot
                    self.visualization_stats["safety_plots"] += 1

                elif plot_type == "social" and hasattr(compound, "social_data"):
                    plot = self._create_social_plot(compound, interactive)
                    visualizations["social_plot"] = plot
                    self.visualization_stats["social_plots"] += 1

            # Update stats
            self.visualization_stats["total_visualizations"] += 1
            self.visualization_stats["visualization_history"].append({
                "timestamp": datetime.now().isoformat(),
                "compound": compound.name,
                "plot_types": plot_types,
                "interactive": interactive,
            })

            # Render template
            html = render_template(
                "compound_visualization.html",
                compound=compound,
                visualizations=visualizations,
                plot_types=plot_types,
                interactive=interactive,
                stats=self.visualization_stats,
            )

            return ViewResult(
                success=True,
                data={
                    "html": html,
                    "visualizations": visualizations,
                    "plot_types": plot_types,
                    "visualization_stats": self.visualization_stats,
                },
            )

        except Exception as e:
            self.logger.error(f"Error rendering visualization: {str(e)}")
            return ViewResult(
                success=False,
                error=str(e),
            )

    def _render_structure(
        self,
        compound: Compound,
        width: int = 400,
        height: int = 400,
    ) -> str:
        """Render compound structure as SVG.
        
        Args:
            compound: Compound to visualize
            width: Image width
            height: Image height
            
        Returns:
            SVG string
        """
        # Parse SMILES
        mol = Chem.MolFromSmiles(compound.smiles)
        if not mol:
            raise ValueError(f"Invalid SMILES: {compound.smiles}")

        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)

        # Draw structure
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        return drawer.GetDrawingText()

    def _create_binding_plot(
        self,
        compound: Compound,
        interactive: bool = True,
    ) -> Dict[str, Any]:
        """Create binding profile plot.
        
        Args:
            compound: Compound to visualize
            interactive: Whether to make plot interactive
            
        Returns:
            Plotly figure dict
        """
        # Extract data
        targets = []
        affinities = []
        confidences = []

        for binding in compound.binding_data:
            targets.append(binding["target"])
            affinities.append(float(binding["affinity"]))
            confidences.append(binding.get("confidence", 1.0))

        # Create plot
        fig = go.Figure()

        # Add bars
        fig.add_trace(
            go.Bar(
                x=targets,
                y=affinities,
                error_y=dict(
                    type="data",
                    array=[1 - c for c in confidences],
                    visible=True,
                ),
                name="Binding Affinity",
            )
        )

        # Update layout
        fig.update_layout(
            title="Receptor Binding Profile",
            xaxis_title="Target",
            yaxis_title="Affinity (Ki)",
            showlegend=True,
        )

        return fig.to_dict() if interactive else fig.to_image(format="svg")

    def _create_activity_plot(
        self,
        compound: Compound,
        interactive: bool = True,
    ) -> Dict[str, Any]:
        """Create activity profile plot.
        
        Args:
            compound: Compound to visualize
            interactive: Whether to make plot interactive
            
        Returns:
            Plotly figure dict
        """
        # Extract data
        activities = compound.predictions.get("activity", {})
        labels = list(activities.keys())
        values = list(activities.values())

        # Create plot
        fig = go.Figure()

        # Add radar plot
        fig.add_trace(
            go.Scatterpolar(
                r=values,
                theta=labels,
                fill="toself",
                name="Activity Profile",
            )
        )

        # Update layout
        fig.update_layout(
            title="Activity Profile",
            polar=dict(
                radialaxis=dict(
                    visible=True,
                    range=[0, 1],
                ),
            ),
            showlegend=True,
        )

        return fig.to_dict() if interactive else fig.to_image(format="svg")

    def _create_safety_plot(
        self,
        compound: Compound,
        interactive: bool = True,
    ) -> Dict[str, Any]:
        """Create safety profile plot.
        
        Args:
            compound: Compound to visualize
            interactive: Whether to make plot interactive
            
        Returns:
            Plotly figure dict
        """
        # Extract data
        safety = compound.predictions.get("safety", {})
        risks = list(safety.keys())
        scores = list(safety.values())

        # Create plot
        fig = go.Figure()

        # Add heatmap
        fig.add_trace(
            go.Heatmap(
                z=[scores],
                x=risks,
                y=["Risk Level"],
                colorscale="RdYlGn_r",
                showscale=True,
            )
        )

        # Update layout
        fig.update_layout(
            title="Safety Profile",
            xaxis_title="Risk Category",
            yaxis_title="",
            showlegend=False,
        )

        return fig.to_dict() if interactive else fig.to_image(format="svg")

    def _create_social_plot(
        self,
        compound: Compound,
        interactive: bool = True,
    ) -> Dict[str, Any]:
        """Create social data visualization.
        
        Args:
            compound: Compound to visualize
            interactive: Whether to make plot interactive
            
        Returns:
            Plotly figure dict
        """
        # Extract data
        social_data = compound.social_data
        dates = []
        counts = []

        # Reddit posts
        if "reddit" in social_data:
            for post in social_data["reddit"].get("posts", []):
                dates.append(datetime.fromisoformat(post["created_utc"]))
                counts.append(1)

        # Twitter mentions
        if "twitter" in social_data:
            for tweet in social_data["twitter"].get("tweets", []):
                dates.append(datetime.fromisoformat(tweet["created_at"]))
                counts.append(1)

        # Create plot
        fig = go.Figure()

        # Add timeline
        fig.add_trace(
            go.Scatter(
                x=dates,
                y=counts,
                mode="markers",
                name="Social Mentions",
            )
        )

        # Update layout
        fig.update_layout(
            title="Social Media Timeline",
            xaxis_title="Date",
            yaxis_title="Mentions",
            showlegend=True,
        )

        return fig.to_dict() if interactive else fig.to_image(format="svg")

    def get_metrics(self) -> Dict[str, Any]:
        """Get component metrics."""
        return {
            "visualization_stats": self.visualization_stats,
        }
