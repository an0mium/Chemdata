"""Enhanced visualization module for abuse potential predictions.

This module provides comprehensive visualization capabilities for:
1. Structural highlighting of abuse-related features
2. Interactive plots of receptor interactions 
3. Risk level visualizations and heatmaps
4. Mechanism contribution networks
5. Uncertainty visualization
6. Web data integration displays
7. 3D structure visualization with highlighted pharmacophores
8. Machine learning model interpretability plots
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from plotly.subplots import make_subplots
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdDepictor
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

from ..predictors.abuse_types import (
    ABUSE_CATEGORIES,
    ABUSE_EFFECTS,
    RECEPTOR_SYSTEMS,
    RISK_LEVELS,
    MECHANISMS,
)
from ...pharmacophore import PharmacophoreGenerator
from ....base import VisualizerBase


class AbusePotentialVisualizer(VisualizerBase):
    """Enhanced visualization tools for abuse potential predictions."""

    # Color schemes for different visualization styles
    COLOR_SCHEMES = {
        "modern": {
            "high_risk": "#FF4B4B",
            "moderate_risk": "#FFA500",
            "low_risk": "#4CAF50",
            "receptor_highlight": "#2196F3",
            "uncertainty": "#9C27B0",
            "background": "#FFFFFF",
            "text": "#000000",
        },
        "dark": {
            "high_risk": "#FF5252",
            "moderate_risk": "#FFB74D",
            "low_risk": "#81C784",
            "receptor_highlight": "#64B5F6",
            "uncertainty": "#BA68C8",
            "background": "#212121",
            "text": "#FFFFFF",
        },
        "classic": {
            "high_risk": "#DC3545",
            "moderate_risk": "#FFC107",
            "low_risk": "#28A745",
            "receptor_highlight": "#007BFF",
            "uncertainty": "#6610F2",
            "background": "#F8F9FA",
            "text": "#212529",
        },
    }

    def __init__(
        self,
        style: str = "modern",
        pharmacophore_gen: Optional[PharmacophoreGenerator] = None,
        interactive: bool = True,
    ):
        """Initialize visualizer.

        Args:
            style: Visualization style ('modern', 'dark', or 'classic')
            pharmacophore_gen: Optional PharmacophoreGenerator instance
            interactive: Whether to use interactive Plotly plots
        """
        super().__init__()
        self.style = style
        self.colors = self.COLOR_SCHEMES[style]
        self.interactive = interactive
        self.pharmacophore_gen = pharmacophore_gen or PharmacophoreGenerator()

        # Set up plotting styles
        if not interactive:
            sns.set_style("darkgrid" if style == "dark" else "whitegrid")
            plt.style.use("dark_background" if style == "dark" else "default")

    def visualize_structure_highlights(
        self,
        mol: Chem.Mol,
        predictions: Dict,
        highlight_type: str = "risk",
        size: Tuple[int, int] = (800, 800),
        show_3d: bool = False,
    ) -> Union[Draw.Image, go.Figure]:
        """Highlight structural features related to abuse potential.

        Args:
            mol: RDKit molecule
            predictions: Prediction results
            highlight_type: Type of highlighting ('risk', 'receptor', or 'mechanism')
            size: Image size (width, height)
            show_3d: Whether to show 3D conformation

        Returns:
            RDKit image or Plotly figure for 3D
        """
        try:
            # Generate 3D conformation if needed
            if show_3d and not mol.GetConformer().Is3D():
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol)
            elif not mol.GetConformer().Is3D():
                rdDepictor.Compute2DCoords(mol)

            # Get atom contributions based on highlight type
            if highlight_type == "risk":
                atom_colors = self._get_risk_based_colors(mol, predictions)
            elif highlight_type == "receptor":
                atom_colors = self._get_receptor_based_colors(mol, predictions)
            else:  # mechanism
                atom_colors = self._get_mechanism_based_colors(mol, predictions)

            if show_3d:
                return self._create_3d_structure_plot(mol, atom_colors)
            else:
                return self._create_2d_structure_plot(mol, atom_colors, size)

        except Exception as e:
            self.logger.error(f"Error in structure visualization: {str(e)}")
            return None

    def plot_abuse_potential(
        self,
        predictions: Dict,
        plot_type: str = "radar",
        show_uncertainty: bool = True,
        show_web_data: bool = True,
    ) -> Union[Figure, go.Figure]:
        """Plot abuse potential predictions.

        Args:
            predictions: Prediction results
            plot_type: Type of plot ('radar', 'bar', or 'sunburst')
            show_uncertainty: Whether to show uncertainty ranges
            show_web_data: Whether to include web data annotations

        Returns:
            Matplotlib or Plotly figure
        """
        try:
            if "abuse_potential" not in predictions:
                raise ValueError("No abuse potential predictions found")

            # Extract data
            categories = []
            probabilities = []
            uncertainties = []
            web_data = []

            for pred in predictions["abuse_potential"]:
                categories.append(pred["category"])
                probabilities.append(pred["probability"])
                if show_uncertainty:
                    uncertainties.append(pred.get("uncertainty", 0))
                if show_web_data:
                    web_data.append(pred.get("web_data", {}))

            if self.interactive:
                return self._create_interactive_abuse_plot(
                    categories,
                    probabilities,
                    uncertainties if show_uncertainty else None,
                    web_data if show_web_data else None,
                    plot_type,
                )
            else:
                return self._create_static_abuse_plot(
                    categories,
                    probabilities,
                    uncertainties if show_uncertainty else None,
                    web_data if show_web_data else None,
                    plot_type,
                )

        except Exception as e:
            self.logger.error(f"Error plotting abuse potential: {str(e)}")
            return None

    def plot_receptor_network(
        self,
        predictions: Dict,
        layout: str = "force",
        min_affinity: float = 0.5,
        show_pharmacophores: bool = True,
    ) -> Union[Figure, go.Figure]:
        """Plot receptor interaction network.

        Args:
            predictions: Prediction results
            layout: Network layout ('force', 'circular', or 'hierarchical')
            min_affinity: Minimum affinity threshold
            show_pharmacophores: Whether to show pharmacophore features

        Returns:
            Matplotlib or Plotly figure
        """
        try:
            if "receptor_affinities" not in predictions:
                raise ValueError("No receptor affinity predictions found")

            # Create network
            G = nx.Graph()

            # Add receptor nodes
            for receptor in predictions["receptor_affinities"]:
                if receptor["affinity"] >= min_affinity:
                    G.add_node(
                        receptor["receptor"],
                        affinity=receptor["affinity"],
                        system=receptor["system"],
                        pharmacophores=(
                            receptor.get("pharmacophores", [])
                            if show_pharmacophores
                            else []
                        ),
                    )

            # Add edges between related receptors
            for r1 in G.nodes():
                for r2 in G.nodes():
                    if r1 < r2:  # Avoid duplicate edges
                        if self._are_receptors_related(r1, r2):
                            weight = (
                                G.nodes[r1]["affinity"] + G.nodes[r2]["affinity"]
                            ) / 2
                            G.add_edge(r1, r2, weight=weight)

            if self.interactive:
                return self._create_interactive_network(G, layout, show_pharmacophores)
            else:
                return self._create_static_network(G, layout, show_pharmacophores)

        except Exception as e:
            self.logger.error(f"Error plotting receptor network: {str(e)}")
            return None

    def plot_mechanism_contributions(
        self,
        predictions: Dict,
        plot_type: str = "sankey",
        min_contribution: float = 0.1,
    ) -> Union[Figure, go.Figure]:
        """Plot mechanism contributions.

        Args:
            predictions: Prediction results
            plot_type: Type of plot ('sankey', 'sunburst', or 'treemap')
            min_contribution: Minimum contribution threshold

        Returns:
            Matplotlib or Plotly figure
        """
        try:
            if "mechanisms" not in predictions:
                raise ValueError("No mechanism predictions found")

            # Extract mechanism data
            mechanisms = []
            contributions = []
            categories = []

            for category, data in predictions["mechanisms"].items():
                if data["contribution"] >= min_contribution:
                    for mechanism in data["mechanisms"]:
                        mechanisms.append(mechanism)
                        contributions.append(data["contribution"])
                        categories.append(category)

            if self.interactive:
                return self._create_interactive_mechanism_plot(
                    mechanisms,
                    contributions,
                    categories,
                    plot_type,
                )
            else:
                return self._create_static_mechanism_plot(
                    mechanisms,
                    contributions,
                    categories,
                    plot_type,
                )

        except Exception as e:
            self.logger.error(f"Error plotting mechanisms: {str(e)}")
            return None

    def plot_model_interpretability(
        self,
        predictions: Dict,
        mol: Optional[Chem.Mol] = None,
        plot_type: str = "attention",
    ) -> Union[Figure, go.Figure]:
        """Plot model interpretability visualizations.

        Args:
            predictions: Prediction results
            mol: Optional RDKit molecule for structural interpretation
            plot_type: Type of plot ('attention', 'feature_importance', or 'embedding')

        Returns:
            Matplotlib or Plotly figure
        """
        try:
            if "model_interpretation" not in predictions:
                raise ValueError("No model interpretation data found")

            interp_data = predictions["model_interpretation"]

            if plot_type == "attention":
                return self._plot_attention_weights(interp_data, mol)
            elif plot_type == "feature_importance":
                return self._plot_feature_importance(interp_data)
            else:  # embedding
                return self._plot_embedding_space(interp_data)

        except Exception as e:
            self.logger.error(f"Error plotting model interpretation: {str(e)}")
            return None

    def _create_3d_structure_plot(
        self,
        mol: Chem.Mol,
        atom_colors: Dict[int, Tuple[float, float, float]],
    ) -> go.Figure:
        """Create 3D structure plot."""
        # Get 3D coordinates
        conf = mol.GetConformer()
        positions = conf.GetPositions()

        # Create figure
        fig = go.Figure()

        # Add atoms
        fig.add_trace(
            go.Scatter3d(
                x=positions[:, 0],
                y=positions[:, 1],
                z=positions[:, 2],
                mode="markers",
                marker=dict(
                    size=10,
                    color=[
                        f"rgb({r*255},{g*255},{b*255})"
                        for r, g, b in [
                            atom_colors.get(i, (0.5, 0.5, 0.5))
                            for i in range(mol.GetNumAtoms())
                        ]
                    ],
                ),
                text=[
                    f"Atom {i}: {mol.GetAtomWithIdx(i).GetSymbol()}"
                    for i in range(mol.GetNumAtoms())
                ],
                hoverinfo="text",
            )
        )

        # Add bonds
        for bond in mol.GetBonds():
            id1 = bond.GetBeginAtomIdx()
            id2 = bond.GetEndAtomIdx()
            fig.add_trace(
                go.Scatter3d(
                    x=[positions[id1, 0], positions[id2, 0]],
                    y=[positions[id1, 1], positions[id2, 1]],
                    z=[positions[id1, 2], positions[id2, 2]],
                    mode="lines",
                    line=dict(color="gray", width=2),
                    hoverinfo="none",
                )
            )

        # Update layout
        fig.update_layout(
            scene=dict(
                xaxis=dict(showticklabels=False),
                yaxis=dict(showticklabels=False),
                zaxis=dict(showticklabels=False),
            ),
            showlegend=False,
            margin=dict(l=0, r=0, t=0, b=0),
        )

        return fig

    def _create_2d_structure_plot(
        self,
        mol: Chem.Mol,
        atom_colors: Dict[int, Tuple[float, float, float]],
        size: Tuple[int, int],
    ) -> Draw.Image:
        """Create 2D structure plot."""
        drawer = Draw.rdMolDraw2DCairo(size[0], size[1])
        drawer.SetFillPolys(True)
        drawer.SetLineWidth(2)

        # Draw molecule with highlights
        Draw.MolToImage(
            mol,
            size=size,
            highlightAtoms=list(range(mol.GetNumAtoms())),
            highlightColor=atom_colors,
            highlightBonds=True,
        )

        return drawer.GetDrawingText()

    def _create_interactive_abuse_plot(
        self,
        categories: List[str],
        probabilities: List[float],
        uncertainties: Optional[List[float]],
        web_data: Optional[List[Dict]],
        plot_type: str,
    ) -> go.Figure:
        """Create interactive abuse potential plot."""
        if plot_type == "radar":
            fig = go.Figure()

            # Add main trace
            fig.add_trace(
                go.Scatterpolar(
                    r=probabilities,
                    theta=categories,
                    fill="toself",
                    name="Risk Level",
                    line_color=self.colors["high_risk"],
                )
            )

            # Add uncertainty range if available
            if uncertainties:
                lower = [max(0, p - u) for p, u in zip(probabilities, uncertainties)]
                upper = [min(1, p + u) for p, u in zip(probabilities, uncertainties)]

                fig.add_trace(
                    go.Scatterpolar(
                        r=upper,
                        theta=categories,
                        fill="tonext",
                        name="Uncertainty",
                        line_color=self.colors["uncertainty"],
                        opacity=0.3,
                    )
                )

                fig.add_trace(
                    go.Scatterpolar(
                        r=lower,
                        theta=categories,
                        fill="tonext",
                        name="Uncertainty",
                        line_color=self.colors["uncertainty"],
                        opacity=0.3,
                    )
                )

        elif plot_type == "sunburst":
            # Create hierarchical data
            labels = []
            parents = []
            values = []
            colors = []

            # Add root
            labels.append("Abuse Potential")
            parents.append("")
            values.append(sum(probabilities))
            colors.append(self.colors["moderate_risk"])

            # Add categories
            for i, (cat, prob) in enumerate(zip(categories, probabilities)):
                labels.append(cat)
                parents.append("Abuse Potential")
                values.append(prob)
                colors.append(
                    self.colors["high_risk"]
                    if prob >= 0.7
                    else (
                        self.colors["moderate_risk"]
                        if prob >= 0.4
                        else self.colors["low_risk"]
                    )
                )

                # Add web data if available
                if web_data and web_data[i]:
                    for key, value in web_data[i].items():
                        labels.append(f"{cat}: {key}")
                        parents.append(cat)
                        values.append(prob * 0.5)  # Scaled down
                        colors.append(self.colors["receptor_highlight"])

            fig = go.Figure(
                go.Sunburst(
                    labels=labels,
                    parents=parents,
                    values=values,
                    marker=dict(colors=colors),
                )
            )

        else:  # bar
            fig = go.Figure()

            # Add bars
            fig.add_trace(
                go.Bar(
                    x=categories,
                    y=probabilities,
                    marker_color=[
                        (
                            self.colors["high_risk"]
                            if p >= 0.7
                            else (
                                self.colors["moderate_risk"]
                                if p >= 0.4
                                else self.colors["low_risk"]
                            )
                        )
                        for p in probabilities
                    ],
                    error_y=dict(
                        type="data",
                        array=uncertainties if uncertainties else None,
                        visible=True,
                        color=self.colors["uncertainty"],
                    ),
                )
            )

            # Add web data annotations if available
            if web_data:
                for i, (cat, data) in enumerate(zip(categories, web_data)):
                    if data:
                        fig.add_annotation(
                            x=cat,
                            y=probabilities[i],
                            text="<br>".join(f"{k}: {v}" for k, v in data.items()),
                            showarrow=True,
                            arrowhead=1,
                        )

        # Update layout
        fig.update_layout(
            template="plotly_dark" if self.style == "dark" else "plotly_white",
            paper_bgcolor=self.colors["background"],
            plot_bgcolor=self.colors["background"],
            font_color=self.colors["text"],
        )

        return fig

    def _create_static_abuse_plot(
        self,
        categories: List[str],
        probabilities: List[float],
        uncertainties: Optional[List[float]],
        web_data: Optional[List[Dict]],
        plot_type: str,
    ) -> Figure:
        """Create static abuse potential plot."""
        fig = plt.figure(figsize=(10, 6))

        if plot_type == "radar":
            ax = plt.subplot(111, projection="polar")
            angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False)

            # Plot data
            ax.plot(angles, probabilities, "o-")
            ax.fill(angles, probabilities, alpha=0.25)

            # Add uncertainty
            if uncertainties:
                lower = [max(0, p - u) for p, u in zip(probabilities, uncertainties)]
                upper = [min(1, p + u) for p, u in zip(probabilities, uncertainties)]
                ax.fill_between(
                    angles,
                    lower,
                    upper,
                    alpha=0.2,
                    color=self.colors["uncertainty"],
                )

            ax.set_xticks(angles)
            ax.set_xticklabels(categories)

        else:  # bar
            x = np.arange(len(categories))
            plt.bar(
                x,
                probabilities,
                color=[
                    (
                        self.colors["high_risk"]
                        if p >= 0.7
                        else (
                            self.colors["moderate_risk"]
                            if p >= 0.4
                            else self.colors["low_risk"]
                        )
                    )
                    for p in probabilities
                ],
            )

            if uncertainties:
                plt.errorbar(
                    x,
                    probabilities,
                    yerr=uncertainties,
                    fmt="none",
                    color=self.colors["uncertainty"],
                )

            plt.xticks(x, categories, rotation=45, ha="right")

            if web_data:
                for i, data in enumerate(web_data):
                    if data:
                        plt.annotate(
                            "\n".join(f"{k}: {v}" for k, v in data.items()),
                            (i, probabilities[i]),
                            xytext=(0, 10),
                            textcoords="offset points",
                            ha="center",
                            va="bottom",
                        )

        plt.title("Abuse Potential Profile")
        plt.tight_layout()
        return fig

    def _create_interactive_network(
        self,
        G: nx.Graph,
        layout: str,
        show_pharmacophores: bool,
    ) -> go.Figure:
        """Create interactive network plot."""
        if layout == "force":
            pos = nx.spring_layout(G)
        elif layout == "circular":
            pos = nx.circular_layout(G)
        else:
            pos = nx.shell_layout(G)

        # Create figure
        fig = go.Figure()

        # Add edges
        edge_x = []
        edge_y = []
        for edge in G.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])

        fig.add_trace(
            go.Scatter(
                x=edge_x,
                y=edge_y,
                line=dict(width=0.5, color="#888"),
                hoverinfo="none",
                mode="lines",
            )
        )

        # Add nodes
        node_x = []
        node_y = []
        node_text = []
        node_color = []
        node_size = []

        for node in G.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)

            # Create hover text
            text = [
                f"Receptor: {node}",
                f"System: {G.nodes[node]['system']}",
                f"Affinity: {G.nodes[node]['affinity']:.2f}",
            ]
            if show_pharmacophores:
                text.extend(
                    [
                        "Pharmacophores:",
                        *[f"- {p}" for p in G.nodes[node]["pharmacophores"]],
                    ]
                )
            node_text.append("<br>".join(text))

            # Set color and size
            node_color.append(G.nodes[node]["affinity"])
            node_size.append(30 + len(G.nodes[node]["pharmacophores"]) * 5)

        fig.add_trace(
            go.Scatter(
                x=node_x,
                y=node_y,
                mode="markers+text",
                hoverinfo="text",
                text=node_text,
                marker=dict(
                    showscale=True,
                    colorscale=[
                        [0, self.colors["low_risk"]],
                        [0.5, self.colors["moderate_risk"]],
                        [1, self.colors["high_risk"]],
                    ],
                    color=node_color,
                    size=node_size,
                    line_width=2,
                ),
            )
        )

        # Update layout
        fig.update_layout(
            showlegend=False,
            hovermode="closest",
            margin=dict(b=20, l=5, r=5, t=40),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            title="Receptor Interaction Network",
        )

        return fig

    def _create_static_network(
        self,
        G: nx.Graph,
        layout: str,
        show_pharmacophores: bool,
    ) -> Figure:
        """Create static network plot."""
        fig = plt.figure(figsize=(10, 10))

        # Get layout
        if layout == "force":
            pos = nx.spring_layout(G)
        elif layout == "circular":
            pos = nx.circular_layout(G)
        else:
            pos = nx.shell_layout(G)

        # Draw edges
        nx.draw_networkx_edges(
            G,
            pos,
            alpha=0.2,
            edge_color="gray",
        )

        # Draw nodes
        node_colors = [G.nodes[node]["affinity"] for node in G.nodes()]
        node_sizes = [
            2000 + len(G.nodes[node]["pharmacophores"]) * 200 for node in G.nodes()
        ]

        nx.draw_networkx_nodes(
            G,
            pos,
            node_color=node_colors,
            node_size=node_sizes,
            cmap=plt.cm.RdYlBu_r,
        )

        # Add labels
        labels = {}
        for node in G.nodes():
            label = node
            if show_pharmacophores and G.nodes[node]["pharmacophores"]:
                label += f"\n({len(G.nodes[node]['pharmacophores'])} features)"
            labels[node] = label

        nx.draw_networkx_labels(G, pos, labels)

        plt.title("Receptor Interaction Network")
        plt.axis("off")
        return fig

    def _plot_attention_weights(
        self,
        interp_data: Dict,
        mol: Optional[Chem.Mol],
    ) -> Union[Figure, go.Figure]:
        """Plot attention weights visualization."""
        if not mol:
            return self._plot_attention_matrix(interp_data)
        else:
            return self._plot_attention_on_structure(interp_data, mol)

    def _plot_attention_matrix(self, interp_data: Dict) -> go.Figure:
        """Plot attention weight matrix."""
        weights = np.array(interp_data["attention_weights"])

        fig = go.Figure(
            data=go.Heatmap(
                z=weights,
                colorscale=[
                    [0, self.colors["low_risk"]],
                    [0.5, self.colors["moderate_risk"]],
                    [1, self.colors["high_risk"]],
                ],
            )
        )

        fig.update_layout(
            title="Attention Weight Matrix",
            xaxis_title="Target",
            yaxis_title="Source",
        )

        return fig

    def _plot_attention_on_structure(
        self,
        interp_data: Dict,
        mol: Chem.Mol,
    ) -> Draw.Image:
        """Plot attention weights on molecular structure."""
        weights = np.array(interp_data["attention_weights"])

        # Average attention weights per atom
        atom_weights = weights.mean(axis=0)

        # Convert to colors
        atom_colors = {}
        max_weight = max(atom_weights)
        for i, weight in enumerate(atom_weights):
            normalized = weight / max_weight
            if normalized >= 0.7:
                atom_colors[i] = self._hex_to_rgb(self.colors["high_risk"])
            elif normalized >= 0.4:
                atom_colors[i] = self._hex_to_rgb(self.colors["moderate_risk"])
            else:
                atom_colors[i] = self._hex_to_rgb(self.colors["low_risk"])

        return self._create_2d_structure_plot(
            mol,
            atom_colors,
            size=(800, 800),
        )

    def _plot_feature_importance(
        self,
        interp_data: Dict,
    ) -> Union[Figure, go.Figure]:
        """Plot feature importance visualization."""
        features = interp_data["feature_names"]
        importance = interp_data["feature_importance"]

        if self.interactive:
            fig = go.Figure(
                data=go.Bar(
                    x=features,
                    y=importance,
                    marker_color=[
                        (
                            self.colors["high_risk"]
                            if imp >= 0.7
                            else (
                                self.colors["moderate_risk"]
                                if imp >= 0.4
                                else self.colors["low_risk"]
                            )
                        )
                        for imp in importance
                    ],
                )
            )

            fig.update_layout(
                title="Feature Importance",
                xaxis_title="Feature",
                yaxis_title="Importance",
                showlegend=False,
            )

            return fig
        else:
            fig = plt.figure(figsize=(10, 6))
            plt.bar(
                features,
                importance,
                color=[
                    (
                        self.colors["high_risk"]
                        if imp >= 0.7
                        else (
                            self.colors["moderate_risk"]
                            if imp >= 0.4
                            else self.colors["low_risk"]
                        )
                    )
                    for imp in importance
                ],
            )
            plt.xticks(rotation=45, ha="right")
            plt.title("Feature Importance")
            plt.tight_layout()
            return fig

    def _plot_embedding_space(
        self,
        interp_data: Dict,
    ) -> Union[Figure, go.Figure]:
        """Plot embedding space visualization."""
        embeddings = np.array(interp_data["embeddings"])
        labels = interp_data["labels"]

        # Reduce dimensionality
        if embeddings.shape[1] > 2:
            embeddings = TSNE(n_components=2).fit_transform(embeddings)

        if self.interactive:
            fig = go.Figure(
                data=go.Scatter(
                    x=embeddings[:, 0],
                    y=embeddings[:, 1],
                    mode="markers",
                    marker=dict(
                        color=labels,
                        colorscale=[
                            [0, self.colors["low_risk"]],
                            [0.5, self.colors["moderate_risk"]],
                            [1, self.colors["high_risk"]],
                        ],
                        size=10,
                    ),
                    text=labels,
                )
            )

            fig.update_layout(
                title="Embedding Space",
                xaxis_title="Dimension 1",
                yaxis_title="Dimension 2",
                showlegend=False,
            )

            return fig
        else:
            fig = plt.figure(figsize=(8, 8))
            plt.scatter(
                embeddings[:, 0],
                embeddings[:, 1],
                c=labels,
                cmap="RdYlBu_r",
                alpha=0.6,
            )
            plt.title("Embedding Space")
            plt.colorbar(label="Risk Level")
            return fig

    def _hex_to_rgb(self, hex_color: str) -> Tuple[float, float, float]:
        """Convert hex color to RGB tuple."""
        hex_color = hex_color.lstrip("#")
        return tuple(int(hex_color[i : i + 2], 16) / 255.0 for i in (0, 2, 4))

    def _get_risk_based_colors(
        self,
        mol: Chem.Mol,
        predictions: Dict,
    ) -> Dict[int, Tuple[float, float, float]]:
        """Get atom colors based on risk contributions."""
        if "atom_contributions" not in predictions:
            return {i: (0.5, 0.5, 0.5) for i in range(mol.GetNumAtoms())}

        contributions = predictions["atom_contributions"]
        max_contrib = max(abs(min(contributions)), abs(max(contributions)))

        colors = {}
        for i, contrib in enumerate(contributions):
            normalized = contrib / max_contrib if max_contrib > 0 else 0
            if normalized >= 0.7:
                colors[i] = self._hex_to_rgb(self.colors["high_risk"])
            elif normalized >= 0.4:
                colors[i] = self._hex_to_rgb(self.colors["moderate_risk"])
            else:
                colors[i] = self._hex_to_rgb(self.colors["low_risk"])

        return colors

    def _get_receptor_based_colors(
        self,
        mol: Chem.Mol,
        predictions: Dict,
    ) -> Dict[int, Tuple[float, float, float]]:
        """Get atom colors based on receptor interactions."""
        if "receptor_contributions" not in predictions:
            return {i: (0.5, 0.5, 0.5) for i in range(mol.GetNumAtoms())}

        contributions = predictions["receptor_contributions"]
        max_contrib = max(abs(min(contributions)), abs(max(contributions)))

        colors = {}
        for i, contrib in enumerate(contributions):
            normalized = contrib / max_contrib if max_contrib > 0 else 0
            colors[i] = (
                self._hex_to_rgb(self.colors["receptor_highlight"])
                if normalized >= 0.5
                else (0.5, 0.5, 0.5)
            )

        return colors

    def _get_mechanism_based_colors(
        self,
        mol: Chem.Mol,
        predictions: Dict,
    ) -> Dict[int, Tuple[float, float, float]]:
        """Get atom colors based on mechanism contributions."""
        if "mechanism_contributions" not in predictions:
            return {i: (0.5, 0.5, 0.5) for i in range(mol.GetNumAtoms())}

        contributions = predictions["mechanism_contributions"]
        max_contrib = max(abs(min(contributions)), abs(max(contributions)))

        colors = {}
        for i, contrib in enumerate(contributions):
            normalized = contrib / max_contrib if max_contrib > 0 else 0
            if normalized >= 0.7:
                colors[i] = self._hex_to_rgb(self.colors["high_risk"])
            elif normalized >= 0.4:
                colors[i] = self._hex_to_rgb(self.colors["moderate_risk"])
            else:
                colors[i] = self._hex_to_rgb(self.colors["low_risk"])

        return colors

    def _are_receptors_related(self, r1: str, r2: str) -> bool:
        """Check if two receptors are functionally related."""
        # Define receptor families
        families = {
            "serotonin": ["5-HT2A", "5-HT2B", "5-HT2C", "5-HT1A"],
            "dopamine": ["D1", "D2", "D3", "D4", "D5"],
            "opioid": ["mu", "kappa", "delta"],
            "cannabinoid": ["CB1", "CB2"],
            "glutamate": ["NMDA", "AMPA", "kainate"],
            "gaba": ["GABA-A", "GABA-B"],
        }

        # Check if receptors belong to same family
        for family in families.values():
            if r1 in family and r2 in family:
                return True

        return False

    def _create_interactive_mechanism_plot(
        self,
        mechanisms: List[str],
        contributions: List[float],
        categories: List[str],
        plot_type: str,
    ) -> go.Figure:
        """Create interactive mechanism contribution plot."""
        if plot_type == "sankey":
            # Create Sankey diagram
            fig = go.Figure(
                data=[
                    go.Sankey(
                        node=dict(
                            pad=15,
                            thickness=20,
                            line=dict(color="black", width=0.5),
                            label=list(set(categories)) + mechanisms,
                            color=[self.colors["high_risk"]] * len(set(categories))
                            + [self.colors["moderate_risk"]] * len(mechanisms),
                        ),
                        link=dict(
                            source=[categories.index(cat) for cat in categories],
                            target=[
                                len(set(categories)) + i for i in range(len(mechanisms))
                            ],
                            value=contributions,
                            color=[self.colors["receptor_highlight"]] * len(mechanisms),
                        ),
                    )
                ]
            )

        elif plot_type == "sunburst":
            # Create hierarchical data
            labels = ["Mechanisms"] + categories + mechanisms
            parents = [""] + ["Mechanisms"] * len(categories) + categories
            values = (
                [sum(contributions)]
                + [
                    sum(
                        c
                        for c, cat in zip(contributions, categories)
                        if cat == category
                    )
                    for category in set(categories)
                ]
                + contributions
            )

            fig = go.Figure(
                go.Sunburst(
                    labels=labels,
                    parents=parents,
                    values=values,
                    branchvalues="total",
                )
            )

        else:  # treemap
            fig = go.Figure(
                go.Treemap(
                    labels=mechanisms,
                    parents=categories,
                    values=contributions,
                )
            )

        # Update layout
        fig.update_layout(
            title="Mechanism Contributions",
            font_size=10,
        )

        return fig

    def _create_static_mechanism_plot(
        self,
        mechanisms: List[str],
        contributions: List[float],
        categories: List[str],
        plot_type: str,
    ) -> Figure:
        """Create static mechanism contribution plot."""
        fig = plt.figure(figsize=(12, 8))

        if plot_type == "treemap":
            # Create treemap using squarify
            import squarify

            # Normalize sizes
            sizes = np.array(contributions) / sum(contributions) * 1000

            # Create color map
            colors = [
                (
                    self.colors["high_risk"]
                    if c >= 0.7
                    else (
                        self.colors["moderate_risk"]
                        if c >= 0.4
                        else self.colors["low_risk"]
                    )
                )
                for c in contributions
            ]

            # Plot treemap
            squarify.plot(
                sizes=sizes,
                label=mechanisms,
                color=colors,
                alpha=0.7,
                text_kwargs={"fontsize": 8},
            )

        else:  # hierarchical
            from matplotlib.patches import Rectangle

            # Group by category
            category_data = {}
            for cat, mech, cont in zip(categories, mechanisms, contributions):
                if cat not in category_data:
                    category_data[cat] = []
                category_data[cat].append((mech, cont))

            # Create hierarchical bar plot
            x = 0
            width = 0.8
            for cat, data in category_data.items():
                # Plot category bar
                cat_height = sum(cont for _, cont in data)
                plt.bar(
                    x,
                    cat_height,
                    width,
                    label=cat,
                    color=self.colors["high_risk"],
                    alpha=0.3,
                )

                # Plot mechanism bars
                y = 0
                for mech, cont in data:
                    plt.bar(
                        x,
                        cont,
                        width * 0.8,
                        bottom=y,
                        label=mech,
                        color=self.colors["moderate_risk"],
                    )
                    y += cont

                x += 1

            plt.xticks(
                range(len(category_data)),
                list(category_data.keys()),
                rotation=45,
                ha="right",
            )

        plt.title("Mechanism Contributions")
        plt.tight_layout()
        return fig
