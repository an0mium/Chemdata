"""Advanced visualization for molecular structure and activity prediction.

This module provides comprehensive visualization capabilities:
1. 3D structure visualization with pharmacophore features
2. Feature distribution and correlation analysis 
3. Dimensionality reduction and clustering
4. Attention and feature importance visualization
5. Uncertainty quantification and visualization
6. Interactive reports and dashboards
"""

import logging
from typing import Dict, List, Optional, Union, Tuple
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN

from .activity import ActivityPredictor
from .base import MLProcessor
from .features import MolecularFeaturizer
from .utils import mol_to_graph
from .predictors.psychoactive import PsychoactivePredictor
from .predictors.nootropic import NootropicPredictor


class MLVisualizer(MLProcessor):
    """Advanced visualization tools for ML predictions."""

    def __init__(
        self,
        activity_predictor: Optional[ActivityPredictor] = None,
        psychoactive_predictor: Optional[PsychoactivePredictor] = None,
        nootropic_predictor: Optional[NootropicPredictor] = None,
        featurizer: Optional[MolecularFeaturizer] = None,
    ):
        """Initialize visualizer.

        Args:
            activity_predictor: Activity prediction model
            psychoactive_predictor: Psychoactive effects predictor
            nootropic_predictor: Nootropic effects predictor
            featurizer: Molecular featurizer
        """
        super().__init__()
        self.activity_predictor = activity_predictor
        self.psychoactive_predictor = psychoactive_predictor
        self.nootropic_predictor = nootropic_predictor
        self.featurizer = featurizer or MolecularFeaturizer()

    def visualize_3d(
        self,
        compound: Union[str, Chem.Mol],
        show_pharmacophore: bool = True,
        show_surface: bool = True,
        highlight_substructures: Optional[List[str]] = None,
        show_predictions: bool = True,
    ) -> py3Dmol.view:
        """Generate enhanced 3D visualization."""
        try:
            mol = (
                compound
                if isinstance(compound, Chem.Mol)
                else Chem.MolFromSmiles(compound)
            )
            if mol is None:
                raise ValueError("Invalid compound input")

            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)

            view = py3Dmol.view(width=800, height=600)
            view.addModel(Chem.MolToMolBlock(mol), "mol")
            view.setStyle({"stick": {}})

            if show_surface:
                view.addSurface(
                    py3Dmol.VDW, {"opacity": 0.7, "colorscheme": {"gradient": "rwb"}}
                )

            if show_pharmacophore and self.predictor:
                features = self.predictor.pharmacophore.generate(mol)
                for feature in features:
                    pos = feature["position"]
                    color = {
                        "donor": "blue",
                        "acceptor": "red",
                        "aromatic": "purple",
                        "hydrophobic": "green",
                        "positive": "cyan",
                        "negative": "magenta",
                    }.get(feature["type"], "gray")

                    view.addSphere(
                        {
                            "center": {"x": pos[0], "y": pos[1], "z": pos[2]},
                            "radius": 0.5,
                            "color": color,
                            "alpha": 0.7,
                        }
                    )

            if highlight_substructures:
                for i, smarts in enumerate(highlight_substructures):
                    pattern = Chem.MolFromSmarts(smarts)
                    if pattern:
                        matches = mol.GetSubstructMatches(pattern)
                        for match in matches:
                            color = px.colors.qualitative.Set3[i % 12]
                            for atom_idx in match:
                                pos = mol.GetConformer().GetAtomPosition(atom_idx)
                                view.addSphere(
                                    {
                                        "center": {"x": pos.x, "y": pos.y, "z": pos.z},
                                        "radius": 0.4,
                                        "color": color,
                                        "alpha": 0.8,
                                    }
                                )

            if show_predictions and self.predictor:
                pred = self.predictor.predict_compound(compound)
                self._add_prediction_labels(view, pred)

            view.zoomTo()
            return view

        except Exception as e:
            self.logger.error(f"Error in visualize_3d: {str(e)}")
            return None

    def plot_feature_distributions(
        self,
        features: Dict[str, Union[np.ndarray, Dict[str, float]]],
        output_path: Optional[str] = None,
        figsize: Tuple[int, int] = (12, 8),
    ) -> None:
        """Plot feature distributions."""
        try:
            fig = plt.figure(figsize=figsize)

            for i, (feat_type, feat_data) in enumerate(features.items(), 1):
                if feat_data is None:
                    continue

                data = (
                    feat_data.flatten()
                    if isinstance(feat_data, np.ndarray)
                    else list(feat_data.values())
                )

                plt.subplot(len(features), 1, i)
                sns.histplot(data=data, kde=True)
                plt.title(f"{feat_type} Distribution")
                plt.xlabel("Value")
                plt.ylabel("Count")

            plt.tight_layout()
            if output_path:
                plt.savefig(output_path)
                plt.close()
            else:
                plt.show()

        except Exception as e:
            self.logger.error(f"Error plotting feature distributions: {str(e)}")

    def plot_feature_correlations(
        self,
        features: Dict[str, Union[np.ndarray, Dict[str, float]]],
        output_path: Optional[str] = None,
        figsize: Tuple[int, int] = (10, 10),
    ) -> None:
        """Plot correlation matrix of features."""
        try:
            feature_matrix = []
            feature_names = []

            for feat_type, feat_data in features.items():
                if feat_data is None:
                    continue

                if isinstance(feat_data, np.ndarray):
                    if feat_data.ndim == 1:
                        feature_matrix.append(feat_data.reshape(-1, 1))
                        feature_names.append(feat_type)
                    else:
                        feature_matrix.append(feat_data)
                        feature_names.extend(
                            [f"{feat_type}_{i}" for i in range(feat_data.shape[1])]
                        )
                elif isinstance(feat_data, dict):
                    values = np.array(list(feat_data.values())).reshape(-1, 1)
                    feature_matrix.append(values)
                    feature_names.extend(feat_data.keys())

            if not feature_matrix:
                return

            feature_matrix = np.hstack(feature_matrix)
            corr_matrix = np.corrcoef(feature_matrix.T)

            plt.figure(figsize=figsize)
            sns.heatmap(
                corr_matrix,
                xticklabels=feature_names,
                yticklabels=feature_names,
                cmap="coolwarm",
                center=0,
                annot=True,
                fmt=".2f",
            )
            plt.title("Feature Correlations")
            plt.xticks(rotation=45, ha="right")
            plt.yticks(rotation=0)
            plt.tight_layout()

            if output_path:
                plt.savefig(output_path)
                plt.close()
            else:
                plt.show()

        except Exception as e:
            self.logger.error(f"Error plotting feature correlations: {str(e)}")

    def plot_dimensionality_reduction(
        self,
        features: Dict[str, Union[np.ndarray, Dict[str, float]]],
        method: str = "umap",
        n_components: int = 2,
        labels: Optional[np.ndarray] = None,
        output_path: Optional[str] = None,
        figsize: Tuple[int, int] = (10, 10),
    ) -> Optional[np.ndarray]:
        """Plot dimensionality reduction."""
        try:
            feature_matrix = []
            for feat_data in features.values():
                if feat_data is None:
                    continue

                if isinstance(feat_data, np.ndarray):
                    feature_matrix.append(feat_data.reshape(feat_data.shape[0], -1))
                else:
                    feature_matrix.append(
                        np.array(list(feat_data.values())).reshape(-1, 1)
                    )

            if not feature_matrix:
                return None

            feature_matrix = np.hstack(feature_matrix)
            scaler = StandardScaler()
            scaled_features = scaler.fit_transform(feature_matrix)

            if method == "pca":
                reducer = PCA(n_components=n_components)
            elif method == "tsne":
                reducer = TSNE(n_components=n_components)
            elif method == "umap":
                reducer = UMAP(n_components=n_components)
            else:
                raise ValueError(f"Unknown reduction method: {method}")

            transformed = reducer.fit_transform(scaled_features)

            plt.figure(figsize=figsize)
            if labels is not None:
                scatter = plt.scatter(
                    transformed[:, 0],
                    transformed[:, 1],
                    c=labels,
                    cmap="viridis",
                )
                plt.colorbar(scatter)
            else:
                plt.scatter(transformed[:, 0], transformed[:, 1])

            plt.title(f"{method.upper()} Projection")
            plt.xlabel("Component 1")
            plt.ylabel("Component 2")
            plt.tight_layout()

            if output_path:
                plt.savefig(output_path)
                plt.close()
            else:
                plt.show()

            return transformed

        except Exception as e:
            self.logger.error(f"Error plotting dimensionality reduction: {str(e)}")
            return None

    def plot_uncertainty_analysis(
        self,
        compounds: List[Union[str, Chem.Mol]],
        n_samples: int = 30,
    ) -> Dict[str, go.Figure]:
        """Generate comprehensive uncertainty visualization."""
        try:
            if not self.predictor:
                return {}

            figures = {}
            results = []
            for compound in compounds:
                pred = self.predictor.predict_with_uncertainty(
                    compound, n_samples=n_samples
                )
                results.append(pred)

            # Binding affinity uncertainty
            fig_binding = go.Figure()
            binding_means = [r["binding"]["mean"] for r in results]
            binding_stds = [r["binding"]["std"] for r in results]

            fig_binding.add_trace(
                go.Scatter(
                    x=list(range(len(compounds))),
                    y=binding_means,
                    error_y=dict(type="data", array=binding_stds, visible=True),
                    mode="markers",
                    name="Binding Affinity",
                )
            )

            fig_binding.update_layout(
                title="Binding Affinity Uncertainty",
                xaxis_title="Compound Index",
                yaxis_title="Predicted Value (pKi)",
            )
            figures["binding"] = fig_binding

            # Activity probabilities
            fig_activity = go.Figure()
            for i, activity in enumerate(self.predictor.ACTIVITY_TYPES):
                means = [r["activity"]["mean"][i] for r in results]
                stds = [r["activity"]["std"][i] for r in results]

                fig_activity.add_trace(
                    go.Box(
                        y=means,
                        name=activity,
                        boxpoints="all",
                        jitter=0.3,
                        pointpos=-1.8,
                    )
                )

            fig_activity.update_layout(
                title="Activity Probability Distributions",
                yaxis_title="Probability",
                showlegend=True,
            )
            figures["activity"] = fig_activity

            return figures

        except Exception as e:
            self.logger.error(f"Error in plot_uncertainty_analysis: {str(e)}")
            return {}

    def create_interactive_report(
        self,
        compound: Union[str, Chem.Mol],
        include_ml: bool = True,
    ) -> Dict:
        """Generate comprehensive visualization report."""
        try:
            report = {}

            # 3D structure
            report["structure_3d"] = self.visualize_3d(
                compound, show_predictions=include_ml
            )

            # Feature analysis
            if self.featurizer:
                mol = (
                    compound
                    if isinstance(compound, Chem.Mol)
                    else Chem.MolFromSmiles(compound)
                )
                features = self.featurizer.generate_features(mol)

                self.plot_feature_distributions(
                    {"molecular": features}, output_path="feature_distributions.png"
                )
                report["distributions"] = "feature_distributions.png"

                self.plot_feature_correlations(
                    {"molecular": features}, output_path="feature_correlations.png"
                )
                report["correlations"] = "feature_correlations.png"

                reduced = self.plot_dimensionality_reduction(
                    {"molecular": features},
                    method="umap",
                    output_path="dim_reduction.png",
                )
                report["embedding"] = "dim_reduction.png"

            # ML analysis
            if include_ml and self.predictor:
                uncertainty_plots = self.plot_uncertainty_analysis([compound])
                report.update(uncertainty_plots)

            return report

        except Exception as e:
            self.logger.error(f"Error in create_interactive_report: {str(e)}")
            return {}

    def _add_prediction_labels(self, view: py3Dmol.view, pred: Dict) -> None:
        """Add prediction labels to 3D view."""
        view.addLabel(
            "Predicted Activities:",
            {
                "position": {"x": -5, "y": -5, "z": 0},
                "backgroundColor": "white",
                "fontColor": "black",
            },
        )

        y_offset = -4
        for activity, prob in pred.get("activities", []):
            if prob > 0.5:
                view.addLabel(
                    f"{activity}: {prob:.2f}",
                    {
                        "position": {"x": -5, "y": y_offset, "z": 0},
                        "backgroundColor": "white",
                        "fontColor": "black",
                    },
                )
                y_offset += 1

    def create_binding_plot(self, data: Dict, colorscale: str = "RdBu") -> go.Figure:
        """Create binding affinity prediction plot.

        Args:
            data: Binding prediction data
            colorscale: Color scale for heatmap

        Returns:
            Plotly figure
        """
        try:
            if not isinstance(data, dict) or "values" not in data:
                raise ValueError("Invalid binding data format")

            values = data["values"]
            receptors = data.get("receptors", list(range(len(values))))

            fig = go.Figure(
                data=go.Heatmap(
                    z=[values],
                    x=receptors,
                    colorscale=colorscale,
                    showscale=True,
                )
            )

            fig.update_layout(
                title="Binding Affinity Predictions",
                xaxis_title="Receptor",
                yaxis_title="Affinity (pKi)",
                showlegend=True,
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating binding plot: {str(e)}")
            return go.Figure()

    def create_abuse_plot(self, data: Dict, colorscale: str = "Reds") -> go.Figure:
        """Create abuse potential prediction plot.

        Args:
            data: Abuse prediction data
            colorscale: Color scale for radar plot

        Returns:
            Plotly figure
        """
        try:
            if not isinstance(data, dict) or "probabilities" not in data:
                raise ValueError("Invalid abuse data format")

            categories = list(data["probabilities"].keys())
            values = list(data["probabilities"].values())
            values.append(values[0])  # Close the polygon
            categories.append(categories[0])

            fig = go.Figure(
                data=go.Scatterpolar(
                    r=values,
                    theta=categories,
                    fill="toself",
                    marker=dict(color=colorscale),
                )
            )

            fig.update_layout(
                title="Abuse Potential Assessment",
                showlegend=False,
                polar=dict(
                    radialaxis=dict(
                        visible=True,
                        range=[0, 1],
                    )
                ),
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating abuse plot: {str(e)}")
            return go.Figure()

    def create_toxicity_plot(
        self, data: Dict, colorscale: str = "RdYlGn_r"
    ) -> go.Figure:
        """Create toxicity prediction plot.

        Args:
            data: Toxicity prediction data
            colorscale: Color scale for bar plot

        Returns:
            Plotly figure
        """
        try:
            if not isinstance(data, dict) or "risks" not in data:
                raise ValueError("Invalid toxicity data format")

            endpoints = list(data["risks"].keys())
            risks = list(data["risks"].values())

            colors = px.colors.sample_colorscale(
                colorscale, np.linspace(0, 1, len(risks))
            )

            fig = go.Figure(
                data=go.Bar(
                    x=endpoints,
                    y=risks,
                    marker_color=colors,
                )
            )

            fig.update_layout(
                title="Toxicity Risk Assessment",
                xaxis_title="Endpoint",
                yaxis_title="Risk Level",
                showlegend=False,
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating toxicity plot: {str(e)}")
            return go.Figure()

    def create_activity_plot(
        self, data: Dict, colorscale: str = "Viridis"
    ) -> go.Figure:
        """Create activity prediction plot.

        Args:
            data: Activity prediction data
            colorscale: Color scale for scatter plot

        Returns:
            Plotly figure
        """
        try:
            if not isinstance(data, dict) or "activities" not in data:
                raise ValueError("Invalid activity data format")

            activities = data["activities"]
            targets = list(activities.keys())
            values = list(activities.values())

            colors = px.colors.sample_colorscale(
                colorscale, np.linspace(0, 1, len(values))
            )

            fig = go.Figure(
                data=go.Scatter(
                    x=targets,
                    y=values,
                    mode="markers+lines",
                    marker=dict(color=colors),
                )
            )

            fig.update_layout(
                title="Activity Profile",
                xaxis_title="Target",
                yaxis_title="Activity",
                showlegend=False,
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating activity plot: {str(e)}")
            return go.Figure()

    def create_psychoactive_plot(
        self, data: Dict, colorscale: str = "Plasma"
    ) -> go.Figure:
        """Create psychoactive effects prediction plot.

        Args:
            data: Psychoactive prediction data
            colorscale: Color scale for radar plot

        Returns:
            Plotly figure
        """
        try:
            if not isinstance(data, dict) or "effects" not in data:
                raise ValueError("Invalid psychoactive data format")

            effects = list(data["effects"].keys())
            values = list(data["effects"].values())
            values.append(values[0])  # Close the polygon
            effects.append(effects[0])

            fig = go.Figure(
                data=go.Scatterpolar(
                    r=values,
                    theta=effects,
                    fill="toself",
                    marker=dict(color=colorscale),
                )
            )

            fig.update_layout(
                title="Psychoactive Effects Profile",
                showlegend=False,
                polar=dict(
                    radialaxis=dict(
                        visible=True,
                        range=[0, 1],
                    )
                ),
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating psychoactive plot: {str(e)}")
            return go.Figure()

    def create_nootropic_plot(self, data: Dict, colorscale: str = "Magma") -> go.Figure:
        """Create nootropic effects prediction plot.

        Args:
            data: Nootropic prediction data
            colorscale: Color scale for radar plot

        Returns:
            Plotly figure
        """
        try:
            if not isinstance(data, dict) or "effects" not in data:
                raise ValueError("Invalid nootropic data format")

            effects = list(data["effects"].keys())
            values = list(data["effects"].values())
            values.append(values[0])  # Close the polygon
            effects.append(effects[0])

            fig = go.Figure(
                data=go.Scatterpolar(
                    r=values,
                    theta=effects,
                    fill="toself",
                    marker=dict(color=colorscale),
                )
            )

            fig.update_layout(
                title="Nootropic Effects Profile",
                showlegend=False,
                polar=dict(
                    radialaxis=dict(
                        visible=True,
                        range=[0, 1],
                    )
                ),
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating nootropic plot: {str(e)}")
            return go.Figure()

    def create_interaction_network(
        self, data: Dict, colorscale: str = "Inferno"
    ) -> go.Figure:
        """Create interaction network plot.

        Args:
            data: Interaction network data
            colorscale: Color scale for network plot

        Returns:
            Plotly figure
        """
        try:
            if not isinstance(data, dict) or "interactions" not in data:
                raise ValueError("Invalid interaction data format")

            # Create network graph
            G = nx.Graph()
            for interaction in data["interactions"]:
                G.add_edge(
                    interaction["source"],
                    interaction["target"],
                    weight=interaction.get("weight", 1),
                )

            # Calculate layout
            pos = nx.spring_layout(G)

            # Create edge traces
            edge_x = []
            edge_y = []
            for edge in G.edges():
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])

            edge_trace = go.Scatter(
                x=edge_x,
                y=edge_y,
                line=dict(width=0.5, color="#888"),
                hoverinfo="none",
                mode="lines",
            )

            # Create node traces
            node_x = []
            node_y = []
            for node in G.nodes():
                x, y = pos[node]
                node_x.append(x)
                node_y.append(y)

            node_trace = go.Scatter(
                x=node_x,
                y=node_y,
                mode="markers+text",
                hoverinfo="text",
                text=list(G.nodes()),
                marker=dict(
                    showscale=True,
                    colorscale=colorscale,
                    size=10,
                    colorbar=dict(
                        thickness=15,
                        title="Node Connections",
                        xanchor="left",
                        titleside="right",
                    ),
                ),
            )

            # Create figure
            fig = go.Figure(
                data=[edge_trace, node_trace],
                layout=go.Layout(
                    title="Interaction Network",
                    showlegend=False,
                    hovermode="closest",
                    margin=dict(b=20, l=5, r=5, t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                ),
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating interaction network: {str(e)}")
            return go.Figure()

    def create_similarity_matrix(
        self, data: Dict, colorscale: str = "Blues"
    ) -> go.Figure:
        """Create similarity matrix plot.

        Args:
            data: Similarity matrix data
            colorscale: Color scale for heatmap

        Returns:
            Plotly figure
        """
        try:
            if not isinstance(data, dict) or "matrix" not in data:
                raise ValueError("Invalid similarity data format")

            matrix = data["matrix"]
            labels = data.get("labels", list(range(len(matrix))))

            fig = go.Figure(
                data=go.Heatmap(
                    z=matrix,
                    x=labels,
                    y=labels,
                    colorscale=colorscale,
                )
            )

            fig.update_layout(
                title="Structural Similarity Matrix",
                xaxis_title="Compound",
                yaxis_title="Compound",
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating similarity matrix: {str(e)}")
            return go.Figure()

    def create_uncertainty_plot(
        self, data: Dict, colorscale: str = "Greys"
    ) -> go.Figure:
        """Create uncertainty visualization plot.

        Args:
            data: Uncertainty data
            colorscale: Color scale for violin plot

        Returns:
            Plotly figure
        """
        try:
            if not isinstance(data, dict) or "distributions" not in data:
                raise ValueError("Invalid uncertainty data format")

            fig = go.Figure()

            for model, dist in data["distributions"].items():
                fig.add_trace(
                    go.Violin(
                        y=dist,
                        name=model,
                        box_visible=True,
                        meanline_visible=True,
                        points="all",
                    )
                )

            fig.update_layout(
                title="Prediction Uncertainty",
                xaxis_title="Model",
                yaxis_title="Prediction Distribution",
                showlegend=True,
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating uncertainty plot: {str(e)}")
            return go.Figure()
