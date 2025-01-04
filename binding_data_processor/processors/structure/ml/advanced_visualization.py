"""Advanced visualization capabilities for compound analysis.

This module provides specialized visualization tools for:
1. Toxicity predictions and risk assessment
2. Abuse potential analysis
3. Psychopharmacological activity profiling
4. Structure-activity relationship visualization
5. Nootropic activity analysis
6. Web data visualization
"""

import logging
from typing import Dict, List, Optional, Union, Tuple
from functools import lru_cache
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from rdkit import Chem
from wordcloud import WordCloud
import networkx as nx
from datetime import datetime, timedelta

from .base import MLProcessor
from .activity import ActivityPredictor
from .predictors.nootropic import NootropicPredictor
from web_enrichment.data_sources.community import CommunityClient
from web_enrichment.data_sources.social import SocialDataHarvester


class AdvancedVisualizer(MLProcessor):
    """Advanced visualization tools for compound analysis."""

    def __init__(
        self,
        predictor: Optional[ActivityPredictor] = None,
        nootropic_predictor: Optional[NootropicPredictor] = None,
        community_client: Optional[CommunityClient] = None,
        social_harvester: Optional[SocialDataHarvester] = None,
    ):
        """Initialize visualizer.

        Args:
            predictor: Activity prediction model
            nootropic_predictor: Nootropic activity predictor
            community_client: Community data client
            social_harvester: Social media data harvester
        """
        super().__init__()
        self.predictor = predictor
        self.nootropic_predictor = nootropic_predictor
        self.community_client = community_client
        self.social_harvester = social_harvester

        # Initialize logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        )
        self.logger.addHandler(handler)

    @lru_cache(maxsize=32)
    def plot_toxicity_predictions(
        self,
        compound: Union[str, Chem.Mol],
        include_confidence: bool = True,
    ) -> go.Figure:
        """Visualize toxicity predictions.

        Args:
            compound: Input compound
            include_confidence: Whether to show confidence intervals

        Returns:
            Plotly figure with toxicity visualization
        """
        try:
            if not self.predictor:
                return None

            pred = self.predictor.predict_toxicity(compound)

            fig = go.Figure()

            # Plot toxicity probabilities
            categories = [
                "Hepatotoxicity",
                "Cardiotoxicity",
                "Nephrotoxicity",
                "Neurotoxicity",
                "Reproductive Toxicity",
                "Mutagenicity",
                "Carcinogenicity",
            ]
            values = [pred[c]["probability"] for c in categories]

            if include_confidence:
                error_plus = [
                    pred[c].get("ci_upper", 0) - pred[c]["probability"]
                    for c in categories
                ]
                error_minus = [
                    pred[c]["probability"] - pred[c].get("ci_lower", 0)
                    for c in categories
                ]

                fig.add_trace(
                    go.Bar(
                        x=categories,
                        y=values,
                        error_y=dict(
                            type="data",
                            symmetric=False,
                            array=error_plus,
                            arrayminus=error_minus,
                        ),
                        marker_color=["red" if v > 0.5 else "green" for v in values],
                    )
                )
            else:
                fig.add_trace(
                    go.Bar(
                        x=categories,
                        y=values,
                        marker_color=["red" if v > 0.5 else "green" for v in values],
                    )
                )

            fig.update_layout(
                title="Predicted Toxicity Risks",
                xaxis_title="Risk Category",
                yaxis_title="Risk Probability",
                yaxis_range=[0, 1],
                showlegend=False,
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error plotting toxicity: {str(e)}")
            return None

    @lru_cache(maxsize=32)
    def plot_abuse_potential(
        self,
        compound: Union[str, Chem.Mol],
        include_mechanisms: bool = True,
    ) -> go.Figure:
        """Visualize abuse potential predictions.

        Args:
            compound: Input compound
            include_mechanisms: Whether to show mechanism predictions

        Returns:
            Plotly figure with abuse potential visualization
        """
        try:
            if not self.predictor:
                return None

            pred = self.predictor.predict_abuse_potential(compound)

            fig = go.Figure()

            # Core abuse metrics
            categories = [
                "Reward Potential",
                "Physical Dependence",
                "Psychological Dependence",
                "Withdrawal Severity",
                "Tolerance Development",
            ]
            values = [pred["core"][c] for c in categories]

            fig.add_trace(
                go.Scatterpolar(
                    r=values,
                    theta=categories,
                    fill="toself",
                    name="Core Metrics",
                )
            )

            if include_mechanisms:
                # Mechanism involvement
                mech_categories = [
                    "Dopamine",
                    "Serotonin",
                    "Norepinephrine",
                    "GABA",
                    "Glutamate",
                    "Opioid",
                    "Cannabinoid",
                ]
                mech_values = [pred["mechanisms"][c] for c in mech_categories]

                fig.add_trace(
                    go.Scatterpolar(
                        r=mech_values,
                        theta=mech_categories,
                        fill="toself",
                        name="Mechanisms",
                    )
                )

            fig.update_layout(
                title="Predicted Abuse Potential Profile",
                polar=dict(radialaxis=dict(visible=True, range=[0, 1])),
                showlegend=True,
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error plotting abuse potential: {str(e)}")
            return None

    @lru_cache(maxsize=32)
    def plot_psychopharm_activity(
        self,
        compound: Union[str, Chem.Mol],
        include_subtypes: bool = True,
    ) -> go.Figure:
        """Visualize predicted psychopharmacological activity.

        Args:
            compound: Input compound
            include_subtypes: Whether to show receptor subtype predictions

        Returns:
            Plotly figure with activity visualization
        """
        try:
            if not self.predictor:
                return None

            pred = self.predictor.predict_psychopharm_activity(compound)

            # Create subplots
            fig = go.Figure()

            # Main activity classes
            categories = [
                "Antidepressant",
                "Anxiolytic",
                "Antipsychotic",
                "Stimulant",
                "Sedative",
                "Hallucinogenic",
                "Dissociative",
                "Entactogenic",
            ]
            values = [pred["main"][c] for c in categories]

            fig.add_trace(
                go.Bar(
                    x=categories,
                    y=values,
                    marker_color=px.colors.qualitative.Set3,
                    name="Main Classes",
                )
            )

            if include_subtypes:
                # Add receptor subtype predictions
                subtypes = {
                    "Serotonin": ["5-HT1A", "5-HT2A", "5-HT2C", "5-HT3", "5-HT7"],
                    "Dopamine": ["D1", "D2", "D3", "D4"],
                    "NMDA": ["NR1", "NR2A", "NR2B"],
                    "GABA": ["GABA-A", "GABA-B"],
                    "Opioid": ["mu", "kappa", "delta"],
                }

                for receptor, subtypes in subtypes.items():
                    values = [pred["subtypes"][receptor][s] for s in subtypes]

                    fig.add_trace(
                        go.Bar(
                            x=subtypes,
                            y=values,
                            name=receptor,
                            visible="legendonly",
                        )
                    )

            fig.update_layout(
                title="Predicted Psychopharmacological Activity Profile",
                xaxis_title="Activity Type",
                yaxis_title="Probability",
                barmode="group",
                showlegend=True,
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error plotting psychopharm activity: {str(e)}")
            return None

    @lru_cache(maxsize=32)
    def plot_nootropic_activity(
        self,
        compound: Union[str, Chem.Mol],
        include_mechanisms: bool = True,
    ) -> go.Figure:
        """Visualize predicted nootropic activity.

        Args:
            compound: Input compound
            include_mechanisms: Whether to show mechanism predictions

        Returns:
            Plotly figure with nootropic activity visualization
        """
        try:
            if not self.nootropic_predictor:
                return None

            pred = self.nootropic_predictor.predict_activity(compound)

            fig = go.Figure()

            # Core cognitive effects
            categories = [
                "Memory Enhancement",
                "Focus/Attention",
                "Learning Speed",
                "Mental Energy",
                "Cognitive Flexibility",
                "Neuroprotection",
            ]
            values = [pred["effects"][c] for c in categories]

            fig.add_trace(
                go.Scatterpolar(
                    r=values,
                    theta=categories,
                    fill="toself",
                    name="Cognitive Effects",
                )
            )

            if include_mechanisms:
                # Mechanism involvement
                mech_categories = [
                    "Cholinergic",
                    "Glutamatergic",
                    "Dopaminergic",
                    "Serotonergic",
                    "Nootropic",
                    "Anti-inflammatory",
                    "Antioxidant",
                ]
                mech_values = [pred["mechanisms"][c] for c in mech_categories]

                fig.add_trace(
                    go.Scatterpolar(
                        r=mech_values,
                        theta=mech_categories,
                        fill="toself",
                        name="Mechanisms",
                    )
                )

            fig.update_layout(
                title="Predicted Nootropic Activity Profile",
                polar=dict(radialaxis=dict(visible=True, range=[0, 1])),
                showlegend=True,
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error plotting nootropic activity: {str(e)}")
            return None

    def plot_community_data(
        self,
        compound_name: str,
        plot_type: str = "timeline",
    ) -> go.Figure:
        """Visualize community data.

        Args:
            compound_name: Name of compound
            plot_type: Type of visualization (timeline, network, heatmap, wordcloud)

        Returns:
            Plotly figure with community data visualization
        """
        try:
            if not self.community_client:
                return None

            data = self.community_client.get_compound_data(compound_name)

            if plot_type == "timeline":
                return self._create_timeline_plot(data)
            elif plot_type == "network":
                return self._create_network_plot(data)
            elif plot_type == "heatmap":
                return self._create_heatmap_plot(data)
            elif plot_type == "wordcloud":
                return self._create_wordcloud_plot(data)
            else:
                raise ValueError(f"Unknown plot type: {plot_type}")

        except Exception as e:
            self.logger.error(f"Error plotting community data: {str(e)}")
            return None

    def plot_social_data(
        self,
        compound_name: str,
        plot_type: str = "timeline",
        time_range: Optional[Tuple[datetime, datetime]] = None,
    ) -> go.Figure:
        """Visualize social media data.

        Args:
            compound_name: Name of compound
            plot_type: Type of visualization (timeline, network, heatmap, wordcloud)
            time_range: Optional time range to filter data

        Returns:
            Plotly figure with social data visualization
        """
        try:
            if not self.social_harvester:
                return None

            data = self.social_harvester.get_compound_mentions(
                compound_name, time_range=time_range
            )

            if plot_type == "timeline":
                return self._create_timeline_plot(data)
            elif plot_type == "network":
                return self._create_network_plot(data)
            elif plot_type == "heatmap":
                return self._create_heatmap_plot(data)
            elif plot_type == "wordcloud":
                return self._create_wordcloud_plot(data)
            else:
                raise ValueError(f"Unknown plot type: {plot_type}")

        except Exception as e:
            self.logger.error(f"Error plotting social data: {str(e)}")
            return None

    def _create_timeline_plot(self, data: Dict) -> go.Figure:
        """Create timeline visualization."""
        try:
            events = []
            for item in data.get("timeline", []):
                events.append(
                    {
                        "date": item["date"],
                        "type": item["type"],
                        "description": item["description"],
                    }
                )

            # Sort events by date
            events.sort(key=lambda x: x["date"])

            fig = go.Figure()

            # Add events as scatter points
            fig.add_trace(
                go.Scatter(
                    x=[e["date"] for e in events],
                    y=[e["type"] for e in events],
                    mode="markers+text",
                    text=[e["description"] for e in events],
                    textposition="top center",
                )
            )

            fig.update_layout(
                title="Timeline of Events",
                xaxis_title="Date",
                yaxis_title="Event Type",
                showlegend=False,
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating timeline plot: {str(e)}")
            return None

    def _create_network_plot(self, data: Dict) -> go.Figure:
        """Create network visualization."""
        try:
            # Create network graph
            G = nx.Graph()

            # Add nodes and edges
            for connection in data.get("connections", []):
                G.add_edge(
                    connection["from"],
                    connection["to"],
                    weight=connection.get("weight", 1),
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
                    colorscale="YlGnBu",
                    size=10,
                ),
            )

            # Create figure
            fig = go.Figure(
                data=[edge_trace, node_trace],
                layout=go.Layout(
                    title="Network Visualization",
                    showlegend=False,
                    hovermode="closest",
                    margin=dict(b=20, l=5, r=5, t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                ),
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating network plot: {str(e)}")
            return None

    def _create_heatmap_plot(self, data: Dict) -> go.Figure:
        """Create heatmap visualization."""
        try:
            # Extract categories and values
            categories = set()
            sources = list(data.keys())
            values = {}

            for source, source_data in data.items():
                if isinstance(source_data, dict):
                    for category, value in source_data.get("metrics", {}).items():
                        categories.add(category)
                        values[(source, category)] = value

            categories = sorted(categories)

            # Create matrix
            matrix = []
            for source in sources:
                row = []
                for category in categories:
                    row.append(values.get((source, category), 0))
                matrix.append(row)

            # Create heatmap
            fig = go.Figure(
                data=go.Heatmap(
                    z=matrix,
                    x=categories,
                    y=sources,
                    colorscale="RdBu",
                )
            )

            fig.update_layout(
                title="Data Heatmap",
                xaxis_title="Category",
                yaxis_title="Source",
            )

            return fig

        except Exception as e:
            self.logger.error(f"Error creating heatmap plot: {str(e)}")
            return None

    def _create_wordcloud_plot(self, data: Dict) -> go.Figure:
        """Create word cloud visualization."""
        try:
            # Extract text data
            text = ""
            for source_data in data.values():
                if isinstance(source_data, dict):
                    # Add mentions
                    for mention in source_data.get("mentions", []):
                        text += f" {mention['text']}"
                    # Add descriptions
                    for desc in source_data.get("descriptions", []):
                        text += f" {desc}"
                    # Add effects
                    for effect in source_data.get("effects", []):
                        text += f" {effect}"

            # Generate word cloud
            if text:
                wordcloud = WordCloud(
                    width=800,
                    height=400,
                    background_color="white",
                ).generate(text)

                # Convert to figure
                fig = go.Figure()
                fig.add_trace(
                    go.Image(
                        z=wordcloud.to_array(),
                    )
                )

                fig.update_layout(
                    title="Word Cloud Visualization",
                    showlegend=False,
                    margin=dict(t=30, b=0, l=0, r=0),
                )

                return fig

            return None

        except Exception as e:
            self.logger.error(f"Error creating wordcloud plot: {str(e)}")
            return None

    def create_advanced_report(
        self,
        compound: Union[str, Chem.Mol],
        compound_name: Optional[str] = None,
    ) -> Dict[str, go.Figure]:
        """Generate comprehensive visualization report.

        Args:
            compound: Input compound
            compound_name: Optional name for web data lookup

        Returns:
            Dictionary of visualization figures
        """
        try:
            report = {}

            # ML predictions
            report["toxicity"] = self.plot_toxicity_predictions(
                compound, include_confidence=True
            )

            report["abuse_potential"] = self.plot_abuse_potential(
                compound, include_mechanisms=True
            )

            report["psychopharm"] = self.plot_psychopharm_activity(
                compound, include_subtypes=True
            )

            report["nootropic"] = self.plot_nootropic_activity(
                compound, include_mechanisms=True
            )

            # Web data visualizations
            if compound_name:
                report["community_timeline"] = self.plot_community_data(
                    compound_name, plot_type="timeline"
                )
                report["community_network"] = self.plot_community_data(
                    compound_name, plot_type="network"
                )
                report["social_timeline"] = self.plot_social_data(
                    compound_name, plot_type="timeline"
                )
                report["social_wordcloud"] = self.plot_social_data(
                    compound_name, plot_type="wordcloud"
                )

            return report

        except Exception as e:
            self.logger.error(f"Error creating advanced report: {str(e)}")
            return {}
