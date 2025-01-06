"""Enhanced dashboard component.

This component provides:
1. Overview statistics
2. Recent activity
3. Data quality metrics
4. System status
5. Quick actions
"""

from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional, Any

import streamlit as st
import plotly.express as px
import plotly.graph_objects as go

from ...models.compound.enhanced import EnhancedCompound


@dataclass
class DashboardConfig:
    """Configuration for dashboard component."""

    # Display options
    show_overview: bool = True
    show_recent: bool = True
    show_quality: bool = True
    show_status: bool = True
    show_actions: bool = True

    # Data options
    max_recent: int = 10
    quality_threshold: float = 0.8
    update_interval: int = 60  # seconds

    # Chart options
    chart_height: int = 300
    chart_width: int = 400


class CompoundDashboardEnhanced:
    """Enhanced dashboard component."""

    def __init__(self, config: Optional[DashboardConfig] = None):
        """Initialize dashboard.

        Args:
            config: Optional dashboard configuration
        """
        self.config = config or DashboardConfig()
        self.last_update = datetime.now()

        # Initialize state
        if "compounds" not in st.session_state:
            st.session_state.compounds = []
        if "quality_metrics" not in st.session_state:
            st.session_state.quality_metrics = {}
        if "system_status" not in st.session_state:
            st.session_state.system_status = {}

    def update_compounds(self, compounds: List[EnhancedCompound]) -> None:
        """Update compound list.

        Args:
            compounds: List of compounds
        """
        st.session_state.compounds = compounds
        self.last_update = datetime.now()

        # Update quality metrics
        self._update_quality_metrics()

        # Update system status
        self._update_system_status()

    def _update_quality_metrics(self) -> None:
        """Update data quality metrics."""
        metrics = {
            "completeness": self._calculate_completeness(),
            "accuracy": self._calculate_accuracy(),
            "consistency": self._calculate_consistency(),
            "timeliness": self._calculate_timeliness(),
        }
        st.session_state.quality_metrics = metrics

    def _update_system_status(self) -> None:
        """Update system status."""
        status = {
            "compounds": len(st.session_state.compounds),
            "last_update": self.last_update,
            "quality_score": self._calculate_quality_score(),
            "processing_rate": self._calculate_processing_rate(),
            "error_rate": self._calculate_error_rate(),
        }
        st.session_state.system_status = status

    def _calculate_completeness(self) -> float:
        """Calculate data completeness score."""
        if not st.session_state.compounds:
            return 0.0

        total_score = 0.0
        for compound in st.session_state.compounds:
            # Check required fields
            score = 0.0
            if compound.name:
                score += 0.2
            if compound.smiles:
                score += 0.2
            if compound.cas:
                score += 0.2

            # Check enhanced data
            if compound.predictions.binding_affinities:
                score += 0.2
            if compound.web_data.community_reports:
                score += 0.2

            total_score += score

        return total_score / len(st.session_state.compounds)

    def _calculate_accuracy(self) -> float:
        """Calculate data accuracy score."""
        if not st.session_state.compounds:
            return 0.0

        total_score = 0.0
        for compound in st.session_state.compounds:
            # Check prediction confidence
            confidence_scores = compound.predictions.confidence_scores
            if confidence_scores:
                total_score += sum(confidence_scores.values()) / len(confidence_scores)

        return total_score / len(st.session_state.compounds)

    def _calculate_consistency(self) -> float:
        """Calculate data consistency score."""
        if not st.session_state.compounds:
            return 0.0

        total_score = 0.0
        for compound in st.session_state.compounds:
            # Check for data inconsistencies
            score = 1.0

            # Check structure consistency
            if compound.smiles and compound.inchi:
                if not compound.mol:
                    score -= 0.5

            # Check prediction consistency
            if compound.predictions.binding_affinities:
                for target, value in compound.predictions.binding_affinities.items():
                    if value < 0:
                        score -= 0.1

            total_score += max(0.0, score)

        return total_score / len(st.session_state.compounds)

    def _calculate_timeliness(self) -> float:
        """Calculate data timeliness score."""
        if not st.session_state.compounds:
            return 0.0

        total_score = 0.0
        now = datetime.now()

        for compound in st.session_state.compounds:
            # Check data freshness
            score = 1.0

            # Check prediction age
            pred_age = (now - compound.predictions.prediction_date).days
            if pred_age > 30:  # Older than 30 days
                score -= 0.3

            # Check web data age
            web_age = (now - compound.web_data.last_updated).days
            if web_age > 7:  # Older than 7 days
                score -= 0.3

            # Check analysis age
            analysis_age = (now - compound.analysis.analysis_date).days
            if analysis_age > 14:  # Older than 14 days
                score -= 0.3

            total_score += max(0.0, score)

        return total_score / len(st.session_state.compounds)

    def _calculate_quality_score(self) -> float:
        """Calculate overall quality score."""
        if not st.session_state.quality_metrics:
            return 0.0

        weights = {
            "completeness": 0.3,
            "accuracy": 0.3,
            "consistency": 0.2,
            "timeliness": 0.2,
        }

        score = 0.0
        for metric, weight in weights.items():
            score += st.session_state.quality_metrics[metric] * weight

        return score

    def _calculate_processing_rate(self) -> float:
        """Calculate compound processing rate (compounds/minute)."""
        if not hasattr(self, "_last_count"):
            self._last_count = 0
            self._last_time = datetime.now()
            return 0.0

        current_count = len(st.session_state.compounds)
        current_time = datetime.now()

        time_diff = (current_time - self._last_time).total_seconds() / 60
        if time_diff == 0:
            return 0.0

        rate = (current_count - self._last_count) / time_diff

        self._last_count = current_count
        self._last_time = current_time

        return rate

    def _calculate_error_rate(self) -> float:
        """Calculate error rate (errors/minute)."""
        if not hasattr(self, "_error_count"):
            self._error_count = 0
            self._error_time = datetime.now()
            return 0.0

        time_diff = (datetime.now() - self._error_time).total_seconds() / 60
        if time_diff == 0:
            return 0.0

        return self._error_count / time_diff

    def render_overview(self) -> None:
        """Render overview section."""
        if not self.config.show_overview:
            return

        st.subheader("Overview")

        # Show key metrics
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric(
                "Total Compounds",
                len(st.session_state.compounds),
            )

        with col2:
            st.metric(
                "Quality Score",
                f"{self._calculate_quality_score():.1%}",
            )

        with col3:
            st.metric(
                "Processing Rate",
                f"{self._calculate_processing_rate():.1f}/min",
            )

        with col4:
            st.metric(
                "Error Rate",
                f"{self._calculate_error_rate():.2f}/min",
            )

    def render_recent(self) -> None:
        """Render recent activity section."""
        if not self.config.show_recent or not st.session_state.compounds:
            return

        st.subheader("Recent Activity")

        # Get recent compounds
        recent = sorted(
            st.session_state.compounds,
            key=lambda x: x.last_updated,
            reverse=True,
        )[: self.config.max_recent]

        # Show recent compounds
        for compound in recent:
            with st.expander(compound.name):
                st.write(f"Updated: {compound.last_updated}")
                st.write(f"SMILES: {compound.smiles}")
                if compound.predictions.binding_affinities:
                    st.write("Binding Predictions:")
                    for target, value in compound.predictions.binding_affinities.items():
                        st.write(f"- {target}: {value:.1f} nM")

    def render_quality(self) -> None:
        """Render data quality section."""
        if not self.config.show_quality:
            return

        st.subheader("Data Quality")

        # Create quality metrics chart
        metrics = st.session_state.quality_metrics
        if not metrics:
            return

        fig = go.Figure()

        # Add radar chart
        fig.add_trace(
            go.Scatterpolar(
                r=[
                    metrics["completeness"],
                    metrics["accuracy"],
                    metrics["consistency"],
                    metrics["timeliness"],
                ],
                theta=[
                    "Completeness",
                    "Accuracy",
                    "Consistency",
                    "Timeliness",
                ],
                fill="toself",
            )
        )

        # Update layout
        fig.update_layout(
            polar=dict(
                radialaxis=dict(
                    visible=True,
                    range=[0, 1],
                ),
            ),
            showlegend=False,
            width=self.config.chart_width,
            height=self.config.chart_height,
        )

        st.plotly_chart(fig)

    def render_status(self) -> None:
        """Render system status section."""
        if not self.config.show_status:
            return

        st.subheader("System Status")

        status = st.session_state.system_status
        if not status:
            return

        # Show status metrics
        col1, col2 = st.columns(2)

        with col1:
            st.write("Last Update:", status["last_update"])
            st.write("Quality Score:", f"{status['quality_score']:.1%}")

        with col2:
            st.write("Processing Rate:", f"{status['processing_rate']:.1f}/min")
            st.write("Error Rate:", f"{status['error_rate']:.2f}/min")

    def render_actions(self) -> None:
        """Render quick actions section."""
        if not self.config.show_actions:
            return

        st.subheader("Quick Actions")

        col1, col2 = st.columns(2)

        with col1:
            if st.button("Update Data"):
                self._update_quality_metrics()
                self._update_system_status()
                st.success("Data updated")

            if st.button("Export Report"):
                # TODO: Implement report export
                st.info("Report export not implemented")

        with col2:
            if st.button("Clear Cache"):
                # TODO: Implement cache clearing
                st.info("Cache clearing not implemented")

            if st.button("System Check"):
                # TODO: Implement system check
                st.info("System check not implemented")

    def render(self) -> None:
        """Render dashboard."""
        st.title("ChemData Dashboard")

        # Render sections
        self.render_overview()
        self.render_recent()
        self.render_quality()
        self.render_status()
        self.render_actions()
