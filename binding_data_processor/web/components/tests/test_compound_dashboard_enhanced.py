"""Tests for enhanced compound dashboard component."""

import pytest
from pathlib import Path
from unittest.mock import Mock, patch
from datetime import datetime, timedelta

from ...components.compound_dashboard_enhanced import CompoundDashboardEnhanced
from ....models.compound import Compound


@pytest.fixture
def dashboard(tmp_path):
    """Create test dashboard."""
    return CompoundDashboardEnhanced(
        template_dir=tmp_path / "templates",
        static_dir=tmp_path / "static",
    )


@pytest.fixture
def compounds():
    """Create test compounds."""
    compounds = []

    # Compound with binding data
    compound1 = Compound(
        name="Test Compound 1",
        smiles="CC1=CC=CC=C1",
        cas_number="100-00-0",
    )
    compound1.binding_data = [
        {
            "target": "5-HT2A",
            "affinity": "1.2",
            "confidence": 0.9,
        },
        {
            "target": "D2",
            "affinity": "2.5",
            "confidence": 0.8,
        },
    ]
    compounds.append(compound1)

    # Compound with predictions
    compound2 = Compound(
        name="Test Compound 2",
        smiles="CC2=CC=CC=C2",
        cas_number="200-00-0",
    )
    compound2.predictions = {
        "activity": {
            "stimulant": 0.7,
            "psychedelic": 0.3,
            "sedative": 0.2,
        },
        "safety": {
            "toxicity": 0.2,
            "addiction": 0.4,
            "overdose": 0.3,
        },
    }
    compounds.append(compound2)

    # Compound with social data
    now = datetime.now()
    compound3 = Compound(
        name="Test Compound 3",
        smiles="CC3=CC=CC=C3",
        cas_number="300-00-0",
    )
    compound3.social_data = {
        "reddit": {
            "posts": [
                {
                    "id": "1",
                    "title": "Test Post 1",
                    "created_utc": (now - timedelta(days=2)).isoformat(),
                },
                {
                    "id": "2",
                    "title": "Test Post 2",
                    "created_utc": (now - timedelta(days=1)).isoformat(),
                },
            ],
        },
        "twitter": {
            "tweets": [
                {
                    "id": "1",
                    "text": "Test Tweet 1",
                    "created_at": (now - timedelta(days=2)).isoformat(),
                },
                {
                    "id": "2",
                    "text": "Test Tweet 2",
                    "created_at": (now - timedelta(days=1)).isoformat(),
                },
            ],
        },
    }
    compounds.append(compound3)

    return compounds


@patch("flask.render_template")
def test_render_dashboard_list(mock_render, dashboard, compounds):
    """Test rendering dashboard list view."""
    # Mock component results
    with patch.multiple(
        dashboard,
        list_view=Mock(),
        search=Mock(),
        export=Mock(),
    ):
        # Configure mocks
        dashboard.list_view.render_list.return_value.success = True
        dashboard.list_view.render_list.return_value.data = {"list": "data"}
        dashboard.search.render_search.return_value.success = True
        dashboard.search.render_search.return_value.data = {
            "compounds": compounds,
            "search": "data",
        }
        dashboard.export.render_export.return_value.success = True
        dashboard.export.render_export.return_value.data = {"export": "data"}
        mock_render.return_value = "<html>Test</html>"

        # Render dashboard
        result = dashboard.render_dashboard(
            compounds=compounds,
            view="list",
            query="test",
            filters={"class": "stimulant"},
            export_format="tsv",
            export_columns=["name", "smiles"],
        )

        # Check result
        assert result.success
        assert result.data["html"] == "<html>Test</html>"
        assert result.data["compounds"] == compounds
        assert result.data["view"] == "list"
        assert "components" in result.data

        # Check components
        components = result.data["components"]
        assert "search" in components
        assert "list" in components
        assert "export" in components

        # Check stats
        stats = dashboard.dashboard_stats
        assert stats["total_views"] == 1
        assert stats["list_views"] == 1
        assert stats["searches"] == 1
        assert stats["exports"] == 1


@patch("flask.render_template")
def test_render_dashboard_detail(mock_render, dashboard, compounds):
    """Test rendering dashboard detail view."""
    # Mock component results
    with patch.multiple(
        dashboard,
        detail_view=Mock(),
        visualization=Mock(),
        analysis=Mock(),
    ):
        # Configure mocks
        dashboard.detail_view.render_detail.return_value.success = True
        dashboard.detail_view.render_detail.return_value.data = {"detail": "data"}
        dashboard.visualization.render_visualization.return_value.success = True
        dashboard.visualization.render_visualization.return_value.data = {"viz": "data"}
        dashboard.analysis.render_analysis.return_value.success = True
        dashboard.analysis.render_analysis.return_value.data = {"analysis": "data"}
        mock_render.return_value = "<html>Test</html>"

        # Render dashboard
        result = dashboard.render_dashboard(
            compounds=compounds,
            selected_compound=compounds[0],
            view="detail",
            analysis_types=["binding", "activity"],
        )

        # Check result
        assert result.success
        assert result.data["html"] == "<html>Test</html>"
        assert result.data["compounds"] == compounds
        assert result.data["view"] == "detail"
        assert result.data["selected_compound"] == compounds[0]

        # Check components
        components = result.data["components"]
        assert "detail" in components
        assert "visualization" in components
        assert "analysis" in components

        # Check stats
        stats = dashboard.dashboard_stats
        assert stats["total_views"] == 1
        assert stats["detail_views"] == 1
        assert stats["analyses"] == 1


def test_component_failures(dashboard, compounds):
    """Test handling of component failures."""
    # Test list view failure
    with patch.object(dashboard.list_view, "render_list") as mock_list:
        mock_list.return_value.success = False
        mock_list.return_value.error = "List error"
        result = dashboard.render_dashboard(compounds, view="list")
        assert not result.success
        assert result.error == "List error"

    # Test detail view failure
    with patch.object(dashboard.detail_view, "render_detail") as mock_detail:
        mock_detail.return_value.success = False
        mock_detail.return_value.error = "Detail error"
        result = dashboard.render_dashboard(
            compounds,
            selected_compound=compounds[0],
            view="detail",
        )
        assert not result.success
        assert result.error == "Detail error"

    # Test search failure
    with patch.object(dashboard.search, "render_search") as mock_search:
        mock_search.return_value.success = False
        mock_search.return_value.error = "Search error"
        result = dashboard.render_dashboard(compounds, query="test")
        assert not result.success
        assert result.error == "Search error"


def test_view_history(dashboard, compounds):
    """Test view history tracking."""
    # View list
    dashboard.render_dashboard(compounds, view="list")

    # View detail
    dashboard.render_dashboard(
        compounds,
        selected_compound=compounds[0],
        view="detail",
    )

    # Search
    dashboard.render_dashboard(
        compounds,
        query="test",
        filters={"class": "stimulant"},
    )

    # Check history
    history = dashboard.dashboard_stats["view_history"]
    assert len(history) == 3
    assert history[0]["view"] == "list"
    assert history[1]["view"] == "detail"
    assert history[1]["compound"] == compounds[0].name
    assert history[2]["query"] == "test"
    assert history[2]["filters"] == {"class": "stimulant"}


def test_metrics(dashboard, compounds):
    """Test metrics collection."""
    # Mock component metrics
    with patch.multiple(
        dashboard,
        list_view=Mock(get_metrics=Mock(return_value={"list": "metrics"})),
        detail_view=Mock(get_metrics=Mock(return_value={"detail": "metrics"})),
        search=Mock(get_metrics=Mock(return_value={"search": "metrics"})),
        visualization=Mock(get_metrics=Mock(return_value={"viz": "metrics"})),
        analysis=Mock(get_metrics=Mock(return_value={"analysis": "metrics"})),
        export=Mock(get_metrics=Mock(return_value={"export": "metrics"})),
    ):
        # Get metrics
        metrics = dashboard.get_metrics()

        # Check metrics
        assert "dashboard_stats" in metrics
        assert metrics["list_stats"] == {"list": "metrics"}
        assert metrics["detail_stats"] == {"detail": "metrics"}
        assert metrics["search_stats"] == {"search": "metrics"}
        assert metrics["visualization_stats"] == {"viz": "metrics"}
        assert metrics["analysis_stats"] == {"analysis": "metrics"}
        assert metrics["export_stats"] == {"export": "metrics"}
