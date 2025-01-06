"""Tests for Bluelight storage system."""

import pytest
from datetime import datetime, timedelta
from pathlib import Path
from unittest.mock import patch

from ..storage.bluelight_storage import BluelightStorage, BluelightAlert


@pytest.fixture
def storage(tmp_path):
    """Create storage instance with temporary directory."""
    return BluelightStorage(storage_dir=str(tmp_path))


@pytest.fixture
def sample_post():
    """Create sample post data."""
    return {
        "title": "Experience Report: Novel Compound XYZ",
        "content": "Detailed report of effects...",
        "author": "researcher123",
        "subforum": "advanced-drug-discussion",
        "url": "https://bluelight.org/xf/threads/123",
        "compounds": ["XYZ", "ABC"],
        "effects": ["stimulation", "focus"],
        "mechanisms": ["dopamine reuptake", "alpha2 agonism"],
        "safety": ["moderate addiction potential"],
        "screenshot": "data:image/png;base64,abc123",
        "javascript_logs": ["Loading complete", "No errors detected"],
    }


def test_store_post(storage, sample_post):
    """Test storing post data."""
    post_id = "123"
    storage.store_post(post_id, sample_post)

    # Verify post storage
    assert post_id in storage.posts
    stored = storage.posts[post_id]
    assert stored["title"] == sample_post["title"]
    assert stored["subforum"] == sample_post["subforum"]
    assert "stored_at" in stored

    # Verify compound storage
    for compound in sample_post["compounds"]:
        assert compound in storage.compounds
        assert post_id in storage.compounds[compound]["posts"]


def test_store_error(storage):
    """Test storing error information."""
    error = "Failed to scrape thread"
    metadata = {"url": "https://bluelight.org/xf/threads/123"}
    storage.store_error(error, metadata)

    assert len(storage.errors) == 1
    stored = storage.errors[0]
    assert stored["error"] == error
    assert stored["metadata"] == metadata
    assert "timestamp" in stored


@patch("binding_data_processor.processors.psychopharm.predictors.bbb.predictors.BBBPredictor")
def test_create_alert(mock_predictor, storage, sample_post):
    """Test creating safety alert."""
    # Mock BBB predictions
    mock_predictor.return_value.predict.return_value = 0.8
    storage.bbb_predictor = mock_predictor.return_value

    # Store sample post
    post_id = "123"
    storage.store_post(post_id, sample_post)

    # Create alert
    alert = storage.create_alert(
        post_id=post_id,
        compounds=sample_post["compounds"],
        safety_notes=["High BBB permeability"],
        severity="high",
        metadata={"source": "automated_analysis"},
    )

    # Verify alert
    assert isinstance(alert, BluelightAlert)
    assert alert.post_id == post_id
    assert alert.subforum == sample_post["subforum"]
    assert alert.severity == "high"
    assert len(alert.bbb_predictions) == len(sample_post["compounds"])
    assert alert.screenshot == sample_post["screenshot"]
    assert alert.javascript_logs == sample_post["javascript_logs"]

    # Verify alert storage
    assert len(storage.alerts) == 1
    stored = storage.alerts[0]
    assert stored["post_id"] == post_id


def test_get_alerts(storage, sample_post):
    """Test alert filtering."""
    # Store posts and create alerts
    posts = [
        {**sample_post, "subforum": "advanced-drug-discussion"},
        {**sample_post, "subforum": "psychedelic-drugs", "compounds": ["LSD"]},
    ]

    for i, post in enumerate(posts):
        post_id = str(i)
        storage.store_post(post_id, post)
        storage.create_alert(
            post_id=post_id,
            compounds=post["compounds"],
            safety_notes=["Test note"],
            severity="medium" if i == 0 else "high",
        )

    # Test filtering
    assert len(storage.get_alerts()) == 2
    assert len(storage.get_alerts(severity="high")) == 1
    assert len(storage.get_alerts(subforum="advanced-drug-discussion")) == 1
    assert len(storage.get_alerts(compound="XYZ")) == 1


def test_get_compound_stats(storage, sample_post):
    """Test compound statistics."""
    # Store multiple posts
    for i in range(3):
        post_id = str(i)
        storage.store_post(post_id, sample_post)
        if i < 2:  # Create alerts for first 2 posts
            storage.create_alert(
                post_id=post_id,
                compounds=sample_post["compounds"],
                safety_notes=["Test note"],
                severity="medium",
            )

    # Get stats for compound
    stats = storage.get_compound_stats("XYZ")
    assert stats["post_count"] == 3
    assert stats["alert_count"] == 2
    assert len(stats["subforums"]) == 1
    assert len(stats["latest_alerts"]) == 2


def test_get_trending_compounds(storage, sample_post):
    """Test trending compounds analysis."""
    # Store posts with different dates
    now = datetime.now()
    posts = [
        {**sample_post, "compounds": ["XYZ"]},
        {**sample_post, "compounds": ["ABC"]},
        {**sample_post, "compounds": ["XYZ"], "stored_at": (now - timedelta(days=10)).isoformat()},
    ]

    for i, post in enumerate(posts):
        storage.store_post(str(i), post)

    # Get trending (last 7 days)
    trending = storage.get_trending_compounds(days=7)
    assert len(trending) == 2
    assert trending[0]["name"] == "XYZ"  # Most recent posts
    assert trending[0]["recent_posts"] == 1


def test_get_safety_summary(storage, sample_post):
    """Test safety summary generation."""
    # Store posts and create alerts
    posts = [
        {**sample_post, "compounds": ["XYZ"]},
        {**sample_post, "compounds": ["ABC"], "subforum": "psychedelic-drugs"},
    ]

    for i, post in enumerate(posts):
        post_id = str(i)
        storage.store_post(post_id, post)
        storage.create_alert(
            post_id=post_id,
            compounds=post["compounds"],
            safety_notes=["Test note"],
            severity="high" if i == 0 else "low",
        )

    # Get summary
    summary = storage.get_safety_summary()
    assert summary["total_alerts"] == 2
    assert summary["unique_compounds"] == 2
    assert len(summary["affected_subforums"]) == 2
    assert summary["severity_distribution"] == {"high": 1, "medium": 0, "low": 1}


def test_get_error_summary(storage):
    """Test error summary generation."""
    # Store errors
    errors = [
        ("Scraping failed", {"url": "url1"}),
        ("Rate limited", {"url": "url2"}),
        ("Old error", {"url": "url3"}, datetime.now() - timedelta(days=10)),
    ]

    for error in errors:
        if len(error) == 3:
            err, meta, timestamp = error
            storage.errors.append(
                {
                    "error": err,
                    "metadata": meta,
                    "timestamp": timestamp.isoformat(),
                }
            )
        else:
            err, meta = error
            storage.store_error(err, meta)

    # Get summary (last 7 days)
    summary = storage.get_error_summary(days=7)
    assert summary["total_errors"] == 2
    assert len(summary["recent_errors"]) == 2
    assert len(summary["error_types"]) == 2


def test_file_persistence(tmp_path):
    """Test data persistence across instances."""
    storage1 = BluelightStorage(storage_dir=str(tmp_path))

    # Store data
    post_id = "123"
    post = {
        "title": "Test Post",
        "subforum": "test-forum",
        "compounds": ["XYZ"],
    }
    storage1.store_post(post_id, post)
    storage1.create_alert(
        post_id=post_id,
        compounds=["XYZ"],
        safety_notes=["Test"],
        severity="medium",
    )

    # Create new instance
    storage2 = BluelightStorage(storage_dir=str(tmp_path))

    # Verify data loaded
    assert post_id in storage2.posts
    assert len(storage2.alerts) == 1
    assert "XYZ" in storage2.compounds


def test_error_handling(storage):
    """Test error handling in storage operations."""
    # Test invalid post ID
    with pytest.raises(ValueError):
        storage.create_alert(
            post_id="invalid",
            compounds=["XYZ"],
            safety_notes=["Test"],
            severity="medium",
        )

    # Test invalid file operations
    storage.posts_file = Path("/invalid/path/file.json")
    storage.store_post("123", {"title": "Test"})  # Should log error but not raise
