"""Tests for Reddit storage system."""

import pytest
from datetime import datetime, timedelta
from unittest.mock import patch

from ..storage.reddit_storage import RedditStorage, RedditAlert


@pytest.fixture
def temp_storage_dir(tmp_path):
    """Create temporary storage directory."""
    return str(tmp_path / "reddit_test")


_BBB_MODULE = "binding_data_processor.processors.psychopharm"
BBB_PREDICTOR_PATH = f"{_BBB_MODULE}.predictors.bbb.predictors.BBBPredictor"


@pytest.fixture
def mock_bbb_predictor():
    """Create mock BBB predictor."""
    with patch(BBB_PREDICTOR_PATH) as mock:
        predictor = mock.return_value
        predictor.predict.return_value = 0.8
        yield predictor


@pytest.fixture
def storage(temp_storage_dir, mock_bbb_predictor):
    """Create storage instance with mocked BBB predictor."""
    return RedditStorage(storage_dir=temp_storage_dir)


@pytest.fixture
def sample_post():
    """Create sample post data."""
    return {
        "title": "Test Compound Analysis",
        "content": "Analysis of test compound shows interesting effects",
        "author": "test_user",
        "subreddit": "DrugNerds",
        "score": 100,
        "num_comments": 10,
        "url": "https://reddit.com/r/DrugNerds/test",
        "compounds": ["TestCompound", "AnotherCompound"],
    }


def test_store_post(storage, sample_post):
    """Test storing post data."""
    post_id = "test123"
    storage.store_post(post_id, sample_post)

    # Verify post storage
    assert post_id in storage.posts
    stored = storage.posts[post_id]
    assert stored["title"] == sample_post["title"]
    assert "stored_at" in stored

    # Verify compound storage
    for compound in sample_post["compounds"]:
        assert compound in storage.compounds
        assert post_id in storage.compounds[compound]["posts"]
        assert "first_seen" in storage.compounds[compound]


def test_create_alert(storage, sample_post, mock_bbb_predictor):
    """Test creating safety alert."""
    post_id = "test123"
    storage.store_post(post_id, sample_post)

    alert = storage.create_alert(
        post_id=post_id,
        compounds=["TestCompound"],
        safety_notes=["Safety concern"],
        severity="high",
    )

    assert isinstance(alert, RedditAlert)
    assert alert.post_id == post_id
    assert alert.compounds == ["TestCompound"]
    assert alert.severity == "high"
    assert alert.bbb_predictions["TestCompound"] == 0.8

    # Verify alert storage
    assert len(storage.alerts) == 1
    stored = storage.alerts[0]
    assert stored["post_id"] == post_id


def test_get_alerts_filtering(storage, sample_post):
    """Test alert filtering."""
    # Create test alerts
    post_id = "test123"
    storage.store_post(post_id, sample_post)

    alerts = [
        {
            "compounds": ["CompoundA"],
            "safety_notes": ["Note A"],
            "severity": "high",
        },
        {
            "compounds": ["CompoundB"],
            "safety_notes": ["Note B"],
            "severity": "medium",
        },
        {
            "compounds": ["CompoundC"],
            "safety_notes": ["Note C"],
            "severity": "low",
        },
    ]

    created = []
    for i, data in enumerate(alerts):
        alert = storage.create_alert(
            post_id=post_id,
            **data,
        )
        created.append(alert)

    # Test filters
    high_alerts = storage.get_alerts(severity="high")
    assert len(high_alerts) == 1
    assert high_alerts[0].severity == "high"

    compound_alerts = storage.get_alerts(compound="CompoundB")
    assert len(compound_alerts) == 1
    assert "CompoundB" in compound_alerts[0].compounds

    subreddit_alerts = storage.get_alerts(subreddit="DrugNerds")
    assert len(subreddit_alerts) == 3


def test_get_compound_stats(storage, sample_post):
    """Test compound statistics."""
    post_id = "test123"
    storage.store_post(post_id, sample_post)

    storage.create_alert(
        post_id=post_id,
        compounds=["TestCompound"],
        safety_notes=["Note"],
        severity="medium",
    )

    stats = storage.get_compound_stats("TestCompound")
    assert stats["name"] == "TestCompound"
    assert stats["post_count"] == 1
    assert stats["alert_count"] == 1
    assert stats["subreddits"] == ["DrugNerds"]
    assert stats["bbb_prediction"] == 0.8
    assert len(stats["latest_alerts"]) == 1


def test_get_trending_compounds(storage):
    """Test trending compounds analysis."""
    # Create posts at different times
    posts = [
        {
            "id": "old_post",
            "data": {
                "title": "Old Post",
                "subreddit": "DrugNerds",
                "compounds": ["OldCompound"],
                "stored_at": (datetime.now() - timedelta(days=10)).isoformat(),
            },
        },
        {
            "id": "new_post",
            "data": {
                "title": "New Post",
                "subreddit": "DrugNerds",
                "compounds": ["TrendingCompound"],
                "stored_at": datetime.now().isoformat(),
            },
        },
    ]

    for post in posts:
        storage.store_post(post["id"], post["data"])

    trending = storage.get_trending_compounds(days=7)
    assert len(trending) == 1
    assert trending[0]["name"] == "TrendingCompound"
    assert trending[0]["recent_posts"] == 1


def test_get_safety_summary(storage, sample_post):
    """Test safety summary generation."""
    post_id = "test123"
    storage.store_post(post_id, sample_post)

    # Create alerts with different severities
    alerts = [
        ("CompoundA", "high"),
        ("CompoundB", "medium"),
        ("CompoundC", "low"),
    ]

    for compound, severity in alerts:
        storage.create_alert(
            post_id=post_id,
            compounds=[compound],
            safety_notes=[f"Note for {compound}"],
            severity=severity,
        )

    summary = storage.get_safety_summary()
    assert summary["total_alerts"] == 3
    assert summary["unique_compounds"] == 3
    assert len(summary["affected_subreddits"]) == 1
    assert summary["severity_distribution"] == {"high": 1, "medium": 1, "low": 1}
    assert len(summary["high_bbb_compounds"]) == 3  # All have BBB prediction of 0.8


def test_error_handling(storage):
    """Test error handling."""
    # Test invalid post ID
    with pytest.raises(ValueError, match="Post .* not found"):
        storage.create_alert(
            post_id="invalid",
            compounds=["Test"],
            safety_notes=["Note"],
        )

    # Test BBB prediction error
    post_id = "test123"
    storage.store_post(post_id, {"title": "Test", "subreddit": "Test"})

    storage.bbb_predictor.predict.side_effect = Exception("Prediction failed")
    alert = storage.create_alert(
        post_id=post_id,
        compounds=["ErrorCompound"],
        safety_notes=["Note"],
    )
    assert not alert.bbb_predictions


def test_file_persistence(temp_storage_dir, sample_post):
    """Test data persistence across instances."""
    # First instance
    storage1 = RedditStorage(storage_dir=temp_storage_dir)
    post_id = "test123"
    storage1.store_post(post_id, sample_post)
    storage1.create_alert(
        post_id=post_id,
        compounds=["TestCompound"],
        safety_notes=["Note"],
    )

    # Second instance should load existing data
    storage2 = RedditStorage(storage_dir=temp_storage_dir)
    assert post_id in storage2.posts
    assert "TestCompound" in storage2.compounds
    assert len(storage2.alerts) == 1
