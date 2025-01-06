"""Integration tests for Bluelight client and storage."""

import pytest
from datetime import datetime, timedelta
from unittest.mock import patch, Mock

from ..clients.bluelight import BluelightCrawl4AIClient
from ..storage.bluelight_storage import BluelightStorage


@pytest.fixture
def storage(tmp_path):
    """Create storage instance with temporary directory."""
    return BluelightStorage(storage_dir=str(tmp_path))


@pytest.fixture
def client():
    """Create Bluelight client instance."""
    return BluelightCrawl4AIClient()


@pytest.fixture
def mock_crawler():
    """Create mock crawler with sample response."""
    mock = Mock()
    mock.arun.return_value = Mock(
        success=True,
        extracted_content=[
            {
                "title": ["<a href='/threads/123'>Experience Report: XYZ</a>"],
                "content": ["Detailed report of effects..."],
                "author": ["researcher123"],
                "subforum": ["advanced-drug-discussion"],
                "date": ["2023-12-01"],
                "compounds": ["XYZ", "ABC"],
                "effects": ["stimulation", "focus"],
                "mechanisms": ["dopamine reuptake", "alpha2 agonism"],
                "safety": ["moderate addiction potential"],
            }
        ],
        screenshot="data:image/png;base64,abc123",
        javascript_logs=["Loading complete", "No errors detected"],
    )
    return mock


@pytest.fixture
def mock_bbb_predictor():
    """Create mock BBB predictor."""
    mock = Mock()
    mock.predict.return_value = 0.8
    return mock


@patch("crawl4ai.AsyncWebCrawler")
async def test_scrape_and_store(mock_crawler_class, mock_crawler, storage, client):
    """Test scraping posts and storing them."""
    # Setup mocks
    mock_crawler_class.return_value.__aenter__.return_value = mock_crawler

    # Scrape posts
    posts = await client.search_posts(
        query="XYZ",
        subforums=["advanced-drug-discussion"],
        max_results=1,
    )

    assert len(posts) == 1
    post = posts[0]

    # Store post
    post_id = post.url.split("/")[-1]
    storage.store_post(post_id, post)

    # Verify storage
    assert post_id in storage.posts
    stored = storage.posts[post_id]
    assert stored["title"] == "Experience Report: XYZ"
    assert stored["subforum"] == "advanced-drug-discussion"
    assert stored["compounds"] == ["XYZ", "ABC"]
    assert stored["screenshot"] == "data:image/png;base64,abc123"
    assert stored["javascript_logs"] == ["Loading complete", "No errors detected"]

    # Verify compounds stored
    for compound in post.compounds:
        assert compound in storage.compounds
        assert post_id in storage.compounds[compound]["posts"]


@patch("crawl4ai.AsyncWebCrawler")
async def test_error_handling(mock_crawler_class, storage, client):
    """Test error handling and storage."""
    # Setup mock to fail
    mock = Mock()
    mock.arun.side_effect = Exception("Rate limited")
    mock_crawler_class.return_value.__aenter__.return_value = mock

    # Attempt scraping
    posts = await client.search_posts(
        query="XYZ",
        subforums=["advanced-drug-discussion"],
    )

    assert len(posts) == 0

    # Verify error stored
    assert len(storage.errors) == 1
    error = storage.errors[0]
    assert "Rate limited" in error["error"]
    assert error["metadata"]["query"] == "XYZ"


@patch("crawl4ai.AsyncWebCrawler")
async def test_bbb_prediction_caching(
    mock_crawler_class,
    mock_crawler,
    mock_bbb_predictor,
    storage,
    client,
):
    """Test BBB prediction caching."""
    # Setup mocks
    mock_crawler_class.return_value.__aenter__.return_value = mock_crawler
    storage.bbb_predictor = mock_bbb_predictor

    # Scrape and store post
    posts = await client.search_posts(query="XYZ", max_results=1)
    post = posts[0]
    post_id = post.url.split("/")[-1]
    storage.store_post(post_id, post)

    # Create first alert
    alert1 = storage.create_alert(
        post_id=post_id,
        compounds=post.compounds,
        safety_notes=["Test"],
        severity="medium",
    )

    # Verify BBB prediction made
    assert mock_bbb_predictor.predict.call_count == len(post.compounds)
    for compound in post.compounds:
        assert alert1.bbb_predictions[compound] == 0.8

    # Reset mock
    mock_bbb_predictor.predict.reset_mock()

    # Create second alert
    alert2 = storage.create_alert(
        post_id=post_id,
        compounds=post.compounds,
        safety_notes=["Test 2"],
        severity="high",
    )

    # Verify cached predictions used
    assert mock_bbb_predictor.predict.call_count == 0
    for compound in post.compounds:
        assert alert2.bbb_predictions[compound] == 0.8


@patch("crawl4ai.AsyncWebCrawler")
async def test_screenshot_handling(mock_crawler_class, mock_crawler, storage, client):
    """Test screenshot and JavaScript log handling."""
    # Setup mock with multiple screenshots
    mock_crawler_class.return_value.__aenter__.return_value = mock_crawler
    mock_crawler.arun.return_value.screenshot = "data:image/png;base64,xyz789"
    mock_crawler.arun.return_value.javascript_logs = [
        "New error detected",
        "Recovery successful",
    ]

    # Scrape and store post
    posts = await client.search_posts(query="XYZ", max_results=1)
    post = posts[0]
    post_id = post.url.split("/")[-1]
    storage.store_post(post_id, post)

    # Create alert
    alert = storage.create_alert(
        post_id=post_id,
        compounds=post.compounds,
        safety_notes=["Test"],
        severity="medium",
    )

    # Verify screenshot and logs stored
    assert alert.screenshot == "data:image/png;base64,xyz789"
    assert len(alert.javascript_logs) == 2
    assert "New error detected" in alert.javascript_logs
    assert "Recovery successful" in alert.javascript_logs


@patch("crawl4ai.AsyncWebCrawler")
async def test_trending_analysis(mock_crawler_class, mock_crawler, storage, client):
    """Test trending compound analysis with stored posts."""
    # Setup mock with different dates
    mock_crawler_class.return_value.__aenter__.return_value = mock_crawler

    # Store posts with different dates
    posts = await client.search_posts(query="XYZ", max_results=1)
    post = posts[0]

    # Recent post
    post_id1 = "123"
    storage.store_post(post_id1, post)

    # Old post
    post_id2 = "456"
    old_post = dict(post)
    old_post["stored_at"] = (datetime.now() - timedelta(days=10)).isoformat()
    storage.store_post(post_id2, old_post)

    # Get trending compounds
    trending = storage.get_trending_compounds(days=7)

    # Verify trending analysis
    assert len(trending) > 0
    assert trending[0]["recent_posts"] == 1  # Only recent post counted
    assert trending[0]["total_posts"] == 2  # Both posts counted
