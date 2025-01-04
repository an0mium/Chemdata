"""Tests for social media monitoring functionality."""

import pytest
from unittest.mock import patch, Mock

from ..base import PsychoactiveClass
from ..compound import PsychoactiveCompound
from ..social_monitoring import SocialMonitoringPipeline


@pytest.fixture
def test_compound():
    """Create a test compound fixture."""
    compound = PsychoactiveCompound(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    compound.psychoactive_class = PsychoactiveClass.STIMULANT
    return compound


@pytest.fixture
def monitoring_pipeline():
    """Create a test social monitoring pipeline fixture."""
    return SocialMonitoringPipeline()


class TestSocialMonitoringPipeline:
    """Tests for SocialMonitoringPipeline class."""

    def test_initialization(self, monitoring_pipeline):
        """Test initialization of SocialMonitoringPipeline."""
        assert monitoring_pipeline.reddit_client is not None
        assert monitoring_pipeline.twitter_client is not None
        assert monitoring_pipeline.stats == {}

    @patch('praw.Reddit')
    def test_reddit_monitoring(self, mock_reddit, test_compound, monitoring_pipeline):
        """Test Reddit data monitoring."""
        # Mock subreddit search results
        mock_submission = Mock()
        mock_submission.title = "Experience with Caffeine"
        mock_submission.selftext = "Great for focus but caused some anxiety"
        mock_submission.created_utc = 1672531200  # 2023-01-01
        mock_submission.score = 50
        mock_submission.num_comments = 10
        
        mock_subreddit = Mock()
        mock_subreddit.search.return_value = [mock_submission]
        
        mock_reddit.return_value.subreddit.return_value = mock_subreddit
        
        # Monitor Reddit
        data = monitoring_pipeline.monitor_reddit(
            test_compound,
            subreddits=["nootropics", "researchchemicals"]
        )
        
        # Check monitored data
        assert len(data["mentions"]) > 0
        assert "focus" in str(data["mentions"][0]["text"]).lower()
        assert "anxiety" in str(data["mentions"][0]["text"]).lower()
        assert data["mentions"][0]["engagement"]["score"] == 50
        assert data["mentions"][0]["engagement"]["comments"] == 10

    @patch('tweepy.Client')
    def test_twitter_monitoring(self, mock_twitter, test_compound, monitoring_pipeline):
        """Test Twitter data monitoring."""
        # Mock tweet search results
        mock_tweet = {
            'id': '1234567890',
            'text': 'Caffeine really helps with productivity #nootropics',
            'created_at': '2023-01-01T00:00:00Z',
            'public_metrics': {
                'retweet_count': 5,
                'like_count': 20,
                'reply_count': 3,
            },
        }
        
        mock_twitter.return_value.search_recent_tweets.return_value.data = [mock_tweet]
        
        # Monitor Twitter
        data = monitoring_pipeline.monitor_twitter(test_compound)
        
        # Check monitored data
        assert len(data["mentions"]) > 0
        assert "productivity" in data["mentions"][0]["text"].lower()
        assert data["mentions"][0]["engagement"]["retweets"] == 5
        assert data["mentions"][0]["engagement"]["likes"] == 20

    def test_sentiment_analysis(self, test_compound, monitoring_pipeline):
        """Test sentiment analysis of social media mentions."""
        # Test data
        mentions = [{
            "text": "Caffeine is great for focus and productivity",
            "source": "reddit",
            "timestamp": "2023-01-01T00:00:00Z",
        }, {
            "text": "Too much caffeine causes anxiety and jitters",
            "source": "twitter",
            "timestamp": "2023-01-01T00:00:00Z",
        }]
        
        # Analyze sentiment
        analysis = monitoring_pipeline.analyze_sentiment(mentions)
        
        # Check analysis results
        assert "positive" in analysis
        assert "negative" in analysis
        assert analysis["positive"]["focus"] > 0.5
        assert analysis["negative"]["anxiety"] > 0.5

    def test_trend_detection(self, test_compound, monitoring_pipeline):
        """Test trend detection in social media mentions."""
        # Test data over time
        historical_data = {
            "2023-01": {"mentions": 100, "sentiment": 0.6},
            "2023-02": {"mentions": 150, "sentiment": 0.7},
            "2023-03": {"mentions": 200, "sentiment": 0.8},
        }
        
        # Detect trends
        trends = monitoring_pipeline.detect_trends(historical_data)
        
        # Check trend analysis
        assert "mention_trend" in trends
        assert "sentiment_trend" in trends
        assert trends["mention_trend"]["direction"] == "increasing"
        assert trends["sentiment_trend"]["direction"] == "improving"

    @patch('praw.Reddit')
    @patch('tweepy.Client')
    def test_full_monitoring(self, mock_twitter, mock_reddit, test_compound, monitoring_pipeline):
        """Test full social media monitoring pipeline."""
        # Mock Reddit data
        mock_submission = Mock()
        mock_submission.title = "Caffeine Experience"
        mock_submission.selftext = "Great for focus"
        mock_submission.created_utc = 1672531200
        
        mock_subreddit = Mock()
        mock_subreddit.search.return_value = [mock_submission]
        mock_reddit.return_value.subreddit.return_value = mock_subreddit
        
        # Mock Twitter data
        mock_tweet = {
            'text': 'Caffeine helps productivity',
            'created_at': '2023-01-01T00:00:00Z',
            'public_metrics': {'like_count': 10},
        }
        mock_twitter.return_value.search_recent_tweets.return_value.data = [mock_tweet]
        
        # Run full monitoring
        enriched_compound = monitoring_pipeline.monitor_compound(test_compound)
        
        # Check monitoring results
        assert len(enriched_compound.social_data["reddit"]) > 0
        assert len(enriched_compound.social_data["twitter"]) > 0
        assert enriched_compound.social_data["trends"] is not None
        assert enriched_compound.last_monitored is not None

    def test_error_handling(self, test_compound, monitoring_pipeline):
        """Test error handling during monitoring."""
        # Mock failed API calls
        with patch('praw.Reddit', side_effect=Exception("API Error")):
            # Attempt monitoring
            data = monitoring_pipeline.monitor_reddit(test_compound)
            
            # Check error handling
            assert data == {}
            assert "reddit_error" in monitoring_pipeline.stats
            assert monitoring_pipeline.stats["failed_requests"] > 0

    def test_rate_limiting(self, monitoring_pipeline):
        """Test rate limiting functionality."""
        # Check rate limit configuration
        assert monitoring_pipeline.rate_limits["reddit"]["requests_per_minute"] > 0
        assert monitoring_pipeline.rate_limits["twitter"]["requests_per_minute"] > 0
        
        # Test rate limiting
        with patch('time.sleep') as mock_sleep:
            for _ in range(10):
                monitoring_pipeline.monitor_reddit(test_compound)
            
            # Check rate limiting was applied
            assert mock_sleep.called


if __name__ == "__main__":
    pytest.main([__file__])
