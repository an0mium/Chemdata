"""Tests for social media monitoring functionality."""

import pytest
from unittest.mock import Mock, patch
from pathlib import Path
import json

from ....models.compound.enhanced import EnhancedCompound
from ..social import SocialMonitor, SocialConfig
from ..config import SocialMedia


@pytest.fixture
def mock_reddit():
    """Create mock Reddit client."""
    mock = Mock()
    
    # Mock subreddit
    mock_subreddit = Mock()
    mock_subreddit.name = "researchchemicals"
    
    # Mock posts
    mock_post = Mock()
    mock_post.title = "New tryptamine synthesis"
    mock_post.selftext = "Synthesized 4-HO-DMT using..."
    
    # Mock comments
    mock_comment = Mock()
    mock_comment.body = "The binding affinity for 5-HT2A is..."
    mock_post.comments.list.return_value = [mock_comment]
    mock_post.comments.replace_more.return_value = None
    
    # Setup mock chain
    mock_subreddit.new.return_value = [mock_post]
    mock.subreddit.return_value = mock_subreddit
    
    return mock


@pytest.fixture
def mock_twitter():
    """Create mock Twitter client."""
    mock = Mock()
    
    # Mock tweets
    mock_tweet = Mock()
    mock_tweet.text = "Novel NMDA antagonist discovered..."
    
    # Setup mock chain
    mock.search_tweets.return_value = [mock_tweet]
    
    return mock


@pytest.fixture
def mock_bluesky():
    """Create mock Bluesky client."""
    mock = Mock()
    
    # Mock posts
    mock_post = Mock()
    mock_post.text = "Interesting paper on ketamine analogs..."
    
    # Setup mock chain
    mock.search_posts.return_value = [mock_post]
    
    return mock


@pytest.fixture
def mock_llm():
    """Create mock LLM processor."""
    mock = Mock()
    
    # Mock compound extraction
    mock.extract_compounds.return_value = [
        EnhancedCompound(
            name="4-HO-DMT",
            smiles="CC1=CC=CC=C1",
            source="reddit",
        )
    ]
    
    return mock


@pytest.fixture
def social_monitor(mock_reddit, mock_twitter, mock_bluesky, mock_llm):
    """Create social monitor with mock components."""
    config = SocialConfig(
        reddit_client_id="test_id",
        reddit_client_secret="test_secret",
        twitter_api_key="test_key",
        twitter_api_secret="test_secret",
        bluesky_handle="test_handle",
        bluesky_password="test_pass",
    )
    
    with patch("praw.Reddit", return_value=mock_reddit), \
         patch("tweepy.API", return_value=mock_twitter), \
         patch("atproto.Client", return_value=mock_bluesky), \
         patch("binding_data_processor.web_enrichment.llm_utils.LLMProcessor", 
               return_value=mock_llm):
        monitor = SocialMonitor(config=config)
        yield monitor


class TestSocialMonitor:
    """Tests for SocialMonitor class."""

    def test_initialization(self, social_monitor):
        """Test initialization of social monitor."""
        assert social_monitor.reddit is not None
        assert social_monitor.twitter is not None
        assert social_monitor.bluesky is not None
        assert social_monitor.llm is not None

    def test_monitor_reddit(self, social_monitor, mock_reddit, mock_llm):
        """Test Reddit monitoring."""
        # Monitor Reddit
        compounds = social_monitor._monitor_reddit()
        
        # Check results
        assert len(compounds) == 1
        assert compounds[0].name == "4-HO-DMT"
        assert compounds[0].source == "reddit"
        
        # Verify calls
        mock_reddit.subreddit.assert_called()
        mock_llm.extract_compounds.assert_called()

    def test_monitor_twitter(self, social_monitor, mock_twitter, mock_llm):
        """Test Twitter monitoring."""
        # Monitor Twitter
        compounds = social_monitor._monitor_twitter()
        
        # Check results
        assert len(compounds) == 1
        assert compounds[0].name == "4-HO-DMT"
        assert compounds[0].source == "twitter"
        
        # Verify calls
        mock_twitter.search_tweets.assert_called()
        mock_llm.extract_compounds.assert_called()

    def test_monitor_bluesky(self, social_monitor, mock_bluesky, mock_llm):
        """Test Bluesky monitoring."""
        # Monitor Bluesky
        compounds = social_monitor._monitor_bluesky()
        
        # Check results
        assert len(compounds) == 1
        assert compounds[0].name == "4-HO-DMT"
        assert compounds[0].source == "bluesky"
        
        # Verify calls
        mock_bluesky.search_posts.assert_called()
        mock_llm.extract_compounds.assert_called()

    def test_monitor_all(self, social_monitor):
        """Test monitoring all sources."""
        # Monitor all sources
        compounds = social_monitor.monitor_all()
        
        # Check results
        assert len(compounds) == 3  # One from each source
        assert all(isinstance(c, EnhancedCompound) for c in compounds)
        
        # Check sources
        sources = {c.source for c in compounds}
        assert sources == {"reddit", "twitter", "bluesky"}

    def test_error_handling(self, social_monitor, mock_reddit):
        """Test error handling."""
        # Make Reddit raise error
        mock_reddit.subreddit.side_effect = Exception("API error")
        
        # Monitor should handle error and continue
        compounds = social_monitor.monitor_all()
        
        # Should still get compounds from other sources
        assert len(compounds) == 2
        assert all(c.source != "reddit" for c in compounds)

    def test_batch_processing(self, social_monitor, mock_llm):
        """Test batch processing of texts."""
        # Create many texts
        texts = ["text"] * 150  # More than batch_size
        
        # Process in monitor
        with patch.object(social_monitor.reddit, "subreddit") as mock_subreddit:
            # Setup mock to return many posts
            mock_post = Mock()
            mock_post.title = "title"
            mock_post.selftext = "text"
            mock_post.comments.list.return_value = []
            mock_post.comments.replace_more.return_value = None
            mock_subreddit.return_value.new.return_value = [mock_post] * 50
            
            # Monitor Reddit
            compounds = social_monitor._monitor_reddit()
            
            # Verify batch processing
            assert len(compounds) > 0
            assert mock_llm.extract_compounds.call_count > 1
            
            # Verify texts were processed
            total_texts = sum(
                len(args[0]) for args, _ in mock_llm.extract_compounds.call_args_list
            )
            assert total_texts >= len(texts)

    def test_config_validation(self):
        """Test configuration validation."""
        # Create invalid config
        config = SocialConfig(
            post_limit=-1,
            comment_limit=0,
            batch_size=0,
        )
        
        # Should raise error
        with pytest.raises(ValueError):
            SocialMonitor(config=config)

    def test_source_filtering(self, social_monitor):
        """Test source filtering."""
        # Monitor with filtered sources
        social_monitor.config.reddit_client_id = None  # Disable Reddit
        compounds = social_monitor.monitor_all()
        
        # Should only get compounds from enabled sources
        assert len(compounds) == 2
        assert all(c.source != "reddit" for c in compounds)

    def test_data_extraction(self, social_monitor, mock_llm):
        """Test compound data extraction."""
        # Create test data
        test_data = {
            "name": "Test Compound",
            "smiles": "CC1=CC=CC=C1",
            "source": "test",
            "data": {"key": "value"},
        }
        
        # Make LLM return test data
        mock_llm.extract_compounds.return_value = [
            EnhancedCompound(**test_data)
        ]
        
        # Monitor sources
        compounds = social_monitor.monitor_all()
        
        # Check extracted data
        assert len(compounds) == 3
        for compound in compounds:
            assert compound.name == test_data["name"]
            assert compound.smiles == test_data["smiles"]
            assert compound.data == test_data["data"]
