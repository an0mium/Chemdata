"""Shared test fixtures and configuration."""

import logging
import pytest

from binding_data_processor.pipeline.infrastructure.monitoring import (
    MetricsCollector,
)


@pytest.fixture
def logger():
    """Create test logger."""
    logger = logging.getLogger("test_web_enrichment")
    logger.setLevel(logging.DEBUG)
    return logger


@pytest.fixture
def metrics_collector(logger):
    """Create test metrics collector."""
    return MetricsCollector(
        namespace="test_web_enrichment",
        logger=logger,
    )


@pytest.fixture
def mock_config():
    """Create test enrichment config."""
    return {
        "reddit_client_id": "test_id",
        "reddit_client_secret": "test_secret",
        "twitter_bearer_token": "test_token",
        "skip_predictions": False,
        "skip_web_data": False,
        "use_cache": True,
        "n_workers": 1,
        "batch_size": 1,
        "model_dir": "test_models",
        "cache_dir": "test_cache",
        "circuit_config": {
            "failure_threshold": 3,
            "recovery_timeout": 60,
        },
        "client_configs": {
            "swiss": {
                "base_url": "https://test.swisstargetprediction.ch/api",
            },
            "community": {
                "base_url": "https://test.psychonautwiki.org",
            },
            "social": {
                "base_url": "https://test.reddit.com",
            },
        },
    }
