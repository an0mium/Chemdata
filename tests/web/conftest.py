"""Shared fixtures and configuration for web tests."""

import pytest
from binding_data_processor.web import (
    request,
    response,
    session,
    database,
    cache,
    auth,
    templates,
    static,
    config,
)


@pytest.fixture
def mock_request():
    """Create a mock request for testing."""
    return request.RequestManager()


@pytest.fixture
def mock_response():
    """Create a mock response for testing."""
    return response.ResponseManager()


@pytest.fixture
def mock_session():
    """Create a mock session for testing."""
    return session.SessionManager()


@pytest.fixture
def mock_database():
    """Create a mock database for testing."""
    return database.DatabaseManager()


@pytest.fixture
def mock_cache():
    """Create a mock cache for testing."""
    return cache.CacheManager()


@pytest.fixture
def mock_auth():
    """Create a mock auth manager for testing."""
    return auth.AuthManager()


@pytest.fixture
def mock_templates():
    """Create a mock template manager for testing."""
    return templates.TemplateManager()


@pytest.fixture
def mock_static():
    """Create a mock static file manager for testing."""
    return static.StaticFileManager()


@pytest.fixture
def mock_config():
    """Create a mock configuration manager for testing."""
    return config.ConfigManager()


@pytest.fixture
def test_app():
    """Create a test application instance."""
    from binding_data_processor.web.app import create_app

    app = create_app({"TESTING": True})
    return app


@pytest.fixture
def test_client(test_app):
    """Create a test client for making requests."""
    return test_app.test_client()


@pytest.fixture
def test_database():
    """Set up a test database."""
    # Set up test database
    yield
    # Clean up test database
