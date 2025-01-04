"""Tests for web server functionality."""

import pytest
import concurrent.futures
from fastapi.testclient import TestClient
from fastapi import status
from pathlib import Path
import json
import logging
from unittest.mock import patch, Mock

from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..web import (
    WebServer,
    ServerConfig,
    CompoundAPI,
    SearchAPI,
    FilterAPI,
    ExportAPI,
    RouteManager,
    APIManager,
    SecurityManager,
    LogManager,
    CacheManager,
    AuthConfig,
    RateLimitConfig,
    ServerResult,
)


@pytest.fixture
def test_compounds():
    """Create test compounds fixture."""
    compounds = []
    
    # Caffeine
    caffeine = PsychoactiveCompound(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    caffeine.psychoactive_class = PsychoactiveClass.STIMULANT
    caffeine.add_receptor_binding(
        "A2A",
        affinity=0.8,
        confidence=0.95,
        activity="antagonist"
    )
    caffeine.effect_profile = {
        "stimulation": (0.8, 0.9),
        "focus": (0.7, 0.8),
    }
    caffeine.safety_alerts = {
        "anxiety": RiskLevel.MODERATE,
        "insomnia": RiskLevel.HIGH,
    }
    compounds.append(caffeine)
    
    # Amphetamine
    amphetamine = PsychoactiveCompound(
        name="Amphetamine",
        smiles="CC(N)CC1=CC=CC=C1",
        cas_number="300-62-9",
    )
    amphetamine.psychoactive_class = PsychoactiveClass.STIMULANT
    amphetamine.add_receptor_binding(
        "DAT",
        affinity=0.05,
        confidence=0.95,
        activity="inhibitor"
    )
    amphetamine.effect_profile = {
        "stimulation": (0.9, 0.95),
        "euphoria": (0.8, 0.9),
    }
    amphetamine.safety_alerts = {
        "addiction": RiskLevel.HIGH,
        "cardiovascular": RiskLevel.HIGH,
    }
    compounds.append(amphetamine)
    
    return compounds


@pytest.fixture
def server_config():
    """Create test server configuration fixture."""
    return ServerConfig(
        host="localhost",
        port=8000,
        debug=True,
        api_prefix="/api/v1",
        cors_origins=["http://localhost:3000"],
        rate_limit=100,
        rate_window=60,
        log_dir=Path("logs"),
        cache_dir=Path("cache"),
        secret_key="test-secret-key",
        auth=AuthConfig(
            enabled=True,
            api_key_header="X-API-Key",
            api_key="test_key_123",
        ),
        rate_limit_config=RateLimitConfig(
            enabled=True,
            requests_per_minute=60,
        ),
        ssl_enabled=False
    )


@pytest.fixture
def web_server(server_config, test_compounds):
    """Create test web server fixture."""
    server = WebServer(config=server_config)
    server.load_compounds(test_compounds)
    return server


@pytest.fixture
def test_client(web_server):
    """Create test client fixture."""
    return TestClient(web_server.app)


@pytest.fixture
def auth_headers():
    """Create test auth headers fixture."""
    return {"X-API-Key": "test_key_123"}


class TestServerInfrastructure:
    """Tests for server infrastructure."""

    def test_initialization(self, web_server, server_config):
        """Test initialization of WebServer."""
        assert isinstance(web_server.routes, RouteManager)
        assert isinstance(web_server.api, APIManager)
        assert isinstance(web_server.security, SecurityManager)
        assert isinstance(web_server.logs, LogManager)
        assert isinstance(web_server.cache, CacheManager)
        assert isinstance(web_server.compound_api, CompoundAPI)
        assert isinstance(web_server.search_api, SearchAPI)
        assert isinstance(web_server.filter_api, FilterAPI)
        assert isinstance(web_server.export_api, ExportAPI)
        assert web_server.config == server_config
        assert web_server.stats == {}

    def test_load_compounds(self, web_server, test_compounds):
        """Test compound loading."""
        # Load compounds
        result = web_server.load_compounds(test_compounds)
        
        # Check result
        assert isinstance(result, ServerResult)
        assert result.success
        assert len(result.compounds) == 2
        assert all(isinstance(c, PsychoactiveCompound) for c in result.compounds)
        assert "data_loading" in web_server.stats

    def test_health_check(self, test_client):
        """Test health check endpoint."""
        response = test_client.get("/api/v1/health")
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert data["status"] == "healthy"
        assert "uptime" in data
        assert "memory_usage" in data
        assert "active_connections" in data

    def test_metrics(self, test_client):
        """Test metrics endpoint."""
        # Generate some traffic
        test_client.get("/api/v1/compounds")
        test_client.get("/api/v1/compounds/58-08-2")
        test_client.get("/api/v1/invalid")
        
        # Check metrics
        response = test_client.get("/api/v1/metrics")
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert data["requests_total"] > 0
        assert data["requests_by_endpoint"]["/api/v1/compounds"] == 1
        assert data["errors_total"] == 1
        assert data["average_response_time"] > 0

    def test_api_docs(self, test_client):
        """Test API documentation endpoints."""
        # Test OpenAPI schema
        response = test_client.get("/openapi.json")
        assert response.status_code == status.HTTP_200_OK
        schema = response.json()
        assert "paths" in schema
        assert "/api/v1/compounds" in schema["paths"]
        assert "/api/v1/search" in schema["paths"]
        assert "/api/v1/export" in schema["paths"]
        
        # Test Swagger UI
        response = test_client.get("/docs")
        assert response.status_code == status.HTTP_200_OK
        assert "text/html" in response.headers["content-type"]


class TestCompoundAPI:
    """Tests for CompoundAPI class."""

    def test_get_compounds(self, test_client, auth_headers):
        """Test getting compound list."""
        response = test_client.get("/api/v1/compounds", headers=auth_headers)
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert len(data["compounds"]) == 2
        assert data["compounds"][0]["name"] == "Caffeine"
        assert data["compounds"][1]["name"] == "Amphetamine"

    def test_get_compound_detail(self, test_client, auth_headers):
        """Test getting compound detail."""
        response = test_client.get(
            "/api/v1/compounds/58-08-2",  # Caffeine CAS
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert data["name"] == "Caffeine"
        assert data["cas_number"] == "58-08-2"
        assert data["psychoactive_class"] == "STIMULANT"
        assert "A2A" in str(data["receptor_bindings"])

    def test_get_nonexistent_compound(self, test_client, auth_headers):
        """Test getting nonexistent compound."""
        response = test_client.get(
            "/api/v1/compounds/invalid-cas",
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_404_NOT_FOUND
        assert "not found" in response.json()["detail"].lower()

    def test_sort_compounds(self, test_client, auth_headers):
        """Test compound sorting."""
        # Sort by name ascending
        response = test_client.get(
            "/api/v1/compounds?sort_by=name&sort_order=asc",
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert data["compounds"][0]["name"] == "Amphetamine"
        assert data["compounds"][1]["name"] == "Caffeine"
        
        # Sort by name descending
        response = test_client.get(
            "/api/v1/compounds?sort_by=name&sort_order=desc",
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert data["compounds"][0]["name"] == "Caffeine"
        assert data["compounds"][1]["name"] == "Amphetamine"

    def test_pagination(self, test_client, auth_headers):
        """Test compound list pagination."""
        # Get first page
        response = test_client.get(
            "/api/v1/compounds?page=1&page_size=1",
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert len(data["compounds"]) == 1
        assert data["total_pages"] == 2
        assert data["current_page"] == 1
        assert data["compounds"][0]["name"] == "Caffeine"
        
        # Get second page
        response = test_client.get(
            "/api/v1/compounds?page=2&page_size=1",
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert len(data["compounds"]) == 1
        assert data["total_pages"] == 2
        assert data["current_page"] == 2
        assert data["compounds"][0]["name"] == "Amphetamine"


class TestSearchAPI:
    """Tests for SearchAPI class."""

    def test_text_search(self, test_client, auth_headers):
        """Test text search endpoint."""
        # Search by name
        response = test_client.get(
            "/api/v1/search?query=caffeine&type=name",
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert len(data["compounds"]) == 1
        assert data["compounds"][0]["name"] == "Caffeine"
        
        # Search by receptor
        response = test_client.get(
            "/api/v1/search?query=DAT&type=receptor",
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert len(data["compounds"]) == 1
        assert data["compounds"][0]["name"] == "Amphetamine"

    def test_structure_search(self, test_client, auth_headers):
        """Test structure search endpoint."""
        params = {
            "smiles": "CN1C=NC2=C1C(=O)N",
            "similarity_threshold": 0.8
        }
        response = test_client.post(
            "/api/v1/search/structure",
            json=params,
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert len(data["compounds"]) == 1
        assert data["compounds"][0]["name"] == "Caffeine"

    def test_invalid_search(self, test_client, auth_headers):
        """Test invalid search request."""
        # Missing query
        response = test_client.get(
            "/api/v1/search",
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_400_BAD_REQUEST
        
        # Invalid search type
        response = test_client.get(
            "/api/v1/search?query=test&type=invalid",
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_400_BAD_REQUEST


class TestFilterAPI:
    """Tests for FilterAPI class."""

    def test_property_filter(self, test_client, auth_headers):
        """Test property filter endpoint."""
        params = {
            "property_name": "binding_affinity",
            "min_value": 0.7,
            "max_value": 1.0
        }
        response = test_client.post(
            "/api/v1/filter/property",
            json=params,
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert len(data["compounds"]) == 1
        assert data["compounds"][0]["name"] == "Caffeine"

    def test_class_filter(self, test_client, auth_headers):
        """Test class filter endpoint."""
        # Filter by class
        params = {
            "classes": ["STIMULANT"]
        }
        response = test_client.post(
            "/api/v1/filter/class",
            json=params,
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert len(data["compounds"]) == 2  # Both are stimulants
        
        # Filter by risk level
        params = {
            "min_risk_level": "HIGH"
        }
        response = test_client.post(
            "/api/v1/filter/risk",
            json=params,
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert len(data["compounds"]) == 2  # Both have high risks

    def test_invalid_filter(self, test_client, auth_headers):
        """Test invalid filter request."""
        params = {
            "property_name": "invalid_property",
            "min_value": 0,
            "max_value": 1
        }
        response = test_client.post(
            "/api/v1/filter/property",
            json=params,
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_400_BAD_REQUEST


class TestExportAPI:
    """Tests for ExportAPI class."""

    def test_export_tsv(self, test_client, auth_headers, tmp_path):
        """Test TSV export endpoint."""
        params = {
            "format": "tsv",
            "columns": ["name", "cas_number", "smiles"],
            "include_metadata": True,
            "output_file": str(tmp_path / "compounds.tsv")
        }
        response = test_client.post(
            "/api/v1/export",
            json=params,
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        assert "text/tab-separated-values" in response.headers["content-type"]
        
        # Check export file
        output_file = tmp_path / "compounds.tsv"
        assert output_file.exists()
        content = output_file.read_text()
        assert "name\tcas_number\tsmiles" in content
        assert "Caffeine\t58-08-2" in content

    def test_export_json(self, test_client, auth_headers):
        """Test JSON export endpoint."""
        params = {
            "format": "json",
            "columns": ["name", "cas_number"],
            "include_metadata": True
        }
        response = test_client.post(
            "/api/v1/export",
            json=params,
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        assert response.headers["content-type"] == "application/json"
        data = response.json()
        assert len(data["compounds"]) == 2
        assert all("name" in c for c in data["compounds"])
        assert all("cas_number" in c for c in data["compounds"])

    def test_invalid_export(self, test_client, auth_headers):
        """Test invalid export request."""
        params = {
            "format": "invalid_format",
            "columns": ["name"]
        }
        response = test_client.post(
            "/api/v1/export",
            json=params,
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_400_BAD_REQUEST


class TestSecurity:
    """Tests for security features."""

    def test_authentication(self, test_client):
        """Test authentication."""
        # Missing API key
        response = test_client.get("/api/v1/compounds")
        assert response.status_code == status.HTTP_401_UNAUTHORIZED
        assert response.json()["detail"] == "Missing API key"
        
        # Invalid API key
        response = test_client.get(
            "/api/v1/compounds",
            headers={"X-API-Key": "wrong_key"}
        )
        assert response.status_code == status.HTTP_401_UNAUTHORIZED
        assert response.json()["detail"] == "Invalid API key"

    def test_rate_limiting(self, test_client, auth_headers):
        """Test rate limiting."""
        # Make requests up to limit
        for _ in range(60):
            response = test_client.get(
                "/api/v1/compounds",
                headers=auth_headers
            )
            assert response.status_code == status.HTTP_200_OK
        
        # Next request should be rate limited
        response = test_client.get(
            "/api/v1/compounds",
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_429_TOO_MANY_REQUESTS
        assert "Rate limit exceeded" in response.json()["detail"]

    def test_security_features(self, test_client, auth_headers):
        """Test security features."""
        # Test XSS protection
        response = test_client.get(
            "/api/v1/compounds?q=<script>alert('xss')</script>",
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert "<script>" not in json.dumps(data)
        
        # Test SQL injection protection
        response = test_client.get(
            "/api/v1/compounds?q='; DROP TABLE compounds;--",
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert all("DROP TABLE" not in str(item) for item in data["compounds"])
        
        # Test CSRF protection
        response = test_client.post(
            "/api/v1/filter/class",
            headers={"X-CSRF-Token": "invalid-token", **auth_headers},
            json={"classes": ["STIMULANT"]}
        )
        assert response.status_code == status.HTTP_403_FORBIDDEN
        assert "csrf" in response.json()["detail"].lower()


class TestInfrastructure:
    """Tests for server infrastructure."""

    def test_logging(self, web_server, test_client, auth_headers):
        """Test logging functionality."""
        # Configure test logger
        test_log = []
        handler = logging.StreamHandler(test_log)
        web_server.logs.logger.addHandler(handler)
        
        # Test access logging
        test_client.get("/api/v1/compounds", headers=auth_headers)
        assert any("GET /api/v1/compounds" in log for log in test_log)
        
        # Test error logging
        test_client.get("/api/v1/invalid", headers=auth_headers)
        assert any("404 Not Found" in log for log in test_log)
        
        # Test warning logging
        test_client.post(
            "/api/v1/filter/property",
            json={"invalid": "data"},
            headers=auth_headers
        )
        assert any("Validation Error" in log for log in test_log)

    def test_caching(self, test_client, auth_headers):
        """Test caching functionality."""
        # Make initial request
        response1 = test_client.get("/api/v1/compounds", headers=auth_headers)
        assert response1.status_code == status.HTTP_200_OK
        assert "X-Cache" not in response1.headers
        
        # Make same request again
        response2 = test_client.get("/api/v1/compounds", headers=auth_headers)
        assert response2.status_code == status.HTTP_200_OK
        assert response2.headers["X-Cache"] == "HIT"
        
        # Make different request
        response3 = test_client.get(
            "/api/v1/compounds/58-08-2",
            headers=auth_headers
        )
        assert response3.status_code == status.HTTP_200_OK
        assert "X-Cache" not in response3.headers

    def test_ssl_configuration(self, server_config):
        """Test SSL configuration."""
        # Enable SSL
        server_config.ssl_enabled = True
        server_config.ssl_cert = Path("cert.pem")
        server_config.ssl_key = Path("key.pem")
        
        # Create server with SSL
        with patch("ssl.SSLContext") as mock_ssl:
            server = WebServer(config=server_config)
            assert server.ssl_context is not None
            mock_ssl.assert_called_once()

    def test_graceful_shutdown(self, web_server):
        """Test graceful shutdown."""
        # Start server
        web_server.start()
        
        # Mock active connections
        active_connections = [Mock(), Mock()]
        web_server.connections = active_connections
        
        # Shutdown server
        web_server.shutdown()
        
        # Check connections were closed
        for conn in active_connections:
            conn.close.assert_called_once()
        
        # Check server stopped
        assert not web_server.is_running
        assert web_server.stats["uptime"] > 0

    def test_concurrent_requests(self, test_client, auth_headers):
        """Test handling concurrent requests."""
        def make_request():
            return test_client.get(
                "/api/v1/compounds",
                headers=auth_headers
            )
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
            futures = [executor.submit(make_request) for _ in range(4)]
            responses = [f.result() for f in futures]
        
        assert all(r.status_code == status.HTTP_200_OK for r in responses)


if __name__ == "__main__":
    pytest.main([__file__])
