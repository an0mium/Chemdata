"""Tests for document processing components."""

import asyncio
from pathlib import Path
from typing import List
from uuid import UUID

import pytest
from fastapi import FastAPI, UploadFile
from fastapi.testclient import TestClient

from ....pipeline.infrastructure.monitoring import Monitor
from ....processors.document.monitor import DirectoryMonitor
from ....processors.document.pdf import PDFProcessor
from ..document_routes import DocumentRouteHandlers, ProcessingStatus, BatchUploadStatus
from ..document_upload import DocumentUploadComponent


@pytest.fixture
def test_storage_dir(tmp_path):
    """Create temporary storage directory."""
    return tmp_path / "documents"


@pytest.fixture
def test_pdf_processor(test_storage_dir):
    """Create test PDF processor."""
    return PDFProcessor(storage_dir=test_storage_dir / "pdf")


@pytest.fixture
def test_directory_monitor(test_pdf_processor):
    """Create test directory monitor."""
    return DirectoryMonitor(processor=test_pdf_processor, watch_dirs={})


@pytest.fixture
def test_monitor():
    """Create test monitor."""
    return Monitor()


@pytest.fixture
def test_route_handlers(test_storage_dir, test_pdf_processor, test_directory_monitor, test_monitor):
    """Create test route handlers."""
    return DocumentRouteHandlers(
        storage_dir=test_storage_dir,
        pdf_processor=test_pdf_processor,
        directory_monitor=test_directory_monitor,
        monitor=test_monitor,
    )


@pytest.fixture
def test_app():
    """Create test FastAPI app."""
    return FastAPI()


@pytest.fixture
def test_component(test_app, test_storage_dir, test_monitor):
    """Create test document upload component."""
    return DocumentUploadComponent(
        app=test_app,
        storage_dir=str(test_storage_dir),
        monitor=test_monitor,
    )


@pytest.fixture
def test_client(test_app, test_component):
    """Create test client."""
    return TestClient(test_app)


class TestDocumentProcessing:
    """Test document processing functionality."""

    async def test_single_file_upload(self, test_client, test_storage_dir):
        """Test single file upload."""
        # Create test file
        test_file = test_storage_dir / "test.pdf"
        test_file.write_bytes(b"test content")

        # Upload file
        with open(test_file, "rb") as f:
            response = test_client.post(
                "/documents/upload",
                files={"files": ("test.pdf", f, "application/pdf")},
                data={"process": "true"},
            )

        assert response.status_code == 200
        results = response.json()
        assert len(results) == 1
        assert results[0]["file_path"].endswith("test.pdf")
        assert results[0]["status"] in ["uploaded", "processed"]

    async def test_batch_upload(self, test_client, test_storage_dir):
        """Test batch file upload."""
        # Create test files
        files = []
        for i in range(3):
            test_file = test_storage_dir / f"test{i}.pdf"
            test_file.write_bytes(f"test content {i}".encode())
            files.append(("files", (f"test{i}.pdf", open(test_file, "rb"), "application/pdf")))

        # Upload files
        response = test_client.post(
            "/documents/batch-upload",
            files=files,
            data={"process": "true"},
        )

        assert response.status_code == 200
        result = response.json()
        assert "batch_id" in result
        assert result["total_files"] == 3
        assert result["status"] == "queued"

        # Check batch status
        batch_id = result["batch_id"]
        response = test_client.get(f"/documents/batch-status/{batch_id}")
        assert response.status_code == 200
        status = response.json()
        assert status["batch_id"] == batch_id
        assert status["total_files"] == 3
        assert status["status"] in ["queued", "processing", "completed"]

        # Clean up
        for file in files:
            file[1][1].close()

    async def test_directory_monitoring(self, test_client, test_storage_dir):
        """Test directory monitoring configuration."""
        # Configure monitoring
        response = test_client.post(
            "/documents/monitor",
            json={
                "path": str(test_storage_dir),
                "patterns": ["*.pdf"],
                "recursive": True,
            },
        )

        assert response.status_code == 200
        result = response.json()
        assert result["status"] == "success"
        assert str(test_storage_dir) in result["message"]

        # Remove monitoring
        response = test_client.delete(
            "/documents/monitor",
            params={"path": str(test_storage_dir)},
        )

        assert response.status_code == 200
        result = response.json()
        assert result["status"] == "success"
        assert str(test_storage_dir) in result["message"]

    async def test_processing_status(self, test_client, test_storage_dir, test_route_handlers):
        """Test processing status retrieval."""
        # Create test status
        test_file = test_storage_dir / "test.pdf"
        status = ProcessingStatus(
            file_path=str(test_file),
            status="processed",
            compounds_found=5,
            processing_time=1.23,
        )
        test_route_handlers.processing_status[str(test_file)] = status

        # Get status
        response = test_client.get(f"/documents/status/{test_file}")
        assert response.status_code == 200
        result = response.json()
        assert result["file_path"] == str(test_file)
        assert result["status"] == "processed"
        assert result["compounds_found"] == 5
        assert result["processing_time"] == 1.23

    async def test_batch_status_tracking(self, test_client, test_storage_dir, test_route_handlers):
        """Test batch status tracking."""
        # Create test batch
        batch_id = UUID("12345678-1234-5678-1234-567812345678")
        batch = BatchUploadStatus(
            batch_id=batch_id,
            total_files=2,
            processed_files=1,
            failed_files=0,
            status="processing",
            files={},
            processing_time=0.5,
        )
        test_route_handlers.batch_status[batch_id] = batch

        # Get batch status
        response = test_client.get(f"/documents/batch-status/{batch_id}")
        assert response.status_code == 200
        result = response.json()
        assert result["batch_id"] == str(batch_id)
        assert result["total_files"] == 2
        assert result["processed_files"] == 1
        assert result["status"] == "processing"

    async def test_error_handling(self, test_client):
        """Test error handling."""
        # Test invalid file upload
        response = test_client.post(
            "/documents/upload",
            files={"files": ("test.txt", b"invalid content", "text/plain")},
        )
        assert response.status_code == 200  # Still returns 200 but with error status
        results = response.json()
        assert len(results) == 1
        assert results[0]["status"] == "failed"
        assert "error_message" in results[0]

        # Test invalid batch ID
        response = test_client.get("/documents/batch-status/invalid-uuid")
        assert response.status_code == 422  # FastAPI validation error

        # Test invalid monitor config
        response = test_client.post(
            "/documents/monitor",
            json={
                "path": "/nonexistent/path",
                "patterns": ["*.pdf"],
            },
        )
        assert response.status_code == 500
        result = response.json()
        assert result["status"] == "error"
        assert "message" in result
