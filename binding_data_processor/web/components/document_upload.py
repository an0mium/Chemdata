"""Web component for document upload and processing."""

import logging
from pathlib import Path
from typing import Optional

from fastapi import FastAPI

from ...pipeline.infrastructure.monitoring import Monitor
from ...processors.document.monitor import DirectoryMonitor
from ...processors.document.pdf import PDFProcessor
from .document_routes import DocumentRouteHandlers


class DocumentUploadComponent:
    """Component for handling document uploads and processing."""

    def __init__(
        self,
        app: FastAPI,
        storage_dir: str = "data/documents",
        monitor: Optional[Monitor] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize component.

        Args:
            app: FastAPI application instance
            storage_dir: Base directory for storing documents
            monitor: Optional monitoring instance
            logger: Optional logger instance
        """
        self.storage_dir = Path(storage_dir)
        self.storage_dir.mkdir(parents=True, exist_ok=True)

        self.monitor = monitor or Monitor()
        self.logger = logger or logging.getLogger(__name__)

        # Initialize processors
        self.pdf_processor = PDFProcessor(
            storage_dir=self.storage_dir / "pdf",
            monitor=self.monitor,
            logger=self.logger,
        )

        # Initialize directory monitor
        self.directory_monitor = DirectoryMonitor(
            processor=self.pdf_processor,
            watch_dirs={},
            monitor=self.monitor,
            logger=self.logger,
        )

        # Initialize route handlers
        self.route_handlers = DocumentRouteHandlers(
            storage_dir=self.storage_dir,
            pdf_processor=self.pdf_processor,
            directory_monitor=self.directory_monitor,
            monitor=self.monitor,
            logger=self.logger,
        )

        # Register routes
        self._register_routes(app)

    def _register_routes(self, app: FastAPI):
        """Register component routes.

        Args:
            app: FastAPI application instance
        """
        app.post("/documents/upload")(self.route_handlers.handle_upload)
        app.post("/documents/monitor")(self.route_handlers.handle_monitor_config)
        app.delete("/documents/monitor")(self.route_handlers.handle_monitor_removal)
        app.get("/documents/status/{file_path:path}")(self.route_handlers.get_processing_status)

    async def __aenter__(self):
        """Enter async context."""
        self.directory_monitor.start()
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Exit async context."""
        self.directory_monitor.stop()
