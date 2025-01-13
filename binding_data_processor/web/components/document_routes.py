"""Route handlers for document processing."""

import asyncio
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set
from uuid import UUID, uuid4

from fastapi import File, Form, UploadFile, BackgroundTasks
from fastapi.responses import JSONResponse
from pydantic import BaseModel


from ...processors.document.pdf import PDFProcessor
from ...processors.document.monitor import DirectoryMonitor
from ...pipeline.infrastructure.monitoring import Monitor


class DirectoryConfig(BaseModel):
    """Configuration for monitored directory."""

    path: str
    patterns: set[str] = {"*.pdf"}
    recursive: bool = True


class ProcessingStatus(BaseModel):
    """Status of document processing."""

    file_path: str
    status: str
    compounds_found: int
    error_message: Optional[str] = None
    processing_time: float


class BatchUploadStatus(BaseModel):
    """Status of batch upload processing."""

    batch_id: UUID
    total_files: int
    processed_files: int
    failed_files: int
    status: str  # "queued", "processing", "completed", "failed"
    files: Dict[str, ProcessingStatus]
    error_message: Optional[str] = None
    processing_time: float = 0.0
    start_time: Optional[float] = None
    completion_time: Optional[float] = None


class DocumentRouteHandlers:
    """Route handlers for document processing endpoints."""

    def __init__(
        self,
        storage_dir: Path,
        pdf_processor: PDFProcessor,
        directory_monitor: Optional[DirectoryMonitor] = None,
        monitor: Optional[Monitor] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize handlers.

        Args:
            storage_dir: Base directory for storing documents
            pdf_processor: PDF processor instance
            directory_monitor: Optional directory monitor instance
            monitor: Optional monitoring instance
            logger: Optional logger instance
        """
        self.storage_dir = storage_dir
        self.pdf_processor = pdf_processor
        self.directory_monitor = directory_monitor
        self.monitor = monitor or Monitor()
        self.logger = logger or logging.getLogger(__name__)
        self.processing_status = {}
        self.batch_status: Dict[UUID, BatchUploadStatus] = {}
        self.processing_queue: Set[str] = set()

    async def handle_upload(
        self,
        files: List[UploadFile] = File(...),
        process: bool = Form(True),
    ) -> List[ProcessingStatus]:
        """Handle document upload request.

        Args:
            files: List of files to upload
            process: Whether to process files after upload

        Returns:
            List of processing status objects
        """
        results = []

        for file in files:
            try:
                # Save file
                file_path = self.storage_dir / file.filename
                content = await file.read()
                file_path.write_bytes(content)

                status = ProcessingStatus(
                    file_path=str(file_path),
                    status="uploaded",
                    compounds_found=0,
                    processing_time=0.0,
                )

                # Process if requested
                if process and file_path.suffix.lower() == ".pdf":
                    try:
                        result = await self.pdf_processor.process_file(file_path)
                        status.status = "processed" if result.success else "failed"
                        status.compounds_found = len(result.compounds)
                        status.processing_time = result.processing_time
                        if result.error_messages:
                            status.error_message = "; ".join(result.error_messages)
                    except Exception as e:
                        status.status = "failed"
                        status.error_message = str(e)

                results.append(status)
                self.processing_status[str(file_path)] = status

            except Exception as e:
                self.logger.error(f"Failed to handle {file.filename}: {str(e)}")
                results.append(
                    ProcessingStatus(
                        file_path=file.filename,
                        status="failed",
                        compounds_found=0,
                        error_message=str(e),
                        processing_time=0.0,
                    )
                )

        return results

    async def handle_batch_upload(
        self,
        files: List[UploadFile] = File(...),
        process: bool = Form(True),
        background_tasks: BackgroundTasks = None,
    ) -> BatchUploadStatus:
        """Handle batch document upload request.

        Args:
            files: List of files to upload
            process: Whether to process files after upload
            background_tasks: FastAPI background tasks

        Returns:
            Batch upload status object
        """
        batch_id = uuid4()
        batch_status = BatchUploadStatus(
            batch_id=batch_id,
            total_files=len(files),
            processed_files=0,
            failed_files=0,
            status="queued",
            files={},
            start_time=None,
            completion_time=None,
        )
        self.batch_status[batch_id] = batch_status

        # Add background processing task
        if background_tasks:
            background_tasks.add_task(
                self._process_batch,
                batch_id=batch_id,
                files=files,
                process=process,
            )

        return batch_status

    async def _process_batch(
        self,
        batch_id: UUID,
        files: List[UploadFile],
        process: bool,
    ):
        """Process batch of files in background.

        Args:
            batch_id: Batch identifier
            files: List of files to process
            process: Whether to process files after upload
        """
        import time

        batch = self.batch_status[batch_id]
        batch.status = "processing"
        batch.start_time = time.time()

        try:
            for file in files:
                if file.filename in self.processing_queue:
                    continue

                self.processing_queue.add(file.filename)
                try:
                    # Save file
                    file_path = self.storage_dir / file.filename
                    content = await file.read()
                    file_path.write_bytes(content)

                    status = ProcessingStatus(
                        file_path=str(file_path),
                        status="uploaded",
                        compounds_found=0,
                        processing_time=0.0,
                    )

                    # Process if requested
                    if process and file_path.suffix.lower() == ".pdf":
                        try:
                            result = await self.pdf_processor.process_file(file_path)
                            status.status = "processed" if result.success else "failed"
                            status.compounds_found = len(result.compounds)
                            status.processing_time = result.processing_time
                            if result.error_messages:
                                status.error_message = "; ".join(result.error_messages)
                            if status.status == "processed":
                                batch.processed_files += 1
                            else:
                                batch.failed_files += 1
                        except Exception as e:
                            status.status = "failed"
                            status.error_message = str(e)
                            batch.failed_files += 1

                    batch.files[file.filename] = status
                    self.processing_status[str(file_path)] = status

                except Exception as e:
                    self.logger.error(f"Failed to handle {file.filename}: {str(e)}")
                    status = ProcessingStatus(
                        file_path=file.filename,
                        status="failed",
                        compounds_found=0,
                        error_message=str(e),
                        processing_time=0.0,
                    )
                    batch.files[file.filename] = status
                    batch.failed_files += 1

                finally:
                    self.processing_queue.remove(file.filename)

            batch.status = "completed"
            batch.completion_time = time.time()
            batch.processing_time = batch.completion_time - batch.start_time

        except Exception as e:
            batch.status = "failed"
            batch.error_message = str(e)
            batch.completion_time = time.time()
            batch.processing_time = batch.completion_time - batch.start_time

    async def handle_monitor_config(self, config: DirectoryConfig) -> JSONResponse:
        """Handle monitor configuration request.

        Args:
            config: Directory monitoring configuration

        Returns:
            JSON response indicating success/failure
        """
        try:
            path = Path(config.path)

            # Update monitoring
            if self.directory_monitor:
                self.directory_monitor.add_directory(path, config.patterns)
            else:
                return JSONResponse(
                    content={
                        "status": "error",
                        "message": "Directory monitor not initialized",
                    },
                    status_code=500,
                )

            return JSONResponse(
                content={
                    "status": "success",
                    "message": f"Now monitoring {path} for {config.patterns}",
                }
            )

        except Exception as e:
            self.logger.error(f"Failed to configure monitoring: {str(e)}")
            return JSONResponse(
                content={
                    "status": "error",
                    "message": str(e),
                },
                status_code=500,
            )

    async def handle_monitor_removal(self, path: str) -> JSONResponse:
        """Handle monitor removal request.

        Args:
            path: Directory path to stop monitoring

        Returns:
            JSON response indicating success/failure
        """
        try:
            if self.directory_monitor:
                self.directory_monitor.remove_directory(path)
                return JSONResponse(
                    content={
                        "status": "success",
                        "message": f"Stopped monitoring {path}",
                    }
                )
            else:
                return JSONResponse(
                    content={
                        "status": "error",
                        "message": "Directory monitor not initialized",
                    },
                    status_code=500,
                )

        except Exception as e:
            self.logger.error(f"Failed to remove monitoring: {str(e)}")
            return JSONResponse(
                content={
                    "status": "error",
                    "message": str(e),
                },
                status_code=500,
            )

    def get_processing_status(self, file_path: str) -> ProcessingStatus:
        """Get processing status for a file.

        Args:
            file_path: Path to file

        Returns:
            Processing status object
        """
        status = self.processing_status.get(file_path)
        if not status:
            return ProcessingStatus(
                file_path=file_path,
                status="unknown",
                compounds_found=0,
                processing_time=0.0,
            )
        return status

    def get_batch_status(self, batch_id: UUID) -> Optional[BatchUploadStatus]:
        """Get status of batch upload.

        Args:
            batch_id: Batch identifier

        Returns:
            Batch upload status object if found, None otherwise
        """
        return self.batch_status.get(batch_id)

    def handle_processed_file(self, path: Path, result):
        """Handle processed file callback.

        Args:
            path: Path to processed file
            result: Processing result
        """
        try:
            status = ProcessingStatus(
                file_path=str(path),
                status="processed" if result.success else "failed",
                compounds_found=len(result.compounds),
                processing_time=result.processing_time,
            )
            if result.error_messages:
                status.error_message = "; ".join(result.error_messages)
            self.processing_status[str(path)] = status
        except Exception as e:
            self.logger.error(f"Failed to handle processed file {path}: {str(e)}")
