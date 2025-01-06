"""Directory monitoring for document processing."""

import asyncio
import logging
from pathlib import Path
from typing import Callable, Dict, Optional, Set, Union
from watchdog.events import FileSystemEvent, FileSystemEventHandler
from watchdog.observers import Observer

from ...pipeline.infrastructure.monitoring import Monitor
from .base import ProcessingResult
from .pdf import PDFProcessor


class DocumentEventHandler(FileSystemEventHandler):
    """Handler for document file system events."""

    def __init__(
        self,
        processor: PDFProcessor,
        watch_patterns: Set[str],
        callback: Optional[Callable[[Path, ProcessingResult], None]] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize handler.

        Args:
            processor: Document processor instance
            watch_patterns: File patterns to watch (e.g. {"*.pdf"})
            callback: Optional callback for processing results
            logger: Optional logger instance
        """
        self.processor = processor
        self.watch_patterns = watch_patterns
        self.callback = callback
        self.logger = logger or logging.getLogger(__name__)
        self.processing_lock = asyncio.Lock()
        self._processing_files: Set[Path] = set()

    def on_created(self, event: FileSystemEvent):
        """Handle file creation event.

        Args:
            event: File system event
        """
        if not event.is_directory:
            path = Path(event.src_path)
            if any(path.match(pattern) for pattern in self.watch_patterns):
                asyncio.create_task(self._process_file(path))

    def on_modified(self, event: FileSystemEvent):
        """Handle file modification event.

        Args:
            event: File system event
        """
        if not event.is_directory:
            path = Path(event.src_path)
            if any(path.match(pattern) for pattern in self.watch_patterns):
                asyncio.create_task(self._process_file(path))

    async def _process_file(self, path: Path):
        """Process a single file.

        Args:
            path: Path to file
        """
        if path in self._processing_files:
            return

        async with self.processing_lock:
            if path in self._processing_files:
                return
            self._processing_files.add(path)

        try:
            # Wait briefly to ensure file is fully written
            await asyncio.sleep(1)

            # Process file
            result = await self.processor.process_file(path)

            # Call callback if provided
            if self.callback:
                try:
                    self.callback(path, result)
                except Exception as e:
                    self.logger.error(f"Callback failed for {path}: {str(e)}")

        except Exception as e:
            self.logger.error(f"Failed to process {path}: {str(e)}")

        finally:
            self._processing_files.remove(path)


class DirectoryMonitor:
    """Monitor directories for new documents."""

    def __init__(
        self,
        processor: PDFProcessor,
        watch_dirs: Dict[Union[str, Path], Set[str]],
        callback: Optional[Callable[[Path, ProcessingResult], None]] = None,
        monitor: Optional[Monitor] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize monitor.

        Args:
            processor: Document processor instance
            watch_dirs: Dictionary mapping directories to sets of file patterns
            callback: Optional callback for processing results
            monitor: Optional monitoring instance
            logger: Optional logger instance
        """
        self.processor = processor
        self.watch_dirs = {Path(d): patterns for d, patterns in watch_dirs.items()}
        self.callback = callback
        self.monitor = monitor or Monitor()
        self.logger = logger or logging.getLogger(__name__)

        self.observer = Observer()
        self.handlers: Dict[Path, DocumentEventHandler] = {}

    def start(self):
        """Start monitoring directories."""
        for directory, patterns in self.watch_dirs.items():
            try:
                # Create directory if it doesn't exist
                directory.mkdir(parents=True, exist_ok=True)

                # Create and register handler
                handler = DocumentEventHandler(
                    processor=self.processor,
                    watch_patterns=patterns,
                    callback=self.callback,
                    logger=self.logger,
                )
                self.observer.schedule(handler, str(directory), recursive=True)
                self.handlers[directory] = handler

                self.logger.info(f"Monitoring directory {directory} for patterns: {patterns}")

            except Exception as e:
                self.logger.error(f"Failed to monitor {directory}: {str(e)}")

        self.observer.start()

    def stop(self):
        """Stop monitoring directories."""
        self.observer.stop()
        self.observer.join()

    async def process_existing(self):
        """Process existing files in watched directories."""
        for directory, patterns in self.watch_dirs.items():
            try:
                for pattern in patterns:
                    files = list(directory.rglob(pattern))
                    if files:
                        self.logger.info(f"Processing {len(files)} existing files in {directory}")
                        for file_path in files:
                            try:
                                result = await self.processor.process_file(file_path)
                                if self.callback:
                                    self.callback(file_path, result)
                            except Exception as e:
                                self.logger.error(f"Failed to process existing file {file_path}: {str(e)}")

            except Exception as e:
                self.logger.error(f"Failed to process existing files in {directory}: {str(e)}")

    def add_directory(self, directory: Union[str, Path], patterns: Set[str]):
        """Add a new directory to monitor.

        Args:
            directory: Directory path to monitor
            patterns: File patterns to watch
        """
        path = Path(directory)
        if path not in self.watch_dirs:
            self.watch_dirs[path] = patterns
            if self.observer.is_alive():
                try:
                    path.mkdir(parents=True, exist_ok=True)
                    handler = DocumentEventHandler(
                        processor=self.processor,
                        watch_patterns=patterns,
                        callback=self.callback,
                        logger=self.logger,
                    )
                    self.observer.schedule(handler, str(path), recursive=True)
                    self.handlers[path] = handler
                    self.logger.info(f"Added monitoring for {path}: {patterns}")
                except Exception as e:
                    self.logger.error(f"Failed to add monitoring for {path}: {str(e)}")

    def remove_directory(self, directory: Union[str, Path]):
        """Remove a directory from monitoring.

        Args:
            directory: Directory path to stop monitoring
        """
        path = Path(directory)
        if path in self.watch_dirs:
            if self.observer.is_alive():
                try:
                    handler = self.handlers.pop(path, None)
                    if handler:
                        for watch in self.observer._watches.copy():  # type: ignore
                            if str(path) in str(watch):
                                self.observer.unschedule(self.observer._watches[watch])  # type: ignore
                    self.logger.info(f"Removed monitoring for {path}")
                except Exception as e:
                    self.logger.error(f"Failed to remove monitoring for {path}: {str(e)}")
            del self.watch_dirs[path]

    async def __aenter__(self):
        """Enter async context."""
        self.start()
        await self.process_existing()
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Exit async context."""
        self.stop()
