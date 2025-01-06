"""Base classes for document processing."""

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Union

from ...models.compound.types import CompoundData
from ...pipeline.infrastructure.monitoring import Monitor


@dataclass
class DocumentMetadata:
    """Metadata extracted from a document."""

    title: str
    authors: List[str]
    date: Optional[str]
    source: str
    document_type: str
    file_path: Path
    file_size: int
    page_count: int
    language: str
    keywords: List[str]
    abstract: Optional[str]
    references: List[str]
    compounds_mentioned: Set[str]
    confidence_score: float


@dataclass
class ProcessingResult:
    """Result of document processing."""

    metadata: DocumentMetadata
    compounds: List[CompoundData]
    extracted_text: str
    error_messages: List[str]
    warning_messages: List[str]
    processing_time: float
    success: bool


class DocumentProcessor(ABC):
    """Base class for document processors."""

    def __init__(
        self,
        storage_dir: Union[str, Path],
        monitor: Optional[Monitor] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize processor.

        Args:
            storage_dir: Directory for storing processed files and data
            monitor: Optional monitoring instance
            logger: Optional logger instance
        """
        self.storage_dir = Path(storage_dir)
        self.storage_dir.mkdir(parents=True, exist_ok=True)

        self.monitor = monitor or Monitor()
        self.logger = logger or logging.getLogger(__name__)

    @abstractmethod
    async def process_file(self, file_path: Union[str, Path]) -> ProcessingResult:
        """Process a single document file.

        Args:
            file_path: Path to document file

        Returns:
            Processing result containing extracted data and metadata
        """
        pass

    @abstractmethod
    async def process_directory(
        self,
        dir_path: Union[str, Path],
        recursive: bool = True,
        file_pattern: str = "*.*",
    ) -> Dict[Path, ProcessingResult]:
        """Process all matching files in a directory.

        Args:
            dir_path: Directory path to process
            recursive: Whether to process subdirectories
            file_pattern: Glob pattern for matching files

        Returns:
            Dictionary mapping file paths to their processing results
        """
        pass

    @abstractmethod
    async def extract_metadata(self, file_path: Union[str, Path]) -> DocumentMetadata:
        """Extract metadata from document.

        Args:
            file_path: Path to document file

        Returns:
            Extracted document metadata
        """
        pass

    @abstractmethod
    async def extract_text(self, file_path: Union[str, Path]) -> str:
        """Extract text content from document.

        Args:
            file_path: Path to document file

        Returns:
            Extracted text content
        """
        pass

    @abstractmethod
    async def extract_compounds(
        self,
        text: str,
        metadata: Optional[DocumentMetadata] = None,
    ) -> List[CompoundData]:
        """Extract compound information from text.

        Args:
            text: Text to extract compounds from
            metadata: Optional document metadata for context

        Returns:
            List of extracted compound data
        """
        pass

    def _validate_file(self, file_path: Union[str, Path]) -> Path:
        """Validate file path and return Path object.

        Args:
            file_path: File path to validate

        Returns:
            Validated Path object

        Raises:
            FileNotFoundError: If file does not exist
            ValueError: If file is not readable or has invalid extension
        """
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")
        if not path.is_file():
            raise ValueError(f"Not a file: {path}")
        if not os.access(path, os.R_OK):
            raise ValueError(f"File not readable: {path}")
        return path

    def _get_matching_files(
        self,
        dir_path: Union[str, Path],
        recursive: bool = True,
        pattern: str = "*.*",
    ) -> List[Path]:
        """Get list of matching files in directory.

        Args:
            dir_path: Directory path to search
            recursive: Whether to search subdirectories
            pattern: Glob pattern for matching files

        Returns:
            List of matching file paths

        Raises:
            NotADirectoryError: If dir_path is not a directory
        """
        path = Path(dir_path)
        if not path.is_dir():
            raise NotADirectoryError(f"Not a directory: {path}")

        if recursive:
            return list(path.rglob(pattern))
        return list(path.glob(pattern))
