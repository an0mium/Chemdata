"""PDF document processor implementation."""

import asyncio
import logging
import os
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Union

import fitz  # PyMuPDF
from rdkit import Chem
from rdkit.Chem import AllChem

from ...models.compound.types import CompoundData
from ...pipeline.infrastructure.monitoring import Monitor
from .base import DocumentMetadata, DocumentProcessor, ProcessingResult


class PDFProcessor(DocumentProcessor):
    """Processor for PDF documents."""

    def __init__(
        self,
        storage_dir: Union[str, Path],
        monitor: Optional[Monitor] = None,
        logger: Optional[logging.Logger] = None,
        **kwargs,
    ):
        """Initialize PDF processor.

        Args:
            storage_dir: Directory for storing processed files and data
            monitor: Optional monitoring instance
            logger: Optional logger instance
            **kwargs: Additional arguments passed to parent
        """
        super().__init__(storage_dir, monitor, logger)
        self.temp_dir = Path(tempfile.mkdtemp(prefix="pdf_processor_"))

    async def process_file(self, file_path: Union[str, Path]) -> ProcessingResult:
        """Process a single PDF file.

        Args:
            file_path: Path to PDF file

        Returns:
            Processing result containing extracted data and metadata

        Raises:
            ValueError: If file is not a PDF
        """
        start_time = datetime.now()
        path = self._validate_file(file_path)

        if path.suffix.lower() != ".pdf":
            raise ValueError(f"Not a PDF file: {path}")

        error_messages = []
        warning_messages = []

        try:
            # Extract metadata and text
            metadata = await self.extract_metadata(path)
            text = await self.extract_text(path)

            # Extract compounds
            compounds = await self.extract_compounds(text, metadata)

            # Update metadata with found compounds
            metadata.compounds_mentioned = {c.name for c in compounds if c.name}

            return ProcessingResult(
                metadata=metadata,
                compounds=compounds,
                extracted_text=text,
                error_messages=error_messages,
                warning_messages=warning_messages,
                processing_time=(datetime.now() - start_time).total_seconds(),
                success=True,
            )

        except Exception as e:
            error_messages.append(f"Error processing PDF: {str(e)}")
            self.logger.error(f"Failed to process PDF {path}: {str(e)}")
            return ProcessingResult(
                metadata=DocumentMetadata(
                    title=path.name,
                    authors=[],
                    date=None,
                    source="pdf",
                    document_type="pdf",
                    file_path=path,
                    file_size=path.stat().st_size,
                    page_count=0,
                    language="unknown",
                    keywords=[],
                    abstract=None,
                    references=[],
                    compounds_mentioned=set(),
                    confidence_score=0.0,
                ),
                compounds=[],
                extracted_text="",
                error_messages=error_messages,
                warning_messages=warning_messages,
                processing_time=(datetime.now() - start_time).total_seconds(),
                success=False,
            )

    async def process_directory(
        self,
        dir_path: Union[str, Path],
        recursive: bool = True,
        file_pattern: str = "*.pdf",
    ) -> Dict[Path, ProcessingResult]:
        """Process all PDF files in directory.

        Args:
            dir_path: Directory path to process
            recursive: Whether to process subdirectories
            file_pattern: Glob pattern for matching files

        Returns:
            Dictionary mapping file paths to their processing results
        """
        results = {}
        files = self._get_matching_files(dir_path, recursive, file_pattern)

        # Process files concurrently
        tasks = [self.process_file(f) for f in files]
        processed = await asyncio.gather(*tasks, return_exceptions=True)

        for file_path, result in zip(files, processed):
            if isinstance(result, Exception):
                self.logger.error(f"Failed to process {file_path}: {str(result)}")
                results[file_path] = ProcessingResult(
                    metadata=DocumentMetadata(
                        title=file_path.name,
                        authors=[],
                        date=None,
                        source="pdf",
                        document_type="pdf",
                        file_path=file_path,
                        file_size=file_path.stat().st_size,
                        page_count=0,
                        language="unknown",
                        keywords=[],
                        abstract=None,
                        references=[],
                        compounds_mentioned=set(),
                        confidence_score=0.0,
                    ),
                    compounds=[],
                    extracted_text="",
                    error_messages=[f"Processing failed: {str(result)}"],
                    warning_messages=[],
                    processing_time=0.0,
                    success=False,
                )
            else:
                results[file_path] = result

        return results

    async def extract_metadata(self, file_path: Union[str, Path]) -> DocumentMetadata:
        """Extract metadata from PDF.

        Args:
            file_path: Path to PDF file

        Returns:
            Extracted document metadata
        """
        path = Path(file_path)
        doc = fitz.open(str(path))

        try:
            # Extract basic metadata
            meta = doc.metadata

            # Get page count
            page_count = len(doc)

            # Extract text from first page for abstract
            first_page = doc[0].get_text() if page_count > 0 else ""
            abstract = first_page[:500] if first_page else None

            # Extract references
            references = []
            for page in doc:
                text = page.get_text()
                # Simple reference extraction - look for lines starting with common patterns
                for line in text.split("\n"):
                    if any(line.strip().startswith(p) for p in ["[", "(", "1.", "1 "]):
                        references.append(line.strip())

            return DocumentMetadata(
                title=meta.get("title", path.stem),
                authors=meta.get("author", "").split(";"),
                date=meta.get("creationDate", None),
                source="pdf",
                document_type="pdf",
                file_path=path,
                file_size=path.stat().st_size,
                page_count=page_count,
                language=meta.get("language", "unknown"),
                keywords=meta.get("keywords", "").split(","),
                abstract=abstract,
                references=references,
                compounds_mentioned=set(),  # Will be populated after compound extraction
                confidence_score=1.0,  # Base confidence for metadata
            )

        finally:
            doc.close()

    async def extract_text(self, file_path: Union[str, Path]) -> str:
        """Extract text content from PDF.

        Args:
            file_path: Path to PDF file

        Returns:
            Extracted text content
        """
        doc = fitz.open(str(file_path))
        try:
            text = ""
            for page in doc:
                text += page.get_text()
            return text
        finally:
            doc.close()

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
        compounds = []

        # Extract SMILES patterns
        smiles_matches = self._extract_smiles(text)
        for smiles in smiles_matches:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Generate 2D coordinates
                    AllChem.Compute2DCoords(mol)

                    compounds.append(
                        CompoundData(
                            smiles=smiles,
                            name=None,  # Name would need to be extracted from context
                            mol=mol,
                            source_doc=str(metadata.file_path) if metadata else None,
                            confidence=0.9,  # High confidence for SMILES
                        )
                    )
            except Exception as e:
                self.logger.warning(f"Failed to process SMILES {smiles}: {str(e)}")

        # Extract IUPAC names
        iupac_matches = self._extract_iupac_names(text)
        for name in iupac_matches:
            try:
                mol = self._name_to_structure(name)
                if mol:
                    compounds.append(
                        CompoundData(
                            smiles=Chem.MolToSmiles(mol),
                            name=name,
                            mol=mol,
                            source_doc=str(metadata.file_path) if metadata else None,
                            confidence=0.7,  # Lower confidence for name conversion
                        )
                    )
            except Exception as e:
                self.logger.warning(f"Failed to process IUPAC name {name}: {str(e)}")

        return compounds

    def _extract_smiles(self, text: str) -> Set[str]:
        """Extract SMILES strings from text.

        Args:
            text: Text to extract from

        Returns:
            Set of unique SMILES strings
        """
        # This is a simplified implementation - would need more robust pattern matching
        smiles = set()

        # Look for common SMILES patterns
        import re

        patterns = [
            r"SMILES[:=]\s*([^\s;]+)",  # SMILES: or SMILES=
            r"\"([CN][^\"]+)\"",  # Quoted strings starting with C or N
            r"\[([^\]]+)\]",  # Bracketed expressions
        ]

        for pattern in patterns:
            matches = re.finditer(pattern, text)
            for match in matches:
                candidate = match.group(1)
                try:
                    # Validate it's a valid SMILES
                    if Chem.MolFromSmiles(candidate):
                        smiles.add(candidate)
                except:
                    continue

        return smiles

    def _extract_iupac_names(self, text: str) -> Set[str]:
        """Extract IUPAC chemical names from text.

        Args:
            text: Text to extract from

        Returns:
            Set of unique IUPAC names
        """
        # This would need a more sophisticated implementation
        # Could use NLP models trained on chemical names
        names = set()

        # Simple pattern matching for common chemical name patterns
        import re

        patterns = [
            r"([0-9]*[A-Z][a-z]*(?:(?:[a-z]|-|[0-9]|\(|\)|\[|\]|{|})+[A-Z][a-z]*)+(?:acid|amine|anol|ane|ate|ene|ide|ine|ol|one|yl|ether|acetate))",
        ]

        for pattern in patterns:
            matches = re.finditer(pattern, text)
            for match in matches:
                names.add(match.group(1))

        return names

    def _name_to_structure(self, name: str) -> Optional[Chem.Mol]:
        """Convert chemical name to structure.

        Args:
            name: Chemical name to convert

        Returns:
            RDKit molecule or None if conversion fails
        """
        # This would need integration with chemical name parsing libraries
        # For now, just try direct SMILES conversion as fallback
        try:
            return Chem.MolFromSmiles(name)
        except:
            return None

    def __del__(self):
        """Clean up temporary files on deletion."""
        import shutil

        try:
            if hasattr(self, "temp_dir") and self.temp_dir.exists():
                shutil.rmtree(self.temp_dir)
        except Exception as e:
            self.logger.error(f"Failed to clean up temp directory: {str(e)}")
