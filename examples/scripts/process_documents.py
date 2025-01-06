"""Example script demonstrating document processing functionality."""

import asyncio
import logging
from pathlib import Path

from fastapi import FastAPI
from uvicorn import Config, Server

from binding_data_processor.processors.document.pdf import PDFProcessor
from binding_data_processor.web.components.document_routes import DirectoryConfig
from binding_data_processor.web.components.document_upload import DocumentUploadComponent


async def main():
    """Run document processing example."""
    # Configure logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # Create FastAPI app
    app = FastAPI(title="Document Processing Example")

    # Setup storage directory
    storage_dir = Path("data/documents")
    storage_dir.mkdir(parents=True, exist_ok=True)

    # Initialize document upload component
    async with DocumentUploadComponent(
        app=app,
        storage_dir=str(storage_dir),
        logger=logger,
    ) as doc_component:
        # Configure initial monitoring
        config = DirectoryConfig(
            path=str(storage_dir),
            patterns={"*.pdf", "*.PDF"},
            recursive=True,
        )
        await doc_component.route_handlers.handle_monitor_config(config)

        # Configure server
        server_config = Config(
            app=app,
            host="127.0.0.1",
            port=8000,
            log_level="info",
        )
        server = Server(server_config)

        # Start server
        logger.info(
            "Server running at http://127.0.0.1:8000\n"
            "Available endpoints:\n"
            "  POST /documents/upload - Upload and process documents\n"
            "  POST /documents/monitor - Configure directory monitoring\n"
            "  DELETE /documents/monitor - Remove directory monitoring\n"
            "  GET /documents/status/{file_path} - Get processing status"
        )
        await server.serve()


def process_single_file():
    """Example of processing a single PDF file."""
    # Configure logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # Initialize processor
    processor = PDFProcessor(
        storage_dir="data/documents/pdf",
        logger=logger,
    )

    # Process file
    async def run():
        try:
            # Process a PDF file
            result = await processor.process_file("path/to/your/file.pdf")

            # Print results
            if result.success:
                logger.info(f"Successfully processed file and found {len(result.compounds)} compounds")
                for compound in result.compounds:
                    logger.info(f"  - {compound.name or compound.smiles}")
            else:
                logger.error(f"Processing failed: {'; '.join(result.error_messages)}")

        except Exception as e:
            logger.error(f"Error processing file: {str(e)}")

    asyncio.run(run())


def monitor_directory():
    """Example of monitoring a directory for new PDFs."""
    # Configure logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # Initialize processor
    processor = PDFProcessor(
        storage_dir="data/documents/pdf",
        logger=logger,
    )

    # Process files
    async def run():
        try:
            # Process all PDFs in a directory
            results = await processor.process_directory(
                "path/to/your/directory",
                recursive=True,
                file_pattern="*.pdf",
            )

            # Print results
            for path, result in results.items():
                if result.success:
                    logger.info(f"Successfully processed {path}")
                    logger.info(f"Found {len(result.compounds)} compounds")
                else:
                    logger.error(f"Failed to process {path}: {'; '.join(result.error_messages)}")

        except Exception as e:
            logger.error(f"Error processing directory: {str(e)}")

    asyncio.run(run())


if __name__ == "__main__":
    # Run web server example
    asyncio.run(main())

    # Uncomment to run file processing example
    # process_single_file()

    # Uncomment to run directory monitoring example
    # monitor_directory()
