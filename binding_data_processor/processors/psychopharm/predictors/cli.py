"""Command-line interface for web application.

This module provides functionality to:
1. Parse command-line arguments
2. Configure web application
3. Load compound data
4. Start services
5. Handle errors
"""

import logging
import argparse
import asyncio
import sys
from pathlib import Path
from typing import List, Optional

from .data_enrichment import EnrichedData
from .web_app import WebApp, WebAppResult


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Start web application for compound data"
    )
    
    # Server options
    parser.add_argument(
        "--host",
        default="localhost",
        help="Server host (default: localhost)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8000,
        help="Server port (default: 8000)",
    )
    
    # Data options
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input data file (TSV format)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory for generated files",
    )
    
    # Logging options
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Logging level (default: INFO)",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        help="Log file path",
    )
    
    # Development options
    parser.add_argument(
        "--dev",
        action="store_true",
        help="Enable development mode",
    )
    parser.add_argument(
        "--reload",
        action="store_true",
        help="Enable auto-reload on file changes",
    )
    
    return parser.parse_args()


def setup_logging(
    level: str,
    log_file: Optional[Path] = None,
) -> None:
    """Set up logging configuration."""
    handlers = []
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        )
    )
    handlers.append(console_handler)
    
    # File handler
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(
            logging.Formatter(
                "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
            )
        )
        handlers.append(file_handler)
    
    # Configure logging
    logging.basicConfig(
        level=getattr(logging, level),
        handlers=handlers,
    )


def load_compounds(
    input_file: Path,
    logger: logging.Logger,
) -> List[EnrichedData]:
    """Load compound data from file."""
    logger.info(f"Loading compounds from {input_file}")
    
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    try:
        # Load compounds from TSV file
        compounds = []
        with open(input_file) as f:
            # Skip header
            header = next(f).strip().split("\t")
            
            # Read compounds
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) != len(header):
                    logger.warning(
                        f"Skipping line with wrong number of fields: {line}"
                    )
                    continue
                
                # Create compound data
                data = dict(zip(header, fields))
                compound = EnrichedData.from_dict(data)
                compounds.append(compound)
        
        logger.info(f"Loaded {len(compounds)} compounds")
        return compounds
        
    except Exception as e:
        logger.error(
            f"Error loading compounds: {str(e)}",
            exc_info=True
        )
        raise


async def start_app(
    args: argparse.Namespace,
    logger: logging.Logger,
) -> WebAppResult:
    """Start web application."""
    try:
        # Load compounds
        compounds = load_compounds(args.input, logger)
        
        # Create web application
        app = WebApp(
            host=args.host,
            port=args.port,
            output_dir=args.output_dir,
            log_level=logger.level,
        )
        
        # Start application
        result = await app.start_app(compounds)
        
        if result.is_valid:
            logger.info(
                f"Web application started at http://{args.host}:{args.port}"
            )
            
            # Print statistics
            logger.info("Application statistics:")
            for key, value in result.stats.items():
                logger.info(f"  {key}: {value}")
            
            return result
            
        else:
            logger.error("Failed to start web application:")
            for issue in result.issues:
                logger.error(f"  {issue}")
            raise RuntimeError("Web application failed to start")
        
    except Exception as e:
        logger.error(
            f"Error starting web application: {str(e)}",
            exc_info=True
        )
        raise


def main() -> None:
    """Main entry point."""
    try:
        # Parse arguments
        args = parse_args()
        
        # Set up logging
        setup_logging(args.log_level, args.log_file)
        logger = logging.getLogger(__name__)
        
        # Start application
        if args.dev:
            # Development mode with auto-reload
            import aiohttp_autoreload
            
            async def run_app():
                result = await start_app(args, logger)
                return result
            
            aiohttp_autoreload.start(
                run_app,
                patterns=["*.py", "*.html", "*.css", "*.js"],
            )
            
        else:
            # Production mode
            asyncio.run(start_app(args, logger))
        
    except KeyboardInterrupt:
        logger.info("Shutting down...")
        sys.exit(0)
        
    except Exception as e:
        logger.error(
            f"Fatal error: {str(e)}",
            exc_info=True
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
