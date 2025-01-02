"""Command-line interface for chemical data collection."""

import argparse
import sys
from pathlib import Path
from typing import List, Optional

from config import DATA_SOURCES, LOG_LEVEL
from data_processor import DataProcessor
from logger import LogManager


def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Args:
        args: Optional list of command line arguments
        
    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Chemical compound data collection tool"
    )
    
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='Input CSV file containing compound data'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output CSV file for processed data'
    )
    
    parser.add_argument(
        '-s', '--sources',
        type=str,
        nargs='+',
        choices=list(DATA_SOURCES.keys()),
        help='Data sources to use (default: all enabled sources)'
    )
    
    parser.add_argument(
        '--log-level',
        type=str,
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default=LOG_LEVEL,
        help='Logging level'
    )
    
    parser.add_argument(
        '--clear-cache',
        action='store_true',
        help='Clear cache before processing'
    )
    
    return parser.parse_args(args)


def validate_files(input_file: str, output_file: str) -> None:
    """
    Validate input and output files.
    
    Args:
        input_file: Input file path
        output_file: Output file path
        
    Raises:
        ValueError: If file validation fails
    """
    # Check input file exists
    input_path = Path(input_file)
    if not input_path.exists():
        raise ValueError(f"Input file does not exist: {input_file}")
    if not input_path.is_file():
        raise ValueError(f"Input path is not a file: {input_file}")
        
    # Check input file is CSV
    if input_path.suffix.lower() != '.csv':
        raise ValueError(f"Input file must be CSV: {input_file}")
        
    # Check output directory exists/can be created
    output_path = Path(output_file)
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        raise ValueError(
            f"Cannot create output directory {output_path.parent}: {str(e)}"
        )
        
    # Check output file can be written
    try:
        output_path.touch()
        output_path.unlink()
    except Exception as e:
        raise ValueError(f"Cannot write to output file {output_file}: {str(e)}")


def main(args: Optional[List[str]] = None) -> int:
    """
    Main entry point.
    
    Args:
        args: Optional list of command line arguments
        
    Returns:
        Exit code (0 for success, non-zero for error)
    """
    try:
        # Parse arguments
        parsed_args = parse_args(args)
        
        # Setup logging
        logger = LogManager().get_logger()
        logger.setLevel(parsed_args.log_level)
        
        # Validate files
        validate_files(parsed_args.input, parsed_args.output)
        
        # Initialize processor
        processor = DataProcessor()
        
        # Clear cache if requested
        if parsed_args.clear_cache:
            logger.info("Clearing cache...")
            processor.cache.clear()
        
        # Process compounds
        logger.info(f"Processing compounds from {parsed_args.input}")
        stats = processor.process_file(
            parsed_args.input,
            parsed_args.output,
            parsed_args.sources
        )
        
        # Log cache statistics
        cache_stats = processor.cache.get_cache_stats()
        logger.info(
            f"Cache statistics:\n"
            f"- Total entries: {cache_stats['total_entries']}\n"
            f"- Total size: {cache_stats['total_size_bytes'] / 1024:.1f} KB"
        )
        
        return 0
        
    except Exception as e:
        logger = LogManager().get_logger()
        logger.exception(f"Error: {str(e)}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
