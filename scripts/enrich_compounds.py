#!/usr/bin/env python3
"""Script to enrich compounds with web data.

This script demonstrates how to:
1. Load compounds from BindingDB
2. Enrich them with web data
3. Export the enriched data

Example usage:
    python enrich_compounds.py \
        --input data/bindingdb_compounds.tsv \
        --output data/enriched_compounds.tsv \
        --reddit-id YOUR_REDDIT_CLIENT_ID \
        --reddit-secret YOUR_REDDIT_CLIENT_SECRET \
        --twitter-token YOUR_TWITTER_BEARER_TOKEN \
        --model-dir models \
        --cache-dir cache \
        --n-workers 4
"""

import argparse
import logging
from pathlib import Path
from typing import List

import pandas as pd
from tqdm import tqdm

from binding_data_processor.data_sources.bindingdb import BindingDBSource
from binding_data_processor.web_enrichment.manager import (
    WebEnrichmentManager,
    EnrichmentConfig,
)
from binding_data_processor.models.compound import Compound


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Enrich compounds with web data"
    )

    # Input/output
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input TSV file with compounds",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output TSV file for enriched compounds",
    )

    # API credentials
    parser.add_argument(
        "--reddit-id",
        help="Reddit API client ID",
    )
    parser.add_argument(
        "--reddit-secret",
        help="Reddit API client secret",
    )
    parser.add_argument(
        "--twitter-token",
        help="Twitter API bearer token",
    )

    # Processing options
    parser.add_argument(
        "--skip-predictions",
        action="store_true",
        help="Skip ML predictions",
    )
    parser.add_argument(
        "--skip-web-data",
        action="store_true",
        help="Skip web data enrichment",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Disable caching",
    )
    parser.add_argument(
        "--n-workers",
        type=int,
        default=4,
        help="Number of worker threads",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="Batch size for processing",
    )

    # Directories
    parser.add_argument(
        "--model-dir",
        type=Path,
        help="Directory containing ML models",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help="Directory for caching",
    )

    # Logging
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level",
    )

    return parser.parse_args()


def setup_logging(level: str) -> logging.Logger:
    """Set up logging configuration.
    
    Args:
        level: Logging level
        
    Returns:
        Configured logger
    """
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger(__name__)


def load_compounds(
    input_file: Path,
    logger: logging.Logger,
) -> List[Compound]:
    """Load compounds from BindingDB file.
    
    Args:
        input_file: Input TSV file
        logger: Logger instance
        
    Returns:
        List of loaded compounds
    """
    logger.info(f"Loading compounds from {input_file}")
    
    bindingdb = BindingDBSource(logger=logger)
    
    # Load compounds with progress bar
    compounds = []
    for chunk in tqdm(
        bindingdb.load_compounds(input_file, chunksize=1000),
        desc="Loading compounds",
        unit="chunk",
    ):
        compounds.extend(chunk)
    
    logger.info(f"Loaded {len(compounds)} compounds")
    return compounds


def save_compounds(
    compounds: List[Compound],
    output_file: Path,
    logger: logging.Logger,
) -> None:
    """Save enriched compounds to TSV file.
    
    Args:
        compounds: List of compounds to save
        output_file: Output TSV file
        logger: Logger instance
    """
    logger.info(f"Saving enriched compounds to {output_file}")
    
    # Convert to DataFrame with progress bar
    data = []
    for compound in tqdm(compounds, desc="Converting compounds"):
        row = {
            # Basic info
            "name": compound.name,
            "smiles": compound.smiles,
            "cas_number": compound.cas_number,
            
            # Swiss data
            "swiss_targets": compound.swiss_data.get("targets", []),
            "swiss_adme": compound.swiss_data.get("adme", {}),
            
            # Community data
            "psychonaut": compound.community_data.get("psychonaut", {}),
            "tripsit": compound.community_data.get("tripsit", {}),
            "erowid": compound.community_data.get("erowid", {}),
            
            # Social data
            "reddit_data": compound.social_data.get("reddit", {}),
            "twitter_data": compound.social_data.get("twitter", {}),
            
            # Metadata
            "enrichment_timestamp": compound.enrichment_metadata.get("timestamp"),
            "enrichment_sources": compound.enrichment_metadata.get("sources", []),
        }
        data.append(row)
    
    df = pd.DataFrame(data)
    df.to_csv(output_file, sep="\t", index=False)
    
    logger.info(f"Saved {len(compounds)} compounds")


def main() -> None:
    """Main entry point."""
    # Parse arguments
    args = parse_args()

    # Set up logging
    logger = setup_logging(args.log_level)

    try:
        # Load compounds
        compounds = load_compounds(args.input, logger)

        # Create enrichment config
        config = EnrichmentConfig(
            reddit_client_id=args.reddit_id,
            reddit_client_secret=args.reddit_secret,
            twitter_bearer_token=args.twitter_token,
            skip_predictions=args.skip_predictions,
            skip_web_data=args.skip_web_data,
            use_cache=not args.no_cache,
            n_workers=args.n_workers,
            batch_size=args.batch_size,
            model_dir=args.model_dir,
            cache_dir=args.cache_dir,
        )

        # Create enrichment manager
        manager = WebEnrichmentManager(config, logger=logger)

        try:
            # Enrich compounds
            logger.info("Enriching compounds with web data")
            manager.enrich_compounds(compounds)

            # Save results
            save_compounds(compounds, args.output, logger)

        finally:
            # Clean up
            manager.close()

    except Exception as e:
        logger.error(f"Error enriching compounds: {str(e)}")
        raise


if __name__ == "__main__":
    main()
