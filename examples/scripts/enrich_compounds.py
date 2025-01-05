#!/usr/bin/env python3
"""Example script demonstrating web enrichment functionality.

This script shows how to:
1. Configure web enrichment
2. Create a manager instance
3. Load compounds from various sources
4. Enrich compounds with web data
5. Export enriched data

The script handles both JSON and TSV input formats:

JSON format:
{
    "compounds": [
        {
            "category": "CNS-Active (BBB Permeable)",
            "compounds": [
                {
                    "name": "Caffeine",
                    "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                    "cas_number": "58-08-2"
                }
            ]
        }
    ]
}

TSV format:
name    smiles    cas_number
# CNS-Active (BBB Permeable)
Caffeine    CN1C=NC2=C1C(=O)N(C(=O)N2C)C    58-08-2
"""

import argparse
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from binding_data_processor.models.compound import Compound
from binding_data_processor.web_enrichment.manager import (
    WebEnrichmentManager,
    EnrichmentConfig,
)
from binding_data_processor.pipeline.infrastructure.circuit_breaker import (
    CircuitConfig,
)
from binding_data_processor.pipeline.infrastructure.monitoring import (
    MetricsCollector,
)


@dataclass
class Category:
    """Category of compounds."""

    name: str
    compounds: List[Compound]
    description: Optional[str] = None


def setup_logging(verbose: bool = False) -> logging.Logger:
    """Set up logging configuration.
    
    Args:
        verbose: Whether to enable debug logging
        
    Returns:
        Configured logger
    """
    # Set up handler
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(
        logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
    )

    # Set up logger
    logger = logging.getLogger("web_enrichment")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.addHandler(handler)

    return logger


def parse_category_from_comment(line: str) -> Optional[str]:
    """Parse category from comment line.
    
    Args:
        line: Comment line from TSV file
        
    Returns:
        Category name if found, None otherwise
    """
    line = line.strip().lstrip("#").strip()
    return line if line else None


def create_compound_from_json(data: dict) -> Compound:
    """Create compound from JSON data.
    
    Args:
        data: JSON data for compound
        
    Returns:
        Created compound
    """
    return Compound(
        name=data["name"],
        smiles=data["smiles"],
        cas_number=data.get("cas_number"),
    )


def create_compound_from_tsv(
    fields: List[str],
    name_idx: int,
    smiles_idx: int,
    cas_idx: Optional[int],
) -> Compound:
    """Create compound from TSV fields.
    
    Args:
        fields: TSV fields
        name_idx: Index of name field
        smiles_idx: Index of SMILES field
        cas_idx: Optional index of CAS number field
        
    Returns:
        Created compound
    """
    return Compound(
        name=fields[name_idx],
        smiles=fields[smiles_idx],
        cas_number=fields[cas_idx] if cas_idx is not None else None,
    )


def load_compounds_from_json(path: Path) -> List[Category]:
    """Load compounds from JSON file.
    
    Args:
        path: Path to JSON file
        
    Returns:
        List of categories
    """
    with open(path) as f:
        data = json.load(f)
        categories = []
        
        for category_data in data["compounds"]:
            compounds = [
                create_compound_from_json(item)
                for item in category_data["compounds"]
            ]
            categories.append(
                Category(
                    name=category_data["category"],
                    compounds=compounds,
                    description=category_data.get("description"),
                )
            )
        
        return categories


def load_compounds_from_tsv(path: Path) -> List[Category]:
    """Load compounds from TSV file.
    
    Args:
        path: Path to TSV file
        
    Returns:
        List of categories
    """
    categories = []
    current_category = Category(
        name="Uncategorized",
        compounds=[],
    )
    
    with open(path) as f:
        # Parse header
        header = f.readline().strip().split("\t")
        name_idx = header.index("name")
        smiles_idx = header.index("smiles")
        cas_idx = header.index("cas_number") if "cas_number" in header else None

        # Parse compounds
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith("#"):
                # Try to parse category from comment
                category_name = parse_category_from_comment(line)
                if category_name:
                    if current_category.compounds:
                        categories.append(current_category)
                    current_category = Category(
                        name=category_name,
                        compounds=[],
                    )
                continue

            # Parse compound data
            fields = line.split("\t")
            if len(fields) >= 2:  # Ensure minimum required fields
                current_category.compounds.append(
                    create_compound_from_tsv(
                        fields,
                        name_idx,
                        smiles_idx,
                        cas_idx,
                    )
                )

    # Add final category
    if current_category.compounds:
        categories.append(current_category)

    return categories


def load_compounds(input_path: Path) -> List[Category]:
    """Load compounds from input file.
    
    Args:
        input_path: Path to input file
        
    Returns:
        List of categories containing compounds
        
    Raises:
        ValueError: If input file format is not supported
    """
    if input_path.suffix == ".json":
        return load_compounds_from_json(input_path)
    elif input_path.suffix == ".tsv":
        return load_compounds_from_tsv(input_path)
    else:
        raise ValueError(f"Unsupported file format: {input_path.suffix}")


def save_enriched_data(
    categories: List[Category],
    output_path: Path,
) -> None:
    """Save enriched compound data.
    
    Args:
        categories: List of categories containing compounds
        output_path: Path to output file
    """
    data = {"compounds": []}
    
    for category in categories:
        category_data = {
            "category": category.name,
            "compounds": [],
        }
        if category.description:
            category_data["description"] = category.description
        
        for compound in category.compounds:
            item = {
                "name": compound.name,
                "smiles": compound.smiles,
                "cas_number": compound.cas_number,
                "swiss_data": compound.swiss_data,
                "community_data": compound.community_data,
                "social_data": compound.social_data,
                "enrichment_metadata": compound.enrichment_metadata,
            }
            category_data["compounds"].append(item)
            
        data["compounds"].append(category_data)

    with open(output_path, "w") as f:
        json.dump(data, f, indent=4)


def main():
    """Run web enrichment example."""
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Enrich compounds with web data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "input_path",
        type=Path,
        help="Path to input file (JSON or TSV)",
    )
    parser.add_argument(
        "output_path",
        type=Path,
        help="Path to output JSON file",
    )
    parser.add_argument(
        "--reddit-id",
        help="Reddit client ID",
    )
    parser.add_argument(
        "--reddit-secret",
        help="Reddit client secret",
    )
    parser.add_argument(
        "--twitter-token",
        help="Twitter bearer token",
    )
    parser.add_argument(
        "--skip-predictions",
        action="store_true",
        help="Skip Swiss predictions",
    )
    parser.add_argument(
        "--skip-web-data",
        action="store_true",
        help="Skip web data collection",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Disable caching",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of worker threads",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10,
        help="Batch size for processing",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging",
    )
    args = parser.parse_args()

    # Set up logging
    logger = setup_logging(args.verbose)

    # Create directories
    cache_dir = Path.home() / ".cache" / "binding_data_processor"
    model_dir = cache_dir / "models"
    cache_dir.mkdir(parents=True, exist_ok=True)
    model_dir.mkdir(parents=True, exist_ok=True)

    # Set up metrics collector
    metrics = MetricsCollector(
        namespace="web_enrichment",
        logger=logger,
    )

    # Configure enrichment
    config = EnrichmentConfig(
        reddit_client_id=args.reddit_id,
        reddit_client_secret=args.reddit_secret,
        twitter_bearer_token=args.twitter_token,
        skip_predictions=args.skip_predictions,
        skip_web_data=args.skip_web_data,
        use_cache=not args.no_cache,
        n_workers=args.workers,
        batch_size=args.batch_size,
        model_dir=model_dir,
        cache_dir=cache_dir,
        circuit_config=CircuitConfig(
            failure_threshold=3,
            recovery_timeout=60,
        ),
        client_configs={
            "swiss": {
                "base_url": "https://www.swisstargetprediction.ch/api",
            },
            "community": {
                "base_url": "https://api.psychonautwiki.org",
            },
            "social": {
                "base_url": "https://api.reddit.com",
            },
        },
    )

    # Create manager
    manager = WebEnrichmentManager(
        config=config,
        logger=logger,
        metrics_collector=metrics,
    )

    try:
        # Load compounds
        logger.info(f"Loading compounds from {args.input_path}")
        categories = load_compounds(args.input_path)
        total_compounds = sum(len(cat.compounds) for cat in categories)
        logger.info(f"Loaded {total_compounds} compounds in {len(categories)} categories")

        # Enrich compounds
        for category in categories:
            logger.info(f"Enriching compounds in category: {category.name}")
            compounds = category.compounds
            manager.enrich_compounds(compounds)

            # Add category to enrichment metadata
            for compound in compounds:
                if compound.enrichment_metadata:
                    compound.enrichment_metadata["category"] = category.name

        # Save results
        logger.info(f"Saving enriched data to {args.output_path}")
        save_enriched_data(categories, args.output_path)

        # Log metrics
        metrics = manager.get_metrics()
        logger.info("Enrichment metrics:")
        logger.info(f"- Processed compounds: {metrics['processed_compounds']}")
        logger.info(f"- Failed compounds: {metrics['failed_compounds']}")
        for client, client_metrics in metrics["clients"].items():
            logger.info(f"- {client} metrics:")
            for key, value in client_metrics.items():
                logger.info(f"  - {key}: {value}")

    except Exception as e:
        logger.error(f"Error enriching compounds: {str(e)}")
        raise

    finally:
        # Clean up
        manager.close()


if __name__ == "__main__":
    main()
