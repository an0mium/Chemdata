#!/usr/bin/env python3
"""Example script demonstrating comprehensive patent searching.

This script shows how to:
1. Search patents across Google Patents, USPTO, and Espacenet
2. Get patent family information from INPADOC
3. Get legal status information
4. Extract chemical compounds
5. Analyze patent classifications
6. Export results

Example usage:
    python search_patents.py --query "5-HT2A antagonist" --structure "CC1=CC=C(C=C1)NC(=O)CN2CCN(CC2)CC3=CC=C(C=C3)F"
"""

import argparse
import logging
import json
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any

from binding_data_processor.web_enrichment.clients.patents import (
    EnhancedPatentClient,
    PatentData,
)


def setup_logging(output_dir: Path) -> None:
    """Setup logging configuration.

    Args:
        output_dir: Directory for log files
    """
    # Create logs directory
    log_dir = output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"patent_search_{timestamp}.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),
        ],
    )


def export_results(
    patents: list[PatentData],
    output_dir: Path,
    timestamp: str,
) -> None:
    """Export patent search results.

    Args:
        patents: List of patent data
        output_dir: Output directory
        timestamp: Timestamp string
    """
    # Create results directory
    results_dir = output_dir / timestamp
    results_dir.mkdir(parents=True, exist_ok=True)

    # Export full JSON results
    with open(results_dir / "patents.json", "w") as f:
        json.dump(
            [patent.dict() for patent in patents],
            f,
            indent=2,
        )

    # Export TSV summary
    with open(results_dir / "patents_summary.tsv", "w") as f:
        # Write header
        f.write("patent_number\ttitle\tassignee\tfiling_date\t")
        f.write("family_size\tcountries\tcompounds\tlegal_status\n")

        # Write data
        for patent in patents:
            f.write(f"{patent.patent_number}\t")
            f.write(f"{patent.title}\t")
            f.write(f"{patent.assignee or ''}\t")
            f.write(f"{patent.filing_date or ''}\t")

            # Family info
            if patent.family:
                f.write(f"{len(patent.family.members)}\t")
                f.write(f"{','.join(sorted(patent.family.countries))}\t")
            else:
                f.write("\t\t")

            # Compounds
            f.write(f"{','.join(patent.compounds)}\t")

            # Legal status
            if patent.legal_status:
                f.write(f"{patent.legal_status.status}")
            f.write("\n")


def print_patent_stats(patents: list[PatentData]) -> None:
    """Print statistics about patent results.

    Args:
        patents: List of patent data
    """
    logger = logging.getLogger(__name__)

    # Basic stats
    logger.info(f"\nFound {len(patents)} patents")

    # Source distribution
    sources = set()
    for patent in patents:
        sources.update(patent.metadata["sources"])

    logger.info("\nSources:")
    for source in sorted(sources):
        count = sum(1 for p in patents if source in p.metadata["sources"])
        logger.info(f"- {source}: {count}")

    # Family coverage
    with_family = sum(1 for p in patents if p.family)
    logger.info(f"\nPatents with family data: {with_family}")

    if with_family:
        # Country distribution
        countries = {}
        for p in patents:
            if p.family:
                for country in p.family.countries:
                    countries[country] = countries.get(country, 0) + 1

        logger.info("\nTop countries:")
        for country, count in sorted(
            countries.items(),
            key=lambda x: x[1],
            reverse=True,
        )[:10]:
            logger.info(f"- {country}: {count}")

    # Legal status
    with_status = sum(1 for p in patents if p.legal_status)
    logger.info(f"\nPatents with legal status: {with_status}")

    if with_status:
        # Status distribution
        statuses = {}
        for p in patents:
            if p.legal_status:
                status = p.legal_status.status
                statuses[status] = statuses.get(status, 0) + 1

        logger.info("\nStatus distribution:")
        for status, count in sorted(
            statuses.items(),
            key=lambda x: x[1],
            reverse=True,
        ):
            logger.info(f"- {status}: {count}")

    # Compound extraction
    with_compounds = sum(1 for p in patents if p.compounds)
    total_compounds = sum(len(p.compounds) for p in patents)
    logger.info(f"\nPatents with compounds: {with_compounds}")
    logger.info(f"Total compounds found: {total_compounds}")


async def main():
    """Run patent search example."""
    # Parse arguments
    parser = argparse.ArgumentParser(description="Search patents with enhanced data")
    parser.add_argument(
        "--query",
        type=str,
        required=True,
        help="Search query",
    )
    parser.add_argument(
        "--structure",
        type=str,
        help="Chemical structure (SMILES/InChI)",
    )
    parser.add_argument(
        "--date-from",
        type=str,
        help="Start date (YYYY-MM-DD)",
    )
    parser.add_argument(
        "--date-to",
        type=str,
        help="End date (YYYY-MM-DD)",
    )
    parser.add_argument(
        "--classification",
        type=str,
        help="Patent classification code",
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=100,
        help="Maximum number of results",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results"),
        help="Output directory",
    )
    parser.add_argument(
        "--espacenet-key",
        type=str,
        help="Espacenet API key",
    )
    parser.add_argument(
        "--uspto-key",
        type=str,
        help="USPTO API key",
    )
    args = parser.parse_args()

    # Setup logging
    setup_logging(args.output)
    logger = logging.getLogger(__name__)

    try:
        # Create patent client
        client = EnhancedPatentClient(
            espacenet_api_key=args.espacenet_key,
            uspto_api_key=args.uspto_key,
        )

        # Prepare date range
        date_range = None
        if args.date_from or args.date_to:
            date_range = (
                args.date_from or "1800-01-01",
                args.date_to or datetime.now().strftime("%Y-%m-%d"),
            )

        # Search patents
        logger.info("Searching patents...")
        patents = await client.search_patents(
            query=args.query,
            structure=args.structure,
            date_range=date_range,
            classification=args.classification,
            max_results=args.max_results,
        )

        if not patents:
            logger.info("No patents found")
            return

        # Print statistics
        print_patent_stats(patents)

        # Export results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        export_results(patents, args.output, timestamp)
        logger.info(f"\nResults exported to {args.output / timestamp}")

    except Exception as e:
        logger.error(f"Search failed: {str(e)}", exc_info=True)
        raise


if __name__ == "__main__":
    import asyncio

    asyncio.run(main())
