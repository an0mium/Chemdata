#!/usr/bin/env python3
"""Example script demonstrating the compound processing pipeline.

This script shows how to:
1. Load compounds from BindingDB
2. Enrich with web data
3. Monitor social media
4. Apply ML predictions
5. Export results

Example usage:
    python process_compounds.py --input bindingdb.tsv --output results/
"""

import argparse
import logging
from pathlib import Path
from datetime import datetime
from collections import Counter

from binding_data_processor.pipeline.processing.pipeline import (
    ProcessingPipeline,
    ProcessingConfig,
)
from binding_data_processor.pipeline.infrastructure import InfrastructureManager


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
    log_file = log_dir / f"processing_{timestamp}.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),
        ],
    )


def main():
    """Run compound processing pipeline."""
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Process compounds from multiple sources"
    )
    parser.add_argument(
        "--input",
        type=str,
        help="Input BindingDB file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results"),
        help="Output directory",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="Batch size for processing",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of worker threads",
    )
    parser.add_argument(
        "--checkpoint-dir",
        type=Path,
        help="Directory for checkpoints",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help="Directory for caching",
    )
    parser.add_argument(
        "--no-ml",
        action="store_true",
        help="Disable ML predictions",
    )
    parser.add_argument(
        "--no-web",
        action="store_true",
        help="Disable web enrichment",
    )
    parser.add_argument(
        "--no-social",
        action="store_true",
        help="Disable social monitoring",
    )
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.output)
    logger = logging.getLogger(__name__)
    
    try:
        # Create configuration
        config = ProcessingConfig(
            data_dir=args.output,
            cache_dir=args.cache_dir,
            batch_size=args.batch_size,
            n_workers=args.workers,
            use_ml_predictions=not args.no_ml,
            use_web_enrichment=not args.no_web,
            use_social_monitoring=not args.no_social,
        )
        
        # Create infrastructure
        infrastructure = InfrastructureManager(
            checkpoint_dir=args.checkpoint_dir,
        )
        
        # Create pipeline
        pipeline = ProcessingPipeline(
            config=config,
            infrastructure=infrastructure,
        )
        
        # Process compounds
        logger.info("Starting compound processing")
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        compounds = pipeline.process_compounds(
            input_file=args.input,
            output_dir=args.output / timestamp,
            checkpoint_key=f"processing_{timestamp}",
        )
        
        # Log basic stats
        logger.info("Processing completed:")
        logger.info(f"- Total compounds: {pipeline.stats['total_compounds']}")
        logger.info(f"- From BindingDB: {pipeline.stats['bindingdb_compounds']}")
        logger.info(f"- From web: {pipeline.stats['web_compounds']}")
        logger.info(f"- From social: {pipeline.stats['social_compounds']}")
        logger.info(f"- With predictions: {pipeline.stats['with_predictions']}")
        logger.info(f"- With web data: {pipeline.stats['with_web_data']}")
        logger.info(f"- Errors: {len(pipeline.stats['errors'])}")
        
        # Log compound statistics
        if compounds:
            # Count sources
            source_counts = Counter(c.source for c in compounds)
            logger.info("\nSource distribution:")
            for source, count in source_counts.most_common():
                logger.info(f"- {source}: {count}")
            
            # Count prediction types
            prediction_counts = Counter()
            for c in compounds:
                if hasattr(c, "activity_predictions") and c.activity_predictions:
                    prediction_counts["activity"] += 1
                if hasattr(c, "toxicity_predictions") and c.toxicity_predictions:
                    prediction_counts["toxicity"] += 1
                if hasattr(c, "abuse_predictions") and c.abuse_predictions:
                    prediction_counts["abuse"] += 1
                if hasattr(c, "bbb_predictions") and c.bbb_predictions:
                    prediction_counts["bbb"] += 1
            
            logger.info("\nPrediction coverage:")
            for pred_type, count in prediction_counts.most_common():
                percentage = count / len(compounds) * 100
                logger.info(f"- {pred_type}: {count} ({percentage:.1f}%)")
        
    except Exception as e:
        logger.error(f"Processing failed: {str(e)}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
