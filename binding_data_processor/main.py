"""Main entry point for the binding data processing pipeline.

This module provides:
1. Command line interface for running the pipeline
2. Configuration management
3. API key and settings handling
4. Progress tracking and reporting
5. Error handling and logging
6. Documentation and help
"""

import os
import sys
import json
import logging
import argparse
from pathlib import Path
from typing import Dict, Any, Optional
from datetime import datetime

import yaml
from tqdm import tqdm

from binding_data_processor.pipeline import PipelineManager
from binding_data_processor.config import load_config, save_config, DEFAULT_CONFIG
from web.app import ChemDataApp
from logger import LogManager

logger = LogManager().get_logger(__name__)


def setup_logging(log_dir: str = "logs") -> None:
    """Set up logging configuration.

    Args:
        log_dir: Directory for log files
    """
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"pipeline_{timestamp}.log")

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout),
        ],
    )


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Process chemical binding data and enrich with web data"
    )

    # Mode selection
    parser.add_argument(
        "--mode",
        choices=["pipeline", "web", "both"],
        default="both",
        help="Mode to run in (default: both)",
    )

    # Input/output options
    parser.add_argument(
        "--bindingdb-file",
        help="Path to BindingDB data file (optional, will download if not provided)",
    )
    parser.add_argument(
        "--output-dir",
        default="output",
        help="Directory for output files (default: output)",
    )
    parser.add_argument(
        "--checkpoint-file",
        help="Path to checkpoint file to resume from",
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
        "--workers",
        type=int,
        default=4,
        help="Number of worker threads (default: 4)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="Batch size for processing (default: 100)",
    )

    # Target filtering
    parser.add_argument(
        "--target-patterns",
        help="JSON file containing target regex patterns",
    )

    # Configuration
    parser.add_argument(
        "--config",
        help="Path to configuration file",
    )
    parser.add_argument(
        "--save-config",
        help="Save current configuration to file",
    )

    # API keys and credentials
    parser.add_argument("--reddit-client-id", help="Reddit API client ID")
    parser.add_argument("--reddit-client-secret", help="Reddit API client secret")
    parser.add_argument("--twitter-api-key", help="Twitter API key")
    parser.add_argument("--twitter-api-secret", help="Twitter API secret")
    parser.add_argument("--discord-token", help="Discord bot token")
    parser.add_argument("--bluesky-handle", help="Bluesky handle")
    parser.add_argument("--bluesky-password", help="Bluesky password")
    parser.add_argument("--llm-api-key", help="API key for LLM analysis")

    # Web options
    parser.add_argument(
        "--host",
        default="0.0.0.0",
        help="Host to run web server on (default: 0.0.0.0)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8050,
        help="Port to run web server on (default: 8050)",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Run in debug mode",
    )

    return parser.parse_args()


def load_target_patterns(patterns_file: str) -> Optional[Dict[str, str]]:
    """Load target regex patterns from file.

    Args:
        patterns_file: Path to JSON file containing patterns

    Returns:
        Dictionary mapping target types to regex patterns
    """
    try:
        with open(patterns_file) as f:
            patterns = json.load(f)
        logger.info(f"Loaded {len(patterns)} target patterns from {patterns_file}")
        return patterns
    except Exception as e:
        logger.error(f"Error loading target patterns: {e}")
        return None


def run_pipeline(args: argparse.Namespace, config: Dict[str, Any]) -> None:
    """Run the data processing pipeline.

    Args:
        args: Command line arguments
        config: Configuration dictionary
    """
    try:
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Load target patterns if provided
        target_patterns = None
        if args.target_patterns:
            target_patterns = load_target_patterns(args.target_patterns)

        # Initialize pipeline
        pipeline = PipelineManager(
            data_dir=config["data_dir"],
            model_dir=config["model_dir"],
            n_workers=args.workers,
            batch_size=args.batch_size,
            checkpoint_interval=config["checkpoint_interval"],
            max_retries=config["max_retries"],
            cache_dir=config["cache_dir"],
            # API keys from args override config
            reddit_client_id=args.reddit_client_id or config.get("reddit_client_id"),
            reddit_client_secret=args.reddit_client_secret
            or config.get("reddit_client_secret"),
            twitter_api_key=args.twitter_api_key or config.get("twitter_api_key"),
            twitter_api_secret=args.twitter_api_secret
            or config.get("twitter_api_secret"),
            discord_token=args.discord_token or config.get("discord_token"),
            bluesky_handle=args.bluesky_handle or config.get("bluesky_handle"),
            bluesky_password=args.bluesky_password or config.get("bluesky_password"),
            llm_api_key=args.llm_api_key or config.get("llm_api_key"),
        )

        # Run pipeline
        stats = pipeline.run_pipeline(
            bindingdb_file=args.bindingdb_file,
            target_patterns=target_patterns,
            skip_predictions=args.skip_predictions,
            skip_web_data=args.skip_web_data,
            output_dir=output_dir,
            checkpoint_file=args.checkpoint_file,
            use_cache=not args.no_cache,
        )

        # Save final stats
        stats_file = output_dir / "pipeline_stats.json"
        with open(stats_file, "w") as f:
            json.dump(stats, f, indent=2)
        logger.info(f"Saved pipeline stats to {stats_file}")

        # Print summary
        print("\nPipeline Summary:")
        print(f"Total compounds processed: {stats['total_compounds']}")
        print(f"BindingDB compounds: {stats['bindingdb_compounds']}")
        if not args.skip_web_data:
            print(f"Web compounds: {stats['web_compounds']}")
            print(f"Social media compounds: {stats['social_compounds']}")
        if not args.skip_predictions:
            print(f"Compounds with predictions: {stats['with_predictions']}")
        print(f"Compounds with web data: {stats['with_web_data']}")
        print(f"Compounds from cache: {stats['from_cache']}")
        if stats["errors"]:
            print(f"\nErrors encountered: {len(stats['errors'])}")
            for error in stats["errors"]:
                print(f"  - {error}")

        return pipeline

    except KeyboardInterrupt:
        logger.info("Pipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        sys.exit(1)


def run_web(
    args: argparse.Namespace,
    config: Dict[str, Any],
    pipeline: Optional[PipelineManager] = None,
) -> None:
    """Run the web application.

    Args:
        args: Command line arguments
        config: Configuration dictionary
        pipeline: Optional pipeline manager instance
    """
    try:
        # Initialize web app
        app = ChemDataApp(
            data_dir=config["data_dir"],
            model_dir=config["model_dir"],
            debug=args.debug,
            n_workers=args.workers,
            batch_size=args.batch_size,
            checkpoint_interval=config["checkpoint_interval"],
            max_retries=config["max_retries"],
            cache_dir=config["cache_dir"],
            reddit_client_id=args.reddit_client_id or config.get("reddit_client_id"),
            reddit_client_secret=args.reddit_client_secret
            or config.get("reddit_client_secret"),
            twitter_api_key=args.twitter_api_key or config.get("twitter_api_key"),
            twitter_api_secret=args.twitter_api_secret
            or config.get("twitter_api_secret"),
            discord_token=args.discord_token or config.get("discord_token"),
            bluesky_handle=args.bluesky_handle or config.get("bluesky_handle"),
            bluesky_password=args.bluesky_password or config.get("bluesky_password"),
        )

        # Use existing pipeline if provided
        if pipeline:
            app.pipeline = pipeline

        # Run web server
        app.run(
            host=args.host,
            port=args.port,
            debug=args.debug,
        )

    except KeyboardInterrupt:
        logger.info("Web server interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Web server failed: {e}")
        sys.exit(1)


def main() -> None:
    """Main entry point."""
    # Parse arguments
    args = parse_args()

    # Set up logging
    setup_logging()

    # Load configuration
    config = DEFAULT_CONFIG
    if args.config:
        config.update(load_config(args.config))

    # Save configuration if requested
    if args.save_config:
        save_config(config, args.save_config)
        logger.info(f"Saved configuration to {args.save_config}")
        return

    # Run in selected mode
    pipeline = None
    if args.mode in ["pipeline", "both"]:
        pipeline = run_pipeline(args, config)

    if args.mode in ["web", "both"]:
        run_web(args, config, pipeline)


if __name__ == "__main__":
    main()
