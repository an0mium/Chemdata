"""Command-line interface for managing ChemData configuration.

This module provides:
1. Configuration initialization and management
2. API credential management
3. Target pattern management
4. Web source configuration
5. ML model configuration
6. Command-line utilities
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Optional, Dict, Any, List

from binding_data_processor.config import (
    Config,
    DEFAULT_CONFIG_FILE,
    TARGET_PATTERNS,
    WEB_SOURCES,
    SOCIAL_SOURCES,
    ML_CONFIG,
)


def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration.

    Args:
        verbose: Whether to enable verbose logging
    """
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(message)s",
    )


def init_config(args: argparse.Namespace) -> None:
    """Initialize configuration.

    Args:
        args: Parsed command line arguments
    """
    config = Config(config_file=args.config)

    # Set directories
    if args.data_dir:
        config.data_dir = Path(args.data_dir)
    if args.model_dir:
        config.model_dir = Path(args.model_dir)
    if args.cache_dir:
        config.cache_dir = Path(args.cache_dir)
    if args.log_dir:
        config.log_dir = Path(args.log_dir)

    # Set processing settings
    if args.workers:
        config.settings["workers"] = args.workers
    if args.batch_size:
        config.settings["batch_size"] = args.batch_size
    if args.checkpoint_interval:
        config.settings["checkpoint_interval"] = args.checkpoint_interval
    if args.max_retries:
        config.settings["max_retries"] = args.max_retries
    if args.similarity_threshold:
        config.settings["similarity_threshold"] = args.similarity_threshold
    if args.confidence_threshold:
        config.settings["confidence_threshold"] = args.confidence_threshold

    # Save config
    config.save()
    logging.info(f"Configuration saved to {config.config_file}")


def set_credentials(args: argparse.Namespace) -> None:
    """Set API credentials.

    Args:
        args: Parsed command line arguments
    """
    config = Config(config_file=args.config)

    # Set credentials
    config.set_api_credentials(
        reddit_client_id=args.reddit_client_id,
        reddit_client_secret=args.reddit_client_secret,
        twitter_api_key=args.twitter_api_key,
        twitter_api_secret=args.twitter_api_secret,
        discord_token=args.discord_token,
        bluesky_handle=args.bluesky_handle,
        bluesky_password=args.bluesky_password,
        llm_api_key=args.llm_api_key,
    )
    logging.info("API credentials updated")


def manage_patterns(args: argparse.Namespace) -> None:
    """Manage target patterns.

    Args:
        args: Parsed command line arguments
    """
    config = Config(config_file=args.config)

    if args.list:
        # List patterns
        patterns = config.get_target_patterns()
        print("\nTarget Patterns:")
        for target, pattern in patterns.items():
            print(f"  {target}: {pattern}")

    elif args.import_file:
        # Import patterns from file
        try:
            with open(args.import_file) as f:
                patterns = json.load(f)
            config.config["target_patterns"] = patterns
            config.save()
            logging.info(f"Imported {len(patterns)} patterns from {args.import_file}")
        except Exception as e:
            logging.error(f"Error importing patterns: {e}")

    elif args.export_file:
        # Export patterns to file
        try:
            patterns = config.get_target_patterns()
            with open(args.export_file, "w") as f:
                json.dump(patterns, f, indent=2)
            logging.info(f"Exported {len(patterns)} patterns to {args.export_file}")
        except Exception as e:
            logging.error(f"Error exporting patterns: {e}")

    elif args.add:
        # Add/update pattern
        if not args.pattern:
            logging.error("Pattern required for add operation")
            return
        config.config.setdefault("target_patterns", {})
        config.config["target_patterns"][args.add] = args.pattern
        config.save()
        logging.info(f"Added pattern for {args.add}")

    elif args.remove:
        # Remove pattern
        if args.remove in config.config.get("target_patterns", {}):
            del config.config["target_patterns"][args.remove]
            config.save()
            logging.info(f"Removed pattern for {args.remove}")
        else:
            logging.error(f"Pattern {args.remove} not found")


def manage_sources(args: argparse.Namespace) -> None:
    """Manage web and social media sources.

    Args:
        args: Parsed command line arguments
    """
    config = Config(config_file=args.config)

    if args.list:
        # List sources
        print("\nWeb Sources:")
        for source in config.get_web_sources():
            print(f"  - {source}")

        print("\nSocial Media Sources:")
        for source_type, sources in config.get_social_sources().items():
            print(f"\n  {source_type}:")
            for source in sources:
                print(f"    - {source}")

    elif args.import_file:
        # Import sources from file
        try:
            with open(args.import_file) as f:
                sources = json.load(f)
            config.config["web_sources"] = sources.get("web_sources", [])
            config.config["social_sources"] = sources.get("social_sources", {})
            config.save()
            logging.info(f"Imported sources from {args.import_file}")
        except Exception as e:
            logging.error(f"Error importing sources: {e}")

    elif args.export_file:
        # Export sources to file
        try:
            sources = {
                "web_sources": config.get_web_sources(),
                "social_sources": config.get_social_sources(),
            }
            with open(args.export_file, "w") as f:
                json.dump(sources, f, indent=2)
            logging.info(f"Exported sources to {args.export_file}")
        except Exception as e:
            logging.error(f"Error exporting sources: {e}")

    elif args.add:
        # Add source
        if args.type == "web":
            config.config.setdefault("web_sources", [])
            if args.add not in config.config["web_sources"]:
                config.config["web_sources"].append(args.add)
                config.save()
                logging.info(f"Added web source {args.add}")
        else:
            config.config.setdefault("social_sources", {})
            config.config["social_sources"].setdefault(args.type, [])
            if args.add not in config.config["social_sources"][args.type]:
                config.config["social_sources"][args.type].append(args.add)
                config.save()
                logging.info(f"Added {args.type} source {args.add}")

    elif args.remove:
        # Remove source
        if args.type == "web":
            if args.remove in config.config.get("web_sources", []):
                config.config["web_sources"].remove(args.remove)
                config.save()
                logging.info(f"Removed web source {args.remove}")
        else:
            if args.remove in config.config.get("social_sources", {}).get(
                args.type, []
            ):
                config.config["social_sources"][args.type].remove(args.remove)
                config.save()
                logging.info(f"Removed {args.type} source {args.remove}")


def manage_ml_config(args: argparse.Namespace) -> None:
    """Manage ML model configurations.

    Args:
        args: Parsed command line arguments
    """
    config = Config(config_file=args.config)

    if args.list:
        # List ML configurations
        print("\nML Model Configurations:")
        for model, model_config in config.ml_config.items():
            print(f"\n  {model}:")
            for key, value in model_config.items():
                print(f"    {key}: {value}")

    elif args.import_file:
        # Import ML config from file
        try:
            with open(args.import_file) as f:
                ml_config = json.load(f)
            config.ml_config.update(ml_config)
            config.save()
            logging.info(f"Imported ML config from {args.import_file}")
        except Exception as e:
            logging.error(f"Error importing ML config: {e}")

    elif args.export_file:
        # Export ML config to file
        try:
            with open(args.export_file, "w") as f:
                json.dump(config.ml_config, f, indent=2)
            logging.info(f"Exported ML config to {args.export_file}")
        except Exception as e:
            logging.error(f"Error exporting ML config: {e}")

    elif args.model:
        # Update model config
        if args.param and args.value:
            if args.model not in config.ml_config:
                config.ml_config[args.model] = {}
            try:
                # Convert value to appropriate type
                if isinstance(ML_CONFIG.get(args.model, {}).get(args.param), int):
                    value = int(args.value)
                elif isinstance(ML_CONFIG.get(args.model, {}).get(args.param), float):
                    value = float(args.value)
                elif isinstance(ML_CONFIG.get(args.model, {}).get(args.param), list):
                    value = json.loads(args.value)
                else:
                    value = args.value
                config.ml_config[args.model][args.param] = value
                config.save()
                logging.info(f"Updated {args.model} {args.param} = {value}")
            except Exception as e:
                logging.error(f"Error updating ML config: {e}")
        else:
            logging.error("Both param and value required for model update")


def show_config(args: argparse.Namespace) -> None:
    """Show current configuration.

    Args:
        args: Parsed command line arguments
    """
    config = Config(config_file=args.config)

    print("\nDirectories:")
    print(f"  Data: {config.data_dir}")
    print(f"  Models: {config.model_dir}")
    print(f"  Cache: {config.cache_dir}")
    print(f"  Logs: {config.log_dir}")

    print("\nProcessing Settings:")
    for key, value in config.get_processing_settings().items():
        print(f"  {key}: {value}")

    print("\nAPI Credentials:")
    creds = config.get_api_credentials()
    for key, value in creds.items():
        # Show first/last 4 chars of credentials
        if value:
            masked = f"{value[:4]}...{value[-4:]}"
        else:
            masked = "Not set"
        print(f"  {key}: {masked}")

    print("\nTarget Patterns:")
    for target, pattern in config.get_target_patterns().items():
        print(f"  {target}: {pattern}")

    print("\nWeb Sources:")
    for source in config.get_web_sources():
        print(f"  - {source}")

    print("\nSocial Media Sources:")
    for source_type, sources in config.get_social_sources().items():
        print(f"\n  {source_type}:")
        for source in sources:
            print(f"    - {source}")

    print("\nML Model Configurations:")
    for model, model_config in config.ml_config.items():
        print(f"\n  {model}:")
        for key, value in model_config.items():
            print(f"    {key}: {value}")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="ChemData configuration management",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        default=DEFAULT_CONFIG_FILE,
        help="Path to config file",
    )

    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Init command
    init_parser = subparsers.add_parser(
        "init",
        help="Initialize configuration",
    )
    init_parser.add_argument(
        "--data-dir",
        type=str,
        help="Data directory path",
    )
    init_parser.add_argument(
        "--model-dir",
        type=str,
        help="Model directory path",
    )
    init_parser.add_argument(
        "--cache-dir",
        type=str,
        help="Cache directory path",
    )
    init_parser.add_argument(
        "--log-dir",
        type=str,
        help="Log directory path",
    )
    init_parser.add_argument(
        "--workers",
        type=int,
        help="Number of worker threads",
    )
    init_parser.add_argument(
        "--batch-size",
        type=int,
        help="Batch size for processing",
    )
    init_parser.add_argument(
        "--checkpoint-interval",
        type=int,
        help="Save checkpoint every N compounds",
    )
    init_parser.add_argument(
        "--max-retries",
        type=int,
        help="Maximum number of retries",
    )
    init_parser.add_argument(
        "--similarity-threshold",
        type=float,
        help="Similarity threshold for deduplication",
    )
    init_parser.add_argument(
        "--confidence-threshold",
        type=float,
        help="Confidence threshold for predictions",
    )

    # Credentials command
    creds_parser = subparsers.add_parser(
        "credentials",
        help="Set API credentials",
    )
    creds_parser.add_argument(
        "--reddit-client-id",
        type=str,
        help="Reddit API client ID",
    )
    creds_parser.add_argument(
        "--reddit-client-secret",
        type=str,
        help="Reddit API client secret",
    )
    creds_parser.add_argument(
        "--twitter-api-key",
        type=str,
        help="Twitter API key",
    )
    creds_parser.add_argument(
        "--twitter-api-secret",
        type=str,
        help="Twitter API secret",
    )
    creds_parser.add_argument(
        "--discord-token",
        type=str,
        help="Discord bot token",
    )
    creds_parser.add_argument(
        "--bluesky-handle",
        type=str,
        help="Bluesky handle",
    )
    creds_parser.add_argument(
        "--bluesky-password",
        type=str,
        help="Bluesky password",
    )
    creds_parser.add_argument(
        "--llm-api-key",
        type=str,
        help="API key for LLM analysis",
    )

    # Patterns command
    patterns_parser = subparsers.add_parser(
        "patterns",
        help="Manage target patterns",
    )
    patterns_group = patterns_parser.add_mutually_exclusive_group(required=True)
    patterns_group.add_argument(
        "--list",
        action="store_true",
        help="List current patterns",
    )
    patterns_group.add_argument(
        "--import-file",
        type=str,
        help="Import patterns from JSON file",
    )
    patterns_group.add_argument(
        "--export-file",
        type=str,
        help="Export patterns to JSON file",
    )
    patterns_group.add_argument(
        "--add",
        type=str,
        help="Add/update pattern",
    )
    patterns_group.add_argument(
        "--remove",
        type=str,
        help="Remove pattern",
    )
    patterns_parser.add_argument(
        "--pattern",
        type=str,
        help="Pattern string for add operation",
    )

    # Sources command
    sources_parser = subparsers.add_parser(
        "sources",
        help="Manage web and social media sources",
    )
    sources_group = sources_parser.add_mutually_exclusive_group(required=True)
    sources_group.add_argument(
        "--list",
        action="store_true",
        help="List current sources",
    )
    sources_group.add_argument(
        "--import-file",
        type=str,
        help="Import sources from JSON file",
    )
    sources_group.add_argument(
        "--export-file",
        type=str,
        help="Export sources to JSON file",
    )
    sources_group.add_argument(
        "--add",
        type=str,
        help="Add source",
    )
    sources_group.add_argument(
        "--remove",
        type=str,
        help="Remove source",
    )
    sources_parser.add_argument(
        "--type",
        choices=[
            "web",
            "subreddits",
            "twitter_queries",
            "discord_servers",
            "bluesky_tags",
        ],
        help="Type of source",
    )

    # ML config command
    ml_parser = subparsers.add_parser(
        "ml",
        help="Manage ML model configurations",
    )
    ml_group = ml_parser.add_mutually_exclusive_group(required=True)
    ml_group.add_argument(
        "--list",
        action="store_true",
        help="List current ML configurations",
    )
    ml_group.add_argument(
        "--import-file",
        type=str,
        help="Import ML config from JSON file",
    )
    ml_group.add_argument(
        "--export-file",
        type=str,
        help="Export ML config to JSON file",
    )
    ml_group.add_argument(
        "--model",
        choices=list(ML_CONFIG.keys()),
        help="Update model configuration",
    )
    ml_parser.add_argument(
        "--param",
        type=str,
        help="Parameter to update",
    )
    ml_parser.add_argument(
        "--value",
        type=str,
        help="New parameter value",
    )

    # Show command
    subparsers.add_parser(
        "show",
        help="Show current configuration",
    )

    return parser.parse_args()


def main() -> None:
    """Main entry point."""
    args = parse_args()
    setup_logging(args.verbose)

    try:
        if args.command == "init":
            init_config(args)
        elif args.command == "credentials":
            set_credentials(args)
        elif args.command == "patterns":
            manage_patterns(args)
        elif args.command == "sources":
            manage_sources(args)
        elif args.command == "ml":
            manage_ml_config(args)
        elif args.command == "show":
            show_config(args)
        else:
            print("No command specified. Use -h for help.")
            sys.exit(1)

    except Exception as e:
        logging.error(f"Error: {str(e)}")
        if args.verbose:
            logging.exception("Detailed error:")
        sys.exit(1)


if __name__ == "__main__":
    main()
