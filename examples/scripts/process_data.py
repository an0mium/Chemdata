"""Example script demonstrating how to use the ChemData pipeline in Python."""

import logging
from pathlib import Path

from binding_data_processor.config import Config
from binding_data_processor.pipeline import PipelineManager
from web.app import ChemDataApp

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def main():
    """Run the ChemData pipeline."""
    try:
        # Initialize configuration
        config = Config()

        # Set up directories
        data_dir = Path("data")
        model_dir = Path("models")
        output_dir = Path("output")
        checkpoint_dir = Path("checkpoints")

        for directory in [data_dir, model_dir, output_dir, checkpoint_dir]:
            directory.mkdir(parents=True, exist_ok=True)

        # Initialize pipeline
        pipeline = PipelineManager(
            data_dir=str(data_dir),
            model_dir=str(model_dir),
            n_workers=4,
            batch_size=100,
            checkpoint_interval=1000,
            max_retries=3,
            cache_dir=str(Path("~/.cache/chemdata").expanduser()),
            # Add your API credentials here
            reddit_client_id="YOUR_REDDIT_CLIENT_ID",
            reddit_client_secret="YOUR_REDDIT_CLIENT_SECRET",
            twitter_api_key="YOUR_TWITTER_API_KEY",
            twitter_api_secret="YOUR_TWITTER_API_SECRET",
            discord_token="YOUR_DISCORD_TOKEN",
            bluesky_handle="YOUR_BLUESKY_HANDLE",
            bluesky_password="YOUR_BLUESKY_PASSWORD",
            llm_api_key="YOUR_LLM_API_KEY",
        )

        # Define target patterns for filtering compounds
        target_patterns = {
            "serotonin": r"5-HT\d*[A-Z]?|serotonin|SLC6A4|tryptamine",
            "dopamine": r"D\d+|dopamine|DAT|SLC6A3|phenethylamine",
            "norepinephrine": r"norepinephrine|NET|SLC6A2|adrenergic",
            "gaba": r"GABA[A-Z]?\d*|SLC6A1|benzodiazepine",
            "glutamate": r"glutamate|NMDA|AMPA|mGluR|GluR|dissociative",
        }

        # Run pipeline
        logger.info("Running pipeline...")
        stats = pipeline.run_pipeline(
            bindingdb_file=str(data_dir / "bindingdb_all.tsv"),
            target_patterns=target_patterns,
            skip_predictions=False,
            skip_web_data=False,
            output_dir=output_dir,
            checkpoint_file=str(checkpoint_dir / "pipeline.pkl"),
            use_cache=True,
        )

        # Print pipeline statistics
        logger.info("Pipeline Statistics:")
        logger.info(f"Total compounds: {stats['total_compounds']}")
        logger.info(f"BindingDB compounds: {stats['bindingdb_compounds']}")
        logger.info(f"Web compounds: {stats['web_compounds']}")
        logger.info(f"Social media compounds: {stats['social_compounds']}")
        logger.info(f"With predictions: {stats['with_predictions']}")
        logger.info(f"With web data: {stats['with_web_data']}")
        logger.info(f"From cache: {stats['from_cache']}")

        if stats["errors"]:
            logger.warning(f"Errors encountered: {len(stats['errors'])}")
            for error in stats["errors"]:
                logger.warning(f"  - {error}")

        # Initialize web app
        logger.info("Starting web interface...")
        app = ChemDataApp(
            data_dir=str(data_dir),
            model_dir=str(model_dir),
            debug=True,
            n_workers=4,
            batch_size=100,
            checkpoint_interval=1000,
            max_retries=3,
            cache_dir=str(Path("~/.cache/chemdata").expanduser()),
            # Use same API credentials as pipeline
            reddit_client_id="YOUR_REDDIT_CLIENT_ID",
            reddit_client_secret="YOUR_REDDIT_CLIENT_SECRET",
            twitter_api_key="YOUR_TWITTER_API_KEY",
            twitter_api_secret="YOUR_TWITTER_API_SECRET",
            discord_token="YOUR_DISCORD_TOKEN",
            bluesky_handle="YOUR_BLUESKY_HANDLE",
            bluesky_password="YOUR_BLUESKY_PASSWORD",
        )

        # Use existing pipeline in web app
        app.pipeline = pipeline

        # Run web server
        app.run(
            host="0.0.0.0",
            port=8050,
            debug=True,
        )

    except KeyboardInterrupt:
        logger.info("Process interrupted by user")
    except Exception as e:
        logger.error(f"Error running pipeline: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
