"""Example script demonstrating Bluelight scraping and analysis workflow."""

import asyncio
import json
import logging
from datetime import datetime
from pathlib import Path

from binding_data_processor.web_enrichment.clients.bluelight import BluelightCrawl4AIClient
from binding_data_processor.web_enrichment.storage.bluelight_storage import BluelightStorage
from binding_data_processor.processors.psychopharm.predictors.bbb.predictors import BBBPredictor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("bluelight_search.log"),
    ],
)
logger = logging.getLogger(__name__)


async def main():
    """Run Bluelight search and analysis workflow."""
    # Initialize components
    storage_dir = Path("data/bluelight")
    storage_dir.mkdir(parents=True, exist_ok=True)

    client = BluelightCrawl4AIClient()
    storage = BluelightStorage(storage_dir=str(storage_dir))
    storage.bbb_predictor = BBBPredictor()

    try:
        # Search for compounds
        compounds = ["XYZ", "ABC"]  # Example compounds
        subforums = ["advanced-drug-discussion", "psychedelic-drugs"]

        for compound in compounds:
            logger.info(f"Searching for compound: {compound}")

            # Search posts
            posts = await client.search_posts(
                query=compound,
                subforums=subforums,
                max_results=5,
            )

            logger.info(f"Found {len(posts)} posts for {compound}")

            # Store posts and create alerts
            for post in posts:
                try:
                    # Extract post ID from URL
                    post_id = post.url.split("/")[-1]

                    # Store post
                    storage.store_post(post_id, post)
                    logger.info(f"Stored post {post_id}")

                    # Create alert if safety notes present
                    if post.safety:
                        # Determine severity based on safety notes
                        safety_text = " ".join(post.safety).lower()
                        severity = "high" if "high risk" in safety_text else "medium"

                        # Create alert
                        alert = storage.create_alert(
                            post_id=post_id,
                            compounds=post.compounds,
                            safety_notes=post.safety,
                            severity=severity,
                        )
                        logger.info(f"Created alert for post {post_id}: {alert.severity}")

                except Exception as e:
                    logger.error(f"Error processing post {post.url}: {e}")
                    storage.store_error(str(e), {"url": post.url})

        # Generate reports
        logger.info("Generating reports...")

        # Get trending compounds
        trending = storage.get_trending_compounds(days=7)
        trending_report = {
            "generated_at": datetime.now().isoformat(),
            "compounds": trending,
        }

        trending_file = storage_dir / "trending_compounds.json"
        with open(trending_file, "w") as f:
            json.dump(trending_report, f, indent=2)
        logger.info(f"Saved trending compounds to {trending_file}")

        # Get safety summary
        safety = storage.get_safety_summary()
        safety_report = {
            "generated_at": datetime.now().isoformat(),
            "summary": safety,
        }

        safety_file = storage_dir / "safety_summary.json"
        with open(safety_file, "w") as f:
            json.dump(safety_report, f, indent=2)
        logger.info(f"Saved safety summary to {safety_file}")

        # Get error summary
        errors = storage.get_error_summary(days=7)
        error_report = {
            "generated_at": datetime.now().isoformat(),
            "summary": errors,
        }

        error_file = storage_dir / "error_summary.json"
        with open(error_file, "w") as f:
            json.dump(error_report, f, indent=2)
        logger.info(f"Saved error summary to {error_file}")

        # Print summary
        print("\nSearch Results Summary:")
        print(f"Total posts processed: {len(storage.posts)}")
        print(f"Total alerts created: {len(storage.alerts)}")
        print(f"Total errors encountered: {len(storage.errors)}")
        print(f"Trending compounds: {len(trending)}")
        print(f"\nReports saved to: {storage_dir}")

    except Exception as e:
        logger.error(f"Fatal error: {e}")
        raise

    finally:
        # Save any remaining data
        storage.save()


if __name__ == "__main__":
    asyncio.run(main())
