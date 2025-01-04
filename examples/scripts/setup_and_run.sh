#!/bin/bash
# Example script demonstrating how to set up and run the ChemData pipeline

# Initialize configuration with default settings
echo "Initializing configuration..."
python -m binding_data_processor.cli init \
    --data-dir data \
    --model-dir models \
    --cache-dir ~/.cache/chemdata \
    --log-dir ~/.local/share/chemdata/logs \
    --workers 4 \
    --batch-size 100

# Set API credentials (replace with your actual credentials)
echo "Setting API credentials..."
python -m binding_data_processor.cli credentials \
    --reddit-client-id YOUR_REDDIT_CLIENT_ID \
    --reddit-client-secret YOUR_REDDIT_CLIENT_SECRET \
    --twitter-api-key YOUR_TWITTER_API_KEY \
    --twitter-api-secret YOUR_TWITTER_API_SECRET \
    --discord-token YOUR_DISCORD_TOKEN \
    --bluesky-handle YOUR_BLUESKY_HANDLE \
    --bluesky-password YOUR_BLUESKY_PASSWORD \
    --llm-api-key YOUR_LLM_API_KEY

# Import target patterns
echo "Importing target patterns..."
python -m binding_data_processor.cli patterns --import-file examples/config/target_patterns.json

# Import web sources
echo "Importing web sources..."
python -m binding_data_processor.cli sources --import-file examples/config/web_sources.json

# Import ML config
echo "Importing ML configuration..."
python -m binding_data_processor.cli ml --import-file examples/config/ml_config.json

# Show current configuration
echo "Current configuration:"
python -m binding_data_processor.cli show

# Run the pipeline
echo "Running pipeline..."
python -m binding_data_processor.main \
    --mode both \
    --bindingdb-file data/bindingdb_all.tsv \
    --output-dir output \
    --checkpoint-file checkpoints/pipeline.pkl \
    --host 0.0.0.0 \
    --port 8050
