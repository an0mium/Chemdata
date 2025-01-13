"""Configuration management for the binding data processor.

This module provides:
1. Configuration file management
2. Default settings and paths
3. API credential management
4. Target pattern management
5. Web source configuration
6. ML model configuration
"""

import os
import json
import logging
from pathlib import Path
from typing import Dict, Any, Optional, List, Union

# Default paths
DEFAULT_DATA_DIR = Path("../data")
DEFAULT_MODEL_DIR = Path("../models")
DEFAULT_CONFIG_FILE = Path("~/.config/chemdata/config.json").expanduser()

# Cache settings
CACHE_DIR = Path("~/.cache/chemdata").expanduser()
CACHE_EXPIRY = 3600  # 1 hour in seconds

# Default processing settings
DEFAULT_SETTINGS = {
    "workers": 4,
    "batch_size": 100,
    "checkpoint_interval": 1000,
    "max_retries": 3,
    "cache_dir": "~/.cache/chemdata",
    "log_dir": "~/.local/share/chemdata/logs",
    "similarity_threshold": 0.7,
    "confidence_threshold": 0.6,
}

# Target patterns for filtering compounds
TARGET_PATTERNS = {
    # Neurotransmitter systems
    "serotonin": r"5-HT\d*[A-Z]?|serotonin|SLC6A4|tryptamine",
    "dopamine": r"D\d+|dopamine|DAT|SLC6A3|phenethylamine",
    "norepinephrine": r"norepinephrine|NET|SLC6A2|adrenergic",
    "gaba": r"GABA[A-Z]?\d*|SLC6A1|benzodiazepine",
    "glutamate": r"glutamate|NMDA|AMPA|mGluR|GluR|dissociative",
    "acetylcholine": r"acetylcholine|nACh|mACh|cholinergic",
    # Receptor systems
    "opioid": r"[μμκδ]?-?opioid|MOR|KOR|DOR|opiate",
    "cannabinoid": r"cannabinoid|CB[12]|endocannabinoid",
    "sigma": r"sigma[12]?|σ[12]?",
    "imidazoline": r"imidazoline|I[123]",
    # Compound classes
    "psychedelic": r"psychedelic|hallucinogen|entheogen|lysergamide|tryptamine",
    "dissociative": r"dissociative|NMDA|ketamine|PCP|arylcyclohexylamine",
    "empathogen": r"empathogen|entactogen|MDxx|phenethylamine",
    "stimulant": r"stimulant|amphetamine|cathinone|cocaine",
    "depressant": r"depressant|sedative|hypnotic|barbiturate",
    "nootropic": r"nootropic|cognitive enhancer|smart drug|racetam",
    "antidepressant": r"antidepressant|SSRI|SNRI|MAOI|tricyclic",
    "anxiolytic": r"anxiolytic|anxiolysis|anti-anxiety",
    "antipsychotic": r"antipsychotic|neuroleptic|major tranquilizer",
}

# Web sources to search
WEB_SOURCES = [
    # Chemical databases
    "pubchem",
    "chembl",
    "drugbank",
    "bindingdb",
    "chemspider",
    "zinc",
    # Scientific sources
    "wikipedia",
    "pubmed",
    "patents",
    "clinicaltrials",
    # Community sources
    "erowid",
    "psychonautwiki",
    "tripsit",
    "bluelight",
    "drugs-forum",
]

# Social media sources
SOCIAL_SOURCES = {
    # Reddit
    "subreddits": [
        "Drugs",
        "researchchemicals",
        "DrugNerds",
        "Nootropics",
        "psychopharmacology",
        "Psychedelics",
        "microdosing",
        "RationalPsychonaut",
        "ReagentTesting",
        "AskDrugNerds",
    ],
    # Twitter/X
    "twitter_queries": [
        "new psychedelic",
        "new antidepressant",
        "research chemical",
        "novel compound",
        "nootropic",
        "receptor binding",
        "pharmacology",
        "drug discovery",
        "medicinal chemistry",
        "psychopharmacology",
    ],
    # Discord servers
    "discord_servers": [
        "PsychonautWiki",
        "DrugNerds",
        "Neuroscience",
        "MedicinalChemistry",
    ],
    # Bluesky
    "bluesky_tags": [
        "psychedelics",
        "pharmacology",
        "chemistry",
        "neuroscience",
        "drugdiscovery",
    ],
}

# ML model configurations
ML_CONFIG = {
    "activity": {
        "model_type": "ensemble",
        "input_size": 2048,
        "hidden_sizes": [1024, 512, 256],
        "output_size": 128,
        "dropout": 0.2,
    },
    "toxicity": {
        "model_type": "ensemble",
        "input_size": 2048,
        "hidden_sizes": [1024, 512, 256],
        "output_size": 64,
        "dropout": 0.3,
    },
    "abuse": {
        "model_type": "ensemble",
        "input_size": 2048,
        "hidden_sizes": [1024, 512, 256],
        "output_size": 32,
        "dropout": 0.3,
    },
    "binding": {
        "model_type": "gnn",
        "node_features": 64,
        "edge_features": 32,
        "hidden_size": 128,
        "num_layers": 4,
        "dropout": 0.2,
    },
}


class Config:
    """Configuration manager."""

    def __init__(
        self,
        config_file: Optional[str] = None,
        data_dir: Optional[str] = None,
        model_dir: Optional[str] = None,
    ):
        """Initialize configuration.

        Args:
            config_file: Optional path to config file
            data_dir: Optional data directory path
            model_dir: Optional model directory path
        """
        # Initialize logger
        self.logger = logging.getLogger(__name__)

        # Load config file
        self.config_file = Path(config_file or DEFAULT_CONFIG_FILE)
        self.config = self._load_config()

        # Set directories
        self.data_dir = Path(data_dir or self.config.get("data_dir") or DEFAULT_DATA_DIR)
        self.model_dir = Path(model_dir or self.config.get("model_dir") or DEFAULT_MODEL_DIR)
        self.cache_dir = Path(self.config.get("cache_dir") or DEFAULT_SETTINGS["cache_dir"]).expanduser()
        self.log_dir = Path(self.config.get("log_dir") or DEFAULT_SETTINGS["log_dir"]).expanduser()

        # Ensure directories exist
        self._ensure_dirs()

        # Load settings
        self.settings = {
            **DEFAULT_SETTINGS,
            **self.config.get("settings", {}),
        }

        # Load credentials
        self.credentials = self.config.get("credentials", {})

        # Load ML config
        self.ml_config = {
            **ML_CONFIG,
            **self.config.get("ml_config", {}),
        }

        # Validate configuration
        self._validate_config()

    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from file.

        Returns:
            Configuration dictionary
        """
        if self.config_file.exists():
            try:
                with open(self.config_file) as f:
                    return json.load(f)
            except Exception as e:
                self.logger.error(f"Error loading config file: {e}")
        return {}

    def _ensure_dirs(self) -> None:
        """Ensure required directories exist."""
        for path in [
            self.data_dir,
            self.model_dir,
            self.cache_dir,
            self.log_dir,
        ]:
            os.makedirs(path, exist_ok=True)

    def _validate_config(self) -> None:
        """Validate configuration settings."""
        # Validate settings
        for key, value in self.settings.items():
            if key in DEFAULT_SETTINGS:
                expected_type = type(DEFAULT_SETTINGS[key])
                if not isinstance(value, expected_type):
                    self.logger.warning(f"Invalid type for setting {key}: expected {expected_type}, got {type(value)}")
                    self.settings[key] = DEFAULT_SETTINGS[key]

        # Validate credentials
        for key in self.credentials:
            if not isinstance(self.credentials[key], str):
                self.logger.warning(f"Invalid credential value for {key}")
                self.credentials[key] = ""

        # Validate ML config
        for model, config in self.ml_config.items():
            if model in ML_CONFIG:
                for key, value in ML_CONFIG[model].items():
                    if key not in config:
                        self.logger.warning(f"Missing ML config key {key} for {model}")
                        config[key] = value
                    elif not isinstance(config[key], type(value)):
                        self.logger.warning(f"Invalid type for ML config {key} in {model}")
                        config[key] = value

    def save(self) -> None:
        """Save configuration to file."""
        try:
            # Ensure config directory exists
            os.makedirs(self.config_file.parent, exist_ok=True)

            # Save config
            with open(self.config_file, "w") as f:
                json.dump(
                    {
                        "data_dir": str(self.data_dir),
                        "model_dir": str(self.model_dir),
                        "cache_dir": str(self.cache_dir),
                        "log_dir": str(self.log_dir),
                        "settings": self.settings,
                        "credentials": self.credentials,
                        "ml_config": self.ml_config,
                    },
                    f,
                    indent=2,
                )

            self.logger.info(f"Saved configuration to {self.config_file}")

        except Exception as e:
            self.logger.error(f"Error saving configuration: {e}")

    def get_api_credentials(self) -> Dict[str, str]:
        """Get API credentials.

        Returns:
            Dictionary containing API credentials
        """
        return {
            "reddit_client_id": self.credentials.get("reddit_client_id", ""),
            "reddit_client_secret": self.credentials.get("reddit_client_secret", ""),
            "twitter_api_key": self.credentials.get("twitter_api_key", ""),
            "twitter_api_secret": self.credentials.get("twitter_api_secret", ""),
            "discord_token": self.credentials.get("discord_token", ""),
            "bluesky_handle": self.credentials.get("bluesky_handle", ""),
            "bluesky_password": self.credentials.get("bluesky_password", ""),
            "llm_api_key": self.credentials.get("llm_api_key", ""),
        }

    def set_api_credentials(
        self,
        reddit_client_id: Optional[str] = None,
        reddit_client_secret: Optional[str] = None,
        twitter_api_key: Optional[str] = None,
        twitter_api_secret: Optional[str] = None,
        discord_token: Optional[str] = None,
        bluesky_handle: Optional[str] = None,
        bluesky_password: Optional[str] = None,
        llm_api_key: Optional[str] = None,
    ) -> None:
        """Set API credentials.

        Args:
            reddit_client_id: Optional Reddit API client ID
            reddit_client_secret: Optional Reddit API client secret
            twitter_api_key: Optional Twitter API key
            twitter_api_secret: Optional Twitter API secret
            discord_token: Optional Discord bot token
            bluesky_handle: Optional Bluesky handle
            bluesky_password: Optional Bluesky password
            llm_api_key: Optional API key for LLM analysis
        """
        if reddit_client_id:
            self.credentials["reddit_client_id"] = reddit_client_id
        if reddit_client_secret:
            self.credentials["reddit_client_secret"] = reddit_client_secret
        if twitter_api_key:
            self.credentials["twitter_api_key"] = twitter_api_key
        if twitter_api_secret:
            self.credentials["twitter_api_secret"] = twitter_api_secret
        if discord_token:
            self.credentials["discord_token"] = discord_token
        if bluesky_handle:
            self.credentials["bluesky_handle"] = bluesky_handle
        if bluesky_password:
            self.credentials["bluesky_password"] = bluesky_password
        if llm_api_key:
            self.credentials["llm_api_key"] = llm_api_key
        self.save()

    def get_target_patterns(self) -> Dict[str, str]:
        """Get target regex patterns.

        Returns:
            Dictionary mapping target types to regex patterns
        """
        return {
            **TARGET_PATTERNS,
            **self.config.get("target_patterns", {}),
        }

    def get_web_sources(self) -> List[str]:
        """Get web sources to search.

        Returns:
            List of web source names
        """
        return list(set(WEB_SOURCES + self.config.get("web_sources", [])))

    def get_social_sources(self) -> Dict[str, List[str]]:
        """Get social media sources.

        Returns:
            Dictionary mapping source types to lists of sources
        """
        config_sources = self.config.get("social_sources", {})
        return {source_type: list(set(sources + config_sources.get(source_type, []))) for source_type, sources in SOCIAL_SOURCES.items()}

    def get_ml_config(self, model_type: str) -> Dict[str, Any]:
        """Get ML model configuration.

        Args:
            model_type: Type of ML model

        Returns:
            Model configuration dictionary
        """
        if model_type in self.ml_config:
            return self.ml_config[model_type]
        self.logger.warning(f"Unknown ML model type: {model_type}")
        return {}

    def get_processing_settings(self) -> Dict[str, Any]:
        """Get processing settings.

        Returns:
            Dictionary containing processing settings
        """
        return self.settings


# Global config instance
config = Config()
