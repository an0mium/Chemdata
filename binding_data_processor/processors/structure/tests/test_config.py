"""Test configuration and settings.

Provides:
1. Test environment configuration
2. Shared test constants
3. Test data paths
4. Mock data settings
5. Test logging setup
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

# Test environment settings
TEST_ENV = {
    "debug": os.environ.get("TEST_DEBUG", "0") == "1",
    "verbose": os.environ.get("TEST_VERBOSE", "0") == "1",
    "skip_slow": os.environ.get("TEST_SKIP_SLOW", "0") == "1",
    "skip_network": os.environ.get("TEST_SKIP_NETWORK", "0") == "1",
}

# Test directories
TEST_ROOT = Path(__file__).parent
TEST_DATA_DIR = TEST_ROOT / "data"
TEST_OUTPUT_DIR = TEST_ROOT / "output"
TEST_TEMP_DIR = TEST_ROOT / "temp"

# Ensure directories exist
for directory in [TEST_DATA_DIR, TEST_OUTPUT_DIR, TEST_TEMP_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# Test molecule settings
MOL_SETTINGS = {
    "max_atoms": 100,  # Maximum atoms for test molecules
    "max_conformers": 10,  # Maximum conformers per molecule
    "embed_timeout": 10,  # Seconds before giving up on 3D embedding
    "optimize_timeout": 30,  # Seconds before giving up on optimization
}

# Mock data settings
MOCK_SETTINGS = {
    "random_seed": 42,  # Fixed seed for reproducibility
    "n_test_mols": 10,  # Number of test molecules to generate
    "n_drug_mols": 5,  # Number of drug-like molecules
    "similarity_threshold": 0.7,  # Threshold for similarity comparisons
    "descriptor_ranges": {  # Ranges for mock descriptors
        "MolWt": (100, 500),
        "LogP": (-2, 5),
        "TPSA": (20, 140),
        "HBA": (0, 10),
        "HBD": (0, 5),
        "RotBonds": (0, 10),
        "AromaticRings": (0, 4),
    },
}

# Test logging setup
LOG_SETTINGS = {
    "level": logging.DEBUG if TEST_ENV["debug"] else logging.INFO,
    "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    "date_format": "%Y-%m-%d %H:%M:%S",
}

# Configure logging
logging.basicConfig(**LOG_SETTINGS)
logger = logging.getLogger("structure_tests")


class TestConfig:
    """Test configuration manager."""

    def __init__(self):
        """Initialize test configuration."""
        self.env = TEST_ENV
        self.mol_settings = MOL_SETTINGS.copy()
        self.mock_settings = MOCK_SETTINGS.copy()
        self.logger = logger

    def get_test_path(self, filename: str) -> Path:
        """Get path in test data directory."""
        return TEST_DATA_DIR / filename

    def get_output_path(self, filename: str) -> Path:
        """Get path in test output directory."""
        return TEST_OUTPUT_DIR / filename

    def get_temp_path(self, filename: str) -> Path:
        """Get path in temporary directory."""
        return TEST_TEMP_DIR / filename

    def should_skip_test(self, test_type: str) -> bool:
        """Check if test should be skipped."""
        if test_type == "slow" and self.env["skip_slow"]:
            return True
        if test_type == "network" and self.env["skip_network"]:
            return True
        return False

    def update_mol_settings(self, **kwargs):
        """Update molecule settings."""
        self.mol_settings.update(kwargs)

    def update_mock_settings(self, **kwargs):
        """Update mock data settings."""
        self.mock_settings.update(kwargs)

    def setup_test_logging(self, name: str):
        """Set up logging for specific test."""
        logger = logging.getLogger(f"structure_tests.{name}")
        if self.env["verbose"]:
            logger.setLevel(logging.DEBUG)
        return logger

    def cleanup_test_files(self):
        """Clean up test output and temp files."""
        import shutil

        for directory in [TEST_OUTPUT_DIR, TEST_TEMP_DIR]:
            if directory.exists():
                shutil.rmtree(directory)
                directory.mkdir(parents=True)


# Global test configuration instance
test_config = TestConfig()


# Decorators for test control
def skip_if_slow(func):
    """Skip test if slow tests are disabled."""

    def wrapper(*args, **kwargs):
        if test_config.should_skip_test("slow"):
            return None
        return func(*args, **kwargs)

    return wrapper


def skip_if_network(func):
    """Skip test if network tests are disabled."""

    def wrapper(*args, **kwargs):
        if test_config.should_skip_test("network"):
            return None
        return func(*args, **kwargs)

    return wrapper


def with_test_logging(func):
    """Set up test-specific logging."""

    def wrapper(*args, **kwargs):
        logger = test_config.setup_test_logging(func.__name__)
        return func(*args, logger=logger, **kwargs)

    return wrapper


# Example usage
if __name__ == "__main__":
    # Get test paths
    data_file = test_config.get_test_path("test.mol2")
    output_file = test_config.get_output_path("result.png")
    temp_file = test_config.get_temp_path("temp.txt")

    # Update settings
    test_config.update_mol_settings(max_atoms=200)
    test_config.update_mock_settings(random_seed=123)

    # Use decorators
    @skip_if_slow
    @skip_if_network
    @with_test_logging
    def example_test(logger):
        """Example test with logging."""
        logger.debug("Running example test")
        logger.info("Test completed")

    # Run test
    example_test()

    # Clean up
    test_config.cleanup_test_files()
