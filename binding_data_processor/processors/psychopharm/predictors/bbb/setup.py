#!/usr/bin/env python
"""Setup script for BBB permeability prediction package.

This script:
1. Installs required dependencies
2. Creates necessary directories
3. Downloads pre-trained models
4. Validates the installation
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

import requests
from tqdm import tqdm


def setup_logging(log_level: str = "INFO") -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )


def install_requirements(requirements_file: Path) -> None:
    """Install required Python packages."""
    logging.info("Installing required packages...")
    try:
        subprocess.check_call([
            sys.executable, "-m", "pip", "install", "-r", str(requirements_file)
        ])
    except subprocess.CalledProcessError as e:
        logging.error(f"Error installing requirements: {str(e)}")
        raise


def create_directories(
    base_dir: Path,
    dirs: List[str] = ["models/bbb", "cache", "output"]
) -> None:
    """Create necessary directories."""
    logging.info("Creating directories...")
    for dir_name in dirs:
        dir_path = base_dir / dir_name
        dir_path.mkdir(parents=True, exist_ok=True)
        logging.debug(f"Created directory: {dir_path}")


def download_file(url: str, output_path: Path) -> None:
    """Download file with progress bar."""
    response = requests.get(url, stream=True)
    total_size = int(response.headers.get("content-length", 0))

    with open(output_path, "wb") as f, tqdm(
        desc=output_path.name,
        total=total_size,
        unit="iB",
        unit_scale=True,
    ) as pbar:
        for data in response.iter_content(chunk_size=1024):
            size = f.write(data)
            pbar.update(size)


def download_models(model_dir: Path) -> None:
    """Download pre-trained models."""
    logging.info("Downloading pre-trained models...")
    
    models = {
        # Base models
        "fingerprint_model.pt": "https://example.com/models/fingerprint_model.pt",
        "descriptor_model.pt": "https://example.com/models/descriptor_model.pt",
        
        # Transporter models
        "pgp_model.pt": "https://example.com/models/pgp_model.pt",
        "bcrp_model.pt": "https://example.com/models/bcrp_model.pt",
        
        # ML models
        "abuse_model.pt": "https://example.com/models/abuse_model.pt",
        "toxicity_model.pt": "https://example.com/models/toxicity_model.pt",
        "receptor_model.pt": "https://example.com/models/receptor_model.pt",
    }
    
    for model_name, model_url in models.items():
        output_path = model_dir / model_name
        if not output_path.exists():
            try:
                download_file(model_url, output_path)
            except Exception as e:
                logging.error(f"Error downloading {model_name}: {str(e)}")
                raise


def validate_installation(
    base_dir: Path,
    test_compound_file: Optional[Path] = None,
) -> None:
    """Validate the installation by running predictions on test compounds."""
    logging.info("Validating installation...")
    
    try:
        # Import the package
        from binding_data_processor.processors.psychopharm.predictors.bbb import (
            BBBPredictorWebEnriched
        )
        
        # Initialize predictor
        predictor = BBBPredictorWebEnriched(
            model_dir=str(base_dir / "models/bbb"),
            cache_dir=str(base_dir / "cache"),
        )
        
        # Run test prediction
        if test_compound_file and test_compound_file.exists():
            from binding_data_processor.models.core import CompoundData
            
            # Load test compound
            with open(test_compound_file) as f:
                next(f)  # Skip header
                line = next(f)
                name, smiles, cas = line.strip().split("\t")
                compound = CompoundData(
                    name=name,
                    smiles=smiles,
                    cas_number=cas,
                )
            
            # Make prediction
            result = predictor.predict(compound)
            logging.info(f"Test prediction successful: {result.value}")
        
        logging.info("Installation validated successfully!")
        
    except Exception as e:
        logging.error(f"Validation failed: {str(e)}")
        raise


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Set up BBB permeability prediction package"
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path.cwd(),
        help="Base directory for installation",
    )
    parser.add_argument(
        "--test-file",
        type=Path,
        default=None,
        help="Path to test compounds file",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging level",
    )
    args = parser.parse_args()

    # Set up logging
    setup_logging(args.log_level)

    try:
        # Get package directory
        package_dir = Path(__file__).parent

        # Install requirements
        install_requirements(package_dir / "requirements.txt")

        # Create directories
        create_directories(args.base_dir)

        # Download models
        download_models(args.base_dir / "models/bbb")

        # Validate installation
        validate_installation(args.base_dir, args.test_file)

        logging.info("Setup completed successfully!")

    except Exception as e:
        logging.error(f"Setup failed: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
