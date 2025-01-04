#!/usr/bin/env python
"""Script to demonstrate BBB permeability prediction.

This script shows how to use the BBB prediction package to:
1. Load compounds from a TSV file
2. Make BBB permeability predictions
3. Enrich predictions with web data
4. Export results to TSV
"""

import argparse
import logging
from pathlib import Path
from typing import List, Optional

from binding_data_processor.models.core import CompoundData
from binding_data_processor.processors.psychopharm.predictors.bbb import (
    BBBPredictorWebEnriched
)


def setup_logging(log_level: str = "INFO") -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )


def load_compounds(input_file: Path) -> List[CompoundData]:
    """Load compounds from TSV file.
    
    Args:
        input_file: Path to TSV file containing compounds
        
    Returns:
        List of CompoundData objects
    """
    compounds = []
    with open(input_file) as f:
        # Skip header
        next(f)
        
        # Parse compounds
        for line in f:
            name, smiles, cas = line.strip().split("\t")
            compound = CompoundData(
                name=name,
                smiles=smiles,
                cas_number=cas,
            )
            compounds.append(compound)
            
    return compounds


def predict_bbb_permeability(
    compounds: List[CompoundData],
    model_dir: Optional[Path] = None,
    cache_dir: Optional[Path] = None,
) -> None:
    """Make BBB permeability predictions for compounds.
    
    Args:
        compounds: List of compounds to predict
        model_dir: Optional directory containing trained models
        cache_dir: Optional directory for caching
    """
    # Initialize predictor
    predictor = BBBPredictorWebEnriched(
        model_dir=str(model_dir) if model_dir else None,
        cache_dir=str(cache_dir) if cache_dir else None,
    )
    
    # Make predictions
    for compound in compounds:
        logging.info(f"Predicting BBB permeability for {compound.name}")
        
        result = predictor.predict(compound)
        
        # Log results
        logging.info(f"BBB Class: {result.value}")
        logging.info(f"Confidence: {result.confidence:.2f}")
        
        # Log supporting data
        logging.debug("Supporting Data:")
        for key, value in result.supporting_data.items():
            logging.debug(f"  {key}: {value}")
            
    return predictor


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Predict BBB permeability for compounds"
    )
    parser.add_argument(
        "input_file",
        type=Path,
        help="Path to TSV file containing compounds",
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        help="Path to save predictions TSV",
        default="predictions.tsv",
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        help="Directory containing trained models",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help="Directory for caching",
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
        # Load compounds
        logging.info(f"Loading compounds from {args.input_file}")
        compounds = load_compounds(args.input_file)
        logging.info(f"Loaded {len(compounds)} compounds")

        # Make predictions
        predictor = predict_bbb_permeability(
            compounds,
            model_dir=args.model_dir,
            cache_dir=args.cache_dir,
        )

        # Export predictions
        logging.info(f"Exporting predictions to {args.output_file}")
        predictor.export_predictions(
            args.output_file,
            include_supporting_data=True,
            include_web_data=True,
        )
        logging.info("Done!")

    except Exception as e:
        logging.error(f"Error: {str(e)}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
