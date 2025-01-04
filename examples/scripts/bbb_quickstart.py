#!/usr/bin/env python
"""Quick start example for BBB permeability prediction.

This script demonstrates how to:
1. Load example compounds
2. Initialize the BBB predictor
3. Make predictions
4. Export results
"""

import logging
from pathlib import Path

from binding_data_processor.models.core import CompoundData
from binding_data_processor.processors.psychopharm.predictors.bbb import (
    BBBPredictorWebEnriched
)


def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )


def load_example_compounds(example_file: Path) -> list[CompoundData]:
    """Load example compounds from TSV file."""
    compounds = []
    
    with open(example_file) as f:
        next(f)  # Skip header
        for line in f:
            if line.strip() and not line.startswith("#"):
                name, smiles, cas = line.strip().split("\t")
                compound = CompoundData(
                    name=name,
                    smiles=smiles,
                    cas_number=cas,
                )
                compounds.append(compound)
    
    return compounds


def main():
    """Run BBB prediction example."""
    # Set up logging
    setup_logging()
    
    try:
        # Get project root directory
        root_dir = Path(__file__).parent.parent.parent
        
        # Load example compounds
        example_file = root_dir / "examples/data/example_compounds.tsv"
        compounds = load_example_compounds(example_file)
        logging.info(f"Loaded {len(compounds)} example compounds")
        
        # Initialize predictor
        predictor = BBBPredictorWebEnriched(
            model_dir=str(root_dir / "models/bbb"),
            cache_dir=str(root_dir / "cache"),
        )
        logging.info("Initialized BBB predictor")
        
        # Make predictions
        for compound in compounds:
            result = predictor.predict(compound)
            
            # Log prediction
            logging.info(
                f"\nCompound: {compound.name}"
                f"\n  BBB Class: {result.value}"
                f"\n  Confidence: {result.confidence:.2f}"
            )
            
            # Log supporting data
            logging.debug("\nSupporting Data:")
            for key, value in result.supporting_data.items():
                logging.debug(f"  {key}: {value}")
        
        # Export predictions
        output_file = root_dir / "output/bbb_predictions.tsv"
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        predictor.export_predictions(
            output_file=output_file,
            include_supporting_data=True,
            include_web_data=True,
        )
        logging.info(f"\nExported predictions to {output_file}")
        
        # Print summary
        print("\nPrediction Summary:")
        print("==================")
        print(f"Total compounds: {len(compounds)}")
        print(f"Results saved to: {output_file}")
        print("\nOutput columns:")
        print("- compound_name")
        print("- smiles")
        print("- cas_number")
        print("- bbb_class")
        print("- confidence")
        print("- supporting_data")
        print("- web_data")
        
    except Exception as e:
        logging.error(f"Error running example: {str(e)}")
        raise


if __name__ == "__main__":
    main()
