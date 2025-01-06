#!/usr/bin/env python
"""Quick start example for Blood-Brain Barrier (BBB) permeability prediction.

This script demonstrates how to:
1. Load and validate compound data from TSV/JSON files
2. Initialize the BBB predictor with web enrichment
3. Make predictions with confidence scores
4. Export detailed results in multiple formats
5. Generate a summary report

The script handles both CNS-active (BBB permeable) and peripherally-selective
(low BBB permeability) compounds, providing predictions and supporting data
from multiple sources.
"""

import json
import logging
import sys
from pathlib import Path
from typing import List, Dict, Union, Optional

from tqdm import tqdm

from binding_data_processor.models.core import CompoundData
from binding_data_processor.processors.psychopharm.predictors.bbb.base import BBBPredictorWebEnriched
from binding_data_processor.processors.structure.properties.descriptors import calculate_molecular_properties


class CompoundLoadError(Exception):
    """Raised when there is an error loading compound data."""

    pass


class PredictionError(Exception):
    """Raised when there is an error making predictions."""

    pass


def setup_logging(log_file: Optional[Path] = None):
    """Set up logging configuration with optional file output.

    Args:
        log_file: Optional path to write logs to file in addition to console
    """
    handlers = [logging.StreamHandler(sys.stdout)]

    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", handlers=handlers
    )


def validate_compound_data(name: str, smiles: str, cas: str) -> None:
    """Validate compound data fields.

    Args:
        name: Compound name
        smiles: SMILES structure
        cas: CAS registry number

    Raises:
        CompoundLoadError: If validation fails
    """
    if not name or not isinstance(name, str):
        raise CompoundLoadError(f"Invalid compound name: {name}")

    if not smiles or not isinstance(smiles, str):
        raise CompoundLoadError(f"Invalid SMILES for compound {name}: {smiles}")

    if not cas or not isinstance(cas, str):
        raise CompoundLoadError(f"Invalid CAS number for compound {name}: {cas}")


def load_compounds(input_file: Path) -> List[CompoundData]:
    """Load and validate compounds from TSV or JSON file.

    Args:
        input_file: Path to input file (.tsv or .json)

    Returns:
        List of validated CompoundData objects

    Raises:
        CompoundLoadError: If file cannot be loaded or data is invalid
    """
    compounds = []

    try:
        if input_file.suffix == ".json":
            with open(input_file) as f:
                data = json.load(f)
                for item in data:
                    validate_compound_data(item["name"], item["smiles"], item["cas_number"])
                    compounds.append(
                        CompoundData(name=item["name"], smiles=item["smiles"], cas_number=item["cas_number"])
                    )

        else:  # TSV
            with open(input_file) as f:
                next(f)  # Skip header
                for line in f:
                    if line.strip() and not line.startswith("#"):
                        try:
                            name, smiles, cas = line.strip().split("\t")
                            validate_compound_data(name, smiles, cas)
                            compounds.append(CompoundData(name=name, smiles=smiles, cas_number=cas))
                        except ValueError as e:
                            raise CompoundLoadError(f"Invalid line format: {line.strip()}") from e

    except (json.JSONDecodeError, OSError) as e:
        raise CompoundLoadError(f"Failed to load compounds from {input_file}: {str(e)}") from e

    if not compounds:
        raise CompoundLoadError(f"No valid compounds found in {input_file}")

    return compounds


def export_results(results: List[Dict], output_file: Path, format: str = "tsv") -> None:
    """Export prediction results in specified format.

    Args:
        results: List of prediction result dictionaries
        output_file: Output file path
        format: Output format ('tsv' or 'json')
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if format == "json":
        with open(output_file, "w") as f:
            json.dump(results, f, indent=2)
    else:
        with open(output_file, "w") as f:
            # Write header
            headers = [
                "compound_name",
                "smiles",
                "cas_number",
                "bbb_class",
                "confidence",
                "molecular_weight",
                "logp",
                "hbd",
                "hba",
                "tpsa",
                "rotatable_bonds",
                "supporting_data",
                "web_data",
            ]
            f.write("\t".join(headers) + "\n")

            # Write data
            for result in results:
                row = [
                    result["compound_name"],
                    result["smiles"],
                    result["cas_number"],
                    result["bbb_class"],
                    f"{result['confidence']:.2f}",
                    f"{result['properties']['molecular_weight']:.1f}",
                    f"{result['properties']['logp']:.1f}",
                    str(result["properties"]["hbd"]),
                    str(result["properties"]["hba"]),
                    f"{result['properties']['tpsa']:.1f}",
                    str(result["properties"]["rotatable_bonds"]),
                    json.dumps(result["supporting_data"]),
                    json.dumps(result["web_data"]),
                ]
                f.write("\t".join(row) + "\n")


def main():
    """Run BBB prediction example with enhanced features."""
    # Set up logging with file output
    root_dir = Path(__file__).parent.parent.parent
    log_file = root_dir / "output/bbb_prediction.log"
    setup_logging(log_file)

    try:
        # Load example compounds
        example_file = root_dir / "examples/data/example_compounds.tsv"
        compounds = load_compounds(example_file)
        logging.info(f"Successfully loaded {len(compounds)} compounds")

        # Initialize predictor
        predictor = BBBPredictorWebEnriched(
            model_dir=str(root_dir / "models/bbb"),
            cache_dir=str(root_dir / "cache"),
        )
        logging.info("Initialized BBB predictor with web enrichment")

        # Make predictions with progress bar
        results = []
        for compound in tqdm(compounds, desc="Making predictions"):
            try:
                # Calculate molecular properties
                properties = calculate_molecular_properties(compound.smiles)

                # Get BBB prediction
                result = predictor.predict(compound)

                # Compile result data
                result_data = {
                    "compound_name": compound.name,
                    "smiles": compound.smiles,
                    "cas_number": compound.cas_number,
                    "bbb_class": result.value,
                    "confidence": result.confidence,
                    "properties": properties,
                    "supporting_data": result.supporting_data,
                    "web_data": result.web_data if hasattr(result, "web_data") else {},
                }
                results.append(result_data)

                # Log detailed prediction
                logging.info(
                    f"\nPrediction for {compound.name}:"
                    f"\n  BBB Class: {result.value}"
                    f"\n  Confidence: {result.confidence:.2f}"
                    f"\n  Molecular Weight: {properties['molecular_weight']:.1f}"
                    f"\n  LogP: {properties['logp']:.1f}"
                    f"\n  TPSA: {properties['tpsa']:.1f}"
                )

            except Exception as e:
                logging.error(f"Error predicting {compound.name}: {str(e)}")
                continue

        # Export results in both formats
        tsv_file = root_dir / "output/bbb_predictions.tsv"
        json_file = root_dir / "output/bbb_predictions.json"

        export_results(results, tsv_file, format="tsv")
        export_results(results, json_file, format="json")

        logging.info(f"\nExported predictions to:")
        logging.info(f"- TSV: {tsv_file}")
        logging.info(f"- JSON: {json_file}")

        # Generate summary statistics
        bbb_positive = sum(1 for r in results if r["bbb_class"] == "BBB_POSITIVE")
        high_confidence = sum(1 for r in results if r["confidence"] > 0.8)

        print("\nPrediction Summary")
        print("==================")
        print(f"Total compounds analyzed: {len(results)}")
        print(f"BBB positive compounds: {bbb_positive}")
        print(f"BBB negative compounds: {len(results) - bbb_positive}")
        print(f"High confidence predictions (>0.8): {high_confidence}")
        print(f"\nResults saved to:")
        print(f"- TSV: {tsv_file}")
        print(f"- JSON: {json_file}")
        print(f"- Log: {log_file}")

        print("\nOutput includes:")
        print("- Basic compound information (name, SMILES, CAS)")
        print("- BBB prediction class and confidence")
        print("- Molecular properties (MW, LogP, TPSA, etc.)")
        print("- Supporting data from prediction model")
        print("- Enrichment data from web sources")

    except Exception as e:
        logging.error(f"Error running BBB prediction: {str(e)}")
        raise

    logging.info("BBB prediction completed successfully")


if __name__ == "__main__":
    main()
