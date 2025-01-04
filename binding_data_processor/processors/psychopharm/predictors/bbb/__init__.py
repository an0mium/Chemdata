"""BBB permeability prediction package.

This package provides classes for predicting blood-brain barrier permeability:

1. BBBPredictorBase - Core BBB prediction functionality
2. BBBPredictor - BBB prediction with transporter analysis
3. BBBPredictorEnhanced - BBB prediction with integrated ML models
4. BBBPredictorWebEnriched - BBB prediction with web data enrichment

The modules are organized in a hierarchy of increasing functionality:
base.py -> predictors.py -> integration.py -> enrichment.py

Example usage:
    ```python
    from binding_data_processor.processors.psychopharm.predictors.bbb import (
        BBBPredictorWebEnriched
    )

    # Initialize predictor with web enrichment
    predictor = BBBPredictorWebEnriched(
        model_dir="models/bbb",
        cache_dir="cache",
    )

    # Make predictions for a compound
    result = predictor.predict(compound)
    print(f"BBB Class: {result.value}")
    print(f"Confidence: {result.confidence:.2f}")
    print("Supporting Data:")
    for key, value in result.supporting_data.items():
        print(f"  {key}: {value}")

    # Export predictions to TSV
    predictor.export_predictions(
        "predictions.tsv",
        include_supporting_data=True,
        include_web_data=True,
    )
    ```
"""

from .base import BBBPredictorBase
from .predictors import BBBPredictor
from .integration import BBBPredictorEnhanced
from .enrichment import BBBPredictorWebEnriched

__all__ = [
    "BBBPredictorBase",
    "BBBPredictor",
    "BBBPredictorEnhanced",
    "BBBPredictorWebEnriched",
]
