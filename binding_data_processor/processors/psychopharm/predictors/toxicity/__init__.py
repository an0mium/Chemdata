"""Toxicity prediction package.

This package provides functionality for predicting toxicological properties
of chemical compounds, including:

1. Core functionality:
   - Toxicity class prediction
   - Mechanism of action prediction
   - Organ-specific toxicity assessment
   - Safety concern evaluation

2. Feature types:
   - Molecular fingerprints
   - Chemical descriptors
   - Enhanced features (interactions & polynomial terms)

3. Model management:
   - Model loading and versioning
   - Feature scaling
   - Prediction history tracking

4. Prediction types:
   - Overall toxicity classification
   - Cellular mechanisms
   - Molecular mechanisms
   - Systemic effects
   - Organ-specific toxicity
   - Safety concerns
"""

from .base import ToxicityPredictorBase as ToxicityPredictor
from .features import (
    extract_fingerprints,
    extract_descriptors,
    extract_features,
    extract_all_features,
    extract_enhanced_features,
)
from .model_loading import (
    load_models,
    initialize_models,
    initialize_scalers,
)
from .prediction import (
    predict_mechanism,
    predict_organ_effects,
    predict_safety_concerns,
    predict_all_effects,
)

__all__ = [
    # Base predictor
    "ToxicityPredictor",
    # Feature extraction
    "extract_fingerprints",
    "extract_descriptors",
    "extract_features",
    "extract_all_features",
    "extract_enhanced_features",
    # Model management
    "load_models",
    "initialize_models",
    "initialize_scalers",
    # Prediction functions
    "predict_mechanism",
    "predict_organ_effects",
    "predict_safety_concerns",
    "predict_all_effects",
]
