"""Pharmacophore-based ML models and utilities.

This module provides comprehensive pharmacophore analysis capabilities:

Core Functionality:
1. Pharmacophore feature detection and extraction
   - Identification of key pharmacophoric points
   - Feature type classification and annotation
   - Spatial relationship analysis
   - Feature vector generation

2. Conformer generation and analysis
   - 3D conformer generation with RDKit
   - Conformer energy minimization
   - Ensemble generation and filtering
   - Conformer clustering and selection

3. Pharmacophore-based similarity calculations
   - Feature-based fingerprint generation
   - 3D pharmacophore alignment
   - Shape-based similarity metrics
   - Pharmacophore matching scores

4. Machine learning model integration
   - Feature extraction and preprocessing
   - Model training on pharmacophore features
   - Prediction using pharmacophore patterns
   - Ensemble methods combining multiple features

Key Components:
- ConformerGenerator: Generates and optimizes 3D conformers
- PharmacophoreGenerator: Extracts and analyzes pharmacophore features
- Feature extraction pipeline for ML model input
- Similarity calculation methods for compound comparison
- Integration with structure-based ML models

The module is designed to work seamlessly with other components:
- Structure processing and analysis
- Machine learning model training
- Binding site prediction
- Activity prediction
- Property calculation
"""

from .conformers import ConformerGenerator
from .generator import PharmacophoreGenerator

__all__ = ["ConformerGenerator", "PharmacophoreGenerator"]
