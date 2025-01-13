"""Chemical structure processing and analysis.

This module provides comprehensive functionality for:

1. Structure Processing and Validation
- SMILES/InChI parsing and conversion 
- Structure cleaning and normalization
- Tautomer/resonance handling
- Salt/solvent removal
- 3D conformation generation
- Force field optimization

2. Molecular Analysis
- Physicochemical descriptors
- Topological indices
- Electronic properties
- Fragment features
- 3D shape descriptors
- Pharmacophore detection

3. Chemical Similarity
- Multiple fingerprint types
- Various similarity metrics
- Substructure matching
- Maximum common substructure
- Shape similarity

4. Machine Learning
- Activity prediction
- Toxicity prediction
- Property prediction
- Structure-activity relationships
- Deep learning on molecular graphs

5. Structure Visualization
- 2D depiction
- 3D conformation viewing
- Property visualization
- Pharmacophore display
- Similarity highlighting

6. Web Integration
- Structure searching
- Property calculation
- Activity prediction
- Data harvesting
- Result visualization
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Type, Union

from .base import BaseStructureProcessor
from .advanced import AdvancedStructureProcessor
from .descriptors import DescriptorCalculator
from .pharmacophore import PharmacophoreDetector
from .similarity import SimilarityProcessor
from .depiction import StructureDepiction
from .ml import (
    ActivityPredictor,
    ToxicityPredictor,
    PropertyPredictor,
    GraphNeuralNetwork,
    DeepChemModel,
)

# Configure logging
logger = logging.getLogger(__name__)

# Core structure processors
__all__ = [
    "BaseStructureProcessor",
    "AdvancedStructureProcessor",
]

# Specialized processors
__all__ += [
    "DescriptorCalculator",
    "PharmacophoreDetector",
    "SimilarityProcessor",
    "StructureDepiction",
]

# ML models
__all__ += [
    "ActivityPredictor",
    "ToxicityPredictor",
    "PropertyPredictor",
    "GraphNeuralNetwork",
    "DeepChemModel",
]

# Version info
__version__ = "1.0.0"

# Default configuration
DEFAULT_CONFIG = {
    # Fingerprint settings
    "fingerprint_bits": 2048,
    "fingerprint_radius": 2,
    "use_features": True,
    # Similarity thresholds
    "similarity_threshold": 0.7,
    "activity_cliff_threshold": 10.0,
    # 3D conformation settings
    "max_conformers": 50,
    "energy_threshold": 10.0,
    "rmsd_threshold": 2.0,
    # ML model settings
    "model_type": "graph",  # graph, deepchem, ensemble
    "batch_size": 32,
    "learning_rate": 0.001,
    "num_epochs": 100,
    # Web scraping settings
    "max_retries": 3,
    "timeout": 30,
    "user_agent": "ChemDataProcessor/1.0",
    # Visualization settings
    "image_size": (400, 400),
    "highlight_color": (0.7, 0.7, 1.0),
    "background_color": (1.0, 1.0, 1.0),
}


def get_processor(processor_type: str = "base", config: Optional[Dict] = None, **kwargs) -> Union[
    BaseStructureProcessor,
    AdvancedStructureProcessor,
    DescriptorCalculator,
    PharmacophoreDetector,
    SimilarityProcessor,
    StructureDepiction,
    ActivityPredictor,
    ToxicityPredictor,
    PropertyPredictor,
    GraphNeuralNetwork,
    DeepChemModel,
]:
    """
    Factory function to get structure processor instance.

    Args:
        processor_type: Type of processor to create
            - base: Basic structure processing
            - advanced: Advanced structure processing
            - descriptors: Molecular descriptor calculation
            - pharmacophore: Pharmacophore detection
            - similarity: Chemical similarity calculation
            - depiction: Structure visualization
            - activity: Activity prediction
            - toxicity: Toxicity prediction
            - property: Property prediction
            - graph: Graph neural networks
            - deepchem: DeepChem models
        config: Configuration dictionary
        **kwargs: Additional arguments passed to processor

    Returns:
        Structure processor instance

    Raises:
        ValueError: If processor_type is unknown
    """
    # Merge config with defaults
    full_config = DEFAULT_CONFIG.copy()
    if config:
        full_config.update(config)

    # Add config to kwargs
    kwargs["config"] = full_config

    # Map processor types to classes
    processors = {
        # Core processors
        "base": BaseStructureProcessor,
        "advanced": AdvancedStructureProcessor,
        # Specialized processors
        "descriptors": DescriptorCalculator,
        "pharmacophore": PharmacophoreDetector,
        "similarity": SimilarityProcessor,
        "depiction": StructureDepiction,
        # ML models
        "activity": ActivityPredictor,
        "toxicity": ToxicityPredictor,
        "property": PropertyPredictor,
        "graph": GraphNeuralNetwork,
        "deepchem": DeepChemModel,
    }

    if processor_type not in processors:
        raise ValueError(f"Unknown processor type: {processor_type}. " f"Available types: {list(processors.keys())}")

    try:
        return processors[processor_type](**kwargs)
    except Exception as e:
        logger.error(f"Error creating {processor_type} processor: {str(e)}")
        raise


def get_all_processors(config: Optional[Dict] = None, **kwargs) -> Dict[
    str,
    Union[
        BaseStructureProcessor,
        AdvancedStructureProcessor,
        DescriptorCalculator,
        PharmacophoreDetector,
        SimilarityProcessor,
        StructureDepiction,
        ActivityPredictor,
        ToxicityPredictor,
        PropertyPredictor,
        GraphNeuralNetwork,
        DeepChemModel,
    ],
]:
    """
    Get instances of all available processors.

    Args:
        config: Configuration dictionary
        **kwargs: Additional arguments passed to processors

    Returns:
        Dictionary mapping processor names to instances
    """
    return {
        name: get_processor(name, config, **kwargs)
        for name in [
            "base",
            "advanced",
            "descriptors",
            "pharmacophore",
            "similarity",
            "depiction",
            "activity",
            "toxicity",
            "property",
            "graph",
            "deepchem",
        ]
    }
