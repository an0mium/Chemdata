"""Chemical data processing pipeline.

This package provides a comprehensive pipeline for:
1. Loading and processing BindingDB data
2. Enriching compounds with web data
3. Running ML predictions
4. Analyzing compounds
5. Exporting results

Key Features:
- Automated BindingDB data processing
- Web data enrichment and social media harvesting
- ML predictions and analysis
- Data validation and deduplication
- Progress tracking and error handling
- Checkpointing and caching
- Parallel processing and error recovery
- Model ensembling and confidence scoring

The pipeline supports both batch processing and individual compound processing,
with extensive configuration options and progress tracking.

Example usage:

```python
from binding_data_processor.pipeline import Pipeline

# Initialize pipeline
pipeline = Pipeline(
    data_dir="data",
    cache_dir="cache",
    checkpoint_dir="checkpoints",
    config={
        "n_workers": 4,
        "batch_size": 100,
        "checkpoint_interval": 1000,
        "max_retries": 3,
    }
)

# Process BindingDB data
pipeline.process_bindingdb(
    input_file="bindingdb_data.tsv",
    output_file="processed_compounds.tsv"
)

# Process individual compound
compound = pipeline.process_compound(
    compound_data,
    enrich=True,
    predict=True,
    analyze=True
)
```
"""

from .base import Pipeline, PipelineConfig, PipelineStats
from .ml import MLManager, ModelConfig, PredictionStats
from .web import WebManager, WebConfig, EnrichmentStats
from .validation import ValidationManager, ValidationConfig, ValidationStats
from .analysis import AnalysisManager, AnalysisConfig, AnalysisStats

__all__ = [
    # Main interface
    "Pipeline",
    "PipelineConfig",
    "PipelineStats",
    
    # ML components
    "MLManager",
    "ModelConfig",
    "PredictionStats",
    
    # Web components
    "WebManager",
    "WebConfig",
    "EnrichmentStats",
    
    # Validation components
    "ValidationManager",
    "ValidationConfig",
    "ValidationStats",
    
    # Analysis components
    "AnalysisManager",
    "AnalysisConfig",
    "AnalysisStats",
]

# Version info
__version__ = "2.1.0"
