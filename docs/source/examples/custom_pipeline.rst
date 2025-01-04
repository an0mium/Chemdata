Custom Pipeline Components
======================

This guide shows how to create custom pipeline components to extend ChemData's functionality.

Custom Stage
----------

Creating a custom pipeline stage:

.. code-block:: python

    from binding_data_processor.pipeline import Stage, StageConfig
    from binding_data_processor.models import CompoundData
    from typing import List, Dict, Any

    class CustomStageConfig(StageConfig):
        """Configuration for custom stage."""
        def __init__(
            self,
            threshold: float = 0.5,
            include_metadata: bool = True,
            **kwargs: Any,
        ):
            super().__init__(**kwargs)
            self.threshold = threshold
            self.include_metadata = include_metadata

    class CustomStage(Stage):
        """Custom pipeline stage for specialized processing."""

        def __init__(self, config: CustomStageConfig):
            super().__init__()
            self.config = config
            self.processed_count = 0

        def process(self, compounds: List[CompoundData]) -> List[CompoundData]:
            """Process compounds with custom logic.
            
            Args:
                compounds: List of compounds to process
                
            Returns:
                Processed compounds
            """
            results = []
            for compound in compounds:
                try:
                    # Custom processing logic
                    processed = self._process_compound(compound)
                    if processed is not None:
                        results.append(processed)
                        self.processed_count += 1
                except Exception as e:
                    self.logger.error(
                        "processing_failed",
                        compound_id=compound.id,
                        error=str(e),
                    )
            return results

        def _process_compound(self, compound: CompoundData) -> Optional[CompoundData]:
            """Process individual compound."""
            # Add custom data
            compound.custom_data = {
                "score": self._calculate_score(compound),
                "metadata": self._get_metadata(compound),
            }
            
            # Filter based on threshold
            if compound.custom_data["score"] < self.config.threshold:
                return None
                
            return compound

        def _calculate_score(self, compound: CompoundData) -> float:
            """Calculate custom score."""
            # Custom scoring logic
            return 0.8

        def _get_metadata(self, compound: CompoundData) -> Dict[str, Any]:
            """Get custom metadata."""
            if not self.config.include_metadata:
                return {}
            
            # Custom metadata logic
            return {
                "processed_at": datetime.now().isoformat(),
                "version": "1.0",
            }

        def get_stats(self) -> Dict[str, Any]:
            """Get stage statistics."""
            return {
                "processed_count": self.processed_count,
                "threshold": self.config.threshold,
            }

Using Custom Stage
---------------

Using the custom stage in a pipeline:

.. code-block:: python

    from binding_data_processor.pipeline import ProcessingPipeline

    # Configure custom stage
    custom_config = CustomStageConfig(
        threshold=0.7,
        include_metadata=True,
    )

    # Create pipeline with custom stage
    pipeline = ProcessingPipeline([
        ValidationStage(),
        EnrichmentStage(),
        CustomStage(custom_config),
        ExportStage(),
    ])

    # Process compounds
    results = pipeline.process_compounds(input_file)

Custom Predictor
-------------

Creating a custom ML predictor:

.. code-block:: python

    from binding_data_processor.pipeline.ml import BasePredictor, PredictorConfig
    import torch
    import torch.nn as nn

    class CustomModel(nn.Module):
        """Custom neural network model."""
        
        def __init__(self, input_size: int, hidden_size: int):
            super().__init__()
            self.layers = nn.Sequential(
                nn.Linear(input_size, hidden_size),
                nn.ReLU(),
                nn.Linear(hidden_size, 1),
                nn.Sigmoid(),
            )

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            return self.layers(x)

    class CustomPredictor(BasePredictor):
        """Custom predictor for specialized predictions."""

        def __init__(
            self,
            config: PredictorConfig,
            model_path: Optional[str] = None,
        ):
            super().__init__(config)
            self.model = self._load_model(model_path)
            self.feature_extractor = self._create_feature_extractor()

        def predict(self, compound: CompoundData) -> PredictionResult:
            """Make prediction for compound."""
            # Extract features
            features = self.feature_extractor.extract_features(compound)
            
            # Convert to tensor
            x = torch.tensor(features, dtype=torch.float32)
            
            # Make prediction
            with torch.no_grad():
                output = self.model(x)
            
            # Process result
            score = float(output.item())
            confidence = self._calculate_confidence(output)
            
            return PredictionResult(
                value=score,
                confidence=confidence,
                metadata=self._get_metadata(compound),
            )

        def _load_model(self, model_path: Optional[str]) -> nn.Module:
            """Load or create model."""
            if model_path:
                return torch.load(model_path)
            
            return CustomModel(
                input_size=self.config.input_size,
                hidden_size=self.config.hidden_size,
            )

        def _create_feature_extractor(self) -> FeatureExtractor:
            """Create feature extractor."""
            return FeatureExtractor(
                fingerprint_size=self.config.input_size,
                use_descriptors=True,
            )

        def _calculate_confidence(self, output: torch.Tensor) -> float:
            """Calculate prediction confidence."""
            # Custom confidence calculation
            return min(1.0, abs(0.5 - float(output.item())) * 2)

Using Custom Predictor
------------------

Using the custom predictor:

.. code-block:: python

    # Configure predictor
    config = PredictorConfig(
        input_size=1024,
        hidden_size=512,
        batch_size=32,
    )

    # Create predictor
    predictor = CustomPredictor(
        config=config,
        model_path="models/custom/model.pt",
    )

    # Make predictions
    for compound in compounds:
        result = predictor.predict(compound)
        print(f"Score: {result.value:.2f}")
        print(f"Confidence: {result.confidence:.2f}")

Custom Analyzer
------------

Creating a custom analyzer:

.. code-block:: python

    from binding_data_processor.pipeline.analysis import BaseAnalyzer, AnalysisConfig
    from typing import Dict, Any

    class CustomAnalyzer(BaseAnalyzer):
        """Custom analyzer for specialized analysis."""

        def __init__(self, config: AnalysisConfig):
            super().__init__(config)
            self.analysis_count = 0

        def analyze(self, compound: CompoundData) -> AnalysisResult:
            """Analyze compound."""
            try:
                # Perform analysis
                data = self._analyze_compound(compound)
                
                # Create result
                result = AnalysisResult(
                    compound_id=compound.id,
                    data=data,
                    metadata=self._get_metadata(compound),
                )
                
                self.analysis_count += 1
                return result
                
            except Exception as e:
                self.logger.error(
                    "analysis_failed",
                    compound_id=compound.id,
                    error=str(e),
                )
                raise

        def _analyze_compound(self, compound: CompoundData) -> Dict[str, Any]:
            """Perform compound analysis."""
            # Custom analysis logic
            return {
                "property_a": self._calculate_property_a(compound),
                "property_b": self._calculate_property_b(compound),
                "score": self._calculate_score(compound),
            }

        def _calculate_property_a(self, compound: CompoundData) -> float:
            """Calculate property A."""
            # Custom calculation
            return 0.5

        def _calculate_property_b(self, compound: CompoundData) -> float:
            """Calculate property B."""
            # Custom calculation
            return 0.8

        def _calculate_score(self, compound: CompoundData) -> float:
            """Calculate analysis score."""
            # Custom scoring
            return 0.7

        def get_stats(self) -> Dict[str, Any]:
            """Get analyzer statistics."""
            return {
                "analysis_count": self.analysis_count,
                "config": self.config.dict(),
            }

Using Custom Analyzer
-----------------

Using the custom analyzer:

.. code-block:: python

    # Configure analyzer
    config = AnalysisConfig(
        include_metadata=True,
        cache_results=True,
    )

    # Create analyzer
    analyzer = CustomAnalyzer(config=config)

    # Analyze compounds
    for compound in compounds:
        result = analyzer.analyze(compound)
        print(f"Compound: {result.compound_id}")
        print(f"Data: {result.data}")
        print(f"Metadata: {result.metadata}")

Integration
---------

Integrating custom components:

.. code-block:: python

    # Create pipeline with custom components
    pipeline = ProcessingPipeline([
        ValidationStage(),
        CustomStage(custom_config),
        PredictionStage([
            CustomPredictor(predictor_config),
        ]),
        AnalysisStage([
            CustomAnalyzer(analyzer_config),
        ]),
        ExportStage(),
    ])

    # Process compounds
    results = pipeline.process_compounds(
        input_file="compounds.tsv",
        output_file="results.tsv",
    )

    # Get statistics
    stats = pipeline.get_stats()
    print(f"Processing Stats: {stats}")
