Analysis Pipelines
=================

This guide shows how to implement comprehensive analysis pipelines for compound data.

Binding Analysis
-------------

Analyzing receptor binding:

.. code-block:: python

    from binding_data_processor.pipeline.analysis import (
        BindingAnalyzer,
        AnalysisConfig,
    )

    # Configure analyzer
    analyzer = BindingAnalyzer(
        config=AnalysisConfig(
            receptors=["5-HT2A", "NMDA", "D2"],
            min_affinity=1e-6,
            include_confidence=True,
        )
    )

    # Analyze compounds
    for compound in compounds:
        # Analyze binding
        results = analyzer.analyze_binding(compound)
        
        # Extract profiles
        binding_profile = results.binding_profile
        selectivity_profile = results.selectivity_profile
        
        # Get strongest targets
        strongest_targets = results.get_strongest_targets(top_n=3)
        
        # Update compound
        compound.binding_analysis = {
            "binding_profile": binding_profile,
            "selectivity_profile": selectivity_profile,
            "strongest_targets": strongest_targets,
        }

Activity Analysis
--------------

Analyzing compound activity:

.. code-block:: python

    from binding_data_processor.pipeline.analysis import ActivityAnalyzer

    # Configure analyzer
    analyzer = ActivityAnalyzer(
        config=AnalysisConfig(
            effect_types=[
                "psychedelic",
                "stimulant",
                "depressant",
                "nootropic",
            ],
            min_confidence=0.7,
        )
    )

    # Analyze compounds
    for compound in compounds:
        # Analyze activity
        results = analyzer.analyze_activity(compound)
        
        # Extract profiles
        effect_profile = results.effect_profile
        mechanism_profile = results.mechanism_profile
        
        # Get primary effects
        primary_effects = results.get_primary_effects(top_n=3)
        
        # Update compound
        compound.activity_analysis = {
            "effect_profile": effect_profile,
            "mechanism_profile": mechanism_profile,
            "primary_effects": primary_effects,
        }

Safety Analysis
------------

Analyzing safety profiles:

.. code-block:: python

    from binding_data_processor.pipeline.analysis import SafetyAnalyzer

    # Configure analyzer
    analyzer = SafetyAnalyzer(
        config=AnalysisConfig(
            risk_types=[
                "toxicity",
                "addiction",
                "interaction",
            ],
            data_sources=[
                "predicted",
                "literature",
                "community",
            ],
        )
    )

    # Analyze compounds
    for compound in compounds:
        # Analyze safety
        results = analyzer.analyze_safety(compound)
        
        # Extract profiles
        risk_profile = results.risk_profile
        interaction_profile = results.interaction_profile
        
        # Get major risks
        major_risks = results.get_major_risks(min_level="HIGH")
        
        # Update compound
        compound.safety_analysis = {
            "risk_profile": risk_profile,
            "interaction_profile": interaction_profile,
            "major_risks": major_risks,
        }

Property Analysis
--------------

Analyzing molecular properties:

.. code-block:: python

    from binding_data_processor.pipeline.analysis import PropertyAnalyzer

    # Configure analyzer
    analyzer = PropertyAnalyzer(
        config=AnalysisConfig(
            properties=[
                "molecular_weight",
                "logp",
                "tpsa",
                "hbd",
                "hba",
            ],
            include_descriptors=True,
        )
    )

    # Analyze compounds
    for compound in compounds:
        # Analyze properties
        results = analyzer.analyze_properties(compound)
        
        # Extract profiles
        property_profile = results.property_profile
        descriptor_profile = results.descriptor_profile
        
        # Get druglikeness
        druglikeness = results.calculate_druglikeness()
        
        # Update compound
        compound.property_analysis = {
            "property_profile": property_profile,
            "descriptor_profile": descriptor_profile,
            "druglikeness": druglikeness,
        }

SAR Analysis
----------

Structure-activity relationship analysis:

.. code-block:: python

    from binding_data_processor.pipeline.analysis import SARAnalyzer

    # Configure analyzer
    analyzer = SARAnalyzer(
        config=AnalysisConfig(
            similarity_threshold=0.7,
            activity_threshold=1.0,
            max_compounds=1000,
        )
    )

    # Analyze compound series
    results = analyzer.analyze_series(compounds)

    # Extract series
    for series in results.series:
        # Get core structure
        core = series.core_structure
        
        # Get activity cliffs
        cliffs = series.activity_cliffs
        
        # Get SAR patterns
        patterns = series.sar_patterns
        
        # Update compounds
        for compound in series.compounds:
            compound.sar_analysis = {
                "series": series.name,
                "core": core,
                "cliffs": cliffs,
                "patterns": patterns,
            }

Analysis Pipeline
--------------

Creating a comprehensive analysis pipeline:

.. code-block:: python

    from binding_data_processor.pipeline import AnalysisPipeline

    # Configure pipeline
    config = AnalysisConfig(
        analyzers=[
            BindingAnalyzer(),
            ActivityAnalyzer(),
            SafetyAnalyzer(),
            PropertyAnalyzer(),
            SARAnalyzer(),
        ],
        cache_enabled=True,
        batch_size=100,
        num_workers=4,
    )

    # Create pipeline
    pipeline = AnalysisPipeline(config)

    # Process compounds
    analyzed_compounds = pipeline.analyze_compounds(
        input_file="compounds.tsv",
        output_file="analyzed_compounds.tsv",
    )

    # Get statistics
    stats = pipeline.get_stats()
    print(f"Analysis Stats: {stats}")

Custom Analysis
------------

Creating custom analyzers:

.. code-block:: python

    from binding_data_processor.pipeline.analysis import BaseAnalyzer

    class CustomAnalyzer(BaseAnalyzer):
        """Custom compound analyzer."""

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
                self.logger.error(f"Analysis error: {str(e)}")
                raise

        def _analyze_compound(self, compound: CompoundData) -> Dict[str, Any]:
            """Perform compound analysis."""
            # Custom analysis logic
            return {
                "metric_a": self._calculate_metric_a(compound),
                "metric_b": self._calculate_metric_b(compound),
                "score": self._calculate_score(compound),
            }

        def get_stats(self) -> Dict[str, Any]:
            """Get analyzer statistics."""
            return {
                "analysis_count": self.analysis_count,
                "config": self.config.dict(),
            }

Analysis Visualization
------------------

Visualizing analysis results:

.. code-block:: python

    from binding_data_processor.pipeline.analysis import AnalysisVisualizer

    # Configure visualizer
    visualizer = AnalysisVisualizer(
        config=AnalysisConfig(
            plot_types=[
                "binding_heatmap",
                "activity_radar",
                "safety_matrix",
                "property_scatter",
            ],
            interactive=True,
        )
    )

    # Create visualizations
    for compound in compounds:
        # Generate plots
        plots = visualizer.create_plots(compound)
        
        # Save plots
        for name, plot in plots.items():
            plot.write_html(f"{compound.id}_{name}.html")

Report Generation
--------------

Generating analysis reports:

.. code-block:: python

    from binding_data_processor.pipeline.analysis import ReportGenerator

    # Configure generator
    generator = ReportGenerator(
        config=AnalysisConfig(
            sections=[
                "binding",
                "activity",
                "safety",
                "properties",
                "sar",
            ],
            include_plots=True,
            template="detailed",
        )
    )

    # Generate reports
    for compound in compounds:
        # Generate report
        report = generator.generate_report(compound)
        
        # Save report
        report.save(f"{compound.id}_report.html")
