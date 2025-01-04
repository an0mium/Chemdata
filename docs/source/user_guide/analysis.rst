Analysis Guide
=============

This guide covers the analysis capabilities of ChemData in detail.

Analysis Overview
--------------

The analysis pipeline includes:

1. Binding Analysis
2. Activity Analysis
3. Safety Analysis
4. Property Analysis
5. SAR Analysis

Basic Usage
---------

Using the analysis pipeline:

.. code-block:: python

    from binding_data_processor.pipeline import ProcessingPipeline
    from binding_data_processor.pipeline.config import ProcessingConfig

    # Create pipeline with analysis
    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            use_analysis=True,
            analysis_dir="analysis/",
        )
    )

    # Process compounds
    compounds = pipeline.process_compounds(
        input_file="compounds.tsv",
        output_dir="results/",
    )

Binding Analysis
-------------

Analyzing receptor binding:

.. code-block:: python

    from binding_data_processor.models.compound.analysis import BindingAnalyzer

    # Configure analyzer
    analyzer = BindingAnalyzer(
        target_receptors=["5-HT2A", "NMDA"],
        min_affinity=1e-6,
    )

    # Analyze binding
    for compound in compounds:
        results = analyzer.analyze_binding(compound)
        print(f"Compound: {compound.name}")
        print(f"Strongest Target: {results.strongest_target}")
        print(f"Binding Profile:")
        for receptor, data in results.binding_profile.items():
            print(f"- {receptor}: {data.affinity} ({data.activity})")

Activity Analysis
--------------

Analyzing compound activity:

.. code-block:: python

    from binding_data_processor.models.compound.analysis import ActivityAnalyzer

    # Configure analyzer
    analyzer = ActivityAnalyzer(
        effect_types=["stimulant", "psychedelic"],
        min_confidence=0.7,
    )

    # Analyze activity
    for compound in compounds:
        results = analyzer.analyze_activity(compound)
        print(f"Compound: {compound.name}")
        print(f"Primary Effect: {results.primary_effect}")
        print(f"Effect Profile:")
        for effect, score in results.effect_profile.items():
            print(f"- {effect}: {score:.2f}")

Safety Analysis
------------

Analyzing safety profiles:

.. code-block:: python

    from binding_data_processor.models.compound.analysis import SafetyAnalyzer

    # Configure analyzer
    analyzer = SafetyAnalyzer(
        risk_types=["toxicity", "addiction"],
        data_sources=["predicted", "literature", "community"],
    )

    # Analyze safety
    for compound in compounds:
        results = analyzer.analyze_safety(compound)
        print(f"Compound: {compound.name}")
        print(f"Overall Risk: {results.risk_level}")
        print(f"Risk Factors:")
        for factor, risk in results.risk_factors.items():
            print(f"- {factor}: {risk.level} ({risk.confidence:.2f})")

Property Analysis
--------------

Analyzing molecular properties:

.. code-block:: python

    from binding_data_processor.models.compound.analysis import PropertyAnalyzer

    # Configure analyzer
    analyzer = PropertyAnalyzer(
        properties=["mw", "logp", "tpsa", "hbd", "hba"],
        include_descriptors=True,
    )

    # Analyze properties
    for compound in compounds:
        results = analyzer.analyze_properties(compound)
        print(f"Compound: {compound.name}")
        print(f"Properties:")
        for name, value in results.properties.items():
            print(f"- {name}: {value}")
        print(f"Drug-likeness: {results.druglikeness_score}")

SAR Analysis
----------

Structure-activity relationship analysis:

.. code-block:: python

    from binding_data_processor.models.compound.analysis import SARAnalyzer

    # Configure analyzer
    analyzer = SARAnalyzer(
        similarity_threshold=0.7,
        activity_threshold=1.0,
    )

    # Analyze SAR
    results = analyzer.analyze_series(compounds)
    for series in results.series:
        print(f"Series: {series.name}")
        print(f"Core: {series.core_structure}")
        print(f"Activity Cliffs:")
        for cliff in series.activity_cliffs:
            print(f"- {cliff.compound_pair}")
            print(f"  Similarity: {cliff.similarity:.2f}")
            print(f"  Activity Difference: {cliff.activity_diff:.2f}")

Advanced Analysis
--------------

Custom Analysis
~~~~~~~~~~~~

Creating custom analyzers:

.. code-block:: python

    from binding_data_processor.models.compound.analysis import BaseAnalyzer

    class CustomAnalyzer(BaseAnalyzer):
        def __init__(self, **kwargs):
            super().__init__()
            self.config = kwargs

        def analyze(self, compound):
            # Custom analysis logic
            results = self._analyze_data(compound)
            return self._format_results(results)

    # Use custom analyzer
    analyzer = CustomAnalyzer(param1="value1")
    results = analyzer.analyze(compound)

Ensemble Analysis
~~~~~~~~~~~~~~

Combining multiple analyses:

.. code-block:: python

    from binding_data_processor.models.compound.analysis import AnalysisEnsemble

    # Configure ensemble
    ensemble = AnalysisEnsemble([
        BindingAnalyzer(),
        ActivityAnalyzer(),
        SafetyAnalyzer(),
    ])

    # Run ensemble analysis
    results = ensemble.analyze(compound)
    print(f"Binding Results: {results.binding}")
    print(f"Activity Results: {results.activity}")
    print(f"Safety Results: {results.safety}")

Statistical Analysis
~~~~~~~~~~~~~~~~

Analyzing compound series:

.. code-block:: python

    from binding_data_processor.models.compound.analysis import StatsAnalyzer

    # Configure analyzer
    analyzer = StatsAnalyzer(
        metrics=["mean", "std", "correlation"],
        grouping="scaffold",
    )

    # Analyze statistics
    stats = analyzer.analyze_series(compounds)
    print(f"Series Statistics:")
    for scaffold, data in stats.items():
        print(f"\nScaffold: {scaffold}")
        print(f"Compounds: {len(data.compounds)}")
        print(f"Mean Activity: {data.mean_activity:.2f}")
        print(f"Activity StdDev: {data.activity_std:.2f}")

Visualization
-----------

Creating analysis visualizations:

.. code-block:: python

    from binding_data_processor.models.compound.analysis import (
        AnalysisVisualizer,
        PlotConfig,
    )

    # Configure visualizer
    visualizer = AnalysisVisualizer(
        config=PlotConfig(
            plot_type="interactive",
            width=800,
            height=600,
        )
    )

    # Create plots
    plots = visualizer.create_plots(results)
    for name, plot in plots.items():
        plot.write_html(f"{name}.html")

Report Generation
--------------

Generating analysis reports:

.. code-block:: python

    from binding_data_processor.models.compound.analysis import ReportGenerator

    # Configure generator
    generator = ReportGenerator(
        sections=[
            "binding",
            "activity",
            "safety",
            "properties",
        ],
        include_plots=True,
    )

    # Generate report
    report = generator.generate_report(
        compounds=compounds,
        analyses=results,
        output_file="report.html",
    )
