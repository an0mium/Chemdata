Quickstart Guide
==============

This guide will help you get started with ChemData quickly. We'll cover:

1. Basic installation
2. Processing compounds from BindingDB
3. Running predictions
4. Using the web interface

Basic Installation
----------------

Install ChemData using pip:

.. code-block:: bash

    # Create virtual environment
    python -m venv .venv
    source .venv/bin/activate  # Linux/macOS
    # or
    .venv\\Scripts\\activate  # Windows

    # Install package
    pip install chemdata

Process Compounds
---------------

Let's process some compounds from BindingDB:

.. code-block:: python

    from binding_data_processor.pipeline import ProcessingPipeline
    from binding_data_processor.pipeline.config import ProcessingConfig

    # Create pipeline
    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            use_ml_predictions=True,
            use_web_enrichment=True,
        )
    )

    # Process compounds
    compounds = pipeline.process_compounds(
        input_file="bindingdb.tsv",
        output_dir="results/",
    )

    # Print results
    for compound in compounds:
        print(f"Name: {compound.name}")
        print(f"SMILES: {compound.smiles}")
        print(f"Binding Data: {compound.binding_data}")
        print()

Run Predictions
-------------

Make predictions for compounds:

.. code-block:: python

    from binding_data_processor.processors.psychopharm.predictors.bbb import (
        BBBPredictorWebEnriched
    )

    # Initialize predictor
    predictor = BBBPredictorWebEnriched(
        model_dir="models/bbb",
        cache_dir="cache",
    )

    # Make predictions
    for compound in compounds:
        result = predictor.predict(compound)
        print(f"Compound: {compound.name}")
        print(f"BBB Class: {result.value}")
        print(f"Confidence: {result.confidence:.2f}")
        print()

Web Interface
-----------

Launch the web interface:

.. code-block:: bash

    # Start web app
    streamlit run examples/web_app/app.py

This will open a browser window with the ChemData web interface where you can:

- Browse compounds
- Search by structure or properties
- View detailed compound information
- Export data

Command Line Interface
-------------------

Process compounds using the CLI:

.. code-block:: bash

    # Process compounds
    python -m binding_data_processor.cli process-compounds \
        --input bindingdb.tsv \
        --output results/ \
        --enable-ml \
        --enable-web \
        --enable-social

    # Run predictions
    python -m binding_data_processor.cli predict-bbb \
        --input results/compounds.tsv \
        --output results/predictions.tsv

    # Generate report
    python -m binding_data_processor.cli generate-report \
        --input results/predictions.tsv \
        --output report.html

Next Steps
---------

- Read the :doc:`user_guide/index` for detailed usage
- Check the :doc:`api_reference/index` for API details
- See :doc:`examples/index` for more examples
- Read :doc:`contributing` to contribute

Example Projects
--------------

1. Process NMDA Antagonists:

.. code-block:: python

    from binding_data_processor.pipeline import ProcessingPipeline
    from pathlib import Path

    # Configure pipeline
    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            target_receptors=["NMDA"],
            min_binding_affinity=1e-6,
            use_ml_predictions=True,
            use_web_enrichment=True,
        )
    )

    # Process compounds
    compounds = pipeline.process_compounds(
        input_file=Path("bindingdb.tsv"),
        output_dir=Path("results/nmda/"),
    )

2. Analyze 5-HT2 Ligands:

.. code-block:: python

    # Configure pipeline
    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            target_receptors=["5-HT2A", "5-HT2B", "5-HT2C"],
            min_binding_affinity=1e-7,
            use_ml_predictions=True,
            use_web_enrichment=True,
            use_social_monitoring=True,
        )
    )

    # Process compounds
    compounds = pipeline.process_compounds(
        input_file=Path("bindingdb.tsv"),
        output_dir=Path("results/5ht2/"),
    )

3. Search for Novel Compounds:

.. code-block:: python

    from binding_data_processor.web_enrichment import WebEnrichmentManager

    # Configure enrichment
    enrichment = WebEnrichmentManager(
        monitor_subreddits=[
            "researchchemicals",
            "nootropics",
            "DrugNerds",
        ],
        monitor_twitter=True,
        monitor_patents=True,
    )

    # Search for compounds
    compounds = enrichment.search_compounds(
        keywords=["novel", "synthesis", "receptor"],
        date_range=("2023-01-01", "2024-01-01"),
    )

Common Tasks
----------

1. Export Data:

.. code-block:: python

    # Export to TSV
    pipeline.export_compounds(
        compounds,
        output_file="compounds.tsv",
        format="tsv",
        columns=[
            "name",
            "smiles",
            "cas_number",
            "binding_affinity",
            "bbb_prediction",
        ],
    )

2. Filter Compounds:

.. code-block:: python

    # Filter by properties
    filtered = [
        c for c in compounds
        if c.molecular_weight < 500
        and c.logp < 5
        and c.bbb_prediction.value == "BBB+"
    ]

3. Visualize Data:

.. code-block:: python

    from binding_data_processor.web.components import CompoundDetails

    # Create visualizer
    visualizer = CompoundDetails()

    # Generate plots
    for compound in compounds:
        visualizer.set_compound(compound)
        plots = visualizer.get_plot_data()
        
        # Save plots
        for name, plot in plots.items():
            plot.write_html(f"{compound.name}_{name}.html")
