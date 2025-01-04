Data Processing Guide
==================

This guide covers the data processing capabilities of ChemData in detail.

Pipeline Overview
--------------

The data processing pipeline consists of several stages:

1. Data Loading
2. Structure Validation
3. Property Calculation
4. Web Enrichment
5. Analysis
6. Export

Basic Usage
---------

The simplest way to process compounds:

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

Data Sources
----------

BindingDB
~~~~~~~~

Loading data from BindingDB:

.. code-block:: python

    from binding_data_processor.data_sources import BindingDBLoader

    # Configure loader
    loader = BindingDBLoader(
        min_binding_affinity=1e-6,
        target_receptors=["5-HT2A", "NMDA"],
    )

    # Load compounds
    compounds = loader.load_compounds("bindingdb.tsv")

ChEMBL
~~~~~

Loading data from ChEMBL:

.. code-block:: python

    from binding_data_processor.data_sources import ChEMBLLoader

    # Configure loader
    loader = ChEMBLLoader(
        target_types=["SINGLE PROTEIN"],
        assay_types=["B"],
    )

    # Load compounds
    compounds = loader.load_compounds_by_target("CHEMBL1983")

PubChem
~~~~~~

Loading data from PubChem:

.. code-block:: python

    from binding_data_processor.data_sources import PubChemLoader

    # Configure loader
    loader = PubChemLoader()

    # Load compounds
    compounds = loader.load_compounds_by_query(
        "dopamine antagonist",
        max_compounds=1000,
    )

Structure Processing
-----------------

Validation
~~~~~~~~~

Validating chemical structures:

.. code-block:: python

    from binding_data_processor.processors.structure import (
        StructureValidator,
        ValidationConfig,
    )

    # Configure validator
    validator = StructureValidator(
        config=ValidationConfig(
            check_valence=True,
            check_aromaticity=True,
            standardize=True,
        )
    )

    # Validate structures
    for compound in compounds:
        result = validator.validate(compound)
        if not result.valid:
            print(f"Invalid structure: {result.errors}")

Property Calculation
~~~~~~~~~~~~~~~~~

Calculating molecular properties:

.. code-block:: python

    from binding_data_processor.processors.structure.properties import (
        PropertyCalculator,
        PropertyConfig,
    )

    # Configure calculator
    calculator = PropertyCalculator(
        config=PropertyConfig(
            calc_descriptors=True,
            calc_fingerprints=True,
        )
    )

    # Calculate properties
    for compound in compounds:
        calculator.calculate_properties(compound)

Similarity Search
~~~~~~~~~~~~~~

Finding similar compounds:

.. code-block:: python

    from binding_data_processor.processors.structure.similarity import (
        SimilaritySearcher,
        SearchConfig,
    )

    # Configure searcher
    searcher = SimilaritySearcher(
        config=SearchConfig(
            similarity_threshold=0.7,
            max_results=100,
        )
    )

    # Search similar compounds
    query = compounds[0]
    similar = searcher.find_similar(query, compounds)

Web Enrichment
------------

Community Data
~~~~~~~~~~~~

Gathering community data:

.. code-block:: python

    from binding_data_processor.web_enrichment import (
        CommunityDataEnricher,
        EnrichmentConfig,
    )

    # Configure enricher
    enricher = CommunityDataEnricher(
        config=EnrichmentConfig(
            use_psychonautwiki=True,
            use_erowid=True,
            use_tripsit=True,
        )
    )

    # Enrich compounds
    for compound in compounds:
        enricher.enrich_compound(compound)

Social Media
~~~~~~~~~~

Monitoring social media:

.. code-block:: python

    from binding_data_processor.web_enrichment import SocialMediaMonitor

    # Configure monitor
    monitor = SocialMediaMonitor(
        subreddits=[
            "researchchemicals",
            "nootropics",
            "DrugNerds",
        ],
        use_twitter=True,
        use_bluesky=True,
    )

    # Monitor compounds
    for compound in compounds:
        mentions = monitor.get_mentions(compound)
        compound.social_data = mentions

Patent Search
~~~~~~~~~~~

Searching patents:

.. code-block:: python

    from binding_data_processor.web_enrichment import PatentSearcher

    # Configure searcher
    searcher = PatentSearcher(
        search_recent=True,
        max_results=1000,
    )

    # Search patents
    for compound in compounds:
        patents = searcher.search_patents(compound)
        compound.patent_data = patents

Data Export
---------

TSV Export
~~~~~~~~

Exporting to TSV:

.. code-block:: python

    from binding_data_processor.pipeline.export import TSVExporter

    # Configure exporter
    exporter = TSVExporter(
        columns=[
            "name",
            "smiles",
            "cas_number",
            "binding_affinity",
            "bbb_prediction",
        ],
    )

    # Export compounds
    exporter.export_compounds(
        compounds,
        "compounds.tsv",
    )

JSON Export
~~~~~~~~~

Exporting to JSON:

.. code-block:: python

    from binding_data_processor.pipeline.export import JSONExporter

    # Configure exporter
    exporter = JSONExporter(
        include_predictions=True,
        include_web_data=True,
    )

    # Export compounds
    exporter.export_compounds(
        compounds,
        "compounds.json",
    )

Advanced Usage
------------

Custom Pipeline
~~~~~~~~~~~~

Creating a custom pipeline:

.. code-block:: python

    from binding_data_processor.pipeline import (
        Pipeline,
        Stage,
        Config,
    )

    # Define custom stage
    class CustomStage(Stage):
        def process(self, compounds):
            # Custom processing
            return compounds

    # Create pipeline
    pipeline = Pipeline([
        StructureValidationStage(),
        PropertyCalculationStage(),
        CustomStage(),
        WebEnrichmentStage(),
        ExportStage(),
    ])

    # Process compounds
    compounds = pipeline.process(compounds)

Batch Processing
~~~~~~~~~~~~~

Processing compounds in batches:

.. code-block:: python

    from binding_data_processor.pipeline import BatchProcessor

    # Configure processor
    processor = BatchProcessor(
        batch_size=100,
        num_workers=4,
    )

    # Process batches
    for batch in processor.process_batches(compounds):
        # Handle batch results
        pass

Error Handling
~~~~~~~~~~~

Handling processing errors:

.. code-block:: python

    from binding_data_processor.pipeline import (
        ErrorHandler,
        ProcessingError,
    )

    # Configure handler
    handler = ErrorHandler(
        retry_count=3,
        ignore_errors=False,
    )

    # Process with error handling
    try:
        compounds = pipeline.process_compounds(
            input_file="bindingdb.tsv",
            error_handler=handler,
        )
    except ProcessingError as e:
        print(f"Processing failed: {e}")

Checkpointing
~~~~~~~~~~~

Using checkpoints:

.. code-block:: python

    from binding_data_processor.pipeline import CheckpointManager

    # Configure checkpoints
    checkpoints = CheckpointManager(
        checkpoint_dir="checkpoints/",
        save_frequency=1000,
    )

    # Process with checkpoints
    compounds = pipeline.process_compounds(
        input_file="bindingdb.tsv",
        checkpoint_manager=checkpoints,
    )
