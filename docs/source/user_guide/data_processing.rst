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

Document Processing
----------------

The document processing system allows you to upload and process PDF documents to extract compound information. The system supports both single file uploads and batch processing with progress tracking.

Single File Upload
~~~~~~~~~~~~~~~~

To upload and process a single PDF file:

.. code-block:: python

    import requests

    # Upload single file
    files = {'files': ('document.pdf', open('document.pdf', 'rb'))}
    response = requests.post('http://localhost:8000/documents/upload', files=files)
    
    # Check result
    result = response.json()[0]
    print(f"Status: {result['status']}")
    print(f"Compounds found: {result['compounds_found']}")

Batch Upload
~~~~~~~~~~

For processing multiple files in a batch:

.. code-block:: python

    import requests
    import time
    from uuid import UUID

    # Upload multiple files
    files = [
        ('files', ('doc1.pdf', open('doc1.pdf', 'rb'))),
        ('files', ('doc2.pdf', open('doc2.pdf', 'rb'))),
        ('files', ('doc3.pdf', open('doc3.pdf', 'rb')))
    ]
    response = requests.post('http://localhost:8000/documents/batch-upload', files=files)
    batch_id = UUID(response.json()['batch_id'])

    # Track progress
    while True:
        status = requests.get(f'http://localhost:8000/documents/batch-status/{batch_id}').json()
        print(f"Progress: {status['processed_files']}/{status['total_files']}")
        
        if status['status'] in ['completed', 'failed']:
            break
            
        time.sleep(1)  # Poll every second

    # Check results
    for filename, file_status in status['files'].items():
        print(f"{filename}: {file_status['status']}")
        if file_status['compounds_found']:
            print(f"Found {file_status['compounds_found']} compounds")

Directory Monitoring
~~~~~~~~~~~~~~~~~

To monitor a directory for new PDF files:

.. code-block:: python

    import requests

    # Configure directory monitoring
    config = {
        'path': '/path/to/documents',
        'patterns': ['*.pdf'],
        'recursive': True
    }
    response = requests.post('http://localhost:8000/documents/monitor', json=config)

    # Stop monitoring when done
    requests.delete('/documents/monitor', params={'path': '/path/to/documents'})

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

Error Handling
~~~~~~~~~~~

The system provides detailed error information:

- For single file uploads, check the ``status`` and ``error_message`` fields in the response
- For batch uploads, check the overall batch status and individual file statuses
- For directory monitoring, check the response status and message

Example error handling:

.. code-block:: python

    import requests

    try:
        response = requests.post('/documents/upload', files={'files': ('doc.pdf', open('doc.pdf', 'rb'))})
        result = response.json()[0]
        
        if result['status'] == 'failed':
            print(f"Processing failed: {result['error_message']}")
        elif 'error_message' in result:
            print(f"Warning: {result['error_message']}")
            
    except requests.exceptions.RequestException as e:
        print(f"Upload failed: {str(e)}")

Best Practices
------------

1. Use batch processing for multiple files to improve performance
2. Monitor batch progress to track processing status
3. Handle errors appropriately at both file and batch levels
4. Use directory monitoring for automated processing
5. Clean up resources (close files, stop monitoring) when done
6. Use checkpoints for long-running processes
7. Implement proper error handling and validation
8. Configure appropriate logging levels
9. Use appropriate data formats for export
10. Follow memory management guidelines

Configuration
-----------

The processing system can be configured through environment variables:

- ``STORAGE_DIR``: Base directory for storing documents (default: "data/documents")
- ``MAX_UPLOAD_SIZE``: Maximum file size in bytes (default: 10MB)
- ``SUPPORTED_FORMATS``: List of supported file formats (default: ["pdf"])
- ``PROCESSING_THREADS``: Number of processing threads (default: 4)
- ``CACHE_DIR``: Directory for caching results (default: "cache/")
- ``LOG_LEVEL``: Logging level (default: "INFO")
- ``BATCH_SIZE``: Processing batch size (default: 1000)
- ``MAX_RETRIES``: Maximum retry attempts (default: 3)

Example configuration:

.. code-block:: bash

    export STORAGE_DIR=/path/to/storage
    export MAX_UPLOAD_SIZE=20971520  # 20MB
    export SUPPORTED_FORMATS='["pdf", "PDF"]'
    export PROCESSING_THREADS=8
    export CACHE_DIR=/path/to/cache
    export LOG_LEVEL=DEBUG
    export BATCH_SIZE=500
    export MAX_RETRIES=5
