Troubleshooting Guide
==================

This guide helps you diagnose and fix common issues you might encounter while using ChemData.

Installation Issues
----------------

ImportError: No module named 'rdkit'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RDKit installation can fail on some systems. Try:

1. Using conda:

.. code-block:: bash

    conda install -c conda-forge rdkit

2. If conda doesn't work, build from source:

.. code-block:: bash

    git clone https://github.com/rdkit/rdkit.git
    cd rdkit
    mkdir build && cd build
    cmake ..
    make -j4
    make install

ImportError: DLL load failed (Windows)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Missing Visual C++ Redistributable. Install from:
https://support.microsoft.com/en-us/help/2977003/the-latest-supported-visual-c-downloads

Data Processing Issues
-------------------

MemoryError During Processing
~~~~~~~~~~~~~~~~~~~~~~~~~

1. Enable batch processing:

.. code-block:: python

    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            batch_size=1000,
            use_disk_cache=True,
        )
    )

2. Reduce memory usage:

.. code-block:: python

    import os
    os.environ["CHEMDATA_LOW_MEMORY"] = "1"

Invalid SMILES Strings
~~~~~~~~~~~~~~~~~~

1. Enable structure validation:

.. code-block:: python

    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            validate_structures=True,
            standardize_structures=True,
        )
    )

2. Handle invalid structures:

.. code-block:: python

    try:
        compounds = pipeline.process_compounds(input_file)
    except InvalidStructureError as e:
        print(f"Invalid structure: {e.smiles}")
        print(f"Reason: {e.reason}")

Missing Data Fields
~~~~~~~~~~~~~~~

1. Set default values:

.. code-block:: python

    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            default_values={
                "binding_affinity": None,
                "confidence": 0.0,
            },
        )
    )

2. Filter incomplete records:

.. code-block:: python

    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            require_complete=True,
            required_fields=["smiles", "cas_number"],
        )
    )

ML Issues
--------

CUDA Out of Memory
~~~~~~~~~~~~~~~

1. Reduce batch size:

.. code-block:: python

    predictor = BindingPredictor(
        config=PredictorConfig(
            batch_size=32,  # Reduce from default
            precision="mixed",  # Use mixed precision
        )
    )

2. Use CPU fallback:

.. code-block:: python

    import os
    os.environ["CUDA_VISIBLE_DEVICES"] = ""  # Disable GPU

Model Loading Failed
~~~~~~~~~~~~~~~~

1. Check model paths:

.. code-block:: python

    predictor = BindingPredictor(
        config=PredictorConfig(
            model_dir="models/binding/",
            version="latest",
        )
    )

2. Download missing models:

.. code-block:: python

    from binding_data_processor.models import download_models

    download_models(
        model_types=["binding", "activity"],
        version="latest",
    )

Web Enrichment Issues
------------------

API Rate Limits
~~~~~~~~~~~~

1. Enable rate limiting:

.. code-block:: python

    enricher = WebEnrichmentManager(
        config=EnrichmentConfig(
            rate_limit=1.0,  # requests per second
            use_cache=True,
        )
    )

2. Use multiple API keys:

.. code-block:: python

    enricher = WebEnrichmentManager(
        config=EnrichmentConfig(
            api_keys={
                "twitter": ["key1", "key2", "key3"],
                "reddit": ["key1", "key2", "key3"],
            },
        )
    )

Connection Timeouts
~~~~~~~~~~~~~~~

1. Adjust timeouts:

.. code-block:: python

    enricher = WebEnrichmentManager(
        config=EnrichmentConfig(
            timeout=30,  # seconds
            retries=3,
        )
    )

2. Enable circuit breaker:

.. code-block:: python

    from binding_data_processor.pipeline.infrastructure import CircuitBreaker

    enricher = WebEnrichmentManager(
        config=EnrichmentConfig(
            circuit_breaker=CircuitBreaker(
                failure_threshold=5,
                recovery_time=60,
            ),
        )
    )

Cache Issues
---------

Cache Corruption
~~~~~~~~~~~~~

1. Clear cache:

.. code-block:: python

    from binding_data_processor.pipeline import clear_cache

    clear_cache(
        cache_types=["web", "predictions", "structures"],
        older_than="7d",
    )

2. Verify cache:

.. code-block:: python

    from binding_data_processor.pipeline import verify_cache

    verify_cache(
        repair=True,
        remove_invalid=True,
    )

Cache Space Issues
~~~~~~~~~~~~~~

1. Set cache limits:

.. code-block:: python

    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            max_cache_size="10GB",
            cache_ttl="30d",
        )
    )

2. Use selective caching:

.. code-block:: python

    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            cache_predictions=True,
            cache_structures=False,
            cache_web_data=True,
        )
    )

Performance Issues
---------------

Slow Processing
~~~~~~~~~~~~

1. Enable parallel processing:

.. code-block:: python

    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            num_workers=4,
            use_threading=True,
        )
    )

2. Profile performance:

.. code-block:: python

    from binding_data_processor.pipeline import profile_pipeline

    profile_pipeline(
        pipeline=pipeline,
        input_file="compounds.tsv",
        output_file="profile.json",
    )

High Memory Usage
~~~~~~~~~~~~~

1. Enable streaming:

.. code-block:: python

    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            stream_processing=True,
            chunk_size=1000,
        )
    )

2. Use memory monitoring:

.. code-block:: python

    from binding_data_processor.pipeline import MemoryMonitor

    with MemoryMonitor(threshold="8GB"):
        pipeline.process_compounds(input_file)

Getting Help
----------

1. Check logs:

.. code-block:: python

    from binding_data_processor import get_logs

    logs = get_logs(
        level="ERROR",
        days=7,
        include_tracebacks=True,
    )

2. Generate diagnostic report:

.. code-block:: python

    from binding_data_processor import generate_report

    report = generate_report(
        include_system_info=True,
        include_logs=True,
        include_config=True,
    )
    report.save("diagnostic_report.txt")
