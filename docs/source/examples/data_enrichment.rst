Data Enrichment Workflows
=====================

This guide shows how to implement data enrichment workflows to gather and integrate data from various sources.

Community Data
------------

Gathering data from community sources:

.. code-block:: python

    from binding_data_processor.web_enrichment import (
        CommunityDataEnricher,
        EnrichmentConfig,
    )

    # Configure enricher
    enricher = CommunityDataEnricher(
        config=EnrichmentConfig(
            sources=[
                "psychonautwiki",
                "erowid",
                "tripsit",
            ],
            cache_dir="cache/community/",
        )
    )

    # Enrich compounds
    for compound in compounds:
        # Get community data
        data = enricher.enrich_compound(compound)
        
        # Extract information
        effects = data.get("effects", [])
        dosage = data.get("dosage", {})
        duration = data.get("duration", {})
        
        # Update compound
        compound.community_data = {
            "effects": effects,
            "dosage": dosage,
            "duration": duration,
        }

Social Media
----------

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
        twitter_keywords=[
            "new compound",
            "synthesis",
            "receptor",
        ],
        cache_dir="cache/social/",
    )

    # Monitor compounds
    for compound in compounds:
        # Get social data
        data = monitor.get_compound_mentions(compound)
        
        # Process Reddit data
        reddit_data = {
            "posts": data["reddit"]["posts"],
            "comments": data["reddit"]["comments"],
            "sentiment": data["reddit"]["sentiment"],
        }
        
        # Process Twitter data
        twitter_data = {
            "tweets": data["twitter"]["tweets"],
            "retweets": data["twitter"]["retweets"],
            "sentiment": data["twitter"]["sentiment"],
        }
        
        # Update compound
        compound.social_data = {
            "reddit": reddit_data,
            "twitter": twitter_data,
        }

Scientific Literature
-----------------

Processing scientific literature:

.. code-block:: python

    from binding_data_processor.web_enrichment import LiteratureProcessor

    # Configure processor
    processor = LiteratureProcessor(
        sources=[
            "pubmed",
            "patents",
            "google_scholar",
        ],
        cache_dir="cache/literature/",
    )

    # Process compounds
    for compound in compounds:
        # Get literature data
        data = processor.process_compound(compound)
        
        # Extract papers
        papers = data.get("papers", [])
        citations = data.get("citations", [])
        
        # Extract patents
        patents = data.get("patents", [])
        claims = data.get("claims", [])
        
        # Update compound
        compound.literature_data = {
            "papers": papers,
            "citations": citations,
            "patents": patents,
            "claims": claims,
        }

Chemical Databases
---------------

Querying chemical databases:

.. code-block:: python

    from binding_data_processor.web_enrichment import DatabaseEnricher

    # Configure enricher
    enricher = DatabaseEnricher(
        databases=[
            "chembl",
            "pubchem",
            "drugbank",
        ],
        cache_dir="cache/databases/",
    )

    # Enrich compounds
    for compound in compounds:
        # Get database data
        data = enricher.enrich_compound(compound)
        
        # Extract ChEMBL data
        chembl_data = {
            "activities": data["chembl"]["activities"],
            "targets": data["chembl"]["targets"],
            "assays": data["chembl"]["assays"],
        }
        
        # Extract PubChem data
        pubchem_data = {
            "properties": data["pubchem"]["properties"],
            "bioactivity": data["pubchem"]["bioactivity"],
            "synonyms": data["pubchem"]["synonyms"],
        }
        
        # Update compound
        compound.database_data = {
            "chembl": chembl_data,
            "pubchem": pubchem_data,
        }

Data Integration
-------------

Integrating data from multiple sources:

.. code-block:: python

    from binding_data_processor.web_enrichment import DataIntegrator

    # Configure integrator
    integrator = DataIntegrator(
        enrichers=[
            CommunityDataEnricher(),
            SocialMediaMonitor(),
            LiteratureProcessor(),
            DatabaseEnricher(),
        ],
        cache_dir="cache/integrated/",
    )

    # Process compounds
    for compound in compounds:
        # Enrich from all sources
        data = integrator.enrich_compound(compound)
        
        # Merge data
        compound.merge_enrichment_data(data)
        
        # Validate data
        integrator.validate_data(compound)
        
        # Save enriched compound
        compound.save()

LLM Processing
-----------

Using LLMs to process text data:

.. code-block:: python

    from binding_data_processor.web_enrichment import LLMProcessor

    # Configure processor
    processor = LLMProcessor(
        model="gpt-4",
        cache_dir="cache/llm/",
    )

    # Process compounds
    for compound in compounds:
        # Extract text data
        text_data = {
            "community": compound.community_data,
            "social": compound.social_data,
            "literature": compound.literature_data,
        }
        
        # Process with LLM
        results = processor.process_text_data(
            text_data,
            tasks=[
                "extract_effects",
                "extract_mechanisms",
                "extract_safety",
            ],
        )
        
        # Update compound
        compound.llm_data = results

Custom Enrichment
--------------

Creating custom enrichment sources:

.. code-block:: python

    from binding_data_processor.web_enrichment import BaseEnricher

    class CustomEnricher(BaseEnricher):
        """Custom data enrichment source."""

        def __init__(self, config: EnrichmentConfig):
            super().__init__(config)
            self.client = self._create_client()

        def enrich_compound(self, compound: CompoundData) -> Dict[str, Any]:
            """Enrich compound with custom data."""
            try:
                # Get custom data
                data = self.client.get_compound_data(compound)
                
                # Process data
                processed = self._process_data(data)
                
                # Validate
                self._validate_data(processed)
                
                return processed
                
            except Exception as e:
                self.logger.error(f"Enrichment error: {str(e)}")
                return {}

        def _process_data(self, data: Dict[str, Any]) -> Dict[str, Any]:
            """Process custom data."""
            # Custom processing logic
            return processed_data

        def _validate_data(self, data: Dict[str, Any]) -> None:
            """Validate custom data."""
            # Custom validation logic
            if not self._is_valid(data):
                raise ValidationError("Invalid data")

Enrichment Pipeline
----------------

Creating an enrichment pipeline:

.. code-block:: python

    from binding_data_processor.pipeline import EnrichmentPipeline

    # Configure pipeline
    config = EnrichmentConfig(
        enrichers=[
            CommunityDataEnricher(),
            SocialMediaMonitor(),
            LiteratureProcessor(),
            DatabaseEnricher(),
            CustomEnricher(),
        ],
        cache_enabled=True,
        batch_size=100,
        num_workers=4,
    )

    # Create pipeline
    pipeline = EnrichmentPipeline(config)

    # Process compounds
    enriched_compounds = pipeline.process_compounds(
        input_file="compounds.tsv",
        output_file="enriched_compounds.tsv",
    )

    # Get statistics
    stats = pipeline.get_stats()
    print(f"Enrichment Stats: {stats}")
