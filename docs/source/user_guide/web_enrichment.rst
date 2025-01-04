Web Enrichment Guide
==================

This guide covers the web enrichment capabilities of ChemData in detail.

Overview
-------

The web enrichment pipeline gathers data from:

1. Community Sources (PsychonautWiki, Erowid, TripSit)
2. Social Media (Reddit, Twitter, Bluesky)
3. Scientific Sources (PubMed, Patents)
4. Chemical Databases (ChEMBL, PubChem)

Basic Usage
---------

Using web enrichment:

.. code-block:: python

    from binding_data_processor.pipeline import ProcessingPipeline
    from binding_data_processor.pipeline.config import ProcessingConfig

    # Create pipeline with web enrichment
    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            use_web_enrichment=True,
            cache_dir="cache/",
        )
    )

    # Process compounds
    compounds = pipeline.process_compounds(
        input_file="compounds.tsv",
        output_dir="results/",
    )

Community Data
------------

PsychonautWiki
~~~~~~~~~~~~

Gathering data from PsychonautWiki:

.. code-block:: python

    from binding_data_processor.web_enrichment import PsychonautWikiClient

    # Configure client
    client = PsychonautWikiClient(
        cache_dir="cache/psychonautwiki",
        rate_limit=1.0,  # requests per second
    )

    # Get compound data
    data = client.get_compound_data("LSD")
    print(f"Effects: {data['effects']}")
    print(f"Duration: {data['duration']}")
    print(f"Dosage: {data['dosage']}")

Erowid
~~~~~

Analyzing Erowid experience reports:

.. code-block:: python

    from binding_data_processor.web_enrichment import ErowidClient

    # Configure client
    client = ErowidClient(
        cache_dir="cache/erowid",
        max_reports=100,
    )

    # Get experience reports
    reports = client.get_experience_reports("MDMA")
    for report in reports:
        print(f"Title: {report.title}")
        print(f"Effects: {report.effects}")
        print(f"Body: {report.text[:200]}...")

TripSit
~~~~~~

Getting factsheet data:

.. code-block:: python

    from binding_data_processor.web_enrichment import TripSitClient

    # Configure client
    client = TripSitClient(
        cache_dir="cache/tripsit",
    )

    # Get factsheet
    data = client.get_factsheet("ketamine")
    print(f"Summary: {data['summary']}")
    print(f"Categories: {data['categories']}")
    print(f"Interactions: {data['interactions']}")

Social Media
----------

Reddit
~~~~~

Monitoring Reddit discussions:

.. code-block:: python

    from binding_data_processor.web_enrichment import RedditMonitor

    # Configure monitor
    monitor = RedditMonitor(
        subreddits=[
            "researchchemicals",
            "nootropics",
            "DrugNerds",
        ],
        client_id="your_client_id",
        client_secret="your_client_secret",
    )

    # Get recent discussions
    posts = monitor.get_recent_posts("2C-B", limit=10)
    for post in posts:
        print(f"Title: {post.title}")
        print(f"Score: {post.score}")
        print(f"Comments: {len(post.comments)}")

Twitter
~~~~~~

Monitoring Twitter mentions:

.. code-block:: python

    from binding_data_processor.web_enrichment import TwitterMonitor

    # Configure monitor
    monitor = TwitterMonitor(
        api_key="your_api_key",
        api_secret="your_api_secret",
        keywords=["new compound", "synthesis", "receptor"],
    )

    # Get recent mentions
    tweets = monitor.get_recent_tweets("5-MeO-DMT")
    for tweet in tweets:
        print(f"Text: {tweet.text}")
        print(f"Likes: {tweet.like_count}")
        print(f"Retweets: {tweet.retweet_count}")

Bluesky
~~~~~~

Monitoring Bluesky posts:

.. code-block:: python

    from binding_data_processor.web_enrichment import BlueskyMonitor

    # Configure monitor
    monitor = BlueskyMonitor(
        handle="your_handle",
        password="your_password",
    )

    # Get recent posts
    posts = monitor.get_recent_posts("phenethylamine")
    for post in posts:
        print(f"Text: {post.text}")
        print(f"Likes: {post.like_count}")
        print(f"Reposts: {post.repost_count}")

Scientific Sources
---------------

PubMed
~~~~~

Searching scientific literature:

.. code-block:: python

    from binding_data_processor.web_enrichment import PubMedClient

    # Configure client
    client = PubMedClient(
        email="your_email",
        api_key="your_api_key",
    )

    # Search papers
    papers = client.search_papers(
        query="5-HT2A antagonist",
        max_results=100,
    )
    for paper in papers:
        print(f"Title: {paper.title}")
        print(f"Abstract: {paper.abstract[:200]}...")
        print(f"DOI: {paper.doi}")

Patents
~~~~~~

Searching patents:

.. code-block:: python

    from binding_data_processor.web_enrichment import PatentClient

    # Configure client
    client = PatentClient(
        cache_dir="cache/patents",
    )

    # Search patents
    patents = client.search_patents(
        query="novel NMDA antagonist",
        date_range=("2020-01-01", "2024-01-01"),
    )
    for patent in patents:
        print(f"Title: {patent.title}")
        print(f"Claims: {patent.claims[:200]}...")
        print(f"Number: {patent.number}")

Chemical Databases
---------------

ChEMBL
~~~~~

Searching ChEMBL:

.. code-block:: python

    from binding_data_processor.web_enrichment import ChEMBLClient

    # Configure client
    client = ChEMBLClient()

    # Search compounds
    compounds = client.search_compounds(
        target="CHEMBL1983",  # 5-HT2A
        min_activity=1e-6,
    )
    for compound in compounds:
        print(f"ID: {compound.chembl_id}")
        print(f"SMILES: {compound.smiles}")
        print(f"Activity: {compound.activity}")

PubChem
~~~~~~

Searching PubChem:

.. code-block:: python

    from binding_data_processor.web_enrichment import PubChemClient

    # Configure client
    client = PubChemClient()

    # Search compounds
    compounds = client.search_compounds(
        query="dopamine antagonist",
        max_results=100,
    )
    for compound in compounds:
        print(f"CID: {compound.pubchem_cid}")
        print(f"SMILES: {compound.smiles}")
        print(f"Properties: {compound.properties}")

Advanced Usage
------------

Custom Enrichment
~~~~~~~~~~~~~~

Creating custom enrichment sources:

.. code-block:: python

    from binding_data_processor.web_enrichment import BaseClient

    class CustomClient(BaseClient):
        def __init__(self, api_key: str):
            super().__init__()
            self.api_key = api_key

        def get_compound_data(self, name: str) -> dict:
            # Custom data gathering logic
            data = self._fetch_data(name)
            return self._process_data(data)

    # Use custom client
    client = CustomClient(api_key="your_key")
    data = client.get_compound_data("compound_name")

Data Integration
~~~~~~~~~~~~~

Integrating data from multiple sources:

.. code-block:: python

    from binding_data_processor.web_enrichment import DataIntegrator

    # Configure integrator
    integrator = DataIntegrator(
        sources=[
            PsychonautWikiClient(),
            ErowidClient(),
            RedditMonitor(),
            PubMedClient(),
        ],
    )

    # Get integrated data
    data = integrator.get_compound_data("ketamine")
    print(f"Community Data: {data['community']}")
    print(f"Social Data: {data['social']}")
    print(f"Scientific Data: {data['scientific']}")

Error Handling
~~~~~~~~~~~

Handling web enrichment errors:

.. code-block:: python

    from binding_data_processor.web_enrichment import (
        EnrichmentManager,
        EnrichmentError,
    )

    # Configure manager
    manager = EnrichmentManager(
        retry_count=3,
        timeout=30,
    )

    # Process with error handling
    try:
        data = manager.enrich_compound(compound)
    except EnrichmentError as e:
        print(f"Enrichment failed: {e}")
        # Handle error or use fallback data
