Installation Guide
==================

System Requirements
-----------------

- Python 3.8 or higher
- Docker and Docker Compose (optional, for containerized deployment)
- PostgreSQL 12 or higher (optional, for persistent storage)
- Redis 6 or higher (optional, for caching)
- 8GB RAM minimum (16GB recommended)
- 20GB disk space

Python Dependencies
-----------------

Core Dependencies
~~~~~~~~~~~~~~

- RDKit: Chemical informatics library
- OpenBabel: Chemical file format conversion
- NumPy: Numerical computing
- Pandas: Data manipulation
- PyTorch: Machine learning (optional)
- Streamlit: Web interface
- Flask: API server
- SQLAlchemy: Database ORM
- Redis-py: Redis client
- Requests: HTTP client
- BeautifulSoup4: Web scraping
- PRAW: Reddit API client
- Tweepy: Twitter API client

Development Dependencies
~~~~~~~~~~~~~~~~~~~~~

- pytest: Testing framework
- mypy: Type checking
- flake8: Code linting
- black: Code formatting
- isort: Import sorting
- bandit: Security checks
- pre-commit: Git hooks
- sphinx: Documentation
- sphinx-rtd-theme: Documentation theme

Installation Methods
------------------

Quick Install
~~~~~~~~~~~

For basic usage without Docker:

.. code-block:: bash

    # Create virtual environment
    python -m venv .venv
    source .venv/bin/activate  # Linux/macOS
    # or
    .venv\\Scripts\\activate  # Windows

    # Install from PyPI
    pip install chemdata

    # Install optional ML dependencies
    pip install chemdata[ml]

Development Install
~~~~~~~~~~~~~~~~

For development and contributing:

.. code-block:: bash

    # Clone repository
    git clone https://github.com/yourusername/chemdata.git
    cd chemdata

    # Create virtual environment
    python -m venv .venv
    source .venv/bin/activate  # Linux/macOS
    # or
    .venv\\Scripts\\activate  # Windows

    # Install dependencies
    ./scripts/setup_dev.sh --dev

    # Set up pre-commit hooks
    pre-commit install

Docker Install
~~~~~~~~~~~~

For containerized deployment:

.. code-block:: bash

    # Clone repository
    git clone https://github.com/yourusername/chemdata.git
    cd chemdata

    # Build and start containers
    docker-compose up -d

Configuration
------------

Environment Variables
~~~~~~~~~~~~~~~~~~

Create a `.env` file with your configuration:

.. code-block:: bash

    # Data directories
    CHEMDATA_DATA_DIR=./data
    CHEMDATA_CACHE_DIR=./cache
    CHEMDATA_LOG_DIR=./logs
    CHEMDATA_OUTPUT_DIR=./output
    CHEMDATA_MODEL_DIR=./models

    # API credentials
    REDDIT_CLIENT_ID=your_client_id
    REDDIT_CLIENT_SECRET=your_client_secret
    TWITTER_API_KEY=your_api_key
    TWITTER_API_SECRET=your_api_secret

    # Database
    POSTGRES_USER=chemdata
    POSTGRES_PASSWORD=chemdata
    POSTGRES_DB=chemdata
    POSTGRES_HOST=postgres

    # Redis
    REDIS_HOST=redis
    REDIS_PORT=6379

    # Web server
    FLASK_APP=binding_data_processor.web.app
    FLASK_ENV=development
    FLASK_DEBUG=1

API Credentials
~~~~~~~~~~~~~

1. Reddit API:
   - Create an app at https://www.reddit.com/prefs/apps
   - Set type to "script"
   - Note the client ID and secret

2. Twitter API:
   - Apply for access at https://developer.twitter.com
   - Create a project and app
   - Note the API key and secret

3. Other APIs:
   - PsychonautWiki: No authentication required
   - Erowid: Contact for API access
   - TripSit: No authentication required
   - ChEMBL: No authentication required
   - PubChem: No authentication required

Verification
-----------

Test Installation
~~~~~~~~~~~~~~~

Run tests to verify installation:

.. code-block:: bash

    # Run all tests
    pytest

    # Run specific test suite
    pytest tests/test_pipeline.py

    # Run with coverage
    pytest --cov=binding_data_processor

Run Example
~~~~~~~~~

Process some example compounds:

.. code-block:: bash

    # Process example data
    python -m binding_data_processor.cli process-compounds \
        --input examples/data/example_compounds.tsv \
        --output results/ \
        --enable-ml \
        --enable-web

Start Web Interface
~~~~~~~~~~~~~~~~

Launch the web interface:

.. code-block:: bash

    # Start web app
    streamlit run examples/web_app/app.py

Troubleshooting
-------------

Common Issues
~~~~~~~~~~~

1. RDKit Installation:
   - Use conda: ``conda install -c conda-forge rdkit``
   - Or build from source following RDKit docs

2. OpenBabel Installation:
   - Linux: ``apt-get install openbabel``
   - macOS: ``brew install open-babel``
   - Windows: Download installer from website

3. Database Connection:
   - Check PostgreSQL is running
   - Verify credentials in .env
   - Check network connectivity

4. Redis Connection:
   - Check Redis is running
   - Verify host/port in .env
   - Check network connectivity

Getting Help
~~~~~~~~~~

- Check the :doc:`troubleshooting` guide
- Search GitHub issues
- Join our Discord server
- Contact maintainers
