Production Deployment
====================

This guide covers deploying ChemData in production environments.

Docker Deployment
--------------

Using Docker Compose:

.. code-block:: yaml

    # docker-compose.yml
    version: '3.8'

    services:
      # Web Application
      web:
        build: .
        image: chemdata-web
        ports:
          - "8000:8000"
        environment:
          - FLASK_ENV=production
          - FLASK_APP=binding_data_processor.web.app
          - DATABASE_URL=postgresql://user:pass@db:5432/chemdata
          - REDIS_URL=redis://cache:6379/0
        depends_on:
          - db
          - cache
          - worker

      # Background Worker
      worker:
        build: .
        image: chemdata-worker
        command: celery -A binding_data_processor.tasks worker
        environment:
          - DATABASE_URL=postgresql://user:pass@db:5432/chemdata
          - REDIS_URL=redis://cache:6379/0
        depends_on:
          - db
          - cache

      # Database
      db:
        image: postgres:13
        volumes:
          - postgres_data:/var/lib/postgresql/data
        environment:
          - POSTGRES_USER=user
          - POSTGRES_PASSWORD=pass
          - POSTGRES_DB=chemdata

      # Cache
      cache:
        image: redis:6
        volumes:
          - redis_data:/data

      # Monitoring
      prometheus:
        image: prom/prometheus
        volumes:
          - ./prometheus.yml:/etc/prometheus/prometheus.yml
        ports:
          - "9090:9090"

      grafana:
        image: grafana/grafana
        ports:
          - "3000:3000"
        depends_on:
          - prometheus

volumes:
  postgres_data:
  redis_data:

Kubernetes Deployment
-----------------

Using Kubernetes manifests:

.. code-block:: yaml

    # kubernetes/web-deployment.yaml
    apiVersion: apps/v1
    kind: Deployment
    metadata:
      name: chemdata-web
    spec:
      replicas: 3
      selector:
        matchLabels:
          app: chemdata-web
      template:
        metadata:
          labels:
            app: chemdata-web
        spec:
          containers:
          - name: web
            image: chemdata-web:latest
            ports:
            - containerPort: 8000
            env:
            - name: FLASK_ENV
              value: production
            - name: DATABASE_URL
              valueFrom:
                secretKeyRef:
                  name: chemdata-secrets
                  key: database-url
            resources:
              requests:
                memory: "512Mi"
                cpu: "250m"
              limits:
                memory: "1Gi"
                cpu: "500m"
            readinessProbe:
              httpGet:
                path: /health
                port: 8000
            livenessProbe:
              httpGet:
                path: /health
                port: 8000

Security Configuration
------------------

SSL/TLS Setup:

.. code-block:: nginx

    # nginx/chemdata.conf
    server {
        listen 443 ssl;
        server_name chemdata.example.com;

        ssl_certificate /etc/letsencrypt/live/chemdata.example.com/fullchain.pem;
        ssl_certificate_key /etc/letsencrypt/live/chemdata.example.com/privkey.pem;

        ssl_protocols TLSv1.2 TLSv1.3;
        ssl_ciphers HIGH:!aNULL:!MD5;

        location / {
            proxy_pass http://web:8000;
            proxy_set_header Host $host;
            proxy_set_header X-Real-IP $remote_addr;
        }
    }

Authentication Configuration:

.. code-block:: python

    # config/production.py
    SECURITY_CONFIG = {
        "require_auth": True,
        "auth_method": "jwt",
        "jwt_secret": os.environ["JWT_SECRET"],
        "jwt_expiry": 3600,
        "rate_limit": {
            "default": "100/hour",
            "api": "1000/hour",
        },
    }

Monitoring Setup
-------------

Prometheus configuration:

.. code-block:: yaml

    # prometheus.yml
    global:
      scrape_interval: 15s

    scrape_configs:
      - job_name: 'chemdata'
        static_configs:
          - targets: ['web:8000']

    rule_files:
      - 'alert.rules'

    alerting:
      alertmanagers:
      - static_configs:
        - targets: ['alertmanager:9093']

Application metrics:

.. code-block:: python

    from prometheus_client import Counter, Histogram

    # Request metrics
    REQUEST_COUNT = Counter(
        'request_count',
        'Number of requests',
        ['method', 'endpoint', 'status']
    )

    REQUEST_LATENCY = Histogram(
        'request_latency_seconds',
        'Request latency',
        ['method', 'endpoint']
    )

    # Processing metrics
    COMPOUNDS_PROCESSED = Counter(
        'compounds_processed_total',
        'Number of compounds processed'
    )

    PROCESSING_TIME = Histogram(
        'compound_processing_seconds',
        'Time spent processing compounds'
    )

Logging Configuration
-----------------

Using structured logging:

.. code-block:: python

    # config/logging.py
    import structlog

    structlog.configure(
        processors=[
            structlog.processors.TimeStamper(fmt="iso"),
            structlog.processors.StackInfoRenderer(),
            structlog.processors.format_exc_info,
            structlog.processors.JSONRenderer(),
        ],
        context_class=dict,
        logger_factory=structlog.PrintLoggerFactory(),
        wrapper_class=structlog.BoundLogger,
        cache_logger_on_first_use=True,
    )

    logger = structlog.get_logger()

    # Example usage
    logger.info(
        "processing_compound",
        compound_id=compound.id,
        stage="validation",
        duration=duration,
    )

Performance Tuning
---------------

Web server configuration:

.. code-block:: python

    # config/gunicorn.py
    import multiprocessing

    # Server socket
    bind = "0.0.0.0:8000"
    backlog = 2048

    # Worker processes
    workers = multiprocessing.cpu_count() * 2 + 1
    worker_class = "uvicorn.workers.UvicornWorker"
    worker_connections = 1000
    timeout = 30
    keepalive = 2

    # Process naming
    proc_name = "chemdata"

    # Logging
    accesslog = "-"
    errorlog = "-"
    loglevel = "info"

Database optimization:

.. code-block:: python

    # config/database.py
    SQLALCHEMY_DATABASE_URI = os.environ["DATABASE_URL"]
    SQLALCHEMY_POOL_SIZE = 5
    SQLALCHEMY_MAX_OVERFLOW = 10
    SQLALCHEMY_POOL_TIMEOUT = 30
    SQLALCHEMY_POOL_RECYCLE = 1800

Cache configuration:

.. code-block:: python

    # config/cache.py
    CACHE_CONFIG = {
        "CACHE_TYPE": "redis",
        "CACHE_REDIS_URL": os.environ["REDIS_URL"],
        "CACHE_DEFAULT_TIMEOUT": 300,
        "CACHE_KEY_PREFIX": "chemdata:",
    }

Backup Strategy
------------

Database backups:

.. code-block:: bash

    #!/bin/bash
    # scripts/backup.sh

    # Set variables
    BACKUP_DIR="/backups/postgres"
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    DB_CONTAINER="chemdata_db_1"

    # Create backup
    docker exec $DB_CONTAINER \
        pg_dump -U user chemdata \
        | gzip > "$BACKUP_DIR/chemdata_$TIMESTAMP.sql.gz"

    # Rotate backups (keep last 7 days)
    find $BACKUP_DIR -type f -mtime +7 -delete

Model backups:

.. code-block:: python

    # scripts/backup_models.py
    from binding_data_processor.models import ModelManager

    manager = ModelManager()

    # Backup models
    manager.backup_models(
        output_dir="/backups/models",
        include_weights=True,
        include_config=True,
        compress=True,
    )

Disaster Recovery
--------------

Recovery procedures:

.. code-block:: bash

    #!/bin/bash
    # scripts/recover.sh

    # Restore database
    BACKUP_FILE="/backups/postgres/latest.sql.gz"
    gunzip -c $BACKUP_FILE | docker exec -i chemdata_db_1 \
        psql -U user -d chemdata

    # Restore models
    python scripts/restore_models.py \
        --backup-dir /backups/models \
        --version latest

    # Verify system
    python scripts/verify_system.py \
        --check-database \
        --check-models \
        --check-cache
