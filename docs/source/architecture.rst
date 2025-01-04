System Architecture
==================

This guide covers the system architecture of ChemData, explaining how components work together.

Overview
-------

High-level architecture:

.. code-block:: text

    ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
    │   Web UI    │     │    API      │     │  Workers    │
    │  (React)    │────▶│  (Flask)    │────▶│  (Celery)   │
    └─────────────┘     └─────────────┘     └─────────────┘
           │                   │                    │
           │                   │                    │
           ▼                   ▼                    ▼
    ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
    │   Cache     │     │  Database   │     │   Storage   │
    │   (Redis)   │     │ (Postgres)  │     │    (S3)     │
    └─────────────┘     └─────────────┘     └─────────────┘

Component Roles
------------

Web Interface:

.. code-block:: python

    # web/app.py
    class ChemDataWeb:
        """Web interface component."""

        def __init__(self):
            self.api_client = APIClient()
            self.cache = CacheClient()
            self.components = {
                "list": CompoundListView(),
                "detail": CompoundDetailView(),
                "search": SearchView(),
            }

        def render_page(self, page_type: str, **kwargs):
            """Render page with data."""
            # Get data from API
            data = self.api_client.get_data(
                endpoint=page_type,
                params=kwargs,
            )

            # Check cache
            cached = self.cache.get(cache_key)
            if cached:
                return cached

            # Render component
            component = self.components[page_type]
            html = component.render(data)

            # Cache result
            self.cache.set(cache_key, html)

            return html

API Layer:

.. code-block:: python

    # api/server.py
    class ChemDataAPI:
        """API server component."""

        def __init__(self):
            self.db = DatabaseClient()
            self.cache = CacheClient()
            self.worker = WorkerClient()
            self.auth = AuthManager()

        def handle_request(self, request: Request) -> Response:
            """Handle API request."""
            # Authenticate
            user = self.auth.authenticate(request)
            if not user:
                return Response(401)

            # Rate limit
            if not self.rate_limiter.allow(user):
                return Response(429)

            # Process request
            try:
                # Check cache
                cached = self.cache.get(cache_key)
                if cached:
                    return Response(200, cached)

                # Handle request
                if request.is_async:
                    # Queue task
                    task_id = self.worker.queue_task(
                        task_type=request.type,
                        params=request.params,
                    )
                    return Response(202, {"task_id": task_id})
                else:
                    # Process synchronously
                    result = self.process_request(request)
                    return Response(200, result)

            except Exception as e:
                return Response(500, str(e))

Worker System:

.. code-block:: python

    # workers/processor.py
    class ChemDataWorker:
        """Background worker component."""

        def __init__(self):
            self.db = DatabaseClient()
            self.storage = StorageClient()
            self.processors = {
                "compounds": CompoundProcessor(),
                "analysis": AnalysisProcessor(),
                "ml": MLProcessor(),
            }

        def process_task(self, task: Task) -> Result:
            """Process background task."""
            try:
                # Get processor
                processor = self.processors[task.type]

                # Process task
                result = processor.process(task.params)

                # Store result
                self.storage.store(
                    task_id=task.id,
                    result=result,
                )

                # Update status
                self.db.update_task_status(
                    task_id=task.id,
                    status="completed",
                    result=result,
                )

                return result

            except Exception as e:
                # Handle failure
                self.db.update_task_status(
                    task_id=task.id,
                    status="failed",
                    error=str(e),
                )
                raise

Data Flow
--------

Request handling:

.. code-block:: text

    1. Client Request
       ├─▶ Load Balancer
       │   ├─▶ Web Server 1
       │   ├─▶ Web Server 2
       │   └─▶ Web Server 3
       │
    2. Authentication
       ├─▶ JWT Validation
       ├─▶ Rate Limiting
       └─▶ Permission Check
       
    3. Request Processing
       ├─▶ Input Validation
       ├─▶ Cache Check
       └─▶ Task Creation
       
    4. Background Processing
       ├─▶ Worker Pool
       │   ├─▶ Worker 1
       │   ├─▶ Worker 2
       │   └─▶ Worker 3
       │
    5. Result Handling
       ├─▶ Result Storage
       ├─▶ Cache Update
       └─▶ Client Notification

Data Storage:

.. code-block:: text

    ┌─────────────┐
    │  Raw Data   │
    │  (Files)    │
    └─────────────┘
          │
          ▼
    ┌─────────────┐
    │  Database   │
    │ (Processed) │
    └─────────────┘
          │
          ▼
    ┌─────────────┐
    │   Cache     │
    │ (Computed)  │
    └─────────────┘

Scalability
---------

Horizontal scaling:

.. code-block:: python

    # config/scaling.py
    SCALING_CONFIG = {
        "web": {
            "min_instances": 3,
            "max_instances": 10,
            "scale_up_threshold": 0.75,  # CPU usage
            "scale_down_threshold": 0.25,
        },
        "worker": {
            "min_instances": 5,
            "max_instances": 20,
            "scale_up_threshold": 100,  # Queue size
            "scale_down_threshold": 10,
        },
        "cache": {
            "shards": 3,
            "replicas": 2,
        },
        "database": {
            "read_replicas": 3,
            "write_nodes": 2,
        },
    }

Load balancing:

.. code-block:: python

    # infrastructure/loadbalancer.py
    class LoadBalancer:
        """Load balancer configuration."""

        def __init__(self):
            self.algorithm = "least_connections"
            self.health_check = {
                "path": "/health",
                "interval": 30,
                "timeout": 5,
                "threshold": 3,
            }
            self.ssl = {
                "enabled": True,
                "redirect_http": True,
            }
            self.session_persistence = {
                "enabled": True,
                "cookie_name": "CHEMDATA_SERVER",
            }

Service Discovery:

.. code-block:: python

    # infrastructure/discovery.py
    class ServiceRegistry:
        """Service discovery component."""

        def __init__(self):
            self.consul = ConsulClient()
            self.services = {}

        def register_service(
            self,
            name: str,
            host: str,
            port: int,
            tags: List[str],
        ):
            """Register service with Consul."""
            service_id = f"{name}-{host}-{port}"
            
            self.consul.register(
                service_id=service_id,
                name=name,
                host=host,
                port=port,
                tags=tags,
                health_check={
                    "http": f"http://{host}:{port}/health",
                    "interval": "30s",
                },
            )
            
            self.services[service_id] = {
                "name": name,
                "host": host,
                "port": port,
                "tags": tags,
            }

Security
-------

Authentication flow:

.. code-block:: python

    # security/auth.py
    class AuthManager:
        """Authentication manager."""

        def __init__(self):
            self.jwt = JWTManager()
            self.oauth = OAuthManager()
            self.mfa = MFAManager()

        def authenticate(self, request: Request) -> Optional[User]:
            """Authenticate request."""
            # Check JWT
            token = request.headers.get("Authorization")
            if token:
                return self.jwt.validate_token(token)

            # Check session
            session = request.cookies.get("session")
            if session:
                return self.validate_session(session)

            return None

        def authorize(self, user: User, resource: str) -> bool:
            """Check user permissions."""
            return self.rbac.check_permission(
                user=user,
                resource=resource,
                action="read",
            )

Monitoring
--------

Health checks:

.. code-block:: python

    # monitoring/health.py
    class HealthChecker:
        """System health monitoring."""

        def __init__(self):
            self.checks = {
                "database": self.check_database,
                "cache": self.check_cache,
                "workers": self.check_workers,
                "storage": self.check_storage,
            }

        def check_health(self) -> Dict[str, bool]:
            """Run all health checks."""
            results = {}
            for name, check in self.checks.items():
                try:
                    check()
                    results[name] = True
                except Exception as e:
                    results[name] = False
                    self.alert(f"Health check failed: {name}", e)
            return results

Metrics collection:

.. code-block:: python

    # monitoring/metrics.py
    class MetricsCollector:
        """System metrics collection."""

        def __init__(self):
            self.prometheus = PrometheusClient()
            self.metrics = {
                "requests": Counter(
                    "request_count",
                    "Number of requests",
                    ["method", "endpoint"],
                ),
                "latency": Histogram(
                    "request_latency",
                    "Request latency in seconds",
                    ["method", "endpoint"],
                ),
                "errors": Counter(
                    "error_count",
                    "Number of errors",
                    ["type", "code"],
                ),
                "processing_time": Histogram(
                    "processing_time",
                    "Task processing time in seconds",
                    ["task_type"],
                ),
            }

        def record_metric(
            self,
            name: str,
            value: float,
            labels: Dict[str, str],
        ):
            """Record metric value."""
            metric = self.metrics[name]
            metric.labels(**labels).observe(value)

Deployment
--------

Service configuration:

.. code-block:: yaml

    # kubernetes/service.yaml
    apiVersion: v1
    kind: Service
    metadata:
      name: chemdata-web
      annotations:
        prometheus.io/scrape: "true"
        prometheus.io/port: "8000"
    spec:
      selector:
        app: chemdata-web
      ports:
        - port: 80
          targetPort: 8000
      type: LoadBalancer

Resource management:

.. code-block:: yaml

    # kubernetes/resources.yaml
    apiVersion: v1
    kind: ResourceQuota
    metadata:
      name: chemdata-quota
    spec:
      hard:
        requests.cpu: "4"
        requests.memory: 8Gi
        limits.cpu: "8"
        limits.memory: 16Gi
        persistentvolumeclaims: "10"
