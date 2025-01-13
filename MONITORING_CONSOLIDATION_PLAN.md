# Monitoring Consolidation Plan

## Current State

We have multiple monitoring implementations:

### 1. Pipeline Monitoring (MonitoringManager)
- Location: binding_data_processor/pipeline/infrastructure/monitoring.py
- Features:
  - Comprehensive metrics
  - Log management
  - Report generation
  - Alert system
  - Performance tracking
  - File-based storage

### 2. Document Monitoring (DirectoryMonitor)
- Location: binding_data_processor/processors/document/monitor.py
- Features:
  - File system events
  - Directory watching
  - Async processing
  - Callback system
  - Basic error tracking

## Consolidation Strategy

### 1. Create Base Monitor Interface

```python
from abc import ABC, abstractmethod
from typing import Any, Dict, Optional
from pathlib import Path

class BaseMonitor(ABC):
    """Base class for all monitoring implementations."""
    
    @abstractmethod
    async def start(self) -> None:
        """Start monitoring."""
        pass
        
    @abstractmethod
    async def stop(self) -> None:
        """Stop monitoring."""
        pass
        
    @abstractmethod
    async def get_stats(self) -> Dict[str, Any]:
        """Get monitoring statistics."""
        pass
        
    @abstractmethod
    async def record_error(self, error: Exception) -> None:
        """Record an error occurrence."""
        pass
```

### 2. Create Monitoring Components

1. Metrics Component:
```python
@dataclass
class MonitoringMetrics:
    """Monitoring metrics."""
    total_operations: int = 0
    completed_operations: int = 0
    failed_operations: int = 0
    error_count: int = 0
    start_time: Optional[datetime] = None
    latencies: List[float] = field(default_factory=list)
```

2. Logging Component:
```python
class LoggingComponent:
    """Logging management component."""
    
    def __init__(
        self,
        log_dir: Path,
        retention_days: int = 30,
        log_format: str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    ):
        self.log_dir = log_dir
        self.retention_days = retention_days
        self._setup_logging(log_format)
```

3. Alert Component:
```python
class AlertComponent:
    """Alert management component."""
    
    def __init__(self, thresholds: Dict[str, float]):
        self.thresholds = thresholds
        self.handlers: Dict[str, Callable] = {}
        
    async def check_threshold(self, metric: str, value: float) -> None:
        """Check if metric exceeds threshold."""
        if metric in self.thresholds and value > self.thresholds[metric]:
            await self._trigger_alert(metric, value)
```

4. Event Component:
```python
class EventComponent:
    """Event handling component."""
    
    def __init__(self):
        self.handlers: Dict[str, List[Callable]] = {}
        self._lock = asyncio.Lock()
        
    async def emit(self, event: str, data: Any) -> None:
        """Emit an event."""
        async with self._lock:
            for handler in self.handlers.get(event, []):
                try:
                    await handler(data)
                except Exception as e:
                    logger.error(f"Event handler error: {str(e)}")
```

### 3. Implement Concrete Monitors

1. Pipeline Monitor:
```python
class PipelineMonitor(BaseMonitor):
    """Pipeline monitoring implementation."""
    
    def __init__(
        self,
        config: MonitoringConfig,
        metrics: Optional[MonitoringMetrics] = None,
        logger: Optional[logging.Logger] = None
    ):
        self.config = config
        self.metrics = metrics or MonitoringMetrics()
        self.logger = logger or logging.getLogger(__name__)
        self.logging = LoggingComponent(config.log_dir)
        self.alerts = AlertComponent(config.alert_thresholds)
        self.events = EventComponent()
```

2. File System Monitor:
```python
class FileSystemMonitor(BaseMonitor):
    """File system monitoring implementation."""
    
    def __init__(
        self,
        watch_dirs: Dict[Path, Set[str]],
        processor: Optional[Callable] = None,
        logger: Optional[logging.Logger] = None
    ):
        self.watch_dirs = watch_dirs
        self.processor = processor
        self.logger = logger or logging.getLogger(__name__)
        self.events = EventComponent()
        self.metrics = MonitoringMetrics()
```

3. Web Monitor:
```python
class WebMonitor(BaseMonitor):
    """Web application monitoring implementation."""
    
    def __init__(
        self,
        metrics: Optional[MonitoringMetrics] = None,
        logger: Optional[logging.Logger] = None
    ):
        self.metrics = metrics or MonitoringMetrics()
        self.logger = logger or logging.getLogger(__name__)
        self.events = EventComponent()
```

### 4. Create Monitor Factory

```python
class MonitorFactory:
    """Factory for creating monitor instances."""
    
    @staticmethod
    def create(
        monitor_type: str,
        config: Optional[Dict[str, Any]] = None
    ) -> BaseMonitor:
        """Create monitor instance."""
        if monitor_type == "pipeline":
            return PipelineMonitor(config or {})
        elif monitor_type == "filesystem":
            return FileSystemMonitor(config or {})
        elif monitor_type == "web":
            return WebMonitor(config or {})
        raise ValueError(f"Unknown monitor type: {monitor_type}")
```

### 5. Migration Steps

1. Create New Structure:
```
binding_data_processor/infrastructure/monitoring/
  ├── __init__.py
  ├── base.py          # Base classes
  ├── components/      # Monitoring components
  │   ├── __init__.py
  │   ├── metrics.py
  │   ├── logging.py
  │   ├── alerts.py
  │   └── events.py
  ├── monitors/        # Monitor implementations
  │   ├── __init__.py
  │   ├── pipeline.py
  │   ├── filesystem.py
  │   └── web.py
  ├── factory.py       # Monitor factory
  └── config.py        # Configuration
```

2. Update Dependencies:
- Update pipeline code to use new monitors
- Update document processor to use new monitors
- Add tests for new implementation
- Update documentation

3. Deprecate Old Implementations:
- Mark old classes as deprecated
- Provide migration guide
- Remove after transition period

### 6. Testing Strategy

1. Unit Tests:
```python
class TestBaseMonitor:
    """Test monitor base implementation."""
    
    async def test_metrics(self, monitor: BaseMonitor):
        """Test metric collection."""
        await monitor.start()
        stats = await monitor.get_stats()
        assert "total_operations" in stats
        
    async def test_error_handling(self, monitor: BaseMonitor):
        """Test error handling."""
        error = Exception("Test error")
        await monitor.record_error(error)
        stats = await monitor.get_stats()
        assert stats["error_count"] > 0
```

2. Integration Tests:
```python
class TestMonitorIntegration:
    """Test monitor integration."""
    
    async def test_pipeline_integration(self, monitor: PipelineMonitor):
        """Test pipeline monitoring."""
        # Test pipeline specific features
        
    async def test_filesystem_integration(self, monitor: FileSystemMonitor):
        """Test filesystem monitoring."""
        # Test filesystem specific features
```

### 7. Documentation

1. API Documentation:
```python
class Monitor:
    """Unified monitoring interface.
    
    This class provides a consistent interface for monitoring across
    the application. It supports metrics collection, logging, alerts,
    and event handling.
    
    Examples:
        >>> monitor = Monitor(monitor_type="pipeline")
        >>> await monitor.start()
        >>> stats = await monitor.get_stats()
    """
```

2. Migration Guide:
- Document changes from old implementations
- Provide examples of updating code
- List breaking changes
- Include troubleshooting tips

### 8. Benefits

1. Consistency:
- Single monitoring interface
- Unified metrics collection
- Standard event handling
- Consistent error tracking

2. Flexibility:
- Modular components
- Pluggable architecture
- Configurable features
- Easy to extend

3. Reliability:
- Better error handling
- Comprehensive testing
- Clear documentation
- Migration support

4. Performance:
- Optimized implementations
- Proper concurrency
- Efficient event handling
- Better resource usage

### 9. Next Steps

1. Implementation:
- Create new package structure
- Implement base classes
- Add monitoring components
- Write tests

2. Migration:
- Update pipeline code
- Update filesystem code
- Run integration tests
- Update docs

3. Deployment:
- Stage rollout
- Monitor performance
- Collect feedback
- Make adjustments

4. Cleanup:
- Remove old code
- Update dependencies
- Final testing
- Release notes
