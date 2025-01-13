"""Monitoring module for structure analysis system.

This module provides:
1. System metrics
2. Health monitoring
3. Performance tracking
4. Monitoring utilities
"""

import logging
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Callable
from pathlib import Path
import json
import threading
from datetime import datetime

from .logging import get_logger
from .errors import ValidationError

logger = get_logger(__name__)


@dataclass
class Metric:
    """Metric data point."""

    name: str
    value: float
    timestamp: float = field(default_factory=time.time)
    labels: Optional[Dict[str, str]] = None


@dataclass
class HealthStatus:
    """Component health status."""

    component: str
    status: str  # "healthy", "degraded", "failed"
    message: Optional[str] = None
    timestamp: float = field(default_factory=time.time)
    details: Optional[Dict[str, Any]] = None


@dataclass
class PerformanceMetrics:
    """Performance metrics collection."""

    operation: str
    start_time: float
    end_time: Optional[float] = None
    metrics: Dict[str, float] = field(default_factory=dict)
    labels: Optional[Dict[str, str]] = None

    @property
    def duration(self) -> Optional[float]:
        """Get operation duration in seconds."""
        if self.end_time is None:
            return None
        return self.end_time - self.start_time


class MetricsCollector:
    """Collector for system metrics."""

    def __init__(self, persist: bool = False, metrics_dir: Optional[Path] = None):
        """Initialize metrics collector.

        Args:
            persist: Whether to persist metrics to disk
            metrics_dir: Directory for persistent metrics
        """
        self.persist = persist
        self.metrics_dir = Path(metrics_dir) if metrics_dir else None
        self._metrics: List[Metric] = []
        self._lock = threading.Lock()

        if self.persist:
            if not self.metrics_dir:
                raise ValueError("Metrics directory required for persistence")
            self.metrics_dir.mkdir(parents=True, exist_ok=True)

    def record(
        self,
        name: str,
        value: float,
        labels: Optional[Dict[str, str]] = None,
    ):
        """Record metric.

        Args:
            name: Metric name
            value: Metric value
            labels: Optional metric labels
        """
        metric = Metric(
            name=name,
            value=value,
            timestamp=time.time(),
            labels=labels,
        )

        with self._lock:
            self._metrics.append(metric)
            if self.persist:
                self._persist_metric(metric)

    def get_metrics(
        self,
        name: Optional[str] = None,
        labels: Optional[Dict[str, str]] = None,
    ) -> List[Metric]:
        """Get recorded metrics.

        Args:
            name: Optional metric name filter
            labels: Optional label filters

        Returns:
            List of matching metrics
        """
        with self._lock:
            metrics = self._metrics.copy()

        if name:
            metrics = [m for m in metrics if m.name == name]

        if labels:
            metrics = [m for m in metrics if m.labels and all(m.labels.get(k) == v for k, v in labels.items())]

        return metrics

    def clear(self):
        """Clear recorded metrics."""
        with self._lock:
            self._metrics.clear()
            if self.persist:
                for path in self.metrics_dir.glob("*.metric"):
                    path.unlink()

    def _persist_metric(self, metric: Metric):
        """Persist metric to disk.

        Args:
            metric: Metric to persist
        """
        try:
            timestamp = datetime.fromtimestamp(metric.timestamp)
            filename = f"{metric.name}_{timestamp.strftime('%Y%m%d_%H%M%S')}.metric"
            path = self.metrics_dir / filename

            data = {
                "name": metric.name,
                "value": metric.value,
                "timestamp": metric.timestamp,
                "labels": metric.labels,
            }

            with open(path, "w") as f:
                json.dump(data, f)

        except Exception as e:
            logger.error(f"Error persisting metric: {str(e)}")


class HealthMonitor:
    """Monitor for component health."""

    def __init__(self):
        """Initialize health monitor."""
        self._status: Dict[str, HealthStatus] = {}
        self._lock = threading.Lock()

    def update(
        self,
        component: str,
        status: str,
        message: Optional[str] = None,
        details: Optional[Dict[str, Any]] = None,
    ):
        """Update component health status.

        Args:
            component: Component name
            status: Health status
            message: Optional status message
            details: Optional status details
        """
        if status not in ["healthy", "degraded", "failed"]:
            raise ValueError(f"Invalid status: {status}")

        health = HealthStatus(
            component=component,
            status=status,
            message=message,
            details=details,
        )

        with self._lock:
            self._status[component] = health

    def get_status(
        self,
        component: Optional[str] = None,
    ) -> Dict[str, HealthStatus]:
        """Get component health status.

        Args:
            component: Optional component name

        Returns:
            Dict mapping component names to health status
        """
        with self._lock:
            if component:
                status = self._status.get(component)
                return {component: status} if status else {}
            return self._status.copy()

    def is_healthy(self, component: str) -> bool:
        """Check if component is healthy.

        Args:
            component: Component name

        Returns:
            True if component is healthy
        """
        with self._lock:
            status = self._status.get(component)
            return status is not None and status.status == "healthy"


class PerformanceMonitor:
    """Monitor for performance metrics."""

    def __init__(self):
        """Initialize performance monitor."""
        self._metrics: Dict[str, List[PerformanceMetrics]] = {}
        self._lock = threading.Lock()

    def start(
        self,
        operation: str,
        labels: Optional[Dict[str, str]] = None,
    ) -> PerformanceMetrics:
        """Start monitoring operation.

        Args:
            operation: Operation name
            labels: Optional operation labels

        Returns:
            Performance metrics instance
        """
        metrics = PerformanceMetrics(
            operation=operation,
            start_time=time.time(),
            labels=labels,
        )

        with self._lock:
            if operation not in self._metrics:
                self._metrics[operation] = []
            self._metrics[operation].append(metrics)

        return metrics

    def stop(
        self,
        metrics: PerformanceMetrics,
        additional_metrics: Optional[Dict[str, float]] = None,
    ):
        """Stop monitoring operation.

        Args:
            metrics: Performance metrics instance
            additional_metrics: Optional additional metrics
        """
        metrics.end_time = time.time()
        if additional_metrics:
            metrics.metrics.update(additional_metrics)

    def get_metrics(
        self,
        operation: Optional[str] = None,
        labels: Optional[Dict[str, str]] = None,
    ) -> List[PerformanceMetrics]:
        """Get performance metrics.

        Args:
            operation: Optional operation name filter
            labels: Optional label filters

        Returns:
            List of matching metrics
        """
        with self._lock:
            if operation:
                metrics = self._metrics.get(operation, []).copy()
            else:
                metrics = [m for metrics in self._metrics.values() for m in metrics]

        if labels:
            metrics = [m for m in metrics if m.labels and all(m.labels.get(k) == v for k, v in labels.items())]

        return metrics


def create_metrics_collector(
    metrics_dir: Optional[Path] = None,
) -> MetricsCollector:
    """Create metrics collector.

    Args:
        metrics_dir: Optional directory for persistent metrics

    Returns:
        MetricsCollector instance
    """
    return MetricsCollector(
        persist=True if metrics_dir else False,
        metrics_dir=metrics_dir,
    )


def monitor_performance(monitor: PerformanceMonitor, operation: str):
    """Decorator for monitoring function performance.

    Args:
        monitor: Performance monitor instance
        operation: Operation name

    Returns:
        Decorated function
    """

    def decorator(func: Callable):
        def wrapper(*args, **kwargs):
            metrics = monitor.start(operation)
            try:
                result = func(*args, **kwargs)
                return result
            finally:
                monitor.stop(metrics)

        return wrapper

    return decorator
