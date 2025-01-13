"""Statistics tracking for infrastructure components.

This module provides metrics collection and statistics tracking for:
1. Cache performance
2. Rate limit usage
3. Circuit breaker states
4. Monitoring metrics

These statistics help monitor system health and performance.
"""

from collections import defaultdict, deque
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from enum import Enum
from typing import Any, Deque, Dict, List, Optional, Set, Union

import time


class MetricType(Enum):
    """Types of metrics that can be tracked."""

    COUNTER = "counter"  # Monotonically increasing value
    GAUGE = "gauge"  # Value that can go up and down
    HISTOGRAM = "histogram"  # Distribution of values
    SUMMARY = "summary"  # Statistical summary of values


@dataclass
class TimeWindow:
    """Time window for tracking metrics."""

    window_size: int  # Window size in seconds
    num_buckets: int  # Number of buckets in window
    bucket_size: int = field(init=False)  # Size of each bucket in seconds

    def __post_init__(self):
        """Initialize derived values."""
        self.bucket_size = self.window_size // self.num_buckets
        if self.window_size % self.num_buckets != 0:
            raise ValueError(f"Window size {self.window_size} must be divisible by " f"number of buckets {self.num_buckets}")


@dataclass
class MetricValue:
    """Value for a single metric."""

    type: MetricType
    value: Union[int, float, Dict[str, float]]
    timestamp: float = field(default_factory=time.time)
    tags: Dict[str, str] = field(default_factory=dict)


class RollingMetric:
    """Metric that maintains values over a rolling time window."""

    def __init__(self, window: TimeWindow, metric_type: MetricType = MetricType.COUNTER):
        """Initialize rolling metric.

        Args:
            window: Time window configuration
            metric_type: Type of metric to track
        """
        self.window = window
        self.type = metric_type
        self._buckets: List[Dict[str, float]] = [{} for _ in range(window.num_buckets)]
        self._current_bucket = 0
        self._last_rotation = time.time()

    def add(self, value: float, timestamp: Optional[float] = None) -> None:
        """Add value to current bucket.

        Args:
            value: Value to add
            timestamp: Optional timestamp (defaults to current time)
        """
        self._rotate_if_needed(timestamp or time.time())
        if self.type == MetricType.COUNTER:
            self._buckets[self._current_bucket]["sum"] = self._buckets[self._current_bucket].get("sum", 0) + value
        elif self.type == MetricType.GAUGE:
            self._buckets[self._current_bucket]["last"] = value
        elif self.type in {MetricType.HISTOGRAM, MetricType.SUMMARY}:
            values = self._buckets[self._current_bucket].setdefault("values", [])
            values.append(value)

    def get_stats(self) -> Dict[str, float]:
        """Get statistics for current window.

        Returns:
            Dictionary of statistics
        """
        self._rotate_if_needed(time.time())

        if self.type == MetricType.COUNTER:
            return {"sum": sum(b.get("sum", 0) for b in self._buckets)}

        elif self.type == MetricType.GAUGE:
            last_value = None
            for b in reversed(self._buckets):
                if "last" in b:
                    last_value = b["last"]
                    break
            return {"value": last_value} if last_value is not None else {}

        elif self.type in {MetricType.HISTOGRAM, MetricType.SUMMARY}:
            values = []
            for b in self._buckets:
                values.extend(b.get("values", []))

            if not values:
                return {}

            values.sort()
            n = len(values)

            stats = {
                "count": n,
                "min": values[0],
                "max": values[-1],
                "mean": sum(values) / n,
                "p50": values[n // 2],
                "p90": values[int(n * 0.9)],
                "p95": values[int(n * 0.95)],
                "p99": values[int(n * 0.99)],
            }

            return stats

    def _rotate_if_needed(self, timestamp: float) -> None:
        """Rotate buckets if needed based on timestamp.

        Args:
            timestamp: Current timestamp
        """
        elapsed = timestamp - self._last_rotation
        if elapsed < self.window.bucket_size:
            return

        num_rotations = int(elapsed / self.window.bucket_size)
        if num_rotations >= self.window.num_buckets:
            # Clear all buckets if more than a full window has elapsed
            self._buckets = [{} for _ in range(self.window.num_buckets)]
            self._current_bucket = 0
        else:
            # Rotate buckets and clear old ones
            for _ in range(num_rotations):
                self._current_bucket = (self._current_bucket + 1) % self.window.num_buckets
                self._buckets[self._current_bucket] = {}

        self._last_rotation = timestamp


class StatsTracker:
    """Tracks statistics for a component."""

    def __init__(self, window_size: int = 60, num_buckets: int = 6):
        """Initialize stats tracker.

        Args:
            window_size: Window size in seconds
            num_buckets: Number of buckets in window
        """
        self.window = TimeWindow(window_size, num_buckets)
        self._metrics: Dict[str, RollingMetric] = {}
        self._tags: Dict[str, str] = {}

    def add_metric(self, name: str, type: MetricType, tags: Optional[Dict[str, str]] = None) -> None:
        """Add a new metric to track.

        Args:
            name: Metric name
            type: Metric type
            tags: Optional metric tags
        """
        if name in self._metrics:
            raise ValueError(f"Metric {name} already exists")

        self._metrics[name] = RollingMetric(self.window, type)
        if tags:
            self._tags.update(tags)

    def record(self, name: str, value: float, timestamp: Optional[float] = None) -> None:
        """Record a metric value.

        Args:
            name: Metric name
            value: Value to record
            timestamp: Optional timestamp
        """
        if name not in self._metrics:
            raise ValueError(f"Unknown metric: {name}")
        self._metrics[name].add(value, timestamp)

    def get_stats(self, names: Optional[List[str]] = None) -> Dict[str, Dict[str, float]]:
        """Get statistics for metrics.

        Args:
            names: Optional list of metric names to get stats for

        Returns:
            Dictionary mapping metric names to their statistics
        """
        if names is None:
            names = list(self._metrics.keys())

        return {name: self._metrics[name].get_stats() for name in names if name in self._metrics}

    def get_tags(self) -> Dict[str, str]:
        """Get metric tags.

        Returns:
            Dictionary of metric tags
        """
        return dict(self._tags)


class StatsManager:
    """Manages statistics for multiple components."""

    def __init__(self):
        """Initialize stats manager."""
        self._trackers: Dict[str, StatsTracker] = {}

    def get_tracker(self, name: str, window_size: int = 60, num_buckets: int = 6) -> StatsTracker:
        """Get or create stats tracker for component.

        Args:
            name: Component name
            window_size: Window size in seconds
            num_buckets: Number of buckets in window

        Returns:
            Stats tracker for component
        """
        if name not in self._trackers:
            self._trackers[name] = StatsTracker(window_size, num_buckets)
        return self._trackers[name]

    def get_all_stats(self) -> Dict[str, Dict[str, Dict[str, float]]]:
        """Get statistics for all components.

        Returns:
            Dictionary mapping component names to their statistics
        """
        return {name: tracker.get_stats() for name, tracker in self._trackers.items()}
