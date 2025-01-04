"""Pipeline monitoring.

This module provides the MonitoringManager class that:
1. Tracks pipeline progress
2. Monitors performance
3. Collects metrics
4. Generates reports
5. Handles alerts
"""

import logging
import threading
from pathlib import Path
from typing import Dict, Optional, Any, List
from dataclasses import dataclass, field
from datetime import datetime, timedelta
import json


@dataclass
class MonitoringConfig:
    """Monitoring configuration."""
    
    # Progress settings
    monitor_progress: bool = True
    progress_interval: int = 5  # seconds
    
    # Performance settings
    monitor_performance: bool = True
    performance_interval: int = 60  # seconds
    
    # Metric settings
    collect_metrics: bool = True
    metric_interval: int = 300  # seconds
    
    # Alert settings
    enable_alerts: bool = True
    alert_thresholds: Dict[str, float] = field(default_factory=lambda: {
        "error_rate": 0.1,  # 10%
        "failure_rate": 0.05,  # 5%
        "latency": 1000.0,  # ms
    })
    
    # Report settings
    generate_reports: bool = True
    report_dir: Optional[Path] = None
    report_format: str = "json"


@dataclass
class MonitoringStats:
    """Monitoring statistics."""
    
    # Progress stats
    total_operations: int = 0
    completed_operations: int = 0
    failed_operations: int = 0
    start_time: Optional[datetime] = None
    
    # Performance stats
    operation_times: List[float] = field(default_factory=list)
    peak_latency: float = 0.0
    avg_latency: float = 0.0
    
    # Metric stats
    total_errors: int = 0
    error_types: Dict[str, int] = field(default_factory=dict)
    component_stats: Dict[str, Dict] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary format."""
        return {
            "progress": {
                "total": self.total_operations,
                "completed": self.completed_operations,
                "failed": self.failed_operations,
                "success_rate": self._get_success_rate(),
                "duration": self._get_duration(),
            },
            "performance": {
                "peak_latency": self.peak_latency,
                "avg_latency": self.avg_latency,
                "operation_count": len(self.operation_times),
            },
            "metrics": {
                "errors": {
                    "total": self.total_errors,
                    "types": self.error_types,
                },
                "components": self.component_stats,
            },
        }
    
    def _get_success_rate(self) -> Optional[float]:
        """Get operation success rate."""
        if not self.total_operations:
            return None
        return (
            self.completed_operations / self.total_operations
            if self.total_operations > 0
            else 0.0
        )
    
    def _get_duration(self) -> Optional[float]:
        """Get total duration in seconds."""
        if not self.start_time:
            return None
        return (datetime.now() - self.start_time).total_seconds()


class MonitoringManager:
    """Manager for pipeline monitoring."""

    def __init__(
        self,
        monitor_resources: bool = True,
        monitor_progress: bool = True,
        log_level: str = "INFO",
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize monitoring manager.
        
        Args:
            monitor_resources: Whether to monitor resources
            monitor_progress: Whether to monitor progress
            log_level: Logging level
            logger: Optional logger instance
        """
        self.config = MonitoringConfig(
            monitor_progress=monitor_progress,
        )
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        
        # Initialize stats
        self.stats = MonitoringStats()
        
        # Initialize state
        self._lock = threading.Lock()
        self._active = False
        self._monitor_thread = None
        
        # Initialize report directory
        self._init_report_dir()

    def _init_report_dir(self) -> None:
        """Initialize report directory."""
        try:
            if self.config.report_dir:
                self.config.report_dir.mkdir(parents=True, exist_ok=True)
                
        except Exception as e:
            self.logger.error(f"Failed to initialize report directory: {str(e)}")
            raise

    def start(self) -> None:
        """Start monitoring manager."""
        with self._lock:
            if self._active:
                return
            
            self._active = True
            self.stats.start_time = datetime.now()
            
            # Start monitoring thread
            if self.config.monitor_progress:
                self._monitor_thread = threading.Thread(
                    target=self._monitor_pipeline,
                    name="pipeline-monitor",
                    daemon=True,
                )
                self._monitor_thread.start()
            
            self.logger.info("Monitoring manager started")

    def stop(self) -> None:
        """Stop monitoring manager."""
        with self._lock:
            if not self._active:
                return
            
            self._active = False
            
            # Wait for monitor thread
            if self._monitor_thread:
                self._monitor_thread.join(timeout=5)
            
            # Generate final report
            if self.config.generate_reports:
                self._generate_report()
            
            self.logger.info("Monitoring manager stopped")

    def update_progress(
        self,
        operations: int = 1,
        failed: bool = False,
    ) -> None:
        """Update progress counters.
        
        Args:
            operations: Number of operations completed
            failed: Whether operations failed
        """
        with self._lock:
            self.stats.total_operations += operations
            if failed:
                self.stats.failed_operations += operations
            else:
                self.stats.completed_operations += operations

    def record_latency(
        self,
        latency: float,
    ) -> None:
        """Record operation latency.
        
        Args:
            latency: Operation latency in milliseconds
        """
        with self._lock:
            self.stats.operation_times.append(latency)
            self.stats.peak_latency = max(
                self.stats.peak_latency,
                latency,
            )
            self.stats.avg_latency = (
                sum(self.stats.operation_times) /
                len(self.stats.operation_times)
            )
            
            # Check alert threshold
            if (
                self.config.enable_alerts and
                latency > self.config.alert_thresholds["latency"]
            ):
                self.logger.warning(
                    f"High latency detected: {latency:.1f}ms"
                )

    def record_error(
        self,
        error_type: str,
        error: Exception,
    ) -> None:
        """Record error occurrence.
        
        Args:
            error_type: Type of error
            error: Exception instance
        """
        with self._lock:
            self.stats.total_errors += 1
            self.stats.error_types[error_type] = (
                self.stats.error_types.get(error_type, 0) + 1
            )
            
            # Check alert threshold
            if self.config.enable_alerts:
                error_rate = (
                    self.stats.total_errors /
                    max(self.stats.total_operations, 1)
                )
                if error_rate > self.config.alert_thresholds["error_rate"]:
                    self.logger.warning(
                        f"High error rate detected: {error_rate:.1%}"
                    )

    def update_component_stats(
        self,
        component: str,
        stats: Dict[str, Any],
    ) -> None:
        """Update component statistics.
        
        Args:
            component: Component name
            stats: Component statistics
        """
        with self._lock:
            self.stats.component_stats[component] = stats

    def _monitor_pipeline(self) -> None:
        """Monitor pipeline progress."""
        last_progress = 0
        last_performance = 0
        last_metric = 0
        
        while self._active:
            try:
                now = datetime.now()
                
                # Check progress interval
                if (
                    self.config.monitor_progress and
                    (now - last_progress).seconds >= self.config.progress_interval
                ):
                    self._log_progress()
                    last_progress = now
                
                # Check performance interval
                if (
                    self.config.monitor_performance and
                    (now - last_performance).seconds >= 
                    self.config.performance_interval
                ):
                    self._log_performance()
                    last_performance = now
                
                # Check metric interval
                if (
                    self.config.collect_metrics and
                    (now - last_metric).seconds >= self.config.metric_interval
                ):
                    self._collect_metrics()
                    last_metric = now
                
                # Sleep until next check
                threading.Event().wait(1)
                
            except Exception as e:
                self.logger.error(f"Monitoring error: {str(e)}")
                threading.Event().wait(5)

    def _log_progress(self) -> None:
        """Log pipeline progress."""
        try:
            # Calculate progress
            total = max(self.stats.total_operations, 1)
            completed = self.stats.completed_operations
            failed = self.stats.failed_operations
            success_rate = completed / total
            
            # Log progress
            self.logger.info(
                f"Progress: {completed}/{total} operations "
                f"({success_rate:.1%} success rate)"
            )
            
            # Check failure threshold
            if self.config.enable_alerts:
                failure_rate = failed / total
                if failure_rate > self.config.alert_thresholds["failure_rate"]:
                    self.logger.warning(
                        f"High failure rate detected: {failure_rate:.1%}"
                    )
            
        except Exception as e:
            self.logger.error(f"Progress logging error: {str(e)}")

    def _log_performance(self) -> None:
        """Log performance metrics."""
        try:
            # Calculate metrics
            if self.stats.operation_times:
                avg_latency = self.stats.avg_latency
                peak_latency = self.stats.peak_latency
                
                # Log metrics
                self.logger.info(
                    f"Performance: {avg_latency:.1f}ms avg, "
                    f"{peak_latency:.1f}ms peak"
                )
            
        except Exception as e:
            self.logger.error(f"Performance logging error: {str(e)}")

    def _collect_metrics(self) -> None:
        """Collect pipeline metrics."""
        try:
            # Generate metrics
            metrics = {
                "timestamp": datetime.now().isoformat(),
                "stats": self.stats.to_dict(),
            }
            
            # Save metrics
            if self.config.report_dir:
                metrics_file = (
                    self.config.report_dir /
                    f"metrics_{self._get_timestamp()}.json"
                )
                with open(metrics_file, "w") as f:
                    json.dump(metrics, f, indent=2)
            
        except Exception as e:
            self.logger.error(f"Metric collection error: {str(e)}")

    def _generate_report(self) -> None:
        """Generate monitoring report."""
        try:
            # Generate report
            report = {
                "timestamp": datetime.now().isoformat(),
                "duration": self.stats._get_duration(),
                "stats": self.stats.to_dict(),
            }
            
            # Save report
            if self.config.report_dir:
                report_file = (
                    self.config.report_dir /
                    f"report_{self._get_timestamp()}.{self.config.report_format}"
                )
                
                if self.config.report_format == "json":
                    with open(report_file, "w") as f:
                        json.dump(report, f, indent=2)
                else:
                    with open(report_file, "w") as f:
                        f.write(self._format_report(report))
            
        except Exception as e:
            self.logger.error(f"Report generation error: {str(e)}")

    def _format_report(
        self,
        report: Dict[str, Any],
    ) -> str:
        """Format report as markdown.
        
        Args:
            report: Report data
            
        Returns:
            Formatted report string
        """
        lines = [
            "# Pipeline Monitoring Report",
            "",
            f"Generated: {report['timestamp']}",
            f"Duration: {timedelta(seconds=int(report['duration']))}",
            "",
            "## Progress",
            "",
            "- Total operations: "
            f"{report['stats']['progress']['total']}",
            "- Completed operations: "
            f"{report['stats']['progress']['completed']}",
            "- Failed operations: "
            f"{report['stats']['progress']['failed']}",
            "- Success rate: "
            f"{report['stats']['progress']['success_rate']:.1%}",
            "",
            "## Performance",
            "",
            "- Peak latency: "
            f"{report['stats']['performance']['peak_latency']:.1f}ms",
            "- Average latency: "
            f"{report['stats']['performance']['avg_latency']:.1f}ms",
            "- Operation count: "
            f"{report['stats']['performance']['operation_count']}",
            "",
            "## Errors",
            "",
            "- Total errors: "
            f"{report['stats']['metrics']['errors']['total']}",
            "",
            "### Error Types",
            "",
        ]
        
        # Add error types
        for error_type, count in report["stats"]["metrics"]["errors"][
            "types"
        ].items():
            lines.append(f"- {error_type}: {count}")
        
        return "\n".join(lines)

    def _get_timestamp(self) -> str:
        """Get formatted timestamp string."""
        return datetime.now().strftime("%Y%m%d_%H%M%S")
