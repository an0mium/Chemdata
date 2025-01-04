"""Infrastructure pipeline coordinator.

This module provides the InfrastructureManager class that:
1. Coordinates infrastructure components
2. Processes compounds in parallel
3. Manages pipeline lifecycle
4. Handles errors and recovery
5. Tracks comprehensive metrics

Example:
    ```python
    # Initialize manager
    manager = InfrastructureManager()
    
    # Process compounds
    results = manager.process_compounds(
        compounds=compounds,
        processor=process_compound,
        checkpoint_key="batch_1"
    )
    
    # Get status
    status = manager.get_status()
    ```
"""

import logging
import threading
from pathlib import Path
from typing import Dict, Optional, Any, List, Union
from dataclasses import dataclass, field
from datetime import datetime
import json
from concurrent.futures import ThreadPoolExecutor

from ...models.compound.enhanced import EnhancedCompound
from ...logger import LogManager
from ...utils.config import PIPELINE_CONFIG
from .circuit_breaker import CircuitBreaker, CircuitBreakerConfig
from .checkpoints import CheckpointManager, CheckpointConfig
from .resources import ResourceManager, ResourceConfig
from .errors import ErrorManager, ErrorConfig
from .monitoring import MonitoringManager, MonitoringConfig


@dataclass
class InfrastructureConfig:
    """Infrastructure configuration."""
    
    # Component settings
    checkpoint_dir: Optional[Path] = None
    error_log: Optional[Path] = None
    metric_dir: Optional[Path] = None
    temp_dir: Optional[Path] = None
    
    # Processing settings
    checkpoint_interval: int = 100  # compounds
    save_checkpoints: bool = True
    
    # Resource limits
    max_workers: int = PIPELINE_CONFIG.get("max_workers", 4)
    max_memory: Optional[int] = PIPELINE_CONFIG.get("max_memory")
    max_temp_size: Optional[int] = PIPELINE_CONFIG.get("max_temp_size")
    
    # Error handling
    raise_errors: bool = False
    max_retries: int = 3
    
    # Monitoring
    track_metrics: bool = True
    track_resources: bool = True
    log_level: str = "INFO"
    
    def __post_init__(self):
        """Initialize directories."""
        for path_attr in ["checkpoint_dir", "error_log", "metric_dir", "temp_dir"]:
            path = getattr(self, path_attr)
            if path:
                path.parent.mkdir(parents=True, exist_ok=True)


@dataclass
class InfrastructureStats:
    """Infrastructure statistics."""
    
    # Component stats
    checkpoint_stats: Dict[str, Any] = field(default_factory=dict)
    resource_stats: Dict[str, Any] = field(default_factory=dict)
    error_stats: Dict[str, Any] = field(default_factory=dict)
    monitoring_stats: Dict[str, Any] = field(default_factory=dict)
    
    # Pipeline stats
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    total_operations: int = 0
    failed_operations: int = 0
    compounds_processed: int = 0
    compounds_failed: int = 0
    total_time: float = 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary format."""
        return {
            "components": {
                "checkpoints": self.checkpoint_stats,
                "resources": self.resource_stats,
                "errors": self.error_stats,
                "monitoring": self.monitoring_stats,
            },
            "pipeline": {
                "start_time": (
                    self.start_time.isoformat()
                    if self.start_time
                    else None
                ),
                "end_time": (
                    self.end_time.isoformat()
                    if self.end_time
                    else None
                ),
                "operations": {
                    "total": self.total_operations,
                    "failed": self.failed_operations,
                    "success_rate": self._get_success_rate(),
                },
                "compounds": {
                    "processed": self.compounds_processed,
                    "failed": self.compounds_failed,
                    "throughput": self._get_throughput(),
                },
            },
        }
    
    def _get_success_rate(self) -> Optional[float]:
        """Get operation success rate."""
        if not self.total_operations:
            return None
        return (
            (self.total_operations - self.failed_operations) /
            self.total_operations
        )
    
    def _get_throughput(self) -> Optional[float]:
        """Get processing throughput."""
        if not self.total_time:
            return None
        return self.compounds_processed / self.total_time


class InfrastructureManager:
    """Manager for pipeline infrastructure."""

    def __init__(
        self,
        config: Optional[InfrastructureConfig] = None,
        logger: Optional[LogManager] = None,
    ):
        """Initialize infrastructure manager.
        
        Args:
            config: Optional configuration
            logger: Optional logger instance
        """
        self.config = config or InfrastructureConfig()
        self.logger = logger or LogManager().get_logger("infrastructure")
        
        # Initialize stats
        self.stats = InfrastructureStats()
        
        # Initialize state
        self._lock = threading.Lock()
        self._active = False
        
        # Initialize components
        self._init_components()

    def _init_components(self) -> None:
        """Initialize infrastructure components."""
        try:
            # Initialize checkpoint manager
            self.checkpoint_manager = CheckpointManager(
                checkpoint_dir=self.config.checkpoint_dir,
                save_interval=self.config.checkpoint_interval,
                logger=self.logger,
            )
            
            # Initialize resource manager
            self.resource_manager = ResourceManager(
                max_workers=self.config.max_workers,
                max_memory=self.config.max_memory,
                logger=self.logger,
            )
            
            # Initialize error manager
            self.error_manager = ErrorManager(
                raise_errors=self.config.raise_errors,
                error_log=self.config.error_log,
                max_retries=self.config.max_retries,
                logger=self.logger,
            )
            
            # Initialize monitoring manager
            self.monitoring_manager = MonitoringManager(
                monitor_resources=self.config.track_resources,
                monitor_progress=self.config.track_metrics,
                log_level=self.config.log_level,
                logger=self.logger,
            )
            
            # Initialize circuit breakers
            self.circuit_breakers: Dict[str, CircuitBreaker] = {}
            
        except Exception as e:
            self.logger.error(f"Failed to initialize components: {str(e)}")
            raise

    def start(self) -> None:
        """Start infrastructure components."""
        with self._lock:
            if self._active:
                return
            
            try:
                # Start components
                self.checkpoint_manager.start()
                self.resource_manager.start()
                self.error_manager.start()
                self.monitoring_manager.start()
                
                # Update stats
                self.stats.start_time = datetime.now()
                self._active = True
                
                self.logger.info("Infrastructure components started")
                
            except Exception as e:
                self.logger.error(f"Failed to start components: {str(e)}")
                self.stop()
                raise

    def stop(self) -> None:
        """Stop infrastructure components."""
        with self._lock:
            if not self._active:
                return
            
            try:
                # Stop components
                self.checkpoint_manager.stop()
                self.resource_manager.stop()
                self.error_manager.stop()
                self.monitoring_manager.stop()
                
                # Update stats
                self.stats.end_time = datetime.now()
                if self.stats.start_time:
                    self.stats.total_time = (
                        datetime.now() - self.stats.start_time
                    ).total_seconds()
                self._update_stats()
                
                self._active = False
                self.logger.info("Infrastructure components stopped")
                
            except Exception as e:
                self.logger.error(f"Failed to stop components: {str(e)}")
                raise

    def process_compounds(
        self,
        compounds: List[EnhancedCompound],
        processor: Any,
        checkpoint_key: Optional[str] = None,
    ) -> List[EnhancedCompound]:
        """Process compounds with infrastructure support.
        
        Args:
            compounds: List of compounds to process
            processor: Processing function or object
            checkpoint_key: Optional key for checkpointing
            
        Returns:
            List of processed compounds
        """
        try:
            # Start pipeline
            self.start()
            
            # Load checkpoint if available
            if checkpoint_key and self.config.save_checkpoints:
                checkpoint = self.checkpoint_manager.load_checkpoint(checkpoint_key)
                if checkpoint:
                    compounds = checkpoint["compounds"]
                    self.stats.checkpoint_stats = checkpoint["stats"]
            
            # Process compounds
            with ThreadPoolExecutor(
                max_workers=self.config.max_workers,
                thread_name_prefix="infra",
            ) as executor:
                # Submit processing tasks
                futures = []
                for compound in compounds:
                    futures.append(
                        executor.submit(
                            self._process_compound,
                            compound,
                            processor,
                        )
                    )
                
                # Process results
                results = []
                for i, future in enumerate(futures):
                    try:
                        result = future.result()
                        if result:
                            results.append(result)
                            self.stats.compounds_processed += 1
                            self.update_progress(1)
                        
                        # Save checkpoint
                        if (
                            checkpoint_key and
                            self.config.save_checkpoints and
                            i > 0 and
                            i % self.config.checkpoint_interval == 0
                        ):
                            self.save_checkpoint(
                                checkpoint_key,
                                {
                                    "compounds": results,
                                    "stats": self.stats.to_dict(),
                                    "timestamp": datetime.now().isoformat(),
                                },
                            )
                            
                    except Exception as e:
                        self.logger.error(
                            f"Error processing compound {i}: {str(e)}"
                        )
                        self.stats.compounds_failed += 1
                        self.update_progress(1, failed=True)
                        self.record_error("compound_processing", e)
            
            return results
            
        except Exception as e:
            self.logger.error(f"Pipeline processing failed: {str(e)}")
            raise
            
        finally:
            # Stop pipeline
            self.stop()

    def _process_compound(
        self,
        compound: EnhancedCompound,
        processor: Any,
    ) -> Optional[EnhancedCompound]:
        """Process single compound with error handling and retries."""
        retries = 0
        while retries <= self.config.max_retries:
            try:
                # Check resources
                self.check_resources()
                
                # Process compound
                start_time = datetime.now()
                result = processor(compound)
                
                # Record latency
                latency = (
                    datetime.now() - start_time
                ).total_seconds() * 1000
                self.record_latency(latency)
                
                return result
                
            except Exception as e:
                retries += 1
                self.stats.total_operations += 1
                self.stats.failed_operations += 1
                
                if retries <= self.config.max_retries:
                    self.logger.warning(
                        f"Retrying compound {compound.name} "
                        f"(attempt {retries}/{self.config.max_retries})"
                    )
                else:
                    self.logger.error(
                        f"Failed to process compound {compound.name} "
                        f"after {retries} attempts: {str(e)}"
                    )
                    self.record_error("compound_processing", e)
                    
                    if self.config.raise_errors:
                        raise
                    
                    return None

    def get_circuit_breaker(
        self,
        name: str,
        config: Optional[CircuitBreakerConfig] = None,
    ) -> CircuitBreaker:
        """Get or create circuit breaker.
        
        Args:
            name: Circuit breaker name
            config: Optional configuration
            
        Returns:
            Circuit breaker instance
        """
        if name not in self.circuit_breakers:
            self.circuit_breakers[name] = CircuitBreaker(
                name=name,
                config=config,
                logger=self.logger,
            )
        return self.circuit_breakers[name]

    def save_checkpoint(
        self,
        key: str,
        data: Dict[str, Any],
    ) -> bool:
        """Save checkpoint data.
        
        Args:
            key: Checkpoint identifier
            data: Data to save
            
        Returns:
            True if successful, False otherwise
        """
        try:
            return self.checkpoint_manager.save_checkpoint(key, data)
        except Exception as e:
            self.error_manager.handle_error(e, "checkpoint_save")
            return False

    def load_checkpoint(
        self,
        key: str,
    ) -> Optional[Dict[str, Any]]:
        """Load checkpoint data.
        
        Args:
            key: Checkpoint identifier
            
        Returns:
            Checkpoint data if successful, None otherwise
        """
        try:
            return self.checkpoint_manager.load_checkpoint(key)
        except Exception as e:
            self.error_manager.handle_error(e, "checkpoint_load")
            return None

    def check_resources(self) -> None:
        """Check resource availability.
        
        Raises:
            RuntimeError: If resources are not available
        """
        try:
            self.resource_manager.check_resources()
        except Exception as e:
            self.error_manager.handle_error(e, "resource_check")
            raise

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
        try:
            self.monitoring_manager.update_progress(operations, failed)
            
            # Update stats
            self.stats.total_operations += operations
            if failed:
                self.stats.failed_operations += operations
                
        except Exception as e:
            self.error_manager.handle_error(e, "progress_update")

    def record_latency(
        self,
        latency: float,
    ) -> None:
        """Record operation latency.
        
        Args:
            latency: Operation latency in milliseconds
        """
        try:
            self.monitoring_manager.record_latency(latency)
        except Exception as e:
            self.error_manager.handle_error(e, "latency_record")

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
        try:
            self.error_manager.record_error(error_type, error)
        except Exception as e:
            self.logger.error(f"Failed to record error: {str(e)}")

    def _update_stats(self) -> None:
        """Update component statistics."""
        try:
            # Get component stats
            self.stats.checkpoint_stats = (
                self.checkpoint_manager.stats.to_dict()
            )
            self.stats.resource_stats = (
                self.resource_manager.stats.to_dict()
            )
            self.stats.error_stats = (
                self.error_manager.stats.to_dict()
            )
            self.stats.monitoring_stats = (
                self.monitoring_manager.stats.to_dict()
            )
            
        except Exception as e:
            self.logger.error(f"Failed to update stats: {str(e)}")

    def get_status(self) -> Dict[str, Any]:
        """Get infrastructure status."""
        self._update_stats()
        return {
            "active": (
                self._active and
                self.checkpoint_manager._active and
                self.resource_manager._active and
                self.error_manager._active and
                self.monitoring_manager._active
            ),
            "uptime": self._get_uptime(),
            "stats": self.stats.to_dict(),
            "circuit_breakers": {
                name: breaker.get_metrics()
                for name, breaker in self.circuit_breakers.items()
            },
        }

    def _get_uptime(self) -> Optional[float]:
        """Get pipeline uptime in seconds."""
        if not self.stats.start_time:
            return None
        return (datetime.now() - self.stats.start_time).total_seconds()
