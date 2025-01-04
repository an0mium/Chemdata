"""Resource management.

This module provides the ResourceManager class that:
1. Monitors system resources
2. Enforces resource limits
3. Manages thread pools
4. Handles resource cleanup
5. Tracks resource usage
"""

import logging
import threading
import psutil
from pathlib import Path
from typing import Dict, Optional, Any
from dataclasses import dataclass, field
from datetime import datetime


@dataclass
class ResourceConfig:
    """Resource configuration."""
    
    # Thread settings
    max_workers: int = 4
    min_workers: int = 1
    worker_timeout: int = 30  # seconds
    
    # Memory settings
    max_memory: Optional[int] = None  # bytes
    memory_threshold: float = 0.9  # percentage
    check_interval: int = 5  # seconds
    
    # Storage settings
    temp_dir: Optional[Path] = None
    max_temp_size: Optional[int] = None  # bytes
    cleanup_temp: bool = True


@dataclass
class ResourceStats:
    """Resource statistics."""
    
    # Thread stats
    active_threads: int = 0
    peak_threads: int = 0
    thread_timeouts: int = 0
    
    # Memory stats
    memory_used: int = 0
    peak_memory: int = 0
    memory_warnings: int = 0
    
    # Storage stats
    temp_size: int = 0
    peak_temp_size: int = 0
    cleanup_count: int = 0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary format."""
        return {
            "threads": {
                "active": self.active_threads,
                "peak": self.peak_threads,
                "timeouts": self.thread_timeouts,
            },
            "memory": {
                "used": self.memory_used,
                "peak": self.peak_memory,
                "warnings": self.memory_warnings,
            },
            "storage": {
                "temp_size": self.temp_size,
                "peak_temp": self.peak_temp_size,
                "cleanups": self.cleanup_count,
            },
        }


class ResourceManager:
    """Manager for system resources."""

    def __init__(
        self,
        max_workers: Optional[int] = None,
        max_memory: Optional[int] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize resource manager.
        
        Args:
            max_workers: Optional worker thread limit
            max_memory: Optional memory limit in bytes
            logger: Optional logger instance
        """
        self.config = ResourceConfig(
            max_workers=max_workers or ResourceConfig.max_workers,
            max_memory=max_memory,
        )
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        
        # Initialize stats
        self.stats = ResourceStats()
        
        # Initialize state
        self._lock = threading.Lock()
        self._active = False
        self._monitor_thread = None
        
        # Initialize process tracking
        self._process = psutil.Process()
        self._threads = set()
        
        # Initialize temp directory
        self._init_temp_dir()

    def _init_temp_dir(self) -> None:
        """Initialize temporary directory."""
        try:
            if self.config.temp_dir:
                self.config.temp_dir.mkdir(parents=True, exist_ok=True)
                
        except Exception as e:
            self.logger.error(f"Failed to initialize temp directory: {str(e)}")
            raise

    def start(self) -> None:
        """Start resource manager."""
        with self._lock:
            if self._active:
                return
            
            self._active = True
            
            # Start monitoring thread
            self._monitor_thread = threading.Thread(
                target=self._monitor_resources,
                name="resource-monitor",
                daemon=True,
            )
            self._monitor_thread.start()
            
            self.logger.info("Resource manager started")

    def stop(self) -> None:
        """Stop resource manager."""
        with self._lock:
            if not self._active:
                return
            
            self._active = False
            
            # Wait for monitor thread
            if self._monitor_thread:
                self._monitor_thread.join(timeout=5)
            
            # Clean up resources
            self._cleanup_resources()
            
            self.logger.info("Resource manager stopped")

    def check_resources(self) -> None:
        """Check resource availability.
        
        Raises:
            RuntimeError: If resources are not available
        """
        if not self._active:
            self.logger.warning("Resource manager not active")
            return
            
        try:
            # Check thread limit
            if len(self._threads) >= self.config.max_workers:
                raise RuntimeError(
                    f"Thread limit exceeded: {len(self._threads)} >= "
                    f"{self.config.max_workers}"
                )
            
            # Check memory limit
            if self.config.max_memory:
                memory = self._process.memory_info().rss
                if memory >= self.config.max_memory:
                    raise RuntimeError(
                        f"Memory limit exceeded: {memory} >= "
                        f"{self.config.max_memory}"
                    )
            
            # Check temp storage
            if self.config.max_temp_size and self.config.temp_dir:
                temp_size = sum(
                    f.stat().st_size
                    for f in self.config.temp_dir.glob("**/*")
                    if f.is_file()
                )
                if temp_size >= self.config.max_temp_size:
                    raise RuntimeError(
                        f"Temp storage limit exceeded: {temp_size} >= "
                        f"{self.config.max_temp_size}"
                    )
            
        except Exception as e:
            self.logger.error(f"Resource check failed: {str(e)}")
            raise

    def register_thread(
        self,
        thread: threading.Thread,
    ) -> None:
        """Register worker thread.
        
        Args:
            thread: Thread to register
        """
        with self._lock:
            self._threads.add(thread)
            self.stats.active_threads = len(self._threads)
            self.stats.peak_threads = max(
                self.stats.peak_threads,
                self.stats.active_threads,
            )

    def unregister_thread(
        self,
        thread: threading.Thread,
    ) -> None:
        """Unregister worker thread.
        
        Args:
            thread: Thread to unregister
        """
        with self._lock:
            self._threads.discard(thread)
            self.stats.active_threads = len(self._threads)

    def _monitor_resources(self) -> None:
        """Monitor system resources."""
        while self._active:
            try:
                # Check memory usage
                memory = self._process.memory_info().rss
                self.stats.memory_used = memory
                self.stats.peak_memory = max(
                    self.stats.peak_memory,
                    memory,
                )
                
                # Check memory threshold
                if self.config.max_memory:
                    usage = memory / self.config.max_memory
                    if usage >= self.config.memory_threshold:
                        self.logger.warning(
                            f"Memory usage high: {usage:.1%}"
                        )
                        self.stats.memory_warnings += 1
                
                # Check thread timeouts
                for thread in list(self._threads):
                    if not thread.is_alive():
                        self.unregister_thread(thread)
                        self.stats.thread_timeouts += 1
                
                # Check temp storage
                if self.config.temp_dir:
                    temp_size = sum(
                        f.stat().st_size
                        for f in self.config.temp_dir.glob("**/*")
                        if f.is_file()
                    )
                    self.stats.temp_size = temp_size
                    self.stats.peak_temp_size = max(
                        self.stats.peak_temp_size,
                        temp_size,
                    )
                
                # Sleep until next check
                threading.Event().wait(self.config.check_interval)
                
            except Exception as e:
                self.logger.error(f"Resource monitoring error: {str(e)}")
                threading.Event().wait(self.config.check_interval)

    def _cleanup_resources(self) -> None:
        """Clean up system resources."""
        try:
            # Clean up threads
            for thread in list(self._threads):
                if thread.is_alive():
                    thread.join(timeout=self.config.worker_timeout)
            self._threads.clear()
            
            # Clean up temp files
            if self.config.cleanup_temp and self.config.temp_dir:
                for path in self.config.temp_dir.glob("**/*"):
                    if path.is_file():
                        path.unlink()
                self.stats.cleanup_count += 1
                
        except Exception as e:
            self.logger.error(f"Resource cleanup error: {str(e)}")

    def get_status(self) -> Dict[str, Any]:
        """Get resource status."""
        return {
            "active": self._active,
            "threads": len(self._threads),
            "memory": self._process.memory_info().rss,
            "temp_size": (
                sum(
                    f.stat().st_size
                    for f in self.config.temp_dir.glob("**/*")
                    if f.is_file()
                )
                if self.config.temp_dir
                else 0
            ),
            "stats": self.stats.to_dict(),
        }
