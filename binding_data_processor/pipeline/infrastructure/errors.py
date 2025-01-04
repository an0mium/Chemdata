"""Error management.

This module provides the ErrorManager class that:
1. Handles pipeline errors
2. Manages error recovery
3. Implements circuit breaking
4. Tracks error patterns
5. Provides error reporting
"""

import logging
import threading
from pathlib import Path
from typing import Dict, Optional, Any, List, Type, Callable
from dataclasses import dataclass, field
from datetime import datetime, timedelta
import json
import traceback
from functools import wraps


@dataclass
class ErrorConfig:
    """Error configuration."""
    
    # Recovery settings
    max_retries: int = 3
    retry_delay: float = 1.0  # seconds
    retry_backoff: float = 2.0
    
    # Circuit breaker settings
    failure_threshold: int = 5
    reset_timeout: int = 60  # seconds
    half_open_timeout: int = 30  # seconds
    
    # Tracking settings
    track_history: bool = True
    history_size: int = 1000
    error_log: Optional[Path] = None
    
    # Alert settings
    raise_errors: bool = False
    alert_threshold: int = 10  # errors per minute


@dataclass
class ErrorStats:
    """Error statistics."""
    
    # Error counts
    total_errors: int = 0
    handled_errors: int = 0
    unhandled_errors: int = 0
    
    # Recovery stats
    total_retries: int = 0
    successful_retries: int = 0
    failed_retries: int = 0
    
    # Circuit stats
    circuit_trips: int = 0
    circuit_resets: int = 0
    current_state: str = "closed"
    
    # Error history
    error_history: List[Dict[str, Any]] = field(default_factory=list)
    error_types: Dict[str, int] = field(default_factory=dict)
    error_components: Dict[str, int] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary format."""
        return {
            "errors": {
                "total": self.total_errors,
                "handled": self.handled_errors,
                "unhandled": self.unhandled_errors,
            },
            "recovery": {
                "retries": self.total_retries,
                "successful": self.successful_retries,
                "failed": self.failed_retries,
                "success_rate": self._get_retry_rate(),
            },
            "circuit": {
                "trips": self.circuit_trips,
                "resets": self.circuit_resets,
                "state": self.current_state,
            },
            "patterns": {
                "types": self.error_types,
                "components": self.error_components,
            },
        }
    
    def _get_retry_rate(self) -> Optional[float]:
        """Get retry success rate."""
        if not self.total_retries:
            return None
        return self.successful_retries / self.total_retries


class CircuitBreaker:
    """Circuit breaker implementation."""
    
    def __init__(
        self,
        failure_threshold: int,
        reset_timeout: int,
        half_open_timeout: int,
    ):
        """Initialize circuit breaker.
        
        Args:
            failure_threshold: Number of failures before opening
            reset_timeout: Seconds before auto-reset
            half_open_timeout: Seconds in half-open state
        """
        self.failure_threshold = failure_threshold
        self.reset_timeout = reset_timeout
        self.half_open_timeout = half_open_timeout
        
        self._failures = 0
        self._state = "closed"
        self._last_failure = None
        self._lock = threading.Lock()

    def record_failure(self) -> None:
        """Record a failure."""
        with self._lock:
            self._failures += 1
            self._last_failure = datetime.now()
            
            if (
                self._state == "closed" and
                self._failures >= self.failure_threshold
            ):
                self._state = "open"

    def record_success(self) -> None:
        """Record a success."""
        with self._lock:
            if self._state == "half-open":
                self._state = "closed"
            self._failures = 0
            self._last_failure = None

    def allow_request(self) -> bool:
        """Check if request is allowed."""
        with self._lock:
            now = datetime.now()
            
            if self._state == "open":
                # Check reset timeout
                if self._last_failure:
                    elapsed = (now - self._last_failure).total_seconds()
                    if elapsed >= self.reset_timeout:
                        self._state = "half-open"
                        return True
                return False
                
            elif self._state == "half-open":
                # Check half-open timeout
                if self._last_failure:
                    elapsed = (now - self._last_failure).total_seconds()
                    return elapsed >= self.half_open_timeout
                return True
                
            return True

    @property
    def state(self) -> str:
        """Get current circuit state."""
        return self._state


class ErrorManager:
    """Manager for pipeline errors."""

    def __init__(
        self,
        raise_errors: bool = False,
        error_log: Optional[Path] = None,
        max_retries: Optional[int] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize error manager.
        
        Args:
            raise_errors: Whether to raise errors
            error_log: Optional error log file
            max_retries: Optional retry limit override
            logger: Optional logger instance
        """
        self.config = ErrorConfig(
            raise_errors=raise_errors,
            error_log=error_log,
            max_retries=max_retries or ErrorConfig.max_retries,
        )
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        
        # Initialize stats
        self.stats = ErrorStats()
        
        # Initialize state
        self._lock = threading.Lock()
        self._active = False
        
        # Initialize circuit breakers
        self._circuits: Dict[str, CircuitBreaker] = {}
        
        # Initialize error log
        self._init_error_log()

    def _init_error_log(self) -> None:
        """Initialize error log."""
        try:
            if self.config.error_log:
                self.config.error_log.parent.mkdir(parents=True, exist_ok=True)
                
        except Exception as e:
            self.logger.error(f"Failed to initialize error log: {str(e)}")
            raise

    def start(self) -> None:
        """Start error manager."""
        with self._lock:
            if self._active:
                return
            
            self._active = True
            self.logger.info("Error manager started")

    def stop(self) -> None:
        """Stop error manager."""
        with self._lock:
            if not self._active:
                return
            
            self._active = False
            self.logger.info("Error manager stopped")

    def handle_error(
        self,
        error: Exception,
        component: str,
        retry_func: Optional[Callable] = None,
    ) -> Optional[Any]:
        """Handle pipeline error.
        
        Args:
            error: Exception instance
            component: Component name
            retry_func: Optional retry function
            
        Returns:
            Retry result if successful, None otherwise
            
        Raises:
            Exception: If raise_errors is True
        """
        try:
            # Update stats
            self.stats.total_errors += 1
            error_type = error.__class__.__name__
            self.stats.error_types[error_type] = (
                self.stats.error_types.get(error_type, 0) + 1
            )
            self.stats.error_components[component] = (
                self.stats.error_components.get(component, 0) + 1
            )
            
            # Check circuit breaker
            circuit = self._get_circuit(component)
            if not circuit.allow_request():
                self.logger.warning(
                    f"Circuit breaker open for {component}"
                )
                if self.config.raise_errors:
                    raise RuntimeError(
                        f"Circuit breaker open for {component}"
                    )
                return None
            
            # Try recovery
            if retry_func:
                result = self._try_recovery(
                    error=error,
                    component=component,
                    retry_func=retry_func,
                    circuit=circuit,
                )
                if result is not None:
                    self.stats.handled_errors += 1
                    return result
            
            # Log error
            self._log_error(error, component)
            
            # Update stats
            self.stats.unhandled_errors += 1
            
            # Check if should raise
            if self.config.raise_errors:
                raise error
            
            return None
            
        except Exception as e:
            self.logger.error(f"Error handling failed: {str(e)}")
            if self.config.raise_errors:
                raise
            return None

    def _try_recovery(
        self,
        error: Exception,
        component: str,
        retry_func: Callable,
        circuit: CircuitBreaker,
    ) -> Optional[Any]:
        """Try error recovery with retries.
        
        Args:
            error: Original error
            component: Component name
            retry_func: Function to retry
            circuit: Circuit breaker instance
            
        Returns:
            Retry result if successful, None otherwise
        """
        retry_count = 0
        retry_delay = self.config.retry_delay
        
        while retry_count < self.config.max_retries:
            try:
                # Wait before retry
                if retry_count > 0:
                    threading.Event().wait(retry_delay)
                    retry_delay *= self.config.retry_backoff
                
                # Try operation
                result = retry_func()
                
                # Update stats
                self.stats.total_retries += 1
                self.stats.successful_retries += 1
                
                # Update circuit
                circuit.record_success()
                
                return result
                
            except Exception as e:
                retry_count += 1
                self.logger.warning(
                    f"Retry {retry_count} failed for {component}: {str(e)}"
                )
                
                # Update stats
                self.stats.total_retries += 1
                self.stats.failed_retries += 1
                
                # Update circuit
                circuit.record_failure()
                if circuit.state == "open":
                    self.stats.circuit_trips += 1
                
                # Check if should continue
                if retry_count >= self.config.max_retries:
                    self.logger.error(
                        f"Max retries ({self.config.max_retries}) "
                        f"exceeded for {component}"
                    )
                    return None
        
        return None

    def _get_circuit(
        self,
        component: str,
    ) -> CircuitBreaker:
        """Get or create circuit breaker.
        
        Args:
            component: Component name
            
        Returns:
            Circuit breaker instance
        """
        with self._lock:
            if component not in self._circuits:
                self._circuits[component] = CircuitBreaker(
                    failure_threshold=self.config.failure_threshold,
                    reset_timeout=self.config.reset_timeout,
                    half_open_timeout=self.config.half_open_timeout,
                )
            return self._circuits[component]

    def _log_error(
        self,
        error: Exception,
        component: str,
    ) -> None:
        """Log error details.
        
        Args:
            error: Exception instance
            component: Component name
        """
        try:
            # Create error entry
            entry = {
                "timestamp": datetime.now().isoformat(),
                "component": component,
                "error_type": error.__class__.__name__,
                "error_message": str(error),
                "traceback": traceback.format_exc(),
            }
            
            # Update history
            if self.config.track_history:
                self.stats.error_history.append(entry)
                if len(self.stats.error_history) > self.config.history_size:
                    self.stats.error_history.pop(0)
            
            # Write to log file
            if self.config.error_log:
                with open(self.config.error_log, "a") as f:
                    json.dump(entry, f)
                    f.write("\n")
            
            # Log error
            self.logger.error(
                f"Error in {component}: {error.__class__.__name__}: {str(error)}"
            )
            
        except Exception as e:
            self.logger.error(f"Error logging failed: {str(e)}")

    def get_error_history(
        self,
        component: Optional[str] = None,
        error_type: Optional[str] = None,
        start_time: Optional[datetime] = None,
        end_time: Optional[datetime] = None,
    ) -> List[Dict[str, Any]]:
        """Get filtered error history.
        
        Args:
            component: Optional component filter
            error_type: Optional error type filter
            start_time: Optional start time filter
            end_time: Optional end time filter
            
        Returns:
            List of matching error entries
        """
        if not self.config.track_history:
            return []
            
        try:
            # Filter history
            history = self.stats.error_history
            
            if component:
                history = [
                    e for e in history
                    if e["component"] == component
                ]
            
            if error_type:
                history = [
                    e for e in history
                    if e["error_type"] == error_type
                ]
            
            if start_time:
                history = [
                    e for e in history
                    if datetime.fromisoformat(e["timestamp"]) >= start_time
                ]
            
            if end_time:
                history = [
                    e for e in history
                    if datetime.fromisoformat(e["timestamp"]) <= end_time
                ]
            
            return history
            
        except Exception as e:
            self.logger.error(f"Error history retrieval failed: {str(e)}")
            return []

    def get_circuit_states(self) -> Dict[str, str]:
        """Get circuit breaker states.
        
        Returns:
            Dictionary of component circuit states
        """
        with self._lock:
            return {
                component: circuit.state
                for component, circuit in self._circuits.items()
            }

    def reset_circuit(
        self,
        component: str,
    ) -> bool:
        """Reset circuit breaker state.
        
        Args:
            component: Component name
            
        Returns:
            True if successful, False otherwise
        """
        try:
            with self._lock:
                if component in self._circuits:
                    circuit = self._circuits[component]
                    if circuit.state == "open":
                        circuit._state = "closed"
                        circuit._failures = 0
                        circuit._last_failure = None
                        self.stats.circuit_resets += 1
                        self.logger.info(
                            f"Circuit breaker reset for {component}"
                        )
                        return True
            return False
            
        except Exception as e:
            self.logger.error(f"Circuit reset failed: {str(e)}")
            return False
