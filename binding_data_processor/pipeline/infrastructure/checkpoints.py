"""Checkpoint management.

This module provides the CheckpointManager class that:
1. Saves pipeline state
2. Loads pipeline state
3. Manages checkpoint files
4. Handles recovery
5. Tracks checkpoint history
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Any
from dataclasses import dataclass, field
from datetime import datetime
import json
import shutil
import threading


@dataclass
class CheckpointConfig:
    """Checkpoint configuration."""
    
    # Directory settings
    checkpoint_dir: Optional[Path] = None
    backup_dir: Optional[Path] = None
    
    # Timing settings
    save_interval: int = 100  # operations
    max_age: Optional[int] = None  # seconds
    
    # Storage settings
    max_checkpoints: Optional[int] = None
    compress: bool = False
    
    # Recovery settings
    auto_recover: bool = True
    validate_checkpoints: bool = True


@dataclass
class CheckpointStats:
    """Checkpoint statistics."""
    
    # Operation counts
    total_saves: int = 0
    total_loads: int = 0
    total_deletes: int = 0
    
    # Storage stats
    total_size: int = 0
    checkpoint_count: int = 0
    oldest_checkpoint: Optional[str] = None
    newest_checkpoint: Optional[str] = None
    
    # Error stats
    save_errors: int = 0
    load_errors: int = 0
    validation_errors: int = 0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary format."""
        return {
            "operations": {
                "saves": self.total_saves,
                "loads": self.total_loads,
                "deletes": self.total_deletes,
            },
            "storage": {
                "size": self.total_size,
                "count": self.checkpoint_count,
                "oldest": self.oldest_checkpoint,
                "newest": self.newest_checkpoint,
            },
            "errors": {
                "saves": self.save_errors,
                "loads": self.load_errors,
                "validation": self.validation_errors,
            },
        }


class CheckpointManager:
    """Manager for pipeline checkpoints."""

    def __init__(
        self,
        checkpoint_dir: Optional[Path] = None,
        save_interval: Optional[int] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize checkpoint manager.
        
        Args:
            checkpoint_dir: Optional checkpoint directory
            save_interval: Optional save interval override
            logger: Optional logger instance
        """
        self.config = CheckpointConfig(
            checkpoint_dir=checkpoint_dir,
            save_interval=save_interval or CheckpointConfig.save_interval,
        )
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        
        # Initialize stats
        self.stats = CheckpointStats()
        
        # Initialize state
        self._lock = threading.Lock()
        self._active = False
        self._operation_count = 0
        
        # Initialize directories
        self._init_directories()

    def _init_directories(self) -> None:
        """Initialize checkpoint directories."""
        try:
            # Create checkpoint directory
            if self.config.checkpoint_dir:
                self.config.checkpoint_dir.mkdir(parents=True, exist_ok=True)
            
            # Create backup directory
            if self.config.backup_dir:
                self.config.backup_dir.mkdir(parents=True, exist_ok=True)
            
            # Load existing checkpoints
            self._load_existing_checkpoints()
            
        except Exception as e:
            self.logger.error(f"Failed to initialize directories: {str(e)}")
            raise

    def _load_existing_checkpoints(self) -> None:
        """Load information about existing checkpoints."""
        if not self.config.checkpoint_dir:
            return
            
        try:
            # Get checkpoint files
            checkpoint_files = list(self.config.checkpoint_dir.glob("*.json"))
            self.stats.checkpoint_count = len(checkpoint_files)
            
            if not checkpoint_files:
                return
            
            # Get timestamps
            timestamps = []
            total_size = 0
            for path in checkpoint_files:
                timestamp = datetime.fromtimestamp(path.stat().st_mtime)
                timestamps.append((timestamp, path))
                total_size += path.stat().st_size
            
            # Update stats
            self.stats.total_size = total_size
            if timestamps:
                oldest = min(timestamps, key=lambda x: x[0])
                newest = max(timestamps, key=lambda x: x[0])
                self.stats.oldest_checkpoint = oldest[0].isoformat()
                self.stats.newest_checkpoint = newest[0].isoformat()
            
        except Exception as e:
            self.logger.error(f"Failed to load existing checkpoints: {str(e)}")

    def start(self) -> None:
        """Start checkpoint manager."""
        with self._lock:
            if self._active:
                return
            
            self._active = True
            self._operation_count = 0
            self.logger.info("Checkpoint manager started")

    def stop(self) -> None:
        """Stop checkpoint manager."""
        with self._lock:
            if not self._active:
                return
            
            self._active = False
            self.logger.info("Checkpoint manager stopped")

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
        if not self._active:
            self.logger.warning("Checkpoint manager not active")
            return False
            
        try:
            # Validate checkpoint directory
            if not self.config.checkpoint_dir:
                self.logger.error("No checkpoint directory configured")
                return False
            
            # Create checkpoint path
            checkpoint_path = (
                self.config.checkpoint_dir / f"{key}_{self._get_timestamp()}.json"
            )
            
            # Create backup if exists
            if checkpoint_path.exists() and self.config.backup_dir:
                backup_path = (
                    self.config.backup_dir / 
                    f"{key}_{self._get_timestamp()}_backup.json"
                )
                shutil.copy2(checkpoint_path, backup_path)
            
            # Save checkpoint
            with open(checkpoint_path, "w") as f:
                json.dump(
                    {
                        "key": key,
                        "timestamp": datetime.now().isoformat(),
                        "data": data,
                    },
                    f,
                    indent=2 if not self.config.compress else None,
                )
            
            # Update stats
            self.stats.total_saves += 1
            self.stats.checkpoint_count += 1
            self.stats.total_size += checkpoint_path.stat().st_size
            self.stats.newest_checkpoint = datetime.now().isoformat()
            
            # Clean old checkpoints
            self._clean_old_checkpoints()
            
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to save checkpoint: {str(e)}")
            self.stats.save_errors += 1
            return False

    def load_checkpoint(
        self,
        key: str,
        validate: Optional[bool] = None,
    ) -> Optional[Dict[str, Any]]:
        """Load checkpoint data.
        
        Args:
            key: Checkpoint identifier
            validate: Optional validation override
            
        Returns:
            Checkpoint data if successful, None otherwise
        """
        if not self._active:
            self.logger.warning("Checkpoint manager not active")
            return None
            
        try:
            # Validate checkpoint directory
            if not self.config.checkpoint_dir:
                self.logger.error("No checkpoint directory configured")
                return None
            
            # Find latest checkpoint
            checkpoint_files = list(
                self.config.checkpoint_dir.glob(f"{key}_*.json")
            )
            if not checkpoint_files:
                return None
            
            latest = max(
                checkpoint_files,
                key=lambda p: datetime.fromtimestamp(p.stat().st_mtime),
            )
            
            # Load checkpoint
            with open(latest) as f:
                checkpoint = json.load(f)
            
            # Validate checkpoint
            if validate or (
                validate is None and self.config.validate_checkpoints
            ):
                if not self._validate_checkpoint(checkpoint):
                    self.logger.error("Invalid checkpoint data")
                    self.stats.validation_errors += 1
                    return None
            
            # Update stats
            self.stats.total_loads += 1
            
            return checkpoint["data"]
            
        except Exception as e:
            self.logger.error(f"Failed to load checkpoint: {str(e)}")
            self.stats.load_errors += 1
            return None

    def delete_checkpoint(
        self,
        key: str,
    ) -> bool:
        """Delete checkpoint data.
        
        Args:
            key: Checkpoint identifier
            
        Returns:
            True if successful, False otherwise
        """
        if not self._active:
            self.logger.warning("Checkpoint manager not active")
            return False
            
        try:
            # Validate checkpoint directory
            if not self.config.checkpoint_dir:
                self.logger.error("No checkpoint directory configured")
                return False
            
            # Find checkpoint files
            checkpoint_files = list(
                self.config.checkpoint_dir.glob(f"{key}_*.json")
            )
            if not checkpoint_files:
                return True
            
            # Delete files
            for path in checkpoint_files:
                # Create backup if configured
                if self.config.backup_dir:
                    backup_path = (
                        self.config.backup_dir / 
                        f"{path.stem}_deleted_{self._get_timestamp()}.json"
                    )
                    shutil.copy2(path, backup_path)
                
                # Delete file
                path.unlink()
                self.stats.total_size -= path.stat().st_size
                self.stats.checkpoint_count -= 1
            
            # Update stats
            self.stats.total_deletes += len(checkpoint_files)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to delete checkpoint: {str(e)}")
            return False

    def _clean_old_checkpoints(self) -> None:
        """Clean old checkpoint files."""
        if not self.config.checkpoint_dir:
            return
            
        try:
            # Check max checkpoints
            if self.config.max_checkpoints:
                checkpoint_files = list(
                    self.config.checkpoint_dir.glob("*.json")
                )
                if len(checkpoint_files) > self.config.max_checkpoints:
                    # Sort by modification time
                    checkpoint_files.sort(
                        key=lambda p: p.stat().st_mtime
                    )
                    
                    # Delete oldest files
                    for path in checkpoint_files[
                        :-(self.config.max_checkpoints)
                    ]:
                        self.delete_checkpoint(path.stem.split("_")[0])
            
            # Check max age
            if self.config.max_age:
                now = datetime.now().timestamp()
                for path in self.config.checkpoint_dir.glob("*.json"):
                    age = now - path.stat().st_mtime
                    if age > self.config.max_age:
                        self.delete_checkpoint(path.stem.split("_")[0])
                        
        except Exception as e:
            self.logger.error(f"Failed to clean old checkpoints: {str(e)}")

    def _validate_checkpoint(
        self,
        checkpoint: Dict[str, Any],
    ) -> bool:
        """Validate checkpoint data."""
        try:
            # Check required fields
            if not all(
                field in checkpoint
                for field in ["key", "timestamp", "data"]
            ):
                return False
            
            # Validate timestamp
            try:
                datetime.fromisoformat(checkpoint["timestamp"])
            except ValueError:
                return False
            
            # Validate data
            if not isinstance(checkpoint["data"], dict):
                return False
            
            return True
            
        except Exception:
            return False

    def _get_timestamp(self) -> str:
        """Get formatted timestamp string."""
        return datetime.now().strftime("%Y%m%d_%H%M%S")
