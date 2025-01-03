"""Checkpoint management for data processing.

This module handles:
1. Saving processing checkpoints
2. Loading previous checkpoints
3. Tracking completed steps
4. Managing temporary data files
"""

import os
import json
import pickle
from typing import Any, Dict, List, Optional, Set
from datetime import datetime
import pandas as pd

from logger import LogManager

logger = LogManager().get_logger("checkpoint_manager")


class CheckpointManager:
    """Manages processing checkpoints and temporary data."""
    
    def __init__(self, base_dir: str = None):
        """
        Initialize checkpoint manager.
        
        Args:
            base_dir: Base directory for checkpoint files
        """
        self.base_dir = base_dir or os.path.join(
            os.path.dirname(__file__),
            'checkpoints'
        )
        os.makedirs(self.base_dir, exist_ok=True)
        
        # Track completed steps
        self.completed_steps: Set[str] = set()
        self.step_data: Dict[str, Any] = {}
        
        # Load existing checkpoints
        self._load_checkpoints()

    def _get_checkpoint_path(self, step_name: str) -> str:
        """Get path for checkpoint file."""
        return os.path.join(
            self.base_dir,
            f"{step_name}.checkpoint"
        )

    def _get_data_path(self, step_name: str) -> str:
        """Get path for step data file."""
        return os.path.join(
            self.base_dir,
            f"{step_name}.data"
        )

    def _load_checkpoints(self) -> None:
        """Load existing checkpoints."""
        try:
            checkpoint_file = os.path.join(self.base_dir, 'completed_steps.json')
            if os.path.exists(checkpoint_file):
                with open(checkpoint_file) as f:
                    checkpoint_data = json.load(f)
                    self.completed_steps = set(checkpoint_data.get('completed_steps', []))
                    
                # Load step data
                for step in self.completed_steps:
                    data_path = self._get_data_path(step)
                    if os.path.exists(data_path):
                        try:
                            if data_path.endswith('.csv') or data_path.endswith('.tsv'):
                                self.step_data[step] = pd.read_csv(
                                    data_path,
                                    sep='\t' if data_path.endswith('.tsv') else ','
                                )
                            else:
                                with open(data_path, 'rb') as f:
                                    self.step_data[step] = pickle.load(f)
                        except Exception as e:
                            logger.error(f"Error loading data for step {step}: {str(e)}")
                            
        except Exception as e:
            logger.error(f"Error loading checkpoints: {str(e)}")

    def save_checkpoint(
        self,
        step_name: str,
        data: Any = None,
        metadata: Optional[Dict] = None
    ) -> None:
        """
        Save processing checkpoint.
        
        Args:
            step_name: Name of processing step
            data: Data to save (optional)
            metadata: Additional metadata (optional)
        """
        try:
            # Save step data if provided
            if data is not None:
                data_path = self._get_data_path(step_name)
                if isinstance(data, pd.DataFrame):
                    # Save DataFrame as TSV/CSV
                    is_tsv = any(col for col in data.columns if '\t' in str(col))
                    data_path = data_path.replace('.data', '.tsv' if is_tsv else '.csv')
                    data.to_csv(
                        data_path,
                        sep='\t' if is_tsv else ',',
                        index=False
                    )
                else:
                    # Save other data types with pickle
                    with open(data_path, 'wb') as f:
                        pickle.dump(data, f)
                        
                self.step_data[step_name] = data
            
            # Mark step as completed
            self.completed_steps.add(step_name)
            
            # Save checkpoint info
            checkpoint_data = {
                'completed_steps': list(self.completed_steps),
                'last_updated': datetime.now().isoformat(),
                'metadata': metadata or {}
            }
            
            checkpoint_file = os.path.join(self.base_dir, 'completed_steps.json')
            with open(checkpoint_file, 'w') as f:
                json.dump(checkpoint_data, f, indent=2)
                
            logger.info(f"Saved checkpoint for step: {step_name}")
            
        except Exception as e:
            logger.error(f"Error saving checkpoint for {step_name}: {str(e)}")

    def load_step_data(self, step_name: str) -> Optional[Any]:
        """
        Load data for a completed step.
        
        Args:
            step_name: Name of processing step
            
        Returns:
            Step data if available, None otherwise
        """
        if step_name in self.completed_steps:
            return self.step_data.get(step_name)
        return None

    def is_step_completed(self, step_name: str) -> bool:
        """
        Check if processing step is completed.
        
        Args:
            step_name: Name of processing step
            
        Returns:
            True if step is completed, False otherwise
        """
        return step_name in self.completed_steps

    def clear_checkpoints(self, steps: Optional[List[str]] = None) -> None:
        """
        Clear checkpoints and temporary data.
        
        Args:
            steps: List of steps to clear (all if None)
        """
        try:
            if steps is None:
                # Clear all checkpoints
                for step in self.completed_steps:
                    data_path = self._get_data_path(step)
                    if os.path.exists(data_path):
                        os.remove(data_path)
                self.completed_steps.clear()
                self.step_data.clear()
                
                # Remove checkpoint file
                checkpoint_file = os.path.join(self.base_dir, 'completed_steps.json')
                if os.path.exists(checkpoint_file):
                    os.remove(checkpoint_file)
            else:
                # Clear specific steps
                for step in steps:
                    if step in self.completed_steps:
                        data_path = self._get_data_path(step)
                        if os.path.exists(data_path):
                            os.remove(data_path)
                        self.completed_steps.remove(step)
                        self.step_data.pop(step, None)
                        
                # Update checkpoint file
                checkpoint_data = {
                    'completed_steps': list(self.completed_steps),
                    'last_updated': datetime.now().isoformat()
                }
                
                checkpoint_file = os.path.join(self.base_dir, 'completed_steps.json')
                with open(checkpoint_file, 'w') as f:
                    json.dump(checkpoint_data, f, indent=2)
                    
            logger.info(
                f"Cleared checkpoints for: {', '.join(steps) if steps else 'all steps'}"
            )
            
        except Exception as e:
            logger.error(f"Error clearing checkpoints: {str(e)}")

    def get_step_metadata(self, step_name: str) -> Optional[Dict]:
        """
        Get metadata for a completed step.
        
        Args:
            step_name: Name of processing step
            
        Returns:
            Step metadata if available, None otherwise
        """
        try:
            checkpoint_file = os.path.join(self.base_dir, 'completed_steps.json')
            if os.path.exists(checkpoint_file):
                with open(checkpoint_file) as f:
                    checkpoint_data = json.load(f)
                    return checkpoint_data.get('metadata', {}).get(step_name)
        except Exception as e:
            logger.error(f"Error getting metadata for {step_name}: {str(e)}")
        return None
