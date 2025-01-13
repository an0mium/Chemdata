"""Utilities module for structure analysis system.

This module provides:
1. Common utilities
2. Helper functions
3. Standard operations
4. Shared tools
"""

import logging
import time
from typing import Dict, List, Optional, Any, TypeVar, Generic, Callable, Union, Tuple
from pathlib import Path
import json
import hashlib
import functools
import threading
from datetime import datetime
import os
import tempfile
import shutil
import contextlib

from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom

from .logging import get_logger
from .errors import ValidationError

logger = get_logger(__name__)

T = TypeVar("T")


def safe_filename(name: str) -> str:
    """Create safe filename from string.

    Args:
        name: Input string

    Returns:
        Safe filename
    """
    # Replace unsafe characters
    safe = name.replace(" ", "_")
    safe = "".join(c for c in safe if c.isalnum() or c in "._-")
    return safe


def hash_object(obj: Any) -> str:
    """Generate hash for object.

    Args:
        obj: Object to hash

    Returns:
        Hash string
    """
    # Convert to JSON-serializable form
    if hasattr(obj, "to_dict"):
        obj = obj.to_dict()
    elif hasattr(obj, "__dict__"):
        obj = obj.__dict__

    # Generate hash
    json_str = json.dumps(obj, sort_keys=True)
    return hashlib.sha256(json_str.encode()).hexdigest()


@contextlib.contextmanager
def temp_directory() -> Path:
    """Create temporary directory.

    Yields:
        Path to temporary directory
    """
    path = Path(tempfile.mkdtemp())
    try:
        yield path
    finally:
        shutil.rmtree(path)


@contextlib.contextmanager
def temp_file(suffix: Optional[str] = None) -> Path:
    """Create temporary file.

    Args:
        suffix: Optional file suffix

    Yields:
        Path to temporary file
    """
    fd, path = tempfile.mkstemp(suffix=suffix)
    try:
        os.close(fd)
        yield Path(path)
    finally:
        os.unlink(path)


def retry(
    max_attempts: int = 3,
    delay: float = 1.0,
    backoff: float = 2.0,
    exceptions: Union[Type[Exception], Tuple[Type[Exception], ...]] = Exception,
):
    """Retry decorator with exponential backoff.

    Args:
        max_attempts: Maximum number of attempts
        delay: Initial delay between attempts
        backoff: Backoff multiplier
        exceptions: Exception types to catch

    Returns:
        Decorated function
    """

    def decorator(func: Callable):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            last_exception = None
            current_delay = delay

            for attempt in range(max_attempts):
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    last_exception = e
                    if attempt < max_attempts - 1:
                        time.sleep(current_delay)
                        current_delay *= backoff
                    continue

            raise last_exception

        return wrapper

    return decorator


def memoize(func: Callable):
    """Memoization decorator.

    Args:
        func: Function to memoize

    Returns:
        Decorated function
    """
    cache = {}
    lock = threading.Lock()

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Create cache key
        key_parts = [str(arg) for arg in args]
        key_parts.extend(f"{k}={v}" for k, v in sorted(kwargs.items()))
        key = "|".join(key_parts)

        with lock:
            if key not in cache:
                cache[key] = func(*args, **kwargs)
            return cache[key]

    return wrapper


def get_structure_info(structure: Structure) -> Dict[str, Any]:
    """Get structure information.

    Args:
        structure: Structure to analyze

    Returns:
        Structure information
    """
    info = {
        "id": structure.id,
        "models": len(structure),
        "chains": {},
        "residues": 0,
        "atoms": 0,
    }

    for model in structure:
        for chain in model:
            chain_info = {
                "id": chain.id,
                "residues": len(chain),
                "atoms": sum(len(residue) for residue in chain),
            }
            info["chains"][chain.id] = chain_info
            info["residues"] += chain_info["residues"]
            info["atoms"] += chain_info["atoms"]

    return info


def get_atom_coordinates(atom: Atom) -> Tuple[float, float, float]:
    """Get atom coordinates.

    Args:
        atom: Atom to get coordinates for

    Returns:
        Tuple of x, y, z coordinates
    """
    return atom.get_coord()


def get_residue_center(residue: Residue) -> Tuple[float, float, float]:
    """Get residue center coordinates.

    Args:
        residue: Residue to get center for

    Returns:
        Tuple of x, y, z coordinates
    """
    coords = [get_atom_coordinates(atom) for atom in residue]
    if not coords:
        raise ValueError(f"Residue {residue.id} has no atoms")

    # Calculate average coordinates
    x = sum(c[0] for c in coords) / len(coords)
    y = sum(c[1] for c in coords) / len(coords)
    z = sum(c[2] for c in coords) / len(coords)

    return (x, y, z)


def get_chain_center(chain: Chain) -> Tuple[float, float, float]:
    """Get chain center coordinates.

    Args:
        chain: Chain to get center for

    Returns:
        Tuple of x, y, z coordinates
    """
    centers = []
    for residue in chain:
        try:
            centers.append(get_residue_center(residue))
        except ValueError:
            continue

    if not centers:
        raise ValueError(f"Chain {chain.id} has no valid residues")

    # Calculate average coordinates
    x = sum(c[0] for c in centers) / len(centers)
    y = sum(c[1] for c in centers) / len(centers)
    z = sum(c[2] for c in centers) / len(centers)

    return (x, y, z)


def get_structure_center(structure: Structure) -> Tuple[float, float, float]:
    """Get structure center coordinates.

    Args:
        structure: Structure to get center for

    Returns:
        Tuple of x, y, z coordinates
    """
    centers = []
    model = structure[0]
    for chain in model:
        try:
            centers.append(get_chain_center(chain))
        except ValueError:
            continue

    if not centers:
        raise ValueError("Structure has no valid chains")

    # Calculate average coordinates
    x = sum(c[0] for c in centers) / len(centers)
    y = sum(c[1] for c in centers) / len(centers)
    z = sum(c[2] for c in centers) / len(centers)

    return (x, y, z)


def distance(p1: Tuple[float, float, float], p2: Tuple[float, float, float]) -> float:
    """Calculate distance between 3D points.

    Args:
        p1: First point coordinates
        p2: Second point coordinates

    Returns:
        Distance between points
    """
    return ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2) ** 0.5


def format_timestamp(timestamp: float) -> str:
    """Format timestamp as string.

    Args:
        timestamp: Unix timestamp

    Returns:
        Formatted timestamp string
    """
    dt = datetime.fromtimestamp(timestamp)
    return dt.strftime("%Y-%m-%d %H:%M:%S")


def parse_timestamp(timestamp_str: str) -> float:
    """Parse timestamp string.

    Args:
        timestamp_str: Timestamp string

    Returns:
        Unix timestamp
    """
    dt = datetime.strptime(timestamp_str, "%Y-%m-%d %H:%M:%S")
    return dt.timestamp()


def format_duration(seconds: float) -> str:
    """Format duration in seconds.

    Args:
        seconds: Duration in seconds

    Returns:
        Formatted duration string
    """
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    seconds = seconds % 60

    parts = []
    if hours > 0:
        parts.append(f"{hours}h")
    if minutes > 0:
        parts.append(f"{minutes}m")
    if seconds > 0 or not parts:
        parts.append(f"{seconds:.1f}s")

    return " ".join(parts)


def format_size(size: int) -> str:
    """Format size in bytes.

    Args:
        size: Size in bytes

    Returns:
        Formatted size string
    """
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024:
            break
        size /= 1024
    return f"{size:.1f} {unit}"


def ensure_directory(path: Union[str, Path]) -> Path:
    """Ensure directory exists.

    Args:
        path: Directory path

    Returns:
        Path object
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def copy_file(src: Union[str, Path], dst: Union[str, Path]):
    """Copy file with directory creation.

    Args:
        src: Source path
        dst: Destination path
    """
    src = Path(src)
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def move_file(src: Union[str, Path], dst: Union[str, Path]):
    """Move file with directory creation.

    Args:
        src: Source path
        dst: Destination path
    """
    src = Path(src)
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(src, dst)


def remove_file(path: Union[str, Path]):
    """Remove file if it exists.

    Args:
        path: File path
    """
    path = Path(path)
    if path.exists():
        path.unlink()


def remove_directory(path: Union[str, Path]):
    """Remove directory if it exists.

    Args:
        path: Directory path
    """
    path = Path(path)
    if path.exists():
        shutil.rmtree(path)


def list_files(
    path: Union[str, Path],
    pattern: str = "*",
    recursive: bool = False,
) -> List[Path]:
    """List files in directory.

    Args:
        path: Directory path
        pattern: Glob pattern
        recursive: Whether to search recursively

    Returns:
        List of file paths
    """
    path = Path(path)
    if recursive:
        return list(path.rglob(pattern))
    return list(path.glob(pattern))


def read_json(path: Union[str, Path]) -> Any:
    """Read JSON file.

    Args:
        path: File path

    Returns:
        Parsed JSON data
    """
    path = Path(path)
    with open(path) as f:
        return json.load(f)


def write_json(path: Union[str, Path], data: Any):
    """Write JSON file.

    Args:
        path: File path
        data: Data to write
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def read_text(path: Union[str, Path]) -> str:
    """Read text file.

    Args:
        path: File path

    Returns:
        File contents
    """
    path = Path(path)
    with open(path) as f:
        return f.read()


def write_text(path: Union[str, Path], text: str):
    """Write text file.

    Args:
        path: File path
        text: Text to write
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(text)
