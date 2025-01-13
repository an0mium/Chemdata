"""Export functionality for compound data.

This module provides classes and utilities for exporting compound data in various formats:
- TSV format with essential data (CAS numbers, structures, properties)
- JSON format with complete data (including all metadata)
- Excel format for easier viewing
- SDF format for chemical structure software
- MOL format for individual structures
"""

from .exporter import CompoundExporter

__all__ = ["CompoundExporter"]
