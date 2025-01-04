"""Export pipeline.

This module provides the ExportPipeline class that:
1. Handles multiple export formats
2. Supports flexible column selection
3. Validates export data
4. Tracks export statistics
5. Generates export reports
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Set
from dataclasses import dataclass, field
from datetime import datetime
import json
import pandas as pd

from ...models.compound.enhanced import EnhancedCompound


@dataclass
class ExportConfig:
    """Export configuration."""
    
    # Format settings
    format: str = "tsv"  # tsv, json, markdown
    delimiter: str = "\t"
    encoding: str = "utf-8"
    
    # Column settings
    required_columns: Set[str] = field(default_factory=lambda: {
        "name", "smiles", "cas_number"
    })
    optional_columns: Set[str] = field(default_factory=lambda: {
        "molecular_weight", "logp", "tpsa",
        "binding_data", "activity_data", "safety_data",
        "ml_predictions", "web_data", "analysis_results",
    })
    
    # Validation settings
    validate_data: bool = True
    allow_missing: bool = False
    
    # Output settings
    output_dir: Optional[Path] = None
    include_metadata: bool = True
    pretty_print: bool = True


@dataclass
class ExportStats:
    """Export statistics."""
    
    # Export counts
    total_exports: int = 0
    successful_exports: int = 0
    failed_exports: int = 0
    
    # Data stats
    total_compounds: int = 0
    valid_compounds: int = 0
    invalid_compounds: int = 0
    
    # Column stats
    column_counts: Dict[str, int] = field(default_factory=dict)
    missing_counts: Dict[str, int] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary format."""
        return {
            "exports": {
                "total": self.total_exports,
                "successful": self.successful_exports,
                "failed": self.failed_exports,
                "success_rate": self._get_success_rate(),
            },
            "compounds": {
                "total": self.total_compounds,
                "valid": self.valid_compounds,
                "invalid": self.invalid_compounds,
                "validity_rate": self._get_validity_rate(),
            },
            "columns": {
                "counts": self.column_counts,
                "missing": self.missing_counts,
            },
        }
    
    def _get_success_rate(self) -> Optional[float]:
        """Get export success rate."""
        if not self.total_exports:
            return None
        return self.successful_exports / self.total_exports
    
    def _get_validity_rate(self) -> Optional[float]:
        """Get compound validity rate."""
        if not self.total_compounds:
            return None
        return self.valid_compounds / self.total_compounds


class ExportPipeline:
    """Pipeline for data export."""

    def __init__(
        self,
        config: Optional[ExportConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize export pipeline.
        
        Args:
            config: Optional export configuration
            logger: Optional logger instance
        """
        self.config = config or ExportConfig()
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        
        # Initialize stats
        self.stats = ExportStats()
        
        # Initialize output directory
        if self.config.output_dir:
            self.config.output_dir.mkdir(parents=True, exist_ok=True)

    def export_compounds(
        self,
        compounds: List[EnhancedCompound],
        columns: Optional[List[str]] = None,
        output_file: Optional[Path] = None,
    ) -> Optional[Path]:
        """Export compound data.
        
        Args:
            compounds: List of compounds to export
            columns: Optional list of columns to export
            output_file: Optional output file path
            
        Returns:
            Path to export file if successful, None otherwise
        """
        try:
            self.stats.total_exports += 1
            self.stats.total_compounds += len(compounds)
            
            # Validate compounds
            if self.config.validate_data:
                compounds = self._validate_compounds(compounds)
            
            # Get export columns
            columns = self._get_export_columns(columns)
            
            # Convert to export format
            if self.config.format == "tsv":
                result = self._export_tsv(compounds, columns, output_file)
            elif self.config.format == "json":
                result = self._export_json(compounds, columns, output_file)
            else:  # markdown
                result = self._export_markdown(compounds, columns, output_file)
            
            self.stats.successful_exports += 1
            return result
            
        except Exception as e:
            self.logger.error(f"Export failed: {str(e)}")
            self.stats.failed_exports += 1
            return None

    def _validate_compounds(
        self,
        compounds: List[EnhancedCompound],
    ) -> List[EnhancedCompound]:
        """Validate compounds for export."""
        valid_compounds = []
        
        for compound in compounds:
            try:
                # Check required fields
                missing_fields = []
                for field in self.config.required_columns:
                    if not hasattr(compound, field):
                        missing_fields.append(field)
                    elif getattr(compound, field) is None:
                        missing_fields.append(field)
                
                if missing_fields and not self.config.allow_missing:
                    self.logger.warning(
                        f"Compound {compound.name} missing required fields: "
                        f"{', '.join(missing_fields)}"
                    )
                    self.stats.invalid_compounds += 1
                    continue
                
                # Update missing counts
                for field in missing_fields:
                    self.stats.missing_counts[field] = (
                        self.stats.missing_counts.get(field, 0) + 1
                    )
                
                valid_compounds.append(compound)
                self.stats.valid_compounds += 1
                
            except Exception as e:
                self.logger.error(
                    f"Error validating compound {compound.name}: {str(e)}"
                )
                self.stats.invalid_compounds += 1
        
        return valid_compounds

    def _get_export_columns(
        self,
        columns: Optional[List[str]] = None,
    ) -> List[str]:
        """Get columns for export."""
        if columns:
            # Validate custom columns
            invalid_columns = []
            for col in columns:
                if col not in self.config.required_columns:
                    if col not in self.config.optional_columns:
                        invalid_columns.append(col)
            
            if invalid_columns:
                self.logger.warning(
                    f"Invalid columns specified: {', '.join(invalid_columns)}"
                )
            
            # Add required columns
            export_columns = list(self.config.required_columns)
            
            # Add valid optional columns
            for col in columns:
                if col not in export_columns:
                    if col in self.config.optional_columns:
                        export_columns.append(col)
        
        else:
            # Use all columns
            export_columns = list(
                self.config.required_columns | self.config.optional_columns
            )
        
        return export_columns

    def _export_tsv(
        self,
        compounds: List[EnhancedCompound],
        columns: List[str],
        output_file: Optional[Path] = None,
    ) -> Optional[Path]:
        """Export compounds to TSV format."""
        try:
            # Convert to DataFrame
            data = []
            for compound in compounds:
                row = {}
                for col in columns:
                    value = getattr(compound, col, None)
                    if isinstance(value, dict):
                        value = json.dumps(value)
                    row[col] = value
                data.append(row)
            
            df = pd.DataFrame(data)
            
            # Add metadata
            if self.config.include_metadata:
                metadata = {
                    "generated": datetime.now().isoformat(),
                    "compounds": len(compounds),
                    "columns": len(columns),
                    "format": "tsv",
                }
                df.attrs["metadata"] = metadata
            
            # Save to file
            if output_file:
                df.to_csv(
                    output_file,
                    sep=self.config.delimiter,
                    encoding=self.config.encoding,
                    index=False,
                )
                return output_file
            
            # Return DataFrame if no output file
            return df
            
        except Exception as e:
            self.logger.error(f"TSV export failed: {str(e)}")
            return None

    def _export_json(
        self,
        compounds: List[EnhancedCompound],
        columns: List[str],
        output_file: Optional[Path] = None,
    ) -> Optional[Path]:
        """Export compounds to JSON format."""
        try:
            # Convert to JSON
            data = []
            for compound in compounds:
                obj = {}
                for col in columns:
                    value = getattr(compound, col, None)
                    obj[col] = value
                data.append(obj)
            
            # Add metadata
            if self.config.include_metadata:
                export_data = {
                    "metadata": {
                        "generated": datetime.now().isoformat(),
                        "compounds": len(compounds),
                        "columns": len(columns),
                        "format": "json",
                    },
                    "compounds": data,
                }
            else:
                export_data = data
            
            # Save to file
            if output_file:
                with open(output_file, "w", encoding=self.config.encoding) as f:
                    if self.config.pretty_print:
                        json.dump(export_data, f, indent=2)
                    else:
                        json.dump(export_data, f)
                return output_file
            
            # Return JSON string if no output file
            if self.config.pretty_print:
                return json.dumps(export_data, indent=2)
            return json.dumps(export_data)
            
        except Exception as e:
            self.logger.error(f"JSON export failed: {str(e)}")
            return None

    def _export_markdown(
        self,
        compounds: List[EnhancedCompound],
        columns: List[str],
        output_file: Optional[Path] = None,
    ) -> Optional[Path]:
        """Export compounds to Markdown format."""
        try:
            # Generate markdown
            lines = []
            
            # Add metadata
            if self.config.include_metadata:
                lines.extend([
                    "# Compound Export",
                    "",
                    "## Metadata",
                    f"- Generated: {datetime.now().isoformat()}",
                    f"- Compounds: {len(compounds)}",
                    f"- Columns: {len(columns)}",
                    f"- Format: markdown",
                    "",
                ])
            
            # Add table header
            lines.extend([
                "## Compounds",
                "",
                "| " + " | ".join(columns) + " |",
                "| " + " | ".join(["---"] * len(columns)) + " |",
            ])
            
            # Add compound rows
            for compound in compounds:
                row = []
                for col in columns:
                    value = getattr(compound, col, "")
                    if isinstance(value, dict):
                        value = json.dumps(value)
                    row.append(str(value))
                lines.append("| " + " | ".join(row) + " |")
            
            markdown = "\n".join(lines)
            
            # Save to file
            if output_file:
                output_file.write_text(markdown, encoding=self.config.encoding)
                return output_file
            
            # Return markdown string if no output file
            return markdown
            
        except Exception as e:
            self.logger.error(f"Markdown export failed: {str(e)}")
            return None

    def get_column_info(self) -> Dict[str, Any]:
        """Get information about export columns."""
        return {
            "required": list(self.config.required_columns),
            "optional": list(self.config.optional_columns),
            "stats": {
                "counts": self.stats.column_counts,
                "missing": self.stats.missing_counts,
            },
        }
