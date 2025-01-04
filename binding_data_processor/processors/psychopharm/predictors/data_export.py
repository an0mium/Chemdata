"""Data export functionality for enriched compound data.

This module provides functionality to:
1. Export enriched compound data in various formats (TSV, JSON, etc.)
2. Filter and select specific data fields for export
3. Format data for different use cases
4. Handle batch exports
5. Generate summary reports
"""

import logging
from typing import Dict, List, Optional, Any, Set
import pandas as pd
import json
from pathlib import Path
from datetime import datetime

from .data_enrichment import EnrichedData


class DataExporter:
    """Exporter for enriched compound data."""

    # Default field groups for export
    FIELD_GROUPS = {
        "basic": {
            "name",
            "smiles",
            "cas",
            "molecular_weight",
            "logp",
        },
        "activity": {
            "bbb_class",
            "bbb_confidence",
            "activity_predictions",
            "target_predictions",
        },
        "safety": {
            "toxicity_class",
            "toxicity_confidence",
            "abuse_potential",
            "warnings",
            "contraindications",
        },
        "sources": {
            "data_sources",
            "prediction_sources",
            "reference_count",
            "last_updated",
        },
    }

    def __init__(
        self,
        output_dir: Optional[str] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize data exporter."""
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        self.output_dir = Path(output_dir) if output_dir else Path.cwd()

    def export_compounds(
        self,
        compounds: List[EnrichedData],
        output_path: str,
        format: str = "tsv",
        fields: Optional[List[str]] = None,
        field_groups: Optional[List[str]] = None,
        filters: Optional[Dict[str, Any]] = None,
        batch_size: int = 1000,
    ) -> None:
        """Export enriched compound data."""
        self.logger.info(
            f"Exporting {len(compounds)} compounds to {output_path}"
        )

        try:
            # Get fields to export
            export_fields = self._get_export_fields(fields, field_groups)
            
            # Filter compounds if needed
            if filters:
                compounds = self._filter_compounds(compounds, filters)
            
            # Convert to DataFrame
            df = self._compounds_to_dataframe(compounds, export_fields)
            
            # Export in batches
            self._export_batches(
                df, output_path, format, batch_size
            )
            
            self.logger.info("Export completed successfully")
            
        except Exception as e:
            self.logger.error(
                f"Error exporting compounds: {str(e)}",
                exc_info=True
            )
            raise

    def _get_export_fields(
        self,
        fields: Optional[List[str]] = None,
        field_groups: Optional[List[str]] = None,
    ) -> Set[str]:
        """Get fields to export."""
        export_fields = set()
        
        # Add specific fields
        if fields:
            export_fields.update(fields)
        
        # Add field groups
        if field_groups:
            for group in field_groups:
                if group in self.FIELD_GROUPS:
                    export_fields.update(self.FIELD_GROUPS[group])
                else:
                    self.logger.warning(f"Unknown field group: {group}")
        
        # Use all fields if none specified
        if not export_fields:
            export_fields = set().union(*self.FIELD_GROUPS.values())
        
        return export_fields

    def _filter_compounds(
        self,
        compounds: List[EnrichedData],
        filters: Dict[str, Any],
    ) -> List[EnrichedData]:
        """Filter compounds based on criteria."""
        filtered = []
        
        for compound in compounds:
            if self._matches_filters(compound, filters):
                filtered.append(compound)
        
        self.logger.debug(
            f"Filtered {len(compounds)} compounds to {len(filtered)}"
        )
        return filtered

    def _matches_filters(
        self,
        compound: EnrichedData,
        filters: Dict[str, Any],
    ) -> bool:
        """Check if compound matches filter criteria."""
        for field, value in filters.items():
            field_value = self._get_field_value(compound, field)
            if field_value is None:
                return False
            
            if not self._matches_filter_value(field_value, value):
                return False
        
        return True

    def _get_field_value(
        self,
        compound: EnrichedData,
        field: str,
    ) -> Any:
        """Get value of a nested field."""
        current = compound
        for part in field.split("."):
            if hasattr(current, part):
                current = getattr(current, part)
            elif isinstance(current, dict) and part in current:
                current = current[part]
            else:
                return None
        return current

    def _matches_filter_value(
        self,
        field_value: Any,
        filter_value: Any,
    ) -> bool:
        """Check if field value matches filter value."""
        if isinstance(filter_value, (list, set, tuple)):
            return field_value in filter_value
        elif isinstance(filter_value, dict):
            return self._matches_range_filter(field_value, filter_value)
        else:
            return field_value == filter_value

    def _matches_range_filter(
        self,
        value: Any,
        range_filter: Dict[str, Any],
    ) -> bool:
        """Check if value matches range filter."""
        if "min" in range_filter and value < range_filter["min"]:
            return False
        if "max" in range_filter and value > range_filter["max"]:
            return False
        return True

    def _compounds_to_dataframe(
        self,
        compounds: List[EnrichedData],
        fields: Set[str],
    ) -> pd.DataFrame:
        """Convert compounds to DataFrame."""
        rows = []
        
        for compound in compounds:
            row = {}
            
            # Add fields by category
            self._add_basic_fields(compound, fields, row)
            self._add_activity_fields(compound, fields, row)
            self._add_safety_fields(compound, fields, row)
            self._add_source_fields(compound, fields, row)
            
            rows.append(row)
        
        return pd.DataFrame(rows)

    def _add_basic_fields(
        self,
        compound: EnrichedData,
        fields: Set[str],
        row: Dict[str, Any],
    ) -> None:
        """Add basic compound fields to row."""
        if "name" in fields:
            row["name"] = compound.compound.name
        if "smiles" in fields:
            row["smiles"] = compound.compound.smiles
        if "cas" in fields:
            row["cas"] = compound.standardized.cas
        if "molecular_weight" in fields:
            row["molecular_weight"] = compound.properties.get(
                "molecular_weight"
            )
        if "logp" in fields:
            row["logp"] = compound.properties.get("logp")

    def _add_activity_fields(
        self,
        compound: EnrichedData,
        fields: Set[str],
        row: Dict[str, Any],
    ) -> None:
        """Add activity-related fields to row."""
        if "bbb_class" in fields:
            row["bbb_class"] = compound.bbb_predictions.get(
                "class"
            )
        if "bbb_confidence" in fields:
            row["bbb_confidence"] = compound.bbb_predictions.get(
                "confidence"
            )
        if "activity_predictions" in fields:
            row["activity_predictions"] = json.dumps(
                compound.activity_predictions.get("predictions", {})
            )
        if "target_predictions" in fields:
            row["target_predictions"] = json.dumps(
                compound.targets
            )

    def _add_safety_fields(
        self,
        compound: EnrichedData,
        fields: Set[str],
        row: Dict[str, Any],
    ) -> None:
        """Add safety-related fields to row."""
        if "toxicity_class" in fields:
            row["toxicity_class"] = compound.toxicity_predictions.get(
                "class"
            )
        if "toxicity_confidence" in fields:
            row["toxicity_confidence"] = compound.toxicity_predictions.get(
                "confidence"
            )
        if "abuse_potential" in fields:
            row["abuse_potential"] = compound.abuse_predictions.get(
                "potential"
            )
        if "warnings" in fields:
            row["warnings"] = "|".join(compound.warnings)
        if "contraindications" in fields:
            row["contraindications"] = "|".join(compound.contraindications)

    def _add_source_fields(
        self,
        compound: EnrichedData,
        fields: Set[str],
        row: Dict[str, Any],
    ) -> None:
        """Add source-related fields to row."""
        if "data_sources" in fields:
            row["data_sources"] = "|".join(compound.sources)
        if "prediction_sources" in fields:
            row["prediction_sources"] = "|".join(
                k for k, v in compound.bbb_predictions.items()
                if v.get("confidence", 0) > 0
            )
        if "reference_count" in fields:
            row["reference_count"] = sum(
                len(refs) for refs in compound.references.values()
            )
        if "last_updated" in fields:
            row["last_updated"] = max(
                compound.timestamps.values(),
                default=datetime.min,
            ).isoformat()

    def _export_batches(
        self,
        df: pd.DataFrame,
        output_path: str,
        format: str,
        batch_size: int,
    ) -> None:
        """Export DataFrame in batches."""
        output_path = Path(output_path)
        
        # Create parent directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Export based on format
        if format == "tsv":
            self._export_tsv(df, output_path, batch_size)
        elif format == "json":
            self._export_json(df, output_path, batch_size)
        elif format == "csv":
            self._export_csv(df, output_path, batch_size)
        else:
            raise ValueError(f"Unsupported format: {format}")

    def _export_tsv(
        self,
        df: pd.DataFrame,
        output_path: Path,
        batch_size: int,
    ) -> None:
        """Export DataFrame to TSV in batches."""
        total_batches = (len(df) + batch_size - 1) // batch_size
        
        for i in range(total_batches):
            start = i * batch_size
            end = min((i + 1) * batch_size, len(df))
            batch = df.iloc[start:end]
            
            # First batch includes header
            header = i == 0
            mode = "w" if header else "a"
            
            batch.to_csv(
                output_path,
                sep="\t",
                index=False,
                header=header,
                mode=mode,
            )
            
            self.logger.debug(
                f"Exported batch {i + 1}/{total_batches} "
                f"({end}/{len(df)} compounds)"
            )

    def _export_json(
        self,
        df: pd.DataFrame,
        output_path: Path,
        batch_size: int,
    ) -> None:
        """Export DataFrame to JSON in batches."""
        total_batches = (len(df) + batch_size - 1) // batch_size
        
        with output_path.open("w") as f:
            # Start JSON array
            f.write("[\n")
            
            for i in range(total_batches):
                start = i * batch_size
                end = min((i + 1) * batch_size, len(df))
                batch = df.iloc[start:end]
                
                # Convert batch to JSON records
                records = batch.to_dict("records")
                json_str = json.dumps(records, indent=2)
                
                # Remove outer brackets and add commas between batches
                json_str = json_str[1:-1]  # Remove [ and ]
                if i > 0:
                    f.write(",\n")
                f.write(json_str)
                
                self.logger.debug(
                    f"Exported batch {i + 1}/{total_batches} "
                    f"({end}/{len(df)} compounds)"
                )
            
            # End JSON array
            f.write("\n]")

    def _export_csv(
        self,
        df: pd.DataFrame,
        output_path: Path,
        batch_size: int,
    ) -> None:
        """Export DataFrame to CSV in batches."""
        total_batches = (len(df) + batch_size - 1) // batch_size
        
        for i in range(total_batches):
            start = i * batch_size
            end = min((i + 1) * batch_size, len(df))
            batch = df.iloc[start:end]
            
            # First batch includes header
            header = i == 0
            mode = "w" if header else "a"
            
            batch.to_csv(
                output_path,
                index=False,
                header=header,
                mode=mode,
            )
            
            self.logger.debug(
                f"Exported batch {i + 1}/{total_batches} "
                f"({end}/{len(df)} compounds)"
            )

    def generate_summary_report(
        self,
        compounds: List[EnrichedData],
        output_path: str,
    ) -> None:
        """Generate summary report of exported compounds."""
        self.logger.info("Generating summary report")
        
        try:
            # Calculate statistics
            stats = {
                "total_compounds": len(compounds),
                "with_bbb_predictions": sum(
                    1 for c in compounds
                    if c.bbb_predictions.get("confidence", 0) > 0
                ),
                "with_activity_data": sum(
                    1 for c in compounds if c.activities
                ),
                "with_toxicity_data": sum(
                    1 for c in compounds
                    if c.toxicity_predictions.get("confidence", 0) > 0
                ),
                "with_abuse_data": sum(
                    1 for c in compounds
                    if c.abuse_predictions.get("confidence", 0) > 0
                ),
                "with_web_data": sum(
                    1 for c in compounds if c.web_data
                ),
                "data_sources": {
                    source: sum(1 for c in compounds if source in c.sources)
                    for source in set().union(
                        *(c.sources for c in compounds)
                    )
                },
                "average_confidence": {
                    "bbb": sum(
                        c.bbb_predictions.get("confidence", 0)
                        for c in compounds
                    ) / len(compounds),
                    "toxicity": sum(
                        c.toxicity_predictions.get("confidence", 0)
                        for c in compounds
                    ) / len(compounds),
                    "abuse": sum(
                        c.abuse_predictions.get("confidence", 0)
                        for c in compounds
                    ) / len(compounds),
                },
                "timestamp": datetime.now().isoformat(),
            }
            
            # Save report
            with open(output_path, "w") as f:
                json.dump(stats, f, indent=2)
            
            self.logger.info(f"Summary report saved to {output_path}")
            
        except Exception as e:
            self.logger.error(
                f"Error generating summary report: {str(e)}",
                exc_info=True
            )
            raise
