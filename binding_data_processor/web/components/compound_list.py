"""Web component for displaying compound lists.

This module provides components for:
1. Displaying compound lists with filtering
2. Showing compound details
3. Visualizing compound data
4. Exporting selected compounds
"""

import logging
from typing import List, Optional, Dict, Any
from dataclasses import dataclass
from pathlib import Path

from ...models.compound.enhanced import EnhancedCompound
from ...pipeline.processing.config import ProcessingConfig


@dataclass
class CompoundListConfig:
    """Configuration for compound list component."""
    
    # Display settings
    page_size: int = 50
    show_structures: bool = True
    show_predictions: bool = True
    show_web_data: bool = True
    
    # Filter settings
    filter_sources: bool = True
    filter_predictions: bool = True
    filter_properties: bool = True
    
    # Export settings
    export_formats: List[str] = None
    export_dir: Optional[Path] = None
    
    def __post_init__(self):
        """Initialize configuration."""
        if self.export_formats is None:
            self.export_formats = ["tsv", "json"]


class CompoundList:
    """Component for displaying compound lists."""

    def __init__(
        self,
        config: Optional[CompoundListConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize compound list component.
        
        Args:
            config: Optional component configuration
            logger: Optional logger instance
        """
        self.config = config or CompoundListConfig()
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        
        # Initialize state
        self.compounds = []
        self.filtered_compounds = []
        self.current_page = 1
        self.filters = {}
        self.sort_by = None
        self.sort_ascending = True

    def update_compounds(self, compounds: List[EnhancedCompound]) -> None:
        """Update compound list.
        
        Args:
            compounds: List of compounds to display
        """
        self.compounds = compounds
        self._apply_filters()
        self._sort_compounds()

    def set_filters(self, filters: Dict[str, Any]) -> None:
        """Set active filters.
        
        Args:
            filters: Dictionary of filter settings
        """
        self.filters = filters
        self._apply_filters()
        self.current_page = 1

    def set_sort(self, field: str, ascending: bool = True) -> None:
        """Set sort settings.
        
        Args:
            field: Field to sort by
            ascending: Sort direction
        """
        self.sort_by = field
        self.sort_ascending = ascending
        self._sort_compounds()

    def get_page(self, page: int) -> List[EnhancedCompound]:
        """Get compounds for specified page.
        
        Args:
            page: Page number (1-based)
            
        Returns:
            List of compounds for page
        """
        start = (page - 1) * self.config.page_size
        end = start + self.config.page_size
        return self.filtered_compounds[start:end]

    def get_total_pages(self) -> int:
        """Get total number of pages.
        
        Returns:
            Total number of pages
        """
        return (len(self.filtered_compounds) + self.config.page_size - 1) // self.config.page_size

    def export_compounds(
        self,
        output_file: Path,
        format: str = "tsv",
        columns: Optional[List[str]] = None,
    ) -> None:
        """Export filtered compounds.
        
        Args:
            output_file: Output file path
            format: Export format (tsv/json)
            columns: Optional list of columns to export
        """
        try:
            # Validate format
            if format not in self.config.export_formats:
                raise ValueError(f"Invalid export format: {format}")
            
            # Create output directory
            output_file.parent.mkdir(parents=True, exist_ok=True)
            
            # Export compounds
            if format == "tsv":
                self._export_tsv(output_file, columns)
            else:
                self._export_json(output_file)
                
            self.logger.info(
                f"Exported {len(self.filtered_compounds)} compounds to {output_file}"
            )
            
        except Exception as e:
            self.logger.error(f"Error exporting compounds: {str(e)}")
            raise

    def _apply_filters(self) -> None:
        """Apply active filters to compounds."""
        filtered = self.compounds
        
        for field, value in self.filters.items():
            if field == "source" and self.config.filter_sources:
                filtered = [c for c in filtered if c.source == value]
                
            elif field == "has_predictions" and self.config.filter_predictions:
                filtered = [
                    c for c in filtered
                    if (
                        hasattr(c, "activity_predictions")
                        and c.activity_predictions is not None
                    )
                ]
                
            elif field == "has_web_data" and self.config.filter_web_data:
                filtered = [
                    c for c in filtered
                    if (
                        hasattr(c, "web_data")
                        and c.web_data is not None
                    )
                ]
                
            elif field.startswith("property_") and self.config.filter_properties:
                prop = field.replace("property_", "")
                min_val = value.get("min")
                max_val = value.get("max")
                
                filtered = [
                    c for c in filtered
                    if (
                        hasattr(c, prop)
                        and getattr(c, prop) is not None
                        and (min_val is None or getattr(c, prop) >= min_val)
                        and (max_val is None or getattr(c, prop) <= max_val)
                    )
                ]
        
        self.filtered_compounds = filtered

    def _sort_compounds(self) -> None:
        """Sort filtered compounds."""
        if not self.sort_by:
            return
            
        def get_sort_value(compound):
            if hasattr(compound, self.sort_by):
                return getattr(compound, self.sort_by) or ""
            return ""
        
        self.filtered_compounds.sort(
            key=get_sort_value,
            reverse=not self.sort_ascending,
        )

    def _export_tsv(self, output_file: Path, columns: Optional[List[str]]) -> None:
        """Export compounds to TSV file.
        
        Args:
            output_file: Output file path
            columns: Optional list of columns to export
        """
        import pandas as pd
        
        # Get compound data
        data = []
        for compound in self.filtered_compounds:
            row = {}
            
            # Add basic fields
            row["name"] = compound.name
            row["smiles"] = compound.smiles
            row["source"] = compound.source
            
            # Add predictions if enabled
            if self.config.show_predictions:
                if hasattr(compound, "activity_predictions"):
                    row["activity_predictions"] = str(compound.activity_predictions)
                if hasattr(compound, "toxicity_predictions"):
                    row["toxicity_predictions"] = str(compound.toxicity_predictions)
                if hasattr(compound, "abuse_predictions"):
                    row["abuse_predictions"] = str(compound.abuse_predictions)
                if hasattr(compound, "bbb_predictions"):
                    row["bbb_predictions"] = str(compound.bbb_predictions)
            
            # Add web data if enabled
            if self.config.show_web_data and hasattr(compound, "web_data"):
                row["web_data"] = str(compound.web_data)
            
            data.append(row)
        
        # Create DataFrame
        df = pd.DataFrame(data)
        
        # Filter columns if specified
        if columns:
            df = df[columns]
        
        # Save TSV
        df.to_csv(output_file, sep="\t", index=False)

    def _export_json(self, output_file: Path) -> None:
        """Export compounds to JSON file.
        
        Args:
            output_file: Output file path
        """
        import json
        
        # Convert compounds to dictionaries
        data = [c.to_dict() for c in self.filtered_compounds]
        
        # Save JSON
        with open(output_file, "w") as f:
            json.dump(data, f, indent=2)
