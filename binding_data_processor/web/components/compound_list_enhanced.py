"""Enhanced compound list view component.

This module provides an enhanced web interface for displaying lists of compounds with:
- Advanced filtering
- Sorting capabilities
- Export functionality
- Interactive visualizations
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
from datetime import datetime

import pandas as pd
import plotly.graph_objects as go
from flask import render_template, request, jsonify

from ...models.compound import Compound
from ..base import BaseComponent, ViewResult


class CompoundListEnhanced(BaseComponent):
    """Enhanced compound list view component."""

    def __init__(
        self,
        template_dir: Optional[Path] = None,
        static_dir: Optional[Path] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize compound list component.
        
        Args:
            template_dir: Optional template directory
            static_dir: Optional static files directory
            logger: Optional logger instance
        """
        super().__init__(template_dir, static_dir, logger)

        # Initialize tracking
        self.view_stats = {
            "total_views": 0,
            "filtered_views": 0,
            "exported_views": 0,
            "visualization_views": 0,
        }

    def render_list(
        self,
        compounds: List[Compound],
        page: int = 1,
        per_page: int = 50,
        sort_by: Optional[str] = None,
        sort_ascending: bool = True,
        filters: Optional[Dict[str, Any]] = None,
    ) -> ViewResult:
        """Render compound list view.
        
        Args:
            compounds: List of compounds to display
            page: Current page number
            per_page: Items per page
            sort_by: Optional field to sort by
            sort_ascending: Sort direction
            filters: Optional filters to apply
            
        Returns:
            ViewResult containing rendered HTML and metadata
        """
        try:
            # Apply filters
            if filters:
                compounds = self._apply_filters(compounds, filters)

            # Sort compounds
            if sort_by:
                compounds = self._sort_compounds(compounds, sort_by, sort_ascending)

            # Paginate
            total_pages = (len(compounds) + per_page - 1) // per_page
            start_idx = (page - 1) * per_page
            end_idx = start_idx + per_page
            page_compounds = compounds[start_idx:end_idx]

            # Generate visualizations
            visualizations = self._generate_visualizations(compounds)

            # Update stats
            self.view_stats["total_views"] += 1
            if filters:
                self.view_stats["filtered_views"] += 1
            if visualizations:
                self.view_stats["visualization_views"] += 1

            # Render template
            html = render_template(
                "compound_list.html",
                compounds=page_compounds,
                total_pages=total_pages,
                current_page=page,
                sort_by=sort_by,
                sort_ascending=sort_ascending,
                filters=filters or {},
                visualizations=visualizations,
                stats=self.view_stats,
            )

            return ViewResult(
                success=True,
                data={
                    "html": html,
                    "compounds": page_compounds,
                    "total_pages": total_pages,
                    "current_page": page,
                    "visualizations": visualizations,
                },
            )

        except Exception as e:
            self.logger.error(f"Error rendering compound list: {str(e)}")
            return ViewResult(
                success=False,
                error=str(e),
            )

    def handle_filter(
        self,
        compounds: List[Compound],
        filters: Dict[str, Any],
    ) -> ViewResult:
        """Handle filter request.
        
        Args:
            compounds: List of compounds to filter
            filters: Filter parameters
            
        Returns:
            ViewResult containing filtered compounds
        """
        try:
            # Apply filters
            filtered = self._apply_filters(compounds, filters)

            # Update stats
            self.view_stats["filtered_views"] += 1

            return ViewResult(
                success=True,
                data={
                    "compounds": filtered,
                    "filter_stats": {
                        "total": len(compounds),
                        "filtered": len(filtered),
                        "filters_applied": list(filters.keys()),
                    },
                },
            )

        except Exception as e:
            self.logger.error(f"Error applying filters: {str(e)}")
            return ViewResult(
                success=False,
                error=str(e),
            )

    def handle_sort(
        self,
        compounds: List[Compound],
        sort_by: str,
        ascending: bool = True,
    ) -> ViewResult:
        """Handle sort request.
        
        Args:
            compounds: List of compounds to sort
            sort_by: Field to sort by
            ascending: Sort direction
            
        Returns:
            ViewResult containing sorted compounds
        """
        try:
            # Sort compounds
            sorted_compounds = self._sort_compounds(compounds, sort_by, ascending)

            return ViewResult(
                success=True,
                data={
                    "compounds": sorted_compounds,
                    "sort_by": sort_by,
                    "ascending": ascending,
                },
            )

        except Exception as e:
            self.logger.error(f"Error sorting compounds: {str(e)}")
            return ViewResult(
                success=False,
                error=str(e),
            )

    def handle_export(
        self,
        compounds: List[Compound],
        format: str = "tsv",
        columns: Optional[List[str]] = None,
    ) -> ViewResult:
        """Handle export request.
        
        Args:
            compounds: List of compounds to export
            format: Export format (tsv/csv/json)
            columns: Optional list of columns to export
            
        Returns:
            ViewResult containing export data
        """
        try:
            # Convert to DataFrame
            data = []
            for compound in compounds:
                row = {
                    "name": compound.name,
                    "smiles": compound.smiles,
                    "cas_number": compound.cas_number,
                }
                if hasattr(compound, "binding_data"):
                    row["binding_data"] = compound.binding_data
                if hasattr(compound, "social_data"):
                    row["social_data"] = compound.social_data
                data.append(row)

            df = pd.DataFrame(data)

            # Filter columns
            if columns:
                df = df[columns]

            # Export
            if format == "tsv":
                export_data = df.to_csv(sep="\t", index=False)
            elif format == "csv":
                export_data = df.to_csv(index=False)
            elif format == "json":
                export_data = df.to_json(orient="records")
            else:
                raise ValueError(f"Unsupported format: {format}")

            # Update stats
            self.view_stats["exported_views"] += 1

            return ViewResult(
                success=True,
                data={
                    "export_data": export_data,
                    "format": format,
                    "columns": list(df.columns),
                    "row_count": len(df),
                },
            )

        except Exception as e:
            self.logger.error(f"Error exporting compounds: {str(e)}")
            return ViewResult(
                success=False,
                error=str(e),
            )

    def _apply_filters(
        self,
        compounds: List[Compound],
        filters: Dict[str, Any],
    ) -> List[Compound]:
        """Apply filters to compound list.
        
        Args:
            compounds: List of compounds to filter
            filters: Filter parameters
            
        Returns:
            Filtered list of compounds
        """
        filtered = compounds

        # Apply name filter
        if "name" in filters:
            filtered = [
                c for c in filtered
                if filters["name"].lower() in c.name.lower()
            ]

        # Apply CAS filter
        if "cas_number" in filters:
            filtered = [
                c for c in filtered
                if filters["cas_number"] in c.cas_number
            ]

        # Apply binding filter
        if "min_binding" in filters:
            filtered = [
                c for c in filtered
                if hasattr(c, "binding_data")
                and any(
                    float(b["affinity"]) >= filters["min_binding"]
                    for b in c.binding_data
                )
            ]

        # Apply social data filter
        if "has_social_data" in filters:
            filtered = [
                c for c in filtered
                if hasattr(c, "social_data") == filters["has_social_data"]
            ]

        return filtered

    def _sort_compounds(
        self,
        compounds: List[Compound],
        sort_by: str,
        ascending: bool = True,
    ) -> List[Compound]:
        """Sort compound list.
        
        Args:
            compounds: List of compounds to sort
            sort_by: Field to sort by
            ascending: Sort direction
            
        Returns:
            Sorted list of compounds
        """
        if sort_by == "name":
            key = lambda c: c.name
        elif sort_by == "cas_number":
            key = lambda c: c.cas_number
        elif sort_by == "binding_count":
            key = lambda c: len(getattr(c, "binding_data", []))
        elif sort_by == "social_count":
            key = lambda c: len(getattr(c, "social_data", {}).get("posts", []))
        else:
            raise ValueError(f"Invalid sort field: {sort_by}")

        return sorted(compounds, key=key, reverse=not ascending)

    def _generate_visualizations(
        self,
        compounds: List[Compound],
    ) -> Dict[str, Any]:
        """Generate visualizations for compound list.
        
        Args:
            compounds: List of compounds to visualize
            
        Returns:
            Dictionary of visualization data
        """
        visualizations = {}

        # Binding distribution
        binding_data = []
        for compound in compounds:
            if hasattr(compound, "binding_data"):
                for binding in compound.binding_data:
                    binding_data.append({
                        "compound": compound.name,
                        "target": binding["target"],
                        "affinity": float(binding["affinity"]),
                    })

        if binding_data:
            df = pd.DataFrame(binding_data)
            fig = go.Figure(data=[
                go.Box(
                    y=df["affinity"],
                    x=df["target"],
                    name="Binding Affinity",
                )
            ])
            fig.update_layout(
                title="Binding Affinity Distribution by Target",
                yaxis_title="Affinity (nM)",
                boxmode="group",
            )
            visualizations["binding_distribution"] = fig.to_json()

        # Social data summary
        social_data = []
        for compound in compounds:
            if hasattr(compound, "social_data"):
                data = compound.social_data
                social_data.append({
                    "compound": compound.name,
                    "reddit_posts": len(data.get("reddit", {}).get("posts", [])),
                    "twitter_mentions": len(data.get("twitter", {}).get("tweets", [])),
                })

        if social_data:
            df = pd.DataFrame(social_data)
            fig = go.Figure(data=[
                go.Bar(
                    name="Reddit Posts",
                    x=df["compound"],
                    y=df["reddit_posts"],
                ),
                go.Bar(
                    name="Twitter Mentions",
                    x=df["compound"],
                    y=df["twitter_mentions"],
                ),
            ])
            fig.update_layout(
                title="Social Media Mentions by Compound",
                barmode="group",
            )
            visualizations["social_summary"] = fig.to_json()

        return visualizations

    def get_metrics(self) -> Dict[str, Any]:
        """Get component metrics."""
        return {
            "view_stats": self.view_stats,
        }
