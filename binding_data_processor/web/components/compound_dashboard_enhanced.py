"""Enhanced compound dashboard component.

This module provides an enhanced web interface that integrates:
- Compound list with search
- Compound details with visualization
- Analysis tools
- Export functionality
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
from datetime import datetime

import pandas as pd
from flask import render_template, request, jsonify

from ...models.compound import Compound
from ..base import BaseComponent, ViewResult
from .compound_list_enhanced import CompoundListEnhanced
from .compound_detail_enhanced import CompoundDetailEnhanced
from .compound_search_enhanced import CompoundSearchEnhanced
from .compound_visualization_enhanced import CompoundVisualizationEnhanced
from .compound_analysis_enhanced import CompoundAnalysisEnhanced
from .compound_export_enhanced import CompoundExportEnhanced


class CompoundDashboardEnhanced(BaseComponent):
    """Enhanced compound dashboard component."""

    def __init__(
        self,
        template_dir: Optional[Path] = None,
        static_dir: Optional[Path] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize dashboard component.
        
        Args:
            template_dir: Optional template directory
            static_dir: Optional static files directory
            logger: Optional logger instance
        """
        super().__init__(template_dir, static_dir, logger)

        # Initialize components
        self.list_view = CompoundListEnhanced(template_dir, static_dir, logger)
        self.detail_view = CompoundDetailEnhanced(template_dir, static_dir, logger)
        self.search = CompoundSearchEnhanced(template_dir, static_dir, logger)
        self.visualization = CompoundVisualizationEnhanced(template_dir, static_dir, logger)
        self.analysis = CompoundAnalysisEnhanced(template_dir, static_dir, logger)
        self.export = CompoundExportEnhanced(template_dir, static_dir, logger)

        # Initialize tracking
        self.dashboard_stats = {
            "total_views": 0,
            "list_views": 0,
            "detail_views": 0,
            "searches": 0,
            "analyses": 0,
            "exports": 0,
            "view_history": [],
        }

    def render_dashboard(
        self,
        compounds: List[Compound],
        selected_compound: Optional[Compound] = None,
        view: str = "list",
        query: Optional[str] = None,
        filters: Optional[Dict[str, Any]] = None,
        analysis_types: Optional[List[str]] = None,
        export_format: Optional[str] = None,
        export_columns: Optional[List[str]] = None,
    ) -> ViewResult:
        """Render dashboard interface.
        
        Args:
            compounds: List of compounds to display
            selected_compound: Optional selected compound for detail view
            view: View type (list/detail)
            query: Optional search query
            filters: Optional filters to apply
            analysis_types: Optional analysis types to perform
            export_format: Optional export format
            export_columns: Optional export columns
            
        Returns:
            ViewResult containing rendered HTML and metadata
        """
        try:
            # Initialize components
            components = {}

            # Handle search/filtering
            if query or filters:
                search_result = self.search.render_search(
                    compounds=compounds,
                    query=query,
                    filters=filters,
                )
                if not search_result.success:
                    return search_result
                compounds = search_result.data["compounds"]
                components["search"] = search_result.data
                self.dashboard_stats["searches"] += 1

            # Handle list view
            if view == "list":
                list_result = self.list_view.render_list(compounds)
                if not list_result.success:
                    return list_result
                components["list"] = list_result.data
                self.dashboard_stats["list_views"] += 1

            # Handle detail view
            elif view == "detail" and selected_compound:
                # Get details
                detail_result = self.detail_view.render_detail(selected_compound)
                if not detail_result.success:
                    return detail_result
                components["detail"] = detail_result.data

                # Get visualizations
                viz_result = self.visualization.render_visualization(selected_compound)
                if not viz_result.success:
                    return viz_result
                components["visualization"] = viz_result.data

                # Get analysis
                if analysis_types:
                    analysis_result = self.analysis.render_analysis(
                        selected_compound,
                        analysis_types=analysis_types,
                    )
                    if not analysis_result.success:
                        return analysis_result
                    components["analysis"] = analysis_result.data
                    self.dashboard_stats["analyses"] += 1

                self.dashboard_stats["detail_views"] += 1

            # Handle export
            if export_format:
                export_result = self.export.render_export(
                    compounds=compounds,
                    format=export_format,
                    columns=export_columns,
                )
                if not export_result.success:
                    return export_result
                components["export"] = export_result.data
                self.dashboard_stats["exports"] += 1

            # Update stats
            self.dashboard_stats["total_views"] += 1
            self.dashboard_stats["view_history"].append({
                "timestamp": datetime.now().isoformat(),
                "view": view,
                "compound": selected_compound.name if selected_compound else None,
                "query": query,
                "filters": filters,
                "analysis_types": analysis_types,
                "export_format": export_format,
            })

            # Render template
            html = render_template(
                "compound_dashboard.html",
                compounds=compounds,
                selected_compound=selected_compound,
                view=view,
                query=query,
                filters=filters or {},
                analysis_types=analysis_types or [],
                export_format=export_format,
                export_columns=export_columns or [],
                components=components,
                stats=self.dashboard_stats,
            )

            return ViewResult(
                success=True,
                data={
                    "html": html,
                    "compounds": compounds,
                    "selected_compound": selected_compound,
                    "view": view,
                    "components": components,
                    "dashboard_stats": self.dashboard_stats,
                },
            )

        except Exception as e:
            self.logger.error(f"Error rendering dashboard: {str(e)}")
            return ViewResult(
                success=False,
                error=str(e),
            )

    def get_metrics(self) -> Dict[str, Any]:
        """Get component metrics."""
        return {
            "dashboard_stats": self.dashboard_stats,
            "list_stats": self.list_view.get_metrics(),
            "detail_stats": self.detail_view.get_metrics(),
            "search_stats": self.search.get_metrics(),
            "visualization_stats": self.visualization.get_metrics(),
            "analysis_stats": self.analysis.get_metrics(),
            "export_stats": self.export.get_metrics(),
        }
