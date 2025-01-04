Web Interface Customization
=======================

This guide shows how to customize the web interface components for specialized visualization and analysis needs.

Custom Components
--------------

Creating custom web components:

.. code-block:: python

    from binding_data_processor.web.components import (
        BaseComponent,
        ComponentConfig,
    )
    from typing import Dict, Any

    class CustomViewer(BaseComponent):
        """Custom compound viewer component."""

        def __init__(self, config: ComponentConfig):
            super().__init__(config)
            self.plot_manager = self._create_plot_manager()

        def render(self, compound: CompoundData) -> str:
            """Render custom view."""
            try:
                # Generate visualizations
                plots = self._generate_plots(compound)
                
                # Build HTML
                html = self.template.render(
                    compound=compound,
                    plots=plots,
                    config=self.config,
                )
                
                return html
                
            except Exception as e:
                self.logger.error(f"Render error: {str(e)}")
                return self._render_error()

        def _generate_plots(self, compound: CompoundData) -> Dict[str, Any]:
            """Generate custom plots."""
            return {
                "structure": self.plot_manager.plot_structure(compound),
                "activity": self.plot_manager.plot_activity(compound),
                "custom": self._create_custom_plot(compound),
            }

        def _create_custom_plot(self, compound: CompoundData) -> Dict[str, Any]:
            """Create custom visualization."""
            # Custom plotting logic
            return {
                "data": [...],
                "layout": {...},
            }

Interactive Components
------------------

Adding interactivity:

.. code-block:: python

    class InteractiveViewer(BaseComponent):
        """Interactive compound viewer."""

        def __init__(self, config: ComponentConfig):
            super().__init__(config)
            self.state = {}

        def render(self, compound: CompoundData) -> str:
            """Render interactive view."""
            return self.template.render(
                compound=compound,
                state=self.state,
                handlers=self._get_handlers(),
            )

        def handle_event(self, event_type: str, data: Dict[str, Any]) -> Dict[str, Any]:
            """Handle user interaction."""
            if event_type == "rotate":
                return self._handle_rotation(data)
            elif event_type == "zoom":
                return self._handle_zoom(data)
            elif event_type == "select":
                return self._handle_selection(data)
            else:
                raise ValueError(f"Unknown event: {event_type}")

        def _handle_rotation(self, data: Dict[str, Any]) -> Dict[str, Any]:
            """Handle structure rotation."""
            angle = data.get("angle", 0)
            self.state["rotation"] = angle
            return {
                "structure": self._update_structure_view(),
            }

        def _handle_zoom(self, data: Dict[str, Any]) -> Dict[str, Any]:
            """Handle zoom action."""
            level = data.get("level", 1.0)
            self.state["zoom"] = level
            return {
                "structure": self._update_structure_view(),
            }

        def _handle_selection(self, data: Dict[str, Any]) -> Dict[str, Any]:
            """Handle atom/bond selection."""
            selection = data.get("selection", [])
            self.state["selection"] = selection
            return {
                "details": self._update_selection_details(),
            }

Search Interface
-------------

Customizing search:

.. code-block:: python

    class AdvancedSearch(BaseComponent):
        """Advanced compound search interface."""

        def __init__(self, config: ComponentConfig):
            super().__init__(config)
            self.search_engine = self._create_search_engine()

        def render(self) -> str:
            """Render search interface."""
            return self.template.render(
                search_fields=self._get_search_fields(),
                recent_searches=self.state.get("recent", []),
            )

        def handle_search(self, query: Dict[str, Any]) -> Dict[str, Any]:
            """Handle search request."""
            try:
                # Perform search
                results = self.search_engine.search(
                    structure=query.get("structure"),
                    properties=query.get("properties"),
                    filters=query.get("filters"),
                )
                
                # Update state
                self._update_recent_searches(query)
                
                return {
                    "results": results,
                    "count": len(results),
                }
                
            except Exception as e:
                self.logger.error(f"Search error: {str(e)}")
                return {"error": str(e)}

        def _get_search_fields(self) -> List[Dict[str, Any]]:
            """Get search field configurations."""
            return [
                {
                    "name": "structure",
                    "type": "structure_editor",
                    "options": {...},
                },
                {
                    "name": "properties",
                    "type": "multi_select",
                    "options": [...],
                },
                {
                    "name": "filters",
                    "type": "filter_builder",
                    "options": {...},
                },
            ]

Export Interface
-------------

Custom export functionality:

.. code-block:: python

    class CustomExport(BaseComponent):
        """Custom export interface."""

        def __init__(self, config: ComponentConfig):
            super().__init__(config)
            self.exporter = self._create_exporter()

        def render(self) -> str:
            """Render export interface."""
            return self.template.render(
                formats=self._get_export_formats(),
                fields=self._get_export_fields(),
            )

        def handle_export(self, request: Dict[str, Any]) -> Dict[str, Any]:
            """Handle export request."""
            try:
                # Prepare export
                result = self.exporter.export(
                    compounds=request["compounds"],
                    format=request["format"],
                    fields=request["fields"],
                    options=request.get("options", {}),
                )
                
                return {
                    "download_url": result.url,
                    "file_size": result.size,
                }
                
            except Exception as e:
                self.logger.error(f"Export error: {str(e)}")
                return {"error": str(e)}

        def _get_export_formats(self) -> List[Dict[str, Any]]:
            """Get available export formats."""
            return [
                {"id": "tsv", "name": "TSV", "description": "..."},
                {"id": "csv", "name": "CSV", "description": "..."},
                {"id": "json", "name": "JSON", "description": "..."},
                {"id": "sdf", "name": "SDF", "description": "..."},
            ]

Dashboard
--------

Creating a custom dashboard:

.. code-block:: python

    class AnalysisDashboard(BaseComponent):
        """Custom analysis dashboard."""

        def __init__(self, config: ComponentConfig):
            super().__init__(config)
            self.components = self._create_components()
            self.layout = self._create_layout()

        def render(self, data: Dict[str, Any]) -> str:
            """Render dashboard."""
            return self.template.render(
                components=self.components,
                layout=self.layout,
                data=data,
            )

        def _create_components(self) -> Dict[str, BaseComponent]:
            """Create dashboard components."""
            return {
                "viewer": CustomViewer(self.config),
                "search": AdvancedSearch(self.config),
                "export": CustomExport(self.config),
                "plots": self._create_plot_components(),
            }

        def _create_plot_components(self) -> Dict[str, BaseComponent]:
            """Create plot components."""
            return {
                "activity": ActivityPlot(self.config),
                "similarity": SimilarityPlot(self.config),
                "property": PropertyPlot(self.config),
            }

        def _create_layout(self) -> Dict[str, Any]:
            """Create dashboard layout."""
            return {
                "type": "grid",
                "rows": [
                    {
                        "height": "60%",
                        "columns": [
                            {"width": "40%", "component": "viewer"},
                            {"width": "60%", "component": "plots"},
                        ],
                    },
                    {
                        "height": "40%",
                        "columns": [
                            {"width": "100%", "component": "search"},
                        ],
                    },
                ],
            }

Integration
---------

Using custom components:

.. code-block:: python

    from binding_data_processor.web import WebApp

    # Configure components
    config = ComponentConfig(
        theme="custom",
        interactive=True,
        cache_enabled=True,
    )

    # Create components
    viewer = CustomViewer(config)
    search = AdvancedSearch(config)
    export = CustomExport(config)
    dashboard = AnalysisDashboard(config)

    # Create web app
    app = WebApp([
        viewer,
        search,
        export,
        dashboard,
    ])

    # Run app
    app.run(
        host="localhost",
        port=8000,
        debug=True,
    )
