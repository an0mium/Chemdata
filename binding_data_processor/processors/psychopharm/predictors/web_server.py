"""Web server functionality for compound data.

This module provides functionality to:
1. Serve web interface
2. Handle API requests
3. Manage WebSocket connections
4. Support real-time updates
5. Handle file uploads/downloads
"""

import logging
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from pathlib import Path
import json
import aiohttp
from aiohttp import web
import aiohttp_cors
from aiohttp_session import setup as setup_session
from aiohttp_session.cookie_storage import EncryptedCookieStorage
import cryptography.fernet

from ....models.validation import ValidationResult
from .data_enrichment import EnrichedData
from .data_analysis import AnalysisResult
from .data_visualization import VisualizationResult
from .web_visualization import WebVisualizationResult
from .web_interface import WebInterfaceResult


@dataclass
class WebServerResult(ValidationResult):
    """Result of web server operation."""
    
    endpoints: Dict[str, str]
    websockets: Dict[str, str]
    sessions: Dict[str, Dict[str, Any]]
    stats: Dict[str, Any]
    issues: List[str]


class WebServer:
    """Web server for compound data."""

    def __init__(
        self,
        host: str = "localhost",
        port: int = 8000,
        static_dir: Optional[str] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize web server."""
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        self.host = host
        self.port = port
        self.static_dir = Path(static_dir) if static_dir else None
        self.app = web.Application()
        self.websockets = {}
        self.sessions = {}

    async def start_server(
        self,
        compounds: List[EnrichedData],
        analysis: AnalysisResult,
        visualization: VisualizationResult,
        web_visualization: WebVisualizationResult,
        web_interface: WebInterfaceResult,
    ) -> WebServerResult:
        """Start web server."""
        self.logger.info(f"Starting server on {self.host}:{self.port}")
        
        try:
            # Set up session storage
            fernet_key = cryptography.fernet.Fernet.generate_key()
            setup_session(
                self.app,
                EncryptedCookieStorage(fernet_key)
            )
            
            # Set up CORS
            aiohttp_cors.setup(self.app, defaults={
                "*": aiohttp_cors.ResourceOptions(
                    allow_credentials=True,
                    expose_headers="*",
                    allow_headers="*",
                ),
            })
            
            # Add routes
            self._add_routes(
                compounds,
                analysis,
                visualization,
                web_visualization,
                web_interface,
            )
            
            # Add static routes
            if self.static_dir:
                self.app.router.add_static(
                    "/static/",
                    self.static_dir,
                )
            
            # Start server
            runner = web.AppRunner(self.app)
            await runner.setup()
            site = web.TCPSite(runner, self.host, self.port)
            await site.start()
            
            # Calculate stats
            stats = self._calculate_stats()
            
            return WebServerResult(
                is_valid=True,
                endpoints=self._get_endpoints(),
                websockets=self.websockets,
                sessions=self.sessions,
                stats=stats,
                issues=[],
            )
            
        except Exception as e:
            self.logger.error(
                f"Error starting server: {str(e)}",
                exc_info=True
            )
            return WebServerResult(
                is_valid=False,
                endpoints={},
                websockets={},
                sessions={},
                stats={},
                issues=[str(e)],
            )

    def _add_routes(
        self,
        compounds: List[EnrichedData],
        analysis: AnalysisResult,
        visualization: VisualizationResult,
        web_visualization: WebVisualizationResult,
        web_interface: WebInterfaceResult,
    ) -> None:
        """Add server routes."""
        # Add page routes
        self.app.router.add_get("/", self._handle_index)
        self.app.router.add_get("/browse", self._handle_browse)
        self.app.router.add_get("/analysis", self._handle_analysis)
        self.app.router.add_get("/settings", self._handle_settings)
        
        # Add API routes
        self.app.router.add_get("/api/compounds", self._handle_compounds)
        self.app.router.add_get(
            "/api/compounds/{id}",
            self._handle_compound_detail
        )
        self.app.router.add_get(
            "/api/analysis/{type}",
            self._handle_analysis_data
        )
        self.app.router.add_get(
            "/api/visualization/{type}",
            self._handle_visualization_data
        )
        
        # Add WebSocket routes
        self.app.router.add_get(
            "/ws/updates",
            self._handle_updates_websocket
        )
        
        # Store data for handlers
        self.app["compounds"] = compounds
        self.app["analysis"] = analysis
        self.app["visualization"] = visualization
        self.app["web_visualization"] = web_visualization
        self.app["web_interface"] = web_interface

    async def _handle_index(
        self,
        request: web.Request,
    ) -> web.Response:
        """Handle index page request."""
        return web.Response(
            text=self.app["web_interface"].templates["base"],
            content_type="text/html",
        )

    async def _handle_browse(
        self,
        request: web.Request,
    ) -> web.Response:
        """Handle browse page request."""
        return web.Response(
            text=self._render_browse_page(),
            content_type="text/html",
        )

    async def _handle_analysis(
        self,
        request: web.Request,
    ) -> web.Response:
        """Handle analysis page request."""
        return web.Response(
            text=self._render_analysis_page(),
            content_type="text/html",
        )

    async def _handle_settings(
        self,
        request: web.Request,
    ) -> web.Response:
        """Handle settings page request."""
        return web.Response(
            text=self._render_settings_page(),
            content_type="text/html",
        )

    async def _handle_compounds(
        self,
        request: web.Request,
    ) -> web.Response:
        """Handle compounds API request."""
        # Get query parameters
        params = request.query
        
        # Filter compounds
        compounds = self._filter_compounds(
            self.app["compounds"],
            params,
        )
        
        # Format response
        response = {
            "compounds": [
                self._format_compound_data(compound)
                for compound in compounds
            ],
            "total": len(compounds),
        }
        
        return web.json_response(response)

    async def _handle_compound_detail(
        self,
        request: web.Request,
    ) -> web.Response:
        """Handle compound detail API request."""
        # Get compound ID
        compound_id = request.match_info["id"]
        
        # Find compound
        compound = next(
            (c for c in self.app["compounds"]
             if c.compound.name == compound_id),
            None,
        )
        
        if not compound:
            raise web.HTTPNotFound()
        
        # Format response
        response = self._format_compound_detail(compound)
        
        return web.json_response(response)

    async def _handle_analysis_data(
        self,
        request: web.Request,
    ) -> web.Response:
        """Handle analysis data API request."""
        # Get analysis type
        analysis_type = request.match_info["type"]
        
        # Get analysis data
        data = self._get_analysis_data(analysis_type)
        
        return web.json_response(data)

    async def _handle_visualization_data(
        self,
        request: web.Request,
    ) -> web.Response:
        """Handle visualization data API request."""
        # Get visualization type
        viz_type = request.match_info["type"]
        
        # Get visualization data
        data = self._get_visualization_data(viz_type)
        
        return web.json_response(data)

    async def _handle_updates_websocket(
        self,
        request: web.Request,
    ) -> web.WebSocketResponse:
        """Handle updates WebSocket connection."""
        ws = web.WebSocketResponse()
        await ws.prepare(request)
        
        # Store WebSocket connection
        session_id = request.cookies.get("session_id")
        if session_id:
            self.websockets[session_id] = ws
        
        try:
            async for msg in ws:
                if msg.type == aiohttp.WSMsgType.TEXT:
                    # Handle message
                    await self._handle_websocket_message(ws, msg.data)
                elif msg.type == aiohttp.WSMsgType.ERROR:
                    self.logger.error(
                        f"WebSocket error: {ws.exception()}"
                    )
        finally:
            # Remove WebSocket connection
            if session_id:
                self.websockets.pop(session_id, None)
        
        return ws

    def _render_browse_page(self) -> str:
        """Render browse page HTML."""
        return self.app["web_interface"].templates["browse"]

    def _render_analysis_page(self) -> str:
        """Render analysis page HTML."""
        return self.app["web_interface"].templates["analysis"]

    def _render_settings_page(self) -> str:
        """Render settings page HTML."""
        return self.app["web_interface"].templates["settings"]

    def _filter_compounds(
        self,
        compounds: List[EnrichedData],
        params: Dict[str, str],
    ) -> List[EnrichedData]:
        """Filter compounds based on parameters."""
        filtered = compounds
        
        # Apply text search
        if "search" in params:
            search = params["search"].lower()
            filtered = [
                c for c in filtered
                if search in c.compound.name.lower()
                or search in c.compound.smiles.lower()
            ]
        
        # Apply target filter
        if "target" in params:
            target = params["target"]
            filtered = [
                c for c in filtered
                if c.targets and c.targets[0] == target
            ]
        
        # Apply activity filter
        if "activity" in params:
            activity = params["activity"]
            filtered = [
                c for c in filtered
                if hasattr(c, "activity_type")
                and c.activity_type == activity
            ]
        
        # Apply property filters
        for prop in {"molecular_weight", "logp", "psa"}:
            if f"{prop}_min" in params:
                min_val = float(params[f"{prop}_min"])
                filtered = [
                    c for c in filtered
                    if c.properties.get(prop, 0) >= min_val
                ]
            if f"{prop}_max" in params:
                max_val = float(params[f"{prop}_max"])
                filtered = [
                    c for c in filtered
                    if c.properties.get(prop, 0) <= max_val
                ]
        
        return filtered

    def _format_compound_data(
        self,
        compound: EnrichedData,
    ) -> Dict[str, Any]:
        """Format compound data for API response."""
        return {
            "name": compound.compound.name,
            "smiles": compound.compound.smiles,
            "target": compound.targets[0] if compound.targets else None,
            "activity": {
                "type": compound.activity_type if hasattr(
                    compound, "activity_type"
                ) else None,
                "value": compound.activity_value if hasattr(
                    compound, "activity_value"
                ) else None,
                "unit": compound.activity_unit if hasattr(
                    compound, "activity_unit"
                ) else None,
            },
            "properties": compound.properties,
            "predictions": {
                pred_type: getattr(compound, f"{pred_type}_predictions", {})
                for pred_type in {
                    "bbb", "activity", "toxicity", "abuse"
                }
            },
        }

    def _format_compound_detail(
        self,
        compound: EnrichedData,
    ) -> Dict[str, Any]:
        """Format compound detail data for API response."""
        return {
            "basic": {
                "name": compound.compound.name,
                "smiles": compound.compound.smiles,
                "properties": compound.properties,
            },
            "activity": {
                "type": compound.activity_type if hasattr(
                    compound, "activity_type"
                ) else None,
                "value": compound.activity_value if hasattr(
                    compound, "activity_value"
                ) else None,
                "unit": compound.activity_unit if hasattr(
                    compound, "activity_unit"
                ) else None,
                "targets": compound.targets,
            },
            "predictions": {
                pred_type: getattr(compound, f"{pred_type}_predictions", {})
                for pred_type in {
                    "bbb", "activity", "toxicity", "abuse"
                }
            },
            "safety": {
                "warnings": compound.warnings if hasattr(
                    compound, "warnings"
                ) else [],
                "risks": compound.risks if hasattr(
                    compound, "risks"
                ) else [],
                "contraindications": compound.contraindications if hasattr(
                    compound, "contraindications"
                ) else [],
            },
            "community": {
                "reports": compound.web_data.get("reports", [])
                if hasattr(compound, "web_data") else [],
                "discussions": compound.web_data.get("discussions", [])
                if hasattr(compound, "web_data") else [],
                "references": compound.web_data.get("references", [])
                if hasattr(compound, "web_data") else [],
            },
        }

    def _get_analysis_data(
        self,
        analysis_type: str,
    ) -> Dict[str, Any]:
        """Get analysis data for API response."""
        analysis = self.app["analysis"]
        
        if analysis_type == "property_stats":
            return analysis.property_stats
        elif analysis_type == "correlations":
            return analysis.correlations
        elif analysis_type == "clusters":
            return analysis.clusters
        elif analysis_type == "outliers":
            return analysis.outliers
        elif analysis_type == "trends":
            return analysis.trends
        else:
            raise web.HTTPBadRequest(
                text=f"Unknown analysis type: {analysis_type}"
            )

    def _get_visualization_data(
        self,
        viz_type: str,
    ) -> Dict[str, Any]:
        """Get visualization data for API response."""
        visualization = self.app["web_visualization"]
        
        if viz_type in visualization.components:
            return visualization.components[viz_type]
        else:
            raise web.HTTPBadRequest(
                text=f"Unknown visualization type: {viz_type}"
            )

    async def _handle_websocket_message(
        self,
        ws: web.WebSocketResponse,
        data: str,
    ) -> None:
        """Handle WebSocket message."""
        try:
            message = json.loads(data)
            
            if message["type"] == "subscribe":
                # Handle subscription request
                await self._handle_subscription(ws, message)
            elif message["type"] == "unsubscribe":
                # Handle unsubscription request
                await self._handle_unsubscription(ws, message)
            else:
                self.logger.warning(
                    f"Unknown message type: {message['type']}"
                )
                
        except json.JSONDecodeError:
            self.logger.error("Invalid JSON message")
        except KeyError:
            self.logger.error("Missing message type")
        except Exception as e:
            self.logger.error(
                f"Error handling message: {str(e)}",
                exc_info=True
            )

    async def _handle_subscription(
        self,
        ws: web.WebSocketResponse,
        message: Dict[str, Any],
    ) -> None:
        """Handle subscription request."""
        # Get subscription details
        subscription = message.get("subscription", {})
        
        # Store subscription
        session_id = message.get("session_id")
        if session_id:
            self.sessions[session_id] = {
                "subscriptions": subscription,
                "websocket": ws,
            }
        
        # Send initial data
        await self._send_subscription_data(ws, subscription)

    async def _handle_unsubscription(
        self,
        ws: web.WebSocketResponse,
        message: Dict[str, Any],
    ) -> None:
        """Handle unsubscription request."""
        # Remove subscription
        session_id = message.get("session_id")
        if session_id:
            self.sessions.pop(session_id, None)

    async def _send_subscription_data(
        self,
        ws: web.WebSocketResponse,
        subscription: Dict[str, Any],
    ) -> None:
        """Send subscription data."""
        data = {}
        
        # Get requested data
        if "compounds" in subscription:
            data["compounds"] = [
                self._format_compound_data(c)
                for c in self.app["compounds"]
            ]
        
        if "analysis" in subscription:
            data["analysis"] = {
                atype: self._get_analysis_data(atype)
                for atype in subscription["analysis"]
            }
        
        if "visualization" in subscription:
            data["visualization"] = {
                vtype: self._get_visualization_data(vtype)
                for vtype in subscription["visualization"]
            }
        
        # Send data
        await ws.send_json({
            "type": "subscription_data",
            "data": data,
        })

    def _get_endpoints(self) -> Dict[str, str]:
        """Get server endpoints."""
        return {
            "index": "/",
            "browse": "/browse",
            "analysis": "/analysis",
            "settings": "/settings",
            "api": {
                "compounds": "/api/compounds",
                "compound_detail": "/api/compounds/{id}",
                "analysis": "/api/analysis/{type}",
                "visualization": "/api/visualization/{type}",
            },
            "websocket": {
                "updates": "/ws/updates",
            },
        }

    def _calculate_stats(self) -> Dict[str, Any]:
        """Calculate server statistics."""
        return {
            "endpoints": {
                "total": len(self._get_endpoints()),
                "types": {
                    "pages": 4,
                    "api": 4,
                    "websocket": 1,
                },
            },
            "websockets": {
                "total": len(self.websockets),
                "active": sum(
                    1 for ws in self.websockets.values()
                    if not ws.closed
                ),
            },
            "sessions": {
                "total": len(self.sessions),
                "subscriptions": sum(
                    len(s.get("subscriptions", {}))
                    for s in self.sessions.values()
                ),
            },
        }
