"""Web application functionality for compound data.

This module provides functionality to:
1. Initialize web application
2. Configure components
3. Start services
4. Handle lifecycle
5. Manage state
"""

import logging
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from pathlib import Path
import json
import tempfile
import shutil

from ....models.validation import ValidationResult
from .data_enrichment import EnrichedData
from .data_analysis import DataAnalyzer
from .data_visualization import DataVisualizer
from .web_visualization import WebVisualizer
from .web_interface import WebInterface
from .web_server import WebServerResult, WebServer


@dataclass
class WebAppResult(ValidationResult):
    """Result of web application operation."""
    
    components: Dict[str, ValidationResult]
    configs: Dict[str, Dict[str, Any]]
    stats: Dict[str, Any]
    issues: List[str]


class WebApp:
    """Web application for compound data."""

    def __init__(
        self,
        host: str = "localhost",
        port: int = 8000,
        output_dir: Optional[str] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize web application."""
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        self.host = host
        self.port = port
        self.output_dir = Path(output_dir) if output_dir else None
        self.temp_dir = None
        self.components = {}
        self.configs = {}

    async def start_app(
        self,
        compounds: List[EnrichedData],
    ) -> WebAppResult:
        """Start web application."""
        self.logger.info("Starting web application")
        
        try:
            # Create temporary directory
            self.temp_dir = Path(tempfile.mkdtemp())
            
            # Initialize components
            await self._initialize_components(compounds)
            
            # Start server
            await self._start_server()
            
            # Calculate stats
            stats = self._calculate_stats()
            
            # Save configs
            if self.output_dir:
                self._save_configs()
            
            return WebAppResult(
                is_valid=all(
                    r.is_valid
                    for r in self.components.values()
                ),
                components=self.components,
                configs=self.configs,
                stats=stats,
                issues=self._collect_issues(),
            )
            
        except Exception as e:
            self.logger.error(
                f"Error starting application: {str(e)}",
                exc_info=True
            )
            return WebAppResult(
                is_valid=False,
                components={},
                configs={},
                stats={},
                issues=[str(e)],
            )
        finally:
            # Clean up temporary directory
            if self.temp_dir:
                shutil.rmtree(self.temp_dir)

    async def _initialize_components(
        self,
        compounds: List[EnrichedData],
    ) -> None:
        """Initialize application components."""
        # Initialize analyzers
        analyzer = DataAnalyzer(
            log_level=self.logger.level
        )
        analysis = analyzer.analyze_compounds(compounds)
        self.components["analysis"] = analysis
        
        # Initialize visualizers
        visualizer = DataVisualizer(
            output_dir=str(self.temp_dir / "visualizations"),
            log_level=self.logger.level,
        )
        visualization = visualizer.visualize_analysis(
            compounds,
            analysis,
        )
        self.components["visualization"] = visualization
        
        # Initialize web visualizers
        web_visualizer = WebVisualizer(
            output_dir=str(self.temp_dir / "web_visualizations"),
            log_level=self.logger.level,
        )
        web_visualization = web_visualizer.create_visualizations(
            compounds,
            analysis,
            visualization,
        )
        self.components["web_visualization"] = web_visualization
        
        # Initialize web interface
        web_interface = WebInterface(
            output_dir=str(self.temp_dir / "web_interface"),
            log_level=self.logger.level,
        )
        interface = web_interface.create_interface(
            compounds,
            analysis,
            visualization,
            web_visualization,
        )
        self.components["web_interface"] = interface
        
        # Store component configs
        self.configs = {
            "analysis": {
                "output_dir": str(self.temp_dir / "analysis"),
            },
            "visualization": {
                "output_dir": str(self.temp_dir / "visualizations"),
            },
            "web_visualization": {
                "output_dir": str(self.temp_dir / "web_visualizations"),
            },
            "web_interface": {
                "output_dir": str(self.temp_dir / "web_interface"),
            },
            "web_server": {
                "host": self.host,
                "port": self.port,
                "static_dir": str(self.temp_dir / "static"),
            },
        }

    async def _start_server(self) -> WebServerResult:
        """Start web server."""
        # Initialize server
        server = WebServer(
            host=self.host,
            port=self.port,
            static_dir=str(self.temp_dir / "static"),
            log_level=self.logger.level,
        )
        
        # Start server
        server_result = await server.start_server(
            self.components["compounds"],
            self.components["analysis"],
            self.components["visualization"],
            self.components["web_visualization"],
            self.components["web_interface"],
        )
        
        self.components["web_server"] = server_result
        return server_result

    def _calculate_stats(self) -> Dict[str, Any]:
        """Calculate application statistics."""
        return {
            "components": {
                "total": len(self.components),
                "valid": sum(
                    1 for r in self.components.values()
                    if r.is_valid
                ),
            },
            "configs": {
                "total": len(self.configs),
                "components": len(self.configs.keys()),
            },
            "storage": {
                "temp_size": sum(
                    f.stat().st_size
                    for f in self.temp_dir.rglob("*")
                    if f.is_file()
                ) if self.temp_dir else 0,
                "output_size": sum(
                    f.stat().st_size
                    for f in self.output_dir.rglob("*")
                    if f.is_file()
                ) if self.output_dir else 0,
            },
        }

    def _collect_issues(self) -> List[str]:
        """Collect component issues."""
        issues = []
        
        for name, result in self.components.items():
            if not result.is_valid:
                issues.extend(
                    f"{name}: {issue}"
                    for issue in result.issues
                )
        
        return issues

    def _save_configs(self) -> None:
        """Save application configurations."""
        if not self.output_dir:
            return
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save configs
        config_file = self.output_dir / "app_config.json"
        with open(config_file, "w") as f:
            json.dump(self.configs, f, indent=2)
        
        # Copy static files
        if self.temp_dir:
            static_src = self.temp_dir / "static"
            static_dst = self.output_dir / "static"
            if static_src.exists():
                if static_dst.exists():
                    shutil.rmtree(static_dst)
                shutil.copytree(static_src, static_dst)
