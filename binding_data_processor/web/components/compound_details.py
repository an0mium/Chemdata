"""Web component for displaying compound details.

This module provides components for:
1. Displaying detailed compound information
2. Visualizing structure and properties
3. Showing predictions and analysis
4. Displaying web data and references
"""

import logging
from typing import Optional, Dict, Any
from dataclasses import dataclass
from pathlib import Path

from ...models.compound.enhanced import EnhancedCompound


@dataclass
class CompoundDetailsConfig:
    """Configuration for compound details component."""
    
    # Display settings
    show_structure: bool = True
    show_predictions: bool = True
    show_web_data: bool = True
    show_references: bool = True
    
    # Visualization settings
    structure_width: int = 400
    structure_height: int = 400
    plot_width: int = 600
    plot_height: int = 400
    
    # Export settings
    allow_export: bool = True
    export_formats: list = None
    
    def __post_init__(self):
        """Initialize configuration."""
        if self.export_formats is None:
            self.export_formats = ["tsv", "json"]


class CompoundDetails:
    """Component for displaying compound details."""

    def __init__(
        self,
        config: Optional[CompoundDetailsConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize compound details component.
        
        Args:
            config: Optional component configuration
            logger: Optional logger instance
        """
        self.config = config or CompoundDetailsConfig()
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        
        # Initialize state
        self.compound = None
        self.active_tab = "info"
        self.plot_data = {}

    def set_compound(self, compound: EnhancedCompound) -> None:
        """Set compound to display.
        
        Args:
            compound: Compound to display
        """
        self.compound = compound
        self._generate_plots()

    def set_active_tab(self, tab: str) -> None:
        """Set active tab.
        
        Args:
            tab: Tab to activate
        """
        self.active_tab = tab

    def get_basic_info(self) -> Dict[str, Any]:
        """Get basic compound information.
        
        Returns:
            Dictionary of basic information
        """
        if not self.compound:
            return {}
            
        info = {
            "name": self.compound.name,
            "smiles": self.compound.smiles,
            "source": self.compound.source,
        }
        
        # Add CAS number if available
        if hasattr(self.compound, "cas_number"):
            info["cas_number"] = self.compound.cas_number
            
        # Add molecular formula if available
        if hasattr(self.compound, "molecular_formula"):
            info["molecular_formula"] = self.compound.molecular_formula
            
        # Add molecular weight if available
        if hasattr(self.compound, "molecular_weight"):
            info["molecular_weight"] = self.compound.molecular_weight
            
        return info

    def get_predictions(self) -> Dict[str, Any]:
        """Get compound predictions.
        
        Returns:
            Dictionary of predictions
        """
        if not self.compound or not self.config.show_predictions:
            return {}
            
        predictions = {}
        
        # Add activity predictions
        if hasattr(self.compound, "activity_predictions"):
            predictions["activity"] = self.compound.activity_predictions
            
        # Add toxicity predictions
        if hasattr(self.compound, "toxicity_predictions"):
            predictions["toxicity"] = self.compound.toxicity_predictions
            
        # Add abuse predictions
        if hasattr(self.compound, "abuse_predictions"):
            predictions["abuse"] = self.compound.abuse_predictions
            
        # Add BBB predictions
        if hasattr(self.compound, "bbb_predictions"):
            predictions["bbb"] = self.compound.bbb_predictions
            
        return predictions

    def get_web_data(self) -> Dict[str, Any]:
        """Get compound web data.
        
        Returns:
            Dictionary of web data
        """
        if not self.compound or not self.config.show_web_data:
            return {}
            
        web_data = {}
        
        # Add community data
        if hasattr(self.compound, "community_data"):
            web_data["community"] = self.compound.community_data
            
        # Add social data
        if hasattr(self.compound, "social_data"):
            web_data["social"] = self.compound.social_data
            
        # Add Swiss data
        if hasattr(self.compound, "swiss_data"):
            web_data["swiss"] = self.compound.swiss_data
            
        return web_data

    def get_references(self) -> Dict[str, Any]:
        """Get compound references.
        
        Returns:
            Dictionary of references
        """
        if not self.compound or not self.config.show_references:
            return {}
            
        references = {}
        
        # Add literature references
        if hasattr(self.compound, "literature_references"):
            references["literature"] = self.compound.literature_references
            
        # Add patent references
        if hasattr(self.compound, "patent_references"):
            references["patents"] = self.compound.patent_references
            
        # Add web references
        if hasattr(self.compound, "web_references"):
            references["web"] = self.compound.web_references
            
        return references

    def get_plot_data(self) -> Dict[str, Any]:
        """Get plot data.
        
        Returns:
            Dictionary of plot data
        """
        return self.plot_data

    def export_details(
        self,
        output_file: Path,
        format: str = "json",
    ) -> None:
        """Export compound details.
        
        Args:
            output_file: Output file path
            format: Export format (json/tsv)
        """
        if not self.compound or not self.config.allow_export:
            return
            
        try:
            # Create output directory
            output_file.parent.mkdir(parents=True, exist_ok=True)
            
            # Get all data
            data = {
                "info": self.get_basic_info(),
                "predictions": self.get_predictions(),
                "web_data": self.get_web_data(),
                "references": self.get_references(),
            }
            
            # Export in requested format
            if format == "json":
                self._export_json(data, output_file)
            else:
                self._export_tsv(data, output_file)
                
            self.logger.info(f"Exported details to {output_file}")
            
        except Exception as e:
            self.logger.error(f"Error exporting details: {str(e)}")
            raise

    def _generate_plots(self) -> None:
        """Generate plots for compound."""
        if not self.compound:
            return
            
        try:
            plots = {}
            
            # Generate structure plot if enabled
            if self.config.show_structure:
                plots["structure"] = self._generate_structure_plot()
            
            # Generate prediction plots if enabled
            if self.config.show_predictions:
                if hasattr(self.compound, "activity_predictions"):
                    plots["activity"] = self._generate_activity_plot()
                if hasattr(self.compound, "toxicity_predictions"):
                    plots["toxicity"] = self._generate_toxicity_plot()
                if hasattr(self.compound, "abuse_predictions"):
                    plots["abuse"] = self._generate_abuse_plot()
            
            self.plot_data = plots
            
        except Exception as e:
            self.logger.error(f"Error generating plots: {str(e)}")
            self.plot_data = {}

    def _generate_structure_plot(self) -> Dict[str, Any]:
        """Generate structure plot.
        
        Returns:
            Plot data dictionary
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw
            
            # Create molecule
            mol = Chem.MolFromSmiles(self.compound.smiles)
            if not mol:
                return {}
            
            # Generate SVG
            drawer = Draw.rdMolDraw2D.MolDraw2DSVG(
                self.config.structure_width,
                self.config.structure_height,
            )
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            
            return {
                "type": "svg",
                "data": svg,
                "width": self.config.structure_width,
                "height": self.config.structure_height,
            }
            
        except Exception as e:
            self.logger.error(f"Error generating structure plot: {str(e)}")
            return {}

    def _generate_activity_plot(self) -> Dict[str, Any]:
        """Generate activity plot.
        
        Returns:
            Plot data dictionary
        """
        try:
            import plotly.graph_objects as go
            
            # Get activity data
            activities = self.compound.activity_predictions
            if not activities:
                return {}
            
            # Create plot
            fig = go.Figure()
            
            # Add bars
            fig.add_trace(go.Bar(
                x=list(activities.keys()),
                y=list(activities.values()),
                name="Activity",
            ))
            
            # Update layout
            fig.update_layout(
                title="Predicted Activities",
                xaxis_title="Activity Type",
                yaxis_title="Score",
                width=self.config.plot_width,
                height=self.config.plot_height,
            )
            
            return {
                "type": "plotly",
                "data": fig.to_dict(),
                "width": self.config.plot_width,
                "height": self.config.plot_height,
            }
            
        except Exception as e:
            self.logger.error(f"Error generating activity plot: {str(e)}")
            return {}

    def _generate_toxicity_plot(self) -> Dict[str, Any]:
        """Generate toxicity plot.
        
        Returns:
            Plot data dictionary
        """
        try:
            import plotly.graph_objects as go
            
            # Get toxicity data
            toxicity = self.compound.toxicity_predictions
            if not toxicity:
                return {}
            
            # Create plot
            fig = go.Figure()
            
            # Add indicators
            fig.add_trace(go.Indicator(
                mode="gauge+number",
                value=toxicity.get("score", 0),
                title={"text": "Toxicity Score"},
                gauge={
                    "axis": {"range": [0, 1]},
                    "bar": {"color": "red"},
                },
            ))
            
            # Update layout
            fig.update_layout(
                title="Predicted Toxicity",
                width=self.config.plot_width,
                height=self.config.plot_height,
            )
            
            return {
                "type": "plotly",
                "data": fig.to_dict(),
                "width": self.config.plot_width,
                "height": self.config.plot_height,
            }
            
        except Exception as e:
            self.logger.error(f"Error generating toxicity plot: {str(e)}")
            return {}

    def _generate_abuse_plot(self) -> Dict[str, Any]:
        """Generate abuse potential plot.
        
        Returns:
            Plot data dictionary
        """
        try:
            import plotly.graph_objects as go
            
            # Get abuse data
            abuse = self.compound.abuse_predictions
            if not abuse:
                return {}
            
            # Create plot
            fig = go.Figure()
            
            # Add indicators
            fig.add_trace(go.Indicator(
                mode="gauge+number",
                value=abuse.get("score", 0),
                title={"text": "Abuse Potential"},
                gauge={
                    "axis": {"range": [0, 1]},
                    "bar": {"color": "orange"},
                },
            ))
            
            # Update layout
            fig.update_layout(
                title="Predicted Abuse Potential",
                width=self.config.plot_width,
                height=self.config.plot_height,
            )
            
            return {
                "type": "plotly",
                "data": fig.to_dict(),
                "width": self.config.plot_width,
                "height": self.config.plot_height,
            }
            
        except Exception as e:
            self.logger.error(f"Error generating abuse plot: {str(e)}")
            return {}

    def _export_json(self, data: Dict[str, Any], output_file: Path) -> None:
        """Export data as JSON.
        
        Args:
            data: Data to export
            output_file: Output file path
        """
        import json
        
        with open(output_file, "w") as f:
            json.dump(data, f, indent=2)

    def _export_tsv(self, data: Dict[str, Any], output_file: Path) -> None:
        """Export data as TSV.
        
        Args:
            data: Data to export
            output_file: Output file path
        """
        import pandas as pd
        
        # Flatten data
        flat_data = {}
        for section, values in data.items():
            for key, value in values.items():
                flat_data[f"{section}_{key}"] = str(value)
        
        # Create DataFrame
        df = pd.DataFrame([flat_data])
        
        # Save TSV
        df.to_csv(output_file, sep="\t", index=False)
