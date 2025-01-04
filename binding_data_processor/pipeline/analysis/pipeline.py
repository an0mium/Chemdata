"""Analysis pipeline.

This module provides the AnalysisPipeline class that:
1. Coordinates analysis components
2. Integrates ML predictions and web data
3. Handles data validation
4. Tracks analysis statistics
5. Generates analysis reports
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from datetime import datetime

from ...models.compound.enhanced import EnhancedCompound
from .base import AnalysisManager, AnalysisConfig


@dataclass
class AnalysisPipelineConfig:
    """Analysis pipeline configuration."""
    
    # Analysis configuration
    analysis_config: AnalysisConfig = field(default_factory=AnalysisConfig)
    
    # Pipeline settings
    batch_size: int = 32
    n_workers: int = 4
    use_cache: bool = True
    
    # Integration settings
    use_ml_predictions: bool = True
    use_web_data: bool = True
    min_confidence: float = 0.5
    
    # Report settings
    generate_reports: bool = True
    report_format: str = "markdown"


@dataclass
class PipelineStats:
    """Analysis pipeline statistics."""
    
    # Processing stats
    total_compounds: int = 0
    processed_compounds: int = 0
    failed_compounds: int = 0
    
    # Analysis stats
    analysis_stats: Dict[str, Any] = field(default_factory=dict)
    
    # Integration stats
    ml_integrations: int = 0
    web_integrations: int = 0
    
    # Report stats
    reports_generated: int = 0
    report_errors: int = 0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary format."""
        return {
            "processing": {
                "total": self.total_compounds,
                "processed": self.processed_compounds,
                "failed": self.failed_compounds,
                "success_rate": self._get_success_rate(),
            },
            "analysis": self.analysis_stats,
            "integration": {
                "ml_integrations": self.ml_integrations,
                "web_integrations": self.web_integrations,
            },
            "reports": {
                "generated": self.reports_generated,
                "errors": self.report_errors,
            },
        }
    
    def _get_success_rate(self) -> Optional[float]:
        """Get processing success rate."""
        if not self.total_compounds:
            return None
        return self.processed_compounds / self.total_compounds


class AnalysisPipeline:
    """Pipeline for compound analysis."""

    def __init__(
        self,
        cache_dir: Optional[Path] = None,
        config: Optional[AnalysisPipelineConfig] = None,
    ):
        """Initialize analysis pipeline.
        
        Args:
            cache_dir: Optional directory for caching
            config: Optional pipeline configuration
        """
        self.logger = logging.getLogger(self.__class__.__name__)
        self.cache_dir = cache_dir
        self.config = config or AnalysisPipelineConfig()
        
        # Initialize analysis manager
        self.analysis = AnalysisManager(
            cache_dir=cache_dir,
            config=self.config.analysis_config,
        )
        
        # Initialize stats
        self.stats = PipelineStats()

    def process_compounds(
        self,
        compounds: List[EnhancedCompound],
        batch_size: Optional[int] = None,
    ) -> List[EnhancedCompound]:
        """Process multiple compounds.
        
        Args:
            compounds: List of compounds to process
            batch_size: Optional batch size override
            
        Returns:
            List of processed compounds
        """
        self.stats.total_compounds += len(compounds)
        batch_size = batch_size or self.config.batch_size
        
        # Process in batches
        for i in range(0, len(compounds), batch_size):
            batch = compounds[i:i + batch_size]
            try:
                self._process_batch(batch)
                self.stats.processed_compounds += len(batch)
            except Exception as e:
                self.logger.error(f"Failed to process batch: {str(e)}")
                self.stats.failed_compounds += len(batch)
        
        return compounds

    def process_compound(
        self,
        compound: EnhancedCompound,
    ) -> EnhancedCompound:
        """Process a single compound.
        
        Args:
            compound: Compound to process
            
        Returns:
            Processed compound
        """
        self.stats.total_compounds += 1
        
        try:
            # Integrate ML predictions
            if self.config.use_ml_predictions:
                self._integrate_ml_predictions(compound)
            
            # Integrate web data
            if self.config.use_web_data:
                self._integrate_web_data(compound)
            
            # Run analysis
            compound = self.analysis.analyze_compound(compound)
            self.stats.analysis_stats = self.analysis.stats.to_dict()
            
            # Generate report
            if self.config.generate_reports:
                self._generate_report(compound)
            
            self.stats.processed_compounds += 1
            return compound
            
        except Exception as e:
            self.logger.error(
                f"Failed to process compound {compound.name}: {str(e)}"
            )
            self.stats.failed_compounds += 1
            raise

    def _process_batch(
        self,
        compounds: List[EnhancedCompound],
    ) -> None:
        """Process a batch of compounds."""
        # Integrate data
        for compound in compounds:
            try:
                # Integrate ML predictions
                if self.config.use_ml_predictions:
                    self._integrate_ml_predictions(compound)
                
                # Integrate web data
                if self.config.use_web_data:
                    self._integrate_web_data(compound)
                    
            except Exception as e:
                self.logger.error(
                    f"Error integrating data for {compound.name}: {str(e)}"
                )
        
        # Run analysis
        for compound in compounds:
            try:
                compound = self.analysis.analyze_compound(compound)
            except Exception as e:
                self.logger.error(
                    f"Error analyzing {compound.name}: {str(e)}"
                )
        
        # Generate reports
        if self.config.generate_reports:
            for compound in compounds:
                try:
                    self._generate_report(compound)
                except Exception as e:
                    self.logger.error(
                        f"Error generating report for {compound.name}: {str(e)}"
                    )

    def _integrate_ml_predictions(
        self,
        compound: EnhancedCompound,
    ) -> None:
        """Integrate ML predictions into analysis."""
        predictions = compound.get_predictions()
        if not predictions:
            return
            
        # Update analysis data
        for pred_type, pred in predictions.items():
            if pred["confidence"] < self.config.min_confidence:
                continue
                
            if pred_type == "binding":
                compound.binding_data.update(pred["data"])
            elif pred_type == "activity":
                compound.activity_data.update(pred["data"])
            elif pred_type == "safety":
                compound.safety_data.update(pred["data"])
            elif pred_type == "properties":
                compound.property_data.update(pred["data"])
        
        self.stats.ml_integrations += 1

    def _integrate_web_data(
        self,
        compound: EnhancedCompound,
    ) -> None:
        """Integrate web data into analysis."""
        # Get web data
        web_data = compound.get_web_data()
        if not web_data:
            return
            
        # Update analysis data
        if "swiss_data" in web_data:
            swiss = web_data["swiss_data"]
            if "targets" in swiss:
                compound.binding_data.update(swiss["targets"])
            if "adme" in swiss:
                compound.property_data.update(swiss["adme"])
        
        if "community_data" in web_data:
            community = web_data["community_data"]
            if "effects" in community:
                compound.activity_data.update(community["effects"])
            if "safety" in community:
                compound.safety_data.update(community["safety"])
        
        if "social_data" in web_data:
            social = web_data["social_data"]
            if "reports" in social:
                compound.activity_data["experience_reports"] = social["reports"]
            if "safety" in social:
                compound.safety_data["community_alerts"] = social["safety"]
        
        self.stats.web_integrations += 1

    def _generate_report(
        self,
        compound: EnhancedCompound,
    ) -> None:
        """Generate analysis report for compound."""
        try:
            # Get analysis results
            results = {
                "binding": compound.binding_data,
                "activity": compound.activity_data,
                "safety": compound.safety_data,
                "properties": compound.property_data,
            }
            
            # Generate report
            if self.config.report_format == "markdown":
                report = self._generate_markdown_report(compound, results)
            else:
                report = self._generate_json_report(compound, results)
            
            # Save report
            report_dir = self.cache_dir / "reports" if self.cache_dir else Path("reports")
            report_dir.mkdir(parents=True, exist_ok=True)
            
            report_path = report_dir / f"{compound.name}.{self.config.report_format}"
            report_path.write_text(report)
            
            self.stats.reports_generated += 1
            
        except Exception as e:
            self.logger.error(
                f"Error generating report for {compound.name}: {str(e)}"
            )
            self.stats.report_errors += 1

    def _generate_markdown_report(
        self,
        compound: EnhancedCompound,
        results: Dict[str, Any],
    ) -> str:
        """Generate markdown format report."""
        lines = [
            f"# Analysis Report: {compound.name}",
            "",
            "## Basic Information",
            f"- Name: {compound.name}",
            f"- SMILES: {compound.smiles}",
            f"- CAS: {compound.cas_number}",
            "",
            "## Binding Analysis",
            self._format_binding_section(results["binding"]),
            "",
            "## Activity Analysis",
            self._format_activity_section(results["activity"]),
            "",
            "## Safety Analysis",
            self._format_safety_section(results["safety"]),
            "",
            "## Property Analysis",
            self._format_property_section(results["properties"]),
            "",
            "## Data Sources",
            "- ML Predictions",
            "- Web Data",
            "- Analysis Results",
            "",
            f"Generated: {datetime.now().isoformat()}",
        ]
        
        return "\n".join(lines)

    def _generate_json_report(
        self,
        compound: EnhancedCompound,
        results: Dict[str, Any],
    ) -> str:
        """Generate JSON format report."""
        import json
        
        report = {
            "name": compound.name,
            "smiles": compound.smiles,
            "cas_number": compound.cas_number,
            "analysis": results,
            "metadata": {
                "generated": datetime.now().isoformat(),
                "sources": [
                    "ML Predictions",
                    "Web Data",
                    "Analysis Results",
                ],
            },
        }
        
        return json.dumps(report, indent=2)

    def _format_binding_section(self, data: Dict[str, Any]) -> str:
        """Format binding analysis section."""
        lines = []
        
        if "strongest_binding" in data:
            binding = data["strongest_binding"]
            lines.extend([
                f"### Strongest Binding",
                f"- Target: {binding['target']}",
                f"- Affinity: {binding['affinity']} {binding['unit']}",
                f"- Activity: {binding['activity']}",
            ])
        
        if "target_selectivity" in data:
            selectivity = data["target_selectivity"]
            lines.extend([
                f"### Target Selectivity",
                f"- Score: {selectivity['score']:.2f}",
                f"- Classification: {selectivity['classification']}",
            ])
        
        return "\n".join(lines)

    def _format_activity_section(self, data: Dict[str, Any]) -> str:
        """Format activity analysis section."""
        lines = []
        
        if "primary_activity" in data:
            activity = data["primary_activity"]
            lines.extend([
                f"### Primary Activity",
                f"- Type: {activity['type']}",
                f"- Score: {activity['score']:.2f}",
                f"- Confidence: {activity['confidence']:.2f}",
            ])
        
        if "mechanisms" in data:
            lines.extend([
                f"### Mechanisms",
                *[f"- {m}" for m in data["mechanisms"]],
            ])
        
        return "\n".join(lines)

    def _format_safety_section(self, data: Dict[str, Any]) -> str:
        """Format safety analysis section."""
        lines = []
        
        if "risk_assessment" in data:
            risk = data["risk_assessment"]
            lines.extend([
                f"### Risk Assessment",
                f"- Level: {risk['level']}",
                f"- Score: {risk['score']:.2f}",
                f"- Confidence: {risk['confidence']:.2f}",
            ])
        
        if "alerts" in data:
            lines.extend([
                f"### Safety Alerts",
                *[f"- {a['type']}: {a['description']}" for a in data["alerts"]],
            ])
        
        return "\n".join(lines)

    def _format_property_section(self, data: Dict[str, Any]) -> str:
        """Format property analysis section."""
        lines = []
        
        if "drug_likeness" in data:
            dl = data["drug_likeness"]
            lines.extend([
                f"### Drug-likeness",
                f"- Score: {dl['score']:.2f}",
                f"- Classification: {dl['classification']}",
            ])
        
        if "properties" in data:
            props = data["properties"]
            lines.extend([
                f"### Properties",
                f"- MW: {props.get('molecular_weight', 'N/A')}",
                f"- LogP: {props.get('logp', 'N/A')}",
                f"- TPSA: {props.get('tpsa', 'N/A')}",
                f"- HBD: {props.get('hbd', 'N/A')}",
                f"- HBA: {props.get('hba', 'N/A')}",
            ])
        
        return "\n".join(lines)
