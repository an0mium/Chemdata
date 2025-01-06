"""Export functionality for compound data.

This module extends AnalyzedCompoundData with export capabilities:
- TSV export with flexible column selection
- JSON export with configurable depth
- Report generation
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set
import json
import csv
from pathlib import Path
from datetime import datetime

from .compound_analysis import AnalyzedCompoundData


@dataclass
class ExportableCompoundData(AnalyzedCompoundData):
    """CompoundData with export capabilities."""

    # Export configuration
    _export_config: Dict = field(default_factory=lambda: {
        "tsv_columns": {
            "basic": [
                "name", "smiles", "cas_number", "inchi", "inchi_key",
                "molecular_weight", "logp", "tpsa",
            ],
            "identifiers": [
                "common_name_1", "common_name_2", "common_name_3",
                "chembl_id", "pubchem_cid", "drugbank_id",
            ],
            "binding": [
                "primary_target", "binding_affinities", "binding_types",
                "binding_confidences",
            ],
            "activity": [
                "activity_types", "activity_scores", "activity_mechanisms",
            ],
            "safety": [
                "safety_scores", "safety_warnings", "contraindications",
                "risk_factors",
            ],
            "predictions": [
                "toxicity_predictions", "abuse_potential", "binding_predictions",
                "activity_predictions",
            ],
            "web": [
                "patent_count", "experience_reports", "safety_profile",
                "approval_status",
            ],
        },
        "json_depth": 3,
        "report_sections": [
            "basic_info",
            "binding_analysis",
            "activity_analysis",
            "safety_analysis",
            "sar_analysis",
            "web_data",
            "predictions",
        ],
    })

    def to_tsv_row(self, columns: Optional[List[str]] = None) -> Dict[str, str]:
        """Convert compound data to TSV row format.
        
        Args:
            columns: Optional list of column names to include
            
        Returns:
            Dict mapping column names to string values
        """
        # Get all data
        data = self.to_dict(include_predictions=True)
        row = {}

        # Use specified columns or all basic columns
        columns = columns or self._export_config["tsv_columns"]["basic"]

        # Convert each field to string representation
        for col in columns:
            value = data.get(col)
            if isinstance(value, (dict, list, set)):
                row[col] = json.dumps(value)
            elif isinstance(value, (int, float)):
                row[col] = str(value)
            elif value is None:
                row[col] = ""
            else:
                row[col] = str(value)

        return row

    def to_json(self, max_depth: Optional[int] = None) -> str:
        """Convert compound data to JSON string with configurable depth.
        
        Args:
            max_depth: Maximum depth for nested objects
            
        Returns:
            JSON string representation
        """
        max_depth = max_depth or self._export_config["json_depth"]
        
        def limit_depth(obj, current_depth=0):
            if current_depth >= max_depth:
                if isinstance(obj, (dict, list, set)):
                    return str(obj)
                return obj
            elif isinstance(obj, dict):
                return {
                    k: limit_depth(v, current_depth + 1)
                    for k, v in obj.items()
                }
            elif isinstance(obj, (list, set)):
                return [
                    limit_depth(item, current_depth + 1)
                    for item in obj
                ]
            return obj

        data = self.to_dict(include_predictions=True)
        limited_data = limit_depth(data)
        return json.dumps(limited_data, indent=2)

    def generate_report(
        self,
        sections: Optional[List[str]] = None,
        output_path: Optional[str] = None,
    ) -> str:
        """Generate detailed report in markdown format.
        
        Args:
            sections: Optional list of sections to include
            output_path: Optional path to save report
            
        Returns:
            Report content as string
        """
        sections = sections or self._export_config["report_sections"]
        
        # Build report content
        content = [
            f"# Compound Report: {self.name}",
            f"Generated: {datetime.now().isoformat()}",
            "",
        ]

        # Add sections
        for section in sections:
            if section == "basic_info":
                content.extend(self._format_basic_info())
            elif section == "binding_analysis":
                content.extend(self._format_binding_analysis())
            elif section == "activity_analysis":
                content.extend(self._format_activity_analysis())
            elif section == "safety_analysis":
                content.extend(self._format_safety_analysis())
            elif section == "sar_analysis":
                content.extend(self._format_sar_analysis())
            elif section == "web_data":
                content.extend(self._format_web_data())
            elif section == "predictions":
                content.extend(self._format_predictions())

        # Join content
        report = "\n".join(content)

        # Save if path provided
        if output_path:
            Path(output_path).write_text(report)

        return report

    def _format_basic_info(self) -> List[str]:
        """Format basic compound information section."""
        return [
            "## Basic Information",
            "",
            f"- Name: {self.name}",
            f"- SMILES: {self.smiles}",
            f"- CAS: {self.cas_number or 'N/A'}",
            f"- InChI: {self.inchi or 'N/A'}",
            f"- Molecular Weight: {self.molecular_weight:.2f}",
            f"- LogP: {self.logp:.2f}",
            f"- TPSA: {self.tpsa:.2f}",
            "",
            "### Common Names",
            "",
            f"1. {self.common_name_1} ({self.common_name_1_results} results)",
            f"2. {self.common_name_2} ({self.common_name_2_results} results)",
            f"3. {self.common_name_3} ({self.common_name_3_results} results)",
            "",
        ]

    def _format_binding_analysis(self) -> List[str]:
        """Format binding analysis section."""
        analysis = self.analyze_binding()
        lines = [
            "## Binding Analysis",
            "",
            f"Total Targets: {analysis['total_targets']}",
            "",
        ]

        if analysis["strongest_binding"]:
            sb = analysis["strongest_binding"]
            lines.extend([
                "### Strongest Binding",
                "",
                f"- Target: {sb['target']}",
                f"- Value: {sb['value']:.2f}",
                f"- Type: {sb['type']}",
                f"- Confidence: {sb['confidence']:.2f}",
                "",
            ])

        if analysis["primary_targets"]:
            lines.extend([
                "### Primary Targets",
                "",
                *[f"- {target}" for target in analysis["primary_targets"]],
                "",
            ])

        return lines

    def _format_activity_analysis(self) -> List[str]:
        """Format activity analysis section."""
        analysis = self.analyze_activity()
        lines = [
            "## Activity Analysis",
            "",
        ]

        if analysis["primary_type"]:
            lines.extend([
                f"Primary Type: {analysis['primary_type']}",
                "",
                "Secondary Types:",
                *[f"- {type_}" for type_ in analysis["secondary_types"]],
                "",
            ])

        if analysis["mechanisms"]:
            lines.extend([
                "### Mechanisms",
                "",
                *[f"- {mech}" for mech in analysis["mechanisms"]],
                "",
            ])

        return lines

    def _format_safety_analysis(self) -> List[str]:
        """Format safety analysis section."""
        analysis = self.analyze_safety()
        lines = [
            "## Safety Analysis",
            "",
        ]

        if analysis["warnings"]:
            lines.extend([
                "### Warnings",
                "",
                *[
                    f"- {warning['type']}: {warning['description']} "
                    f"(Severity: {warning['severity']:.2f})"
                    for warning in analysis["warnings"]
                ],
                "",
            ])

        if analysis["contraindications"]:
            lines.extend([
                "### Contraindications",
                "",
                *[f"- {contra}" for contra in analysis["contraindications"]],
                "",
            ])

        return lines

    def _format_sar_analysis(self) -> List[str]:
        """Format SAR analysis section."""
        analysis = self.analyze_sar()
        lines = [
            "## Structure-Activity Analysis",
            "",
        ]

        if analysis["pharmacophores"]:
            lines.extend([
                "### Pharmacophores",
                "",
                *[
                    f"- {pharm['type']}: {pharm['description']} "
                    f"(Score: {pharm['score']:.2f})"
                    for pharm in analysis["pharmacophores"]
                ],
                "",
            ])

        if analysis["mechanism_predictions"]:
            lines.extend([
                "### Predicted Mechanisms",
                "",
                *[
                    f"- {pred['mechanism']}: {pred['probability']:.2f} "
                    f"(Confidence: {pred['confidence']:.2f})"
                    for pred in analysis["mechanism_predictions"]
                ],
                "",
            ])

        return lines

    def _format_web_data(self) -> List[str]:
        """Format web data section."""
        return [
            "## Web Data",
            "",
            f"Patent Count: {self.patent_count}",
            f"Experience Reports: {len(self.experience_reports)}",
            "",
            "### URLs",
            "",
            f"- PubChem: {self.pubchem_url}",
            f"- ChEMBL: {self.chembl_url}",
            f"- Erowid: {self.erowid_url}",
            f"- Wikipedia: {self.wikipedia_url}",
            "",
        ]

    def _format_predictions(self) -> List[str]:
        """Format predictions section."""
        lines = [
            "## ML Predictions",
            "",
        ]

        # Add prediction results
        for pred_type, result in self._prediction_cache.items():
            lines.extend([
                f"### {pred_type}",
                "",
                f"Value: {result.value}",
                f"Confidence: {result.confidence:.2f}",
                "",
            ])

        return lines
