"""Enhanced compound list component.

This module provides an enhanced web interface for displaying lists of compounds with:
1. Advanced filtering and sorting
2. Pagination and batch operations
3. Interactive visualizations
4. Export functionality
5. Comprehensive statistics tracking
"""

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Callable

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors, Crippen

from ...models.compound.enhanced import EnhancedCompound


@dataclass
class ListConfig:
    """Configuration for list component."""

    # Display options
    page_size: int = 10
    show_structures: bool = True
    show_predictions: bool = True
    show_analysis: bool = True
    show_visualizations: bool = True

    # Sort options
    default_sort: str = "name"
    default_ascending: bool = True

    # Filter options
    enable_filtering: bool = True
    filter_columns: Optional[List[str]] = None

    # Selection options
    enable_selection: bool = True
    enable_batch_ops: bool = True

    # Export options
    export_formats: List[str] = field(default_factory=lambda: ["tsv", "csv", "json", "sdf"])
    default_export_columns: List[str] = field(
        default_factory=lambda: [
            "name",
            "smiles",
            "cas",
            "binding_affinities",
            "toxicity_score",
            "abuse_potential",
            "bbb_permeability",
        ]
    )

    # Visualization options
    chart_height: int = 300
    chart_width: int = 400
    structure_width: int = 200
    structure_height: int = 200


class CompoundListEnhanced:
    """Enhanced compound list component."""

    def __init__(
        self,
        config: Optional[ListConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize list.

        Args:
            config: Optional list configuration
            logger: Optional logger instance
        """
        self.config = config or ListConfig()
        self.logger = logger or logging.getLogger(__name__)

        # Initialize state
        if "compounds" not in st.session_state:
            st.session_state.compounds = []
        if "filtered_compounds" not in st.session_state:
            st.session_state.filtered_compounds = []
        if "selected_compounds" not in st.session_state:
            st.session_state.selected_compounds = set()
        if "sort_by" not in st.session_state:
            st.session_state.sort_by = self.config.default_sort
        if "sort_ascending" not in st.session_state:
            st.session_state.sort_ascending = self.config.default_ascending
        if "current_page" not in st.session_state:
            st.session_state.current_page = 1
        if "filters" not in st.session_state:
            st.session_state.filters = {}
        if "view_stats" not in st.session_state:
            st.session_state.view_stats = {
                "total_views": 0,
                "filtered_views": 0,
                "exported_views": 0,
                "visualization_views": 0,
                "last_update": datetime.now(),
            }

    def update_compounds(self, compounds: List[EnhancedCompound]) -> None:
        """Update compound list.

        Args:
            compounds: List of compounds
        """
        try:
            st.session_state.compounds = compounds
            self._apply_filters()
            self._apply_sort()

            # Update stats
            st.session_state.view_stats["total_views"] += 1
            st.session_state.view_stats["last_update"] = datetime.now()

            self.logger.info(
                f"Updated compounds: {len(compounds)} total, " f"{len(st.session_state.filtered_compounds)} filtered"
            )

        except Exception as e:
            self.logger.error(f"Error updating compounds: {str(e)}")
            st.error(f"Error updating compounds: {str(e)}")

    def set_sort(self, sort_by: str, ascending: bool = True) -> None:
        """Set sort order.

        Args:
            sort_by: Column to sort by
            ascending: Sort ascending if True
        """
        try:
            st.session_state.sort_by = sort_by
            st.session_state.sort_ascending = ascending
            self._apply_sort()

            self.logger.info(f"Applied sort: {sort_by} ({'ascending' if ascending else 'descending'})")

        except Exception as e:
            self.logger.error(f"Error setting sort: {str(e)}")
            st.error(f"Error setting sort: {str(e)}")

    def set_filters(self, filters: Dict[str, Any]) -> None:
        """Set filters.

        Args:
            filters: Dictionary of column filters
        """
        try:
            st.session_state.filters = filters
            self._apply_filters()
            self._apply_sort()

            # Update stats
            st.session_state.view_stats["filtered_views"] += 1

            self.logger.info(
                f"Applied filters: {list(filters.keys())}, "
                f"matched {len(st.session_state.filtered_compounds)} compounds"
            )

        except Exception as e:
            self.logger.error(f"Error setting filters: {str(e)}")
            st.error(f"Error setting filters: {str(e)}")

    def get_page(self, page: int) -> List[EnhancedCompound]:
        """Get compounds for specified page.

        Args:
            page: Page number (1-based)

        Returns:
            List of compounds for page
        """
        if not st.session_state.filtered_compounds:
            return []

        start = (page - 1) * self.config.page_size
        end = start + self.config.page_size
        return st.session_state.filtered_compounds[start:end]

    def get_total_pages(self) -> int:
        """Get total number of pages.

        Returns:
            Total pages
        """
        if not st.session_state.filtered_compounds:
            return 1

        return (len(st.session_state.filtered_compounds) + self.config.page_size - 1) // self.config.page_size

    def get_selected(self) -> List[EnhancedCompound]:
        """Get selected compounds.

        Returns:
            List of selected compounds
        """
        return [c for c in st.session_state.filtered_compounds if c.name in st.session_state.selected_compounds]

    def export_compounds(
        self,
        compounds: Optional[List[EnhancedCompound]] = None,
        format: str = "tsv",
        columns: Optional[List[str]] = None,
        output_file: Optional[Path] = None,
    ) -> Optional[str]:
        """Export compounds to file or string.

        Args:
            compounds: Optional list of compounds to export (defaults to filtered compounds)
            format: Export format (tsv/csv/json/sdf)
            columns: Optional list of columns to export
            output_file: Optional output file path

        Returns:
            Export data as string if no output file provided
        """
        try:
            compounds = compounds or st.session_state.filtered_compounds
            columns = columns or self.config.default_export_columns

            # Convert to DataFrame
            data = []
            for compound in compounds:
                row = {
                    "name": compound.name,
                    "smiles": compound.smiles,
                    "cas": compound.cas,
                    "compound_type": compound.compound_type.value,
                    "legal_status": compound.legal_status.value,
                }

                # Add predictions
                if "binding_affinities" in columns:
                    row["binding_affinities"] = compound.predictions.binding_affinities
                if "toxicity_score" in columns:
                    row["toxicity_score"] = compound.predictions.toxicity_score
                if "abuse_potential" in columns:
                    row["abuse_potential"] = compound.predictions.abuse_potential
                if "bbb_permeability" in columns:
                    row["bbb_permeability"] = compound.predictions.bbb_permeability.value

                # Add web data
                if "community_reports" in columns:
                    row["community_reports"] = compound.web_data.community_reports
                if "pubmed_articles" in columns:
                    row["pubmed_articles"] = compound.web_data.pubmed_articles
                if "patents" in columns:
                    row["patents"] = compound.web_data.patents

                # Add analysis
                if "structural_alerts" in columns:
                    row["structural_alerts"] = compound.analysis.structural_alerts
                if "similar_compounds" in columns:
                    row["similar_compounds"] = compound.analysis.similar_compounds
                if "target_interactions" in columns:
                    row["target_interactions"] = compound.analysis.target_interactions

                data.append(row)

            df = pd.DataFrame(data)

            # Filter columns
            if columns:
                df = df[[c for c in columns if c in df.columns]]

            # Export
            if format == "tsv":
                export_data = df.to_csv(sep="\t", index=False)
            elif format == "csv":
                export_data = df.to_csv(index=False)
            elif format == "json":
                export_data = df.to_json(orient="records", indent=2)
            elif format == "sdf":
                # Create SDF file
                writer = Chem.SDWriter(str(output_file) if output_file else None)
                for compound in compounds:
                    if not compound.mol:
                        continue
                    # Add properties
                    for key, value in df.loc[df["name"] == compound.name].iloc[0].items():
                        if isinstance(value, (str, int, float)):
                            compound.mol.SetProp(key, str(value))
                    writer.write(compound.mol)
                writer.close()
                export_data = None
            else:
                raise ValueError(f"Unsupported format: {format}")

            # Write to file or return string
            if output_file and export_data:
                with open(output_file, "w") as f:
                    f.write(export_data)

            # Update stats
            st.session_state.view_stats["exported_views"] += 1

            self.logger.info(
                f"Exported {len(compounds)} compounds to {format} format " f"with {len(df.columns)} columns"
            )

            return export_data

        except Exception as e:
            self.logger.error(f"Error exporting compounds: {str(e)}")
            st.error(f"Error exporting compounds: {str(e)}")
            return None

    def _apply_filters(self) -> None:
        """Apply filters to compounds."""
        filtered = st.session_state.compounds

        for column, filter_value in st.session_state.filters.items():
            if not filter_value:
                continue

            if column == "name":
                filtered = [c for c in filtered if filter_value.lower() in c.name.lower()]

            elif column == "smiles":
                filtered = [c for c in filtered if filter_value.lower() in c.smiles.lower()]

            elif column == "cas":
                filtered = [c for c in filtered if filter_value in c.cas]

            elif column == "molecular_weight":
                min_mw = filter_value.get("min")
                max_mw = filter_value.get("max")
                if min_mw is not None:
                    filtered = [c for c in filtered if c.mol and Chem.Descriptors.ExactMolWt(c.mol) >= min_mw]
                if max_mw is not None:
                    filtered = [c for c in filtered if c.mol and Chem.Descriptors.ExactMolWt(c.mol) <= max_mw]

            elif column == "logp":
                min_logp = filter_value.get("min")
                max_logp = filter_value.get("max")
                if min_logp is not None:
                    filtered = [c for c in filtered if c.mol and Chem.Crippen.MolLogP(c.mol) >= min_logp]
                if max_logp is not None:
                    filtered = [c for c in filtered if c.mol and Chem.Crippen.MolLogP(c.mol) <= max_logp]

            elif column == "binding_affinity":
                target = filter_value.get("target")
                min_value = filter_value.get("min")
                max_value = filter_value.get("max")
                if target:
                    filtered = [
                        c
                        for c in filtered
                        if target in c.predictions.binding_affinities
                        and (min_value is None or c.predictions.binding_affinities[target] >= min_value)
                        and (max_value is None or c.predictions.binding_affinities[target] <= max_value)
                    ]

            elif column == "toxicity":
                max_score = filter_value.get("max")
                if max_score is not None:
                    filtered = [c for c in filtered if c.predictions.toxicity_score <= max_score]

            elif column == "abuse_potential":
                max_score = filter_value.get("max")
                if max_score is not None:
                    filtered = [c for c in filtered if c.predictions.abuse_potential <= max_score]

            elif column == "bbb_permeability":
                if filter_value:
                    filtered = [c for c in filtered if c.predictions.bbb_permeability.value == filter_value]

        st.session_state.filtered_compounds = filtered

    def _apply_sort(self) -> None:
        """Apply sort to filtered compounds."""
        if not st.session_state.filtered_compounds:
            return

        sort_by = st.session_state.sort_by
        ascending = st.session_state.sort_ascending

        if sort_by == "name":
            st.session_state.filtered_compounds.sort(
                key=lambda x: x.name.lower(),
                reverse=not ascending,
            )

        elif sort_by == "molecular_weight":
            st.session_state.filtered_compounds.sort(
                key=lambda x: Chem.Descriptors.ExactMolWt(x.mol) if x.mol else 0,
                reverse=not ascending,
            )

        elif sort_by == "logp":
            st.session_state.filtered_compounds.sort(
                key=lambda x: Chem.Crippen.MolLogP(x.mol) if x.mol else 0,
                reverse=not ascending,
            )

        elif sort_by == "toxicity_score":
            st.session_state.filtered_compounds.sort(
                key=lambda x: x.predictions.toxicity_score,
                reverse=not ascending,
            )

        elif sort_by == "abuse_potential":
            st.session_state.filtered_compounds.sort(
                key=lambda x: x.predictions.abuse_potential,
                reverse=not ascending,
            )

        elif sort_by == "quality_score":
            st.session_state.filtered_compounds.sort(
                key=lambda x: self._calculate_quality_score(x),
                reverse=not ascending,
            )

    def _calculate_quality_score(self, compound: EnhancedCompound) -> float:
        """Calculate quality score for compound.

        Args:
            compound: Compound to score

        Returns:
            Quality score (0-1)
        """
        score = 0.0
        total = 0.0

        # Check required fields
        if compound.name:
            score += 1
            total += 1
        if compound.smiles:
            score += 1
            total += 1
        if compound.cas:
            score += 1
            total += 1

        # Check predictions
        if compound.predictions.binding_affinities:
            score += 1
            total += 1
        if compound.predictions.toxicity_score is not None:
            score += 1
            total += 1
        if compound.predictions.abuse_potential is not None:
            score += 1
            total += 1
        if compound.predictions.bbb_permeability is not None:
            score += 1
            total += 1

        # Check web data
        if compound.web_data.community_reports:
            score += 1
            total += 1
        if compound.web_data.pubmed_articles:
            score += 1
            total += 1
        if compound.web_data.patents:
            score += 1
            total += 1

        # Check analysis
        if compound.analysis.structural_alerts:
            score += 1
            total += 1
        if compound.analysis.similar_compounds:
            score += 1
            total += 1
        if compound.analysis.target_interactions:
            score += 1
            total += 1

        return score / total if total > 0 else 0.0

    def _generate_visualizations(self) -> Dict[str, Any]:
        """Generate visualizations for current compounds.

        Returns:
            Dictionary of visualization data
        """
        if not self.config.show_visualizations:
            return {}

        visualizations = {}
        compounds = st.session_state.filtered_compounds

        # Binding distribution
        binding_data = []
        for compound in compounds:
            if compound.predictions.binding_affinities:
                for target, value in compound.predictions.binding_affinities.items():
                    binding_data.append(
                        {
                            "compound": compound.name,
                            "target": target,
                            "affinity": value,
                        }
                    )

        if binding_data:
            df = pd.DataFrame(binding_data)
            fig = go.Figure(
                data=[
                    go.Box(
                        y=df["affinity"],
                        x=df["target"],
                        name="Binding Affinity",
                    )
                ]
            )
            fig.update_layout(
                title="Binding Affinity Distribution by Target",
                yaxis_title="Affinity (nM)",
                boxmode="group",
                height=self.config.chart_height,
                width=self.config.chart_width,
            )
            visualizations["binding_distribution"] = fig

        # Safety profile
        safety_data = []
        for compound in compounds:
            safety_data.append(
                {
                    "compound": compound.name,
                    "toxicity_score": compound.predictions.toxicity_score,
                    "abuse_potential": compound.predictions.abuse_potential,
                }
            )

        if safety_data:
            df = pd.DataFrame(safety_data)
            fig = go.Figure(
                data=[
                    go.Scatter(
                        x=df["toxicity_score"],
                        y=df["abuse_potential"],
                        mode="markers+text",
                        text=df["compound"],
                        textposition="top center",
                    )
                ]
            )
            fig.update_layout(
                title="Safety Profile Distribution",
                xaxis_title="Toxicity Score",
                yaxis_title="Abuse Potential",
                height=self.config.chart_height,
                width=self.config.chart_width,
            )
            visualizations["safety_profile"] = fig

        # Update stats
        if visualizations:
            st.session_state.view_stats["visualization_views"] += 1

        return visualizations

    def render_structure(self, smiles: str) -> str:
        """Render chemical structure as SVG.

        Args:
            smiles: SMILES string

        Returns:
            SVG string
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return ""

        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(
            self.config.structure_width,
            self.config.structure_height,
        )
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    def render_predictions(self, compound: EnhancedCompound) -> None:
        """Render predictions section.

        Args:
            compound: Compound to render predictions for
        """
        if not self.config.show_predictions:
            return

        st.write("Predictions:")

        # Binding predictions
        if compound.predictions.binding_affinities:
            st.write("Binding Affinities:")
            for target, value in compound.predictions.binding_affinities.items():
                st.write(f"- {target}: {value:.1f} nM")

        # Safety predictions
        st.write("Safety Profile:")
        st.write(f"- Toxicity Score: {compound.predictions.toxicity_score:.2f}")
        st.write(f"- Abuse Potential: {compound.predictions.abuse_potential:.2f}")
        st.write(f"- BBB Permeability: {compound.predictions.bbb_permeability.value}")

    def render_analysis(self, compound: EnhancedCompound) -> None:
        """Render analysis section.

        Args:
            compound: Compound to render analysis for
        """
        if not self.config.show_analysis:
            return

        st.write("Analysis:")

        # Structure analysis
        if compound.analysis.structural_alerts:
            st.write("Structural Alerts:")
            for alert in compound.analysis.structural_alerts:
                st.write(f"- {alert}")

        # Similar compounds
        if compound.analysis.similar_compounds:
            st.write("Similar Compounds:")
            for similar in compound.analysis.similar_compounds[:3]:
                st.write(f"- {similar}")

        # Target interactions
        if compound.analysis.target_interactions:
            st.write("Target Interactions:")
            for target, value in compound.analysis.target_interactions.items():
                st.write(f"- {target}: {value:.2f}")

    def render(self) -> None:
        """Render compound list."""
        # Get current page compounds
        compounds = self.get_page(st.session_state.current_page)
        if not compounds:
            st.info("No compounds found")
            return

        # Generate visualizations
        if self.config.show_visualizations:
            visualizations = self._generate_visualizations()
            for name, fig in visualizations.items():
                st.plotly_chart(fig, use_container_width=True)

        # Show compounds
        for compound in compounds:
            with st.expander(compound.name):
                col1, col2 = st.columns(2)

                # Show structure
                with col1:
                    if self.config.show_structures and compound.smiles:
                        svg = self.render_structure(compound.smiles)
                        st.image(svg)

                # Show details
                with col2:
                    st.write("Basic Information:")
                    st.write(f"- Name: {compound.name}")
                    st.write(f"- SMILES: {compound.smiles}")
                    st.write(f"- CAS: {compound.cas}")
                    st.write(f"- Type: {compound.compound_type.value}")
                    st.write(f"- Legal Status: {compound.legal_status.value}")

                    # Show predictions
                    self.render_predictions(compound)

                    # Show analysis
                    self.render_analysis(compound)

                # Selection checkbox
                if self.config.enable_selection:
                    selected = st.checkbox(
                        "Select",
                        key=f"select_{compound.name}",
                        value=compound.name in st.session_state.selected_compounds,
                    )
                    if selected:
                        st.session_state.selected_compounds.add(compound.name)
                    else:
                        st.session_state.selected_compounds.discard(compound.name)

        # Show pagination
        st.write(f"Page {st.session_state.current_page} of {self.get_total_pages()}")
        col1, col2 = st.columns(2)
        with col1:
            if st.session_state.current_page > 1:
                if st.button("Previous Page"):
                    st.session_state.current_page -= 1
        with col2:
            if st.session_state.current_page < self.get_total_pages():
                if st.button("Next Page"):
                    st.session_state.current_page += 1

        # Show batch operations
        if self.config.enable_batch_ops and st.session_state.selected_compounds:
            st.subheader("Batch Operations")
            st.write(f"Selected: {len(st.session_state.selected_compounds)} compounds")

            col1, col2 = st.columns(2)
            with col1:
                if st.button("Export Selected"):
                    selected = self.get_selected()
                    export_data = self.export_compounds(
                        compounds=selected,
                        format="tsv",
                    )
                    if export_data:
                        st.download_button(
                            "Download TSV",
                            export_data,
                            file_name="selected_compounds.tsv",
                            mime="text/tab-separated-values",
                        )

            with col2:
                if st.button("Clear Selection"):
                    st.session_state.selected_compounds.clear()

        # Show stats
        if st.session_state.view_stats:
            st.subheader("View Statistics")
            col1, col2 = st.columns(2)
            with col1:
                st.write("Total Views:", st.session_state.view_stats["total_views"])
                st.write("Filtered Views:", st.session_state.view_stats["filtered_views"])
            with col2:
                st.write("Exported Views:", st.session_state.view_stats["exported_views"])
                st.write("Visualization Views:", st.session_state.view_stats["visualization_views"])
            st.write("Last Update:", st.session_state.view_stats["last_update"])
