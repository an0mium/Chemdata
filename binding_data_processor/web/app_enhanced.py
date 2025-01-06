"""Enhanced web application for compound data.

This app provides:
1. BindingDB data processing
2. Web data enrichment (community, social, patents, literature)
3. ML predictions (activity, toxicity, abuse potential, BBB)
4. Rich visualization and analysis
5. Flexible export options
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List

import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

from ..models.compound.enhanced import EnhancedCompound
from ..pipeline.processing.pipeline import ProcessingPipeline, ProcessingConfig
from .components.compound_dashboard_enhanced import CompoundDashboardEnhanced
from .components.compound_list_enhanced import CompoundListEnhanced
from .components.compound_detail_enhanced import CompoundDetailEnhanced
from .components.compound_search_enhanced import CompoundSearchEnhanced
from .components.compound_visualization_enhanced import CompoundVisualizationEnhanced
from .components.compound_analysis_enhanced import CompoundAnalysisEnhanced
from .components.compound_export_enhanced import CompoundExportEnhanced


def setup_logging() -> None:
    """Setup logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )


def load_compounds(input_file: Optional[Path] = None) -> List[EnhancedCompound]:
    """Load and process compounds.

    Args:
        input_file: Optional input file path

    Returns:
        List of processed compounds
    """
    # Create pipeline with all features enabled
    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            use_ml_predictions=True,
            use_web_enrichment=True,
            use_social_monitoring=True,
            use_patent_search=True,
            use_literature_search=True,
            use_community_data=True,
            use_swiss_data=True,
            use_enhanced_analysis=True,
        )
    )

    # Process compounds
    compounds = pipeline.process_compounds(
        input_file=input_file,
        output_dir=None,
    )

    return compounds


def render_structure(smiles: str, width: int = 400, height: int = 400) -> str:
    """Render chemical structure as SVG.

    Args:
        smiles: SMILES string
        width: Image width
        height: Image height

    Returns:
        SVG string
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return ""

    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def main():
    """Run enhanced web application."""
    st.set_page_config(
        page_title="ChemData Explorer",
        page_icon="🧪",
        layout="wide",
    )

    # Setup logging
    setup_logging()
    logger = logging.getLogger(__name__)

    try:
        # Initialize enhanced components
        dashboard = CompoundDashboardEnhanced()
        compound_list = CompoundListEnhanced()
        compound_details = CompoundDetailEnhanced()
        compound_search = CompoundSearchEnhanced()
        compound_visualization = CompoundVisualizationEnhanced()
        compound_analysis = CompoundAnalysisEnhanced()
        compound_export = CompoundExportEnhanced()

        # Add sidebar
        with st.sidebar:
            st.title("ChemData Explorer")

            # File upload
            input_file = st.file_uploader(
                "Upload BindingDB file",
                type=["tsv"],
            )

            # Search options
            st.subheader("Search")
            search_type = st.selectbox(
                "Search type",
                ["Text", "Structure", "Property", "Activity", "Safety"],
            )

            if search_type == "Text":
                query = st.text_input("Search query")
                if query:
                    results = compound_search.text_search(query)
                    compound_list.update_compounds(results)

            elif search_type == "Structure":
                smiles = st.text_input("SMILES")
                threshold = st.slider(
                    "Similarity threshold",
                    min_value=0.0,
                    max_value=1.0,
                    value=0.7,
                )
                if smiles:
                    results = compound_search.structure_search(
                        smiles,
                        threshold=threshold,
                    )
                    compound_list.update_compounds(results)

            elif search_type == "Property":
                st.write("Property filters")
                mw_min = st.number_input("Min MW", value=0)
                mw_max = st.number_input("Max MW", value=1000)
                logp_min = st.number_input("Min LogP", value=-5)
                logp_max = st.number_input("Max LogP", value=10)
                psa_min = st.number_input("Min PSA", value=0)
                psa_max = st.number_input("Max PSA", value=200)

                if st.button("Apply property filters"):
                    results = compound_search.property_search(
                        {
                            "molecular_weight": {
                                "min": mw_min,
                                "max": mw_max,
                            },
                            "logp": {
                                "min": logp_min,
                                "max": logp_max,
                            },
                            "psa": {
                                "min": psa_min,
                                "max": psa_max,
                            },
                        }
                    )
                    compound_list.update_compounds(results)

            elif search_type == "Activity":
                st.write("Activity filters")
                target = st.selectbox(
                    "Target",
                    ["5-HT2A", "NMDA", "D2", "SERT", "NET", "DAT"],
                )
                activity_min = st.number_input("Min Activity (nM)", value=0.0)
                activity_max = st.number_input("Max Activity (nM)", value=1000.0)

                if st.button("Apply activity filters"):
                    results = compound_search.activity_search(
                        {
                            "target": target,
                            "activity": {
                                "min": activity_min,
                                "max": activity_max,
                            },
                        }
                    )
                    compound_list.update_compounds(results)

            else:  # Safety search
                st.write("Safety filters")
                bbb = st.checkbox("BBB permeable")
                toxicity = st.selectbox(
                    "Max toxicity",
                    ["Low", "Medium", "High"],
                )
                abuse = st.selectbox(
                    "Max abuse potential",
                    ["Low", "Medium", "High"],
                )

                if st.button("Apply safety filters"):
                    results = compound_search.safety_search(
                        {
                            "bbb_permeable": bbb,
                            "max_toxicity": toxicity.lower(),
                            "max_abuse_potential": abuse.lower(),
                        }
                    )
                    compound_list.update_compounds(results)

            # Analysis options
            st.subheader("Analysis")
            analysis_types = st.multiselect(
                "Analysis types",
                [
                    "Activity",
                    "Structure-Activity",
                    "Safety",
                    "Community",
                    "Literature",
                ],
            )

            # Export options
            st.subheader("Export")
            export_format = st.selectbox(
                "Export format",
                ["TSV", "JSON", "SDF", "PDF"],
            )
            export_columns = st.multiselect(
                "Export columns",
                [
                    "Name",
                    "SMILES",
                    "CAS",
                    "Target",
                    "Activity",
                    "Properties",
                    "Predictions",
                    "Safety",
                    "Community",
                    "Literature",
                ],
                default=["Name", "SMILES", "CAS"],
            )
            if st.button("Export"):
                if len(compound_list.compounds) > 0:
                    output_file = Path(f"export.{export_format.lower()}")
                    compound_export.export_compounds(
                        compound_list.compounds,
                        output_file,
                        format=export_format.lower(),
                        columns=export_columns,
                    )
                    st.success(f"Exported to {output_file}")

        # Main content
        if input_file:
            # Load compounds
            compounds = load_compounds(Path(input_file.name))

            # Update components
            compound_search.update_compounds(compounds)
            compound_list.update_compounds(compounds)

            # Show compound list
            st.subheader("Compounds")

            # Add sorting options
            col1, col2 = st.columns(2)
            with col1:
                sort_by = st.selectbox(
                    "Sort by",
                    [
                        "Name",
                        "Molecular Weight",
                        "LogP",
                        "Activity",
                        "Safety Score",
                        "Community Score",
                    ],
                )
            with col2:
                ascending = st.checkbox("Ascending", value=True)

            compound_list.set_sort(
                sort_by.lower().replace(" ", "_"),
                ascending=ascending,
            )

            # Show current page
            current_page = st.number_input(
                "Page",
                min_value=1,
                max_value=compound_list.get_total_pages(),
                value=1,
            )

            # Display compounds
            compounds = compound_list.get_page(current_page)
            for compound in compounds:
                with st.expander(compound.name):
                    col1, col2 = st.columns(2)

                    # Show structure
                    with col1:
                        svg = render_structure(compound.smiles)
                        st.image(svg)

                    # Show details
                    with col2:
                        compound_details.set_compound(compound)
                        info = compound_details.get_basic_info()

                        st.write("Basic Information:")
                        for key, value in info.items():
                            st.write(f"- {key}: {value}")

                        if compound_details.config.show_predictions:
                            st.write("Predictions:")
                            predictions = compound_details.get_predictions()
                            for key, value in predictions.items():
                                st.write(f"- {key}: {value}")

                    # Show analysis if selected
                    if analysis_types:
                        st.write("Analysis:")
                        analysis_results = compound_analysis.analyze_compound(
                            compound,
                            analysis_types=analysis_types,
                        )
                        for analysis_type, result in analysis_results.items():
                            st.write(f"\n{analysis_type} Analysis:")
                            compound_visualization.render_analysis(
                                result,
                                analysis_type=analysis_type,
                            )

            # Show pagination
            st.write(f"Page {current_page} of {compound_list.get_total_pages()}")

        else:
            st.info("Please upload a BindingDB file to start")

    except Exception as e:
        logger.error(f"Application error: {str(e)}")
        st.error(f"Error: {str(e)}")


if __name__ == "__main__":
    main()
