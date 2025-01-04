#!/usr/bin/env python3
"""Example web application demonstrating compound data interface.

This app shows how to:
1. Load and process compounds
2. Display compound lists with filtering
3. Show compound details with visualizations
4. Enable compound search
5. Export compound data
"""

import logging
from pathlib import Path
from typing import Optional

import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

from binding_data_processor.models.compound.enhanced import EnhancedCompound
from binding_data_processor.pipeline.processing.pipeline import (
    ProcessingPipeline,
    ProcessingConfig,
)
from binding_data_processor.web.components.compound_list import (
    CompoundList,
    CompoundListConfig,
)
from binding_data_processor.web.components.compound_details import (
    CompoundDetails,
    CompoundDetailsConfig,
)
from binding_data_processor.web.components.compound_search import (
    CompoundSearch,
    CompoundSearchConfig,
)


def setup_logging() -> None:
    """Setup logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )


def load_compounds(input_file: Optional[Path] = None) -> list[EnhancedCompound]:
    """Load compounds from input file.
    
    Args:
        input_file: Optional input file path
        
    Returns:
        List of loaded compounds
    """
    # Create pipeline
    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            use_ml_predictions=True,
            use_web_enrichment=True,
            use_social_monitoring=True,
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
    """Run web application."""
    st.set_page_config(
        page_title="ChemData Explorer",
        page_icon="🧪",
        layout="wide",
    )
    
    # Setup logging
    setup_logging()
    logger = logging.getLogger(__name__)
    
    try:
        # Initialize components
        compound_list = CompoundList(
            config=CompoundListConfig(
                page_size=25,
                show_structures=True,
                show_predictions=True,
            )
        )
        
        compound_details = CompoundDetails(
            config=CompoundDetailsConfig(
                show_structure=True,
                show_predictions=True,
                structure_width=400,
                structure_height=400,
            )
        )
        
        compound_search = CompoundSearch(
            config=CompoundSearchConfig(
                enable_text_search=True,
                enable_structure_search=True,
                enable_property_search=True,
            )
        )
        
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
                ["Text", "Structure", "Property"],
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
                    
            else:  # Property search
                st.write("Property filters")
                mw_min = st.number_input("Min MW", value=0)
                mw_max = st.number_input("Max MW", value=1000)
                logp_min = st.number_input("Min LogP", value=-5)
                logp_max = st.number_input("Max LogP", value=10)
                
                if st.button("Apply filters"):
                    results = compound_search.property_search({
                        "molecular_weight": {
                            "min": mw_min,
                            "max": mw_max,
                        },
                        "logp": {
                            "min": logp_min,
                            "max": logp_max,
                        },
                    })
                    compound_list.update_compounds(results)
            
            # Export options
            st.subheader("Export")
            export_format = st.selectbox(
                "Export format",
                ["TSV", "JSON"],
            )
            if st.button("Export"):
                if len(compound_list.compounds) > 0:
                    output_file = Path(f"export.{export_format.lower()}")
                    compound_list.export_compounds(
                        output_file,
                        format=export_format.lower(),
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
                    ["Name", "Molecular Weight", "LogP"],
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
            
            # Show pagination
            st.write(f"Page {current_page} of {compound_list.get_total_pages()}")
            
        else:
            st.info("Please upload a BindingDB file to start")
        
    except Exception as e:
        logger.error(f"Application error: {str(e)}")
        st.error(f"Error: {str(e)}")


if __name__ == "__main__":
    main()
