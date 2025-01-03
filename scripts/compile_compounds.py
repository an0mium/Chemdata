#!/usr/bin/env python3
"""
Script to search for and compile compound data from patents and other sources.
"""

import json
import pandas as pd
from typing import List, Dict, Any
import asyncio
from datetime import datetime
import os

from binding_data_processor import BindingDataProcessor

from api_client import PubMedClient
from web_enrichment import WebEnrichment


class CompoundCompiler:
    def __init__(self, output_dir: str = "data"):
        """Initialize the compiler with output directory."""
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

        # Initialize clients
        self.pubmed_client = PubMedClient()
        self.web_client = WebEnrichment()
        self.processor = BindingDataProcessor(self.pubmed_client, self.web_client)

        # Define search queries
        self.queries = [
            {
                "name": "5-HT2_ligands",
                "query": "5-HT2 receptor ligand OR serotonin 2 receptor ligand",
                "target_type": "5-HT2",
            },
            {
                "name": "NMDA_antagonists",
                "query": "NMDA receptor antagonist OR N-methyl-D-aspartate receptor antagonist",
                "target_type": "NMDA",
            },
        ]

        # Define source URL fields to include
        self.url_fields = [
            "pubmed_url",
            "chembl_url",
            "pubchem_url",
            "bindingdb_url",
            "patents_url",
            "wikipedia_url",
            "drugbank_url",
            "zinc_url",
            "psychonaut_url",
            "erowid_url",
        ]

    async def search_patents(self, query: str) -> List[Dict[str, Any]]:
        """Search patents using MCP server."""
        try:
            result = await self.web_client.use_mcp_tool(
                "patent-search", "search_patents", {"query": query}
            )
            if isinstance(result, dict) and "patents" in result:
                return result["patents"]
            return []
        except Exception as e:
            print(f"Error searching patents: {str(e)}")
            return []

    async def extract_compounds(self, patent_number: str) -> List[Dict[str, Any]]:
        """Extract compounds from a patent using MCP server."""
        try:
            result = await self.web_client.use_mcp_tool(
                "patent-search", "extract_compounds", {"patent_number": patent_number}
            )
            if isinstance(result, dict) and "compounds" in result:
                return result["compounds"]
            return []
        except Exception as e:
            print(f"Error extracting compounds: {str(e)}")
            return []

    async def compile_compounds(self):
        """Search patents and compile compound data."""
        all_compounds = []

        for search in self.queries:
            print(f"\nSearching for {search['name']}...")

            # Search patents
            patents = await self.search_patents(search["query"])
            print(f"Found {len(patents)} relevant patents")

            # Extract compounds from each patent
            for patent in patents:
                try:
                    compounds = await self.extract_compounds(patent["patent_number"])

                    # Add source information
                    for compound in compounds:
                        compound["source_patent"] = patent["patent_number"]
                        compound["source_patent_title"] = patent["patent_title"]
                        compound["target_type"] = search["target_type"]

                        # Get additional data from web sources
                        urls = self.web_client.get_reference_urls(compound["name"])
                        compound["reference_urls"] = urls

                        all_compounds.append(compound)

                except Exception as e:
                    print(
                        f"Error processing patent {patent['patent_number']}: {str(e)}"
                    )
                    continue

        print(f"\nFound total of {len(all_compounds)} unique compounds")
        return all_compounds

    def create_tsv(self, compounds: List[Dict[str, Any]], filename: str):
        """Create TSV file with compound data."""
        if not compounds:
            # Create empty DataFrame with all expected columns
            columns = [
                "name",
                "smiles",
                "inchi",
                "target_type",
                "source_patent",
                "source_patent_title",
                "relevance_score",
                "reference_urls",
                "mentions",
                "cas_number",
                "molecular_formula",
                "molecular_weight",
                "logp",
                "hbd",
                "hba",
                "tpsa",
                "rotatable_bonds",
                "activity_type",
                "activity_value",
                "activity_unit",
                "mechanism",
                "primary_target",
                "secondary_targets",
                "legal_status",
                "scheduling",
                "controlled_status",
            ] + self.url_fields
            df = pd.DataFrame(columns=columns)
        else:
            # Convert to DataFrame
            df = pd.DataFrame(compounds)

            # Add URL columns
            for url_field in self.url_fields:
                df[url_field] = df["reference_urls"].apply(
                    lambda x: x.get(url_field, "") if isinstance(x, dict) else ""
                )

        # Define column order
        columns = [
            "name",
            "smiles",
            "inchi",
            "target_type",
            "source_patent",
            "source_patent_title",
            "relevance_score",
        ] + self.url_fields

        # Add activity columns
        activity_columns = []
        for mention in df["mentions"].explode():
            if isinstance(mention, dict):
                if mention.get("has_activity"):
                    activity_columns.append(f"Activity: {mention['context'][:100]}...")
                if mention.get("has_binding"):
                    activity_columns.append(f"Binding: {mention['context'][:100]}...")
                if mention.get("has_affinity"):
                    activity_columns.append(f"Affinity: {mention['context'][:100]}...")

        df["activity_data"] = df["mentions"].apply(
            lambda x: (
                "\n".join(
                    f"{m['context']}"
                    for m in x
                    if isinstance(m, dict)
                    and (
                        m.get("has_activity")
                        or m.get("has_binding")
                        or m.get("has_affinity")
                    )
                )
                if isinstance(x, list)
                else ""
            )
        )

        columns.append("activity_data")

        # Save to TSV
        output_path = os.path.join(self.output_dir, filename)
        df[columns].to_csv(output_path, sep="\t", index=False)
        print(f"\nSaved compound data to {output_path}")

        # Return stats
        return {
            "total_compounds": len(df),
            "with_activity": len(df[df["activity_data"].str.len() > 0]),
            "by_target": df["target_type"].value_counts().to_dict(),
        }


async def main():
    compiler = CompoundCompiler()

    # Search and compile compounds
    compounds = await compiler.compile_compounds()

    # Create TSV file
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    stats = compiler.create_tsv(compounds, f"receptor_compounds_{timestamp}.tsv")

    print("\nCompilation complete!")
    print(f"Total compounds: {stats['total_compounds']}")
    print(f"Compounds with activity data: {stats['with_activity']}")
    print("\nCompounds by target type:")
    for target, count in stats["by_target"].items():
        print(f"  {target}: {count}")


if __name__ == "__main__":
    asyncio.run(main())
