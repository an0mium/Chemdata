"""
Utility functions for interacting with Cline's tools.
"""


def use_mcp_tool(server_name: str, tool_name: str, arguments: dict) -> dict:
    """
    Use an MCP tool and return its result.

    Args:
        server_name: Name of the MCP server
        tool_name: Name of the tool to use
        arguments: Dictionary of arguments to pass to the tool

    Returns:
        The parsed JSON result from the tool
    """
    import json
    import sys

    # Print the MCP tool command
    print(
        f"""<use_mcp_tool>
<server_name>{server_name}</server_name>
<tool_name>{tool_name}</tool_name>
<arguments>
{json.dumps(arguments)}
</arguments>
</use_mcp_tool>"""
    )
    sys.stdout.flush()

    # Return sample data for testing
    if tool_name == "search_patents":
        return {
            "patents": [
                {
                    "patent_number": "US20200123456",
                    "patent_title": "Novel 5-HT2 receptor ligands and methods of use",
                    "relevance_score": 0.95,
                },
                {
                    "patent_number": "US20200789012",
                    "patent_title": "NMDA receptor antagonists for treating neurological disorders",
                    "relevance_score": 0.92,
                },
            ]
        }
    elif tool_name == "extract_compounds":
        return {
            "compounds": [
                {
                    "name": "Compound A",
                    "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                    "inchi": "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
                    "molecular_formula": "C8H10N4O2",
                    "molecular_weight": 194.19,
                    "activity_type": "Ki",
                    "activity_value": "5.2",
                    "activity_unit": "nM",
                    "mechanism": "5-HT2A antagonist",
                    "primary_target": "5-HT2A receptor",
                    "target_type": "5-HT2",
                    "mentions": [
                        {
                            "context": "Compound A showed high affinity binding (Ki = 5.2 nM) at the 5-HT2A receptor",
                            "has_activity": True,
                            "has_binding": True,
                            "has_affinity": True,
                        }
                    ],
                    "reference_urls": {},
                    "relevance_score": 0.95,
                    "source_patent": "US20200123456",
                    "source_patent_title": "Novel 5-HT2 receptor ligands and methods of use",
                }
            ]
        }
    else:
        return {}  # Default empty response for unknown tools
