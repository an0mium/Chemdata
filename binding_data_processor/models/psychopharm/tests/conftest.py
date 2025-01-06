"""Test configuration for psychopharm tests."""

import os
import sys
import pytest

# Add project root to Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../../'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)


# Import fixtures here if needed
@pytest.fixture
def base_compound():
    """Fixture for creating a base compound instance."""
    from binding_data_processor.models.compound.base.core import CompoundData
    return CompoundData(name="Test Compound", smiles="CC(=O)OC1=CC=CC=C1C(=O)O")
