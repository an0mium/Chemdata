#!/usr/bin/env python3
"""Script to run the web application.

This script provides a convenient way to:
1. Start the web application
2. Configure settings
3. Load data
4. Handle errors
"""

import os
import sys
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))
os.environ.setdefault("PYTHONPATH", str(project_root))
os.environ.setdefault("PYTHONUNBUFFERED", "1")

# Import after path setup
from binding_data_processor.processors.psychopharm.predictors.cli import main  # noqa: E402

if __name__ == "__main__":
    # Run application
    main()
