#!/usr/bin/env python3
"""Script to run the enhanced web application.

This script provides a convenient way to:
1. Start the enhanced web application
2. Configure settings
3. Load data
4. Handle errors

The enhanced app provides:
1. BindingDB data processing
2. Web data enrichment (community, social, patents, literature)
3. ML predictions (activity, toxicity, abuse potential, BBB)
4. Rich visualization and analysis
5. Flexible export options
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
from binding_data_processor.web.app_enhanced import main  # noqa: E402

if __name__ == "__main__":
    # Run enhanced application
    main()
