#!/bin/bash

# Setup script for model consolidation

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}Starting model consolidation setup...${NC}"

# Create directory structure
echo -e "\n${BLUE}Creating directory structure...${NC}"
mkdir -p binding_data_processor/models/compound/{base,ml,enrichment,analysis,export}
mkdir -p binding_data_processor/models/compound/{base,ml,enrichment,analysis,export}/tests

# Create __init__.py files
echo -e "\n${BLUE}Creating __init__.py files...${NC}"
for dir in binding_data_processor/models/compound/{base,ml,enrichment,analysis,export}; do
    touch "$dir/__init__.py"
    echo "Created $dir/__init__.py"
done

# Move base files
echo -e "\n${BLUE}Moving base files...${NC}"
if [ -f binding_data_processor/models/compound.py ]; then
    mv binding_data_processor/models/compound.py binding_data_processor/models/compound/base/core.py
    echo "Moved compound.py to base/core.py"
fi

if [ -f binding_data_processor/models/mixins.py ]; then
    mv binding_data_processor/models/mixins.py binding_data_processor/models/compound/base/mixins.py
    echo "Moved mixins.py to base/mixins.py"
fi

if [ -f binding_data_processor/models/types.py ]; then
    mv binding_data_processor/models/types.py binding_data_processor/models/compound/base/types.py
    echo "Moved types.py to base/types.py"
fi

# Move ML files
echo -e "\n${BLUE}Moving ML files...${NC}"
if [ -f binding_data_processor/models/compound_ml.py ]; then
    mv binding_data_processor/models/compound_ml.py binding_data_processor/models/compound/ml/predictors.py
    echo "Moved compound_ml.py to ml/predictors.py"
fi

# Move enrichment files
echo -e "\n${BLUE}Moving enrichment files...${NC}"
if [ -f binding_data_processor/models/compound_enrichment.py ]; then
    mv binding_data_processor/models/compound_enrichment.py binding_data_processor/models/compound/enrichment/web.py
    echo "Moved compound_enrichment.py to enrichment/web.py"
fi

# Move analysis files
echo -e "\n${BLUE}Moving analysis files...${NC}"
if [ -f binding_data_processor/models/compound_analysis.py ]; then
    mv binding_data_processor/models/compound_analysis.py binding_data_processor/models/compound/analysis/base.py
    echo "Moved compound_analysis.py to analysis/base.py"
fi

# Move export files
echo -e "\n${BLUE}Moving export files...${NC}"
if [ -f binding_data_processor/models/compound_export.py ]; then
    mv binding_data_processor/models/compound_export.py binding_data_processor/models/compound/export/formats.py
    echo "Moved compound_export.py to export/formats.py"
fi

# Update imports
echo -e "\n${BLUE}Updating imports...${NC}"
find binding_data_processor -name "*.py" -exec sed -i '' 's/from models\./from models.compound./g' {} +
echo "Updated imports in Python files"

# Create test files
echo -e "\n${BLUE}Creating test files...${NC}"
for module in base ml enrichment analysis export; do
    mkdir -p "tests/models/compound/$module"
    touch "tests/models/compound/$module/test_${module}.py"
    echo "Created test_${module}.py"
done

# Run tests
echo -e "\n${BLUE}Running tests...${NC}"
if pytest tests/models/compound/; then
    echo -e "${GREEN}Tests passed successfully${NC}"
else
    echo -e "${RED}Tests failed - manual intervention needed${NC}"
fi

# Check coverage
echo -e "\n${BLUE}Checking test coverage...${NC}"
pytest --cov=binding_data_processor/models/compound/ tests/models/compound/

# Run linters
echo -e "\n${BLUE}Running linters...${NC}"
flake8 binding_data_processor/models/compound/
mypy binding_data_processor/models/compound/

echo -e "\n${GREEN}Setup complete!${NC}"
echo -e "Next steps:"
echo -e "1. Review moved files"
echo -e "2. Fix any broken imports"
echo -e "3. Run tests and fix failures"
echo -e "4. Update documentation"
