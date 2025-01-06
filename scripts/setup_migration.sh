#!/bin/bash

# Setup script for model consolidation migration

# Create feature branch (handle existing branch)
echo "Setting up feature branch..."
if git rev-parse --verify feature/model-consolidation >/dev/null 2>&1; then
    echo "Branch exists, checking out..."
    git checkout feature/model-consolidation
else
    echo "Creating new branch..."
    git checkout -b feature/model-consolidation
fi

# Create directories
echo "Creating directories..."
mkdir -p backups/models
mkdir -p backups/tests
mkdir -p reports/coverage
mkdir -p reports/test-results

# Backup critical files
echo "Backing up critical files..."
cp binding_data_processor/models/compound_base.py backups/models/
cp binding_data_processor/models/compound.py backups/models/
cp binding_data_processor/models/psychopharm/base.py backups/models/
cp binding_data_processor/models/psychopharm/compound.py backups/models/
cp binding_data_processor/models/psychopharm/tests/test_base.py backups/tests/
cp binding_data_processor/models/psychopharm/tests/test_compound.py backups/tests/

# Run baseline tests
echo "Running baseline tests..."
pytest binding_data_processor/models/psychopharm/tests/test_base.py binding_data_processor/models/psychopharm/tests/test_compound.py -v > backups/baseline_test_results.txt

# Run tests with coverage
echo "Running tests with coverage..."
pytest --cov=binding_data_processor binding_data_processor/models/psychopharm/tests/test_base.py binding_data_processor/models/psychopharm/tests/test_compound.py > backups/coverage_test_results.txt

# Document current coverage
echo "Documenting test coverage..."
coverage report > backups/baseline_coverage.txt

# Create migration log
echo "Creating migration log..."
cat > backups/migration.log << EOL
Model Consolidation Migration Log
===============================

Started: $(date)

Baseline Test Results: See baseline_test_results.txt
Coverage Test Results: See coverage_test_results.txt
Baseline Coverage: See baseline_coverage.txt

Critical Files Backed Up:
- compound_base.py
- compound.py
- psychopharm/base.py
- psychopharm/compound.py
- psychopharm/tests/test_base.py
- psychopharm/tests/test_compound.py

Next Steps:
1. Review compound_base.py and psychopharm/base.py
2. Map overlapping functionality
3. Begin migration of unique features
4. Update tests
5. Verify coverage

Branch Status:
$(git status)

Test Summary:
$(tail -n 10 backups/baseline_test_results.txt)

Coverage Summary:
$(tail -n 10 backups/baseline_coverage.txt)
EOL

echo "Setup complete. See backups/migration.log for details."
