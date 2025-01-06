# Model Structure Analysis

## Directory Overview

### Legacy Model Files (binding_data_processor/models/)
1. compound.py
   - Purpose: Main compound model
   - Status: To be migrated
   - Action: Merge into psychopharm/compound.py
   - Notes: Contains valuable features that need preservation
   - Integration: Check test coverage in psychopharm/tests/test_compound.py

2. compound_base.py
   - Purpose: Base compound functionality
   - Status: To be migrated
   - Action: Merge into psychopharm/base.py
   - Notes: Core functionality to be preserved
   - Integration: Check test coverage in psychopharm/tests/test_base.py

3. compound_ml.py
   - Purpose: ML prediction functionality
   - Status: To be migrated
   - Action: Merge into psychopharm/binding.py
   - Notes: ML features to be integrated
   - Integration: Check test coverage in psychopharm/tests/test_ml_models.py

4. compound_enrichment.py
   - Purpose: Web data enrichment
   - Status: To be migrated
   - Action: Merge into psychopharm/enrichment.py
   - Notes: Web features to be preserved
   - Integration: Check test coverage in psychopharm/tests/test_enrichment.py

5. compound_analysis.py
   - Purpose: Analysis functionality
   - Status: To be migrated
   - Action: Split between activity.py and safety.py
   - Notes: Analysis features to be preserved
   - Integration: Check test coverage in respective test files

### New Model Structure (binding_data_processor/models/psychopharm/)
1. base.py
   - Purpose: Core psychopharm model
   - Status: Well integrated
   - Keep: Yes
   - Notes: Primary base class
   - Tests: test_base.py shows good coverage

2. binding.py
   - Purpose: Binding predictions
   - Status: Well integrated
   - Keep: Yes
   - Notes: Needs ML integration from compound_ml.py
   - Tests: test_binding.py shows good coverage

3. activity.py
   - Purpose: Activity analysis
   - Status: Well integrated
   - Keep: Yes
   - Notes: Will receive analysis features
   - Tests: test_activity.py shows good coverage

4. safety.py
   - Purpose: Safety assessment
   - Status: Well integrated
   - Keep: Yes
   - Notes: Will receive analysis features
   - Tests: test_safety.py shows good coverage

5. enrichment.py
   - Purpose: Web data enrichment
   - Status: Well integrated
   - Keep: Yes
   - Notes: Will receive web features
   - Tests: test_enrichment.py shows good coverage

### Enhanced Components (binding_data_processor/models/compound/)
1. enhanced.py
   - Purpose: Enhanced compound features
   - Status: New features
   - Action: Merge with base components
   - Notes: Contains improvements
   - Integration: Preserve enhancements during merge

### Test Coverage
1. psychopharm/tests/
   - Good coverage of core functionality
   - Comprehensive test suite
   - Integration tests present
   - Performance tests needed

2. compound/tests/
   - Legacy test coverage
   - Some unique test cases
   - Migration needed
   - Coverage gaps identified

## Integration Status

### Well Integrated Components
1. Core Models
   - psychopharm/base.py
   - psychopharm/binding.py
   - psychopharm/activity.py
   - psychopharm/safety.py
   - psychopharm/enrichment.py

2. Test Suite
   - test_base.py
   - test_binding.py
   - test_activity.py
   - test_safety.py
   - test_enrichment.py

### Needs Migration
1. Legacy Models
   - compound.py -> psychopharm/compound.py
   - compound_base.py -> psychopharm/base.py
   - compound_ml.py -> psychopharm/binding.py
   - compound_enrichment.py -> psychopharm/enrichment.py
   - compound_analysis.py -> activity.py/safety.py

2. Legacy Tests
   - Migrate unique test cases
   - Preserve edge case coverage
   - Update test dependencies
   - Add missing coverage

## Migration Steps

1. Code Analysis
   - Review all legacy files
   - Identify unique features
   - Map dependencies
   - Plan test migration

2. Feature Preservation
   - Document core functionality
   - Identify enhancements
   - Map test coverage
   - Plan integration tests

3. Integration Process
   - Migrate core features first
   - Preserve enhancements
   - Update dependencies
   - Run test suite
   - Fix integration issues
   - Add missing tests

4. Validation
   - Run full test suite
   - Check coverage metrics
   - Verify functionality
   - Document changes

## Next Steps

1. Immediate Actions
   - Start with compound.py migration
   - Focus on core functionality
   - Preserve test coverage
   - Document progress

2. Short-term Goals
   - Complete base class migration
   - Integrate ML features
   - Update documentation
   - Expand test suite

3. Long-term Goals
   - Full integration
   - Enhanced features
   - Complete coverage
   - Updated docs
