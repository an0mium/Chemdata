# Migration Inventory

## Overview

This document outlines the specific migration steps needed to consolidate and improve the codebase. It identifies files that need to be merged, moved, or updated.

## Migration Tasks

### 1. Legacy Root Files Migration

These files in the root directory need to be moved to appropriate modules:

1. API Clients
- api_client.py -> web_enrichment/http_client.py
- chembl_client.py -> data_sources/chembl.py
- pubmed_processor.py -> data_sources/pubmed.py

2. Infrastructure
- cache_manager.py -> pipeline/infrastructure/cache.py
- checkpoint_manager.py -> pipeline/infrastructure/checkpoints.py
- logger.py -> pipeline/infrastructure/logging.py

3. Core Processing
- binding_data_processor.py -> binding_data_processor/core.py
- chemical_properties.py -> processors/structure/properties/core.py
- structure_utils.py -> processors/structure/utils.py

4. Entry Points
- cli.py -> binding_data_processor/cli.py
- run.py -> scripts/run.py

### 2. Enhanced Components Consolidation

These enhanced components need to be merged with their base versions:

1. Web Components
- compound_dashboard_enhanced.py -> compound_dashboard.py
- compound_list_enhanced.py -> compound_list.py
- compound_detail_enhanced.py -> compound_detail.py
- compound_search_enhanced.py -> compound_search.py
- compound_visualization_enhanced.py -> compound_visualization.py

2. Web Enrichment Clients
- http_client_enhanced.py -> http_client.py
- manager_enhanced.py -> manager.py
- community_client_enhanced.py -> community_client.py
- social_client_enhanced.py -> social_client.py
- swiss_client_enhanced.py -> swiss_client.py

3. BBB Predictor Components
- bbb_enhanced.py -> bbb/base.py
- bbb_web_enrichment.py -> bbb/enrichment.py

4. Nootropic Predictor Components
- nootropic_enhanced.py -> nootropic/base.py

5. Toxicity Predictor Components
- toxicity_enhanced.py -> toxicity/base.py

### 3. Documentation Consolidation

1. Analysis Documents
- WEB_ENRICHMENT_ANALYSIS.md -> docs/analysis/web_enrichment.rst
- WEB_COMPONENTS_ANALYSIS.md -> docs/analysis/web_components.rst
- MODEL_STRUCTURE_ANALYSIS.md -> docs/analysis/model_structure.rst
- FINDINGS_SUMMARY.md -> docs/analysis/findings.rst

2. Strategy Documents
- TESTING_STRATEGY.md -> docs/strategy/testing.rst
- DEPLOYMENT_STRATEGY.md -> docs/strategy/deployment.rst
- MAINTENANCE_STRATEGY.md -> docs/strategy/maintenance.rst
- DOCUMENTATION_STRATEGY.md -> docs/strategy/documentation.rst
- RELEASE_STRATEGY.md -> docs/strategy/release.rst
- INTEGRATION_STRATEGY.md -> docs/strategy/integration.rst
- OPTIMIZATION_STRATEGY.md -> docs/strategy/optimization.rst
- MONITORING_STRATEGY.md -> docs/strategy/monitoring.rst
- SCALING_STRATEGY.md -> docs/strategy/scaling.rst
- BACKUP_STRATEGY.md -> docs/strategy/backup.rst
- SECURITY_STRATEGY.md -> docs/strategy/security.rst

3. Implementation Documents
- MODEL_CONSOLIDATION_STEPS.md -> docs/implementation/model_consolidation.rst
- DATA_SOURCE_INTEGRATION_STEPS.md -> docs/implementation/data_source_integration.rst
- ML_PIPELINE_ENHANCEMENT_STEPS.md -> docs/implementation/ml_pipeline_enhancement.rst
- WEB_INTERFACE_ENHANCEMENT_STEPS.md -> docs/implementation/web_interface_enhancement.rst
- WEB_ENRICHMENT_CONSOLIDATION_STEPS.md -> docs/implementation/web_enrichment_consolidation.rst

### 4. Test Coverage Enhancement

1. Model Tests
- Add tests for migrated model functionality
- Update existing model tests
- Add integration tests between models

2. Web Component Tests
- Add tests for enhanced web components
- Update existing web component tests
- Add integration tests between components

3. Web Enrichment Tests
- Add tests for enhanced clients
- Update existing client tests
- Add integration tests between clients

4. Infrastructure Tests
- Add tests for cache functionality
- Add tests for checkpoint functionality
- Add tests for monitoring functionality

## Migration Steps

### Phase 1: Infrastructure (Week 1)
1. Move infrastructure files
2. Update imports
3. Add tests
4. Update documentation

### Phase 2: API Clients (Week 2)
1. Move API client files
2. Update imports
3. Add tests
4. Update documentation

### Phase 3: Core Processing (Week 3)
1. Move processing files
2. Update imports
3. Add tests
4. Update documentation

### Phase 4: Entry Points (Week 4)
1. Move entry point files
2. Update imports
3. Add tests
4. Update documentation

### Phase 5: Enhanced Components (Weeks 5-6)
1. Merge enhanced web components
2. Merge enhanced clients
3. Merge enhanced predictors
4. Update tests
5. Update documentation

### Phase 6: Documentation (Weeks 7-8)
1. Convert markdown to RST
2. Organize documentation structure
3. Update cross-references
4. Add new sections

## Success Criteria

1. Code Quality
- All tests passing
- No duplicate code
- Clear documentation
- Type hints complete

2. Functionality
- All features preserved
- Enhanced capabilities
- Good performance
- Clean interfaces

3. Documentation
- Updated API docs
- Clear examples
- Migration guide
- Best practices

4. Testing
- High coverage
- Integration tests
- Performance tests
- Error cases

## Monitoring and Validation

1. Code Quality Metrics
- Test coverage reports
- Linting reports
- Type checking reports
- Complexity metrics

2. Performance Metrics
- Response times
- Memory usage
- CPU usage
- I/O operations

3. Documentation Metrics
- Documentation coverage
- Example coverage
- API documentation completeness
- User guide completeness

4. Integration Metrics
- Build success rate
- Test pass rate
- Integration test coverage
- End-to-end test coverage

This inventory will be updated as migration progresses to track completion status and any issues encountered.
