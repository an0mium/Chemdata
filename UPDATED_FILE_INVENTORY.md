# Updated File Inventory

## Overview

This document provides a comprehensive inventory of all files in the codebase, along with their status, integration level, and relationships to other files. The inventory is organized by major components and includes both existing files and planned files.

## Core Components

### Root Directory Structure
```
/
├── binding_data_processor/    # Main package directory
├── docs/                      # Documentation
├── examples/                  # Example code
├── scripts/                   # Utility scripts
├── tests/                     # Test files
├── web/                       # Web interface
└── web_enrichment/           # Web data enrichment
```

### Key Configuration Files
1. Package Configuration
- pyproject.toml (exists, active)
- setup.cfg (exists, active)
- setup.py (exists, active)
- requirements.txt (exists, active)

2. Development Configuration
- .pre-commit-config.yaml (exists, active)
- .bandit.yaml (exists, active)
- .flake8 (exists, active)
- .gitignore (exists, active)
- pytest.ini (exists, active)
- .coveragerc (exists, active)

3. Container Configuration
- Dockerfile (exists, active)
- docker-compose.yml (exists, active)

## Binding Data Processor Package

### Core Package Files
1. Base Files
- binding_data_processor/__init__.py (exists, active)
- binding_data_processor/config.py (exists, active)
- binding_data_processor/main.py (exists, active)
- binding_data_processor/pipeline.py (exists, active)
- binding_data_processor/cli.py (exists, active)

### Models Module
1. Base Models
- binding_data_processor/models/__init__.py (exists, active)
- binding_data_processor/models/core.py (exists, active)
- binding_data_processor/models/validation.py (exists, active)
- binding_data_processor/models/mixins.py (exists, active)

2. Compound Models
- binding_data_processor/models/compound/__init__.py (exists, active)
- binding_data_processor/models/compound/base.py (exists, active)
- binding_data_processor/models/compound/types.py (exists, active)
- binding_data_processor/models/compound/ml.py (exists, active)
- binding_data_processor/models/compound/enrichment.py (exists, active)
- binding_data_processor/models/compound/analysis.py (exists, active)
- binding_data_processor/models/compound/enhanced.py (exists, active)

3. Psychopharm Models
- binding_data_processor/models/psychopharm/__init__.py (exists, active)
- binding_data_processor/models/psychopharm/base.py (exists, active)
- binding_data_processor/models/psychopharm/types.py (exists, active)
- binding_data_processor/models/psychopharm/binding.py (exists, active)
- binding_data_processor/models/psychopharm/activity.py (exists, active)
- binding_data_processor/models/psychopharm/safety.py (exists, active)
- binding_data_processor/models/psychopharm/enrichment.py (exists, active)
- binding_data_processor/models/psychopharm/analysis.py (exists, active)
- binding_data_processor/models/psychopharm/community.py (exists, active)

### Pipeline Module
1. Core Pipeline
- binding_data_processor/pipeline/__init__.py (exists, active)
- binding_data_processor/pipeline/base.py (exists, active)
- binding_data_processor/pipeline/ml.py (exists, active)
- binding_data_processor/pipeline/web.py (exists, active)
- binding_data_processor/pipeline/validation.py (exists, active)

2. Analysis Pipeline
- binding_data_processor/pipeline/analysis/__init__.py (exists, active)
- binding_data_processor/pipeline/analysis/base.py (exists, active)
- binding_data_processor/pipeline/analysis/binding.py (exists, active)
- binding_data_processor/pipeline/analysis/activity.py (exists, active)
- binding_data_processor/pipeline/analysis/properties.py (exists, active)
- binding_data_processor/pipeline/analysis/safety.py (exists, active)
- binding_data_processor/pipeline/analysis/sar.py (exists, active)

3. Infrastructure
- binding_data_processor/pipeline/infrastructure/__init__.py (exists, active)
- binding_data_processor/pipeline/infrastructure/cache.py (exists, active)
- binding_data_processor/pipeline/infrastructure/checkpoints.py (exists, active)
- binding_data_processor/pipeline/infrastructure/circuit_breaker.py (exists, active)
- binding_data_processor/pipeline/infrastructure/errors.py (exists, active)
- binding_data_processor/pipeline/infrastructure/monitoring.py (exists, active)
- binding_data_processor/pipeline/infrastructure/resources.py (exists, active)

### Processors Module
1. Structure Processing
- binding_data_processor/processors/structure/properties/base.py (exists, active)
- binding_data_processor/processors/structure/properties/descriptors.py (exists, active)
- binding_data_processor/processors/structure/properties/similarity.py (exists, active)
- binding_data_processor/processors/structure/properties/standardization.py (exists, active)

2. Psychopharm Predictors
- binding_data_processor/processors/psychopharm/predictors/base.py (exists, active)
- binding_data_processor/processors/psychopharm/predictors/bbb/ (directory exists, active)
- binding_data_processor/processors/psychopharm/predictors/nootropic/ (directory exists, active)
- binding_data_processor/processors/psychopharm/predictors/toxicity/ (directory exists, active)
- binding_data_processor/processors/psychopharm/predictors/abuse.py (exists, active)
- binding_data_processor/processors/psychopharm/predictors/receptors.py (exists, active)
- binding_data_processor/processors/psychopharm/predictors/psychoactive.py (exists, active)

### Web Module
1. Core Web
- binding_data_processor/web/app.py (exists, active)
- binding_data_processor/web/app_enhanced.py (exists, active)
- binding_data_processor/web/report_generator.py (exists, active)

2. Components
- binding_data_processor/web/components/compound_list.py (exists, active)
- binding_data_processor/web/components/compound_details.py (exists, active)
- binding_data_processor/web/components/compound_search.py (exists, active)
- binding_data_processor/web/components/compound_dashboard_enhanced.py (exists, active)
- binding_data_processor/web/components/compound_list_enhanced.py (exists, active)
- binding_data_processor/web/components/compound_detail_enhanced.py (exists, active)
- binding_data_processor/web/components/compound_search_enhanced.py (exists, active)
- binding_data_processor/web/components/compound_visualization_enhanced.py (exists, active)

3. Templates
- binding_data_processor/web/templates/base.html (exists, active)
- binding_data_processor/web/templates/compound_dashboard.html (exists, active)
- binding_data_processor/web/templates/modals/search_modal.html (exists, active)
- binding_data_processor/web/templates/modals/filter_modal.html (exists, active)
- binding_data_processor/web/templates/modals/export_modal.html (exists, active)

### Web Enrichment Module
1. Core Clients
- binding_data_processor/web_enrichment/base_client.py (exists, active)
- binding_data_processor/web_enrichment/http_client.py (exists, active)
- binding_data_processor/web_enrichment/manager.py (exists, active)
- binding_data_processor/web_enrichment/llm_utils.py (exists, active)

2. Data Source Clients
- binding_data_processor/web_enrichment/clients/reddit.py (exists, active)
- binding_data_processor/web_enrichment/clients/bluelight.py (exists, active)
- binding_data_processor/web_enrichment/clients/pubmed.py (exists, active)
- binding_data_processor/web_enrichment/clients/patents.py (exists, active)
- binding_data_processor/web_enrichment/clients/scholar.py (exists, active)
- binding_data_processor/web_enrichment/clients/sciencedirect.py (exists, active)

3. Enhanced Clients
- binding_data_processor/web_enrichment/http_client_enhanced.py (exists, active)
- binding_data_processor/web_enrichment/manager_enhanced.py (exists, active)
- binding_data_processor/web_enrichment/community_client_enhanced.py (exists, active)
- binding_data_processor/web_enrichment/social_client_enhanced.py (exists, active)

## Integration Status

### Components to Merge
1. Enhanced Web Components
- Merge enhanced versions into base components
- Preserve new features
- Update tests
- Update documentation

2. Enhanced Web Enrichment Clients
- Merge enhanced versions into base clients
- Preserve new features
- Update tests
- Update documentation

3. Legacy Root Files
- Move functionality to appropriate modules
- Update imports
- Add tests
- Remove legacy files

### Next Steps

1. Code Migration
- Start with core utilities
- Then move API clients
- Then move processors
- Finally move entry points

2. Documentation Update
- Merge analysis documents
- Update API docs
- Add migration guides
- Update examples

3. Test Coverage
- Add tests for migrated code
- Update existing tests
- Verify coverage
- Add missing cases

4. Configuration Cleanup
- Consolidate settings
- Update build process
- Update test process
- Update linting rules

This inventory will be updated as files are migrated and consolidated.
