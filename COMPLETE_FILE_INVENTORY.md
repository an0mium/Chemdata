# Complete File Inventory

## Core Modules

### Models
```
binding_data_processor/models/
├── __init__.py
├── core.py
├── mixins.py
├── validation.py
├── enrichment.py
├── analysis.py
├── predictions.py
├── compound/
│   ├── __init__.py
│   ├── base/
│   │   ├── __init__.py
│   │   ├── core.py
│   │   ├── types.py
│   │   ├── validation.py
│   │   └── mixins.py
│   ├── ml/
│   │   ├── __init__.py
│   │   ├── predictors.py
│   │   ├── features.py
│   │   ├── training.py
│   │   └── ensemble.py
│   ├── analysis/
│   │   ├── __init__.py
│   │   ├── base.py
│   │   ├── binding_analysis.py
│   │   ├── activity_analysis.py
│   │   ├── safety_analysis.py
│   │   ├── property_analysis.py
│   │   └── sar_analysis.py
│   ├── enrichment/
│   │   ├── __init__.py
│   │   ├── validation/
│   │   │   ├── __init__.py
│   │   │   ├── schema.py
│   │   │   └── data.py
│   │   └── clients/
│   │       ├── __init__.py
│   │       ├── base.py
│   │       ├── http.py
│   │       ├── community.py
│   │       ├── social.py
│   │       └── swiss.py
│   └── export/
│       ├── __init__.py
│       ├── formats.py
│       └── validation.py
└── psychopharm/
    ├── __init__.py
    ├── base.py
    ├── binding.py
    ├── activity.py
    ├── safety.py
    ├── enrichment.py
    ├── community.py
    ├── compound.py
    ├── types.py
    └── analysis.py
```

### Pipeline
```
binding_data_processor/pipeline/
├── __init__.py
├── base.py
├── ml/
│   ├── __init__.py
│   ├── pipeline.py
│   └── core/
│       └── base.py
├── web/
│   ├── __init__.py
│   └── pipeline.py
├── analysis/
│   ├── __init__.py
│   ├── base.py
│   ├── binding.py
│   ├── activity.py
│   ├── safety.py
│   ├── properties.py
│   └── sar.py
├── export/
│   ├── __init__.py
│   └── pipeline.py
├── infrastructure/
│   ├── __init__.py
│   ├── resources.py
│   ├── errors.py
│   ├── cache.py
│   ├── checkpoints.py
│   ├── monitoring.py
│   ├── circuit_breaker.py
│   └── pipeline.py
└── processing/
    ├── __init__.py
    ├── pipeline.py
    ├── config.py
    └── social.py
```

### Data Sources
```
binding_data_processor/data_sources/
├── bindingdb.py
├── chembl.py
├── pubchem.py
└── pubmed.py
```

### Web Enrichment
```
binding_data_processor/web_enrichment/
├── __init__.py
├── base_client.py
├── http_client.py
├── http_client_enhanced.py
├── social_client.py
├── social_client_enhanced.py
├── community_client.py
├── community_client_enhanced.py
├── swiss_client.py
├── swiss_client_enhanced.py
├── manager.py
├── manager_enhanced.py
├── crawl4ai_client.py
├── llm_utils.py
├── validation/
│   ├── __init__.py
│   ├── data.py
│   ├── schema.py
│   └── enhanced.py
└── clients/
    ├── __init__.py
    ├── base.py
    ├── pubmed.py
    ├── scholar.py
    ├── sciencedirect.py
    ├── patents.py
    ├── reddit.py
    ├── bluelight.py
    └── swiss.py
```

### Web Interface
```
binding_data_processor/web/
├── __init__.py
├── app_enhanced.py
├── components/
│   ├── __init__.py
│   ├── compound_list.py
│   ├── compound_list_enhanced.py
│   ├── compound_details.py
│   ├── compound_detail_enhanced.py
│   ├── compound_search.py
│   ├── compound_search_enhanced.py
│   ├── compound_export.py
│   ├── compound_export_enhanced.py
│   ├── compound_visualization.py
│   ├── compound_visualization_enhanced.py
│   ├── compound_analysis.py
│   ├── compound_analysis_enhanced.py
│   └── compound_dashboard_enhanced.py
├── templates/
│   ├── base.html
│   ├── compound_dashboard.html
│   └── modals/
│       ├── search_modal.html
│       ├── filter_modal.html
│       └── export_modal.html
└── static/
    ├── css/
    │   └── style.css
    └── js/
        ├── app.js
        ├── package.json
        ├── .babelrc
        ├── .eslintrc.js
        ├── .prettierrc.js
        └── tests/
            ├── app.test.js
            ├── setup.js
            ├── globalSetup.js
            ├── globalTeardown.js
            └── fixtures/
                ├── compounds.json
                ├── predictions.json
                ├── web_data.json
                ├── literature_data.json
                └── analysis_results.json
```

### Scripts
```
scripts/
├── setup_project.sh
├── setup_dev.sh
├── setup_models.sh
├── setup_migration.sh
├── install_special_deps.sh
├── process_bindingdb.sh
├── enrich_compounds.sh
├── analyze_compounds.sh
├── generate_report.sh
├── setup_and_run.sh
├── run_pipeline.sh
├── run_web_app.py
├── predict_bbb.py
├── test_bbb_predictor.sh
└── manage_*.sh
```

### Documentation
```
docs/
└── source/
    ├── conf.py
    ├── index.rst
    ├── installation.rst
    ├── quickstart.rst
    ├── architecture.rst
    ├── deployment.rst
    ├── troubleshooting.rst
    ├── best_practices.rst
    ├── api_reference/
    │   ├── pipeline.rst
    │   ├── models.rst
    │   ├── predictors.rst
    │   └── web.rst
    ├── user_guide/
    │   ├── data_processing.rst
    │   ├── machine_learning.rst
    │   ├── web_enrichment.rst
    │   └── analysis.rst
    └── examples/
        ├── custom_pipeline.rst
        ├── ml_training.rst
        ├── web_interface.rst
        ├── data_enrichment.rst
        └── analysis_pipelines.rst
```

### Tests
```
tests/
├── conftest.py
├── test_enrichment.py
├── test_ml_predictions.py
├── test_pipeline.py
├── test_web_enrichment.py
├── test_web_interface.py
├── models/
│   └── compound/
│       ├── base/
│       │   └── test_core.py
│       ├── ml/
│       │   ├── test_predictors.py
│       │   ├── test_ensemble.py
│       │   ├── test_features.py
│       │   └── test_pipeline.py
│       ├── enrichment/
│       │   ├── test_web.py
│       │   ├── test_base.py
│       │   ├── test_validation.py
│       │   └── test_swiss.py
│       ├── analysis/
│       │   ├── test_base.py
│       │   ├── test_binding_analysis.py
│       │   ├── test_activity_analysis.py
│       │   ├── test_safety_analysis.py
│       │   ├── test_property_analysis.py
│       │   └── test_sar_analysis.py
│       └── export/
│           ├── test_formats.py
│           └── test_validation.py
└── web/
    ├── conftest.py
    ├── test_app.py
    ├── test_config.py
    ├── test_utils.py
    ├── test_fixtures.py
    ├── test_cache.py
    ├── test_database.py
    ├── test_templates.py
    ├── test_static.py
    ├── test_routes.py
    ├── test_api.py
    ├── test_server.py
    ├── test_monitoring.py
    ├── test_logging.py
    ├── test_auth.py
    ├── test_session.py
    ├── test_security.py
    ├── test_validation.py
    ├── test_request.py
    ├── test_response.py
    ├── test_errors.py
    ├── test_middleware.py
    └── components/
        ├── test_compound_dashboard.py
        ├── test_compound_search.py
        └── test_compound_details.py
```

### Configuration Files
```
./
├── .bandit.yaml
├── .coveragerc
├── .flake8
├── .gitignore
├── .pre-commit-config.yaml
├── docker-compose.yml
├── Dockerfile
├── pyproject.toml
├── pytest.ini
├── requirements.txt
└── setup.py
```

### Documentation Files
```
./
├── README.md
├── CHANGELOG.md
├── CONTRIBUTING.md
├── IMMEDIATE_STEPS.md
├── project_overview.md
├── project_analysis.md
├── project_plan.md
├── project_roadmap.md
├── codebase_status.md
├── codebase_analysis.md
├── file_analysis.md
├── findings_summary.md
├── implementation_files.md
├── action_plan.md
├── SUMMARY.md
├── CHECKLIST.md
├── COMPLETE_FILE_ANALYSIS.md
├── COMPLETE_FILE_INVENTORY.md
├── UPDATED_FILE_INVENTORY.md
├── MODELS_INVENTORY.md
├── PROCESSORS_INVENTORY.md
├── SCRIPTS_INVENTORY.md
├── TESTS_INVENTORY.md
├── INFRASTRUCTURE_INVENTORY.md
├── WEB_COMPONENTS_INVENTORY.md
├── MIGRATION_INVENTORY.md
├── DATA_SOURCE_INTEGRATION_STEPS.md
├── ML_PIPELINE_ENHANCEMENT_STEPS.md
├── WEB_INTERFACE_ENHANCEMENT_STEPS.md
├── WEB_ENRICHMENT_CONSOLIDATION_STEPS.md
├── MODEL_CONSOLIDATION_STEPS.md
├── MODEL_STRUCTURE_ANALYSIS.md
├── TESTING_STRATEGY.md
├── DEPLOYMENT_STRATEGY.md
├── MAINTENANCE_STRATEGY.md
├── DOCUMENTATION_STRATEGY.md
├── RELEASE_STRATEGY.md
├── INTEGRATION_STRATEGY.md
├── OPTIMIZATION_STRATEGY.md
├── MONITORING_STRATEGY.md
├── SCALING_STRATEGY.md
├── BACKUP_STRATEGY.md
└── SECURITY_STRATEGY.md
```

### Examples
```
examples/
├── README.md
├── pytest.ini
├── requirements-dev.txt
├── setup.py
├── setup_and_run.sh
├── data/
│   ├── example_compounds.tsv
│   └── example_compounds.json
├── scripts/
│   ├── enrich_compounds.py
│   ├── process_compounds.py
│   ├── run_bbb_prediction.sh
│   └── run_bbb_quickstart.sh
└── tests/
    ├── conftest.py
    └── test_enrich_compounds.py
