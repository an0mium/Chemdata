-- Main schema file that handles proper ordering and initialization
-- This file ensures tables are created in the correct order to handle dependencies

-- Create base reference tables that other tables depend on
CREATE TABLE IF NOT EXISTS toxicity_endpoints (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    description TEXT,
    unit TEXT,
    endpoint_type TEXT
);
COMMENT ON TABLE toxicity_endpoints IS 'Reference data for toxicity measurement endpoints';

CREATE TABLE IF NOT EXISTS receptor_families (
    id SERIAL PRIMARY KEY,
    family_name TEXT NOT NULL,
    description TEXT,
    protein_type TEXT,
    signaling_type TEXT,
    receptor_class TEXT,
    primary_effects TEXT[],
    therapeutic_areas TEXT[]
);
COMMENT ON TABLE receptor_families IS 'Classification of receptor families and their properties';

CREATE TABLE IF NOT EXISTS therapeutic_classes (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    description TEXT
);
COMMENT ON TABLE therapeutic_classes IS 'Classification of therapeutic uses for compounds';

CREATE TABLE IF NOT EXISTS subjective_effect_categories (
    id SERIAL PRIMARY KEY,
    category_name TEXT NOT NULL,
    description TEXT,
    variability_factors TEXT[]
);
COMMENT ON TABLE subjective_effect_categories IS 'Categories of subjective effects and their variability';

CREATE TABLE IF NOT EXISTS web_data_sources (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    source_type TEXT,
    reliability_score FLOAT,
    update_frequency INTERVAL
);
COMMENT ON TABLE web_data_sources IS 'External data sources and their reliability metrics';

CREATE TABLE IF NOT EXISTS toxicity_mechanisms (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    description TEXT,
    mechanism_type TEXT
);
COMMENT ON TABLE toxicity_mechanisms IS 'Mechanisms by which compounds can cause toxicity';

CREATE TABLE IF NOT EXISTS organ_toxicity_types (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    description TEXT,
    assessment_methods TEXT[]
);
COMMENT ON TABLE organ_toxicity_types IS 'Organ-specific toxicity types and assessment methods';

CREATE TABLE IF NOT EXISTS safety_thresholds (
    id SERIAL PRIMARY KEY,
    parameter_name TEXT NOT NULL,
    threshold_value FLOAT,
    unit TEXT,
    severity_level TEXT,
    description TEXT
);
COMMENT ON TABLE safety_thresholds IS 'Safety threshold values for various parameters';

CREATE TABLE IF NOT EXISTS monitoring_parameters (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    category TEXT,
    frequency TEXT,
    monitoring_method TEXT,
    alert_conditions TEXT[]
);
COMMENT ON TABLE monitoring_parameters IS 'Parameters to monitor for safety assessment';

CREATE TABLE IF NOT EXISTS intervention_thresholds (
    id SERIAL PRIMARY KEY,
    parameter_name TEXT NOT NULL,
    threshold_level TEXT,
    intervention_type TEXT,
    urgency_level TEXT,
    intervention_protocol TEXT
);
COMMENT ON TABLE intervention_thresholds IS 'Thresholds for clinical interventions';

-- Load consolidated core schema first
\i consolidated/00_core.sql

-- Load base compound and receptor tables
\i compounds/01_base.sql
\i compounds/02_binding.sql
\i compounds/03_molecular_descriptors.sql
\i compounds/05_quantum_criticality.sql
\i receptors/01_base.sql
\i receptors/02_subtypes.sql

-- Load reference data
\i reference_data/01_toxicity_endpoints.sql
\i reference_data/02_receptor_families.sql
\i reference_data/03_therapeutic_classes.sql
\i reference_data/04_subjective_effects.sql
\i reference_data/05_web_sources.sql
\i reference_data/06_toxicity_mechanisms.sql
\i reference_data/07_organ_toxicity.sql
\i reference_data/08_safety_thresholds.sql
\i reference_data/09_monitoring_parameters.sql
\i reference_data/10_intervention_thresholds.sql
\i reference_data/11_quantum_parameters.sql

-- Load core data tables
\i safety/01_toxicity_data.sql
\i safety/02_risk_assessment.sql
\i pharmacology/01_pharmacology_data.sql
\i pharmacology/02_profiles.sql
\i clinical/01_clinical_data.sql
\i clinical/02_experience.sql
\i regulatory/01_regulatory_data.sql

-- Load analysis and machine learning tables
\i analysis/01_literature.sql
\i analysis/02_machine_learning.sql
\i analysis/03_sar_data.sql
\i research/01_quantum_research.sql

-- Load literature and documentation tables
\i literature/01_literature_data.sql
\i literature/02_analysis.sql

-- Load web interface tables
\i web/web_templates.sql
\i web/web_settings.sql
\i web/03_data_integration.sql

-- Load social media and community data tables
\i social/01_community_data.sql
\i social/02_media_data.sql
\i social/03_community_platforms.sql
\i social/04_additional_forums.sql
\i social/05_more_forums.sql
\i social/06_erowid.sql
\i social/07_psychonautwiki.sql
\i social/08_tripsit.sql
\i social/09_longecity.sql
\i social/10_reddit.sql
\i social/11_twitter.sql
\i social/12_cross_platform_analytics.sql
\i social/13_alerts.sql

-- Load monitoring and ML model tables last
\i ml/models.sql
\i monitoring/alerts.sql

-- Create GiST index for quantum data
CREATE INDEX IF NOT EXISTS idx_density_grid_points ON electronic_structure USING gist (density_grid_points);
CREATE INDEX IF NOT EXISTS idx_quantum_calculations_compound ON quantum_calculations(compound_id);
CREATE INDEX IF NOT EXISTS idx_quantum_findings_type ON quantum_research_findings(finding_type);
CREATE INDEX IF NOT EXISTS idx_quantum_correlations_property ON quantum_structure_correlations(property_type);

-- Verify all audit triggers are properly set up
DO $$
DECLARE
    r RECORD;
    trigger_exists BOOLEAN;
BEGIN
    FOR r IN (SELECT tablename FROM pg_tables WHERE schemaname = 'public')
    LOOP
        -- Check if trigger already exists
        SELECT EXISTS (
            SELECT 1 
            FROM pg_trigger 
            WHERE tgname = 'audit_' || r.tablename || '_trigger'
        ) INTO trigger_exists;
        
        -- Create trigger if it doesn't exist
        IF NOT trigger_exists THEN
            EXECUTE format(
                'CREATE TRIGGER audit_%I_trigger 
                 AFTER INSERT OR UPDATE OR DELETE ON %I 
                 FOR EACH ROW EXECUTE FUNCTION audit_trigger_func()',
                r.tablename, r.tablename
            );
        END IF;
    END LOOP;
END
$$;

-- Verify schema integrity
DO $$
DECLARE
    missing_tables TEXT[];
    table_record RECORD;
BEGIN
    -- Check core tables exist
    FOR table_record IN (
        SELECT tablename 
        FROM pg_tables 
        WHERE schemaname = 'public'
        AND tablename IN (
            -- Core tables
            'compounds', 'receptor_families', 'binding_data',
            'toxicity_endpoints', 'therapeutic_classes',
            'subjective_effect_categories', 'web_data_sources',
            'toxicity_mechanisms', 'organ_toxicity_types',
            -- Data quality tables
            'data_quality_metrics', 'validation_rules',
            'validation_results', 'confidence_scores',
            'data_quality_issues', 'data_quality_thresholds',
            -- Social media tables
            'community_data', 'reddit_data', 'twitter_data',
            'bluelight_data',
            -- Quantum tables
            'electronic_structure', 'quantum_critical_params',
            'quantum_dynamics', 'phase_transitions',
            'scaling_analysis', 'quantum_observables',
            'quantum_basis_sets', 'quantum_functionals',
            'quantum_observables_ref', 'phase_transition_types',
            'quantum_research_projects', 'quantum_calculations',
            'quantum_research_findings', 'quantum_structure_correlations'
        )
    )
    LOOP
        IF NOT EXISTS (SELECT 1 FROM pg_tables WHERE tablename = table_record.tablename) THEN
            missing_tables := array_append(missing_tables, table_record.tablename);
        END IF;
    END LOOP;

    -- Raise error if any core tables are missing
    IF array_length(missing_tables, 1) > 0 THEN
        RAISE EXCEPTION 'Schema integrity check failed. Missing tables: %', array_to_string(missing_tables, ', ');
    END IF;
END
$$;
