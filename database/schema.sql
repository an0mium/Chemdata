-- Main schema file for Chemdata database
-- This file loads all schema components in the correct dependency order

-- Set configuration
SET client_min_messages TO warning;
SET timezone TO 'UTC';

-- Enable required extensions first
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";
CREATE EXTENSION IF NOT EXISTS "unaccent";
CREATE EXTENSION IF NOT EXISTS "pg_trgm";

-- Create audit log table first since it's referenced by the trigger function
CREATE TABLE audit_log (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    table_name text NOT NULL,
    record_id uuid NOT NULL,
    action text NOT NULL,
    old_data jsonb,
    new_data jsonb,
    changed_by text,
    changed_at timestamptz NOT NULL DEFAULT now()
);

-- Create common functions
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = now();
    RETURN NEW;
END;
$$ language 'plpgsql';

CREATE OR REPLACE FUNCTION audit_trigger_func()
RETURNS TRIGGER AS $$
BEGIN
    IF TG_OP = 'INSERT' THEN
        INSERT INTO audit_log (table_name, record_id, action, new_data)
        VALUES (TG_TABLE_NAME, NEW.id, 'INSERT', row_to_json(NEW));
    ELSIF TG_OP = 'UPDATE' THEN
        INSERT INTO audit_log (table_name, record_id, action, old_data, new_data)
        VALUES (TG_TABLE_NAME, NEW.id, 'UPDATE', row_to_json(OLD), row_to_json(NEW));
    ELSIF TG_OP = 'DELETE' THEN
        INSERT INTO audit_log (table_name, record_id, action, old_data)
        VALUES (TG_TABLE_NAME, OLD.id, 'DELETE', row_to_json(OLD));
    END IF;
    RETURN NULL;
END;
$$ LANGUAGE plpgsql;

-- Load core functionality
\i database/schema/core/00_extensions.sql
\i database/schema/core/01_functions.sql

-- Load reference data tables first
-- These need to be created before tables that reference them
\i database/schema/reference_data/01_toxicity_endpoints.sql
\i database/schema/reference_data/02_receptor_families.sql
\i database/schema/reference_data/03_therapeutic_classes.sql
\i database/schema/reference_data/04_subjective_effects.sql
\i database/schema/reference_data/05_web_sources.sql
\i database/schema/reference_data/06_toxicity_mechanisms.sql
\i database/schema/reference_data/07_organ_toxicity.sql
\i database/schema/reference_data/08_safety_thresholds.sql
\i database/schema/reference_data/09_monitoring_parameters.sql
\i database/schema/reference_data/10_intervention_thresholds.sql

-- Load base compound tables
-- These need to exist before binding data tables
\i database/schema/compounds/01_base.sql

-- Load base receptor tables
-- These need to exist before receptor subtypes and binding data
\i database/schema/receptors/01_base.sql

-- Load binding data tables
-- These reference both compounds and receptors
\i database/schema/compounds/02_binding.sql

-- Load receptor subtype tables
-- These reference receptor families
\i database/schema/receptors/02_subtypes.sql

-- Load safety and pharmacology data tables
-- These reference compounds and toxicity endpoints
\i database/schema/safety/01_toxicity_data.sql
\i database/schema/pharmacology/01_pharmacology_data.sql

-- Load clinical and regulatory data tables
-- These reference compounds and various reference data
\i database/schema/clinical/01_clinical_data.sql
\i database/schema/regulatory/01_regulatory_data.sql

-- Load social and literature data tables
-- These reference compounds and web sources
\i database/schema/social/01_community_data.sql
\i database/schema/literature/01_literature_data.sql

-- Load analysis data tables
-- These reference compounds and various other tables
\i database/schema/analysis/01_literature.sql
\i database/schema/analysis/02_machine_learning.sql
\i database/schema/analysis/03_sar_data.sql

-- Verify schema installation
DO $$
BEGIN
    -- Check that key tables exist
    ASSERT EXISTS (
        SELECT FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename = 'compounds'
    ), 'compounds table not found';

    ASSERT EXISTS (
        SELECT FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename = 'receptor_families'
    ), 'receptor_families table not found';

    ASSERT EXISTS (
        SELECT FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename = 'audit_log'
    ), 'audit_log table not found';

    -- Check that extensions are installed
    ASSERT EXISTS (
        SELECT FROM pg_extension 
        WHERE extname = 'uuid-ossp'
    ), 'uuid-ossp extension not found';

    ASSERT EXISTS (
        SELECT FROM pg_extension 
        WHERE extname = 'unaccent'
    ), 'unaccent extension not found';

    ASSERT EXISTS (
        SELECT FROM pg_extension 
        WHERE extname = 'pg_trgm'
    ), 'pg_trgm extension not found';

    -- Check that functions exist
    ASSERT EXISTS (
        SELECT FROM pg_proc 
        WHERE proname = 'update_updated_at_column'
    ), 'update_updated_at_column function not found';

    ASSERT EXISTS (
        SELECT FROM pg_proc 
        WHERE proname = 'audit_trigger_func'
    ), 'audit_trigger_func function not found';

    -- Check that triggers are created
    ASSERT EXISTS (
        SELECT FROM pg_trigger
        WHERE tgname = 'audit_compounds_trigger'
    ), 'compounds audit trigger not found';

    -- Check that indexes are created
    ASSERT EXISTS (
        SELECT FROM pg_indexes
        WHERE tablename = 'compounds' 
        AND indexname = 'idx_compounds_name'
    ), 'compounds name index not found';

    -- Log success
    RAISE NOTICE 'Schema verification complete - all required objects exist';
END $$;

-- Set search path
SET search_path TO public;

-- Final notice
DO $$
BEGIN
    RAISE NOTICE 'Schema installation complete';
    RAISE NOTICE 'Database is ready for use';
    RAISE NOTICE 'Installed extensions: uuid-ossp, unaccent, pg_trgm';
    RAISE NOTICE 'Created audit logging system';
    RAISE NOTICE 'Verified all required tables, functions, triggers and indexes';
END $$;
