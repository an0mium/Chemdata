-- Main schema file for Chemdata database
-- This file loads all schema components in the correct dependency order

-- Set configuration
SET client_min_messages TO warning;
SET timezone TO 'UTC';

-- Enable required extensions first
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";
CREATE EXTENSION IF NOT EXISTS "unaccent";
CREATE EXTENSION IF NOT EXISTS "pg_trgm";

-- Create audit log table first since it's referenced by triggers
CREATE TABLE IF NOT EXISTS audit_log (
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

-- SMILES validation function
CREATE OR REPLACE FUNCTION validate_smiles(smiles text)
RETURNS boolean AS $$
BEGIN
    -- Basic SMILES validation
    -- Check for balanced parentheses and brackets
    IF (
        LENGTH(REGEXP_REPLACE(smiles, '[^\(\)]', '', 'g')) % 2 != 0 OR
        LENGTH(REGEXP_REPLACE(smiles, '[^\[\]]', '', 'g')) % 2 != 0
    ) THEN
        RETURN false;
    END IF;

    -- Check for valid atoms and special characters
    IF NOT REGEXP_MATCH(smiles, '^[A-Za-z0-9\(\)\[\]\+\-\=\#\$\:\\/\.\@\*\{\}]+$') THEN
        RETURN false;
    END IF;

    -- Check for valid atom symbols
    IF NOT REGEXP_MATCH(smiles, '^[A-Z][a-z]?|[a-z]|[\d\(\)\[\]\+\-\=\#\$\:\\/\.\@\*\{\}]') THEN
        RETURN false;
    END IF;

    -- Check for valid bond symbols
    IF NOT REGEXP_MATCH(smiles, '^[\-\=\#\:\~\.]') THEN
        RETURN false;
    END IF;

    -- Check for valid ring numbers
    IF NOT REGEXP_MATCH(smiles, '^\%?\d{1,2}') THEN
        RETURN false;
    END IF;

    RETURN true;
END;
$$ LANGUAGE plpgsql;

-- Load schema components in dependency order

-- 1. Core functionality
\i database/schema/consolidated/00_core.sql

-- 2. Reference data (needs to be loaded first as other tables reference it)
\i database/schema/consolidated/02_reference_data.sql

-- 3. Compounds (core tables with molecular descriptors)
\i database/schema/consolidated/01_compounds.sql

-- 4. Safety and toxicity data
\i database/schema/consolidated/04_safety.sql

-- 5. Social and community data
\i database/schema/consolidated/03_social.sql

-- 6. Analysis and machine learning
\i database/schema/consolidated/05_analysis.sql

-- 7. Machine learning models
\i database/schema/consolidated/06_ml.sql

-- 8. Web interface data
\i database/schema/consolidated/07_web.sql

-- 9. Quantum chemistry and computational data
\i database/schema/consolidated/08_quantum.sql

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
        AND tablename = 'descriptors_2d'
    ), 'descriptors_2d table not found';

    ASSERT EXISTS (
        SELECT FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename = 'molecular_fingerprints'
    ), 'molecular_fingerprints table not found';

    ASSERT EXISTS (
        SELECT FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename = 'audit_log'
    ), 'audit_log table not found';

    ASSERT EXISTS (
        SELECT FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename = 'quantum_properties'
    ), 'quantum_properties table not found';

    ASSERT EXISTS (
        SELECT FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename = 'electronic_structure'
    ), 'electronic_structure table not found';

    ASSERT EXISTS (
        SELECT FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename = 'quantum_critical_params'
    ), 'quantum_critical_params table not found';

    ASSERT EXISTS (
        SELECT FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename = 'quantum_dynamics'
    ), 'quantum_dynamics table not found';

    ASSERT EXISTS (
        SELECT FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename = 'quantum_basis_sets'
    ), 'quantum_basis_sets table not found';

    ASSERT EXISTS (
        SELECT FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename = 'quantum_functionals'
    ), 'quantum_functionals table not found';

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
        WHERE tgname = 'update_compounds_modtime'
    ), 'compounds update trigger not found';

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
