-- Core PostgreSQL schema
-- This file contains:
-- 1. Required extensions
-- 2. Schema version tracking
-- 3. Core functions
-- 4. Data quality tables and functions
-- 5. Social validation rules and functions
-- 6. Schema verification functions

-- Enable required extensions
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";      -- For UUID generation
CREATE EXTENSION IF NOT EXISTS "unaccent";       -- For text search without accents
CREATE EXTENSION IF NOT EXISTS "pg_trgm";        -- For fuzzy text matching
CREATE EXTENSION IF NOT EXISTS "hstore";         -- For key-value storage
CREATE EXTENSION IF NOT EXISTS "btree_gin";      -- For GIN indexes on B-tree-indexable columns
CREATE EXTENSION IF NOT EXISTS "pg_stat_statements"; -- For query performance monitoring
CREATE EXTENSION IF NOT EXISTS "pgcrypto";       -- For cryptographic functions

-- Set configuration parameters
ALTER DATABASE chemdata SET timezone TO 'UTC';
ALTER DATABASE chemdata SET search_path TO public;

-- Create schema version tracking
CREATE TABLE IF NOT EXISTS schema_versions (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    version text NOT NULL,
    description text,
    applied_at timestamptz NOT NULL DEFAULT now(),
    applied_by text,
    script_name text,
    checksum text,
    execution_time interval,
    status text NOT NULL,
    error_message text
);

-- Create function to update schema version
CREATE OR REPLACE FUNCTION record_schema_version(
    p_version text,
    p_description text,
    p_script_name text,
    p_checksum text
) RETURNS uuid AS $$
DECLARE
    v_start_time timestamptz;
    v_schema_version_id uuid;
BEGIN
    v_start_time := clock_timestamp();
    
    INSERT INTO schema_versions (
        version,
        description,
        script_name,
        checksum,
        applied_by,
        execution_time,
        status
    ) VALUES (
        p_version,
        p_description,
        p_script_name,
        p_checksum,
        current_user,
        clock_timestamp() - v_start_time,
        'SUCCESS'
    ) RETURNING id INTO v_schema_version_id;
    
    RETURN v_schema_version_id;
EXCEPTION WHEN OTHERS THEN
    INSERT INTO schema_versions (
        version,
        description,
        script_name,
        checksum,
        applied_by,
        execution_time,
        status,
        error_message
    ) VALUES (
        p_version,
        p_description,
        p_script_name,
        p_checksum,
        current_user,
        clock_timestamp() - v_start_time,
        'ERROR',
        SQLERRM
    );
    RAISE;
END;
$$ LANGUAGE plpgsql;

-- Create function to check if a schema version exists
CREATE OR REPLACE FUNCTION check_schema_version(p_version text)
RETURNS boolean AS $$
BEGIN
    RETURN EXISTS (
        SELECT 1 
        FROM schema_versions 
        WHERE version = p_version 
        AND status = 'SUCCESS'
    );
END;
$$ LANGUAGE plpgsql;

-- Create function to get current schema version
CREATE OR REPLACE FUNCTION get_current_schema_version()
RETURNS text AS $$
BEGIN
    RETURN version 
    FROM schema_versions 
    WHERE status = 'SUCCESS' 
    ORDER BY applied_at DESC 
    LIMIT 1;
END;
$$ LANGUAGE plpgsql;

-- Create indexes for schema_versions
CREATE INDEX IF NOT EXISTS idx_schema_versions_version ON schema_versions(version);
CREATE INDEX IF NOT EXISTS idx_schema_versions_status ON schema_versions(status);
CREATE INDEX IF NOT EXISTS idx_schema_versions_applied ON schema_versions(applied_at);

-- Create audit log table
CREATE TABLE IF NOT EXISTS audit_log (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    table_name text NOT NULL,
    record_id uuid NOT NULL,
    action text NOT NULL,
    old_data jsonb,
    new_data jsonb,
    changed_by text NOT NULL,
    changed_at timestamptz NOT NULL DEFAULT now()
);

-- Create indexes for audit log
CREATE INDEX IF NOT EXISTS idx_audit_log_table ON audit_log(table_name);
CREATE INDEX IF NOT EXISTS idx_audit_log_record ON audit_log(record_id);
CREATE INDEX IF NOT EXISTS idx_audit_log_action ON audit_log(action);
CREATE INDEX IF NOT EXISTS idx_audit_log_changed ON audit_log(changed_at);

-- Create function for timestamps
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Create function for array validation
CREATE OR REPLACE FUNCTION validate_float_array_min(arr FLOAT[], min_val FLOAT)
RETURNS BOOLEAN AS $$
BEGIN
    -- Handle null array
    IF arr IS NULL THEN
        RETURN TRUE;
    END IF;
    
    -- Check if all values are >= min_val
    RETURN NOT EXISTS (
        SELECT 1
        FROM unnest(arr) AS val
        WHERE val < min_val
    );
END;
$$ LANGUAGE plpgsql;

-- Create function for SMILES validation
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

-- Audit logging function
CREATE OR REPLACE FUNCTION audit_trigger_func()
RETURNS TRIGGER AS $$
DECLARE
    audit_row audit_log;
    excluded_cols text[] = ARRAY[]::text[];
BEGIN
    -- Skip audit logging for the audit_log table itself to prevent recursion
    IF TG_TABLE_NAME = 'audit_log' THEN
        RETURN NULL;
    END IF;

    IF TG_OP = 'INSERT' THEN
        audit_row = ROW(
            uuid_generate_v4(),          -- id
            TG_TABLE_NAME::text,         -- table_name
            NEW.id,                      -- record_id
            'INSERT',                    -- action
            NULL,                        -- old_data
            to_jsonb(NEW),              -- new_data
            current_user,                -- changed_by
            CURRENT_TIMESTAMP           -- changed_at
        );
    ELSIF TG_OP = 'UPDATE' THEN
        audit_row = ROW(
            uuid_generate_v4(),          -- id
            TG_TABLE_NAME::text,         -- table_name
            NEW.id,                      -- record_id
            'UPDATE',                    -- action
            to_jsonb(OLD),              -- old_data
            to_jsonb(NEW),              -- new_data
            current_user,                -- changed_by
            CURRENT_TIMESTAMP           -- changed_at
        );
    ELSIF TG_OP = 'DELETE' THEN
        audit_row = ROW(
            uuid_generate_v4(),          -- id
            TG_TABLE_NAME::text,         -- table_name
            OLD.id,                      -- record_id
            'DELETE',                    -- action
            to_jsonb(OLD),              -- old_data
            NULL,                        -- new_data
            current_user,                -- changed_by
            CURRENT_TIMESTAMP           -- changed_at
        );
    END IF;

    INSERT INTO audit_log VALUES (audit_row.*);
    RETURN NULL;
END;
$$ LANGUAGE plpgsql;

-- Create audit triggers
CREATE TRIGGER audit_schema_versions_trigger
AFTER INSERT OR UPDATE OR DELETE ON schema_versions
FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_audit_log_trigger
AFTER INSERT OR UPDATE OR DELETE ON audit_log
FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Create text search configurations for different content types
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM pg_ts_config WHERE cfgname = 'social_media_search') THEN
        CREATE TEXT SEARCH CONFIGURATION social_media_search (COPY = english);
        ALTER TEXT SEARCH CONFIGURATION social_media_search
            ALTER MAPPING FOR hword, hword_part, word
            WITH english_stem;
    END IF;
END
$$;

-- Schema verification functions
CREATE OR REPLACE FUNCTION verify_required_indexes()
RETURNS TABLE (
    table_name text,
    missing_indexes text[]
) AS $$
DECLARE
    required_indexes jsonb;
    table_record record;
    index_name text;
    missing text[];
BEGIN
    -- Define required indexes for each table
    required_indexes := '{
        "compounds": ["idx_compounds_name", "idx_compounds_smiles", "idx_compounds_cas"],
        "binding_data": ["idx_binding_data_compound", "idx_binding_data_receptor"],
        "social_posts": [
            "idx_social_posts_content_search",
            "idx_reddit_content_search",
            "idx_twitter_content_search",
            "idx_community_content_search",
            "idx_reddit_score",
            "idx_twitter_engagement",
            "idx_community_spam",
            "idx_community_toxicity"
        ],
        "data_quality_metrics": ["idx_quality_metrics_table", "idx_quality_metrics_status"],
        "validation_results": ["idx_validation_results_rule", "idx_validation_results_record"],
        "confidence_scores": ["idx_confidence_scores_table", "idx_confidence_scores_value"]
    }'::jsonb;

    -- Check each table
    FOR table_record IN 
        SELECT key as tname, value as indexes
        FROM jsonb_each(required_indexes)
    LOOP
        missing := ARRAY[]::text[];
        
        -- Check each required index
        FOR index_name IN SELECT jsonb_array_elements_text(table_record.indexes)
        LOOP
            IF NOT EXISTS (
                SELECT 1
                FROM pg_indexes
                WHERE schemaname = 'public'
                AND tablename = table_record.tname
                AND indexname = index_name
            ) THEN
                missing := array_append(missing, index_name);
            END IF;
        END LOOP;

        IF array_length(missing, 1) > 0 THEN
            table_name := table_record.tname;
            missing_indexes := missing;
            RETURN NEXT;
        END IF;
    END LOOP;
END;
$$ LANGUAGE plpgsql;

-- Function to verify text search configurations
CREATE OR REPLACE FUNCTION verify_text_search_configs()
RETURNS TABLE (
    config_name text,
    status text
) AS $$
DECLARE
    required_configs text[];
    config text;
BEGIN
    required_configs := ARRAY['social_media_search'];
    
    FOREACH config IN ARRAY required_configs
    LOOP
        IF EXISTS (
            SELECT 1 
            FROM pg_ts_config 
            WHERE cfgname = config
        ) THEN
            config_name := config;
            status := 'OK';
        ELSE
            config_name := config;
            status := 'MISSING';
        END IF;
        RETURN NEXT;
    END LOOP;
END;
$$ LANGUAGE plpgsql;

-- Function to verify audit triggers
CREATE OR REPLACE FUNCTION verify_audit_triggers()
RETURNS TABLE (
    table_name text,
    trigger_status text
) AS $$
DECLARE
    r RECORD;
BEGIN
    FOR r IN (
        SELECT tablename 
        FROM pg_tables 
        WHERE schemaname = 'public'
    )
    LOOP
        IF EXISTS (
            SELECT 1 
            FROM pg_trigger 
            WHERE tgrelid = (r.tablename::regclass)
            AND tgname LIKE 'audit_%'
        ) THEN
            table_name := r.tablename;
            trigger_status := 'OK';
        ELSE
            table_name := r.tablename;
            trigger_status := 'MISSING AUDIT TRIGGER';
        END IF;
        RETURN NEXT;
    END LOOP;
END;
$$ LANGUAGE plpgsql;

-- Function to verify data quality constraints
CREATE OR REPLACE FUNCTION verify_data_quality_constraints()
RETURNS TABLE (
    result_table_name text,
    constraint_type text,
    status text
) AS $$
DECLARE
    tables_to_check text[];
    t text;
BEGIN
    tables_to_check := ARRAY[
        'social_posts'
    ];
    
    FOREACH t IN ARRAY tables_to_check
    LOOP
        -- Check for NOT NULL constraints
        IF EXISTS (
            SELECT 1
            FROM information_schema.columns c
            WHERE c.table_schema = 'public'
            AND c.table_name = t
            AND c.column_name IN ('content', 'url', 'platform', 'external_id')
            AND c.is_nullable = 'NO'
        ) THEN
            result_table_name := t;
            constraint_type := 'NOT NULL';
            status := 'OK';
        ELSE
            result_table_name := t;
            constraint_type := 'NOT NULL';
            status := 'MISSING';
        END IF;
        RETURN NEXT;

        -- Check for CHECK constraints
        IF EXISTS (
            SELECT 1
            FROM information_schema.check_constraints cc
            JOIN information_schema.constraint_column_usage cu
            ON cc.constraint_name = cu.constraint_name
            WHERE cu.table_schema = 'public'
                AND cu.table_name = 'social_posts'
                AND (
                    cc.check_clause LIKE '%check_engagement_metrics%' OR
                    cc.check_clause LIKE '%check_classification_scores%' OR
                    cc.check_clause LIKE '%check_sentiment_scores%'
                )
        ) THEN
            result_table_name := t;
            constraint_type := 'CHECK';
            status := 'OK';
        ELSE
            result_table_name := t;
            constraint_type := 'CHECK';
            status := 'MISSING';
        END IF;
        RETURN NEXT;
    END LOOP;
END;
$$ LANGUAGE plpgsql;

-- Run all verifications
DO $$
DECLARE
    missing_idx record;
    ts_config record;
    audit_status record;
    constraint_status record;
    has_errors boolean := false;
BEGIN
    RAISE NOTICE 'Starting schema verification...';
    
    -- Check indexes
    FOR missing_idx IN SELECT * FROM verify_required_indexes()
    LOOP
        RAISE WARNING 'Missing indexes for table %: %', 
            missing_idx.table_name, 
            array_to_string(missing_idx.missing_indexes, ', ');
        has_errors := true;
    END LOOP;

    -- Check text search configs
    FOR ts_config IN SELECT * FROM verify_text_search_configs()
    LOOP
        IF ts_config.status = 'MISSING' THEN
            RAISE WARNING 'Missing text search configuration: %', ts_config.config_name;
            has_errors := true;
        END IF;
    END LOOP;

    -- Check audit triggers
    FOR audit_status IN SELECT * FROM verify_audit_triggers()
    LOOP
        IF audit_status.trigger_status = 'MISSING AUDIT TRIGGER' THEN
            RAISE WARNING 'Missing audit trigger for table: %', audit_status.table_name;
            has_errors := true;
        END IF;
    END LOOP;

    -- Check data quality constraints
    FOR constraint_status IN SELECT * FROM verify_data_quality_constraints()
    LOOP
        IF constraint_status.status = 'MISSING' THEN
            RAISE WARNING 'Missing % constraint for table: %', 
                constraint_status.constraint_type, 
                constraint_status.result_table_name;
            has_errors := true;
        END IF;
    END LOOP;

    IF has_errors THEN
        RAISE NOTICE 'Schema verification completed with warnings. See above for details.';
    ELSE
        RAISE NOTICE 'Schema verification completed successfully.';
    END IF;
END;
$$;

-- Record initial schema version
SELECT record_schema_version(
    '0.1.0',
    'Initial core schema setup',
    '00_core.sql',
    NULL
);
