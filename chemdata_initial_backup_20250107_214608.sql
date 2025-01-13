--
-- PostgreSQL database dump
--

-- Dumped from database version 14.15 (Homebrew)
-- Dumped by pg_dump version 14.15 (Homebrew)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: btree_gin; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS btree_gin WITH SCHEMA public;


--
-- Name: EXTENSION btree_gin; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION btree_gin IS 'support for indexing common datatypes in GIN';


--
-- Name: hstore; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS hstore WITH SCHEMA public;


--
-- Name: EXTENSION hstore; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION hstore IS 'data type for storing sets of (key, value) pairs';


--
-- Name: pg_stat_statements; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS pg_stat_statements WITH SCHEMA public;


--
-- Name: EXTENSION pg_stat_statements; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION pg_stat_statements IS 'track planning and execution statistics of all SQL statements executed';


--
-- Name: pg_trgm; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS pg_trgm WITH SCHEMA public;


--
-- Name: EXTENSION pg_trgm; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION pg_trgm IS 'text similarity measurement and index searching based on trigrams';


--
-- Name: pgcrypto; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS pgcrypto WITH SCHEMA public;


--
-- Name: EXTENSION pgcrypto; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION pgcrypto IS 'cryptographic functions';


--
-- Name: unaccent; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS unaccent WITH SCHEMA public;


--
-- Name: EXTENSION unaccent; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION unaccent IS 'text search dictionary that removes accents';


--
-- Name: uuid-ossp; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS "uuid-ossp" WITH SCHEMA public;


--
-- Name: EXTENSION "uuid-ossp"; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION "uuid-ossp" IS 'generate universally unique identifiers (UUIDs)';


--
-- Name: analyze_enaqt_properties(uuid); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.analyze_enaqt_properties(compound_id uuid) RETURNS TABLE(transport_efficiency double precision, coherence_time double precision, decoherence_rate double precision, anti_zeno_factor double precision, is_optimal boolean)
    LANGUAGE plpgsql
    AS $_$
BEGIN
    RETURN QUERY
    WITH dynamics_data AS (
        SELECT 
            qd.coherence_times[1] as coh_time,
            qd.dephasing_rates[1] as deph_rate,
            qd.conductivity_tensor[1] as conductivity,
            els.distribution_type,
            wa.correlation_dimension as d2,
            qh.disorder_strength
        FROM quantum_dynamics qd
        LEFT JOIN energy_level_statistics els ON qd.compound_id = els.compound_id
        LEFT JOIN wavefunction_analysis wa ON qd.compound_id = wa.compound_id
        LEFT JOIN quantum_hamiltonians qh ON qd.compound_id = qh.compound_id
        WHERE qd.compound_id = $1
    )
    SELECT
        conductivity as transport_efficiency,
        coh_time as coherence_time,
        deph_rate as decoherence_rate,
        CASE 
            WHEN distribution_type = 'Semi-Poisson' 
                 AND abs(disorder_strength - 16.5) < 0.5  -- Critical disorder strength
            THEN 1.0
            ELSE deph_rate * coh_time
        END as anti_zeno_factor,
        CASE
            -- Optimal transport conditions from paper:
            -- 1. Semi-Poissonian level statistics (critical)
            -- 2. D2 ≈ 0.5 (multifractal)
            -- 3. W ≈ 16.5 (Anderson transition)
            -- 4. Algebraic coherence decay
            WHEN distribution_type = 'Semi-Poisson' 
                 AND abs(d2 - 0.5) < 0.05
                 AND abs(disorder_strength - 16.5) < 0.5
                 AND deph_rate * coh_time < 1.0  -- Coherent regime
            THEN true
            ELSE false
        END as is_optimal
    FROM dynamics_data;
END;
$_$;


ALTER FUNCTION public.analyze_enaqt_properties(compound_id uuid) OWNER TO armand;

--
-- Name: analyze_level_spacing(uuid, integer); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.analyze_level_spacing(compound_id uuid, min_levels integer DEFAULT 1000) RETURNS TABLE(distribution_type text, gamma_value double precision, confidence_score double precision, chi_squared double precision)
    LANGUAGE plpgsql
    AS $_$
DECLARE
    spacing_count INTEGER;
    mean_spacing FLOAT;
    spacing_variance FLOAT;
BEGIN
    -- Get normalized level spacings
    WITH spacings AS (
        SELECT level_spacings 
        FROM energy_level_statistics 
        WHERE compound_id = $1
    )
    SELECT 
        array_length(level_spacings, 1),
        avg(unnest),
        variance(unnest)
    INTO spacing_count, mean_spacing, spacing_variance
    FROM spacings, unnest(level_spacings);

    -- Check if we have enough statistics
    IF spacing_count < min_levels THEN
        RETURN QUERY SELECT 
            'insufficient_data'::TEXT,
            NULL::FLOAT,
            0::FLOAT,
            NULL::FLOAT;
        RETURN;
    END IF;

    -- Calculate fit scores for each distribution
    RETURN QUERY
    WITH distribution_fits AS (
        SELECT
            -- Poisson: P(s) = exp(-s)
            sum((hist_count - exp(-bin_center))^2) as poisson_chi2,
            -- Wigner: P(s) = (π/2)s*exp(-πs²/4)
            sum((hist_count - (pi()/2)*bin_center*exp(-pi()*bin_center^2/4))^2) as wigner_chi2,
            -- Semi-Poisson: P(s) = 4s*exp(-2s)
            sum((hist_count - 4*bin_center*exp(-2*bin_center))^2) as semi_poisson_chi2
        FROM (
            SELECT 
                width_bucket(unnest(level_spacings), 0, 5, 50) as bin,
                count(*)::float/spacing_count as hist_count,
                (width_bucket(unnest(level_spacings), 0, 5, 50) - 0.5)*0.1 as bin_center
            FROM energy_level_statistics
            WHERE compound_id = $1
            GROUP BY bin
        ) as histogram
    )
    SELECT
        CASE
            WHEN poisson_chi2 < LEAST(wigner_chi2, semi_poisson_chi2) THEN 'Poisson'
            WHEN wigner_chi2 < LEAST(poisson_chi2, semi_poisson_chi2) THEN 'Wigner'
            ELSE 'Semi-Poisson'
        END,
        CASE 
            WHEN semi_poisson_chi2 < LEAST(poisson_chi2, wigner_chi2) THEN 0.0  -- γ ≈ 0 for critical
            ELSE NULL
        END,
        1.0 - LEAST(poisson_chi2, wigner_chi2, semi_poisson_chi2)/
            (poisson_chi2 + wigner_chi2 + semi_poisson_chi2),
        LEAST(poisson_chi2, wigner_chi2, semi_poisson_chi2)
    FROM distribution_fits;
END;
$_$;


ALTER FUNCTION public.analyze_level_spacing(compound_id uuid, min_levels integer) OWNER TO armand;

--
-- Name: analyze_transport_mechanism(uuid); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.analyze_transport_mechanism(compound_id uuid) RETURNS TABLE(mechanism_type text, transport_regime text, efficiency_score double precision, mechanism_confidence double precision)
    LANGUAGE plpgsql
    AS $_$
BEGIN
    RETURN QUERY
    WITH transport_data AS (
        SELECT 
            enaqt.*,
            els.distribution_type,
            wa.correlation_dimension as d2,
            qh.disorder_strength
        FROM analyze_enaqt_properties($1) enaqt
        LEFT JOIN energy_level_statistics els ON els.compound_id = $1
        LEFT JOIN wavefunction_analysis wa ON wa.compound_id = $1
        LEFT JOIN quantum_hamiltonians qh ON qh.compound_id = $1
    )
    SELECT
        -- Classify transport mechanism
        CASE
            WHEN is_optimal THEN 'ENAQT'
            WHEN distribution_type = 'Wigner' THEN 'Quantum'
            ELSE 'Classical'
        END as mechanism_type,
        
        -- Determine transport regime
        CASE
            WHEN distribution_type = 'Semi-Poisson' 
                 AND abs(d2 - 0.5) < 0.05 THEN 'Critical'
            WHEN decoherence_rate * coherence_time < 1.0 THEN 'Coherent'
            ELSE 'Localized'
        END as transport_regime,
        
        -- Calculate efficiency score
        transport_efficiency * 
        CASE
            WHEN is_optimal THEN 2.0  -- Boost score for optimal ENAQT
            ELSE 1.0
        END as efficiency_score,
        
        -- Calculate confidence
        CASE
            WHEN is_optimal THEN 1.0
            WHEN distribution_type = 'Semi-Poisson' 
                 OR abs(d2 - 0.5) < 0.05 THEN 0.8
            ELSE 0.6
        END as mechanism_confidence
    FROM transport_data;
END;
$_$;


ALTER FUNCTION public.analyze_transport_mechanism(compound_id uuid) OWNER TO armand;

--
-- Name: analyze_wavefunction_multifractality(uuid, integer[], integer, integer); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.analyze_wavefunction_multifractality(compound_id uuid, q_values integer[] DEFAULT ARRAY['-10'::integer, '-5'::integer, '-2'::integer, '-1'::integer, 0, 1, 2, 5, 10], min_box_size integer DEFAULT 10, max_box_size integer DEFAULT 1000) RETURNS TABLE(q integer, tau double precision, dimension double precision, r_squared double precision, is_critical boolean)
    LANGUAGE plpgsql
    AS $_$
BEGIN
    RETURN QUERY
    WITH box_counting AS (
        -- Calculate box probabilities μ_k(l) at different scales
        SELECT 
            q_val,
            box_size,
            sum(power(prob, q_val)) as moment
        FROM (
            SELECT 
                q.q as q_val,
                l.l as box_size,
                sum(power(abs(coeff), 2)) as prob
            FROM unnest(q_values) as q(q),
                 generate_series(min_box_size, max_box_size, min_box_size) as l(l),
                 (SELECT wavefunction FROM electronic_structure WHERE compound_id = $1) as wf,
                 jsonb_array_elements(wf) as coeff
            GROUP BY q.q, l.l
        ) as boxes
        GROUP BY q_val, box_size
    ),
    scaling_fits AS (
        -- Fit scaling exponents τ(q)
        SELECT
            q_val,
            regr_slope(ln(moment), ln(box_size)) as tau,
            regr_r2(ln(moment), ln(box_size)) as r2
        FROM box_counting
        GROUP BY q_val
    )
    SELECT 
        q_val,
        tau,
        CASE 
            WHEN q_val = 1 THEN tau
            ELSE tau/(q_val - 1)
        END as dq,
        r2,
        CASE
            WHEN q_val = 2 THEN abs(tau/(q_val - 1) - 0.5) < 0.05  -- D_2 ≈ 0.5 criterion
            ELSE NULL
        END
    FROM scaling_fits
    ORDER BY q_val;
END;
$_$;


ALTER FUNCTION public.analyze_wavefunction_multifractality(compound_id uuid, q_values integer[], min_box_size integer, max_box_size integer) OWNER TO armand;

--
-- Name: audit_trigger_func(); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.audit_trigger_func() RETURNS trigger
    LANGUAGE plpgsql
    AS $$
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
$$;


ALTER FUNCTION public.audit_trigger_func() OWNER TO armand;

--
-- Name: calculate_quantum_criticality_indicators(uuid); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.calculate_quantum_criticality_indicators(compound_id uuid) RETURNS TABLE(indicator_name text, indicator_value double precision, confidence_score double precision)
    LANGUAGE plpgsql
    AS $_$
BEGIN
    RETURN QUERY
    WITH electronic_data AS (
        SELECT 
            homo_lumo_gap,
            total_electronic_energy,
            electron_correlation_energy
        FROM electronic_structure
        WHERE compound_id = $1
    ),
    critical_data AS (
        SELECT 
            correlation_length,
            entanglement_entropy
        FROM quantum_critical_params
        WHERE compound_id = $1
    ),
    level_stats AS (
        SELECT * FROM analyze_level_spacing($1)
    ),
    fractal_analysis AS (
        SELECT * FROM analyze_wavefunction_multifractality($1)
        WHERE q = 2  -- Get correlation dimension D_2
    )
    SELECT 'level_spacing'::TEXT,
           CASE 
               WHEN distribution_type = 'Semi-Poisson' THEN 1.0
               WHEN distribution_type = 'Poisson' THEN 0.0
               ELSE 0.5
           END,
           confidence_score
    FROM level_stats
    UNION ALL
    SELECT 'correlation_dimension',
           dimension,
           r_squared
    FROM fractal_analysis
    UNION ALL
    SELECT 'energy_gap',
           homo_lumo_gap,
           CASE 
               WHEN homo_lumo_gap > 0 THEN 1.0
               ELSE 0.5
           END
    FROM electronic_data
    UNION ALL
    SELECT 'correlation_strength',
           correlation_length,
           CASE 
               WHEN correlation_length > 10 THEN 1.0
               WHEN correlation_length > 5 THEN 0.8
               ELSE 0.6
           END
    FROM critical_data;
END;
$_$;


ALTER FUNCTION public.calculate_quantum_criticality_indicators(compound_id uuid) OWNER TO armand;

--
-- Name: check_schema_version(text); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.check_schema_version(p_version text) RETURNS boolean
    LANGUAGE plpgsql
    AS $$
BEGIN
    RETURN EXISTS (
        SELECT 1 
        FROM schema_versions 
        WHERE version = p_version 
        AND status = 'SUCCESS'
    );
END;
$$;


ALTER FUNCTION public.check_schema_version(p_version text) OWNER TO armand;

--
-- Name: get_current_schema_version(); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.get_current_schema_version() RETURNS text
    LANGUAGE plpgsql
    AS $$
BEGIN
    RETURN version 
    FROM schema_versions 
    WHERE status = 'SUCCESS' 
    ORDER BY applied_at DESC 
    LIMIT 1;
END;
$$;


ALTER FUNCTION public.get_current_schema_version() OWNER TO armand;

--
-- Name: record_schema_version(text, text, text, text); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.record_schema_version(p_version text, p_description text, p_script_name text, p_checksum text) RETURNS uuid
    LANGUAGE plpgsql
    AS $$
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
$$;


ALTER FUNCTION public.record_schema_version(p_version text, p_description text, p_script_name text, p_checksum text) OWNER TO armand;

--
-- Name: update_quantum_properties(); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.update_quantum_properties() RETURNS trigger
    LANGUAGE plpgsql
    AS $$
BEGIN
    -- Update electronic structure when critical parameters change
    IF TG_TABLE_NAME = 'quantum_critical_params' THEN
        UPDATE electronic_structure
        SET updated_at = CURRENT_TIMESTAMP
        WHERE compound_id = NEW.compound_id;
    END IF;
    
    -- Update scaling analysis when dynamics change
    IF TG_TABLE_NAME = 'quantum_dynamics' THEN
        UPDATE scaling_analysis
        SET updated_at = CURRENT_TIMESTAMP
        WHERE compound_id = NEW.compound_id;
    END IF;
    
    RETURN NEW;
END;
$$;


ALTER FUNCTION public.update_quantum_properties() OWNER TO armand;

--
-- Name: update_updated_at_column(); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.update_updated_at_column() RETURNS trigger
    LANGUAGE plpgsql
    AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$;


ALTER FUNCTION public.update_updated_at_column() OWNER TO armand;

--
-- Name: validate_float_array_min(double precision[], double precision); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.validate_float_array_min(arr double precision[], min_val double precision) RETURNS boolean
    LANGUAGE plpgsql
    AS $$
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
$$;


ALTER FUNCTION public.validate_float_array_min(arr double precision[], min_val double precision) OWNER TO armand;

--
-- Name: validate_smiles(text); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.validate_smiles(smiles text) RETURNS boolean
    LANGUAGE plpgsql
    AS $_$
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
$_$;


ALTER FUNCTION public.validate_smiles(smiles text) OWNER TO armand;

--
-- Name: verify_audit_triggers(); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.verify_audit_triggers() RETURNS TABLE(table_name text, trigger_status text)
    LANGUAGE plpgsql
    AS $$
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
$$;


ALTER FUNCTION public.verify_audit_triggers() OWNER TO armand;

--
-- Name: verify_data_quality_constraints(); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.verify_data_quality_constraints() RETURNS TABLE(result_table_name text, constraint_type text, status text)
    LANGUAGE plpgsql
    AS $$
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
$$;


ALTER FUNCTION public.verify_data_quality_constraints() OWNER TO armand;

--
-- Name: verify_required_indexes(); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.verify_required_indexes() RETURNS TABLE(table_name text, missing_indexes text[])
    LANGUAGE plpgsql
    AS $$
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
$$;


ALTER FUNCTION public.verify_required_indexes() OWNER TO armand;

--
-- Name: verify_text_search_configs(); Type: FUNCTION; Schema: public; Owner: armand
--

CREATE FUNCTION public.verify_text_search_configs() RETURNS TABLE(config_name text, status text)
    LANGUAGE plpgsql
    AS $$
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
$$;


ALTER FUNCTION public.verify_text_search_configs() OWNER TO armand;

--
-- Name: social_media_search; Type: TEXT SEARCH CONFIGURATION; Schema: public; Owner: armand
--

CREATE TEXT SEARCH CONFIGURATION public.social_media_search (
    PARSER = pg_catalog."default" );

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR asciiword WITH english_stem;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR word WITH english_stem;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR numword WITH simple;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR email WITH simple;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR url WITH simple;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR host WITH simple;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR sfloat WITH simple;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR version WITH simple;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR hword_numpart WITH simple;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR hword_part WITH english_stem;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR hword_asciipart WITH english_stem;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR numhword WITH simple;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR asciihword WITH english_stem;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR hword WITH english_stem;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR url_path WITH simple;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR file WITH simple;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR "float" WITH simple;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR "int" WITH simple;

ALTER TEXT SEARCH CONFIGURATION public.social_media_search
    ADD MAPPING FOR uint WITH simple;


ALTER TEXT SEARCH CONFIGURATION public.social_media_search OWNER TO armand;

SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- Name: activity_cliffs; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.activity_cliffs (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_pair uuid[] NOT NULL,
    activity_difference double precision,
    structural_similarity double precision,
    cliff_magnitude double precision,
    activity_type text,
    detection_method text,
    significance_score double precision,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.activity_cliffs OWNER TO armand;

--
-- Name: TABLE activity_cliffs; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.activity_cliffs IS 'Activity cliff analysis between compound pairs';


--
-- Name: activity_correlations; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.activity_correlations (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    feature_id uuid NOT NULL,
    activity_type text NOT NULL,
    correlation_coefficient double precision,
    statistical_significance double precision,
    analysis_method text,
    sample_size integer,
    confidence_interval jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.activity_correlations OWNER TO armand;

--
-- Name: TABLE activity_correlations; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.activity_correlations IS 'Correlations between structural features and activities';


--
-- Name: alert_triggers; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.alert_triggers (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    trigger_type text NOT NULL,
    target_type text NOT NULL,
    target_id uuid NOT NULL,
    conditions jsonb NOT NULL,
    actions jsonb NOT NULL,
    severity text,
    priority text,
    notification_channels text[],
    escalation_rules jsonb,
    cooldown_period interval,
    is_active boolean DEFAULT true NOT NULL,
    last_triggered_at timestamp with time zone,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.alert_triggers OWNER TO armand;

--
-- Name: analysis_parameters; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.analysis_parameters (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    parameter_name text NOT NULL,
    description text,
    data_type text NOT NULL,
    validation_rules jsonb,
    default_value jsonb,
    allowed_range jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.analysis_parameters OWNER TO armand;

--
-- Name: api_endpoints; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.api_endpoints (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    path text NOT NULL,
    method text NOT NULL,
    description text,
    parameters jsonb,
    response_schema jsonb,
    auth_required boolean DEFAULT true,
    rate_limit integer,
    cache_ttl interval,
    version text NOT NULL,
    is_deprecated boolean DEFAULT false,
    deprecated_reason text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.api_endpoints OWNER TO armand;

--
-- Name: api_keys; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.api_keys (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    key_hash text NOT NULL,
    user_id text NOT NULL,
    name text NOT NULL,
    permissions text[],
    rate_limit integer,
    expires_at timestamp with time zone,
    last_used_at timestamp with time zone,
    is_active boolean DEFAULT true,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.api_keys OWNER TO armand;

--
-- Name: audit_log; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.audit_log (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    table_name text NOT NULL,
    record_id uuid NOT NULL,
    action text NOT NULL,
    old_data jsonb,
    new_data jsonb,
    changed_by text,
    changed_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.audit_log OWNER TO armand;

--
-- Name: binding_assay_protocols; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.binding_assay_protocols (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    assay_type text NOT NULL,
    protocol_steps jsonb,
    reagents jsonb,
    equipment text[],
    controls jsonb,
    validation_criteria jsonb,
    limitations text[],
    reference_list text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.binding_assay_protocols OWNER TO armand;

--
-- Name: TABLE binding_assay_protocols; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.binding_assay_protocols IS 'Standardized protocols for binding assays';


--
-- Name: binding_assay_types; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.binding_assay_types (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    method_type text,
    detection_type text,
    typical_unit text,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL
);


ALTER TABLE public.binding_assay_types OWNER TO armand;

--
-- Name: TABLE binding_assay_types; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.binding_assay_types IS 'Types of binding assays and their characteristics';


--
-- Name: binding_data; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.binding_data (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    receptor_family_id uuid NOT NULL,
    assay_type_id uuid NOT NULL,
    value double precision NOT NULL,
    unit text NOT NULL,
    confidence_score double precision,
    data_source text,
    publication_doi text,
    experimental_conditions jsonb,
    measurement_date date,
    binding_site text,
    binding_mode text,
    binding_kinetics jsonb,
    activity_type text,
    activity_value double precision,
    activity_unit text,
    efficacy double precision,
    potency double precision,
    assay_description text,
    assay_organism text,
    assay_type text,
    assay_conditions jsonb,
    assay_cell_line text,
    assay_tissue_type text,
    sar_data jsonb,
    pharmacophore_features jsonb,
    binding_pocket_residues text[],
    experimental_method text,
    validation_method text,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    CONSTRAINT binding_data_confidence_score_check CHECK (((confidence_score >= (0)::double precision) AND (confidence_score <= (1)::double precision))),
    CONSTRAINT valid_value CHECK ((value > (0)::double precision))
);


ALTER TABLE public.binding_data OWNER TO armand;

--
-- Name: TABLE binding_data; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.binding_data IS 'Core table for receptor binding measurements';


--
-- Name: binding_data_quality; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.binding_data_quality (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    binding_data_id uuid NOT NULL,
    replicate_count integer,
    standard_deviation double precision,
    confidence_interval double precision,
    quality_score double precision,
    validation_notes text,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    CONSTRAINT binding_data_quality_quality_score_check CHECK (((quality_score >= (0)::double precision) AND (quality_score <= (1)::double precision)))
);


ALTER TABLE public.binding_data_quality OWNER TO armand;

--
-- Name: TABLE binding_data_quality; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.binding_data_quality IS 'Quality metrics for binding measurements';


--
-- Name: binding_kinetics; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.binding_kinetics (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    binding_data_id uuid NOT NULL,
    kon_value double precision,
    kon_unit text,
    koff_value double precision,
    koff_unit text,
    residence_time double precision,
    residence_time_unit text,
    temperature double precision,
    ph double precision,
    ionic_strength text,
    buffer_conditions jsonb,
    method_details text,
    equipment_used text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.binding_kinetics OWNER TO armand;

--
-- Name: TABLE binding_kinetics; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.binding_kinetics IS 'Detailed binding kinetics measurements';


--
-- Name: binding_sar; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.binding_sar (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    receptor_family_id uuid NOT NULL,
    structural_feature text,
    effect_type text,
    effect_magnitude double precision,
    confidence_score double precision,
    evidence_type text,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    CONSTRAINT binding_sar_confidence_score_check CHECK (((confidence_score >= (0)::double precision) AND (confidence_score <= (1)::double precision)))
);


ALTER TABLE public.binding_sar OWNER TO armand;

--
-- Name: TABLE binding_sar; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.binding_sar IS 'Structure-activity relationships for binding data';


--
-- Name: binding_site_mapping; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.binding_site_mapping (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    binding_data_id uuid NOT NULL,
    site_name text NOT NULL,
    residues text[],
    interaction_types text[],
    binding_pocket_volume double precision,
    surface_accessibility double precision,
    conservation_score double precision,
    mutation_effects jsonb,
    structural_features jsonb,
    modeling_method text,
    confidence_score double precision,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.binding_site_mapping OWNER TO armand;

--
-- Name: TABLE binding_site_mapping; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.binding_site_mapping IS 'Binding site characterization and mapping';


--
-- Name: carcinogenicity_data; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.carcinogenicity_data (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    study_type text NOT NULL,
    species text,
    duration interval,
    dose_levels jsonb,
    tumor_types jsonb,
    mechanism text[],
    histopathology jsonb,
    reference_doi text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.carcinogenicity_data OWNER TO armand;

--
-- Name: TABLE carcinogenicity_data; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.carcinogenicity_data IS 'Cancer-related toxicity data';


--
-- Name: cardiotoxicity_data; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.cardiotoxicity_data (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    study_type text NOT NULL,
    cardiac_effects text[],
    mechanism text[],
    herg_ic50 double precision,
    ecg_changes jsonb,
    hemodynamic_effects jsonb,
    reference_doi text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.cardiotoxicity_data OWNER TO armand;

--
-- Name: TABLE cardiotoxicity_data; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.cardiotoxicity_data IS 'Cardiovascular system toxicity data';


--
-- Name: compounds; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.compounds (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    smiles text NOT NULL,
    inchi text,
    inchi_key text,
    cas_number text,
    pubchem_cid text,
    chembl_id text,
    drugbank_id text,
    unii text,
    kegg_id text,
    chemspider_id text,
    zinc_id text,
    chebi_id text,
    iupac_name text,
    preferred_iupac_name text,
    common_names text[],
    einecs_number text,
    rtecs_number text,
    hsdb_number text,
    ccdc_number text,
    reaxys_id text,
    lipidmaps_id text,
    nsc_number text,
    qcarchive_id text,
    nomad_id text,
    materials_project_id text,
    basis_set_id text,
    method_id text,
    calculation_id text,
    molecular_weight double precision,
    molecular_formula text,
    stereochemistry jsonb,
    crystal_structure jsonb,
    solubility_data jsonb,
    pka_values jsonb,
    partition_coefficients jsonb,
    surface_properties jsonb,
    conformational_analysis jsonb,
    chirality_info jsonb,
    isomer_details jsonb,
    melting_point double precision,
    boiling_point double precision,
    density double precision,
    refractive_index double precision,
    optical_rotation double precision,
    drug_class text[],
    mechanism_categories text[],
    therapeutic_categories text[],
    pharmacological_effects text[],
    administration_routes text[],
    bioavailability_data jsonb,
    metabolism_data jsonb,
    distribution_data jsonb,
    legal_status jsonb,
    scheduling_info jsonb,
    approval_status jsonb,
    clinical_trial_status jsonb,
    patent_status jsonb,
    registration_numbers jsonb,
    control_status jsonb,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    CONSTRAINT valid_basis_set CHECK ((basis_set_id ~ '^[A-Z0-9\-]+$'::text)),
    CONSTRAINT valid_calculation CHECK ((calculation_id ~ '^CALC[0-9]+$'::text)),
    CONSTRAINT valid_cas_number CHECK ((cas_number ~ '^[0-9]{1,7}-[0-9]{2}-[0-9]$'::text)),
    CONSTRAINT valid_ccdc CHECK ((ccdc_number ~ '^[0-9]{6}$'::text)),
    CONSTRAINT valid_chebi_id CHECK ((chebi_id ~ '^CHEBI:[0-9]+$'::text)),
    CONSTRAINT valid_chembl_id CHECK ((chembl_id ~ '^CHEMBL[0-9]+$'::text)),
    CONSTRAINT valid_chemspider_id CHECK ((chemspider_id ~ '^[0-9]+$'::text)),
    CONSTRAINT valid_drugbank_id CHECK ((drugbank_id ~ '^DB[0-9]{5}$'::text)),
    CONSTRAINT valid_einecs CHECK ((einecs_number ~ '^[0-9]{3}-[0-9]{3}-[0-9]$'::text)),
    CONSTRAINT valid_hsdb CHECK ((hsdb_number ~ '^[0-9]+$'::text)),
    CONSTRAINT valid_kegg_id CHECK ((kegg_id ~ '^C[0-9]{5}$|^D[0-9]{5}$'::text)),
    CONSTRAINT valid_lipidmaps CHECK ((lipidmaps_id ~ '^LM[A-Z]{2}[0-9]{8}$'::text)),
    CONSTRAINT valid_materials_project CHECK ((materials_project_id ~ '^mp-[0-9]+$'::text)),
    CONSTRAINT valid_method CHECK ((method_id ~ '^[A-Z0-9\-]+$'::text)),
    CONSTRAINT valid_nomad CHECK ((nomad_id ~ '^NOMAD[0-9]+$'::text)),
    CONSTRAINT valid_nsc CHECK ((nsc_number ~ '^NSC[0-9]+$'::text)),
    CONSTRAINT valid_pubchem_cid CHECK ((pubchem_cid ~ '^[0-9]+$'::text)),
    CONSTRAINT valid_qcarchive CHECK ((qcarchive_id ~ '^QCA[0-9]+$'::text)),
    CONSTRAINT valid_reaxys CHECK ((reaxys_id ~ '^RX[0-9]+$'::text)),
    CONSTRAINT valid_rtecs CHECK ((rtecs_number ~ '^[A-Z]{2}[0-9]{4,7}$'::text)),
    CONSTRAINT valid_smiles CHECK (public.validate_smiles(smiles)),
    CONSTRAINT valid_unii CHECK ((unii ~ '^[A-Z0-9]{10}$'::text)),
    CONSTRAINT valid_zinc_id CHECK ((zinc_id ~ '^ZINC[0-9]{12}$'::text))
);


ALTER TABLE public.compounds OWNER TO armand;

--
-- Name: TABLE compounds; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.compounds IS 'Core compound information and properties';


--
-- Name: COLUMN compounds.id; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.id IS 'Unique identifier for each compound';


--
-- Name: COLUMN compounds.name; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.name IS 'Common or systematic name of the compound';


--
-- Name: COLUMN compounds.smiles; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.smiles IS 'SMILES notation representing chemical structure';


--
-- Name: COLUMN compounds.inchi; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.inchi IS 'InChI identifier for the compound';


--
-- Name: COLUMN compounds.inchi_key; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.inchi_key IS 'InChIKey for efficient compound lookup';


--
-- Name: COLUMN compounds.cas_number; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.cas_number IS 'CAS Registry Number';


--
-- Name: COLUMN compounds.pubchem_cid; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.pubchem_cid IS 'PubChem Compound ID';


--
-- Name: COLUMN compounds.chembl_id; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.chembl_id IS 'ChEMBL database identifier';


--
-- Name: COLUMN compounds.drugbank_id; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.drugbank_id IS 'DrugBank database identifier';


--
-- Name: COLUMN compounds.unii; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.unii IS 'FDA Unique Ingredient Identifier';


--
-- Name: COLUMN compounds.kegg_id; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.kegg_id IS 'KEGG database identifier';


--
-- Name: COLUMN compounds.chemspider_id; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.chemspider_id IS 'ChemSpider database identifier';


--
-- Name: COLUMN compounds.zinc_id; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.zinc_id IS 'ZINC database identifier';


--
-- Name: COLUMN compounds.chebi_id; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.chebi_id IS 'ChEBI identifier';


--
-- Name: COLUMN compounds.iupac_name; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.iupac_name IS 'IUPAC systematic name';


--
-- Name: COLUMN compounds.preferred_iupac_name; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.preferred_iupac_name IS 'Preferred IUPAC Name (PIN)';


--
-- Name: COLUMN compounds.common_names; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.common_names IS 'Array of common/trivial names';


--
-- Name: COLUMN compounds.molecular_weight; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.molecular_weight IS 'Molecular weight in g/mol';


--
-- Name: COLUMN compounds.molecular_formula; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.molecular_formula IS 'Molecular formula';


--
-- Name: COLUMN compounds.stereochemistry; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.stereochemistry IS 'Stereochemical configuration details';


--
-- Name: COLUMN compounds.crystal_structure; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.crystal_structure IS 'Crystal structure parameters';


--
-- Name: COLUMN compounds.solubility_data; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.solubility_data IS 'Solubility in various solvents';


--
-- Name: COLUMN compounds.pka_values; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.pka_values IS 'Acid dissociation constants';


--
-- Name: COLUMN compounds.partition_coefficients; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.partition_coefficients IS 'Various partition coefficients';


--
-- Name: COLUMN compounds.surface_properties; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.surface_properties IS 'Surface tension, etc.';


--
-- Name: COLUMN compounds.conformational_analysis; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.conformational_analysis IS 'Conformational states';


--
-- Name: COLUMN compounds.chirality_info; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.chirality_info IS 'Chirality details';


--
-- Name: COLUMN compounds.isomer_details; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.isomer_details IS 'Isomer information';


--
-- Name: COLUMN compounds.drug_class; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.drug_class IS 'Therapeutic/pharmacological classes';


--
-- Name: COLUMN compounds.mechanism_categories; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.mechanism_categories IS 'Mechanism of action categories';


--
-- Name: COLUMN compounds.therapeutic_categories; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.therapeutic_categories IS 'Therapeutic use categories';


--
-- Name: COLUMN compounds.pharmacological_effects; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.pharmacological_effects IS 'Known pharmacological effects';


--
-- Name: COLUMN compounds.administration_routes; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.administration_routes IS 'Routes of administration';


--
-- Name: COLUMN compounds.bioavailability_data; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.bioavailability_data IS 'Bioavailability information';


--
-- Name: COLUMN compounds.metabolism_data; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.metabolism_data IS 'Metabolic pathway information';


--
-- Name: COLUMN compounds.distribution_data; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.distribution_data IS 'Distribution information';


--
-- Name: COLUMN compounds.legal_status; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.legal_status IS 'Legal status by jurisdiction';


--
-- Name: COLUMN compounds.scheduling_info; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.scheduling_info IS 'Drug scheduling information';


--
-- Name: COLUMN compounds.approval_status; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.approval_status IS 'Approval status by region';


--
-- Name: COLUMN compounds.clinical_trial_status; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.clinical_trial_status IS 'Clinical trial information';


--
-- Name: COLUMN compounds.patent_status; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.patent_status IS 'Patent status information';


--
-- Name: COLUMN compounds.registration_numbers; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.registration_numbers IS 'Various registration numbers';


--
-- Name: COLUMN compounds.control_status; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.compounds.control_status IS 'Control status by jurisdiction';


--
-- Name: cytotoxicity_data; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.cytotoxicity_data (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    cell_line text NOT NULL,
    tissue_type text,
    assay_method text,
    exposure_time interval,
    ic50_value double precision,
    ic50_unit text,
    cell_viability double precision,
    cytotoxicity_mechanism text[],
    morphological_changes text[],
    reference_doi text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.cytotoxicity_data OWNER TO armand;

--
-- Name: TABLE cytotoxicity_data; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.cytotoxicity_data IS 'Cell-based toxicity data';


--
-- Name: emergency_procedures; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.emergency_procedures (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    spill_response text[],
    fire_fighting_measures text[],
    first_aid_procedures jsonb,
    evacuation_criteria text[],
    emergency_contacts jsonb,
    special_hazards text[],
    cleanup_procedures text[],
    disposal_procedures text[],
    reporting_requirements text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.emergency_procedures OWNER TO armand;

--
-- Name: TABLE emergency_procedures; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.emergency_procedures IS 'Emergency response and spill control procedures';


--
-- Name: genotoxicity_data; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.genotoxicity_data (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    test_type text NOT NULL,
    organism text,
    metabolic_activation boolean,
    result text NOT NULL,
    mutation_type text[],
    dna_damage_type text[],
    mechanism text[],
    reference_doi text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.genotoxicity_data OWNER TO armand;

--
-- Name: TABLE genotoxicity_data; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.genotoxicity_data IS 'DNA and chromosome damage data';


--
-- Name: hazard_classifications; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.hazard_classifications (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    ghs_classifications text[],
    signal_word text,
    pictograms text[],
    hazard_statements text[],
    precautionary_statements text[],
    nfpa_health integer,
    nfpa_fire integer,
    nfpa_reactivity integer,
    nfpa_special text,
    hmis_health integer,
    hmis_fire integer,
    hmis_physical integer,
    hmis_ppe text,
    classification_source text,
    classification_date date,
    review_date date,
    notes text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL,
    CONSTRAINT valid_hmis_ratings CHECK ((((hmis_health >= 0) AND (hmis_health <= 4)) AND ((hmis_fire >= 0) AND (hmis_fire <= 4)) AND ((hmis_physical >= 0) AND (hmis_physical <= 4)))),
    CONSTRAINT valid_nfpa_ratings CHECK ((((nfpa_health >= 0) AND (nfpa_health <= 4)) AND ((nfpa_fire >= 0) AND (nfpa_fire <= 4)) AND ((nfpa_reactivity >= 0) AND (nfpa_reactivity <= 4))))
);


ALTER TABLE public.hazard_classifications OWNER TO armand;

--
-- Name: TABLE hazard_classifications; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.hazard_classifications IS 'GHS and other hazard classification systems for compounds';


--
-- Name: ppe_requirements; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.ppe_requirements (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    eye_protection text,
    skin_protection text,
    respiratory_protection text,
    hand_protection text,
    body_protection text,
    minimum_ppe_rating text,
    special_requirements text[],
    exposure_limits jsonb,
    monitoring_requirements text[],
    decontamination_procedures text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.ppe_requirements OWNER TO armand;

--
-- Name: TABLE ppe_requirements; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.ppe_requirements IS 'Personal protective equipment and exposure control requirements';


--
-- Name: storage_requirements; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.storage_requirements (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    storage_temp_min numeric,
    storage_temp_max numeric,
    temp_unit text,
    humidity_requirements text,
    light_sensitivity boolean,
    air_sensitivity boolean,
    storage_conditions text[],
    container_type text[],
    incompatible_materials text[],
    segregation_requirements text[],
    ventilation_requirements text,
    static_protection boolean,
    max_storage_time interval,
    storage_precautions text[],
    handling_precautions text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.storage_requirements OWNER TO armand;

--
-- Name: TABLE storage_requirements; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.storage_requirements IS 'Storage conditions and safety requirements';


--
-- Name: toxicity_assays; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.toxicity_assays (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    assay_type text NOT NULL,
    cell_line text,
    organism text,
    endpoint text NOT NULL,
    concentration double precision,
    concentration_unit text,
    exposure_time interval,
    result_value double precision,
    result_unit text,
    result_type text,
    confidence_score double precision,
    protocol_details jsonb,
    reference_doi text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.toxicity_assays OWNER TO armand;

--
-- Name: TABLE toxicity_assays; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.toxicity_assays IS 'General toxicity assay data for compounds';


--
-- Name: compound_safety_overview; Type: VIEW; Schema: public; Owner: armand
--

CREATE VIEW public.compound_safety_overview AS
 SELECT c.id AS compound_id,
    c.name AS compound_name,
    hc.ghs_classifications,
    hc.signal_word,
    hc.pictograms,
    hc.nfpa_health,
    hc.nfpa_fire,
    hc.nfpa_reactivity,
    hc.nfpa_special,
    sr.storage_conditions,
    sr.incompatible_materials,
    pr.minimum_ppe_rating,
    ep.special_hazards,
    count(DISTINCT ta.id) AS toxicity_assay_count,
    count(DISTINCT cd.id) AS cytotoxicity_data_count,
    count(DISTINCT gd.id) AS genotoxicity_data_count
   FROM (((((((public.compounds c
     LEFT JOIN public.hazard_classifications hc ON ((c.id = hc.compound_id)))
     LEFT JOIN public.storage_requirements sr ON ((c.id = sr.compound_id)))
     LEFT JOIN public.ppe_requirements pr ON ((c.id = pr.compound_id)))
     LEFT JOIN public.emergency_procedures ep ON ((c.id = ep.compound_id)))
     LEFT JOIN public.toxicity_assays ta ON ((c.id = ta.compound_id)))
     LEFT JOIN public.cytotoxicity_data cd ON ((c.id = cd.compound_id)))
     LEFT JOIN public.genotoxicity_data gd ON ((c.id = gd.compound_id)))
  GROUP BY c.id, c.name, hc.ghs_classifications, hc.signal_word, hc.pictograms, hc.nfpa_health, hc.nfpa_fire, hc.nfpa_reactivity, hc.nfpa_special, sr.storage_conditions, sr.incompatible_materials, pr.minimum_ppe_rating, ep.special_hazards;


ALTER TABLE public.compound_safety_overview OWNER TO armand;

--
-- Name: descriptors_2d; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.descriptors_2d (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    wiener_index double precision,
    balaban_j_index double precision,
    bertz_ct double precision,
    schultz_molecular_topological_index double precision,
    chi_0_index double precision,
    chi_1_index double precision,
    chi_2_index double precision,
    chi_3_index double precision,
    sum_estate_indices double precision,
    mean_estate_indices double precision,
    kappa_1 double precision,
    kappa_2 double precision,
    kappa_3 double precision,
    molar_refractivity double precision,
    van_der_waals_volume double precision,
    polarizability double precision,
    formal_charge integer,
    ring_count integer,
    aromatic_ring_count integer,
    aliphatic_ring_count integer,
    ring_fusion_degree integer,
    rotatable_bond_count integer,
    rigid_bond_count integer,
    chain_atom_count integer,
    chain_bond_count integer,
    carbon_count integer,
    nitrogen_count integer,
    oxygen_count integer,
    sulfur_count integer,
    phosphorus_count integer,
    halogen_count integer,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    CONSTRAINT valid_counts CHECK (((ring_count >= 0) AND (aromatic_ring_count >= 0) AND (aliphatic_ring_count >= 0) AND (rotatable_bond_count >= 0) AND (chain_atom_count >= 0)))
);


ALTER TABLE public.descriptors_2d OWNER TO armand;

--
-- Name: descriptors_3d; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.descriptors_3d (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    conformer_id integer,
    radius_of_gyration double precision,
    molecular_volume double precision,
    molecular_surface_area double precision,
    solvent_accessible_surface_area double precision,
    polar_surface_area_3d double precision,
    spherosity double precision,
    asphericity double precision,
    eccentricity double precision,
    inertial_shape_factor double precision,
    principal_moment_1 double precision,
    principal_moment_2 double precision,
    principal_moment_3 double precision,
    gravitational_index double precision,
    radius_of_distribution double precision,
    molecular_surface_potential double precision,
    average_surface_charge double precision,
    total_energy double precision,
    strain_energy double precision,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    CONSTRAINT valid_3d_ranges CHECK (((molecular_volume > (0)::double precision) AND (molecular_surface_area > (0)::double precision) AND (solvent_accessible_surface_area > (0)::double precision)))
);


ALTER TABLE public.descriptors_3d OWNER TO armand;

--
-- Name: developmental_toxicity; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.developmental_toxicity (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    study_type text NOT NULL,
    species text,
    exposure_period text,
    dose_levels jsonb,
    developmental_effects jsonb,
    mechanism text[],
    teratogenicity boolean,
    reference_doi text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.developmental_toxicity OWNER TO armand;

--
-- Name: TABLE developmental_toxicity; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.developmental_toxicity IS 'Developmental and reproductive toxicity data';


--
-- Name: device_settings; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.device_settings (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    device_id text NOT NULL,
    device_type text NOT NULL,
    os_version text,
    app_version text,
    screen_size text,
    capabilities jsonb,
    settings jsonb,
    performance_profile jsonb,
    security_settings jsonb,
    network_preferences jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.device_settings OWNER TO armand;

--
-- Name: effect_relationships; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.effect_relationships (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    effect_id uuid NOT NULL,
    related_effect_id uuid NOT NULL,
    relationship_type text NOT NULL,
    strength double precision,
    evidence_level text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.effect_relationships OWNER TO armand;

--
-- Name: electronic_structure; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.electronic_structure (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    total_electronic_energy double precision,
    electron_correlation_energy double precision,
    exchange_energy double precision,
    kinetic_energy double precision,
    homo_lumo_gap double precision,
    orbital_energies double precision[],
    orbital_occupancies integer[],
    density_matrix jsonb,
    density_grid_points jsonb,
    density_values double precision[],
    wavefunction jsonb,
    state_type character varying(50),
    band_structure jsonb,
    method text,
    basis_set text,
    convergence_criteria jsonb,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    metadata jsonb,
    CONSTRAINT valid_energies CHECK (((total_electronic_energy < (0)::double precision) AND (homo_lumo_gap >= (0)::double precision)))
);


ALTER TABLE public.electronic_structure OWNER TO armand;

--
-- Name: TABLE electronic_structure; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.electronic_structure IS 'Quantum mechanical electronic structure data';


--
-- Name: COLUMN electronic_structure.density_matrix; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.electronic_structure.density_matrix IS 'Electron density matrix stored as sparse matrix in JSONB format';


--
-- Name: COLUMN electronic_structure.wavefunction; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.electronic_structure.wavefunction IS 'Wavefunction data stored as complex numbers in JSONB format';


--
-- Name: energy_level_statistics; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.energy_level_statistics (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    energy_levels double precision[],
    level_spacings double precision[],
    spacing_distribution jsonb,
    cumulative_distribution jsonb,
    distribution_type text,
    gamma_parameter double precision,
    confidence_score double precision,
    analysis_parameters jsonb,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    CONSTRAINT valid_gamma CHECK ((gamma_parameter >= (0)::double precision))
);


ALTER TABLE public.energy_level_statistics OWNER TO armand;

--
-- Name: experience_categories; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.experience_categories (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    category_name text NOT NULL,
    parent_category text,
    description text,
    report_count integer DEFAULT 0 NOT NULL,
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.experience_categories OWNER TO armand;

--
-- Name: external_services; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.external_services (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    service_type text NOT NULL,
    base_url text NOT NULL,
    auth_config jsonb,
    rate_limit jsonb,
    timeout interval,
    retry_config jsonb,
    is_active boolean DEFAULT true,
    health_check_path text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.external_services OWNER TO armand;

--
-- Name: feature_definitions; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.feature_definitions (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    feature_type text NOT NULL,
    data_type text NOT NULL,
    calculation_method text,
    dependencies text[],
    validation_rules jsonb,
    metadata jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.feature_definitions OWNER TO armand;

--
-- Name: feature_pipelines; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.feature_pipelines (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    steps jsonb NOT NULL,
    input_features text[],
    output_features text[],
    parameters jsonb,
    validation_rules jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.feature_pipelines OWNER TO armand;

--
-- Name: feature_selection; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.feature_selection (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    model_id uuid NOT NULL,
    selection_method text NOT NULL,
    selected_features text[] NOT NULL,
    importance_scores jsonb,
    selection_criteria jsonb,
    validation_metrics jsonb,
    selection_date timestamp with time zone NOT NULL,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.feature_selection OWNER TO armand;

--
-- Name: feature_values; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.feature_values (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    feature_id uuid NOT NULL,
    value jsonb NOT NULL,
    calculation_date timestamp with time zone NOT NULL,
    confidence_score double precision,
    metadata jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.feature_values OWNER TO armand;

--
-- Name: genes; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.genes (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    symbol text NOT NULL,
    description text,
    organism text NOT NULL,
    chromosome text,
    location text,
    gene_type text,
    sequence text,
    alternative_symbols text[],
    ensembl_id text,
    entrez_id text,
    hgnc_id text,
    mgi_id text,
    rgd_id text,
    uniprot_ids text[],
    refseq_ids text[],
    regulatory_elements jsonb,
    expression_data jsonb,
    pathway_involvement text[],
    disease_associations jsonb,
    gene_ontology jsonb,
    evolutionary_conservation jsonb,
    variants jsonb,
    interactions text[],
    literature_references text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL,
    CONSTRAINT valid_gene_symbol CHECK ((symbol ~ '^[A-Za-z0-9-]+$'::text))
);


ALTER TABLE public.genes OWNER TO armand;

--
-- Name: TABLE genes; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.genes IS 'Gene information and annotations';


--
-- Name: hepatotoxicity_data; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.hepatotoxicity_data (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    study_type text NOT NULL,
    liver_effects text[],
    mechanism text[],
    enzyme_elevations jsonb,
    histopathology jsonb,
    clinical_significance text,
    reference_doi text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.hepatotoxicity_data OWNER TO armand;

--
-- Name: TABLE hepatotoxicity_data; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.hepatotoxicity_data IS 'Liver toxicity data';


--
-- Name: immunotoxicity_data; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.immunotoxicity_data (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    study_type text NOT NULL,
    immune_parameters text[],
    effect_type text,
    mechanism text[],
    clinical_relevance text,
    reference_doi text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.immunotoxicity_data OWNER TO armand;

--
-- Name: TABLE immunotoxicity_data; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.immunotoxicity_data IS 'Immune system toxicity data';


--
-- Name: literature_findings; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.literature_findings (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    paper_id uuid NOT NULL,
    compound_id uuid NOT NULL,
    finding_type text NOT NULL,
    finding_details text NOT NULL,
    statistical_significance double precision,
    confidence_interval jsonb,
    methodology_notes text,
    limitations text[],
    replication_status text,
    validation_method text[],
    supporting_evidence jsonb,
    contradicting_evidence jsonb,
    clinical_relevance text,
    research_implications text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.literature_findings OWNER TO armand;

--
-- Name: TABLE literature_findings; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.literature_findings IS 'Specific findings from scientific literature';


--
-- Name: meta_analyses; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.meta_analyses (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    title text NOT NULL,
    topic text NOT NULL,
    included_papers uuid[] NOT NULL,
    methodology text NOT NULL,
    total_sample_size integer,
    pooled_effect_size double precision,
    heterogeneity_metrics jsonb,
    subgroup_analyses jsonb,
    sensitivity_analyses jsonb,
    publication_bias_assessment jsonb,
    quality_assessment_method text,
    evidence_strength text,
    clinical_implications text[],
    research_gaps text[],
    conclusions text,
    limitations text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.meta_analyses OWNER TO armand;

--
-- Name: TABLE meta_analyses; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.meta_analyses IS 'Meta-analyses of multiple research papers';


--
-- Name: metabolic_toxicity; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.metabolic_toxicity (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    enzyme_affected text[],
    inhibition_type text,
    ki_value double precision,
    ki_unit text,
    metabolites jsonb,
    pathway_disruption text[],
    clinical_significance text,
    reference_doi text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.metabolic_toxicity OWNER TO armand;

--
-- Name: TABLE metabolic_toxicity; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.metabolic_toxicity IS 'Metabolic enzyme interactions and toxicity';


--
-- Name: ml_models; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.ml_models (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    version text NOT NULL,
    model_type text NOT NULL,
    description text,
    target_variable text NOT NULL,
    features text[] NOT NULL,
    hyperparameters jsonb NOT NULL,
    architecture jsonb,
    preprocessing_steps jsonb,
    training_config jsonb,
    metadata jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.ml_models OWNER TO armand;

--
-- Name: model_deployments; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.model_deployments (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    model_id uuid NOT NULL,
    deployment_name text NOT NULL,
    deployment_environment text NOT NULL,
    deployment_date timestamp with time zone NOT NULL,
    status text NOT NULL,
    version_tag text NOT NULL,
    configuration jsonb,
    performance_metrics jsonb,
    monitoring_config jsonb,
    rollback_info jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.model_deployments OWNER TO armand;

--
-- Name: model_metrics; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.model_metrics (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    model_id uuid NOT NULL,
    version_id uuid NOT NULL,
    metric_type text NOT NULL,
    metric_value double precision NOT NULL,
    metric_date timestamp with time zone NOT NULL,
    dataset_info jsonb,
    calculation_method text,
    confidence_interval jsonb,
    metadata jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.model_metrics OWNER TO armand;

--
-- Name: model_monitoring; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.model_monitoring (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    model_id uuid NOT NULL,
    monitoring_date timestamp with time zone NOT NULL,
    metric_name text NOT NULL,
    metric_value double precision NOT NULL,
    threshold_value double precision,
    alert_status text,
    data_drift_metrics jsonb,
    performance_metrics jsonb,
    resource_usage jsonb,
    alert_history jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.model_monitoring OWNER TO armand;

--
-- Name: model_predictions; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.model_predictions (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    model_id uuid NOT NULL,
    version_id uuid NOT NULL,
    compound_id uuid NOT NULL,
    prediction_type text NOT NULL,
    predicted_value jsonb NOT NULL,
    confidence_score double precision,
    prediction_date timestamp with time zone NOT NULL,
    input_features jsonb,
    explanation jsonb,
    metadata jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.model_predictions OWNER TO armand;

--
-- Name: model_validation_results; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.model_validation_results (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    model_id uuid NOT NULL,
    validation_type text NOT NULL,
    validation_date timestamp with time zone NOT NULL,
    test_dataset_info jsonb,
    validation_metrics jsonb,
    test_set_performance jsonb,
    cross_validation_results jsonb,
    error_analysis jsonb,
    validation_plots jsonb,
    validation_notes text,
    recommendations text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.model_validation_results OWNER TO armand;

--
-- Name: model_versions; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.model_versions (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    model_id uuid NOT NULL,
    version_number text NOT NULL,
    changes_description text,
    performance_metrics jsonb,
    validation_results jsonb,
    deployment_status text,
    deployed_at timestamp with time zone,
    deprecated_at timestamp with time zone,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.model_versions OWNER TO armand;

--
-- Name: molecular_descriptors; Type: VIEW; Schema: public; Owner: armand
--

CREATE VIEW public.molecular_descriptors AS
 SELECT descriptors_2d.id,
    descriptors_2d.compound_id,
    descriptors_2d.wiener_index,
    descriptors_2d.balaban_j_index,
    descriptors_2d.bertz_ct,
    descriptors_2d.schultz_molecular_topological_index,
    descriptors_2d.chi_0_index,
    descriptors_2d.chi_1_index,
    descriptors_2d.chi_2_index,
    descriptors_2d.chi_3_index,
    descriptors_2d.sum_estate_indices,
    descriptors_2d.mean_estate_indices,
    descriptors_2d.kappa_1,
    descriptors_2d.kappa_2,
    descriptors_2d.kappa_3,
    descriptors_2d.molar_refractivity,
    descriptors_2d.van_der_waals_volume,
    descriptors_2d.polarizability,
    descriptors_2d.formal_charge,
    descriptors_2d.ring_count,
    descriptors_2d.aromatic_ring_count,
    descriptors_2d.aliphatic_ring_count,
    descriptors_2d.ring_fusion_degree,
    descriptors_2d.rotatable_bond_count,
    descriptors_2d.rigid_bond_count,
    descriptors_2d.chain_atom_count,
    descriptors_2d.chain_bond_count,
    descriptors_2d.carbon_count,
    descriptors_2d.nitrogen_count,
    descriptors_2d.oxygen_count,
    descriptors_2d.sulfur_count,
    descriptors_2d.phosphorus_count,
    descriptors_2d.halogen_count,
    descriptors_2d.created_at,
    descriptors_2d.updated_at
   FROM public.descriptors_2d;


ALTER TABLE public.molecular_descriptors OWNER TO armand;

--
-- Name: molecular_fingerprints; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.molecular_fingerprints (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    maccs_keys bit varying(166),
    pubchem_bits bit varying(881),
    ecfp4_bits bit varying(1024),
    ecfp6_bits bit varying(1024),
    daylight_bits bit varying(1024),
    pharma_bits bit varying(512),
    atom_pairs bit varying(1024),
    torsion_bits bit varying(1024),
    morgan_bits bit varying(1024),
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL
);


ALTER TABLE public.molecular_fingerprints OWNER TO armand;

--
-- Name: monitoring_parameters; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.monitoring_parameters (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    category text NOT NULL,
    frequency text NOT NULL,
    monitoring_method text NOT NULL,
    alert_conditions text[],
    normal_range jsonb,
    data_type text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.monitoring_parameters OWNER TO armand;

--
-- Name: TABLE monitoring_parameters; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.monitoring_parameters IS 'Parameters to monitor for safety assessment';


--
-- Name: COLUMN monitoring_parameters.name; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.monitoring_parameters.name IS 'Name of the parameter to monitor';


--
-- Name: COLUMN monitoring_parameters.category; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.monitoring_parameters.category IS 'Category of the monitoring parameter';


--
-- Name: COLUMN monitoring_parameters.frequency; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.monitoring_parameters.frequency IS 'How often the parameter should be monitored';


--
-- Name: COLUMN monitoring_parameters.monitoring_method; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.monitoring_parameters.monitoring_method IS 'Method used to monitor this parameter';


--
-- Name: COLUMN monitoring_parameters.alert_conditions; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.monitoring_parameters.alert_conditions IS 'Conditions that should trigger alerts';


--
-- Name: neurotoxicity_data; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.neurotoxicity_data (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    study_type text NOT NULL,
    brain_regions_affected text[],
    behavioral_effects text[],
    cellular_effects text[],
    mechanism text[],
    reversibility text,
    long_term_effects jsonb,
    reference_doi text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.neurotoxicity_data OWNER TO armand;

--
-- Name: TABLE neurotoxicity_data; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.neurotoxicity_data IS 'Nervous system toxicity data';


--
-- Name: offline_data; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.offline_data (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    device_id text NOT NULL,
    data_type text NOT NULL,
    data_id uuid NOT NULL,
    content jsonb NOT NULL,
    version integer NOT NULL,
    priority integer DEFAULT 0,
    compression_type text,
    encryption_type text,
    validation_hash text,
    expires_at timestamp with time zone,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.offline_data OWNER TO armand;

--
-- Name: organ_systems; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.organ_systems (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    major_components text[],
    key_functions text[],
    vulnerability_factors text[],
    assessment_methods text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.organ_systems OWNER TO armand;

--
-- Name: organ_toxicity_patterns; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.organ_toxicity_patterns (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    organ_system_id uuid NOT NULL,
    name text NOT NULL,
    description text,
    cellular_targets text[],
    molecular_mechanisms text[],
    histological_changes text[],
    functional_impacts text[],
    early_biomarkers text[],
    diagnostic_markers text[],
    progression_pattern text,
    reversibility_potential text,
    risk_factors text[],
    protective_factors text[],
    monitoring_parameters jsonb,
    intervention_thresholds jsonb,
    treatment_approaches text[],
    prevention_strategies text[],
    research_status text,
    evidence_level text,
    notes text,
    text_search_vector tsvector GENERATED ALWAYS AS ((setweight(to_tsvector('english'::regconfig, COALESCE(name, ''::text)), 'A'::"char") || setweight(to_tsvector('english'::regconfig, COALESCE(description, ''::text)), 'B'::"char"))) STORED,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.organ_toxicity_patterns OWNER TO armand;

--
-- Name: pathway_analysis; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.pathway_analysis (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    pathway_name text NOT NULL,
    pathway_type text NOT NULL,
    affected_proteins text[],
    regulation_effects jsonb,
    downstream_effects jsonb,
    feedback_mechanisms jsonb,
    pathway_crosstalk jsonb,
    temporal_dynamics jsonb,
    tissue_specificity jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.pathway_analysis OWNER TO armand;

--
-- Name: TABLE pathway_analysis; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.pathway_analysis IS 'Analysis of pathway effects and regulation';


--
-- Name: performance_metrics; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.performance_metrics (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    metric_type text NOT NULL,
    component text NOT NULL,
    value double precision NOT NULL,
    unit text,
    threshold double precision,
    status text,
    metadata jsonb,
    tags text[],
    alert_triggered boolean DEFAULT false,
    resolution_steps jsonb,
    measured_at timestamp with time zone DEFAULT now() NOT NULL,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.performance_metrics OWNER TO armand;

--
-- Name: pharmacological_classes; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.pharmacological_classes (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    mechanism_type text NOT NULL,
    target_systems text[],
    typical_effects text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.pharmacological_classes OWNER TO armand;

--
-- Name: pharmacophore_features; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.pharmacophore_features (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    h_bond_donors_count integer,
    h_bond_acceptors_count integer,
    pos_charge_groups_count integer,
    neg_charge_groups_count integer,
    aromatic_rings_count integer,
    hydrophobic_groups_count integer,
    donor_positions jsonb,
    acceptor_positions jsonb,
    charge_positions jsonb,
    aromatic_positions jsonb,
    hydrophobic_positions jsonb,
    donor_strengths double precision[],
    acceptor_strengths double precision[],
    charge_strengths double precision[],
    feature_type text NOT NULL,
    coordinates jsonb,
    strength double precision,
    interaction_radius double precision,
    optional boolean,
    detection_method text,
    confidence_score double precision,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    CONSTRAINT valid_feature_counts CHECK (((h_bond_donors_count >= 0) AND (h_bond_acceptors_count >= 0) AND (pos_charge_groups_count >= 0) AND (neg_charge_groups_count >= 0) AND (aromatic_rings_count >= 0) AND (hydrophobic_groups_count >= 0)))
);


ALTER TABLE public.pharmacophore_features OWNER TO armand;

--
-- Name: phase_transition_types; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.phase_transition_types (
    id integer NOT NULL,
    name text NOT NULL,
    description text,
    transition_order integer,
    critical_exponents jsonb,
    universality_class text,
    characteristic_properties text[]
);


ALTER TABLE public.phase_transition_types OWNER TO armand;

--
-- Name: TABLE phase_transition_types; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.phase_transition_types IS 'Classification of phase transitions and their properties';


--
-- Name: phase_transition_types_id_seq; Type: SEQUENCE; Schema: public; Owner: armand
--

CREATE SEQUENCE public.phase_transition_types_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.phase_transition_types_id_seq OWNER TO armand;

--
-- Name: phase_transition_types_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: armand
--

ALTER SEQUENCE public.phase_transition_types_id_seq OWNED BY public.phase_transition_types.id;


--
-- Name: phase_transitions; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.phase_transitions (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    transition_type character varying(50) NOT NULL,
    critical_temperature double precision,
    critical_pressure double precision,
    order_parameter jsonb,
    correlation_length double precision,
    transition_order integer,
    hysteresis_data jsonb,
    fluctuation_data jsonb,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    metadata jsonb,
    CONSTRAINT valid_transition_params CHECK (((critical_temperature >= (0)::double precision) AND (critical_pressure >= (0)::double precision) AND (correlation_length > (0)::double precision) AND (transition_order > 0)))
);


ALTER TABLE public.phase_transitions OWNER TO armand;

--
-- Name: TABLE phase_transitions; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.phase_transitions IS 'Phase transition data and properties';


--
-- Name: COLUMN phase_transitions.order_parameter; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.phase_transitions.order_parameter IS 'Order parameter data stored in JSONB format';


--
-- Name: proteins; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.proteins (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    gene_id uuid,
    name text NOT NULL,
    symbol text NOT NULL,
    description text,
    organism text NOT NULL,
    sequence text,
    length integer,
    molecular_weight double precision,
    uniprot_id text,
    pdb_ids text[],
    refseq_ids text[],
    alternative_names text[],
    protein_family text,
    domains jsonb,
    motifs jsonb,
    subcellular_location text[],
    post_translational_modifications jsonb,
    structure_data jsonb,
    function_data jsonb,
    interactions jsonb,
    expression_pattern jsonb,
    regulatory_mechanisms jsonb,
    disease_associations jsonb,
    drug_interactions text[],
    pathway_involvement text[],
    literature_references text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.proteins OWNER TO armand;

--
-- Name: TABLE proteins; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.proteins IS 'Protein structures and annotations';


--
-- Name: quantum_basis_sets; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.quantum_basis_sets (
    id integer NOT NULL,
    name text NOT NULL,
    description text,
    basis_type text,
    elements text[],
    accuracy_level text,
    computational_cost text,
    reference_citation text
);


ALTER TABLE public.quantum_basis_sets OWNER TO armand;

--
-- Name: TABLE quantum_basis_sets; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.quantum_basis_sets IS 'Standard basis sets for quantum calculations';


--
-- Name: quantum_basis_sets_id_seq; Type: SEQUENCE; Schema: public; Owner: armand
--

CREATE SEQUENCE public.quantum_basis_sets_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.quantum_basis_sets_id_seq OWNER TO armand;

--
-- Name: quantum_basis_sets_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: armand
--

ALTER SEQUENCE public.quantum_basis_sets_id_seq OWNED BY public.quantum_basis_sets.id;


--
-- Name: quantum_calculation_summary; Type: VIEW; Schema: public; Owner: armand
--

CREATE VIEW public.quantum_calculation_summary AS
SELECT
    NULL::uuid AS id,
    NULL::uuid AS compound_id,
    NULL::text AS calculation_type,
    NULL::text AS basis_set,
    NULL::text AS functional,
    NULL::double precision AS energy_hartree,
    NULL::boolean AS convergence_achieved,
    NULL::bigint AS finding_count,
    NULL::text AS finding_types;


ALTER TABLE public.quantum_calculation_summary OWNER TO armand;

--
-- Name: quantum_calculations; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.quantum_calculations (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    project_id uuid,
    compound_id uuid,
    calculation_type text,
    basis_set text,
    functional text,
    calculation_status text,
    start_timestamp timestamp without time zone,
    end_timestamp timestamp without time zone,
    cpu_hours double precision,
    memory_gb double precision,
    convergence_achieved boolean,
    energy_hartree double precision,
    energy_gradient_norm double precision,
    calculation_parameters jsonb,
    output_files text[],
    metadata jsonb,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    CONSTRAINT check_energy_gradient CHECK ((energy_gradient_norm >= (0)::double precision)),
    CONSTRAINT check_energy_hartree CHECK ((energy_hartree < (0)::double precision)),
    CONSTRAINT check_resources CHECK (((cpu_hours > (0)::double precision) AND (memory_gb > (0)::double precision)))
);


ALTER TABLE public.quantum_calculations OWNER TO armand;

--
-- Name: TABLE quantum_calculations; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.quantum_calculations IS 'Quantum mechanical calculations and their parameters';


--
-- Name: quantum_critical_params; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.quantum_critical_params (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    critical_temperature double precision,
    critical_pressure double precision,
    critical_field double precision,
    primary_order_parameter text,
    order_parameter_values jsonb,
    alpha double precision,
    beta double precision,
    gamma double precision,
    delta double precision,
    nu double precision,
    eta double precision,
    correlation_length double precision,
    correlation_function jsonb,
    dynamic_exponent_z double precision,
    phase_boundaries jsonb,
    multicriticality_type text,
    coherence_length double precision,
    entanglement_entropy double precision,
    quantum_fluctuations jsonb,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    metadata jsonb,
    CONSTRAINT valid_critical_params CHECK (((critical_temperature >= (0)::double precision) AND (critical_pressure >= (0)::double precision) AND (correlation_length > (0)::double precision)))
);


ALTER TABLE public.quantum_critical_params OWNER TO armand;

--
-- Name: TABLE quantum_critical_params; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.quantum_critical_params IS 'Parameters characterizing quantum critical behavior';


--
-- Name: COLUMN quantum_critical_params.order_parameter_values; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.quantum_critical_params.order_parameter_values IS 'Order parameter data at different conditions stored in JSONB format';


--
-- Name: quantum_hamiltonians; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.quantum_hamiltonians (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    eh_hamiltonian jsonb,
    overlap_matrix jsonb,
    transformed_hamiltonian jsonb,
    disorder_strength double precision,
    is_critical boolean,
    calculation_method text,
    calculation_parameters jsonb,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL
);


ALTER TABLE public.quantum_hamiltonians OWNER TO armand;

--
-- Name: wavefunction_analysis; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.wavefunction_analysis (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    box_probabilities jsonb,
    scaling_exponents jsonb,
    fractal_dimensions jsonb,
    correlation_dimension double precision,
    localization_length double precision,
    participation_ratio double precision,
    analysis_parameters jsonb,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    CONSTRAINT valid_dimensions CHECK (((correlation_dimension >= (0)::double precision) AND (correlation_dimension <= (1)::double precision)))
);


ALTER TABLE public.wavefunction_analysis OWNER TO armand;

--
-- Name: quantum_criticality_summary; Type: VIEW; Schema: public; Owner: armand
--

CREATE VIEW public.quantum_criticality_summary AS
 SELECT c.id AS compound_id,
    c.name AS compound_name,
    els.distribution_type,
    els.gamma_parameter,
    els.confidence_score AS spacing_confidence,
    wa.correlation_dimension AS d2,
    wa.participation_ratio,
    qh.disorder_strength,
    qh.is_critical AS hamiltonian_critical,
        CASE
            WHEN ((els.distribution_type = 'Semi-Poisson'::text) AND (abs((wa.correlation_dimension - (0.5)::double precision)) < (0.05)::double precision) AND (qh.is_critical = true)) THEN true
            ELSE false
        END AS is_critical,
    (((els.confidence_score + (
        CASE
            WHEN (abs((wa.correlation_dimension - (0.5)::double precision)) < (0.05)::double precision) THEN 1.0
            ELSE 0.0
        END)::double precision) + (
        CASE
            WHEN qh.is_critical THEN 1.0
            ELSE 0.0
        END)::double precision) / (3.0)::double precision) AS criticality_confidence
   FROM (((public.compounds c
     LEFT JOIN public.energy_level_statistics els ON ((c.id = els.compound_id)))
     LEFT JOIN public.wavefunction_analysis wa ON ((c.id = wa.compound_id)))
     LEFT JOIN public.quantum_hamiltonians qh ON ((c.id = qh.compound_id)));


ALTER TABLE public.quantum_criticality_summary OWNER TO armand;

--
-- Name: VIEW quantum_criticality_summary; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON VIEW public.quantum_criticality_summary IS 'Consolidated view of quantum criticality indicators for compounds based on:
1. Level spacing distribution (should be Semi-Poissonian)
2. Correlation dimension D2 ≈ 0.5 (multifractal criterion)
3. Critical disorder strength at metal-insulator transition';


--
-- Name: quantum_dynamics; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.quantum_dynamics (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    time_points double precision[],
    wavefunction_evolution jsonb,
    density_matrix_evolution jsonb,
    coherence_times double precision[],
    relaxation_rates double precision[],
    dephasing_rates double precision[],
    conductivity_tensor double precision[],
    hall_conductance double precision,
    thermal_conductivity double precision,
    spectral_function jsonb,
    optical_conductivity jsonb,
    entanglement_spectrum double precision[],
    mutual_information double precision,
    dissipation_kernel jsonb,
    noise_spectrum jsonb,
    metadata jsonb,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL
);


ALTER TABLE public.quantum_dynamics OWNER TO armand;

--
-- Name: TABLE quantum_dynamics; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.quantum_dynamics IS 'Time-dependent quantum properties and dynamics';


--
-- Name: quantum_functionals; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.quantum_functionals (
    id integer NOT NULL,
    name text NOT NULL,
    description text,
    functional_type text,
    properties_handled text[],
    accuracy_metrics jsonb,
    computational_cost text,
    reference_citation text
);


ALTER TABLE public.quantum_functionals OWNER TO armand;

--
-- Name: TABLE quantum_functionals; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.quantum_functionals IS 'Density functionals for electronic structure calculations';


--
-- Name: quantum_functionals_id_seq; Type: SEQUENCE; Schema: public; Owner: armand
--

CREATE SEQUENCE public.quantum_functionals_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.quantum_functionals_id_seq OWNER TO armand;

--
-- Name: quantum_functionals_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: armand
--

ALTER SEQUENCE public.quantum_functionals_id_seq OWNED BY public.quantum_functionals.id;


--
-- Name: quantum_observables; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.quantum_observables (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    observable_type character varying(50) NOT NULL,
    value double precision,
    uncertainty double precision,
    measurement_basis text,
    operator_type text,
    expectation_value double precision,
    variance double precision,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    metadata jsonb,
    CONSTRAINT valid_uncertainty CHECK ((uncertainty >= (0)::double precision))
);


ALTER TABLE public.quantum_observables OWNER TO armand;

--
-- Name: TABLE quantum_observables; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.quantum_observables IS 'Quantum observable measurements and uncertainties';


--
-- Name: quantum_observables_ref; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.quantum_observables_ref (
    id integer NOT NULL,
    name text NOT NULL,
    description text,
    operator_type text,
    measurement_units text,
    uncertainty_type text,
    standard_deviation_typical double precision,
    measurement_protocol text
);


ALTER TABLE public.quantum_observables_ref OWNER TO armand;

--
-- Name: TABLE quantum_observables_ref; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.quantum_observables_ref IS 'Reference data for quantum mechanical observables';


--
-- Name: quantum_observables_ref_id_seq; Type: SEQUENCE; Schema: public; Owner: armand
--

CREATE SEQUENCE public.quantum_observables_ref_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.quantum_observables_ref_id_seq OWNER TO armand;

--
-- Name: quantum_observables_ref_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: armand
--

ALTER SEQUENCE public.quantum_observables_ref_id_seq OWNED BY public.quantum_observables_ref.id;


--
-- Name: quantum_parameters; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.quantum_parameters (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    parameter_name text NOT NULL,
    description text,
    unit text,
    calculation_method text NOT NULL,
    typical_range jsonb,
    accuracy_metrics jsonb,
    validation_criteria jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.quantum_parameters OWNER TO armand;

--
-- Name: quantum_properties; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.quantum_properties (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    total_electronic_energy double precision,
    homo_lumo_gap double precision,
    electron_density jsonb,
    orbital_energies double precision[],
    critical_temperature double precision,
    critical_pressure double precision,
    correlation_length double precision,
    order_parameter jsonb,
    coherence_time double precision,
    relaxation_rate double precision,
    dephasing_rate double precision,
    transition_type text,
    transition_order integer,
    universality_class text,
    basis_set text,
    method text,
    calculation_parameters jsonb,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    CONSTRAINT valid_quantum_ranges CHECK (((total_electronic_energy < (0)::double precision) AND (homo_lumo_gap >= (0)::double precision) AND (critical_temperature >= (0)::double precision) AND (critical_pressure >= (0)::double precision) AND (correlation_length > (0)::double precision) AND (coherence_time >= (0)::double precision) AND (relaxation_rate >= (0)::double precision) AND (dephasing_rate >= (0)::double precision) AND (transition_order > 0)))
);


ALTER TABLE public.quantum_properties OWNER TO armand;

--
-- Name: TABLE quantum_properties; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.quantum_properties IS 'Consolidated quantum mechanical properties for compounds';


--
-- Name: quantum_research_findings; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.quantum_research_findings (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    calculation_id uuid,
    finding_type text,
    description text,
    numerical_value double precision,
    units text,
    confidence_level double precision,
    methodology text,
    validation_method text,
    publication_reference text,
    discovery_date date,
    significance_level text,
    metadata jsonb,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    CONSTRAINT check_confidence CHECK (((confidence_level >= (0)::double precision) AND (confidence_level <= (1)::double precision))),
    CONSTRAINT check_significance CHECK ((significance_level = ANY (ARRAY['high'::text, 'medium'::text, 'low'::text])))
);


ALTER TABLE public.quantum_research_findings OWNER TO armand;

--
-- Name: TABLE quantum_research_findings; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.quantum_research_findings IS 'Research findings from quantum calculations';


--
-- Name: quantum_research_projects; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.quantum_research_projects (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    start_date date,
    end_date date,
    status text,
    principal_investigator text,
    research_type text,
    funding_source text,
    budget numeric,
    objectives text[],
    metadata jsonb,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL
);


ALTER TABLE public.quantum_research_projects OWNER TO armand;

--
-- Name: TABLE quantum_research_projects; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.quantum_research_projects IS 'Research projects involving quantum mechanical studies';


--
-- Name: quantum_structure_correlations; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.quantum_structure_correlations (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid,
    property_type text,
    property_value double precision,
    correlation_type text,
    correlation_coefficient double precision,
    statistical_significance double precision,
    sample_size integer,
    methodology text,
    validation_metrics jsonb,
    metadata jsonb,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    CONSTRAINT check_correlation CHECK (((correlation_coefficient >= ('-1'::integer)::double precision) AND (correlation_coefficient <= (1)::double precision))),
    CONSTRAINT check_significance_value CHECK (((statistical_significance >= (0)::double precision) AND (statistical_significance <= (1)::double precision)))
);


ALTER TABLE public.quantum_structure_correlations OWNER TO armand;

--
-- Name: TABLE quantum_structure_correlations; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.quantum_structure_correlations IS 'Structure-property correlations from quantum studies';


--
-- Name: quantum_structure_insights; Type: VIEW; Schema: public; Owner: armand
--

CREATE VIEW public.quantum_structure_insights AS
 SELECT c.compound_id,
    count(DISTINCT c.id) AS calculation_count,
    count(DISTINCT f.id) AS finding_count,
    count(DISTINCT corr.id) AS correlation_count,
    string_agg(DISTINCT c.calculation_type, ', '::text) AS calculation_types,
    string_agg(DISTINCT f.finding_type, ', '::text) AS finding_types
   FROM ((public.quantum_calculations c
     LEFT JOIN public.quantum_research_findings f ON ((c.id = f.calculation_id)))
     LEFT JOIN public.quantum_structure_correlations corr ON ((c.compound_id = corr.compound_id)))
  GROUP BY c.compound_id;


ALTER TABLE public.quantum_structure_insights OWNER TO armand;

--
-- Name: quantum_transport_mechanism_summary; Type: VIEW; Schema: public; Owner: armand
--

CREATE VIEW public.quantum_transport_mechanism_summary AS
 SELECT c.id AS compound_id,
    c.name AS compound_name,
    tm.mechanism_type,
    tm.transport_regime,
    tm.efficiency_score,
    tm.mechanism_confidence,
    qcs.distribution_type,
    qcs.d2 AS correlation_dimension,
    qcs.disorder_strength,
    qcs.is_critical,
    enaqt.transport_efficiency,
    enaqt.coherence_time,
    enaqt.decoherence_rate,
    enaqt.anti_zeno_factor,
    enaqt.is_optimal AS transport_optimal,
        CASE
            WHEN ((tm.mechanism_type = 'ENAQT'::text) AND qcs.is_critical) THEN 'Optimal quantum transport at criticality'::text
            WHEN (tm.mechanism_type = 'ENAQT'::text) THEN 'Environment-assisted transport'::text
            WHEN qcs.is_critical THEN 'Critical but not environment-assisted'::text
            WHEN (tm.mechanism_type = 'Quantum'::text) THEN 'Quantum transport'::text
            ELSE 'Classical transport'::text
        END AS transport_classification,
    ((tm.mechanism_confidence + qcs.criticality_confidence) / (2.0)::double precision) AS overall_confidence
   FROM (((public.compounds c
     LEFT JOIN LATERAL ( SELECT analyze_transport_mechanism.mechanism_type,
            analyze_transport_mechanism.transport_regime,
            analyze_transport_mechanism.efficiency_score,
            analyze_transport_mechanism.mechanism_confidence
           FROM public.analyze_transport_mechanism(c.id) analyze_transport_mechanism(mechanism_type, transport_regime, efficiency_score, mechanism_confidence)) tm ON (true))
     LEFT JOIN public.quantum_criticality_summary qcs ON ((c.id = qcs.compound_id)))
     LEFT JOIN LATERAL ( SELECT analyze_enaqt_properties.transport_efficiency,
            analyze_enaqt_properties.coherence_time,
            analyze_enaqt_properties.decoherence_rate,
            analyze_enaqt_properties.anti_zeno_factor,
            analyze_enaqt_properties.is_optimal
           FROM public.analyze_enaqt_properties(c.id) analyze_enaqt_properties(transport_efficiency, coherence_time, decoherence_rate, anti_zeno_factor, is_optimal)) enaqt ON (true));


ALTER TABLE public.quantum_transport_mechanism_summary OWNER TO armand;

--
-- Name: VIEW quantum_transport_mechanism_summary; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON VIEW public.quantum_transport_mechanism_summary IS 'Comprehensive view of quantum transport properties combining:
1. Transport mechanism classification (ENAQT, Quantum, Classical)
2. Quantum criticality indicators
3. Environment-assisted transport metrics
4. Overall transport efficiency assessment';


--
-- Name: quantum_transport_summary; Type: VIEW; Schema: public; Owner: armand
--

CREATE VIEW public.quantum_transport_summary AS
 SELECT c.id AS compound_id,
    c.name AS compound_name,
    qcs.distribution_type,
    qcs.d2 AS correlation_dimension,
    qcs.is_critical,
    enaqt.transport_efficiency,
    enaqt.coherence_time,
    enaqt.decoherence_rate,
    enaqt.anti_zeno_factor,
    enaqt.is_optimal AS transport_optimal,
        CASE
            WHEN (qcs.is_critical AND enaqt.is_optimal) THEN 'Optimal quantum transport'::text
            WHEN qcs.is_critical THEN 'Critical but suboptimal transport'::text
            WHEN enaqt.is_optimal THEN 'Optimal transport but not critical'::text
            ELSE 'Neither critical nor optimal'::text
        END AS transport_classification,
    qcs.criticality_confidence,
    ((qcs.criticality_confidence +
        CASE
            WHEN enaqt.is_optimal THEN (1.0)::double precision
            ELSE enaqt.transport_efficiency
        END) / (2.0)::double precision) AS overall_confidence
   FROM ((public.compounds c
     LEFT JOIN public.quantum_criticality_summary qcs ON ((c.id = qcs.compound_id)))
     LEFT JOIN LATERAL ( SELECT analyze_enaqt_properties.transport_efficiency,
            analyze_enaqt_properties.coherence_time,
            analyze_enaqt_properties.decoherence_rate,
            analyze_enaqt_properties.anti_zeno_factor,
            analyze_enaqt_properties.is_optimal
           FROM public.analyze_enaqt_properties(c.id) analyze_enaqt_properties(transport_efficiency, coherence_time, decoherence_rate, anti_zeno_factor, is_optimal)) enaqt ON (true));


ALTER TABLE public.quantum_transport_summary OWNER TO armand;

--
-- Name: VIEW quantum_transport_summary; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON VIEW public.quantum_transport_summary IS 'Consolidated view of quantum transport properties combining:
1. Quantum criticality indicators (level spacing, multifractality)
2. Environment-assisted Quantum Transport (ENAQT) metrics
3. Overall transport efficiency assessment';


--
-- Name: rate_limits; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.rate_limits (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    api_key_id uuid,
    endpoint_id uuid,
    requests_count integer DEFAULT 0 NOT NULL,
    window_start timestamp with time zone NOT NULL,
    window_end timestamp with time zone NOT NULL,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.rate_limits OWNER TO armand;

--
-- Name: receptor_families; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.receptor_families (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    category_id uuid NOT NULL,
    name text NOT NULL,
    abbreviation text,
    description text,
    protein_type text NOT NULL,
    signaling_type text NOT NULL,
    primary_endogenous_ligands text[],
    primary_effects text[],
    therapeutic_areas text[],
    expression_pattern jsonb,
    signaling_pathways jsonb,
    pharmacological_properties jsonb,
    clinical_significance text,
    research_status text,
    notes text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.receptor_families OWNER TO armand;

--
-- Name: receptor_family_categories; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.receptor_family_categories (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    receptor_type text NOT NULL,
    signaling_mechanism text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.receptor_family_categories OWNER TO armand;

--
-- Name: receptor_subtypes; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.receptor_subtypes (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    family_id uuid NOT NULL,
    subtype_name text NOT NULL,
    description text,
    protein_sequence text,
    species text,
    expression_pattern jsonb,
    signaling_pathways text[],
    pharmacological_profile jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.receptor_subtypes OWNER TO armand;

--
-- Name: research_findings; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.research_findings (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    finding_type text NOT NULL,
    description text NOT NULL,
    methodology text[],
    experimental_data jsonb,
    statistical_analysis jsonb,
    conclusions text[],
    limitations text[],
    future_directions text[],
    reference_dois text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.research_findings OWNER TO armand;

--
-- Name: TABLE research_findings; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.research_findings IS 'Research findings and conclusions';


--
-- Name: safety_documents; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.safety_documents (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    sds_url text,
    msds_url text,
    sds_last_updated date,
    sds_provider text,
    sds_version text,
    safety_data_sheet jsonb,
    handling_precautions text[],
    storage_precautions text[],
    disposal_instructions text[],
    first_aid_measures jsonb,
    firefighting_measures jsonb,
    accidental_release_measures jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.safety_documents OWNER TO armand;

--
-- Name: TABLE safety_documents; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.safety_documents IS 'Safety documentation including SDS and handling instructions';


--
-- Name: safety_thresholds; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.safety_thresholds (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    parameter_name text NOT NULL,
    threshold_value double precision NOT NULL,
    unit text,
    severity_level text NOT NULL,
    description text,
    intervention_required boolean DEFAULT false,
    validation_method text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.safety_thresholds OWNER TO armand;

--
-- Name: TABLE safety_thresholds; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.safety_thresholds IS 'Safety threshold values for various clinical and laboratory parameters';


--
-- Name: COLUMN safety_thresholds.parameter_name; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.safety_thresholds.parameter_name IS 'Name of the safety parameter being measured';


--
-- Name: COLUMN safety_thresholds.threshold_value; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.safety_thresholds.threshold_value IS 'Numerical threshold value';


--
-- Name: COLUMN safety_thresholds.unit; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.safety_thresholds.unit IS 'Unit of measurement';


--
-- Name: COLUMN safety_thresholds.severity_level; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.safety_thresholds.severity_level IS 'Severity level when threshold is exceeded';


--
-- Name: COLUMN safety_thresholds.description; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.safety_thresholds.description IS 'Description of the threshold and its significance';


--
-- Name: sar_analysis; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.sar_analysis (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    analysis_type text NOT NULL,
    structural_features text[],
    activity_correlations jsonb,
    pharmacophore_model jsonb,
    binding_patterns jsonb,
    selectivity_patterns jsonb,
    structure_modifications text[],
    predicted_effects jsonb,
    confidence_metrics jsonb,
    validation_results jsonb,
    reference_compounds text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.sar_analysis OWNER TO armand;

--
-- Name: TABLE sar_analysis; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.sar_analysis IS 'Structure-activity relationship analysis results';


--
-- Name: sar_patterns; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.sar_patterns (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    pattern_type text NOT NULL,
    structural_elements text[],
    activity_impact jsonb,
    confidence_score double precision,
    supporting_compounds text[],
    detection_method text,
    validation_status text,
    notes text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.sar_patterns OWNER TO armand;

--
-- Name: TABLE sar_patterns; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.sar_patterns IS 'Identified structure-activity relationship patterns';


--
-- Name: scaffold_analysis; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.scaffold_analysis (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    scaffold_smiles text NOT NULL,
    compound_count integer,
    average_activity double precision,
    activity_range jsonb,
    diversity_score double precision,
    important_substitutions text[],
    analysis_method text,
    notes text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.scaffold_analysis OWNER TO armand;

--
-- Name: TABLE scaffold_analysis; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.scaffold_analysis IS 'Analysis of molecular scaffolds and their properties';


--
-- Name: scaling_analysis; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.scaling_analysis (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    phase_transition_id uuid,
    scaling_function_type text,
    scaling_variables jsonb,
    scaling_dimensions double precision[],
    rg_flow_equations jsonb,
    fixed_points jsonb,
    relevant_operators jsonb,
    universality_class text,
    central_charge double precision,
    operator_spectrum jsonb,
    size_scaling_exponents double precision[],
    correction_exponents double precision[],
    crossover_scales jsonb,
    crossover_functions jsonb,
    analysis_method text,
    confidence_metrics jsonb,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP NOT NULL,
    metadata jsonb
);


ALTER TABLE public.scaling_analysis OWNER TO armand;

--
-- Name: TABLE scaling_analysis; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.scaling_analysis IS 'Scaling analysis and renormalization group results';


--
-- Name: COLUMN scaling_analysis.rg_flow_equations; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON COLUMN public.scaling_analysis.rg_flow_equations IS 'Renormalization group flow equations stored in JSONB format';


--
-- Name: schema_versions; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.schema_versions (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    version text NOT NULL,
    description text,
    applied_at timestamp with time zone DEFAULT now() NOT NULL,
    applied_by text,
    script_name text,
    checksum text,
    execution_time interval,
    status text NOT NULL,
    error_message text
);


ALTER TABLE public.schema_versions OWNER TO armand;

--
-- Name: scientific_papers; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.scientific_papers (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    doi text,
    title text NOT NULL,
    authors text[] NOT NULL,
    journal text,
    publication_date date,
    abstract text,
    full_text text,
    methodology text[],
    study_type text,
    sample_size integer,
    study_duration interval,
    quality_metrics jsonb,
    evidence_level text,
    key_findings text[],
    limitations text[],
    compounds_studied uuid[],
    validation_status text,
    peer_review_status text,
    citation_count integer,
    impact_factor double precision,
    external_links jsonb,
    supplementary_data jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.scientific_papers OWNER TO armand;

--
-- Name: TABLE scientific_papers; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.scientific_papers IS 'Scientific literature and research papers';


--
-- Name: service_status; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.service_status (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    service_id uuid NOT NULL,
    status text NOT NULL,
    response_time interval,
    error_message text,
    check_timestamp timestamp with time zone NOT NULL,
    metrics jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.service_status OWNER TO armand;

--
-- Name: social_alert_rules; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_alert_rules (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    alert_type text NOT NULL,
    severity text NOT NULL,
    platforms text[] NOT NULL,
    platform_subdivisions text[],
    compounds uuid[],
    trigger_conditions jsonb NOT NULL,
    required_metrics text[] NOT NULL,
    threshold_values jsonb NOT NULL,
    cooldown_period interval,
    is_active boolean DEFAULT true NOT NULL,
    last_triggered_at timestamp with time zone,
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_alert_rules OWNER TO armand;

--
-- Name: social_comments; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_comments (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    post_id uuid NOT NULL,
    platform text NOT NULL,
    platform_subdivision text,
    external_id text NOT NULL,
    parent_id text,
    content text NOT NULL,
    author_id text NOT NULL,
    author_username text NOT NULL,
    created_at timestamp with time zone NOT NULL,
    engagement_metrics jsonb DEFAULT '{}'::jsonb NOT NULL,
    classification_data jsonb DEFAULT '{}'::jsonb NOT NULL,
    sentiment_data jsonb DEFAULT '{}'::jsonb NOT NULL,
    is_scientific boolean DEFAULT false NOT NULL,
    has_citations boolean DEFAULT false NOT NULL,
    reported_effects text[] DEFAULT '{}'::text[] NOT NULL,
    reported_side_effects text[] DEFAULT '{}'::text[] NOT NULL,
    platform_specific_data jsonb DEFAULT '{}'::jsonb NOT NULL,
    created_at_internal timestamp with time zone DEFAULT now() NOT NULL,
    CONSTRAINT check_comment_classification CHECK ((((((classification_data ->> 'spam_score'::text))::double precision >= (0)::double precision) AND (((classification_data ->> 'spam_score'::text))::double precision <= (1)::double precision)) AND ((((classification_data ->> 'toxicity_score'::text))::double precision >= (0)::double precision) AND (((classification_data ->> 'toxicity_score'::text))::double precision <= (1)::double precision)))),
    CONSTRAINT check_comment_engagement CHECK (
CASE platform
    WHEN 'reddit'::text THEN (((engagement_metrics ->> 'score'::text))::integer >= 0)
    WHEN 'twitter'::text THEN (((engagement_metrics ->> 'likes'::text))::integer >= 0)
    ELSE true
END),
    CONSTRAINT check_comment_sentiment CHECK ((((((sentiment_data ->> 'score'::text))::double precision >= ('-1'::integer)::double precision) AND (((sentiment_data ->> 'score'::text))::double precision <= (1)::double precision)) AND ((((sentiment_data ->> 'positive_ratio'::text))::double precision >= (0)::double precision) AND (((sentiment_data ->> 'positive_ratio'::text))::double precision <= (1)::double precision)) AND ((((sentiment_data ->> 'negative_ratio'::text))::double precision >= (0)::double precision) AND (((sentiment_data ->> 'negative_ratio'::text))::double precision <= (1)::double precision)) AND ((((sentiment_data ->> 'neutral_ratio'::text))::double precision >= (0)::double precision) AND (((sentiment_data ->> 'neutral_ratio'::text))::double precision <= (1)::double precision))))
);


ALTER TABLE public.social_comments OWNER TO armand;

--
-- Name: social_compound_combinations; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_compound_combinations (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id_1 uuid NOT NULL,
    compound_id_2 uuid NOT NULL,
    interaction_type text NOT NULL,
    risk_level text NOT NULL,
    description text,
    mechanism text,
    evidence_level text,
    sources text[],
    platforms text[] NOT NULL,
    reported_count integer DEFAULT 0 NOT NULL,
    reports uuid[],
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_compound_combinations OWNER TO armand;

--
-- Name: social_compound_mentions; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_compound_mentions (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    platform text NOT NULL,
    platform_subdivision text,
    analysis_date date NOT NULL,
    total_mentions integer DEFAULT 0 NOT NULL,
    unique_authors integer DEFAULT 0 NOT NULL,
    total_engagement integer DEFAULT 0 NOT NULL,
    scientific_mentions integer DEFAULT 0 NOT NULL,
    experience_reports integer DEFAULT 0 NOT NULL,
    harm_reduction_mentions integer DEFAULT 0 NOT NULL,
    sentiment_distribution jsonb,
    topic_distribution jsonb,
    user_demographics jsonb,
    geographic_distribution jsonb,
    temporal_patterns jsonb,
    common_contexts text[],
    related_compounds text[],
    platform_specific_metrics jsonb,
    cross_platform_engagement_flow jsonb,
    user_influence_metrics jsonb,
    content_propagation_patterns jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_compound_mentions OWNER TO armand;

--
-- Name: social_compound_trends; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_compound_trends (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    trend_start_date date NOT NULL,
    trend_end_date date NOT NULL,
    platforms text[] NOT NULL,
    total_mentions integer DEFAULT 0 NOT NULL,
    unique_authors integer DEFAULT 0 NOT NULL,
    total_engagement integer DEFAULT 0 NOT NULL,
    trend_velocity numeric,
    trend_acceleration numeric,
    platform_distribution jsonb,
    topic_evolution jsonb,
    sentiment_evolution jsonb,
    influencer_impact jsonb,
    news_correlation jsonb,
    research_correlation jsonb,
    correlation_strength numeric,
    context_similarity numeric,
    user_overlap_patterns jsonb,
    temporal_similarity numeric,
    platform_specific_metrics jsonb,
    trend_confidence_score numeric,
    trend_validation_metrics jsonb,
    cross_platform_correlations jsonb,
    trend_seasonality jsonb,
    trend_anomaly_scores jsonb,
    trend_prediction_metrics jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_compound_trends OWNER TO armand;

--
-- Name: social_content_quality; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_content_quality (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    platform text NOT NULL,
    platform_subdivision text,
    analysis_date date NOT NULL,
    scientific_accuracy_score numeric,
    harm_reduction_quality_score numeric,
    experience_report_quality_score numeric,
    information_completeness_score numeric,
    citation_quality_score numeric,
    misinformation_prevalence_score numeric,
    content_depth_distribution jsonb,
    quality_trends jsonb,
    improvement_suggestions text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_content_quality OWNER TO armand;

--
-- Name: social_dosage_data; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_dosage_data (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    platform text NOT NULL,
    platform_subdivision text,
    route_of_administration text NOT NULL,
    threshold_dose text,
    light_dose text,
    common_dose text,
    strong_dose text,
    heavy_dose text,
    warning_dose text,
    duration_total interval,
    duration_onset interval,
    duration_peak interval,
    duration_offset interval,
    bioavailability numeric,
    dosage_notes text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_dosage_data OWNER TO armand;

--
-- Name: social_effect_reports; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_effect_reports (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    effect_id uuid NOT NULL,
    report_text text NOT NULL,
    intensity integer,
    duration interval,
    onset_time interval,
    conditions jsonb,
    platform text NOT NULL,
    platform_subdivision text,
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_effect_reports OWNER TO armand;

--
-- Name: social_effects; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_effects (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    effect_name text NOT NULL,
    description text,
    type text,
    url text,
    analysis_data jsonb,
    related_effects text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_effects OWNER TO armand;

--
-- Name: social_experience_reports; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_experience_reports (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    post_id uuid,
    compound_id uuid,
    stack_id uuid,
    protocol_id uuid,
    title text NOT NULL,
    content text NOT NULL,
    author text,
    platform text NOT NULL,
    platform_subdivision text,
    experience_date date,
    report_date timestamp with time zone NOT NULL,
    substance_data jsonb[],
    duration interval,
    setting_data jsonb,
    intention text,
    effects_timeline jsonb[],
    timeline_data jsonb,
    reported_effects text[],
    side_effects text[],
    interactions text[],
    test_kit_info jsonb,
    detection_time_data jsonb,
    after_effects_data jsonb,
    body_weight numeric,
    weight_unit text,
    gender text,
    age integer,
    experience_level text,
    harm_reduction_notes text[],
    classification_data jsonb,
    sentiment_data jsonb,
    report_version integer,
    experience_category text[],
    total_views integer,
    report_quality_score double precision,
    medical_conditions text[],
    medications text[],
    baseline_metrics jsonb,
    outcome_metrics jsonb,
    testing_methods text[],
    overall_rating integer,
    platform_specific_data jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_experience_reports OWNER TO armand;

--
-- Name: social_harm_reduction; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_harm_reduction (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    post_id uuid,
    compound_id uuid,
    title text NOT NULL,
    content text NOT NULL,
    author text,
    platform text,
    platform_subdivision text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    category text NOT NULL,
    importance_level text,
    safety_notes text[],
    warnings text[],
    contraindications text[],
    emergency_procedures text[],
    sources text[],
    verification_status text DEFAULT 'unverified'::text NOT NULL,
    verification_notes text,
    last_reviewed_at timestamp with time zone,
    created_at_internal timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_harm_reduction OWNER TO armand;

--
-- Name: social_influence_networks; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_influence_networks (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    platform text NOT NULL,
    platform_subdivision text,
    analysis_date date NOT NULL,
    top_influencers jsonb[],
    influence_connections jsonb,
    community_clusters jsonb,
    information_flow_patterns jsonb,
    key_opinion_leaders text[],
    emerging_voices text[],
    platform_specific_metrics jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_influence_networks OWNER TO armand;

--
-- Name: social_platform_stats; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_platform_stats (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    platform text NOT NULL,
    platform_subdivision text,
    compound_id uuid,
    total_posts integer DEFAULT 0 NOT NULL,
    total_comments integer DEFAULT 0 NOT NULL,
    unique_authors integer DEFAULT 0 NOT NULL,
    scientific_post_ratio numeric,
    experience_report_ratio numeric,
    harm_reduction_ratio numeric,
    top_compounds text[],
    topic_distribution jsonb,
    sentiment_distribution jsonb,
    quality_metrics jsonb,
    engagement_metrics jsonb,
    last_updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_platform_stats OWNER TO armand;

--
-- Name: social_posts; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_posts (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    platform text NOT NULL,
    platform_subdivision text,
    external_id text NOT NULL,
    url text NOT NULL,
    title text,
    content text NOT NULL,
    author_id text NOT NULL,
    author_username text NOT NULL,
    post_type text NOT NULL,
    created_at timestamp with time zone NOT NULL,
    engagement_metrics jsonb DEFAULT '{}'::jsonb NOT NULL,
    classification_data jsonb DEFAULT '{}'::jsonb NOT NULL,
    sentiment_data jsonb DEFAULT '{}'::jsonb NOT NULL,
    metadata jsonb DEFAULT '{}'::jsonb NOT NULL,
    is_scientific boolean DEFAULT false NOT NULL,
    is_experience_report boolean DEFAULT false NOT NULL,
    is_harm_reduction boolean DEFAULT false NOT NULL,
    platform_specific_data jsonb DEFAULT '{}'::jsonb NOT NULL,
    tags text[] DEFAULT '{}'::text[] NOT NULL,
    research_citations text[] DEFAULT '{}'::text[] NOT NULL,
    view_count integer DEFAULT 0 NOT NULL,
    created_at_internal timestamp with time zone DEFAULT now() NOT NULL,
    CONSTRAINT check_classification_scores CHECK ((((((classification_data ->> 'spam_score'::text))::double precision >= (0)::double precision) AND (((classification_data ->> 'spam_score'::text))::double precision <= (1)::double precision)) AND ((((classification_data ->> 'toxicity_score'::text))::double precision >= (0)::double precision) AND (((classification_data ->> 'toxicity_score'::text))::double precision <= (1)::double precision)))),
    CONSTRAINT check_engagement_metrics CHECK (
CASE platform
    WHEN 'reddit'::text THEN (((engagement_metrics ->> 'score'::text))::integer >= 0)
    WHEN 'twitter'::text THEN ((((engagement_metrics ->> 'retweet_count'::text))::integer >= 0) AND (((engagement_metrics ->> 'favorite_count'::text))::integer >= 0))
    ELSE true
END),
    CONSTRAINT check_sentiment_scores CHECK ((((((sentiment_data ->> 'score'::text))::double precision >= ('-1'::integer)::double precision) AND (((sentiment_data ->> 'score'::text))::double precision <= (1)::double precision)) AND ((((sentiment_data ->> 'positive_ratio'::text))::double precision >= (0)::double precision) AND (((sentiment_data ->> 'positive_ratio'::text))::double precision <= (1)::double precision)) AND ((((sentiment_data ->> 'negative_ratio'::text))::double precision >= (0)::double precision) AND (((sentiment_data ->> 'negative_ratio'::text))::double precision <= (1)::double precision)) AND ((((sentiment_data ->> 'neutral_ratio'::text))::double precision >= (0)::double precision) AND (((sentiment_data ->> 'neutral_ratio'::text))::double precision <= (1)::double precision))))
);


ALTER TABLE public.social_posts OWNER TO armand;

--
-- Name: social_protocols; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_protocols (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    protocol_name text NOT NULL,
    creator text,
    description text NOT NULL,
    target_outcome text,
    compounds uuid[] NOT NULL,
    protocol_steps jsonb[],
    duration text,
    frequency text,
    monitoring_parameters text[],
    success_metrics text[],
    warnings text[],
    contraindications text[],
    reference_citations text[],
    review_score numeric,
    platform text NOT NULL,
    platform_subdivision text,
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_protocols OWNER TO armand;

--
-- Name: social_research_reviews; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_research_reviews (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    title text NOT NULL,
    content text NOT NULL,
    author text,
    publication_date timestamp with time zone NOT NULL,
    study_type text[],
    methodology text,
    findings text[],
    limitations text[],
    research_quality_score numeric,
    citations text[],
    peer_review_notes text[],
    platform text NOT NULL,
    platform_subdivision text,
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_research_reviews OWNER TO armand;

--
-- Name: social_safety_incidents; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_safety_incidents (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    incident_type text NOT NULL,
    severity text NOT NULL,
    description text NOT NULL,
    reported_effects text[],
    reported_causes text[],
    platforms text[] NOT NULL,
    platform_subdivisions text[],
    content_urls text[],
    verification_status text DEFAULT 'unverified'::text NOT NULL,
    verification_notes text,
    response_actions text[],
    resolution_status text DEFAULT 'open'::text NOT NULL,
    resolution_notes text,
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_safety_incidents OWNER TO armand;

--
-- Name: social_scientific_content; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_scientific_content (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    post_id uuid,
    compound_id uuid,
    title text NOT NULL,
    content text NOT NULL,
    author text,
    platform text,
    platform_subdivision text,
    publication_date timestamp with time zone NOT NULL,
    content_type text NOT NULL,
    research_topics text[],
    methodology text,
    findings text[],
    limitations text[],
    citations text[],
    peer_review_notes text[],
    quality_score numeric,
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_scientific_content OWNER TO armand;

--
-- Name: social_stacks; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_stacks (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    stack_name text NOT NULL,
    creator text,
    description text,
    purpose text,
    compounds uuid[] NOT NULL,
    dosages jsonb,
    timing_schedule jsonb,
    duration text,
    reported_effects text[],
    side_effects text[],
    interactions text[],
    warnings text[],
    rating numeric,
    review_count integer,
    platform text NOT NULL,
    platform_subdivision text,
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_stacks OWNER TO armand;

--
-- Name: social_substance_data; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.social_substance_data (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    platform text NOT NULL,
    platform_subdivision text,
    external_id text NOT NULL,
    common_names text[],
    chemical_class text[],
    psychoactive_class text[],
    summary text,
    tolerance_data jsonb,
    roa_data jsonb,
    effect_data jsonb,
    onset text,
    duration text,
    after_effects text,
    detection_time text,
    test_kits jsonb,
    experiences text[],
    effects text[],
    aliases text[],
    avoid text[],
    warning_message text,
    dangerous_interactions text[],
    uncertain_interactions text[],
    unsafe_interactions text[],
    risk_potential jsonb,
    toxicity_data text,
    addiction_potential text,
    cross_tolerance text[],
    chemistry_data jsonb,
    dosage_data jsonb,
    duration_data jsonb,
    health_effects jsonb,
    risk_factors text[],
    contraindications text[],
    interactions jsonb,
    legal_status jsonb,
    research_status text,
    history text,
    traditional_use text,
    harm_reduction_notes text[],
    platform_specific_data jsonb,
    last_updated_at timestamp with time zone,
    created_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.social_substance_data OWNER TO armand;

--
-- Name: structure_similarity; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.structure_similarity (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id_1 uuid NOT NULL,
    compound_id_2 uuid NOT NULL,
    similarity_metric text NOT NULL,
    similarity_score double precision,
    comparison_method text,
    fingerprint_type text,
    calculation_parameters jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.structure_similarity OWNER TO armand;

--
-- Name: TABLE structure_similarity; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.structure_similarity IS 'Pairwise structural similarity between compounds';


--
-- Name: subjective_effect_categories; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.subjective_effect_categories (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    domain text NOT NULL,
    level integer NOT NULL,
    parent_category_id uuid,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL,
    CONSTRAINT valid_level CHECK ((level > 0))
);


ALTER TABLE public.subjective_effect_categories OWNER TO armand;

--
-- Name: subjective_effects; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.subjective_effects (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    category_id uuid NOT NULL,
    name text NOT NULL,
    description text,
    onset_characteristics jsonb,
    duration_characteristics jsonb,
    intensity_characteristics jsonb,
    common_variations text[],
    contributing_factors text[],
    risk_factors text[],
    management_strategies text[],
    research_status text,
    notes text,
    text_search_vector tsvector GENERATED ALWAYS AS ((setweight(to_tsvector('english'::regconfig, COALESCE(name, ''::text)), 'A'::"char") || setweight(to_tsvector('english'::regconfig, COALESCE(description, ''::text)), 'B'::"char"))) STORED,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.subjective_effects OWNER TO armand;

--
-- Name: sync_status; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.sync_status (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    device_id text NOT NULL,
    data_type text NOT NULL,
    last_sync_at timestamp with time zone NOT NULL,
    sync_version integer NOT NULL,
    status text NOT NULL,
    conflict_resolution jsonb,
    retry_count integer DEFAULT 0,
    next_retry_at timestamp with time zone,
    error_details text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.sync_status OWNER TO armand;

--
-- Name: systems_biology_analysis; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.systems_biology_analysis (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    compound_id uuid NOT NULL,
    analysis_level text NOT NULL,
    network_effects jsonb,
    cellular_responses jsonb,
    metabolic_impact jsonb,
    signaling_cascades jsonb,
    regulatory_networks jsonb,
    adaptation_mechanisms jsonb,
    system_robustness jsonb,
    emergent_properties jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.systems_biology_analysis OWNER TO armand;

--
-- Name: TABLE systems_biology_analysis; Type: COMMENT; Schema: public; Owner: armand
--

COMMENT ON TABLE public.systems_biology_analysis IS 'Systems-level biological analysis';


--
-- Name: therapeutic_class_categories; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.therapeutic_class_categories (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    level integer NOT NULL,
    parent_category_id uuid,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL,
    CONSTRAINT valid_level CHECK ((level > 0))
);


ALTER TABLE public.therapeutic_class_categories OWNER TO armand;

--
-- Name: therapeutic_classes; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.therapeutic_classes (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    category_id uuid NOT NULL,
    name text NOT NULL,
    abbreviation text,
    description text,
    mechanism_of_action text,
    primary_targets text[],
    therapeutic_uses text[],
    contraindications text[],
    typical_dosing jsonb,
    side_effects jsonb,
    drug_interactions jsonb,
    regulatory_status text,
    clinical_guidelines text[],
    research_status text,
    notes text,
    text_search_vector tsvector GENERATED ALWAYS AS (((setweight(to_tsvector('english'::regconfig, COALESCE(name, ''::text)), 'A'::"char") || setweight(to_tsvector('english'::regconfig, COALESCE(description, ''::text)), 'B'::"char")) || setweight(to_tsvector('english'::regconfig, COALESCE(mechanism_of_action, ''::text)), 'C'::"char"))) STORED,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.therapeutic_classes OWNER TO armand;

--
-- Name: toxicity_endpoint_categories; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.toxicity_endpoint_categories (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    measurement_type text NOT NULL,
    units text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.toxicity_endpoint_categories OWNER TO armand;

--
-- Name: toxicity_endpoints; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.toxicity_endpoints (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    category_id uuid NOT NULL,
    name text NOT NULL,
    description text,
    standard_unit text NOT NULL,
    conversion_factors jsonb,
    detection_methods text[],
    validation_criteria jsonb,
    reference_ranges jsonb,
    severity_thresholds jsonb,
    regulatory_limits jsonb,
    notes text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.toxicity_endpoints OWNER TO armand;

--
-- Name: toxicity_mechanism_categories; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.toxicity_mechanism_categories (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    level integer NOT NULL,
    mechanism_type text NOT NULL,
    parent_category_id uuid,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL,
    CONSTRAINT valid_level CHECK ((level > 0))
);


ALTER TABLE public.toxicity_mechanism_categories OWNER TO armand;

--
-- Name: toxicity_mechanisms; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.toxicity_mechanisms (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    category_id uuid NOT NULL,
    name text NOT NULL,
    description text,
    molecular_targets text[],
    cellular_effects text[],
    tissue_effects text[],
    systemic_effects text[],
    biomarkers text[],
    detection_methods text[],
    time_course jsonb,
    dose_response_characteristics jsonb,
    reversibility text,
    risk_factors text[],
    preventive_measures text[],
    treatment_approaches text[],
    research_status text,
    evidence_level text,
    notes text,
    text_search_vector tsvector GENERATED ALWAYS AS ((setweight(to_tsvector('english'::regconfig, COALESCE(name, ''::text)), 'A'::"char") || setweight(to_tsvector('english'::regconfig, COALESCE(description, ''::text)), 'B'::"char"))) STORED,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.toxicity_mechanisms OWNER TO armand;

--
-- Name: training_datasets; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.training_datasets (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    dataset_name text NOT NULL,
    description text,
    source text NOT NULL,
    version text NOT NULL,
    compound_count integer,
    feature_count integer,
    feature_names text[],
    feature_types text[],
    target_names text[],
    target_types text[],
    preprocessing_steps jsonb,
    split_strategy text,
    validation_method text,
    data_statistics jsonb,
    quality_metrics jsonb,
    metadata jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.training_datasets OWNER TO armand;

--
-- Name: training_history; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.training_history (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    model_id uuid NOT NULL,
    version_id uuid NOT NULL,
    training_run_id text NOT NULL,
    start_time timestamp with time zone NOT NULL,
    end_time timestamp with time zone,
    parameters jsonb,
    metrics jsonb,
    loss_history jsonb,
    validation_history jsonb,
    hardware_metrics jsonb,
    status text NOT NULL,
    error_logs text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.training_history OWNER TO armand;

--
-- Name: trend_analysis; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.trend_analysis (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    analysis_type text NOT NULL,
    target_type text NOT NULL,
    target_id uuid NOT NULL,
    time_period text NOT NULL,
    metrics jsonb NOT NULL,
    insights jsonb,
    recommendations jsonb,
    priority text,
    status text,
    assigned_to uuid,
    resolution_notes text,
    next_review_date timestamp with time zone,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.trend_analysis OWNER TO armand;

--
-- Name: ui_settings; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.ui_settings (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    user_id text NOT NULL,
    theme text,
    layout_preferences jsonb,
    display_options jsonb,
    notification_settings jsonb,
    accessibility_settings jsonb,
    custom_views jsonb[],
    dashboard_config jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.ui_settings OWNER TO armand;

--
-- Name: usage_statistics; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.usage_statistics (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    feature text NOT NULL,
    action text NOT NULL,
    user_agent text,
    device_type text,
    session_id uuid,
    user_id uuid,
    duration_ms integer,
    success boolean,
    error_details text,
    performance_metrics jsonb,
    user_feedback jsonb,
    metadata jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.usage_statistics OWNER TO armand;

--
-- Name: user_preferences; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.user_preferences (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    user_id uuid NOT NULL,
    device_id text NOT NULL,
    theme text,
    layout text,
    notifications jsonb,
    accessibility jsonb,
    data_preferences jsonb,
    sync_settings jsonb,
    privacy_settings jsonb,
    feature_flags jsonb,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.user_preferences OWNER TO armand;

--
-- Name: web_components; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.web_components (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    component_type text NOT NULL,
    configuration jsonb NOT NULL,
    dependencies text[],
    styling jsonb,
    client_scripts jsonb,
    server_scripts jsonb,
    version text NOT NULL,
    is_active boolean DEFAULT true,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.web_components OWNER TO armand;

--
-- Name: web_data_sources; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.web_data_sources (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    category_id uuid NOT NULL,
    name text NOT NULL,
    base_url text NOT NULL,
    description text,
    api_endpoint text,
    access_method text NOT NULL,
    authentication_type text,
    rate_limits jsonb,
    data_format text,
    update_frequency text,
    last_validated timestamp with time zone,
    validation_status text,
    data_quality_metrics jsonb,
    coverage_areas text[],
    known_limitations text[],
    usage_requirements text,
    citation_format text,
    notes text,
    active boolean DEFAULT true NOT NULL,
    text_search_vector tsvector GENERATED ALWAYS AS ((setweight(to_tsvector('english'::regconfig, COALESCE(name, ''::text)), 'A'::"char") || setweight(to_tsvector('english'::regconfig, COALESCE(description, ''::text)), 'B'::"char"))) STORED,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.web_data_sources OWNER TO armand;

--
-- Name: web_hook_deliveries; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.web_hook_deliveries (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    hook_id uuid NOT NULL,
    event_type text NOT NULL,
    payload jsonb NOT NULL,
    response_status integer,
    response_body text,
    delivery_status text NOT NULL,
    attempt_count integer DEFAULT 0,
    next_retry_at timestamp with time zone,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.web_hook_deliveries OWNER TO armand;

--
-- Name: web_hooks; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.web_hooks (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    url text NOT NULL,
    event_types text[] NOT NULL,
    headers jsonb,
    is_active boolean DEFAULT true,
    secret_key text,
    retry_config jsonb,
    timeout interval,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.web_hooks OWNER TO armand;

--
-- Name: web_source_categories; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.web_source_categories (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    source_type text NOT NULL,
    reliability_rating integer NOT NULL,
    validation_requirements text[],
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL,
    CONSTRAINT web_source_categories_reliability_rating_check CHECK (((reliability_rating >= 1) AND (reliability_rating <= 5)))
);


ALTER TABLE public.web_source_categories OWNER TO armand;

--
-- Name: web_templates; Type: TABLE; Schema: public; Owner: armand
--

CREATE TABLE public.web_templates (
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    name text NOT NULL,
    description text,
    template_type text NOT NULL,
    content text NOT NULL,
    parameters jsonb,
    styling jsonb,
    scripts jsonb,
    version text NOT NULL,
    is_active boolean DEFAULT true,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL
);


ALTER TABLE public.web_templates OWNER TO armand;

--
-- Name: phase_transition_types id; Type: DEFAULT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.phase_transition_types ALTER COLUMN id SET DEFAULT nextval('public.phase_transition_types_id_seq'::regclass);


--
-- Name: quantum_basis_sets id; Type: DEFAULT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_basis_sets ALTER COLUMN id SET DEFAULT nextval('public.quantum_basis_sets_id_seq'::regclass);


--
-- Name: quantum_functionals id; Type: DEFAULT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_functionals ALTER COLUMN id SET DEFAULT nextval('public.quantum_functionals_id_seq'::regclass);


--
-- Name: quantum_observables_ref id; Type: DEFAULT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_observables_ref ALTER COLUMN id SET DEFAULT nextval('public.quantum_observables_ref_id_seq'::regclass);


--
-- Data for Name: activity_cliffs; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.activity_cliffs (id, compound_pair, activity_difference, structural_similarity, cliff_magnitude, activity_type, detection_method, significance_score, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: activity_correlations; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.activity_correlations (id, compound_id, feature_id, activity_type, correlation_coefficient, statistical_significance, analysis_method, sample_size, confidence_interval, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: alert_triggers; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.alert_triggers (id, trigger_type, target_type, target_id, conditions, actions, severity, priority, notification_channels, escalation_rules, cooldown_period, is_active, last_triggered_at, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: analysis_parameters; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.analysis_parameters (id, parameter_name, description, data_type, validation_rules, default_value, allowed_range, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: api_endpoints; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.api_endpoints (id, path, method, description, parameters, response_schema, auth_required, rate_limit, cache_ttl, version, is_deprecated, deprecated_reason, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: api_keys; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.api_keys (id, key_hash, user_id, name, permissions, rate_limit, expires_at, last_used_at, is_active, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: audit_log; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.audit_log (id, table_name, record_id, action, old_data, new_data, changed_by, changed_at) FROM stdin;
c6c5a6e5-4ace-476e-b3cd-fbd069f7e766	schema_versions	6cd8cb4b-94b3-4ba1-b242-a2dfa5c69304	INSERT	\N	{"id": "6cd8cb4b-94b3-4ba1-b242-a2dfa5c69304", "status": "SUCCESS", "version": "0.1.0", "checksum": null, "applied_at": "2025-01-08T03:46:07.874039+00:00", "applied_by": "armand", "description": "Initial core schema setup", "script_name": "00_core.sql", "error_message": null, "execution_time": "00:00:00.000088"}	armand	2025-01-08 03:46:07.874039+00
84830754-e1c5-49bd-91fb-439d39170856	binding_assay_types	2a55ba47-c256-4db7-858e-3a0371ea6703	INSERT	\N	{"id": "2a55ba47-c256-4db7-858e-3a0371ea6703", "name": "Radioligand Binding", "created_at": "2025-01-08T03:46:07.950642+00:00", "updated_at": "2025-01-08T03:46:07.950642+00:00", "description": "Direct measurement using radioactive ligands", "method_type": "Competition", "typical_unit": "Ki (nM)", "detection_type": "Scintillation"}	armand	2025-01-08 03:46:07.950642+00
26b7c461-ca3e-4759-9d3d-b47ce1d7b6c2	binding_assay_types	fad9a424-29a5-4a18-96bd-ffa430676e08	INSERT	\N	{"id": "fad9a424-29a5-4a18-96bd-ffa430676e08", "name": "Fluorescence Binding", "created_at": "2025-01-08T03:46:07.950642+00:00", "updated_at": "2025-01-08T03:46:07.950642+00:00", "description": "Fluorescence-based binding assays", "method_type": "Direct", "typical_unit": "Kd (nM)", "detection_type": "Fluorescence"}	armand	2025-01-08 03:46:07.950642+00
36ebe976-251e-4cf1-88bf-25d3b953b08e	binding_assay_types	27c3d2f1-9821-4f7c-a311-44d1a8d805d8	INSERT	\N	{"id": "27c3d2f1-9821-4f7c-a311-44d1a8d805d8", "name": "SPR", "created_at": "2025-01-08T03:46:07.950642+00:00", "updated_at": "2025-01-08T03:46:07.950642+00:00", "description": "Surface plasmon resonance binding", "method_type": "Direct", "typical_unit": "KD (M)", "detection_type": "Optical"}	armand	2025-01-08 03:46:07.950642+00
c89994f6-5b31-49e7-abd9-9ff2df888b18	binding_assay_types	3f05260d-5289-4355-a482-ea256f410fa8	INSERT	\N	{"id": "3f05260d-5289-4355-a482-ea256f410fa8", "name": "FRET", "created_at": "2025-01-08T03:46:07.950642+00:00", "updated_at": "2025-01-08T03:46:07.950642+00:00", "description": "Förster resonance energy transfer", "method_type": "Proximity", "typical_unit": "EC50 (nM)", "detection_type": "Fluorescence"}	armand	2025-01-08 03:46:07.950642+00
783033b4-f886-4b55-9ba9-90cb0ebdaa68	binding_assay_types	1d7a8df7-584a-4022-96a5-13d5416016e5	INSERT	\N	{"id": "1d7a8df7-584a-4022-96a5-13d5416016e5", "name": "BRET", "created_at": "2025-01-08T03:46:07.950642+00:00", "updated_at": "2025-01-08T03:46:07.950642+00:00", "description": "Bioluminescence resonance energy transfer", "method_type": "Proximity", "typical_unit": "EC50 (nM)", "detection_type": "Luminescence"}	armand	2025-01-08 03:46:07.950642+00
346b22cc-23a1-487d-9d94-c6edde60453d	social_effects	3fdb5d16-59c9-42d7-86ad-8b993b227ded	INSERT	\N	{"id": "3fdb5d16-59c9-42d7-86ad-8b993b227ded", "url": null, "type": "Visual", "created_at": "2025-01-08T03:46:08.024246+00:00", "description": "The experience of textures, surfaces, and objects appearing to move or flow", "effect_name": "Visual Drifting", "analysis_data": null, "related_effects": null}	armand	2025-01-08 03:46:08.024246+00
0b595057-468f-4e85-bd85-87b0946093c4	social_effects	060a4b57-d1df-4d46-b673-9723c3e37d69	INSERT	\N	{"id": "060a4b57-d1df-4d46-b673-9723c3e37d69", "url": null, "type": "Visual", "created_at": "2025-01-08T03:46:08.024246+00:00", "description": "The experience of seeing various geometric patterns and forms", "effect_name": "Geometric Patterns", "analysis_data": null, "related_effects": null}	armand	2025-01-08 03:46:08.024246+00
e4c13546-77b6-415f-a2c8-4915473a0293	social_effects	8c5c3eb6-1487-49b4-aae3-16e37188e79f	INSERT	\N	{"id": "8c5c3eb6-1487-49b4-aae3-16e37188e79f", "url": null, "type": "Cognitive", "created_at": "2025-01-08T03:46:08.024246+00:00", "description": "Alterations in the perception of time passing", "effect_name": "Time Distortion", "analysis_data": null, "related_effects": null}	armand	2025-01-08 03:46:08.024246+00
7d9f68a2-1106-4754-ae4e-fb59067a0102	social_effects	2bf0d809-7ce8-455e-a0db-d519d6110c91	INSERT	\N	{"id": "2bf0d809-7ce8-455e-a0db-d519d6110c91", "url": null, "type": "Physical/Cognitive", "created_at": "2025-01-08T03:46:08.024246+00:00", "description": "A state of intense happiness and well-being", "effect_name": "Euphoria", "analysis_data": null, "related_effects": null}	armand	2025-01-08 03:46:08.024246+00
13689d67-3d15-4fba-9b15-8158c1cf1201	social_effects	8637b2a6-910a-4608-95cb-c5cac98e17ca	INSERT	\N	{"id": "8637b2a6-910a-4608-95cb-c5cac98e17ca", "url": null, "type": "Auditory", "created_at": "2025-01-08T03:46:08.024246+00:00", "description": "Music sounds more detailed, meaningful, or emotionally impactful", "effect_name": "Enhanced Music Appreciation", "analysis_data": null, "related_effects": null}	armand	2025-01-08 03:46:08.024246+00
507ba983-6094-456b-9da3-c072494afc6d	social_effects	1d636ef8-02f0-4331-b7f6-6503ecbd26b9	INSERT	\N	{"id": "1d636ef8-02f0-4331-b7f6-6503ecbd26b9", "url": null, "type": "Cognitive", "created_at": "2025-01-08T03:46:08.024246+00:00", "description": "The experience of a decreased sense of self-identity", "effect_name": "Ego Dissolution", "analysis_data": null, "related_effects": null}	armand	2025-01-08 03:46:08.024246+00
137810a3-bf9e-452e-8c9d-a3d8470bc052	social_effects	4b403726-f6e7-4338-bf85-adff0603e994	INSERT	\N	{"id": "4b403726-f6e7-4338-bf85-adff0603e994", "url": null, "type": "Cognitive", "created_at": "2025-01-08T03:46:08.024246+00:00", "description": "The mixing of sensory modalities", "effect_name": "Synesthesia", "analysis_data": null, "related_effects": null}	armand	2025-01-08 03:46:08.024246+00
796692af-6b79-4964-8ac9-79751dcec5b7	social_effects	637a8080-1718-4ed2-89c0-a47ae0860245	INSERT	\N	{"id": "637a8080-1718-4ed2-89c0-a47ae0860245", "url": null, "type": "Physical", "created_at": "2025-01-08T03:46:08.024246+00:00", "description": "Increased sensitivity to physical touch and textures", "effect_name": "Enhanced Tactile Sensation", "analysis_data": null, "related_effects": null}	armand	2025-01-08 03:46:08.024246+00
d559214a-34c8-42f6-b729-3a12c5570e08	social_effects	b79560db-90a6-418c-aec7-35343c873cb9	INSERT	\N	{"id": "b79560db-90a6-418c-aec7-35343c873cb9", "url": null, "type": "Cognitive", "created_at": "2025-01-08T03:46:08.024246+00:00", "description": "Abstract thoughts become more vivid and meaningful", "effect_name": "Conceptual Thinking", "analysis_data": null, "related_effects": null}	armand	2025-01-08 03:46:08.024246+00
29da0961-e51b-4b99-98e7-ffcd4f0f29fb	social_effects	5f4f7070-5d23-4c92-971c-969a78b8a42c	INSERT	\N	{"id": "5f4f7070-5d23-4c92-971c-969a78b8a42c", "url": null, "type": "Visual", "created_at": "2025-01-08T03:46:08.024246+00:00", "description": "Improved clarity and sharpness of vision", "effect_name": "Visual Acuity Enhancement", "analysis_data": null, "related_effects": null}	armand	2025-01-08 03:46:08.024246+00
327ef275-064d-480b-8874-ae83c0d2facf	social_harm_reduction	3012fd35-b302-401e-bda9-ed43a8a9d665	INSERT	\N	{"id": "3012fd35-b302-401e-bda9-ed43a8a9d665", "title": "When to Call Emergency Services", "author": null, "content": "Call emergency services immediately if someone experiences: severe confusion, unconsciousness, difficulty breathing, seizures, severe overheating, or chest pain.", "post_id": null, "sources": null, "category": "Emergency", "platform": null, "warnings": null, "created_at": "2025-01-08T03:46:08.024795+00:00", "compound_id": null, "safety_notes": null, "importance_level": "Critical", "last_reviewed_at": null, "contraindications": null, "verification_notes": null, "created_at_internal": "2025-01-08T03:46:08.024795+00:00", "verification_status": "unverified", "emergency_procedures": null, "platform_subdivision": null}	armand	2025-01-08 03:46:08.024795+00
ba948f4b-8912-4c75-bc56-31f511a45bcb	social_harm_reduction	5297da67-67e5-4a08-9990-a53b46ebedda	INSERT	\N	{"id": "5297da67-67e5-4a08-9990-a53b46ebedda", "title": "Test Your Substances", "author": null, "content": "Always test your substances with multiple reagent tests. Never consume unidentified substances.", "post_id": null, "sources": null, "category": "General", "platform": null, "warnings": null, "created_at": "2025-01-08T03:46:08.024795+00:00", "compound_id": null, "safety_notes": null, "importance_level": "Critical", "last_reviewed_at": null, "contraindications": null, "verification_notes": null, "created_at_internal": "2025-01-08T03:46:08.024795+00:00", "verification_status": "unverified", "emergency_procedures": null, "platform_subdivision": null}	armand	2025-01-08 03:46:08.024795+00
8e2c7e49-5a81-45c5-8da4-96b08ab638bf	social_harm_reduction	779d2433-774f-41a5-8611-89ae6782ee8b	INSERT	\N	{"id": "779d2433-774f-41a5-8611-89ae6782ee8b", "title": "Start Low, Go Slow", "author": null, "content": "Always start with a low dose, especially with new substances or batches. Wait sufficient time before considering redosing.", "post_id": null, "sources": null, "category": "Dosage", "platform": null, "warnings": null, "created_at": "2025-01-08T03:46:08.024795+00:00", "compound_id": null, "safety_notes": null, "importance_level": "Critical", "last_reviewed_at": null, "contraindications": null, "verification_notes": null, "created_at_internal": "2025-01-08T03:46:08.024795+00:00", "verification_status": "unverified", "emergency_procedures": null, "platform_subdivision": null}	armand	2025-01-08 03:46:08.024795+00
82e963ad-e515-4a81-8c69-1d6b6e18c8b8	social_harm_reduction	be83b6c1-7ff6-497f-9015-35fde83f44b5	INSERT	\N	{"id": "be83b6c1-7ff6-497f-9015-35fde83f44b5", "title": "Safe Injection Practices", "author": null, "content": "Use clean equipment, never share needles, and practice proper hygiene. Know the proper injection techniques for harm reduction.", "post_id": null, "sources": null, "category": "ROA", "platform": null, "warnings": null, "created_at": "2025-01-08T03:46:08.024795+00:00", "compound_id": null, "safety_notes": null, "importance_level": "Critical", "last_reviewed_at": null, "contraindications": null, "verification_notes": null, "created_at_internal": "2025-01-08T03:46:08.024795+00:00", "verification_status": "unverified", "emergency_procedures": null, "platform_subdivision": null}	armand	2025-01-08 03:46:08.024795+00
4f070ee6-c921-4e98-a844-5824031f70bf	social_harm_reduction	b9279197-5d0b-45b9-9b47-1754c3257166	INSERT	\N	{"id": "b9279197-5d0b-45b9-9b47-1754c3257166", "title": "Avoid Dangerous Combinations", "author": null, "content": "Research interactions before combining substances. Many combinations can be unexpectedly dangerous.", "post_id": null, "sources": null, "category": "Combinations", "platform": null, "warnings": null, "created_at": "2025-01-08T03:46:08.024795+00:00", "compound_id": null, "safety_notes": null, "importance_level": "Critical", "last_reviewed_at": null, "contraindications": null, "verification_notes": null, "created_at_internal": "2025-01-08T03:46:08.024795+00:00", "verification_status": "unverified", "emergency_procedures": null, "platform_subdivision": null}	armand	2025-01-08 03:46:08.024795+00
d7e62765-7c5c-4ac4-b03f-dc44e4db1403	quantum_research_projects	1df598d0-cb75-4d9e-914f-536d0bfa5824	INSERT	\N	{"id": "1df598d0-cb75-4d9e-914f-536d0bfa5824", "name": "QM Structure Analysis", "budget": null, "status": "active", "end_date": null, "metadata": null, "created_at": "2025-01-08T03:46:08.143179+00:00", "objectives": null, "start_date": null, "updated_at": "2025-01-08T03:46:08.143179+00:00", "description": "Quantum mechanical analysis of molecular structures", "research_type": "computational", "funding_source": null, "principal_investigator": null}	armand	2025-01-08 03:46:08.143179+00
80decb0a-358b-4e96-b22e-95ed7e126d64	quantum_research_projects	15da953f-2dde-4cbb-8360-2185bfbd7bc8	INSERT	\N	{"id": "15da953f-2dde-4cbb-8360-2185bfbd7bc8", "name": "Electronic Properties", "budget": null, "status": "active", "end_date": null, "metadata": null, "created_at": "2025-01-08T03:46:08.143179+00:00", "objectives": null, "start_date": null, "updated_at": "2025-01-08T03:46:08.143179+00:00", "description": "Investigation of electronic structure properties", "research_type": "theoretical", "funding_source": null, "principal_investigator": null}	armand	2025-01-08 03:46:08.143179+00
e7e5783b-4de8-4a9e-b696-3d16e64e4953	quantum_research_projects	49fa405c-53e9-4b79-9e16-ddabfa42d504	INSERT	\N	{"id": "49fa405c-53e9-4b79-9e16-ddabfa42d504", "name": "Reaction Mechanisms", "budget": null, "status": "planned", "end_date": null, "metadata": null, "created_at": "2025-01-08T03:46:08.143179+00:00", "objectives": null, "start_date": null, "updated_at": "2025-01-08T03:46:08.143179+00:00", "description": "Quantum study of reaction pathways", "research_type": "computational", "funding_source": null, "principal_investigator": null}	armand	2025-01-08 03:46:08.143179+00
\.


--
-- Data for Name: binding_assay_protocols; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.binding_assay_protocols (id, name, description, assay_type, protocol_steps, reagents, equipment, controls, validation_criteria, limitations, reference_list, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: binding_assay_types; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.binding_assay_types (id, name, description, method_type, detection_type, typical_unit, created_at, updated_at) FROM stdin;
2a55ba47-c256-4db7-858e-3a0371ea6703	Radioligand Binding	Direct measurement using radioactive ligands	Competition	Scintillation	Ki (nM)	2025-01-08 03:46:07.950642+00	2025-01-08 03:46:07.950642+00
fad9a424-29a5-4a18-96bd-ffa430676e08	Fluorescence Binding	Fluorescence-based binding assays	Direct	Fluorescence	Kd (nM)	2025-01-08 03:46:07.950642+00	2025-01-08 03:46:07.950642+00
27c3d2f1-9821-4f7c-a311-44d1a8d805d8	SPR	Surface plasmon resonance binding	Direct	Optical	KD (M)	2025-01-08 03:46:07.950642+00	2025-01-08 03:46:07.950642+00
3f05260d-5289-4355-a482-ea256f410fa8	FRET	Förster resonance energy transfer	Proximity	Fluorescence	EC50 (nM)	2025-01-08 03:46:07.950642+00	2025-01-08 03:46:07.950642+00
1d7a8df7-584a-4022-96a5-13d5416016e5	BRET	Bioluminescence resonance energy transfer	Proximity	Luminescence	EC50 (nM)	2025-01-08 03:46:07.950642+00	2025-01-08 03:46:07.950642+00
\.


--
-- Data for Name: binding_data; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.binding_data (id, compound_id, receptor_family_id, assay_type_id, value, unit, confidence_score, data_source, publication_doi, experimental_conditions, measurement_date, binding_site, binding_mode, binding_kinetics, activity_type, activity_value, activity_unit, efficacy, potency, assay_description, assay_organism, assay_type, assay_conditions, assay_cell_line, assay_tissue_type, sar_data, pharmacophore_features, binding_pocket_residues, experimental_method, validation_method, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: binding_data_quality; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.binding_data_quality (id, binding_data_id, replicate_count, standard_deviation, confidence_interval, quality_score, validation_notes, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: binding_kinetics; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.binding_kinetics (id, binding_data_id, kon_value, kon_unit, koff_value, koff_unit, residence_time, residence_time_unit, temperature, ph, ionic_strength, buffer_conditions, method_details, equipment_used, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: binding_sar; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.binding_sar (id, compound_id, receptor_family_id, structural_feature, effect_type, effect_magnitude, confidence_score, evidence_type, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: binding_site_mapping; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.binding_site_mapping (id, binding_data_id, site_name, residues, interaction_types, binding_pocket_volume, surface_accessibility, conservation_score, mutation_effects, structural_features, modeling_method, confidence_score, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: carcinogenicity_data; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.carcinogenicity_data (id, compound_id, study_type, species, duration, dose_levels, tumor_types, mechanism, histopathology, reference_doi, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: cardiotoxicity_data; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.cardiotoxicity_data (id, compound_id, study_type, cardiac_effects, mechanism, herg_ic50, ecg_changes, hemodynamic_effects, reference_doi, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: compounds; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.compounds (id, name, smiles, inchi, inchi_key, cas_number, pubchem_cid, chembl_id, drugbank_id, unii, kegg_id, chemspider_id, zinc_id, chebi_id, iupac_name, preferred_iupac_name, common_names, einecs_number, rtecs_number, hsdb_number, ccdc_number, reaxys_id, lipidmaps_id, nsc_number, qcarchive_id, nomad_id, materials_project_id, basis_set_id, method_id, calculation_id, molecular_weight, molecular_formula, stereochemistry, crystal_structure, solubility_data, pka_values, partition_coefficients, surface_properties, conformational_analysis, chirality_info, isomer_details, melting_point, boiling_point, density, refractive_index, optical_rotation, drug_class, mechanism_categories, therapeutic_categories, pharmacological_effects, administration_routes, bioavailability_data, metabolism_data, distribution_data, legal_status, scheduling_info, approval_status, clinical_trial_status, patent_status, registration_numbers, control_status, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: cytotoxicity_data; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.cytotoxicity_data (id, compound_id, cell_line, tissue_type, assay_method, exposure_time, ic50_value, ic50_unit, cell_viability, cytotoxicity_mechanism, morphological_changes, reference_doi, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: descriptors_2d; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.descriptors_2d (id, compound_id, wiener_index, balaban_j_index, bertz_ct, schultz_molecular_topological_index, chi_0_index, chi_1_index, chi_2_index, chi_3_index, sum_estate_indices, mean_estate_indices, kappa_1, kappa_2, kappa_3, molar_refractivity, van_der_waals_volume, polarizability, formal_charge, ring_count, aromatic_ring_count, aliphatic_ring_count, ring_fusion_degree, rotatable_bond_count, rigid_bond_count, chain_atom_count, chain_bond_count, carbon_count, nitrogen_count, oxygen_count, sulfur_count, phosphorus_count, halogen_count, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: descriptors_3d; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.descriptors_3d (id, compound_id, conformer_id, radius_of_gyration, molecular_volume, molecular_surface_area, solvent_accessible_surface_area, polar_surface_area_3d, spherosity, asphericity, eccentricity, inertial_shape_factor, principal_moment_1, principal_moment_2, principal_moment_3, gravitational_index, radius_of_distribution, molecular_surface_potential, average_surface_charge, total_energy, strain_energy, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: developmental_toxicity; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.developmental_toxicity (id, compound_id, study_type, species, exposure_period, dose_levels, developmental_effects, mechanism, teratogenicity, reference_doi, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: device_settings; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.device_settings (id, device_id, device_type, os_version, app_version, screen_size, capabilities, settings, performance_profile, security_settings, network_preferences, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: effect_relationships; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.effect_relationships (id, effect_id, related_effect_id, relationship_type, strength, evidence_level, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: electronic_structure; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.electronic_structure (id, compound_id, total_electronic_energy, electron_correlation_energy, exchange_energy, kinetic_energy, homo_lumo_gap, orbital_energies, orbital_occupancies, density_matrix, density_grid_points, density_values, wavefunction, state_type, band_structure, method, basis_set, convergence_criteria, created_at, updated_at, metadata) FROM stdin;
\.


--
-- Data for Name: emergency_procedures; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.emergency_procedures (id, compound_id, spill_response, fire_fighting_measures, first_aid_procedures, evacuation_criteria, emergency_contacts, special_hazards, cleanup_procedures, disposal_procedures, reporting_requirements, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: energy_level_statistics; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.energy_level_statistics (id, compound_id, energy_levels, level_spacings, spacing_distribution, cumulative_distribution, distribution_type, gamma_parameter, confidence_score, analysis_parameters, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: experience_categories; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.experience_categories (id, category_name, parent_category, description, report_count, created_at) FROM stdin;
ad435445-a5a7-4461-85fd-e99472597a04	First Time	\N	First experiences with a substance	0	2025-01-08 03:46:08.023842+00
410fb348-8060-4e23-a377-3638405fc095	Difficult Experience	\N	Challenging or negative experiences	0	2025-01-08 03:46:08.023842+00
9636613c-349b-48a0-8afa-535c196b990f	Medical Use	\N	Use for medical or therapeutic purposes	0	2025-01-08 03:46:08.023842+00
0f85d5cf-24b8-4a06-913e-d8f2a726a5d8	Spiritual	\N	Experiences with spiritual or mystical elements	0	2025-01-08 03:46:08.023842+00
c742e24a-a8ed-424b-a4ac-8ff251334d89	Train Wrecks & Trip Disasters	\N	Particularly difficult or dangerous experiences	0	2025-01-08 03:46:08.023842+00
4ad4b80b-5a3f-4652-9ea3-b54621c20039	Health Problems	\N	Experiences involving health issues	0	2025-01-08 03:46:08.023842+00
362a6433-f937-4ecd-ad66-370c1fece302	Retrospective / Summary	\N	Looking back on past experiences	0	2025-01-08 03:46:08.023842+00
30ce678c-0e08-4dda-a643-027f5dc351d1	Small Collection	\N	Brief or collected experiences	0	2025-01-08 03:46:08.023842+00
0c7e5ecf-5056-4443-8260-e8d3e456ee03	Master/Teacher Plants	\N	Experiences with traditional plant medicines	0	2025-01-08 03:46:08.023842+00
c4109e4d-5f72-4024-b184-017b9f01aad2	Combinations	\N	Experiences with multiple substances	0	2025-01-08 03:46:08.023842+00
871d2fd7-f580-4a2d-b097-82181214c2be	Preparation / Recipes	\N	Methods of preparation or consumption	0	2025-01-08 03:46:08.023842+00
5b1007da-1ebb-492b-b3bc-7db8ace6b5f9	Bad Trips	\N	Specifically difficult psychological experiences	0	2025-01-08 03:46:08.023842+00
e7528b0c-ac26-4761-8320-2fa0f375ce98	Glowing Experiences	\N	Particularly positive experiences	0	2025-01-08 03:46:08.023842+00
420540e2-a693-4998-9a43-af5a570225f0	Clinical Research	\N	Experiences in research settings	0	2025-01-08 03:46:08.023842+00
\.


--
-- Data for Name: external_services; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.external_services (id, name, service_type, base_url, auth_config, rate_limit, timeout, retry_config, is_active, health_check_path, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: feature_definitions; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.feature_definitions (id, name, description, feature_type, data_type, calculation_method, dependencies, validation_rules, metadata, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: feature_pipelines; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.feature_pipelines (id, name, description, steps, input_features, output_features, parameters, validation_rules, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: feature_selection; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.feature_selection (id, model_id, selection_method, selected_features, importance_scores, selection_criteria, validation_metrics, selection_date, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: feature_values; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.feature_values (id, compound_id, feature_id, value, calculation_date, confidence_score, metadata, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: genes; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.genes (id, name, symbol, description, organism, chromosome, location, gene_type, sequence, alternative_symbols, ensembl_id, entrez_id, hgnc_id, mgi_id, rgd_id, uniprot_ids, refseq_ids, regulatory_elements, expression_data, pathway_involvement, disease_associations, gene_ontology, evolutionary_conservation, variants, interactions, literature_references, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: genotoxicity_data; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.genotoxicity_data (id, compound_id, test_type, organism, metabolic_activation, result, mutation_type, dna_damage_type, mechanism, reference_doi, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: hazard_classifications; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.hazard_classifications (id, compound_id, ghs_classifications, signal_word, pictograms, hazard_statements, precautionary_statements, nfpa_health, nfpa_fire, nfpa_reactivity, nfpa_special, hmis_health, hmis_fire, hmis_physical, hmis_ppe, classification_source, classification_date, review_date, notes, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: hepatotoxicity_data; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.hepatotoxicity_data (id, compound_id, study_type, liver_effects, mechanism, enzyme_elevations, histopathology, clinical_significance, reference_doi, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: immunotoxicity_data; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.immunotoxicity_data (id, compound_id, study_type, immune_parameters, effect_type, mechanism, clinical_relevance, reference_doi, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: literature_findings; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.literature_findings (id, paper_id, compound_id, finding_type, finding_details, statistical_significance, confidence_interval, methodology_notes, limitations, replication_status, validation_method, supporting_evidence, contradicting_evidence, clinical_relevance, research_implications, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: meta_analyses; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.meta_analyses (id, title, topic, included_papers, methodology, total_sample_size, pooled_effect_size, heterogeneity_metrics, subgroup_analyses, sensitivity_analyses, publication_bias_assessment, quality_assessment_method, evidence_strength, clinical_implications, research_gaps, conclusions, limitations, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: metabolic_toxicity; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.metabolic_toxicity (id, compound_id, enzyme_affected, inhibition_type, ki_value, ki_unit, metabolites, pathway_disruption, clinical_significance, reference_doi, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: ml_models; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.ml_models (id, name, version, model_type, description, target_variable, features, hyperparameters, architecture, preprocessing_steps, training_config, metadata, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: model_deployments; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.model_deployments (id, model_id, deployment_name, deployment_environment, deployment_date, status, version_tag, configuration, performance_metrics, monitoring_config, rollback_info, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: model_metrics; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.model_metrics (id, model_id, version_id, metric_type, metric_value, metric_date, dataset_info, calculation_method, confidence_interval, metadata, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: model_monitoring; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.model_monitoring (id, model_id, monitoring_date, metric_name, metric_value, threshold_value, alert_status, data_drift_metrics, performance_metrics, resource_usage, alert_history, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: model_predictions; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.model_predictions (id, model_id, version_id, compound_id, prediction_type, predicted_value, confidence_score, prediction_date, input_features, explanation, metadata, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: model_validation_results; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.model_validation_results (id, model_id, validation_type, validation_date, test_dataset_info, validation_metrics, test_set_performance, cross_validation_results, error_analysis, validation_plots, validation_notes, recommendations, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: model_versions; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.model_versions (id, model_id, version_number, changes_description, performance_metrics, validation_results, deployment_status, deployed_at, deprecated_at, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: molecular_fingerprints; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.molecular_fingerprints (id, compound_id, maccs_keys, pubchem_bits, ecfp4_bits, ecfp6_bits, daylight_bits, pharma_bits, atom_pairs, torsion_bits, morgan_bits, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: monitoring_parameters; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.monitoring_parameters (id, name, category, frequency, monitoring_method, alert_conditions, normal_range, data_type, created_at, updated_at) FROM stdin;
4383f91a-068f-45d2-aced-9e51a8f18d00	Liver Function	biochemical	1 week	Blood test	{"ALT > 3x ULN","AST > 3x ULN","ALP > 2x ULN"}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
bd7b361c-45af-4160-bab1-498345c872c7	Kidney Function	biochemical	1 week	Blood test	{"Creatinine > 1.5x baseline","eGFR decrease > 25%"}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
38d2856e-16bd-4718-a585-d07865b665aa	Cardiac Function	physiological	1 day	ECG	{"QTc > 500ms","QTc increase > 60ms"}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
c40398a7-ee8f-43b7-91db-579438add951	Blood Pressure	physiological	6 hours	Automated measurement	{"Systolic > 160","Diastolic > 100"}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
1b315f7b-1063-4873-992c-ea504ff5fd56	Blood Count	hematological	1 week	Blood test	{"WBC < 3000/µL","Platelets < 100k/µL"}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
d948d10d-0ae3-4f75-b0a3-c41332848de9	Mental Status	neurological	6 hours	Clinical assessment	{Confusion,Agitation,Drowsiness}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
0537362d-e6f3-424c-b692-8aafe2ef9d23	Body Temperature	physiological	6 hours	Temperature measurement	{"> 38.5°C","< 35.5°C"}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
f1b972ce-a75e-49f3-bd31-e8c2f82627c1	Respiratory Rate	physiological	6 hours	Clinical measurement	{"> 24/min","< 8/min"}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
1228783f-7e7c-4080-a706-9fddb5f4440e	Oxygen Saturation	physiological	6 hours	Pulse oximetry	{"< 92%","Drop > 4% from baseline"}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
23620d18-ee65-4892-923a-a10f0d7ffb25	Cognitive Function	neurological	1 day	Cognitive assessment	{"MMSE decrease > 2 points","New onset confusion"}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
9dc34c63-f29b-405d-8475-9fce11f6f6f9	Sleep Pattern	behavioral	1 day	Sleep diary	{"Insomnia > 2 hours","Excessive drowsiness"}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
50c33a74-e1d8-4e68-9d2e-6009771f1b31	Appetite	behavioral	1 day	Food intake log	{"Decrease > 50%","Complete loss of appetite"}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
a712742c-7c4a-4268-a368-34b0f0d28b07	Mood	psychological	1 day	Mood scale	{"Severe depression",Mania,Anxiety}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
76adec8a-5eb7-40a3-8798-2484633eb230	Movement	neurological	1 day	Clinical observation	{Tremor,Ataxia,Dystonia}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
b2ef6aef-9456-4263-815c-be1311724e3e	Pain Level	subjective	6 hours	Pain scale	{"Score > 7/10","Acute increase > 3 points"}	\N	\N	2025-01-08 03:46:07.899239+00	2025-01-08 03:46:07.899239+00
\.


--
-- Data for Name: neurotoxicity_data; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.neurotoxicity_data (id, compound_id, study_type, brain_regions_affected, behavioral_effects, cellular_effects, mechanism, reversibility, long_term_effects, reference_doi, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: offline_data; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.offline_data (id, device_id, data_type, data_id, content, version, priority, compression_type, encryption_type, validation_hash, expires_at, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: organ_systems; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.organ_systems (id, name, description, major_components, key_functions, vulnerability_factors, assessment_methods, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: organ_toxicity_patterns; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.organ_toxicity_patterns (id, organ_system_id, name, description, cellular_targets, molecular_mechanisms, histological_changes, functional_impacts, early_biomarkers, diagnostic_markers, progression_pattern, reversibility_potential, risk_factors, protective_factors, monitoring_parameters, intervention_thresholds, treatment_approaches, prevention_strategies, research_status, evidence_level, notes, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: pathway_analysis; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.pathway_analysis (id, compound_id, pathway_name, pathway_type, affected_proteins, regulation_effects, downstream_effects, feedback_mechanisms, pathway_crosstalk, temporal_dynamics, tissue_specificity, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: performance_metrics; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.performance_metrics (id, metric_type, component, value, unit, threshold, status, metadata, tags, alert_triggered, resolution_steps, measured_at, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: pharmacological_classes; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.pharmacological_classes (id, name, description, mechanism_type, target_systems, typical_effects, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: pharmacophore_features; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.pharmacophore_features (id, compound_id, h_bond_donors_count, h_bond_acceptors_count, pos_charge_groups_count, neg_charge_groups_count, aromatic_rings_count, hydrophobic_groups_count, donor_positions, acceptor_positions, charge_positions, aromatic_positions, hydrophobic_positions, donor_strengths, acceptor_strengths, charge_strengths, feature_type, coordinates, strength, interaction_radius, optional, detection_method, confidence_score, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: phase_transition_types; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.phase_transition_types (id, name, description, transition_order, critical_exponents, universality_class, characteristic_properties) FROM stdin;
1	Continuous	Second-order phase transition	2	{"nu": 0.63, "beta": 0.325, "alpha": 0.110, "delta": 4.82, "gamma": 1.24}	Ising	\N
2	Discontinuous	First-order phase transition	1	\N	\N	\N
3	BCS	Superconducting transition	2	{"nu": 0.5, "beta": 0.5, "alpha": 0, "delta": 3.0, "gamma": 1.0}	Mean Field	\N
4	BEC	Bose-Einstein condensation	2	{"nu": 0.5, "beta": 0.5, "alpha": -1, "delta": 3.0, "gamma": 1.0}	Gaussian	\N
\.


--
-- Data for Name: phase_transitions; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.phase_transitions (id, compound_id, transition_type, critical_temperature, critical_pressure, order_parameter, correlation_length, transition_order, hysteresis_data, fluctuation_data, created_at, updated_at, metadata) FROM stdin;
\.


--
-- Data for Name: ppe_requirements; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.ppe_requirements (id, compound_id, eye_protection, skin_protection, respiratory_protection, hand_protection, body_protection, minimum_ppe_rating, special_requirements, exposure_limits, monitoring_requirements, decontamination_procedures, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: proteins; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.proteins (id, gene_id, name, symbol, description, organism, sequence, length, molecular_weight, uniprot_id, pdb_ids, refseq_ids, alternative_names, protein_family, domains, motifs, subcellular_location, post_translational_modifications, structure_data, function_data, interactions, expression_pattern, regulatory_mechanisms, disease_associations, drug_interactions, pathway_involvement, literature_references, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: quantum_basis_sets; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.quantum_basis_sets (id, name, description, basis_type, elements, accuracy_level, computational_cost, reference_citation) FROM stdin;
1	STO-3G	Minimal basis set using 3 Gaussian functions	minimal	{H,C,N,O,F,P,S,Cl}	low	very low	\N
2	6-31G	Split valence basis set	double-zeta	{H,C,N,O,F,P,S,Cl}	medium	medium	\N
3	cc-pVDZ	Correlation consistent double-zeta basis	double-zeta	{H,C,N,O,F,P,S,Cl}	high	high	\N
4	cc-pVTZ	Correlation consistent triple-zeta basis	triple-zeta	{H,C,N,O,F,P,S,Cl}	very high	very high	\N
\.


--
-- Data for Name: quantum_calculations; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.quantum_calculations (id, project_id, compound_id, calculation_type, basis_set, functional, calculation_status, start_timestamp, end_timestamp, cpu_hours, memory_gb, convergence_achieved, energy_hartree, energy_gradient_norm, calculation_parameters, output_files, metadata, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: quantum_critical_params; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.quantum_critical_params (id, compound_id, critical_temperature, critical_pressure, critical_field, primary_order_parameter, order_parameter_values, alpha, beta, gamma, delta, nu, eta, correlation_length, correlation_function, dynamic_exponent_z, phase_boundaries, multicriticality_type, coherence_length, entanglement_entropy, quantum_fluctuations, created_at, updated_at, metadata) FROM stdin;
\.


--
-- Data for Name: quantum_dynamics; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.quantum_dynamics (id, compound_id, time_points, wavefunction_evolution, density_matrix_evolution, coherence_times, relaxation_rates, dephasing_rates, conductivity_tensor, hall_conductance, thermal_conductivity, spectral_function, optical_conductivity, entanglement_spectrum, mutual_information, dissipation_kernel, noise_spectrum, metadata, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: quantum_functionals; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.quantum_functionals (id, name, description, functional_type, properties_handled, accuracy_metrics, computational_cost, reference_citation) FROM stdin;
1	LDA	Local Density Approximation	LDA	{exchange,correlation}	\N	low	\N
2	PBE	Perdew-Burke-Ernzerhof	GGA	{exchange,correlation}	\N	medium	\N
3	B3LYP	Becke 3-parameter Lee-Yang-Parr	hybrid	{exchange,correlation}	\N	high	\N
4	M06-2X	Minnesota 06 functional with 2X exchange	hybrid-meta-GGA	{exchange,correlation,dispersion}	\N	very high	\N
\.


--
-- Data for Name: quantum_hamiltonians; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.quantum_hamiltonians (id, compound_id, eh_hamiltonian, overlap_matrix, transformed_hamiltonian, disorder_strength, is_critical, calculation_method, calculation_parameters, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: quantum_observables; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.quantum_observables (id, compound_id, observable_type, value, uncertainty, measurement_basis, operator_type, expectation_value, variance, created_at, updated_at, metadata) FROM stdin;
\.


--
-- Data for Name: quantum_observables_ref; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.quantum_observables_ref (id, name, description, operator_type, measurement_units, uncertainty_type, standard_deviation_typical, measurement_protocol) FROM stdin;
1	Energy	Total electronic energy	energy	Hartree	absolute	\N	\N
2	Dipole Moment	Electric dipole moment	electromagnetic	Debye	relative	\N	\N
3	Spin	Electron spin	angular momentum	ℏ/2	discrete	\N	\N
4	Electron Density	Probability density of electrons	density	e/Å³	statistical	\N	\N
\.


--
-- Data for Name: quantum_parameters; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.quantum_parameters (id, parameter_name, description, unit, calculation_method, typical_range, accuracy_metrics, validation_criteria, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: quantum_properties; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.quantum_properties (id, compound_id, total_electronic_energy, homo_lumo_gap, electron_density, orbital_energies, critical_temperature, critical_pressure, correlation_length, order_parameter, coherence_time, relaxation_rate, dephasing_rate, transition_type, transition_order, universality_class, basis_set, method, calculation_parameters, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: quantum_research_findings; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.quantum_research_findings (id, calculation_id, finding_type, description, numerical_value, units, confidence_level, methodology, validation_method, publication_reference, discovery_date, significance_level, metadata, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: quantum_research_projects; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.quantum_research_projects (id, name, description, start_date, end_date, status, principal_investigator, research_type, funding_source, budget, objectives, metadata, created_at, updated_at) FROM stdin;
1df598d0-cb75-4d9e-914f-536d0bfa5824	QM Structure Analysis	Quantum mechanical analysis of molecular structures	\N	\N	active	\N	computational	\N	\N	\N	\N	2025-01-08 03:46:08.143179+00	2025-01-08 03:46:08.143179+00
15da953f-2dde-4cbb-8360-2185bfbd7bc8	Electronic Properties	Investigation of electronic structure properties	\N	\N	active	\N	theoretical	\N	\N	\N	\N	2025-01-08 03:46:08.143179+00	2025-01-08 03:46:08.143179+00
49fa405c-53e9-4b79-9e16-ddabfa42d504	Reaction Mechanisms	Quantum study of reaction pathways	\N	\N	planned	\N	computational	\N	\N	\N	\N	2025-01-08 03:46:08.143179+00	2025-01-08 03:46:08.143179+00
\.


--
-- Data for Name: quantum_structure_correlations; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.quantum_structure_correlations (id, compound_id, property_type, property_value, correlation_type, correlation_coefficient, statistical_significance, sample_size, methodology, validation_metrics, metadata, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: rate_limits; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.rate_limits (id, api_key_id, endpoint_id, requests_count, window_start, window_end, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: receptor_families; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.receptor_families (id, category_id, name, abbreviation, description, protein_type, signaling_type, primary_endogenous_ligands, primary_effects, therapeutic_areas, expression_pattern, signaling_pathways, pharmacological_properties, clinical_significance, research_status, notes, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: receptor_family_categories; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.receptor_family_categories (id, name, description, receptor_type, signaling_mechanism, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: receptor_subtypes; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.receptor_subtypes (id, family_id, subtype_name, description, protein_sequence, species, expression_pattern, signaling_pathways, pharmacological_profile, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: research_findings; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.research_findings (id, compound_id, finding_type, description, methodology, experimental_data, statistical_analysis, conclusions, limitations, future_directions, reference_dois, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: safety_documents; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.safety_documents (id, compound_id, sds_url, msds_url, sds_last_updated, sds_provider, sds_version, safety_data_sheet, handling_precautions, storage_precautions, disposal_instructions, first_aid_measures, firefighting_measures, accidental_release_measures, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: safety_thresholds; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.safety_thresholds (id, parameter_name, threshold_value, unit, severity_level, description, intervention_required, validation_method, created_at, updated_at) FROM stdin;
3bf67437-2d90-46dd-bffc-385d3f557b65	ALT	40	U/L	warning	Upper limit for alanine aminotransferase	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
ce3b5d75-6812-4828-a84d-042e4e37af34	AST	40	U/L	warning	Upper limit for aspartate aminotransferase	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
429b2523-f4d3-4887-8582-a4cbdaef2ad1	Creatinine	1.2	mg/dL	warning	Upper limit for serum creatinine	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
3ad82c61-2416-4cd2-94a4-fc0b3595c477	QTc	450	ms	warning	Upper limit for corrected QT interval	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
caafe287-1d2d-4974-87c5-138fd2ea69ad	Neutrophils	1500	cells/µL	warning	Lower limit for neutrophil count	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
0fcde142-ffd8-4aed-905b-74243b81a24f	Platelets	150000	cells/µL	warning	Lower limit for platelet count	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
c328858b-5397-404c-b595-5b2d748109ab	Heart Rate	100	bpm	warning	Upper limit for resting heart rate	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
228faec7-17b0-468b-9c43-26ce498b5952	Blood Pressure	140	mmHg	warning	Upper limit for systolic blood pressure	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
a02ee0bf-0a9c-4a4e-83f0-43d6e7172b76	Body Temperature	38.3	°C	warning	Upper limit for body temperature	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
2dd5f5fe-c66e-4f6c-9fd5-bb33460ceca0	Respiratory Rate	20	breaths/min	warning	Upper limit for respiratory rate	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
d491766b-d67e-4ba3-8415-5dc43c907fc2	Glucose	126	mg/dL	warning	Upper limit for fasting glucose	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
aaad45fe-3795-46b9-9cad-80492282728a	Total Bilirubin	1.2	mg/dL	warning	Upper limit for total bilirubin	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
ec7d6b88-dabe-403f-946f-bed13915cb82	Albumin	3.5	g/dL	warning	Lower limit for serum albumin	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
21156daf-7b2b-4109-8fc0-28bd5ab0421b	eGFR	60	mL/min	warning	Lower limit for estimated glomerular filtration rate	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
c2ab03ac-ab3a-4dfb-a376-b2b833387a3a	Oxygen Saturation	95	%	warning	Lower limit for oxygen saturation	f	\N	2025-01-08 03:46:07.896988+00	2025-01-08 03:46:07.896988+00
\.


--
-- Data for Name: sar_analysis; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.sar_analysis (id, compound_id, analysis_type, structural_features, activity_correlations, pharmacophore_model, binding_patterns, selectivity_patterns, structure_modifications, predicted_effects, confidence_metrics, validation_results, reference_compounds, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: sar_patterns; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.sar_patterns (id, pattern_type, structural_elements, activity_impact, confidence_score, supporting_compounds, detection_method, validation_status, notes, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: scaffold_analysis; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.scaffold_analysis (id, scaffold_smiles, compound_count, average_activity, activity_range, diversity_score, important_substitutions, analysis_method, notes, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: scaling_analysis; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.scaling_analysis (id, compound_id, phase_transition_id, scaling_function_type, scaling_variables, scaling_dimensions, rg_flow_equations, fixed_points, relevant_operators, universality_class, central_charge, operator_spectrum, size_scaling_exponents, correction_exponents, crossover_scales, crossover_functions, analysis_method, confidence_metrics, created_at, updated_at, metadata) FROM stdin;
\.


--
-- Data for Name: schema_versions; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.schema_versions (id, version, description, applied_at, applied_by, script_name, checksum, execution_time, status, error_message) FROM stdin;
6cd8cb4b-94b3-4ba1-b242-a2dfa5c69304	0.1.0	Initial core schema setup	2025-01-08 03:46:07.874039+00	armand	00_core.sql	\N	00:00:00.000088	SUCCESS	\N
\.


--
-- Data for Name: scientific_papers; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.scientific_papers (id, doi, title, authors, journal, publication_date, abstract, full_text, methodology, study_type, sample_size, study_duration, quality_metrics, evidence_level, key_findings, limitations, compounds_studied, validation_status, peer_review_status, citation_count, impact_factor, external_links, supplementary_data, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: service_status; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.service_status (id, service_id, status, response_time, error_message, check_timestamp, metrics, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: social_alert_rules; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_alert_rules (id, name, description, alert_type, severity, platforms, platform_subdivisions, compounds, trigger_conditions, required_metrics, threshold_values, cooldown_period, is_active, last_triggered_at, created_at) FROM stdin;
\.


--
-- Data for Name: social_comments; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_comments (id, post_id, platform, platform_subdivision, external_id, parent_id, content, author_id, author_username, created_at, engagement_metrics, classification_data, sentiment_data, is_scientific, has_citations, reported_effects, reported_side_effects, platform_specific_data, created_at_internal) FROM stdin;
\.


--
-- Data for Name: social_compound_combinations; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_compound_combinations (id, compound_id_1, compound_id_2, interaction_type, risk_level, description, mechanism, evidence_level, sources, platforms, reported_count, reports, created_at) FROM stdin;
\.


--
-- Data for Name: social_compound_mentions; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_compound_mentions (id, compound_id, platform, platform_subdivision, analysis_date, total_mentions, unique_authors, total_engagement, scientific_mentions, experience_reports, harm_reduction_mentions, sentiment_distribution, topic_distribution, user_demographics, geographic_distribution, temporal_patterns, common_contexts, related_compounds, platform_specific_metrics, cross_platform_engagement_flow, user_influence_metrics, content_propagation_patterns, created_at) FROM stdin;
\.


--
-- Data for Name: social_compound_trends; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_compound_trends (id, compound_id, trend_start_date, trend_end_date, platforms, total_mentions, unique_authors, total_engagement, trend_velocity, trend_acceleration, platform_distribution, topic_evolution, sentiment_evolution, influencer_impact, news_correlation, research_correlation, correlation_strength, context_similarity, user_overlap_patterns, temporal_similarity, platform_specific_metrics, trend_confidence_score, trend_validation_metrics, cross_platform_correlations, trend_seasonality, trend_anomaly_scores, trend_prediction_metrics, created_at) FROM stdin;
\.


--
-- Data for Name: social_content_quality; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_content_quality (id, compound_id, platform, platform_subdivision, analysis_date, scientific_accuracy_score, harm_reduction_quality_score, experience_report_quality_score, information_completeness_score, citation_quality_score, misinformation_prevalence_score, content_depth_distribution, quality_trends, improvement_suggestions, created_at) FROM stdin;
\.


--
-- Data for Name: social_dosage_data; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_dosage_data (id, compound_id, platform, platform_subdivision, route_of_administration, threshold_dose, light_dose, common_dose, strong_dose, heavy_dose, warning_dose, duration_total, duration_onset, duration_peak, duration_offset, bioavailability, dosage_notes, created_at) FROM stdin;
\.


--
-- Data for Name: social_effect_reports; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_effect_reports (id, compound_id, effect_id, report_text, intensity, duration, onset_time, conditions, platform, platform_subdivision, created_at) FROM stdin;
\.


--
-- Data for Name: social_effects; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_effects (id, effect_name, description, type, url, analysis_data, related_effects, created_at) FROM stdin;
3fdb5d16-59c9-42d7-86ad-8b993b227ded	Visual Drifting	The experience of textures, surfaces, and objects appearing to move or flow	Visual	\N	\N	\N	2025-01-08 03:46:08.024246+00
060a4b57-d1df-4d46-b673-9723c3e37d69	Geometric Patterns	The experience of seeing various geometric patterns and forms	Visual	\N	\N	\N	2025-01-08 03:46:08.024246+00
8c5c3eb6-1487-49b4-aae3-16e37188e79f	Time Distortion	Alterations in the perception of time passing	Cognitive	\N	\N	\N	2025-01-08 03:46:08.024246+00
2bf0d809-7ce8-455e-a0db-d519d6110c91	Euphoria	A state of intense happiness and well-being	Physical/Cognitive	\N	\N	\N	2025-01-08 03:46:08.024246+00
8637b2a6-910a-4608-95cb-c5cac98e17ca	Enhanced Music Appreciation	Music sounds more detailed, meaningful, or emotionally impactful	Auditory	\N	\N	\N	2025-01-08 03:46:08.024246+00
1d636ef8-02f0-4331-b7f6-6503ecbd26b9	Ego Dissolution	The experience of a decreased sense of self-identity	Cognitive	\N	\N	\N	2025-01-08 03:46:08.024246+00
4b403726-f6e7-4338-bf85-adff0603e994	Synesthesia	The mixing of sensory modalities	Cognitive	\N	\N	\N	2025-01-08 03:46:08.024246+00
637a8080-1718-4ed2-89c0-a47ae0860245	Enhanced Tactile Sensation	Increased sensitivity to physical touch and textures	Physical	\N	\N	\N	2025-01-08 03:46:08.024246+00
b79560db-90a6-418c-aec7-35343c873cb9	Conceptual Thinking	Abstract thoughts become more vivid and meaningful	Cognitive	\N	\N	\N	2025-01-08 03:46:08.024246+00
5f4f7070-5d23-4c92-971c-969a78b8a42c	Visual Acuity Enhancement	Improved clarity and sharpness of vision	Visual	\N	\N	\N	2025-01-08 03:46:08.024246+00
\.


--
-- Data for Name: social_experience_reports; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_experience_reports (id, post_id, compound_id, stack_id, protocol_id, title, content, author, platform, platform_subdivision, experience_date, report_date, substance_data, duration, setting_data, intention, effects_timeline, timeline_data, reported_effects, side_effects, interactions, test_kit_info, detection_time_data, after_effects_data, body_weight, weight_unit, gender, age, experience_level, harm_reduction_notes, classification_data, sentiment_data, report_version, experience_category, total_views, report_quality_score, medical_conditions, medications, baseline_metrics, outcome_metrics, testing_methods, overall_rating, platform_specific_data, created_at) FROM stdin;
\.


--
-- Data for Name: social_harm_reduction; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_harm_reduction (id, post_id, compound_id, title, content, author, platform, platform_subdivision, created_at, category, importance_level, safety_notes, warnings, contraindications, emergency_procedures, sources, verification_status, verification_notes, last_reviewed_at, created_at_internal) FROM stdin;
3012fd35-b302-401e-bda9-ed43a8a9d665	\N	\N	When to Call Emergency Services	Call emergency services immediately if someone experiences: severe confusion, unconsciousness, difficulty breathing, seizures, severe overheating, or chest pain.	\N	\N	\N	2025-01-08 03:46:08.024795+00	Emergency	Critical	\N	\N	\N	\N	\N	unverified	\N	\N	2025-01-08 03:46:08.024795+00
5297da67-67e5-4a08-9990-a53b46ebedda	\N	\N	Test Your Substances	Always test your substances with multiple reagent tests. Never consume unidentified substances.	\N	\N	\N	2025-01-08 03:46:08.024795+00	General	Critical	\N	\N	\N	\N	\N	unverified	\N	\N	2025-01-08 03:46:08.024795+00
779d2433-774f-41a5-8611-89ae6782ee8b	\N	\N	Start Low, Go Slow	Always start with a low dose, especially with new substances or batches. Wait sufficient time before considering redosing.	\N	\N	\N	2025-01-08 03:46:08.024795+00	Dosage	Critical	\N	\N	\N	\N	\N	unverified	\N	\N	2025-01-08 03:46:08.024795+00
be83b6c1-7ff6-497f-9015-35fde83f44b5	\N	\N	Safe Injection Practices	Use clean equipment, never share needles, and practice proper hygiene. Know the proper injection techniques for harm reduction.	\N	\N	\N	2025-01-08 03:46:08.024795+00	ROA	Critical	\N	\N	\N	\N	\N	unverified	\N	\N	2025-01-08 03:46:08.024795+00
b9279197-5d0b-45b9-9b47-1754c3257166	\N	\N	Avoid Dangerous Combinations	Research interactions before combining substances. Many combinations can be unexpectedly dangerous.	\N	\N	\N	2025-01-08 03:46:08.024795+00	Combinations	Critical	\N	\N	\N	\N	\N	unverified	\N	\N	2025-01-08 03:46:08.024795+00
\.


--
-- Data for Name: social_influence_networks; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_influence_networks (id, compound_id, platform, platform_subdivision, analysis_date, top_influencers, influence_connections, community_clusters, information_flow_patterns, key_opinion_leaders, emerging_voices, platform_specific_metrics, created_at) FROM stdin;
\.


--
-- Data for Name: social_platform_stats; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_platform_stats (id, platform, platform_subdivision, compound_id, total_posts, total_comments, unique_authors, scientific_post_ratio, experience_report_ratio, harm_reduction_ratio, top_compounds, topic_distribution, sentiment_distribution, quality_metrics, engagement_metrics, last_updated_at) FROM stdin;
\.


--
-- Data for Name: social_posts; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_posts (id, compound_id, platform, platform_subdivision, external_id, url, title, content, author_id, author_username, post_type, created_at, engagement_metrics, classification_data, sentiment_data, metadata, is_scientific, is_experience_report, is_harm_reduction, platform_specific_data, tags, research_citations, view_count, created_at_internal) FROM stdin;
\.


--
-- Data for Name: social_protocols; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_protocols (id, protocol_name, creator, description, target_outcome, compounds, protocol_steps, duration, frequency, monitoring_parameters, success_metrics, warnings, contraindications, reference_citations, review_score, platform, platform_subdivision, created_at) FROM stdin;
\.


--
-- Data for Name: social_research_reviews; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_research_reviews (id, compound_id, title, content, author, publication_date, study_type, methodology, findings, limitations, research_quality_score, citations, peer_review_notes, platform, platform_subdivision, created_at) FROM stdin;
\.


--
-- Data for Name: social_safety_incidents; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_safety_incidents (id, compound_id, incident_type, severity, description, reported_effects, reported_causes, platforms, platform_subdivisions, content_urls, verification_status, verification_notes, response_actions, resolution_status, resolution_notes, created_at) FROM stdin;
\.


--
-- Data for Name: social_scientific_content; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_scientific_content (id, post_id, compound_id, title, content, author, platform, platform_subdivision, publication_date, content_type, research_topics, methodology, findings, limitations, citations, peer_review_notes, quality_score, created_at) FROM stdin;
\.


--
-- Data for Name: social_stacks; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_stacks (id, stack_name, creator, description, purpose, compounds, dosages, timing_schedule, duration, reported_effects, side_effects, interactions, warnings, rating, review_count, platform, platform_subdivision, created_at) FROM stdin;
\.


--
-- Data for Name: social_substance_data; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.social_substance_data (id, compound_id, platform, platform_subdivision, external_id, common_names, chemical_class, psychoactive_class, summary, tolerance_data, roa_data, effect_data, onset, duration, after_effects, detection_time, test_kits, experiences, effects, aliases, avoid, warning_message, dangerous_interactions, uncertain_interactions, unsafe_interactions, risk_potential, toxicity_data, addiction_potential, cross_tolerance, chemistry_data, dosage_data, duration_data, health_effects, risk_factors, contraindications, interactions, legal_status, research_status, history, traditional_use, harm_reduction_notes, platform_specific_data, last_updated_at, created_at) FROM stdin;
\.


--
-- Data for Name: storage_requirements; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.storage_requirements (id, compound_id, storage_temp_min, storage_temp_max, temp_unit, humidity_requirements, light_sensitivity, air_sensitivity, storage_conditions, container_type, incompatible_materials, segregation_requirements, ventilation_requirements, static_protection, max_storage_time, storage_precautions, handling_precautions, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: structure_similarity; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.structure_similarity (id, compound_id_1, compound_id_2, similarity_metric, similarity_score, comparison_method, fingerprint_type, calculation_parameters, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: subjective_effect_categories; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.subjective_effect_categories (id, name, description, domain, level, parent_category_id, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: subjective_effects; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.subjective_effects (id, category_id, name, description, onset_characteristics, duration_characteristics, intensity_characteristics, common_variations, contributing_factors, risk_factors, management_strategies, research_status, notes, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: sync_status; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.sync_status (id, device_id, data_type, last_sync_at, sync_version, status, conflict_resolution, retry_count, next_retry_at, error_details, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: systems_biology_analysis; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.systems_biology_analysis (id, compound_id, analysis_level, network_effects, cellular_responses, metabolic_impact, signaling_cascades, regulatory_networks, adaptation_mechanisms, system_robustness, emergent_properties, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: therapeutic_class_categories; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.therapeutic_class_categories (id, name, description, level, parent_category_id, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: therapeutic_classes; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.therapeutic_classes (id, category_id, name, abbreviation, description, mechanism_of_action, primary_targets, therapeutic_uses, contraindications, typical_dosing, side_effects, drug_interactions, regulatory_status, clinical_guidelines, research_status, notes, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: toxicity_assays; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.toxicity_assays (id, compound_id, assay_type, cell_line, organism, endpoint, concentration, concentration_unit, exposure_time, result_value, result_unit, result_type, confidence_score, protocol_details, reference_doi, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: toxicity_endpoint_categories; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.toxicity_endpoint_categories (id, name, description, measurement_type, units, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: toxicity_endpoints; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.toxicity_endpoints (id, category_id, name, description, standard_unit, conversion_factors, detection_methods, validation_criteria, reference_ranges, severity_thresholds, regulatory_limits, notes, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: toxicity_mechanism_categories; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.toxicity_mechanism_categories (id, name, description, level, mechanism_type, parent_category_id, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: toxicity_mechanisms; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.toxicity_mechanisms (id, category_id, name, description, molecular_targets, cellular_effects, tissue_effects, systemic_effects, biomarkers, detection_methods, time_course, dose_response_characteristics, reversibility, risk_factors, preventive_measures, treatment_approaches, research_status, evidence_level, notes, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: training_datasets; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.training_datasets (id, dataset_name, description, source, version, compound_count, feature_count, feature_names, feature_types, target_names, target_types, preprocessing_steps, split_strategy, validation_method, data_statistics, quality_metrics, metadata, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: training_history; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.training_history (id, model_id, version_id, training_run_id, start_time, end_time, parameters, metrics, loss_history, validation_history, hardware_metrics, status, error_logs, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: trend_analysis; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.trend_analysis (id, analysis_type, target_type, target_id, time_period, metrics, insights, recommendations, priority, status, assigned_to, resolution_notes, next_review_date, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: ui_settings; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.ui_settings (id, user_id, theme, layout_preferences, display_options, notification_settings, accessibility_settings, custom_views, dashboard_config, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: usage_statistics; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.usage_statistics (id, feature, action, user_agent, device_type, session_id, user_id, duration_ms, success, error_details, performance_metrics, user_feedback, metadata, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: user_preferences; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.user_preferences (id, user_id, device_id, theme, layout, notifications, accessibility, data_preferences, sync_settings, privacy_settings, feature_flags, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: wavefunction_analysis; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.wavefunction_analysis (id, compound_id, box_probabilities, scaling_exponents, fractal_dimensions, correlation_dimension, localization_length, participation_ratio, analysis_parameters, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: web_components; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.web_components (id, name, description, component_type, configuration, dependencies, styling, client_scripts, server_scripts, version, is_active, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: web_data_sources; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.web_data_sources (id, category_id, name, base_url, description, api_endpoint, access_method, authentication_type, rate_limits, data_format, update_frequency, last_validated, validation_status, data_quality_metrics, coverage_areas, known_limitations, usage_requirements, citation_format, notes, active, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: web_hook_deliveries; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.web_hook_deliveries (id, hook_id, event_type, payload, response_status, response_body, delivery_status, attempt_count, next_retry_at, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: web_hooks; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.web_hooks (id, name, url, event_types, headers, is_active, secret_key, retry_config, timeout, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: web_source_categories; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.web_source_categories (id, name, description, source_type, reliability_rating, validation_requirements, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: web_templates; Type: TABLE DATA; Schema: public; Owner: armand
--

COPY public.web_templates (id, name, description, template_type, content, parameters, styling, scripts, version, is_active, created_at, updated_at) FROM stdin;
\.


--
-- Name: phase_transition_types_id_seq; Type: SEQUENCE SET; Schema: public; Owner: armand
--

SELECT pg_catalog.setval('public.phase_transition_types_id_seq', 4, true);


--
-- Name: quantum_basis_sets_id_seq; Type: SEQUENCE SET; Schema: public; Owner: armand
--

SELECT pg_catalog.setval('public.quantum_basis_sets_id_seq', 4, true);


--
-- Name: quantum_functionals_id_seq; Type: SEQUENCE SET; Schema: public; Owner: armand
--

SELECT pg_catalog.setval('public.quantum_functionals_id_seq', 4, true);


--
-- Name: quantum_observables_ref_id_seq; Type: SEQUENCE SET; Schema: public; Owner: armand
--

SELECT pg_catalog.setval('public.quantum_observables_ref_id_seq', 4, true);


--
-- Name: activity_cliffs activity_cliffs_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.activity_cliffs
    ADD CONSTRAINT activity_cliffs_pkey PRIMARY KEY (id);


--
-- Name: activity_correlations activity_correlations_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.activity_correlations
    ADD CONSTRAINT activity_correlations_pkey PRIMARY KEY (id);


--
-- Name: alert_triggers alert_triggers_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.alert_triggers
    ADD CONSTRAINT alert_triggers_pkey PRIMARY KEY (id);


--
-- Name: analysis_parameters analysis_parameters_parameter_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.analysis_parameters
    ADD CONSTRAINT analysis_parameters_parameter_name_key UNIQUE (parameter_name);


--
-- Name: analysis_parameters analysis_parameters_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.analysis_parameters
    ADD CONSTRAINT analysis_parameters_pkey PRIMARY KEY (id);


--
-- Name: api_endpoints api_endpoints_path_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.api_endpoints
    ADD CONSTRAINT api_endpoints_path_key UNIQUE (path);


--
-- Name: api_endpoints api_endpoints_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.api_endpoints
    ADD CONSTRAINT api_endpoints_pkey PRIMARY KEY (id);


--
-- Name: api_keys api_keys_key_hash_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.api_keys
    ADD CONSTRAINT api_keys_key_hash_key UNIQUE (key_hash);


--
-- Name: api_keys api_keys_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.api_keys
    ADD CONSTRAINT api_keys_pkey PRIMARY KEY (id);


--
-- Name: audit_log audit_log_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.audit_log
    ADD CONSTRAINT audit_log_pkey PRIMARY KEY (id);


--
-- Name: binding_assay_protocols binding_assay_protocols_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_assay_protocols
    ADD CONSTRAINT binding_assay_protocols_pkey PRIMARY KEY (id);


--
-- Name: binding_assay_types binding_assay_types_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_assay_types
    ADD CONSTRAINT binding_assay_types_pkey PRIMARY KEY (id);


--
-- Name: binding_data binding_data_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_data
    ADD CONSTRAINT binding_data_pkey PRIMARY KEY (id);


--
-- Name: binding_data_quality binding_data_quality_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_data_quality
    ADD CONSTRAINT binding_data_quality_pkey PRIMARY KEY (id);


--
-- Name: binding_kinetics binding_kinetics_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_kinetics
    ADD CONSTRAINT binding_kinetics_pkey PRIMARY KEY (id);


--
-- Name: binding_sar binding_sar_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_sar
    ADD CONSTRAINT binding_sar_pkey PRIMARY KEY (id);


--
-- Name: binding_site_mapping binding_site_mapping_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_site_mapping
    ADD CONSTRAINT binding_site_mapping_pkey PRIMARY KEY (id);


--
-- Name: carcinogenicity_data carcinogenicity_data_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.carcinogenicity_data
    ADD CONSTRAINT carcinogenicity_data_pkey PRIMARY KEY (id);


--
-- Name: cardiotoxicity_data cardiotoxicity_data_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.cardiotoxicity_data
    ADD CONSTRAINT cardiotoxicity_data_pkey PRIMARY KEY (id);


--
-- Name: compounds compounds_cas_number_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.compounds
    ADD CONSTRAINT compounds_cas_number_key UNIQUE (cas_number);


--
-- Name: compounds compounds_inchi_key_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.compounds
    ADD CONSTRAINT compounds_inchi_key_key UNIQUE (inchi_key);


--
-- Name: compounds compounds_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.compounds
    ADD CONSTRAINT compounds_pkey PRIMARY KEY (id);


--
-- Name: cytotoxicity_data cytotoxicity_data_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.cytotoxicity_data
    ADD CONSTRAINT cytotoxicity_data_pkey PRIMARY KEY (id);


--
-- Name: descriptors_2d descriptors_2d_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.descriptors_2d
    ADD CONSTRAINT descriptors_2d_pkey PRIMARY KEY (id);


--
-- Name: descriptors_3d descriptors_3d_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.descriptors_3d
    ADD CONSTRAINT descriptors_3d_pkey PRIMARY KEY (id);


--
-- Name: developmental_toxicity developmental_toxicity_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.developmental_toxicity
    ADD CONSTRAINT developmental_toxicity_pkey PRIMARY KEY (id);


--
-- Name: device_settings device_settings_device_id_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.device_settings
    ADD CONSTRAINT device_settings_device_id_key UNIQUE (device_id);


--
-- Name: device_settings device_settings_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.device_settings
    ADD CONSTRAINT device_settings_pkey PRIMARY KEY (id);


--
-- Name: effect_relationships effect_relationships_effect_id_related_effect_id_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.effect_relationships
    ADD CONSTRAINT effect_relationships_effect_id_related_effect_id_key UNIQUE (effect_id, related_effect_id);


--
-- Name: effect_relationships effect_relationships_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.effect_relationships
    ADD CONSTRAINT effect_relationships_pkey PRIMARY KEY (id);


--
-- Name: electronic_structure electronic_structure_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.electronic_structure
    ADD CONSTRAINT electronic_structure_pkey PRIMARY KEY (id);


--
-- Name: emergency_procedures emergency_procedures_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.emergency_procedures
    ADD CONSTRAINT emergency_procedures_pkey PRIMARY KEY (id);


--
-- Name: energy_level_statistics energy_level_statistics_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.energy_level_statistics
    ADD CONSTRAINT energy_level_statistics_pkey PRIMARY KEY (id);


--
-- Name: experience_categories experience_categories_category_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.experience_categories
    ADD CONSTRAINT experience_categories_category_name_key UNIQUE (category_name);


--
-- Name: experience_categories experience_categories_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.experience_categories
    ADD CONSTRAINT experience_categories_pkey PRIMARY KEY (id);


--
-- Name: external_services external_services_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.external_services
    ADD CONSTRAINT external_services_name_key UNIQUE (name);


--
-- Name: external_services external_services_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.external_services
    ADD CONSTRAINT external_services_pkey PRIMARY KEY (id);


--
-- Name: feature_definitions feature_definitions_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.feature_definitions
    ADD CONSTRAINT feature_definitions_name_key UNIQUE (name);


--
-- Name: feature_definitions feature_definitions_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.feature_definitions
    ADD CONSTRAINT feature_definitions_pkey PRIMARY KEY (id);


--
-- Name: feature_pipelines feature_pipelines_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.feature_pipelines
    ADD CONSTRAINT feature_pipelines_name_key UNIQUE (name);


--
-- Name: feature_pipelines feature_pipelines_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.feature_pipelines
    ADD CONSTRAINT feature_pipelines_pkey PRIMARY KEY (id);


--
-- Name: feature_selection feature_selection_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.feature_selection
    ADD CONSTRAINT feature_selection_pkey PRIMARY KEY (id);


--
-- Name: feature_values feature_values_compound_id_feature_id_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.feature_values
    ADD CONSTRAINT feature_values_compound_id_feature_id_key UNIQUE (compound_id, feature_id);


--
-- Name: feature_values feature_values_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.feature_values
    ADD CONSTRAINT feature_values_pkey PRIMARY KEY (id);


--
-- Name: genes genes_ensembl_id_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.genes
    ADD CONSTRAINT genes_ensembl_id_key UNIQUE (ensembl_id);


--
-- Name: genes genes_entrez_id_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.genes
    ADD CONSTRAINT genes_entrez_id_key UNIQUE (entrez_id);


--
-- Name: genes genes_hgnc_id_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.genes
    ADD CONSTRAINT genes_hgnc_id_key UNIQUE (hgnc_id);


--
-- Name: genes genes_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.genes
    ADD CONSTRAINT genes_pkey PRIMARY KEY (id);


--
-- Name: genes genes_symbol_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.genes
    ADD CONSTRAINT genes_symbol_key UNIQUE (symbol);


--
-- Name: genotoxicity_data genotoxicity_data_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.genotoxicity_data
    ADD CONSTRAINT genotoxicity_data_pkey PRIMARY KEY (id);


--
-- Name: hazard_classifications hazard_classifications_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.hazard_classifications
    ADD CONSTRAINT hazard_classifications_pkey PRIMARY KEY (id);


--
-- Name: hepatotoxicity_data hepatotoxicity_data_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.hepatotoxicity_data
    ADD CONSTRAINT hepatotoxicity_data_pkey PRIMARY KEY (id);


--
-- Name: immunotoxicity_data immunotoxicity_data_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.immunotoxicity_data
    ADD CONSTRAINT immunotoxicity_data_pkey PRIMARY KEY (id);


--
-- Name: literature_findings literature_findings_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.literature_findings
    ADD CONSTRAINT literature_findings_pkey PRIMARY KEY (id);


--
-- Name: meta_analyses meta_analyses_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.meta_analyses
    ADD CONSTRAINT meta_analyses_pkey PRIMARY KEY (id);


--
-- Name: metabolic_toxicity metabolic_toxicity_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.metabolic_toxicity
    ADD CONSTRAINT metabolic_toxicity_pkey PRIMARY KEY (id);


--
-- Name: ml_models ml_models_name_version_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.ml_models
    ADD CONSTRAINT ml_models_name_version_key UNIQUE (name, version);


--
-- Name: ml_models ml_models_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.ml_models
    ADD CONSTRAINT ml_models_pkey PRIMARY KEY (id);


--
-- Name: model_deployments model_deployments_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_deployments
    ADD CONSTRAINT model_deployments_pkey PRIMARY KEY (id);


--
-- Name: model_metrics model_metrics_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_metrics
    ADD CONSTRAINT model_metrics_pkey PRIMARY KEY (id);


--
-- Name: model_monitoring model_monitoring_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_monitoring
    ADD CONSTRAINT model_monitoring_pkey PRIMARY KEY (id);


--
-- Name: model_predictions model_predictions_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_predictions
    ADD CONSTRAINT model_predictions_pkey PRIMARY KEY (id);


--
-- Name: model_validation_results model_validation_results_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_validation_results
    ADD CONSTRAINT model_validation_results_pkey PRIMARY KEY (id);


--
-- Name: model_versions model_versions_model_id_version_number_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_versions
    ADD CONSTRAINT model_versions_model_id_version_number_key UNIQUE (model_id, version_number);


--
-- Name: model_versions model_versions_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_versions
    ADD CONSTRAINT model_versions_pkey PRIMARY KEY (id);


--
-- Name: molecular_fingerprints molecular_fingerprints_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.molecular_fingerprints
    ADD CONSTRAINT molecular_fingerprints_pkey PRIMARY KEY (id);


--
-- Name: monitoring_parameters monitoring_parameters_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.monitoring_parameters
    ADD CONSTRAINT monitoring_parameters_name_key UNIQUE (name);


--
-- Name: monitoring_parameters monitoring_parameters_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.monitoring_parameters
    ADD CONSTRAINT monitoring_parameters_pkey PRIMARY KEY (id);


--
-- Name: neurotoxicity_data neurotoxicity_data_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.neurotoxicity_data
    ADD CONSTRAINT neurotoxicity_data_pkey PRIMARY KEY (id);


--
-- Name: offline_data offline_data_device_id_data_type_data_id_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.offline_data
    ADD CONSTRAINT offline_data_device_id_data_type_data_id_key UNIQUE (device_id, data_type, data_id);


--
-- Name: offline_data offline_data_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.offline_data
    ADD CONSTRAINT offline_data_pkey PRIMARY KEY (id);


--
-- Name: organ_systems organ_systems_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.organ_systems
    ADD CONSTRAINT organ_systems_name_key UNIQUE (name);


--
-- Name: organ_systems organ_systems_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.organ_systems
    ADD CONSTRAINT organ_systems_pkey PRIMARY KEY (id);


--
-- Name: organ_toxicity_patterns organ_toxicity_patterns_organ_system_id_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.organ_toxicity_patterns
    ADD CONSTRAINT organ_toxicity_patterns_organ_system_id_name_key UNIQUE (organ_system_id, name);


--
-- Name: organ_toxicity_patterns organ_toxicity_patterns_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.organ_toxicity_patterns
    ADD CONSTRAINT organ_toxicity_patterns_pkey PRIMARY KEY (id);


--
-- Name: pathway_analysis pathway_analysis_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.pathway_analysis
    ADD CONSTRAINT pathway_analysis_pkey PRIMARY KEY (id);


--
-- Name: performance_metrics performance_metrics_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.performance_metrics
    ADD CONSTRAINT performance_metrics_pkey PRIMARY KEY (id);


--
-- Name: pharmacological_classes pharmacological_classes_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.pharmacological_classes
    ADD CONSTRAINT pharmacological_classes_name_key UNIQUE (name);


--
-- Name: pharmacological_classes pharmacological_classes_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.pharmacological_classes
    ADD CONSTRAINT pharmacological_classes_pkey PRIMARY KEY (id);


--
-- Name: pharmacophore_features pharmacophore_features_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.pharmacophore_features
    ADD CONSTRAINT pharmacophore_features_pkey PRIMARY KEY (id);


--
-- Name: phase_transition_types phase_transition_types_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.phase_transition_types
    ADD CONSTRAINT phase_transition_types_name_key UNIQUE (name);


--
-- Name: phase_transition_types phase_transition_types_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.phase_transition_types
    ADD CONSTRAINT phase_transition_types_pkey PRIMARY KEY (id);


--
-- Name: phase_transitions phase_transitions_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.phase_transitions
    ADD CONSTRAINT phase_transitions_pkey PRIMARY KEY (id);


--
-- Name: ppe_requirements ppe_requirements_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.ppe_requirements
    ADD CONSTRAINT ppe_requirements_pkey PRIMARY KEY (id);


--
-- Name: proteins proteins_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.proteins
    ADD CONSTRAINT proteins_pkey PRIMARY KEY (id);


--
-- Name: proteins proteins_uniprot_id_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.proteins
    ADD CONSTRAINT proteins_uniprot_id_key UNIQUE (uniprot_id);


--
-- Name: quantum_basis_sets quantum_basis_sets_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_basis_sets
    ADD CONSTRAINT quantum_basis_sets_name_key UNIQUE (name);


--
-- Name: quantum_basis_sets quantum_basis_sets_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_basis_sets
    ADD CONSTRAINT quantum_basis_sets_pkey PRIMARY KEY (id);


--
-- Name: quantum_calculations quantum_calculations_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_calculations
    ADD CONSTRAINT quantum_calculations_pkey PRIMARY KEY (id);


--
-- Name: quantum_critical_params quantum_critical_params_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_critical_params
    ADD CONSTRAINT quantum_critical_params_pkey PRIMARY KEY (id);


--
-- Name: quantum_dynamics quantum_dynamics_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_dynamics
    ADD CONSTRAINT quantum_dynamics_pkey PRIMARY KEY (id);


--
-- Name: quantum_functionals quantum_functionals_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_functionals
    ADD CONSTRAINT quantum_functionals_name_key UNIQUE (name);


--
-- Name: quantum_functionals quantum_functionals_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_functionals
    ADD CONSTRAINT quantum_functionals_pkey PRIMARY KEY (id);


--
-- Name: quantum_hamiltonians quantum_hamiltonians_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_hamiltonians
    ADD CONSTRAINT quantum_hamiltonians_pkey PRIMARY KEY (id);


--
-- Name: quantum_observables quantum_observables_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_observables
    ADD CONSTRAINT quantum_observables_pkey PRIMARY KEY (id);


--
-- Name: quantum_observables_ref quantum_observables_ref_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_observables_ref
    ADD CONSTRAINT quantum_observables_ref_name_key UNIQUE (name);


--
-- Name: quantum_observables_ref quantum_observables_ref_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_observables_ref
    ADD CONSTRAINT quantum_observables_ref_pkey PRIMARY KEY (id);


--
-- Name: quantum_parameters quantum_parameters_parameter_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_parameters
    ADD CONSTRAINT quantum_parameters_parameter_name_key UNIQUE (parameter_name);


--
-- Name: quantum_parameters quantum_parameters_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_parameters
    ADD CONSTRAINT quantum_parameters_pkey PRIMARY KEY (id);


--
-- Name: quantum_properties quantum_properties_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_properties
    ADD CONSTRAINT quantum_properties_pkey PRIMARY KEY (id);


--
-- Name: quantum_research_findings quantum_research_findings_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_research_findings
    ADD CONSTRAINT quantum_research_findings_pkey PRIMARY KEY (id);


--
-- Name: quantum_research_projects quantum_research_projects_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_research_projects
    ADD CONSTRAINT quantum_research_projects_pkey PRIMARY KEY (id);


--
-- Name: quantum_structure_correlations quantum_structure_correlations_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_structure_correlations
    ADD CONSTRAINT quantum_structure_correlations_pkey PRIMARY KEY (id);


--
-- Name: rate_limits rate_limits_api_key_id_endpoint_id_window_start_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.rate_limits
    ADD CONSTRAINT rate_limits_api_key_id_endpoint_id_window_start_key UNIQUE (api_key_id, endpoint_id, window_start);


--
-- Name: rate_limits rate_limits_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.rate_limits
    ADD CONSTRAINT rate_limits_pkey PRIMARY KEY (id);


--
-- Name: receptor_families receptor_families_category_id_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.receptor_families
    ADD CONSTRAINT receptor_families_category_id_name_key UNIQUE (category_id, name);


--
-- Name: receptor_families receptor_families_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.receptor_families
    ADD CONSTRAINT receptor_families_pkey PRIMARY KEY (id);


--
-- Name: receptor_family_categories receptor_family_categories_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.receptor_family_categories
    ADD CONSTRAINT receptor_family_categories_name_key UNIQUE (name);


--
-- Name: receptor_family_categories receptor_family_categories_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.receptor_family_categories
    ADD CONSTRAINT receptor_family_categories_pkey PRIMARY KEY (id);


--
-- Name: receptor_subtypes receptor_subtypes_family_id_subtype_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.receptor_subtypes
    ADD CONSTRAINT receptor_subtypes_family_id_subtype_name_key UNIQUE (family_id, subtype_name);


--
-- Name: receptor_subtypes receptor_subtypes_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.receptor_subtypes
    ADD CONSTRAINT receptor_subtypes_pkey PRIMARY KEY (id);


--
-- Name: research_findings research_findings_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.research_findings
    ADD CONSTRAINT research_findings_pkey PRIMARY KEY (id);


--
-- Name: safety_documents safety_documents_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.safety_documents
    ADD CONSTRAINT safety_documents_pkey PRIMARY KEY (id);


--
-- Name: safety_thresholds safety_thresholds_parameter_name_severity_level_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.safety_thresholds
    ADD CONSTRAINT safety_thresholds_parameter_name_severity_level_key UNIQUE (parameter_name, severity_level);


--
-- Name: safety_thresholds safety_thresholds_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.safety_thresholds
    ADD CONSTRAINT safety_thresholds_pkey PRIMARY KEY (id);


--
-- Name: sar_analysis sar_analysis_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.sar_analysis
    ADD CONSTRAINT sar_analysis_pkey PRIMARY KEY (id);


--
-- Name: sar_patterns sar_patterns_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.sar_patterns
    ADD CONSTRAINT sar_patterns_pkey PRIMARY KEY (id);


--
-- Name: scaffold_analysis scaffold_analysis_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.scaffold_analysis
    ADD CONSTRAINT scaffold_analysis_pkey PRIMARY KEY (id);


--
-- Name: scaling_analysis scaling_analysis_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.scaling_analysis
    ADD CONSTRAINT scaling_analysis_pkey PRIMARY KEY (id);


--
-- Name: schema_versions schema_versions_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.schema_versions
    ADD CONSTRAINT schema_versions_pkey PRIMARY KEY (id);


--
-- Name: scientific_papers scientific_papers_doi_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.scientific_papers
    ADD CONSTRAINT scientific_papers_doi_key UNIQUE (doi);


--
-- Name: scientific_papers scientific_papers_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.scientific_papers
    ADD CONSTRAINT scientific_papers_pkey PRIMARY KEY (id);


--
-- Name: service_status service_status_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.service_status
    ADD CONSTRAINT service_status_pkey PRIMARY KEY (id);


--
-- Name: social_alert_rules social_alert_rules_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_alert_rules
    ADD CONSTRAINT social_alert_rules_pkey PRIMARY KEY (id);


--
-- Name: social_comments social_comments_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_comments
    ADD CONSTRAINT social_comments_pkey PRIMARY KEY (id);


--
-- Name: social_comments social_comments_platform_external_id_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_comments
    ADD CONSTRAINT social_comments_platform_external_id_key UNIQUE (platform, external_id);


--
-- Name: social_compound_combinations social_compound_combinations_compound_id_1_compound_id_2_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_compound_combinations
    ADD CONSTRAINT social_compound_combinations_compound_id_1_compound_id_2_key UNIQUE (compound_id_1, compound_id_2);


--
-- Name: social_compound_combinations social_compound_combinations_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_compound_combinations
    ADD CONSTRAINT social_compound_combinations_pkey PRIMARY KEY (id);


--
-- Name: social_compound_mentions social_compound_mentions_compound_id_platform_platform_subd_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_compound_mentions
    ADD CONSTRAINT social_compound_mentions_compound_id_platform_platform_subd_key UNIQUE (compound_id, platform, platform_subdivision, analysis_date);


--
-- Name: social_compound_mentions social_compound_mentions_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_compound_mentions
    ADD CONSTRAINT social_compound_mentions_pkey PRIMARY KEY (id);


--
-- Name: social_compound_trends social_compound_trends_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_compound_trends
    ADD CONSTRAINT social_compound_trends_pkey PRIMARY KEY (id);


--
-- Name: social_content_quality social_content_quality_compound_id_platform_platform_subdiv_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_content_quality
    ADD CONSTRAINT social_content_quality_compound_id_platform_platform_subdiv_key UNIQUE (compound_id, platform, platform_subdivision, analysis_date);


--
-- Name: social_content_quality social_content_quality_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_content_quality
    ADD CONSTRAINT social_content_quality_pkey PRIMARY KEY (id);


--
-- Name: social_dosage_data social_dosage_data_compound_id_platform_route_of_administra_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_dosage_data
    ADD CONSTRAINT social_dosage_data_compound_id_platform_route_of_administra_key UNIQUE (compound_id, platform, route_of_administration);


--
-- Name: social_dosage_data social_dosage_data_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_dosage_data
    ADD CONSTRAINT social_dosage_data_pkey PRIMARY KEY (id);


--
-- Name: social_effect_reports social_effect_reports_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_effect_reports
    ADD CONSTRAINT social_effect_reports_pkey PRIMARY KEY (id);


--
-- Name: social_effects social_effects_effect_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_effects
    ADD CONSTRAINT social_effects_effect_name_key UNIQUE (effect_name);


--
-- Name: social_effects social_effects_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_effects
    ADD CONSTRAINT social_effects_pkey PRIMARY KEY (id);


--
-- Name: social_experience_reports social_experience_reports_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_experience_reports
    ADD CONSTRAINT social_experience_reports_pkey PRIMARY KEY (id);


--
-- Name: social_harm_reduction social_harm_reduction_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_harm_reduction
    ADD CONSTRAINT social_harm_reduction_pkey PRIMARY KEY (id);


--
-- Name: social_influence_networks social_influence_networks_compound_id_platform_platform_sub_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_influence_networks
    ADD CONSTRAINT social_influence_networks_compound_id_platform_platform_sub_key UNIQUE (compound_id, platform, platform_subdivision, analysis_date);


--
-- Name: social_influence_networks social_influence_networks_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_influence_networks
    ADD CONSTRAINT social_influence_networks_pkey PRIMARY KEY (id);


--
-- Name: social_platform_stats social_platform_stats_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_platform_stats
    ADD CONSTRAINT social_platform_stats_pkey PRIMARY KEY (id);


--
-- Name: social_platform_stats social_platform_stats_platform_platform_subdivision_compoun_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_platform_stats
    ADD CONSTRAINT social_platform_stats_platform_platform_subdivision_compoun_key UNIQUE (platform, platform_subdivision, compound_id);


--
-- Name: social_posts social_posts_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_posts
    ADD CONSTRAINT social_posts_pkey PRIMARY KEY (id);


--
-- Name: social_posts social_posts_platform_external_id_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_posts
    ADD CONSTRAINT social_posts_platform_external_id_key UNIQUE (platform, external_id);


--
-- Name: social_protocols social_protocols_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_protocols
    ADD CONSTRAINT social_protocols_pkey PRIMARY KEY (id);


--
-- Name: social_research_reviews social_research_reviews_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_research_reviews
    ADD CONSTRAINT social_research_reviews_pkey PRIMARY KEY (id);


--
-- Name: social_safety_incidents social_safety_incidents_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_safety_incidents
    ADD CONSTRAINT social_safety_incidents_pkey PRIMARY KEY (id);


--
-- Name: social_scientific_content social_scientific_content_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_scientific_content
    ADD CONSTRAINT social_scientific_content_pkey PRIMARY KEY (id);


--
-- Name: social_stacks social_stacks_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_stacks
    ADD CONSTRAINT social_stacks_pkey PRIMARY KEY (id);


--
-- Name: social_substance_data social_substance_data_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_substance_data
    ADD CONSTRAINT social_substance_data_pkey PRIMARY KEY (id);


--
-- Name: social_substance_data social_substance_data_platform_external_id_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_substance_data
    ADD CONSTRAINT social_substance_data_platform_external_id_key UNIQUE (platform, external_id);


--
-- Name: storage_requirements storage_requirements_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.storage_requirements
    ADD CONSTRAINT storage_requirements_pkey PRIMARY KEY (id);


--
-- Name: structure_similarity structure_similarity_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.structure_similarity
    ADD CONSTRAINT structure_similarity_pkey PRIMARY KEY (id);


--
-- Name: subjective_effect_categories subjective_effect_categories_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.subjective_effect_categories
    ADD CONSTRAINT subjective_effect_categories_name_key UNIQUE (name);


--
-- Name: subjective_effect_categories subjective_effect_categories_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.subjective_effect_categories
    ADD CONSTRAINT subjective_effect_categories_pkey PRIMARY KEY (id);


--
-- Name: subjective_effects subjective_effects_category_id_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.subjective_effects
    ADD CONSTRAINT subjective_effects_category_id_name_key UNIQUE (category_id, name);


--
-- Name: subjective_effects subjective_effects_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.subjective_effects
    ADD CONSTRAINT subjective_effects_pkey PRIMARY KEY (id);


--
-- Name: sync_status sync_status_device_id_data_type_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.sync_status
    ADD CONSTRAINT sync_status_device_id_data_type_key UNIQUE (device_id, data_type);


--
-- Name: sync_status sync_status_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.sync_status
    ADD CONSTRAINT sync_status_pkey PRIMARY KEY (id);


--
-- Name: systems_biology_analysis systems_biology_analysis_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.systems_biology_analysis
    ADD CONSTRAINT systems_biology_analysis_pkey PRIMARY KEY (id);


--
-- Name: therapeutic_class_categories therapeutic_class_categories_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.therapeutic_class_categories
    ADD CONSTRAINT therapeutic_class_categories_name_key UNIQUE (name);


--
-- Name: therapeutic_class_categories therapeutic_class_categories_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.therapeutic_class_categories
    ADD CONSTRAINT therapeutic_class_categories_pkey PRIMARY KEY (id);


--
-- Name: therapeutic_classes therapeutic_classes_category_id_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.therapeutic_classes
    ADD CONSTRAINT therapeutic_classes_category_id_name_key UNIQUE (category_id, name);


--
-- Name: therapeutic_classes therapeutic_classes_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.therapeutic_classes
    ADD CONSTRAINT therapeutic_classes_pkey PRIMARY KEY (id);


--
-- Name: toxicity_assays toxicity_assays_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.toxicity_assays
    ADD CONSTRAINT toxicity_assays_pkey PRIMARY KEY (id);


--
-- Name: toxicity_endpoint_categories toxicity_endpoint_categories_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.toxicity_endpoint_categories
    ADD CONSTRAINT toxicity_endpoint_categories_name_key UNIQUE (name);


--
-- Name: toxicity_endpoint_categories toxicity_endpoint_categories_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.toxicity_endpoint_categories
    ADD CONSTRAINT toxicity_endpoint_categories_pkey PRIMARY KEY (id);


--
-- Name: toxicity_endpoints toxicity_endpoints_category_id_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.toxicity_endpoints
    ADD CONSTRAINT toxicity_endpoints_category_id_name_key UNIQUE (category_id, name);


--
-- Name: toxicity_endpoints toxicity_endpoints_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.toxicity_endpoints
    ADD CONSTRAINT toxicity_endpoints_pkey PRIMARY KEY (id);


--
-- Name: toxicity_mechanism_categories toxicity_mechanism_categories_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.toxicity_mechanism_categories
    ADD CONSTRAINT toxicity_mechanism_categories_name_key UNIQUE (name);


--
-- Name: toxicity_mechanism_categories toxicity_mechanism_categories_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.toxicity_mechanism_categories
    ADD CONSTRAINT toxicity_mechanism_categories_pkey PRIMARY KEY (id);


--
-- Name: toxicity_mechanisms toxicity_mechanisms_category_id_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.toxicity_mechanisms
    ADD CONSTRAINT toxicity_mechanisms_category_id_name_key UNIQUE (category_id, name);


--
-- Name: toxicity_mechanisms toxicity_mechanisms_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.toxicity_mechanisms
    ADD CONSTRAINT toxicity_mechanisms_pkey PRIMARY KEY (id);


--
-- Name: training_datasets training_datasets_dataset_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.training_datasets
    ADD CONSTRAINT training_datasets_dataset_name_key UNIQUE (dataset_name);


--
-- Name: training_datasets training_datasets_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.training_datasets
    ADD CONSTRAINT training_datasets_pkey PRIMARY KEY (id);


--
-- Name: training_history training_history_model_id_training_run_id_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.training_history
    ADD CONSTRAINT training_history_model_id_training_run_id_key UNIQUE (model_id, training_run_id);


--
-- Name: training_history training_history_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.training_history
    ADD CONSTRAINT training_history_pkey PRIMARY KEY (id);


--
-- Name: trend_analysis trend_analysis_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.trend_analysis
    ADD CONSTRAINT trend_analysis_pkey PRIMARY KEY (id);


--
-- Name: ui_settings ui_settings_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.ui_settings
    ADD CONSTRAINT ui_settings_pkey PRIMARY KEY (id);


--
-- Name: ui_settings ui_settings_user_id_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.ui_settings
    ADD CONSTRAINT ui_settings_user_id_key UNIQUE (user_id);


--
-- Name: usage_statistics usage_statistics_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.usage_statistics
    ADD CONSTRAINT usage_statistics_pkey PRIMARY KEY (id);


--
-- Name: user_preferences user_preferences_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.user_preferences
    ADD CONSTRAINT user_preferences_pkey PRIMARY KEY (id);


--
-- Name: user_preferences user_preferences_user_id_device_id_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.user_preferences
    ADD CONSTRAINT user_preferences_user_id_device_id_key UNIQUE (user_id, device_id);


--
-- Name: wavefunction_analysis wavefunction_analysis_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.wavefunction_analysis
    ADD CONSTRAINT wavefunction_analysis_pkey PRIMARY KEY (id);


--
-- Name: web_components web_components_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.web_components
    ADD CONSTRAINT web_components_name_key UNIQUE (name);


--
-- Name: web_components web_components_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.web_components
    ADD CONSTRAINT web_components_pkey PRIMARY KEY (id);


--
-- Name: web_data_sources web_data_sources_category_id_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.web_data_sources
    ADD CONSTRAINT web_data_sources_category_id_name_key UNIQUE (category_id, name);


--
-- Name: web_data_sources web_data_sources_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.web_data_sources
    ADD CONSTRAINT web_data_sources_pkey PRIMARY KEY (id);


--
-- Name: web_hook_deliveries web_hook_deliveries_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.web_hook_deliveries
    ADD CONSTRAINT web_hook_deliveries_pkey PRIMARY KEY (id);


--
-- Name: web_hooks web_hooks_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.web_hooks
    ADD CONSTRAINT web_hooks_pkey PRIMARY KEY (id);


--
-- Name: web_source_categories web_source_categories_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.web_source_categories
    ADD CONSTRAINT web_source_categories_name_key UNIQUE (name);


--
-- Name: web_source_categories web_source_categories_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.web_source_categories
    ADD CONSTRAINT web_source_categories_pkey PRIMARY KEY (id);


--
-- Name: web_templates web_templates_name_key; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.web_templates
    ADD CONSTRAINT web_templates_name_key UNIQUE (name);


--
-- Name: web_templates web_templates_pkey; Type: CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.web_templates
    ADD CONSTRAINT web_templates_pkey PRIMARY KEY (id);


--
-- Name: idx_activity_cliffs_compounds; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_activity_cliffs_compounds ON public.activity_cliffs USING gin (compound_pair);


--
-- Name: idx_activity_cliffs_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_activity_cliffs_type ON public.activity_cliffs USING btree (activity_type);


--
-- Name: idx_activity_correlations_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_activity_correlations_compound ON public.activity_correlations USING btree (compound_id);


--
-- Name: idx_activity_correlations_feature; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_activity_correlations_feature ON public.activity_correlations USING btree (feature_id);


--
-- Name: idx_activity_correlations_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_activity_correlations_type ON public.activity_correlations USING btree (activity_type);


--
-- Name: idx_alert_rules_active; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_alert_rules_active ON public.social_alert_rules USING btree (is_active);


--
-- Name: idx_alert_rules_compounds; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_alert_rules_compounds ON public.social_alert_rules USING gin (compounds);


--
-- Name: idx_alert_rules_platforms; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_alert_rules_platforms ON public.social_alert_rules USING gin (platforms);


--
-- Name: idx_alert_rules_severity; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_alert_rules_severity ON public.social_alert_rules USING btree (severity);


--
-- Name: idx_alert_rules_subdivisions; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_alert_rules_subdivisions ON public.social_alert_rules USING gin (platform_subdivisions);


--
-- Name: idx_alert_rules_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_alert_rules_type ON public.social_alert_rules USING btree (alert_type);


--
-- Name: idx_alert_triggers_active; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_alert_triggers_active ON public.alert_triggers USING btree (is_active);


--
-- Name: idx_alert_triggers_priority; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_alert_triggers_priority ON public.alert_triggers USING btree (priority);


--
-- Name: idx_alert_triggers_severity; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_alert_triggers_severity ON public.alert_triggers USING btree (severity);


--
-- Name: idx_alert_triggers_target; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_alert_triggers_target ON public.alert_triggers USING btree (target_type, target_id);


--
-- Name: idx_alert_triggers_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_alert_triggers_type ON public.alert_triggers USING btree (trigger_type);


--
-- Name: idx_analysis_params_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_analysis_params_name ON public.analysis_parameters USING btree (parameter_name);


--
-- Name: idx_analysis_params_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_analysis_params_type ON public.analysis_parameters USING btree (data_type);


--
-- Name: idx_api_keys_active; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_api_keys_active ON public.api_keys USING btree (is_active);


--
-- Name: idx_api_keys_expires; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_api_keys_expires ON public.api_keys USING btree (expires_at);


--
-- Name: idx_api_keys_hash; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_api_keys_hash ON public.api_keys USING btree (key_hash);


--
-- Name: idx_api_keys_user; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_api_keys_user ON public.api_keys USING btree (user_id);


--
-- Name: idx_audit_log_action; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_audit_log_action ON public.audit_log USING btree (action);


--
-- Name: idx_audit_log_changed; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_audit_log_changed ON public.audit_log USING btree (changed_at);


--
-- Name: idx_audit_log_record; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_audit_log_record ON public.audit_log USING btree (record_id);


--
-- Name: idx_audit_log_table; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_audit_log_table ON public.audit_log USING btree (table_name);


--
-- Name: idx_combinations_compounds; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_combinations_compounds ON public.social_compound_combinations USING btree (compound_id_1, compound_id_2);


--
-- Name: idx_combinations_platforms; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_combinations_platforms ON public.social_compound_combinations USING gin (platforms);


--
-- Name: idx_combinations_risk; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_combinations_risk ON public.social_compound_combinations USING btree (risk_level);


--
-- Name: idx_combinations_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_combinations_type ON public.social_compound_combinations USING btree (interaction_type);


--
-- Name: idx_community_content_search; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_community_content_search ON public.social_posts USING gin (to_tsvector('english'::regconfig, content));


--
-- Name: idx_community_spam; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_community_spam ON public.social_posts USING btree (((classification_data ->> 'spam_score'::text)));


--
-- Name: idx_community_toxicity; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_community_toxicity ON public.social_posts USING btree (((classification_data ->> 'toxicity_score'::text)));


--
-- Name: idx_components_active; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_components_active ON public.web_components USING btree (is_active);


--
-- Name: idx_components_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_components_name ON public.web_components USING btree (name);


--
-- Name: idx_components_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_components_type ON public.web_components USING btree (component_type);


--
-- Name: idx_compounds_cas; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_compounds_cas ON public.compounds USING btree (cas_number);


--
-- Name: idx_compounds_chembl; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_compounds_chembl ON public.compounds USING btree (chembl_id);


--
-- Name: idx_compounds_created; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_compounds_created ON public.compounds USING btree (created_at);


--
-- Name: idx_compounds_drugbank; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_compounds_drugbank ON public.compounds USING btree (drugbank_id);


--
-- Name: idx_compounds_inchi; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_compounds_inchi ON public.compounds USING btree (inchi);


--
-- Name: idx_compounds_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_compounds_name ON public.compounds USING btree (name);


--
-- Name: idx_compounds_pubchem; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_compounds_pubchem ON public.compounds USING btree (pubchem_cid);


--
-- Name: idx_compounds_smiles; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_compounds_smiles ON public.compounds USING btree (smiles);


--
-- Name: idx_compounds_updated; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_compounds_updated ON public.compounds USING btree (updated_at);


--
-- Name: idx_device_settings_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_device_settings_type ON public.device_settings USING btree (device_type);


--
-- Name: idx_device_settings_version; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_device_settings_version ON public.device_settings USING btree (app_version);


--
-- Name: idx_dosage_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_dosage_compound ON public.social_dosage_data USING btree (compound_id);


--
-- Name: idx_dosage_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_dosage_platform ON public.social_dosage_data USING btree (platform);


--
-- Name: idx_dosage_roa; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_dosage_roa ON public.social_dosage_data USING btree (route_of_administration);


--
-- Name: idx_dosage_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_dosage_subdivision ON public.social_dosage_data USING btree (platform_subdivision);


--
-- Name: idx_effect_categories_domain; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_effect_categories_domain ON public.subjective_effect_categories USING btree (domain);


--
-- Name: idx_effect_categories_level; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_effect_categories_level ON public.subjective_effect_categories USING btree (level);


--
-- Name: idx_effect_categories_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_effect_categories_name ON public.subjective_effect_categories USING btree (name);


--
-- Name: idx_effect_categories_parent; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_effect_categories_parent ON public.subjective_effect_categories USING btree (parent_category_id);


--
-- Name: idx_effect_relationships_effect; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_effect_relationships_effect ON public.effect_relationships USING btree (effect_id);


--
-- Name: idx_effect_relationships_related; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_effect_relationships_related ON public.effect_relationships USING btree (related_effect_id);


--
-- Name: idx_effect_relationships_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_effect_relationships_type ON public.effect_relationships USING btree (relationship_type);


--
-- Name: idx_effect_reports_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_effect_reports_compound ON public.social_effect_reports USING btree (compound_id);


--
-- Name: idx_effect_reports_effect; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_effect_reports_effect ON public.social_effect_reports USING btree (effect_id);


--
-- Name: idx_effect_reports_intensity; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_effect_reports_intensity ON public.social_effect_reports USING btree (intensity);


--
-- Name: idx_effect_reports_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_effect_reports_platform ON public.social_effect_reports USING btree (platform);


--
-- Name: idx_effect_reports_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_effect_reports_subdivision ON public.social_effect_reports USING btree (platform_subdivision);


--
-- Name: idx_effects_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_effects_name ON public.social_effects USING btree (effect_name);


--
-- Name: idx_effects_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_effects_type ON public.social_effects USING btree (type);


--
-- Name: idx_electronic_structure_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_electronic_structure_compound ON public.electronic_structure USING btree (compound_id);


--
-- Name: idx_electronic_structure_density; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_electronic_structure_density ON public.electronic_structure USING gin (density_matrix);


--
-- Name: idx_endpoints_deprecated; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_endpoints_deprecated ON public.api_endpoints USING btree (is_deprecated);


--
-- Name: idx_endpoints_method; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_endpoints_method ON public.api_endpoints USING btree (method);


--
-- Name: idx_endpoints_path; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_endpoints_path ON public.api_endpoints USING btree (path);


--
-- Name: idx_endpoints_version; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_endpoints_version ON public.api_endpoints USING btree (version);


--
-- Name: idx_experience_categories_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_experience_categories_name ON public.experience_categories USING btree (category_name);


--
-- Name: idx_experience_categories_parent; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_experience_categories_parent ON public.experience_categories USING btree (parent_category);


--
-- Name: idx_experience_reports_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_experience_reports_compound ON public.social_experience_reports USING btree (compound_id);


--
-- Name: idx_experience_reports_date; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_experience_reports_date ON public.social_experience_reports USING btree (experience_date);


--
-- Name: idx_experience_reports_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_experience_reports_platform ON public.social_experience_reports USING btree (platform);


--
-- Name: idx_experience_reports_post; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_experience_reports_post ON public.social_experience_reports USING btree (post_id);


--
-- Name: idx_experience_reports_protocol; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_experience_reports_protocol ON public.social_experience_reports USING btree (protocol_id);


--
-- Name: idx_experience_reports_stack; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_experience_reports_stack ON public.social_experience_reports USING btree (stack_id);


--
-- Name: idx_experience_reports_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_experience_reports_subdivision ON public.social_experience_reports USING btree (platform_subdivision);


--
-- Name: idx_feature_values_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_feature_values_compound ON public.feature_values USING btree (compound_id);


--
-- Name: idx_feature_values_date; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_feature_values_date ON public.feature_values USING btree (calculation_date);


--
-- Name: idx_feature_values_feature; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_feature_values_feature ON public.feature_values USING btree (feature_id);


--
-- Name: idx_features_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_features_name ON public.feature_definitions USING btree (name);


--
-- Name: idx_features_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_features_type ON public.feature_definitions USING btree (feature_type);


--
-- Name: idx_findings_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_findings_compound ON public.literature_findings USING btree (compound_id);


--
-- Name: idx_findings_paper; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_findings_paper ON public.literature_findings USING btree (paper_id);


--
-- Name: idx_findings_replication; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_findings_replication ON public.literature_findings USING btree (replication_status);


--
-- Name: idx_findings_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_findings_type ON public.literature_findings USING btree (finding_type);


--
-- Name: idx_genes_ensembl; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_genes_ensembl ON public.genes USING btree (ensembl_id);


--
-- Name: idx_genes_entrez; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_genes_entrez ON public.genes USING btree (entrez_id);


--
-- Name: idx_genes_hgnc; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_genes_hgnc ON public.genes USING btree (hgnc_id);


--
-- Name: idx_genes_organism; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_genes_organism ON public.genes USING btree (organism);


--
-- Name: idx_genes_symbol; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_genes_symbol ON public.genes USING btree (symbol);


--
-- Name: idx_harm_reduction_category; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_harm_reduction_category ON public.social_harm_reduction USING btree (category);


--
-- Name: idx_harm_reduction_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_harm_reduction_compound ON public.social_harm_reduction USING btree (compound_id);


--
-- Name: idx_harm_reduction_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_harm_reduction_platform ON public.social_harm_reduction USING btree (platform);


--
-- Name: idx_harm_reduction_post; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_harm_reduction_post ON public.social_harm_reduction USING btree (post_id);


--
-- Name: idx_harm_reduction_status; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_harm_reduction_status ON public.social_harm_reduction USING btree (verification_status);


--
-- Name: idx_harm_reduction_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_harm_reduction_subdivision ON public.social_harm_reduction USING btree (platform_subdivision);


--
-- Name: idx_hook_deliveries_event; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_hook_deliveries_event ON public.web_hook_deliveries USING btree (event_type);


--
-- Name: idx_hook_deliveries_hook; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_hook_deliveries_hook ON public.web_hook_deliveries USING btree (hook_id);


--
-- Name: idx_hook_deliveries_retry; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_hook_deliveries_retry ON public.web_hook_deliveries USING btree (next_retry_at);


--
-- Name: idx_hook_deliveries_status; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_hook_deliveries_status ON public.web_hook_deliveries USING btree (delivery_status);


--
-- Name: idx_hooks_active; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_hooks_active ON public.web_hooks USING btree (is_active);


--
-- Name: idx_hooks_events; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_hooks_events ON public.web_hooks USING gin (event_types);


--
-- Name: idx_hooks_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_hooks_name ON public.web_hooks USING btree (name);


--
-- Name: idx_influence_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_influence_compound ON public.social_influence_networks USING btree (compound_id);


--
-- Name: idx_influence_date; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_influence_date ON public.social_influence_networks USING btree (analysis_date);


--
-- Name: idx_influence_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_influence_platform ON public.social_influence_networks USING btree (platform);


--
-- Name: idx_influence_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_influence_subdivision ON public.social_influence_networks USING btree (platform_subdivision);


--
-- Name: idx_mentions_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_mentions_compound ON public.social_compound_mentions USING btree (compound_id);


--
-- Name: idx_mentions_date; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_mentions_date ON public.social_compound_mentions USING btree (analysis_date);


--
-- Name: idx_mentions_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_mentions_platform ON public.social_compound_mentions USING btree (platform);


--
-- Name: idx_mentions_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_mentions_subdivision ON public.social_compound_mentions USING btree (platform_subdivision);


--
-- Name: idx_meta_analyses_evidence; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_meta_analyses_evidence ON public.meta_analyses USING btree (evidence_strength);


--
-- Name: idx_meta_analyses_papers; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_meta_analyses_papers ON public.meta_analyses USING gin (included_papers);


--
-- Name: idx_meta_analyses_topic; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_meta_analyses_topic ON public.meta_analyses USING btree (topic);


--
-- Name: idx_metrics_date; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_metrics_date ON public.model_metrics USING btree (metric_date);


--
-- Name: idx_metrics_model; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_metrics_model ON public.model_metrics USING btree (model_id);


--
-- Name: idx_metrics_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_metrics_type ON public.model_metrics USING btree (metric_type);


--
-- Name: idx_metrics_version; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_metrics_version ON public.model_metrics USING btree (version_id);


--
-- Name: idx_model_deployments_env; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_model_deployments_env ON public.model_deployments USING btree (deployment_environment);


--
-- Name: idx_model_deployments_model; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_model_deployments_model ON public.model_deployments USING btree (model_id);


--
-- Name: idx_model_deployments_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_model_deployments_name ON public.model_deployments USING btree (deployment_name);


--
-- Name: idx_model_deployments_status; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_model_deployments_status ON public.model_deployments USING btree (status);


--
-- Name: idx_model_monitoring_date; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_model_monitoring_date ON public.model_monitoring USING btree (monitoring_date);


--
-- Name: idx_model_monitoring_metric; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_model_monitoring_metric ON public.model_monitoring USING btree (metric_name);


--
-- Name: idx_model_monitoring_model; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_model_monitoring_model ON public.model_monitoring USING btree (model_id);


--
-- Name: idx_model_monitoring_status; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_model_monitoring_status ON public.model_monitoring USING btree (alert_status);


--
-- Name: idx_models_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_models_name ON public.ml_models USING btree (name);


--
-- Name: idx_models_target; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_models_target ON public.ml_models USING btree (target_variable);


--
-- Name: idx_models_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_models_type ON public.ml_models USING btree (model_type);


--
-- Name: idx_monitoring_params_category; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_monitoring_params_category ON public.monitoring_parameters USING btree (category);


--
-- Name: idx_monitoring_params_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_monitoring_params_name ON public.monitoring_parameters USING btree (name);


--
-- Name: idx_offline_data_device; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_offline_data_device ON public.offline_data USING btree (device_id);


--
-- Name: idx_offline_data_expires; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_offline_data_expires ON public.offline_data USING btree (expires_at);


--
-- Name: idx_offline_data_priority; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_offline_data_priority ON public.offline_data USING btree (priority);


--
-- Name: idx_offline_data_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_offline_data_type ON public.offline_data USING btree (data_type);


--
-- Name: idx_organ_systems_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_organ_systems_name ON public.organ_systems USING btree (name);


--
-- Name: idx_organ_toxicity_patterns_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_organ_toxicity_patterns_name ON public.organ_toxicity_patterns USING btree (name);


--
-- Name: idx_organ_toxicity_patterns_organ; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_organ_toxicity_patterns_organ ON public.organ_toxicity_patterns USING btree (organ_system_id);


--
-- Name: idx_organ_toxicity_patterns_text_search; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_organ_toxicity_patterns_text_search ON public.organ_toxicity_patterns USING gin (text_search_vector);


--
-- Name: idx_papers_compounds; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_papers_compounds ON public.scientific_papers USING gin (compounds_studied);


--
-- Name: idx_papers_date; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_papers_date ON public.scientific_papers USING btree (publication_date);


--
-- Name: idx_papers_doi; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_papers_doi ON public.scientific_papers USING btree (doi);


--
-- Name: idx_papers_evidence; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_papers_evidence ON public.scientific_papers USING btree (evidence_level);


--
-- Name: idx_papers_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_papers_type ON public.scientific_papers USING btree (study_type);


--
-- Name: idx_papers_validation; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_papers_validation ON public.scientific_papers USING btree (validation_status);


--
-- Name: idx_pathway_analysis_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_pathway_analysis_compound ON public.pathway_analysis USING btree (compound_id);


--
-- Name: idx_pathway_analysis_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_pathway_analysis_name ON public.pathway_analysis USING btree (pathway_name);


--
-- Name: idx_pathway_analysis_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_pathway_analysis_type ON public.pathway_analysis USING btree (pathway_type);


--
-- Name: idx_performance_metrics_component; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_performance_metrics_component ON public.performance_metrics USING btree (component);


--
-- Name: idx_performance_metrics_measured; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_performance_metrics_measured ON public.performance_metrics USING btree (measured_at);


--
-- Name: idx_performance_metrics_status; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_performance_metrics_status ON public.performance_metrics USING btree (status);


--
-- Name: idx_performance_metrics_tags; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_performance_metrics_tags ON public.performance_metrics USING gin (tags);


--
-- Name: idx_performance_metrics_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_performance_metrics_type ON public.performance_metrics USING btree (metric_type);


--
-- Name: idx_pharmacological_classes_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_pharmacological_classes_name ON public.pharmacological_classes USING btree (name);


--
-- Name: idx_pharmacological_classes_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_pharmacological_classes_type ON public.pharmacological_classes USING btree (mechanism_type);


--
-- Name: idx_phase_transitions_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_phase_transitions_compound ON public.phase_transitions USING btree (compound_id);


--
-- Name: idx_phase_transitions_order; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_phase_transitions_order ON public.phase_transition_types USING btree (transition_order);


--
-- Name: idx_pipelines_input; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_pipelines_input ON public.feature_pipelines USING gin (input_features);


--
-- Name: idx_pipelines_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_pipelines_name ON public.feature_pipelines USING btree (name);


--
-- Name: idx_pipelines_output; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_pipelines_output ON public.feature_pipelines USING gin (output_features);


--
-- Name: idx_platform_stats_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_platform_stats_compound ON public.social_platform_stats USING btree (compound_id);


--
-- Name: idx_platform_stats_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_platform_stats_platform ON public.social_platform_stats USING btree (platform);


--
-- Name: idx_platform_stats_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_platform_stats_subdivision ON public.social_platform_stats USING btree (platform_subdivision);


--
-- Name: idx_platform_stats_updated; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_platform_stats_updated ON public.social_platform_stats USING btree (last_updated_at);


--
-- Name: idx_predictions_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_predictions_compound ON public.model_predictions USING btree (compound_id);


--
-- Name: idx_predictions_date; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_predictions_date ON public.model_predictions USING btree (prediction_date);


--
-- Name: idx_predictions_model; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_predictions_model ON public.model_predictions USING btree (model_id);


--
-- Name: idx_predictions_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_predictions_type ON public.model_predictions USING btree (prediction_type);


--
-- Name: idx_predictions_version; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_predictions_version ON public.model_predictions USING btree (version_id);


--
-- Name: idx_proteins_family; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_proteins_family ON public.proteins USING btree (protein_family);


--
-- Name: idx_proteins_gene; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_proteins_gene ON public.proteins USING btree (gene_id);


--
-- Name: idx_proteins_organism; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_proteins_organism ON public.proteins USING btree (organism);


--
-- Name: idx_proteins_symbol; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_proteins_symbol ON public.proteins USING btree (symbol);


--
-- Name: idx_proteins_uniprot; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_proteins_uniprot ON public.proteins USING btree (uniprot_id);


--
-- Name: idx_quality_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quality_compound ON public.social_content_quality USING btree (compound_id);


--
-- Name: idx_quality_date; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quality_date ON public.social_content_quality USING btree (analysis_date);


--
-- Name: idx_quality_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quality_platform ON public.social_content_quality USING btree (platform);


--
-- Name: idx_quality_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quality_subdivision ON public.social_content_quality USING btree (platform_subdivision);


--
-- Name: idx_quantum_basis_sets_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_basis_sets_type ON public.quantum_basis_sets USING btree (basis_type);


--
-- Name: idx_quantum_calculations_basis; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_calculations_basis ON public.quantum_calculations USING btree (basis_set);


--
-- Name: idx_quantum_calculations_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_calculations_compound ON public.quantum_calculations USING btree (compound_id);


--
-- Name: idx_quantum_calculations_functional; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_calculations_functional ON public.quantum_calculations USING btree (functional);


--
-- Name: idx_quantum_calculations_project; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_calculations_project ON public.quantum_calculations USING btree (project_id);


--
-- Name: idx_quantum_correlations_property; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_correlations_property ON public.quantum_structure_correlations USING btree (property_type);


--
-- Name: idx_quantum_critical_params_boundaries; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_critical_params_boundaries ON public.quantum_critical_params USING gin (phase_boundaries);


--
-- Name: idx_quantum_critical_params_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_critical_params_compound ON public.quantum_critical_params USING btree (compound_id);


--
-- Name: idx_quantum_dynamics_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_dynamics_compound ON public.quantum_dynamics USING btree (compound_id);


--
-- Name: idx_quantum_dynamics_evolution; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_dynamics_evolution ON public.quantum_dynamics USING gin (wavefunction_evolution);


--
-- Name: idx_quantum_findings_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_findings_type ON public.quantum_research_findings USING btree (finding_type);


--
-- Name: idx_quantum_functionals_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_functionals_type ON public.quantum_functionals USING btree (functional_type);


--
-- Name: idx_quantum_observables_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_observables_compound ON public.quantum_observables USING btree (compound_id);


--
-- Name: idx_quantum_observables_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_observables_type ON public.quantum_observables_ref USING btree (operator_type);


--
-- Name: idx_quantum_params_method; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_params_method ON public.quantum_parameters USING btree (calculation_method);


--
-- Name: idx_quantum_params_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_params_name ON public.quantum_parameters USING btree (parameter_name);


--
-- Name: idx_quantum_properties_basis; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_properties_basis ON public.quantum_properties USING btree (basis_set);


--
-- Name: idx_quantum_properties_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_properties_compound ON public.quantum_properties USING btree (compound_id);


--
-- Name: idx_quantum_properties_method; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_properties_method ON public.quantum_properties USING btree (method);


--
-- Name: idx_quantum_properties_transition; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_quantum_properties_transition ON public.quantum_properties USING btree (transition_type);


--
-- Name: idx_rate_limits_endpoint; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_rate_limits_endpoint ON public.rate_limits USING btree (endpoint_id);


--
-- Name: idx_rate_limits_key; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_rate_limits_key ON public.rate_limits USING btree (api_key_id);


--
-- Name: idx_rate_limits_window; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_rate_limits_window ON public.rate_limits USING btree (window_start, window_end);


--
-- Name: idx_receptor_families_category; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_receptor_families_category ON public.receptor_families USING btree (category_id);


--
-- Name: idx_receptor_families_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_receptor_families_name ON public.receptor_families USING btree (name);


--
-- Name: idx_receptor_families_signaling; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_receptor_families_signaling ON public.receptor_families USING btree (signaling_type);


--
-- Name: idx_receptor_families_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_receptor_families_type ON public.receptor_families USING btree (protein_type);


--
-- Name: idx_receptor_family_categories_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_receptor_family_categories_name ON public.receptor_family_categories USING btree (name);


--
-- Name: idx_receptor_family_categories_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_receptor_family_categories_type ON public.receptor_family_categories USING btree (receptor_type);


--
-- Name: idx_receptor_subtypes_family; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_receptor_subtypes_family ON public.receptor_subtypes USING btree (family_id);


--
-- Name: idx_receptor_subtypes_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_receptor_subtypes_name ON public.receptor_subtypes USING btree (subtype_name);


--
-- Name: idx_reddit_content_search; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_reddit_content_search ON public.social_posts USING gin (to_tsvector('english'::regconfig, content)) WHERE (platform = 'reddit'::text);


--
-- Name: idx_reddit_score; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_reddit_score ON public.social_posts USING btree (((engagement_metrics ->> 'score'::text))) WHERE (platform = 'reddit'::text);


--
-- Name: idx_research_findings_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_research_findings_compound ON public.research_findings USING btree (compound_id);


--
-- Name: idx_research_findings_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_research_findings_type ON public.research_findings USING btree (finding_type);


--
-- Name: idx_safety_incidents_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_safety_incidents_compound ON public.social_safety_incidents USING btree (compound_id);


--
-- Name: idx_safety_incidents_platforms; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_safety_incidents_platforms ON public.social_safety_incidents USING gin (platforms);


--
-- Name: idx_safety_incidents_resolution; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_safety_incidents_resolution ON public.social_safety_incidents USING btree (resolution_status);


--
-- Name: idx_safety_incidents_severity; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_safety_incidents_severity ON public.social_safety_incidents USING btree (severity);


--
-- Name: idx_safety_incidents_status; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_safety_incidents_status ON public.social_safety_incidents USING btree (verification_status);


--
-- Name: idx_safety_incidents_subdivisions; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_safety_incidents_subdivisions ON public.social_safety_incidents USING gin (platform_subdivisions);


--
-- Name: idx_safety_incidents_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_safety_incidents_type ON public.social_safety_incidents USING btree (incident_type);


--
-- Name: idx_safety_thresholds_param; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_safety_thresholds_param ON public.safety_thresholds USING btree (parameter_name);


--
-- Name: idx_safety_thresholds_severity; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_safety_thresholds_severity ON public.safety_thresholds USING btree (severity_level);


--
-- Name: idx_sar_analysis_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_sar_analysis_compound ON public.sar_analysis USING btree (compound_id);


--
-- Name: idx_sar_analysis_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_sar_analysis_type ON public.sar_analysis USING btree (analysis_type);


--
-- Name: idx_sar_patterns_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_sar_patterns_type ON public.sar_patterns USING btree (pattern_type);


--
-- Name: idx_scaffold_analysis_smiles; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_scaffold_analysis_smiles ON public.scaffold_analysis USING btree (scaffold_smiles);


--
-- Name: idx_scaling_analysis_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_scaling_analysis_compound ON public.scaling_analysis USING btree (compound_id);


--
-- Name: idx_scaling_analysis_flow; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_scaling_analysis_flow ON public.scaling_analysis USING gin (rg_flow_equations);


--
-- Name: idx_scaling_analysis_transition; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_scaling_analysis_transition ON public.scaling_analysis USING btree (phase_transition_id);


--
-- Name: idx_schema_versions_applied; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_schema_versions_applied ON public.schema_versions USING btree (applied_at);


--
-- Name: idx_schema_versions_status; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_schema_versions_status ON public.schema_versions USING btree (status);


--
-- Name: idx_schema_versions_version; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_schema_versions_version ON public.schema_versions USING btree (version);


--
-- Name: idx_scientific_content_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_scientific_content_compound ON public.social_scientific_content USING btree (compound_id);


--
-- Name: idx_scientific_content_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_scientific_content_platform ON public.social_scientific_content USING btree (platform);


--
-- Name: idx_scientific_content_post; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_scientific_content_post ON public.social_scientific_content USING btree (post_id);


--
-- Name: idx_scientific_content_quality; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_scientific_content_quality ON public.social_scientific_content USING btree (quality_score);


--
-- Name: idx_scientific_content_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_scientific_content_subdivision ON public.social_scientific_content USING btree (platform_subdivision);


--
-- Name: idx_scientific_content_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_scientific_content_type ON public.social_scientific_content USING btree (content_type);


--
-- Name: idx_selection_date; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_selection_date ON public.feature_selection USING btree (selection_date);


--
-- Name: idx_selection_method; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_selection_method ON public.feature_selection USING btree (selection_method);


--
-- Name: idx_selection_model; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_selection_model ON public.feature_selection USING btree (model_id);


--
-- Name: idx_service_status_service; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_service_status_service ON public.service_status USING btree (service_id);


--
-- Name: idx_service_status_status; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_service_status_status ON public.service_status USING btree (status);


--
-- Name: idx_service_status_timestamp; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_service_status_timestamp ON public.service_status USING btree (check_timestamp);


--
-- Name: idx_services_active; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_services_active ON public.external_services USING btree (is_active);


--
-- Name: idx_services_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_services_name ON public.external_services USING btree (name);


--
-- Name: idx_services_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_services_type ON public.external_services USING btree (service_type);


--
-- Name: idx_settings_user; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_settings_user ON public.ui_settings USING btree (user_id);


--
-- Name: idx_social_comments_created; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_comments_created ON public.social_comments USING btree (created_at);


--
-- Name: idx_social_comments_external; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_comments_external ON public.social_comments USING btree (external_id);


--
-- Name: idx_social_comments_parent; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_comments_parent ON public.social_comments USING btree (parent_id);


--
-- Name: idx_social_comments_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_comments_platform ON public.social_comments USING btree (platform);


--
-- Name: idx_social_comments_post; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_comments_post ON public.social_comments USING btree (post_id);


--
-- Name: idx_social_comments_scientific; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_comments_scientific ON public.social_comments USING btree (is_scientific);


--
-- Name: idx_social_comments_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_comments_subdivision ON public.social_comments USING btree (platform_subdivision);


--
-- Name: idx_social_posts_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_posts_compound ON public.social_posts USING btree (compound_id);


--
-- Name: idx_social_posts_content_search; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_posts_content_search ON public.social_posts USING gin (to_tsvector('english'::regconfig, content));


--
-- Name: idx_social_posts_created; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_posts_created ON public.social_posts USING btree (created_at);


--
-- Name: idx_social_posts_external_id; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_posts_external_id ON public.social_posts USING btree (external_id);


--
-- Name: idx_social_posts_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_posts_platform ON public.social_posts USING btree (platform);


--
-- Name: idx_social_posts_scientific; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_posts_scientific ON public.social_posts USING btree (is_scientific);


--
-- Name: idx_social_posts_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_posts_subdivision ON public.social_posts USING btree (platform_subdivision);


--
-- Name: idx_social_posts_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_posts_type ON public.social_posts USING btree (post_type);


--
-- Name: idx_social_protocols_compounds; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_protocols_compounds ON public.social_protocols USING gin (compounds);


--
-- Name: idx_social_protocols_outcome; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_protocols_outcome ON public.social_protocols USING btree (target_outcome);


--
-- Name: idx_social_protocols_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_protocols_platform ON public.social_protocols USING btree (platform);


--
-- Name: idx_social_protocols_score; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_protocols_score ON public.social_protocols USING btree (review_score);


--
-- Name: idx_social_protocols_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_protocols_subdivision ON public.social_protocols USING btree (platform_subdivision);


--
-- Name: idx_social_research_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_research_compound ON public.social_research_reviews USING btree (compound_id);


--
-- Name: idx_social_research_date; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_research_date ON public.social_research_reviews USING btree (publication_date);


--
-- Name: idx_social_research_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_research_platform ON public.social_research_reviews USING btree (platform);


--
-- Name: idx_social_research_quality; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_research_quality ON public.social_research_reviews USING btree (research_quality_score);


--
-- Name: idx_social_research_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_research_subdivision ON public.social_research_reviews USING btree (platform_subdivision);


--
-- Name: idx_social_research_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_research_type ON public.social_research_reviews USING gin (study_type);


--
-- Name: idx_social_stacks_compounds; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_stacks_compounds ON public.social_stacks USING gin (compounds);


--
-- Name: idx_social_stacks_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_stacks_platform ON public.social_stacks USING btree (platform);


--
-- Name: idx_social_stacks_purpose; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_stacks_purpose ON public.social_stacks USING btree (purpose);


--
-- Name: idx_social_stacks_rating; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_stacks_rating ON public.social_stacks USING btree (rating);


--
-- Name: idx_social_stacks_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_social_stacks_subdivision ON public.social_stacks USING btree (platform_subdivision);


--
-- Name: idx_structure_similarity_compound1; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_structure_similarity_compound1 ON public.structure_similarity USING btree (compound_id_1);


--
-- Name: idx_structure_similarity_compound2; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_structure_similarity_compound2 ON public.structure_similarity USING btree (compound_id_2);


--
-- Name: idx_structure_similarity_metric; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_structure_similarity_metric ON public.structure_similarity USING btree (similarity_metric);


--
-- Name: idx_subjective_effects_category; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_subjective_effects_category ON public.subjective_effects USING btree (category_id);


--
-- Name: idx_subjective_effects_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_subjective_effects_name ON public.subjective_effects USING btree (name);


--
-- Name: idx_subjective_effects_text_search; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_subjective_effects_text_search ON public.subjective_effects USING gin (text_search_vector);


--
-- Name: idx_substance_data_class; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_substance_data_class ON public.social_substance_data USING gin (chemical_class, psychoactive_class);


--
-- Name: idx_substance_data_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_substance_data_compound ON public.social_substance_data USING btree (compound_id);


--
-- Name: idx_substance_data_external; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_substance_data_external ON public.social_substance_data USING btree (external_id);


--
-- Name: idx_substance_data_names; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_substance_data_names ON public.social_substance_data USING gin (common_names);


--
-- Name: idx_substance_data_platform; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_substance_data_platform ON public.social_substance_data USING btree (platform);


--
-- Name: idx_substance_data_subdivision; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_substance_data_subdivision ON public.social_substance_data USING btree (platform_subdivision);


--
-- Name: idx_sync_status_device; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_sync_status_device ON public.sync_status USING btree (device_id);


--
-- Name: idx_sync_status_status; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_sync_status_status ON public.sync_status USING btree (status);


--
-- Name: idx_sync_status_sync; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_sync_status_sync ON public.sync_status USING btree (last_sync_at);


--
-- Name: idx_sync_status_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_sync_status_type ON public.sync_status USING btree (data_type);


--
-- Name: idx_systems_biology_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_systems_biology_compound ON public.systems_biology_analysis USING btree (compound_id);


--
-- Name: idx_systems_biology_level; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_systems_biology_level ON public.systems_biology_analysis USING btree (analysis_level);


--
-- Name: idx_templates_active; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_templates_active ON public.web_templates USING btree (is_active);


--
-- Name: idx_templates_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_templates_name ON public.web_templates USING btree (name);


--
-- Name: idx_templates_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_templates_type ON public.web_templates USING btree (template_type);


--
-- Name: idx_therapeutic_class_categories_level; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_therapeutic_class_categories_level ON public.therapeutic_class_categories USING btree (level);


--
-- Name: idx_therapeutic_class_categories_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_therapeutic_class_categories_name ON public.therapeutic_class_categories USING btree (name);


--
-- Name: idx_therapeutic_class_categories_parent; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_therapeutic_class_categories_parent ON public.therapeutic_class_categories USING btree (parent_category_id);


--
-- Name: idx_therapeutic_classes_category; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_therapeutic_classes_category ON public.therapeutic_classes USING btree (category_id);


--
-- Name: idx_therapeutic_classes_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_therapeutic_classes_name ON public.therapeutic_classes USING btree (name);


--
-- Name: idx_therapeutic_classes_text_search; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_therapeutic_classes_text_search ON public.therapeutic_classes USING gin (text_search_vector);


--
-- Name: idx_toxicity_endpoint_categories_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_toxicity_endpoint_categories_name ON public.toxicity_endpoint_categories USING btree (name);


--
-- Name: idx_toxicity_endpoint_categories_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_toxicity_endpoint_categories_type ON public.toxicity_endpoint_categories USING btree (measurement_type);


--
-- Name: idx_toxicity_endpoints_category; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_toxicity_endpoints_category ON public.toxicity_endpoints USING btree (category_id);


--
-- Name: idx_toxicity_endpoints_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_toxicity_endpoints_name ON public.toxicity_endpoints USING btree (name);


--
-- Name: idx_toxicity_endpoints_unit; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_toxicity_endpoints_unit ON public.toxicity_endpoints USING btree (standard_unit);


--
-- Name: idx_toxicity_mechanism_categories_level; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_toxicity_mechanism_categories_level ON public.toxicity_mechanism_categories USING btree (level);


--
-- Name: idx_toxicity_mechanism_categories_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_toxicity_mechanism_categories_name ON public.toxicity_mechanism_categories USING btree (name);


--
-- Name: idx_toxicity_mechanism_categories_parent; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_toxicity_mechanism_categories_parent ON public.toxicity_mechanism_categories USING btree (parent_category_id);


--
-- Name: idx_toxicity_mechanism_categories_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_toxicity_mechanism_categories_type ON public.toxicity_mechanism_categories USING btree (mechanism_type);


--
-- Name: idx_toxicity_mechanisms_category; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_toxicity_mechanisms_category ON public.toxicity_mechanisms USING btree (category_id);


--
-- Name: idx_toxicity_mechanisms_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_toxicity_mechanisms_name ON public.toxicity_mechanisms USING btree (name);


--
-- Name: idx_toxicity_mechanisms_text_search; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_toxicity_mechanisms_text_search ON public.toxicity_mechanisms USING gin (text_search_vector);


--
-- Name: idx_training_datasets_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_training_datasets_name ON public.training_datasets USING btree (dataset_name);


--
-- Name: idx_training_datasets_source; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_training_datasets_source ON public.training_datasets USING btree (source);


--
-- Name: idx_training_datasets_version; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_training_datasets_version ON public.training_datasets USING btree (version);


--
-- Name: idx_training_date; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_training_date ON public.training_history USING btree (start_time);


--
-- Name: idx_training_model; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_training_model ON public.training_history USING btree (model_id);


--
-- Name: idx_training_status; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_training_status ON public.training_history USING btree (status);


--
-- Name: idx_training_version; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_training_version ON public.training_history USING btree (version_id);


--
-- Name: idx_trend_analysis_period; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_trend_analysis_period ON public.trend_analysis USING btree (time_period);


--
-- Name: idx_trend_analysis_priority; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_trend_analysis_priority ON public.trend_analysis USING btree (priority);


--
-- Name: idx_trend_analysis_status; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_trend_analysis_status ON public.trend_analysis USING btree (status);


--
-- Name: idx_trend_analysis_target; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_trend_analysis_target ON public.trend_analysis USING btree (target_type, target_id);


--
-- Name: idx_trend_analysis_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_trend_analysis_type ON public.trend_analysis USING btree (analysis_type);


--
-- Name: idx_trends_compound; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_trends_compound ON public.social_compound_trends USING btree (compound_id);


--
-- Name: idx_trends_correlation; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_trends_correlation ON public.social_compound_trends USING btree (correlation_strength);


--
-- Name: idx_trends_dates; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_trends_dates ON public.social_compound_trends USING btree (trend_start_date, trend_end_date);


--
-- Name: idx_trends_platforms; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_trends_platforms ON public.social_compound_trends USING gin (platforms);


--
-- Name: idx_trends_similarity; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_trends_similarity ON public.social_compound_trends USING btree (temporal_similarity);


--
-- Name: idx_twitter_content_search; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_twitter_content_search ON public.social_posts USING gin (to_tsvector('english'::regconfig, content)) WHERE (platform = 'twitter'::text);


--
-- Name: idx_twitter_engagement; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_twitter_engagement ON public.social_posts USING btree (((engagement_metrics ->> 'retweet_count'::text))) WHERE (platform = 'twitter'::text);


--
-- Name: idx_usage_statistics_action; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_usage_statistics_action ON public.usage_statistics USING btree (action);


--
-- Name: idx_usage_statistics_created; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_usage_statistics_created ON public.usage_statistics USING btree (created_at);


--
-- Name: idx_usage_statistics_device; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_usage_statistics_device ON public.usage_statistics USING btree (device_type);


--
-- Name: idx_usage_statistics_feature; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_usage_statistics_feature ON public.usage_statistics USING btree (feature);


--
-- Name: idx_usage_statistics_session; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_usage_statistics_session ON public.usage_statistics USING btree (session_id);


--
-- Name: idx_usage_statistics_user; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_usage_statistics_user ON public.usage_statistics USING btree (user_id);


--
-- Name: idx_user_preferences_device; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_user_preferences_device ON public.user_preferences USING btree (device_id);


--
-- Name: idx_user_preferences_user; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_user_preferences_user ON public.user_preferences USING btree (user_id);


--
-- Name: idx_validation_results_date; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_validation_results_date ON public.model_validation_results USING btree (validation_date);


--
-- Name: idx_validation_results_model; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_validation_results_model ON public.model_validation_results USING btree (model_id);


--
-- Name: idx_validation_results_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_validation_results_type ON public.model_validation_results USING btree (validation_type);


--
-- Name: idx_versions_deployed; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_versions_deployed ON public.model_versions USING btree (deployed_at);


--
-- Name: idx_versions_model; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_versions_model ON public.model_versions USING btree (model_id);


--
-- Name: idx_versions_status; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_versions_status ON public.model_versions USING btree (deployment_status);


--
-- Name: idx_web_source_categories_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_web_source_categories_name ON public.web_source_categories USING btree (name);


--
-- Name: idx_web_source_categories_type; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_web_source_categories_type ON public.web_source_categories USING btree (source_type);


--
-- Name: idx_web_sources_active; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_web_sources_active ON public.web_data_sources USING btree (active);


--
-- Name: idx_web_sources_category; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_web_sources_category ON public.web_data_sources USING btree (category_id);


--
-- Name: idx_web_sources_name; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_web_sources_name ON public.web_data_sources USING btree (name);


--
-- Name: idx_web_sources_text_search; Type: INDEX; Schema: public; Owner: armand
--

CREATE INDEX idx_web_sources_text_search ON public.web_data_sources USING gin (text_search_vector);


--
-- Name: quantum_calculation_summary _RETURN; Type: RULE; Schema: public; Owner: armand
--

CREATE OR REPLACE VIEW public.quantum_calculation_summary AS
 SELECT c.id,
    c.compound_id,
    c.calculation_type,
    c.basis_set,
    c.functional,
    c.energy_hartree,
    c.convergence_achieved,
    count(f.id) AS finding_count,
    string_agg(DISTINCT f.finding_type, ', '::text) AS finding_types
   FROM (public.quantum_calculations c
     LEFT JOIN public.quantum_research_findings f ON ((c.id = f.calculation_id)))
  GROUP BY c.id;


--
-- Name: activity_cliffs audit_activity_cliffs_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_activity_cliffs_trigger AFTER INSERT OR DELETE OR UPDATE ON public.activity_cliffs FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: activity_correlations audit_activity_correlations_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_activity_correlations_trigger AFTER INSERT OR DELETE OR UPDATE ON public.activity_correlations FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: audit_log audit_audit_log_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_audit_log_trigger AFTER INSERT OR DELETE OR UPDATE ON public.audit_log FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: binding_assay_protocols audit_binding_assay_protocols_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_binding_assay_protocols_trigger AFTER INSERT OR DELETE OR UPDATE ON public.binding_assay_protocols FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: binding_assay_types audit_binding_assay_types_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_binding_assay_types_trigger AFTER INSERT OR DELETE OR UPDATE ON public.binding_assay_types FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: binding_data_quality audit_binding_data_quality_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_binding_data_quality_trigger AFTER INSERT OR DELETE OR UPDATE ON public.binding_data_quality FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: binding_data audit_binding_data_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_binding_data_trigger AFTER INSERT OR DELETE OR UPDATE ON public.binding_data FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: binding_kinetics audit_binding_kinetics_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_binding_kinetics_trigger AFTER INSERT OR DELETE OR UPDATE ON public.binding_kinetics FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: binding_sar audit_binding_sar_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_binding_sar_trigger AFTER INSERT OR DELETE OR UPDATE ON public.binding_sar FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: binding_site_mapping audit_binding_site_mapping_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_binding_site_mapping_trigger AFTER INSERT OR DELETE OR UPDATE ON public.binding_site_mapping FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: carcinogenicity_data audit_carcinogenicity_data_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_carcinogenicity_data_trigger AFTER INSERT OR DELETE OR UPDATE ON public.carcinogenicity_data FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: cardiotoxicity_data audit_cardiotoxicity_data_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_cardiotoxicity_data_trigger AFTER INSERT OR DELETE OR UPDATE ON public.cardiotoxicity_data FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: compounds audit_compounds_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_compounds_trigger AFTER INSERT OR DELETE OR UPDATE ON public.compounds FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: cytotoxicity_data audit_cytotoxicity_data_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_cytotoxicity_data_trigger AFTER INSERT OR DELETE OR UPDATE ON public.cytotoxicity_data FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: descriptors_2d audit_descriptors_2d_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_descriptors_2d_trigger AFTER INSERT OR DELETE OR UPDATE ON public.descriptors_2d FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: descriptors_3d audit_descriptors_3d_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_descriptors_3d_trigger AFTER INSERT OR DELETE OR UPDATE ON public.descriptors_3d FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: developmental_toxicity audit_developmental_toxicity_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_developmental_toxicity_trigger AFTER INSERT OR DELETE OR UPDATE ON public.developmental_toxicity FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: electronic_structure audit_electronic_structure_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_electronic_structure_trigger AFTER INSERT OR DELETE OR UPDATE ON public.electronic_structure FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: emergency_procedures audit_emergency_procedures_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_emergency_procedures_trigger AFTER INSERT OR DELETE OR UPDATE ON public.emergency_procedures FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: genes audit_genes_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_genes_trigger AFTER INSERT OR DELETE OR UPDATE ON public.genes FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: genotoxicity_data audit_genotoxicity_data_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_genotoxicity_data_trigger AFTER INSERT OR DELETE OR UPDATE ON public.genotoxicity_data FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: hazard_classifications audit_hazard_classifications_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_hazard_classifications_trigger AFTER INSERT OR DELETE OR UPDATE ON public.hazard_classifications FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: hepatotoxicity_data audit_hepatotoxicity_data_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_hepatotoxicity_data_trigger AFTER INSERT OR DELETE OR UPDATE ON public.hepatotoxicity_data FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: immunotoxicity_data audit_immunotoxicity_data_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_immunotoxicity_data_trigger AFTER INSERT OR DELETE OR UPDATE ON public.immunotoxicity_data FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: literature_findings audit_literature_findings_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_literature_findings_trigger AFTER INSERT OR DELETE OR UPDATE ON public.literature_findings FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: meta_analyses audit_meta_analyses_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_meta_analyses_trigger AFTER INSERT OR DELETE OR UPDATE ON public.meta_analyses FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: metabolic_toxicity audit_metabolic_toxicity_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_metabolic_toxicity_trigger AFTER INSERT OR DELETE OR UPDATE ON public.metabolic_toxicity FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: molecular_fingerprints audit_molecular_fingerprints_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_molecular_fingerprints_trigger AFTER INSERT OR DELETE OR UPDATE ON public.molecular_fingerprints FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: monitoring_parameters audit_monitoring_params_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_monitoring_params_trigger AFTER INSERT OR DELETE OR UPDATE ON public.monitoring_parameters FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: neurotoxicity_data audit_neurotoxicity_data_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_neurotoxicity_data_trigger AFTER INSERT OR DELETE OR UPDATE ON public.neurotoxicity_data FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: organ_systems audit_organ_systems_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_organ_systems_trigger AFTER INSERT OR DELETE OR UPDATE ON public.organ_systems FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: organ_toxicity_patterns audit_organ_toxicity_patterns_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_organ_toxicity_patterns_trigger AFTER INSERT OR DELETE OR UPDATE ON public.organ_toxicity_patterns FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: pathway_analysis audit_pathway_analysis_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_pathway_analysis_trigger AFTER INSERT OR DELETE OR UPDATE ON public.pathway_analysis FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: pharmacological_classes audit_pharmacological_classes_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_pharmacological_classes_trigger AFTER INSERT OR DELETE OR UPDATE ON public.pharmacological_classes FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: pharmacophore_features audit_pharmacophore_features_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_pharmacophore_features_trigger AFTER INSERT OR DELETE OR UPDATE ON public.pharmacophore_features FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: phase_transitions audit_phase_transitions_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_phase_transitions_trigger AFTER INSERT OR DELETE OR UPDATE ON public.phase_transitions FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: ppe_requirements audit_ppe_requirements_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_ppe_requirements_trigger AFTER INSERT OR DELETE OR UPDATE ON public.ppe_requirements FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: proteins audit_proteins_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_proteins_trigger AFTER INSERT OR DELETE OR UPDATE ON public.proteins FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: quantum_calculations audit_quantum_calculations_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_quantum_calculations_trigger AFTER INSERT OR DELETE OR UPDATE ON public.quantum_calculations FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: quantum_critical_params audit_quantum_critical_params_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_quantum_critical_params_trigger AFTER INSERT OR DELETE OR UPDATE ON public.quantum_critical_params FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: quantum_dynamics audit_quantum_dynamics_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_quantum_dynamics_trigger AFTER INSERT OR DELETE OR UPDATE ON public.quantum_dynamics FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: quantum_observables audit_quantum_observables_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_quantum_observables_trigger AFTER INSERT OR DELETE OR UPDATE ON public.quantum_observables FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: quantum_properties audit_quantum_properties_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_quantum_properties_trigger AFTER INSERT OR DELETE OR UPDATE ON public.quantum_properties FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: quantum_research_findings audit_quantum_research_findings_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_quantum_research_findings_trigger AFTER INSERT OR DELETE OR UPDATE ON public.quantum_research_findings FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: quantum_research_projects audit_quantum_research_projects_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_quantum_research_projects_trigger AFTER INSERT OR DELETE OR UPDATE ON public.quantum_research_projects FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: quantum_structure_correlations audit_quantum_structure_correlations_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_quantum_structure_correlations_trigger AFTER INSERT OR DELETE OR UPDATE ON public.quantum_structure_correlations FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: receptor_families audit_receptor_families_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_receptor_families_trigger AFTER INSERT OR DELETE OR UPDATE ON public.receptor_families FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: receptor_family_categories audit_receptor_family_categories_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_receptor_family_categories_trigger AFTER INSERT OR DELETE OR UPDATE ON public.receptor_family_categories FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: receptor_subtypes audit_receptor_subtypes_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_receptor_subtypes_trigger AFTER INSERT OR DELETE OR UPDATE ON public.receptor_subtypes FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: research_findings audit_research_findings_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_research_findings_trigger AFTER INSERT OR DELETE OR UPDATE ON public.research_findings FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: safety_documents audit_safety_documents_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_safety_documents_trigger AFTER INSERT OR DELETE OR UPDATE ON public.safety_documents FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: safety_thresholds audit_safety_thresholds_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_safety_thresholds_trigger AFTER INSERT OR DELETE OR UPDATE ON public.safety_thresholds FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: sar_analysis audit_sar_analysis_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_sar_analysis_trigger AFTER INSERT OR DELETE OR UPDATE ON public.sar_analysis FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: sar_patterns audit_sar_patterns_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_sar_patterns_trigger AFTER INSERT OR DELETE OR UPDATE ON public.sar_patterns FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: scaffold_analysis audit_scaffold_analysis_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_scaffold_analysis_trigger AFTER INSERT OR DELETE OR UPDATE ON public.scaffold_analysis FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: scaling_analysis audit_scaling_analysis_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_scaling_analysis_trigger AFTER INSERT OR DELETE OR UPDATE ON public.scaling_analysis FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: schema_versions audit_schema_versions_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_schema_versions_trigger AFTER INSERT OR DELETE OR UPDATE ON public.schema_versions FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: scientific_papers audit_scientific_papers_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_scientific_papers_trigger AFTER INSERT OR DELETE OR UPDATE ON public.scientific_papers FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_alert_rules audit_social_alert_rules_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_alert_rules_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_alert_rules FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_comments audit_social_comments_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_comments_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_comments FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_compound_combinations audit_social_compound_combinations_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_compound_combinations_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_compound_combinations FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_compound_mentions audit_social_compound_mentions_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_compound_mentions_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_compound_mentions FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_compound_trends audit_social_compound_trends_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_compound_trends_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_compound_trends FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_content_quality audit_social_content_quality_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_content_quality_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_content_quality FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_dosage_data audit_social_dosage_data_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_dosage_data_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_dosage_data FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_effect_reports audit_social_effect_reports_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_effect_reports_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_effect_reports FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_effects audit_social_effects_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_effects_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_effects FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_experience_reports audit_social_experience_reports_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_experience_reports_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_experience_reports FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_harm_reduction audit_social_harm_reduction_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_harm_reduction_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_harm_reduction FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_influence_networks audit_social_influence_networks_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_influence_networks_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_influence_networks FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_platform_stats audit_social_platform_stats_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_platform_stats_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_platform_stats FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_posts audit_social_posts_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_posts_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_posts FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_protocols audit_social_protocols_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_protocols_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_protocols FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_research_reviews audit_social_research_reviews_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_research_reviews_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_research_reviews FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_safety_incidents audit_social_safety_incidents_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_safety_incidents_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_safety_incidents FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_scientific_content audit_social_scientific_content_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_scientific_content_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_scientific_content FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_stacks audit_social_stacks_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_stacks_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_stacks FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: social_substance_data audit_social_substance_data_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_social_substance_data_trigger AFTER INSERT OR DELETE OR UPDATE ON public.social_substance_data FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: storage_requirements audit_storage_requirements_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_storage_requirements_trigger AFTER INSERT OR DELETE OR UPDATE ON public.storage_requirements FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: structure_similarity audit_structure_similarity_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_structure_similarity_trigger AFTER INSERT OR DELETE OR UPDATE ON public.structure_similarity FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: systems_biology_analysis audit_systems_biology_analysis_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_systems_biology_analysis_trigger AFTER INSERT OR DELETE OR UPDATE ON public.systems_biology_analysis FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: therapeutic_class_categories audit_therapeutic_class_categories_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_therapeutic_class_categories_trigger AFTER INSERT OR DELETE OR UPDATE ON public.therapeutic_class_categories FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: therapeutic_classes audit_therapeutic_classes_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_therapeutic_classes_trigger AFTER INSERT OR DELETE OR UPDATE ON public.therapeutic_classes FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: toxicity_assays audit_toxicity_assays_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_toxicity_assays_trigger AFTER INSERT OR DELETE OR UPDATE ON public.toxicity_assays FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: toxicity_endpoint_categories audit_toxicity_endpoint_categories_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_toxicity_endpoint_categories_trigger AFTER INSERT OR DELETE OR UPDATE ON public.toxicity_endpoint_categories FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: toxicity_endpoints audit_toxicity_endpoints_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_toxicity_endpoints_trigger AFTER INSERT OR DELETE OR UPDATE ON public.toxicity_endpoints FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: toxicity_mechanism_categories audit_toxicity_mechanism_categories_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_toxicity_mechanism_categories_trigger AFTER INSERT OR DELETE OR UPDATE ON public.toxicity_mechanism_categories FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: toxicity_mechanisms audit_toxicity_mechanisms_trigger; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER audit_toxicity_mechanisms_trigger AFTER INSERT OR DELETE OR UPDATE ON public.toxicity_mechanisms FOR EACH ROW EXECUTE FUNCTION public.audit_trigger_func();


--
-- Name: quantum_critical_params trigger_quantum_critical_updates; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER trigger_quantum_critical_updates AFTER INSERT OR UPDATE ON public.quantum_critical_params FOR EACH ROW EXECUTE FUNCTION public.update_quantum_properties();


--
-- Name: quantum_dynamics trigger_quantum_dynamics_updates; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER trigger_quantum_dynamics_updates AFTER INSERT OR UPDATE ON public.quantum_dynamics FOR EACH ROW EXECUTE FUNCTION public.update_quantum_properties();


--
-- Name: activity_cliffs update_activity_cliffs_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_activity_cliffs_modtime BEFORE UPDATE ON public.activity_cliffs FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: activity_correlations update_activity_correlations_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_activity_correlations_modtime BEFORE UPDATE ON public.activity_correlations FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: alert_triggers update_alert_triggers_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_alert_triggers_modtime BEFORE UPDATE ON public.alert_triggers FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: api_keys update_api_keys_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_api_keys_modtime BEFORE UPDATE ON public.api_keys FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: binding_assay_protocols update_binding_assay_protocols_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_binding_assay_protocols_modtime BEFORE UPDATE ON public.binding_assay_protocols FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: binding_assay_types update_binding_assay_types_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_binding_assay_types_modtime BEFORE UPDATE ON public.binding_assay_types FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: binding_data update_binding_data_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_binding_data_modtime BEFORE UPDATE ON public.binding_data FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: binding_data_quality update_binding_data_quality_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_binding_data_quality_modtime BEFORE UPDATE ON public.binding_data_quality FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: binding_kinetics update_binding_kinetics_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_binding_kinetics_modtime BEFORE UPDATE ON public.binding_kinetics FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: binding_sar update_binding_sar_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_binding_sar_modtime BEFORE UPDATE ON public.binding_sar FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: binding_site_mapping update_binding_site_mapping_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_binding_site_mapping_modtime BEFORE UPDATE ON public.binding_site_mapping FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: web_components update_components_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_components_modtime BEFORE UPDATE ON public.web_components FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: compounds update_compounds_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_compounds_modtime BEFORE UPDATE ON public.compounds FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: descriptors_2d update_descriptors_2d_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_descriptors_2d_modtime BEFORE UPDATE ON public.descriptors_2d FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: descriptors_3d update_descriptors_3d_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_descriptors_3d_modtime BEFORE UPDATE ON public.descriptors_3d FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: device_settings update_device_settings_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_device_settings_modtime BEFORE UPDATE ON public.device_settings FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: electronic_structure update_electronic_structure_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_electronic_structure_modtime BEFORE UPDATE ON public.electronic_structure FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: api_endpoints update_endpoints_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_endpoints_modtime BEFORE UPDATE ON public.api_endpoints FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: feature_values update_feature_values_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_feature_values_modtime BEFORE UPDATE ON public.feature_values FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: feature_definitions update_features_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_features_modtime BEFORE UPDATE ON public.feature_definitions FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: genes update_genes_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_genes_modtime BEFORE UPDATE ON public.genes FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: web_hook_deliveries update_hook_deliveries_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_hook_deliveries_modtime BEFORE UPDATE ON public.web_hook_deliveries FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: web_hooks update_hooks_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_hooks_modtime BEFORE UPDATE ON public.web_hooks FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: literature_findings update_literature_findings_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_literature_findings_modtime BEFORE UPDATE ON public.literature_findings FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: meta_analyses update_meta_analyses_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_meta_analyses_modtime BEFORE UPDATE ON public.meta_analyses FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: model_metrics update_metrics_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_metrics_modtime BEFORE UPDATE ON public.model_metrics FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: ml_models update_models_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_models_modtime BEFORE UPDATE ON public.ml_models FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: molecular_fingerprints update_molecular_fingerprints_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_molecular_fingerprints_modtime BEFORE UPDATE ON public.molecular_fingerprints FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: monitoring_parameters update_monitoring_params_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_monitoring_params_modtime BEFORE UPDATE ON public.monitoring_parameters FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: offline_data update_offline_data_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_offline_data_modtime BEFORE UPDATE ON public.offline_data FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: organ_systems update_organ_systems_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_organ_systems_modtime BEFORE UPDATE ON public.organ_systems FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: organ_toxicity_patterns update_organ_toxicity_patterns_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_organ_toxicity_patterns_modtime BEFORE UPDATE ON public.organ_toxicity_patterns FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: pathway_analysis update_pathway_analysis_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_pathway_analysis_modtime BEFORE UPDATE ON public.pathway_analysis FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: performance_metrics update_performance_metrics_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_performance_metrics_modtime BEFORE UPDATE ON public.performance_metrics FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: pharmacological_classes update_pharmacological_classes_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_pharmacological_classes_modtime BEFORE UPDATE ON public.pharmacological_classes FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: pharmacophore_features update_pharmacophore_features_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_pharmacophore_features_modtime BEFORE UPDATE ON public.pharmacophore_features FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: phase_transitions update_phase_transitions_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_phase_transitions_modtime BEFORE UPDATE ON public.phase_transitions FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: feature_pipelines update_pipelines_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_pipelines_modtime BEFORE UPDATE ON public.feature_pipelines FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: model_predictions update_predictions_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_predictions_modtime BEFORE UPDATE ON public.model_predictions FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: proteins update_proteins_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_proteins_modtime BEFORE UPDATE ON public.proteins FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: quantum_calculations update_quantum_calculations_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_quantum_calculations_modtime BEFORE UPDATE ON public.quantum_calculations FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: quantum_critical_params update_quantum_critical_params_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_quantum_critical_params_modtime BEFORE UPDATE ON public.quantum_critical_params FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: quantum_dynamics update_quantum_dynamics_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_quantum_dynamics_modtime BEFORE UPDATE ON public.quantum_dynamics FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: quantum_observables update_quantum_observables_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_quantum_observables_modtime BEFORE UPDATE ON public.quantum_observables FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: quantum_properties update_quantum_properties_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_quantum_properties_modtime BEFORE UPDATE ON public.quantum_properties FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: quantum_research_findings update_quantum_research_findings_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_quantum_research_findings_modtime BEFORE UPDATE ON public.quantum_research_findings FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: quantum_research_projects update_quantum_research_projects_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_quantum_research_projects_modtime BEFORE UPDATE ON public.quantum_research_projects FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: quantum_structure_correlations update_quantum_structure_correlations_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_quantum_structure_correlations_modtime BEFORE UPDATE ON public.quantum_structure_correlations FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: rate_limits update_rate_limits_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_rate_limits_modtime BEFORE UPDATE ON public.rate_limits FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: receptor_families update_receptor_families_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_receptor_families_modtime BEFORE UPDATE ON public.receptor_families FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: receptor_family_categories update_receptor_family_categories_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_receptor_family_categories_modtime BEFORE UPDATE ON public.receptor_family_categories FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: receptor_subtypes update_receptor_subtypes_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_receptor_subtypes_modtime BEFORE UPDATE ON public.receptor_subtypes FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: research_findings update_research_findings_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_research_findings_modtime BEFORE UPDATE ON public.research_findings FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: safety_thresholds update_safety_thresholds_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_safety_thresholds_modtime BEFORE UPDATE ON public.safety_thresholds FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: sar_analysis update_sar_analysis_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_sar_analysis_modtime BEFORE UPDATE ON public.sar_analysis FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: sar_patterns update_sar_patterns_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_sar_patterns_modtime BEFORE UPDATE ON public.sar_patterns FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: scaffold_analysis update_scaffold_analysis_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_scaffold_analysis_modtime BEFORE UPDATE ON public.scaffold_analysis FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: scaling_analysis update_scaling_analysis_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_scaling_analysis_modtime BEFORE UPDATE ON public.scaling_analysis FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: scientific_papers update_scientific_papers_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_scientific_papers_modtime BEFORE UPDATE ON public.scientific_papers FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: feature_selection update_selection_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_selection_modtime BEFORE UPDATE ON public.feature_selection FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: service_status update_service_status_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_service_status_modtime BEFORE UPDATE ON public.service_status FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: external_services update_services_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_services_modtime BEFORE UPDATE ON public.external_services FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: ui_settings update_settings_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_settings_modtime BEFORE UPDATE ON public.ui_settings FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: structure_similarity update_structure_similarity_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_structure_similarity_modtime BEFORE UPDATE ON public.structure_similarity FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: sync_status update_sync_status_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_sync_status_modtime BEFORE UPDATE ON public.sync_status FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: systems_biology_analysis update_systems_biology_analysis_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_systems_biology_analysis_modtime BEFORE UPDATE ON public.systems_biology_analysis FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: web_templates update_templates_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_templates_modtime BEFORE UPDATE ON public.web_templates FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: therapeutic_class_categories update_therapeutic_class_categories_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_therapeutic_class_categories_modtime BEFORE UPDATE ON public.therapeutic_class_categories FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: therapeutic_classes update_therapeutic_classes_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_therapeutic_classes_modtime BEFORE UPDATE ON public.therapeutic_classes FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: toxicity_endpoint_categories update_toxicity_endpoint_categories_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_toxicity_endpoint_categories_modtime BEFORE UPDATE ON public.toxicity_endpoint_categories FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: toxicity_endpoints update_toxicity_endpoints_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_toxicity_endpoints_modtime BEFORE UPDATE ON public.toxicity_endpoints FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: toxicity_mechanism_categories update_toxicity_mechanism_categories_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_toxicity_mechanism_categories_modtime BEFORE UPDATE ON public.toxicity_mechanism_categories FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: toxicity_mechanisms update_toxicity_mechanisms_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_toxicity_mechanisms_modtime BEFORE UPDATE ON public.toxicity_mechanisms FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: training_history update_training_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_training_modtime BEFORE UPDATE ON public.training_history FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: trend_analysis update_trend_analysis_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_trend_analysis_modtime BEFORE UPDATE ON public.trend_analysis FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: usage_statistics update_usage_statistics_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_usage_statistics_modtime BEFORE UPDATE ON public.usage_statistics FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: user_preferences update_user_preferences_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_user_preferences_modtime BEFORE UPDATE ON public.user_preferences FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: model_versions update_versions_modtime; Type: TRIGGER; Schema: public; Owner: armand
--

CREATE TRIGGER update_versions_modtime BEFORE UPDATE ON public.model_versions FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();


--
-- Name: activity_correlations activity_correlations_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.activity_correlations
    ADD CONSTRAINT activity_correlations_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: binding_data binding_data_assay_type_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_data
    ADD CONSTRAINT binding_data_assay_type_id_fkey FOREIGN KEY (assay_type_id) REFERENCES public.binding_assay_types(id);


--
-- Name: binding_data binding_data_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_data
    ADD CONSTRAINT binding_data_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: binding_data_quality binding_data_quality_binding_data_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_data_quality
    ADD CONSTRAINT binding_data_quality_binding_data_id_fkey FOREIGN KEY (binding_data_id) REFERENCES public.binding_data(id) ON DELETE CASCADE;


--
-- Name: binding_data binding_data_receptor_family_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_data
    ADD CONSTRAINT binding_data_receptor_family_id_fkey FOREIGN KEY (receptor_family_id) REFERENCES public.receptor_families(id);


--
-- Name: binding_kinetics binding_kinetics_binding_data_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_kinetics
    ADD CONSTRAINT binding_kinetics_binding_data_id_fkey FOREIGN KEY (binding_data_id) REFERENCES public.binding_data(id) ON DELETE CASCADE;


--
-- Name: binding_sar binding_sar_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_sar
    ADD CONSTRAINT binding_sar_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: binding_sar binding_sar_receptor_family_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_sar
    ADD CONSTRAINT binding_sar_receptor_family_id_fkey FOREIGN KEY (receptor_family_id) REFERENCES public.receptor_families(id);


--
-- Name: binding_site_mapping binding_site_mapping_binding_data_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.binding_site_mapping
    ADD CONSTRAINT binding_site_mapping_binding_data_id_fkey FOREIGN KEY (binding_data_id) REFERENCES public.binding_data(id) ON DELETE CASCADE;


--
-- Name: carcinogenicity_data carcinogenicity_data_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.carcinogenicity_data
    ADD CONSTRAINT carcinogenicity_data_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: cardiotoxicity_data cardiotoxicity_data_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.cardiotoxicity_data
    ADD CONSTRAINT cardiotoxicity_data_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: cytotoxicity_data cytotoxicity_data_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.cytotoxicity_data
    ADD CONSTRAINT cytotoxicity_data_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: descriptors_2d descriptors_2d_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.descriptors_2d
    ADD CONSTRAINT descriptors_2d_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: descriptors_3d descriptors_3d_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.descriptors_3d
    ADD CONSTRAINT descriptors_3d_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: developmental_toxicity developmental_toxicity_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.developmental_toxicity
    ADD CONSTRAINT developmental_toxicity_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: effect_relationships effect_relationships_effect_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.effect_relationships
    ADD CONSTRAINT effect_relationships_effect_id_fkey FOREIGN KEY (effect_id) REFERENCES public.subjective_effects(id);


--
-- Name: effect_relationships effect_relationships_related_effect_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.effect_relationships
    ADD CONSTRAINT effect_relationships_related_effect_id_fkey FOREIGN KEY (related_effect_id) REFERENCES public.subjective_effects(id);


--
-- Name: electronic_structure electronic_structure_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.electronic_structure
    ADD CONSTRAINT electronic_structure_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: emergency_procedures emergency_procedures_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.emergency_procedures
    ADD CONSTRAINT emergency_procedures_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: energy_level_statistics energy_level_statistics_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.energy_level_statistics
    ADD CONSTRAINT energy_level_statistics_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: experience_categories experience_categories_parent_category_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.experience_categories
    ADD CONSTRAINT experience_categories_parent_category_fkey FOREIGN KEY (parent_category) REFERENCES public.experience_categories(category_name);


--
-- Name: feature_selection feature_selection_model_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.feature_selection
    ADD CONSTRAINT feature_selection_model_id_fkey FOREIGN KEY (model_id) REFERENCES public.ml_models(id);


--
-- Name: feature_values feature_values_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.feature_values
    ADD CONSTRAINT feature_values_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: feature_values feature_values_feature_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.feature_values
    ADD CONSTRAINT feature_values_feature_id_fkey FOREIGN KEY (feature_id) REFERENCES public.feature_definitions(id);


--
-- Name: electronic_structure fk_basis_set; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.electronic_structure
    ADD CONSTRAINT fk_basis_set FOREIGN KEY (basis_set) REFERENCES public.quantum_basis_sets(name);


--
-- Name: electronic_structure fk_functional; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.electronic_structure
    ADD CONSTRAINT fk_functional FOREIGN KEY (method) REFERENCES public.quantum_functionals(name);


--
-- Name: quantum_observables fk_observable_type; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_observables
    ADD CONSTRAINT fk_observable_type FOREIGN KEY (observable_type) REFERENCES public.quantum_observables_ref(name);


--
-- Name: phase_transitions fk_transition_type; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.phase_transitions
    ADD CONSTRAINT fk_transition_type FOREIGN KEY (transition_type) REFERENCES public.phase_transition_types(name);


--
-- Name: genotoxicity_data genotoxicity_data_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.genotoxicity_data
    ADD CONSTRAINT genotoxicity_data_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: hazard_classifications hazard_classifications_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.hazard_classifications
    ADD CONSTRAINT hazard_classifications_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: hepatotoxicity_data hepatotoxicity_data_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.hepatotoxicity_data
    ADD CONSTRAINT hepatotoxicity_data_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: immunotoxicity_data immunotoxicity_data_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.immunotoxicity_data
    ADD CONSTRAINT immunotoxicity_data_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: literature_findings literature_findings_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.literature_findings
    ADD CONSTRAINT literature_findings_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: literature_findings literature_findings_paper_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.literature_findings
    ADD CONSTRAINT literature_findings_paper_id_fkey FOREIGN KEY (paper_id) REFERENCES public.scientific_papers(id) ON DELETE CASCADE;


--
-- Name: metabolic_toxicity metabolic_toxicity_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.metabolic_toxicity
    ADD CONSTRAINT metabolic_toxicity_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: model_deployments model_deployments_model_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_deployments
    ADD CONSTRAINT model_deployments_model_id_fkey FOREIGN KEY (model_id) REFERENCES public.ml_models(id) ON DELETE CASCADE;


--
-- Name: model_metrics model_metrics_model_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_metrics
    ADD CONSTRAINT model_metrics_model_id_fkey FOREIGN KEY (model_id) REFERENCES public.ml_models(id);


--
-- Name: model_metrics model_metrics_version_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_metrics
    ADD CONSTRAINT model_metrics_version_id_fkey FOREIGN KEY (version_id) REFERENCES public.model_versions(id);


--
-- Name: model_monitoring model_monitoring_model_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_monitoring
    ADD CONSTRAINT model_monitoring_model_id_fkey FOREIGN KEY (model_id) REFERENCES public.ml_models(id) ON DELETE CASCADE;


--
-- Name: model_predictions model_predictions_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_predictions
    ADD CONSTRAINT model_predictions_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: model_predictions model_predictions_model_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_predictions
    ADD CONSTRAINT model_predictions_model_id_fkey FOREIGN KEY (model_id) REFERENCES public.ml_models(id);


--
-- Name: model_predictions model_predictions_version_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_predictions
    ADD CONSTRAINT model_predictions_version_id_fkey FOREIGN KEY (version_id) REFERENCES public.model_versions(id);


--
-- Name: model_validation_results model_validation_results_model_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_validation_results
    ADD CONSTRAINT model_validation_results_model_id_fkey FOREIGN KEY (model_id) REFERENCES public.ml_models(id) ON DELETE CASCADE;


--
-- Name: model_versions model_versions_model_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.model_versions
    ADD CONSTRAINT model_versions_model_id_fkey FOREIGN KEY (model_id) REFERENCES public.ml_models(id) ON DELETE CASCADE;


--
-- Name: molecular_fingerprints molecular_fingerprints_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.molecular_fingerprints
    ADD CONSTRAINT molecular_fingerprints_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: neurotoxicity_data neurotoxicity_data_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.neurotoxicity_data
    ADD CONSTRAINT neurotoxicity_data_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: organ_toxicity_patterns organ_toxicity_patterns_organ_system_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.organ_toxicity_patterns
    ADD CONSTRAINT organ_toxicity_patterns_organ_system_id_fkey FOREIGN KEY (organ_system_id) REFERENCES public.organ_systems(id);


--
-- Name: pathway_analysis pathway_analysis_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.pathway_analysis
    ADD CONSTRAINT pathway_analysis_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: pharmacophore_features pharmacophore_features_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.pharmacophore_features
    ADD CONSTRAINT pharmacophore_features_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: phase_transitions phase_transitions_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.phase_transitions
    ADD CONSTRAINT phase_transitions_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: ppe_requirements ppe_requirements_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.ppe_requirements
    ADD CONSTRAINT ppe_requirements_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: proteins proteins_gene_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.proteins
    ADD CONSTRAINT proteins_gene_id_fkey FOREIGN KEY (gene_id) REFERENCES public.genes(id);


--
-- Name: quantum_calculations quantum_calculations_basis_set_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_calculations
    ADD CONSTRAINT quantum_calculations_basis_set_fkey FOREIGN KEY (basis_set) REFERENCES public.quantum_basis_sets(name);


--
-- Name: quantum_calculations quantum_calculations_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_calculations
    ADD CONSTRAINT quantum_calculations_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id);


--
-- Name: quantum_calculations quantum_calculations_functional_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_calculations
    ADD CONSTRAINT quantum_calculations_functional_fkey FOREIGN KEY (functional) REFERENCES public.quantum_functionals(name);


--
-- Name: quantum_calculations quantum_calculations_project_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_calculations
    ADD CONSTRAINT quantum_calculations_project_id_fkey FOREIGN KEY (project_id) REFERENCES public.quantum_research_projects(id);


--
-- Name: quantum_critical_params quantum_critical_params_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_critical_params
    ADD CONSTRAINT quantum_critical_params_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: quantum_dynamics quantum_dynamics_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_dynamics
    ADD CONSTRAINT quantum_dynamics_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: quantum_hamiltonians quantum_hamiltonians_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_hamiltonians
    ADD CONSTRAINT quantum_hamiltonians_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: quantum_observables quantum_observables_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_observables
    ADD CONSTRAINT quantum_observables_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: quantum_properties quantum_properties_basis_set_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_properties
    ADD CONSTRAINT quantum_properties_basis_set_fkey FOREIGN KEY (basis_set) REFERENCES public.quantum_basis_sets(name);


--
-- Name: quantum_properties quantum_properties_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_properties
    ADD CONSTRAINT quantum_properties_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: quantum_properties quantum_properties_method_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_properties
    ADD CONSTRAINT quantum_properties_method_fkey FOREIGN KEY (method) REFERENCES public.quantum_functionals(name);


--
-- Name: quantum_research_findings quantum_research_findings_calculation_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_research_findings
    ADD CONSTRAINT quantum_research_findings_calculation_id_fkey FOREIGN KEY (calculation_id) REFERENCES public.quantum_calculations(id);


--
-- Name: quantum_structure_correlations quantum_structure_correlations_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.quantum_structure_correlations
    ADD CONSTRAINT quantum_structure_correlations_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id);


--
-- Name: rate_limits rate_limits_api_key_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.rate_limits
    ADD CONSTRAINT rate_limits_api_key_id_fkey FOREIGN KEY (api_key_id) REFERENCES public.api_keys(id) ON DELETE CASCADE;


--
-- Name: rate_limits rate_limits_endpoint_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.rate_limits
    ADD CONSTRAINT rate_limits_endpoint_id_fkey FOREIGN KEY (endpoint_id) REFERENCES public.api_endpoints(id) ON DELETE CASCADE;


--
-- Name: receptor_families receptor_families_category_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.receptor_families
    ADD CONSTRAINT receptor_families_category_id_fkey FOREIGN KEY (category_id) REFERENCES public.receptor_family_categories(id);


--
-- Name: receptor_subtypes receptor_subtypes_family_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.receptor_subtypes
    ADD CONSTRAINT receptor_subtypes_family_id_fkey FOREIGN KEY (family_id) REFERENCES public.receptor_families(id) ON DELETE CASCADE;


--
-- Name: research_findings research_findings_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.research_findings
    ADD CONSTRAINT research_findings_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: safety_documents safety_documents_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.safety_documents
    ADD CONSTRAINT safety_documents_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: sar_analysis sar_analysis_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.sar_analysis
    ADD CONSTRAINT sar_analysis_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: scaling_analysis scaling_analysis_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.scaling_analysis
    ADD CONSTRAINT scaling_analysis_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: scaling_analysis scaling_analysis_phase_transition_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.scaling_analysis
    ADD CONSTRAINT scaling_analysis_phase_transition_id_fkey FOREIGN KEY (phase_transition_id) REFERENCES public.phase_transitions(id);


--
-- Name: service_status service_status_service_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.service_status
    ADD CONSTRAINT service_status_service_id_fkey FOREIGN KEY (service_id) REFERENCES public.external_services(id) ON DELETE CASCADE;


--
-- Name: social_comments social_comments_post_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_comments
    ADD CONSTRAINT social_comments_post_id_fkey FOREIGN KEY (post_id) REFERENCES public.social_posts(id) ON DELETE CASCADE;


--
-- Name: social_compound_combinations social_compound_combinations_compound_id_1_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_compound_combinations
    ADD CONSTRAINT social_compound_combinations_compound_id_1_fkey FOREIGN KEY (compound_id_1) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_compound_combinations social_compound_combinations_compound_id_2_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_compound_combinations
    ADD CONSTRAINT social_compound_combinations_compound_id_2_fkey FOREIGN KEY (compound_id_2) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_compound_mentions social_compound_mentions_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_compound_mentions
    ADD CONSTRAINT social_compound_mentions_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_compound_trends social_compound_trends_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_compound_trends
    ADD CONSTRAINT social_compound_trends_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_content_quality social_content_quality_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_content_quality
    ADD CONSTRAINT social_content_quality_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_dosage_data social_dosage_data_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_dosage_data
    ADD CONSTRAINT social_dosage_data_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_effect_reports social_effect_reports_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_effect_reports
    ADD CONSTRAINT social_effect_reports_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_effect_reports social_effect_reports_effect_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_effect_reports
    ADD CONSTRAINT social_effect_reports_effect_id_fkey FOREIGN KEY (effect_id) REFERENCES public.social_effects(id) ON DELETE CASCADE;


--
-- Name: social_experience_reports social_experience_reports_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_experience_reports
    ADD CONSTRAINT social_experience_reports_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_experience_reports social_experience_reports_post_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_experience_reports
    ADD CONSTRAINT social_experience_reports_post_id_fkey FOREIGN KEY (post_id) REFERENCES public.social_posts(id) ON DELETE CASCADE;


--
-- Name: social_experience_reports social_experience_reports_protocol_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_experience_reports
    ADD CONSTRAINT social_experience_reports_protocol_id_fkey FOREIGN KEY (protocol_id) REFERENCES public.social_protocols(id) ON DELETE SET NULL;


--
-- Name: social_experience_reports social_experience_reports_stack_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_experience_reports
    ADD CONSTRAINT social_experience_reports_stack_id_fkey FOREIGN KEY (stack_id) REFERENCES public.social_stacks(id) ON DELETE SET NULL;


--
-- Name: social_harm_reduction social_harm_reduction_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_harm_reduction
    ADD CONSTRAINT social_harm_reduction_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_harm_reduction social_harm_reduction_post_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_harm_reduction
    ADD CONSTRAINT social_harm_reduction_post_id_fkey FOREIGN KEY (post_id) REFERENCES public.social_posts(id) ON DELETE CASCADE;


--
-- Name: social_influence_networks social_influence_networks_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_influence_networks
    ADD CONSTRAINT social_influence_networks_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_platform_stats social_platform_stats_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_platform_stats
    ADD CONSTRAINT social_platform_stats_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_posts social_posts_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_posts
    ADD CONSTRAINT social_posts_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_research_reviews social_research_reviews_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_research_reviews
    ADD CONSTRAINT social_research_reviews_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_safety_incidents social_safety_incidents_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_safety_incidents
    ADD CONSTRAINT social_safety_incidents_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_scientific_content social_scientific_content_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_scientific_content
    ADD CONSTRAINT social_scientific_content_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: social_scientific_content social_scientific_content_post_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_scientific_content
    ADD CONSTRAINT social_scientific_content_post_id_fkey FOREIGN KEY (post_id) REFERENCES public.social_posts(id) ON DELETE CASCADE;


--
-- Name: social_substance_data social_substance_data_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.social_substance_data
    ADD CONSTRAINT social_substance_data_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: storage_requirements storage_requirements_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.storage_requirements
    ADD CONSTRAINT storage_requirements_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: structure_similarity structure_similarity_compound_id_1_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.structure_similarity
    ADD CONSTRAINT structure_similarity_compound_id_1_fkey FOREIGN KEY (compound_id_1) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: structure_similarity structure_similarity_compound_id_2_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.structure_similarity
    ADD CONSTRAINT structure_similarity_compound_id_2_fkey FOREIGN KEY (compound_id_2) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: subjective_effect_categories subjective_effect_categories_parent_category_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.subjective_effect_categories
    ADD CONSTRAINT subjective_effect_categories_parent_category_id_fkey FOREIGN KEY (parent_category_id) REFERENCES public.subjective_effect_categories(id);


--
-- Name: subjective_effects subjective_effects_category_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.subjective_effects
    ADD CONSTRAINT subjective_effects_category_id_fkey FOREIGN KEY (category_id) REFERENCES public.subjective_effect_categories(id);


--
-- Name: systems_biology_analysis systems_biology_analysis_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.systems_biology_analysis
    ADD CONSTRAINT systems_biology_analysis_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: therapeutic_class_categories therapeutic_class_categories_parent_category_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.therapeutic_class_categories
    ADD CONSTRAINT therapeutic_class_categories_parent_category_id_fkey FOREIGN KEY (parent_category_id) REFERENCES public.therapeutic_class_categories(id);


--
-- Name: therapeutic_classes therapeutic_classes_category_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.therapeutic_classes
    ADD CONSTRAINT therapeutic_classes_category_id_fkey FOREIGN KEY (category_id) REFERENCES public.therapeutic_class_categories(id);


--
-- Name: toxicity_assays toxicity_assays_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.toxicity_assays
    ADD CONSTRAINT toxicity_assays_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: toxicity_endpoints toxicity_endpoints_category_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.toxicity_endpoints
    ADD CONSTRAINT toxicity_endpoints_category_id_fkey FOREIGN KEY (category_id) REFERENCES public.toxicity_endpoint_categories(id);


--
-- Name: toxicity_mechanism_categories toxicity_mechanism_categories_parent_category_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.toxicity_mechanism_categories
    ADD CONSTRAINT toxicity_mechanism_categories_parent_category_id_fkey FOREIGN KEY (parent_category_id) REFERENCES public.toxicity_mechanism_categories(id);


--
-- Name: toxicity_mechanisms toxicity_mechanisms_category_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.toxicity_mechanisms
    ADD CONSTRAINT toxicity_mechanisms_category_id_fkey FOREIGN KEY (category_id) REFERENCES public.toxicity_mechanism_categories(id);


--
-- Name: training_history training_history_model_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.training_history
    ADD CONSTRAINT training_history_model_id_fkey FOREIGN KEY (model_id) REFERENCES public.ml_models(id) ON DELETE CASCADE;


--
-- Name: training_history training_history_version_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.training_history
    ADD CONSTRAINT training_history_version_id_fkey FOREIGN KEY (version_id) REFERENCES public.model_versions(id);


--
-- Name: wavefunction_analysis wavefunction_analysis_compound_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.wavefunction_analysis
    ADD CONSTRAINT wavefunction_analysis_compound_id_fkey FOREIGN KEY (compound_id) REFERENCES public.compounds(id) ON DELETE CASCADE;


--
-- Name: web_data_sources web_data_sources_category_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.web_data_sources
    ADD CONSTRAINT web_data_sources_category_id_fkey FOREIGN KEY (category_id) REFERENCES public.web_source_categories(id);


--
-- Name: web_hook_deliveries web_hook_deliveries_hook_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: armand
--

ALTER TABLE ONLY public.web_hook_deliveries
    ADD CONSTRAINT web_hook_deliveries_hook_id_fkey FOREIGN KEY (hook_id) REFERENCES public.web_hooks(id) ON DELETE CASCADE;


--
-- PostgreSQL database dump complete
--

