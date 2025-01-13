# Chemdata Database Schema

## Overview
This schema defines the database structure for the Chemdata project, organizing chemical compound data, analysis results, and social media monitoring capabilities. The schema is designed to be modular, maintainable, and ensures data integrity through comprehensive validation and verification mechanisms.

## Directory Structure
```
schema/
├── 00_main.sql                    # Main schema file that handles initialization and proper ordering
├── core/                          # Core database components
│   ├── 00_extensions.sql         # PostgreSQL extensions
│   ├── 01_functions.sql          # Common database functions
│   ├── 02_data_quality.sql       # Data quality tracking and validation
│   ├── 03_social_validation.sql  # Social media content validation
│   └── 04_schema_verification.sql # Schema integrity verification
├── compounds/                     # Compound-related tables
│   ├── 01_base.sql              # Base compound tables
│   └── 02_binding.sql           # Binding data tables
├── receptors/                    # Receptor-related tables
│   ├── 01_base.sql             # Base receptor tables
│   └── 02_subtypes.sql         # Receptor subtype tables
├── reference_data/              # Reference and lookup tables
│   ├── 01_toxicity_endpoints.sql
│   ├── 02_receptor_families.sql
│   ├── 03_therapeutic_classes.sql
│   ├── 04_subjective_effects.sql
│   ├── 05_web_sources.sql
│   ├── 06_toxicity_mechanisms.sql
│   ├── 07_organ_toxicity.sql
│   ├── 08_safety_thresholds.sql
│   ├── 09_monitoring_parameters.sql
│   └── 10_intervention_thresholds.sql
├── safety/                      # Safety and toxicity tables
│   ├── 01_toxicity_data.sql
│   └── 02_risk_assessment.sql
├── pharmacology/                # Pharmacological data tables
│   ├── 01_pharmacology_data.sql
│   └── 02_profiles.sql
├── clinical/                    # Clinical data tables
│   ├── 01_clinical_data.sql
│   └── 02_experience.sql
├── regulatory/                  # Regulatory information tables
│   └── 01_regulatory_data.sql
├── social/                      # Social media and community data tables
│   ├── 01_community_data.sql
│   ├── 02_media_data.sql
│   ├── 03_community_platforms.sql
│   ├── 04_additional_forums.sql
│   ├── 05_more_forums.sql
│   ├── 06_erowid.sql
│   ├── 07_psychonautwiki.sql
│   ├── 08_tripsit.sql
│   ├── 09_longecity.sql
│   ├── 10_reddit.sql
│   ├── 11_twitter.sql
│   ├── 12_cross_platform_analytics.sql
│   └── 13_alerts.sql
├── literature/                  # Literature and documentation tables
│   ├── 01_literature_data.sql
│   └── 02_analysis.sql
├── analysis/                    # Analysis and machine learning tables
│   ├── 01_literature.sql
│   ├── 02_machine_learning.sql
│   └── 03_sar_data.sql
├── web/                        # Web interface related tables
│   ├── web_templates.sql
│   ├── web_settings.sql
│   └── 03_data_integration.sql
├── ml/                         # Machine learning model tables
│   └── models.sql
└── monitoring/                 # System monitoring tables
    └── alerts.sql
```

## Core Components

### 1. Extensions and Base Setup (00_extensions.sql)
- Required PostgreSQL extensions
- Text search configurations
- Array and JSON support
- Base data types and domains

### 2. Common Functions (01_functions.sql)
- Utility functions
  * `format_chemical_name(text)`: Standardize chemical nomenclature
  * `calculate_molecular_weight(text)`: Calculate from SMILES
  * `validate_smiles(text)`: Verify SMILES notation
  * `standardize_units(text, text)`: Convert between unit systems
- Trigger functions
  * `audit_trigger_func()`: Track all data changes
  * `update_timestamp_func()`: Maintain timestamps
  * `validate_before_insert_func()`: Pre-insert validation
  * `validate_before_update_func()`: Pre-update validation
- Validation functions
  * `validate_structure(text)`: Validate chemical structure
  * `validate_binding_data(json)`: Validate binding measurements
  * `validate_toxicity_data(json)`: Validate toxicity data
  * `validate_clinical_data(json)`: Validate clinical results
- Data processing functions
  * `process_raw_binding_data(json)`: Process binding data
  * `calculate_statistics(float[])`: Statistical calculations
  * `generate_fingerprints(text)`: Generate molecular fingerprints
  * `calculate_similarity(text, text)`: Calculate structural similarity

### 3. Data Quality (02_data_quality.sql)
- Data quality metrics
- Validation rules
- Quality tracking
- Monitoring functions

### 4. Social Validation (03_social_validation.sql)
- Content validation rules
- Platform-specific validation
- Community guidelines enforcement
- Content quality metrics

### 5. Schema Verification (04_schema_verification.sql)
- Schema integrity checks
- Index verification
- Constraint validation
- Dependency validation

## Key Tables

### Core Data
- `compounds`: Base compound information
- `receptor_families`: Receptor classification
- `binding_data`: Receptor binding measurements
- `genes`: Genetic information
- `proteins`: Protein structures and properties

### Reference Data
- `toxicity_endpoints`: Standard toxicity measurements
- `therapeutic_classes`: Drug classification system
- `safety_thresholds`: Safety monitoring thresholds
- `monitoring_parameters`: Parameters to track
- `intervention_thresholds`: Clinical intervention criteria

### Social Media Components
- Community data tracking
  * User activity monitoring
  * Content engagement metrics
  * Reputation scoring
  * Trust level assessment
- Platform-specific data collection
  * Reddit: Subreddit monitoring, post/comment tracking
  * Twitter: Tweet analysis, hashtag tracking
  * Bluelight: Forum discussions, experience reports
  * PsychonautWiki: Documentation and guidelines
  * Erowid: Experience reports and safety information
- Cross-platform analytics
  * Sentiment analysis
  * Topic modeling
  * Trend detection
  * Pattern recognition
  * Anomaly detection
- Real-time monitoring and alerts
  * Safety concerns detection
  * Misinformation tracking
  * Abuse pattern detection
  * Emergency response triggers
- Content moderation and quality control
  * Automated content filtering
  * Manual review queues
  * Quality scoring
  * Source reliability assessment
  * Content verification workflows

## Common Patterns

### Timestamps and Auditing
All tables include:
- `created_at`: Creation timestamp
- `updated_at`: Last modification timestamp
- `created_by`: User who created the record
- `updated_by`: User who last modified the record

### Data Quality Tracking
- Validation status
- Quality metrics
- Confidence scores
- Verification flags

### JSON/JSONB Fields
Used for:
- Complex structured data
- Flexible property storage
- Analytics results
- Configuration data

### Indexing Strategy
- B-tree indexes on foreign keys
- GIN indexes on JSONB/array fields
- Text search indexes on content
- Composite indexes for common queries

## Schema Initialization

The schema is initialized through `00_main.sql`, which:
1. Creates required extensions
2. Loads core components in dependency order
3. Creates base tables
4. Loads reference data
5. Creates indexes and constraints
6. Verifies schema integrity

## Validation and Verification

### Data Quality
- Input validation rules
- Data quality metrics
- Content validation
- Cross-reference checks

### Schema Integrity
- Table existence checks
  ```sql
  -- Check if required tables exist
  SELECT tablename 
  FROM pg_tables 
  WHERE schemaname = 'public' 
  AND tablename IN ('compounds', 'binding_data', 'toxicity_data');
  ```
- Foreign key validation
  ```sql
  -- Verify foreign key constraints
  SELECT conname, conrelid::regclass, confrelid::regclass
  FROM pg_constraint
  WHERE contype = 'f'
  AND connamespace = 'public'::regnamespace;
  ```
- Index verification
  ```sql
  -- Check required indexes
  SELECT schemaname, tablename, indexname, indexdef
  FROM pg_indexes
  WHERE schemaname = 'public'
  AND tablename IN ('compounds', 'binding_data');
  ```
- Constraint validation
  ```sql
  -- Verify table constraints
  SELECT conname, contype, pg_get_constraintdef(oid)
  FROM pg_constraint
  WHERE connamespace = 'public'::regnamespace;
  ```

### Performance Monitoring
- Query performance tracking
- Index usage statistics
- Table statistics
- Resource utilization

## Usage

### Installation
```bash
createdb chemdata
psql -d chemdata -f database/schema/00_main.sql
```

### Verification
```sql
-- Check core tables
SELECT tablename FROM pg_tables WHERE schemaname = 'public';

-- Verify extensions
SELECT extname FROM pg_extension;

-- Check triggers
SELECT tgname FROM pg_trigger;
```

### Maintenance
```sql
-- Analyze tables
ANALYZE verbose;

-- Check indexes
SELECT schemaname, tablename, indexname FROM pg_indexes;

-- View audit logs
SELECT * FROM audit_log ORDER BY created_at DESC LIMIT 10;
```

## Best Practices

### Database Design
1. Use explicit foreign key constraints
2. Create appropriate indexes
3. Include table and column comments
4. Maintain proper dependency ordering
5. Use consistent naming conventions
6. Include audit triggers
7. Validate data integrity
8. Implement proper security controls

### Error Handling and Recovery

1. Chemical Data Validation
   ```sql
   -- Validate chemical structure before insert
   CREATE OR REPLACE FUNCTION validate_chemical_structure()
   RETURNS TRIGGER AS $$
   BEGIN
     -- Validate SMILES notation
     IF NOT is_valid_smiles(NEW.smiles) THEN
       RAISE EXCEPTION 'Invalid SMILES notation: %', NEW.smiles
         USING HINT = 'Check chemical structure format';
     END IF;
     
     -- Validate molecular weight
     IF NEW.molecular_weight <= 0 OR NEW.molecular_weight > 2000 THEN
       RAISE EXCEPTION 'Invalid molecular weight: %', NEW.molecular_weight
         USING HINT = 'Weight must be between 0 and 2000';
     END IF;
     
     RETURN NEW;
   END;
   $$ LANGUAGE plpgsql;
   ```

2. Transaction Management
   ```sql
   -- Example: Batch compound import with validation
   DO $$
   BEGIN
     -- Start transaction
     BEGIN;
     
     -- Create savepoint before import
     SAVEPOINT before_import;
     
     -- Import compounds
     INSERT INTO compounds (name, smiles, molecular_weight)
     SELECT name, smiles, calculate_molecular_weight(smiles)
     FROM import_staging_table;
     
     -- Validate imported data
     PERFORM validate_imported_compounds();
     
     -- Calculate additional properties
     UPDATE compounds
     SET logp = calculate_logp(smiles),
         tpsa = calculate_tpsa(smiles)
     WHERE created_at >= current_timestamp - interval '1 hour';
     
     COMMIT;
   EXCEPTION WHEN OTHERS THEN
     -- Log error details
     INSERT INTO error_log (
       operation_type,
       error_code,
       error_message,
       affected_records,
       stack_trace
     ) VALUES (
       'compound_import',
       SQLSTATE,
       SQLERRM,
       (SELECT count(*) FROM import_staging_table),
       format('%s %s %s', SQLSTATE, SQLERRM, PG_CONTEXT)
     );
     
     -- Rollback to savepoint
     ROLLBACK TO before_import;
     
     -- Raise error to caller
     RAISE EXCEPTION 'Import failed: %', SQLERRM;
   END;
   $$;
   ```

3. Data Recovery and Backup
   ```sql
   -- Create backup with chemical structure validation
   CREATE OR REPLACE PROCEDURE backup_chemical_database(
     backup_path text,
     validate boolean DEFAULT true
   ) AS $$
   DECLARE
     validation_errors text[];
   BEGIN
     -- Validate structures if requested
     IF validate THEN
       SELECT array_agg(id || ': ' || error_message)
       INTO validation_errors
       FROM validate_all_structures();
       
       IF array_length(validation_errors, 1) > 0 THEN
         RAISE EXCEPTION 'Structure validation failed: %',
           array_to_string(validation_errors, E'\n');
       END IF;
     END IF;
     
     -- Create backup point
     PERFORM pg_create_restore_point(
       'chemdata_backup_' || to_char(current_timestamp, 'YYYY_MM_DD_HH24_MI_SS')
     );
     
     -- Perform backup
     PERFORM pg_dump_database(
       'chemdata',
       backup_path,
       array['--clean', '--if-exists']
     );
     
     -- Log backup
     INSERT INTO backup_log (
       backup_path,
       validated,
       record_counts
     ) VALUES (
       backup_path,
       validate,
       jsonb_build_object(
         'compounds', (SELECT count(*) FROM compounds),
         'binding_data', (SELECT count(*) FROM binding_data),
         'toxicity_data', (SELECT count(*) FROM toxicity_data)
       )
     );
   END;
   $$ LANGUAGE plpgsql;
   ```

4. Error Monitoring and Alerts
   ```sql
   -- Monitor data quality issues
   CREATE OR REPLACE FUNCTION monitor_data_quality()
   RETURNS TABLE (
     issue_type text,
     severity text,
     record_count bigint,
     details jsonb
   ) AS $$
   BEGIN
     RETURN QUERY
     
     -- Check for invalid structures
     SELECT 
       'invalid_structure' as issue_type,
       'high' as severity,
       count(*),
       jsonb_agg(jsonb_build_object(
         'id', id,
         'name', name,
         'smiles', smiles,
         'error', error_message
       ))
     FROM validate_chemical_structures()
     WHERE is_valid = false
     
     UNION ALL
     
     -- Check for missing binding data
     SELECT
       'missing_binding_data',
       'medium',
       count(*),
       jsonb_agg(jsonb_build_object(
         'compound_id', c.id,
         'name', c.name
       ))
     FROM compounds c
     LEFT JOIN binding_data b ON c.id = b.compound_id
     WHERE b.id IS NULL
     
     UNION ALL
     
     -- Check for outlier values
     SELECT
       'property_outliers',
       'low',
       count(*),
       jsonb_agg(jsonb_build_object(
         'compound_id', id,
         'name', name,
         'property', outlier_property,
         'value', outlier_value
       ))
     FROM detect_property_outliers();
   END;
   $$ LANGUAGE plpgsql;
   ```

### Data Quality Metrics
1. Completeness Checks
   ```sql
   -- Monitor data completeness
   SELECT table_name, 
          count(*) as total_rows,
          sum(case when critical_fields_complete then 1 else 0 end) as complete_rows,
          sum(case when critical_fields_complete then 1 else 0 end)::float / count(*) as completion_rate
   FROM data_quality_metrics
   GROUP BY table_name;
   ```

2. Accuracy Metrics
   ```sql
   -- Track data accuracy
   SELECT metric_name,
          avg(accuracy_score) as avg_accuracy,
          min(accuracy_score) as min_accuracy,
          max(accuracy_score) as max_accuracy
   FROM quality_metrics
   WHERE metric_type = 'accuracy'
   GROUP BY metric_name;
   ```

3. Validation Thresholds
   ```sql
   -- Define quality thresholds
   INSERT INTO quality_thresholds (
     metric_name, min_value, max_value, severity
   ) VALUES
     ('completeness_rate', 0.95, 1.0, 'high'),
     ('accuracy_score', 0.90, 1.0, 'high'),
     ('validation_rate', 0.98, 1.0, 'critical');
   ```

### Security Controls
1. Access Control
   ```sql
   -- Role-based access
   CREATE ROLE chemdata_reader;
   CREATE ROLE chemdata_writer;
   CREATE ROLE chemdata_admin;
   
   -- Grant permissions
   GRANT SELECT ON ALL TABLES IN SCHEMA public TO chemdata_reader;
   GRANT INSERT, UPDATE ON allowed_tables TO chemdata_writer;
   GRANT ALL ON ALL TABLES IN SCHEMA public TO chemdata_admin;
   ```

2. Data Encryption
   ```sql
   -- Enable encryption
   ALTER TABLE sensitive_data 
   ALTER COLUMN confidential_info 
   SET DATA TYPE bytea 
   USING encrypt_sensitive_data(confidential_info::text);
   ```

3. Audit Logging
   ```sql
   -- Enhanced audit logging
   CREATE TRIGGER enhanced_audit_trigger
   AFTER INSERT OR UPDATE OR DELETE ON sensitive_tables
   FOR EACH ROW EXECUTE FUNCTION log_data_changes();
   ```

4. Row-Level Security
   ```sql
   -- Enable RLS
   ALTER TABLE compounds ENABLE ROW LEVEL SECURITY;
   
   -- Create access policy
   CREATE POLICY compound_access_policy ON compounds
   FOR ALL
   TO chemdata_users
   USING (created_by = current_user OR is_public = true);
   ```

### Data Types
- Use `uuid` for primary keys
- Use `timestamptz` for timestamps
- Use `jsonb` for complex data
- Use `text[]` for arrays
- Use `double precision` for measurements

### Social Media Data
1. Implement rate limiting
2. Handle API quotas
3. Validate data quality
4. Monitor content appropriateness
5. Track data provenance
6. Implement privacy controls

## Maintenance

### Regular Tasks
1. Run schema verification
   ```sql
   -- Verify schema integrity
   SELECT verify_schema_integrity();
   
   -- Check for missing indexes
   SELECT * FROM verify_required_indexes();
   
   -- Validate foreign keys
   SELECT verify_foreign_keys();
   ```
2. Check index usage
   ```sql
   -- Monitor index usage
   SELECT schemaname, tablename, indexname, idx_scan, idx_tup_read
   FROM pg_stat_user_indexes
   ORDER BY idx_scan DESC;
   ```
3. Verify data quality metrics
   ```sql
   -- Check data quality scores
   SELECT * FROM data_quality_metrics
   WHERE score < threshold
   ORDER BY severity DESC;
   ```
4. Update validation rules
   ```sql
   -- Review and update rules
   SELECT * FROM validation_rules
   WHERE last_updated < current_date - interval '30 days';
   ```
5. Monitor performance
   ```sql
   -- Check query performance
   SELECT * FROM pg_stat_statements
   ORDER BY total_time DESC
   LIMIT 10;
   ```
6. Review audit logs
   ```sql
   -- Check recent changes
   SELECT * FROM audit_log
   WHERE created_at > current_timestamp - interval '24 hours'
   ORDER BY created_at DESC;
   ```
7. Update documentation
   - Review and update README
   - Document new features
   - Update maintenance procedures
   - Track schema changes

### Adding New Components
1. Create component in appropriate directory
2. Add to `00_main.sql` in correct order
3. Update schema verification
4. Add required indexes
5. Document changes
6. Test thoroughly

### Modifying Existing Components
1. Maintain referential integrity
2. Update affected constraints
3. Adjust indexes as needed
4. Test changes thoroughly
5. Update documentation
6. Consider backward compatibility

## Documentation
- Keep this README updated
- Document all schema changes
- Maintain change logs
- Update dependency diagrams
- Document best practices
- Track schema versions

## Testing and Continuous Integration

### Automated Schema Tests
```sql
-- Test suite for schema validation
CREATE OR REPLACE FUNCTION test_schema_integrity()
RETURNS SETOF text AS $$
BEGIN
  -- Test required tables exist
  RETURN NEXT ok(
    EXISTS (
      SELECT 1 FROM pg_tables 
      WHERE schemaname = 'public' 
      AND tablename = 'compounds'
    ),
    'compounds table exists'
  );

  -- Test chemical structure validation
  RETURN NEXT lives_ok(
    $$INSERT INTO compounds (name, smiles) 
      VALUES ('Caffeine', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C')$$,
    'Valid SMILES notation accepted'
  );

  RETURN NEXT throws_ok(
    $$INSERT INTO compounds (name, smiles) 
      VALUES ('Invalid', 'NOT-A-SMILES')$$,
    'Invalid SMILES notation',
    'Invalid chemical structure rejected'
  );

  -- Test foreign key constraints
  RETURN NEXT throws_ok(
    $$INSERT INTO binding_data (compound_id, receptor_id, ki_value) 
      VALUES (999999, 1, 1.0)$$,
    '23503',
    'Foreign key violation caught'
  );
END;
$$ LANGUAGE plpgsql;

-- Run as part of CI pipeline
SELECT * FROM runtests();
```

### Integration Tests
```sql
-- Test data import pipeline
CREATE OR REPLACE PROCEDURE test_compound_import()
LANGUAGE plpgsql AS $$
DECLARE
  test_data jsonb;
BEGIN
  -- Prepare test data
  test_data := '[
    {"name": "Test Compound", "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"},
    {"name": "Another Test", "smiles": "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"}
  ]'::jsonb;
  
  -- Test full import pipeline
  BEGIN
    -- Stage data
    CREATE TEMP TABLE import_staging (data jsonb);
    INSERT INTO import_staging VALUES (test_data);
    
    -- Run import
    CALL import_compounds_from_staging();
    
    -- Verify results
    ASSERT (
      SELECT count(*) = 2 
      FROM compounds 
      WHERE source = 'test_import'
    ), 'Import failed - wrong number of records';
    
    -- Verify computed properties
    ASSERT (
      SELECT count(*) = 2 
      FROM compounds 
      WHERE molecular_weight IS NOT NULL 
      AND logp IS NOT NULL
    ), 'Property calculation failed';
    
    -- Cleanup
    DROP TABLE import_staging;
  EXCEPTION WHEN OTHERS THEN
    RAISE NOTICE 'Test failed: %', SQLERRM;
    RAISE;
  END;
END;
$$;

-- CI/CD Pipeline Integration
CREATE OR REPLACE PROCEDURE run_ci_tests()
LANGUAGE plpgsql AS $$
BEGIN
  -- Run schema tests
  PERFORM test_schema_integrity();
  
  -- Run data validation tests
  PERFORM test_data_validation();
  
  -- Run integration tests
  CALL test_compound_import();
  
  -- Run performance tests
  PERFORM test_query_performance();
  
  -- Generate test report
  INSERT INTO test_results (
    run_id,
    timestamp,
    tests_passed,
    tests_failed,
    performance_metrics
  ) VALUES (
    gen_random_uuid(),
    current_timestamp,
    (SELECT count(*) FROM test_results WHERE status = 'pass'),
    (SELECT count(*) FROM test_results WHERE status = 'fail'),
    jsonb_build_object(
      'avg_query_time', (SELECT avg(execution_time) FROM query_stats),
      'max_memory_used', (SELECT max(memory_used) FROM resource_stats)
    )
  );
END;
$$;
```

### Performance Benchmarks
```sql
-- Benchmark critical operations
CREATE OR REPLACE FUNCTION benchmark_critical_operations()
RETURNS TABLE (
  operation_name text,
  avg_duration interval,
  p95_duration interval,
  success_rate float
) AS $$
BEGIN
  RETURN QUERY
  WITH benchmarks AS (
    -- Benchmark structure validation
    SELECT
      'structure_validation' as op,
      clock_timestamp() - start_time as duration,
      success
    FROM (
      SELECT clock_timestamp() as start_time,
      validate_chemical_structure(smiles) as success
      FROM compounds
      LIMIT 1000
    ) v
    
    UNION ALL
    
    -- Benchmark similarity search
    SELECT
      'similarity_search' as op,
      clock_timestamp() - start_time as duration,
      true as success
    FROM (
      SELECT 
        clock_timestamp() as start_time,
        find_similar_compounds(smiles, 0.7)
      FROM compounds
      LIMIT 100
    ) s
  )
  SELECT 
    op as operation_name,
    avg(duration) as avg_duration,
    percentile_cont(0.95) 
      WITHIN GROUP (ORDER BY duration) as p95_duration,
    sum(case when success then 1 else 0 end)::float / count(*) as success_rate
  FROM benchmarks
  GROUP BY op;
END;
$$ LANGUAGE plpgsql;
```

## High Availability and Disaster Recovery

### Replication Setup
```sql
-- Configure streaming replication
ALTER SYSTEM SET wal_level = replica;
ALTER SYSTEM SET max_wal_senders = 10;
ALTER SYSTEM SET max_replication_slots = 10;

-- Create replication slot
SELECT * FROM pg_create_physical_replication_slot('chemdata_replica');

-- Monitor replication status
SELECT slot_name, active, restart_lsn
FROM pg_replication_slots;

SELECT client_addr, state, sent_lsn, write_lsn, flush_lsn
FROM pg_stat_replication;
```

### Backup Strategies
1. Physical Backups
   ```bash
   # Full backup with validation
   pg_basebackup -D /backup/chemdata/base -Ft -Xs -P -U replicator
   
   # Validate backup
   pg_verifybackup /backup/chemdata/base
   ```

2. Logical Backups
   ```bash
   # Schema-only backup
   pg_dump -s chemdata > schema_backup.sql
   
   # Full database backup with compression
   pg_dump -Fc chemdata > chemdata_full.dump
   
   # Custom-format backup with parallel jobs
   pg_dump -Fd chemdata -j 4 -f backup/chemdata
   ```

3. Point-in-Time Recovery
   ```sql
   -- Create recovery target
   SELECT pg_create_restore_point('before_major_update');
   
   -- Configure recovery
   ALTER SYSTEM SET restore_command = 'cp /path/to/archive/%f %p';
   ALTER SYSTEM SET recovery_target_name = 'before_major_update';
   ALTER SYSTEM SET recovery_target_action = 'promote';
   ```

### Monitoring and Alerts
```sql
-- Monitor replication lag
CREATE OR REPLACE FUNCTION check_replication_lag()
RETURNS TABLE (
  replica_host text,
  lag_bytes bigint,
  lag_time interval
) AS $$
BEGIN
  RETURN QUERY
  SELECT 
    client_addr::text,
    pg_wal_lsn_diff(pg_current_wal_lsn(), replay_lsn),
    clock_timestamp() - reply_time
  FROM pg_stat_replication;
END;
$$ LANGUAGE plpgsql;

-- Set up lag alerts
CREATE OR REPLACE FUNCTION alert_on_replication_lag()
RETURNS trigger AS $$
BEGIN
  IF NEW.lag_bytes > 1024*1024*100 OR  -- 100MB
     NEW.lag_time > interval '5 minutes' THEN
    
    INSERT INTO alert_log (
      alert_type,
      severity,
      message
    ) VALUES (
      'replication_lag',
      'high',
      format(
        'Replica %s is lagging: %s bytes, %s time',
        NEW.replica_host,
        NEW.lag_bytes,
        NEW.lag_time
      )
    );
  END IF;
  RETURN NEW;
END;
$$ LANGUAGE plpgsql;
```

### Failover Procedures
```sql
-- Check if failover is needed
CREATE OR REPLACE FUNCTION check_failover_conditions()
RETURNS boolean AS $$
DECLARE
  primary_healthy boolean;
  replica_ready boolean;
BEGIN
  -- Check primary health
  SELECT pg_is_in_recovery() = false
  INTO primary_healthy
  FROM pg_stat_database
  WHERE datname = 'chemdata'
  LIMIT 1;
  
  -- Check replica readiness
  SELECT 
    pg_wal_lsn_diff(pg_last_wal_receive_lsn(), pg_last_wal_replay_lsn()) = 0
    AND pg_is_in_recovery() = true
  INTO replica_ready
  FROM pg_stat_wal_receiver;
  
  RETURN NOT primary_healthy AND replica_ready;
END;
$$ LANGUAGE plpgsql;

-- Perform automated failover
CREATE OR REPLACE PROCEDURE execute_failover()
LANGUAGE plpgsql AS $$
BEGIN
  -- Log failover start
  INSERT INTO failover_log (
    timestamp,
    initiated_by,
    reason
  ) VALUES (
    current_timestamp,
    current_user,
    'Primary server unresponsive'
  );
  
  -- Promote replica
  PERFORM pg_promote(wait => true, wait_seconds => 300);
  
  -- Update connection info
  UPDATE connection_endpoints 
  SET is_primary = true 
  WHERE host = current_setting('server_name');
  
  -- Notify applications
  NOTIFY 'failover_complete', '';
  
  -- Log completion
  UPDATE failover_log 
  SET completed_at = current_timestamp,
      status = 'success'
  WHERE id = currval('failover_log_id_seq');
END;
$$;
```
