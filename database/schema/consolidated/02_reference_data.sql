-- Consolidated reference data schema
-- Combines all reference tables used across the system

-- Proteins and genes
CREATE TABLE IF NOT EXISTS genes (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL,
    symbol text NOT NULL UNIQUE,
    description text,
    organism text NOT NULL,
    chromosome text,
    location text,
    gene_type text,
    sequence text,
    -- Enhanced fields
    alternative_symbols text[],
    ensembl_id text UNIQUE,
    entrez_id text UNIQUE,
    hgnc_id text UNIQUE,
    mgi_id text,  -- For mouse genes
    rgd_id text,  -- For rat genes
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
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    CONSTRAINT valid_gene_symbol CHECK (symbol ~ '^[A-Za-z0-9-]+$')
);

CREATE TABLE IF NOT EXISTS proteins (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    gene_id uuid REFERENCES genes(id),
    name text NOT NULL,
    symbol text NOT NULL,
    description text,
    organism text NOT NULL,
    sequence text,
    length integer,
    molecular_weight double precision,
    -- Enhanced fields
    uniprot_id text UNIQUE,
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
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Create indexes for genes and proteins
CREATE INDEX IF NOT EXISTS idx_genes_symbol ON genes(symbol);
CREATE INDEX IF NOT EXISTS idx_genes_ensembl ON genes(ensembl_id);
CREATE INDEX IF NOT EXISTS idx_genes_entrez ON genes(entrez_id);
CREATE INDEX IF NOT EXISTS idx_genes_hgnc ON genes(hgnc_id);
CREATE INDEX IF NOT EXISTS idx_genes_organism ON genes(organism);

CREATE INDEX IF NOT EXISTS idx_proteins_symbol ON proteins(symbol);
CREATE INDEX IF NOT EXISTS idx_proteins_uniprot ON proteins(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_proteins_gene ON proteins(gene_id);
CREATE INDEX IF NOT EXISTS idx_proteins_organism ON proteins(organism);
CREATE INDEX IF NOT EXISTS idx_proteins_family ON proteins(protein_family);

-- Add triggers for genes and proteins
CREATE TRIGGER update_genes_modtime
    BEFORE UPDATE ON genes
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_proteins_modtime
    BEFORE UPDATE ON proteins
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER audit_genes_trigger
    AFTER INSERT OR UPDATE OR DELETE ON genes
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_proteins_trigger
    AFTER INSERT OR UPDATE OR DELETE ON proteins
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Add comments
COMMENT ON TABLE genes IS 'Gene information and annotations';
COMMENT ON TABLE proteins IS 'Protein structures and annotations';

-- Receptor classification and properties
CREATE TABLE IF NOT EXISTS receptor_family_categories (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    description text,
    receptor_type text NOT NULL,
    signaling_mechanism text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

CREATE TABLE IF NOT EXISTS receptor_families (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    category_id uuid NOT NULL REFERENCES receptor_family_categories(id),
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
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(category_id, name)
);

CREATE TABLE IF NOT EXISTS receptor_subtypes (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    family_id uuid NOT NULL REFERENCES receptor_families(id) ON DELETE CASCADE,
    subtype_name text NOT NULL,
    description text,
    protein_sequence text,
    species text,
    expression_pattern jsonb,
    signaling_pathways text[],
    pharmacological_profile jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(family_id, subtype_name)
);

-- Add triggers for receptor tables
CREATE TRIGGER update_receptor_family_categories_modtime
    BEFORE UPDATE ON receptor_family_categories
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_receptor_families_modtime
    BEFORE UPDATE ON receptor_families
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_receptor_subtypes_modtime
    BEFORE UPDATE ON receptor_subtypes
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER audit_receptor_family_categories_trigger
    AFTER INSERT OR UPDATE OR DELETE ON receptor_family_categories
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_receptor_families_trigger
    AFTER INSERT OR UPDATE OR DELETE ON receptor_families
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_receptor_subtypes_trigger
    AFTER INSERT OR UPDATE OR DELETE ON receptor_subtypes
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Therapeutic and pharmacological classifications
CREATE TABLE IF NOT EXISTS therapeutic_class_categories (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    description text,
    level integer NOT NULL, -- Hierarchical level (1=highest, e.g., CNS drugs)
    parent_category_id uuid REFERENCES therapeutic_class_categories(id),
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    CONSTRAINT valid_level CHECK (level > 0)
);

CREATE TABLE IF NOT EXISTS therapeutic_classes (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    category_id uuid NOT NULL REFERENCES therapeutic_class_categories(id),
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
    text_search_vector tsvector GENERATED ALWAYS AS (
        setweight(to_tsvector('english', coalesce(name, '')), 'A') ||
        setweight(to_tsvector('english', coalesce(description, '')), 'B') ||
        setweight(to_tsvector('english', coalesce(mechanism_of_action, '')), 'C')
    ) STORED,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(category_id, name)
);

CREATE TABLE IF NOT EXISTS pharmacological_classes (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    description text,
    mechanism_type text NOT NULL,
    target_systems text[],
    typical_effects text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Add triggers for therapeutic classes
CREATE TRIGGER update_therapeutic_class_categories_modtime
    BEFORE UPDATE ON therapeutic_class_categories
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_therapeutic_classes_modtime
    BEFORE UPDATE ON therapeutic_classes
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_pharmacological_classes_modtime
    BEFORE UPDATE ON pharmacological_classes
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER audit_therapeutic_class_categories_trigger
    AFTER INSERT OR UPDATE OR DELETE ON therapeutic_class_categories
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_therapeutic_classes_trigger
    AFTER INSERT OR UPDATE OR DELETE ON therapeutic_classes
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_pharmacological_classes_trigger
    AFTER INSERT OR UPDATE OR DELETE ON pharmacological_classes
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Toxicity measurement and classification
CREATE TABLE IF NOT EXISTS toxicity_endpoint_categories (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    description text,
    measurement_type text NOT NULL,
    units text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

CREATE TABLE IF NOT EXISTS toxicity_endpoints (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    category_id uuid NOT NULL REFERENCES toxicity_endpoint_categories(id),
    name text NOT NULL,
    description text,
    standard_unit text NOT NULL,
    conversion_factors jsonb, -- For unit conversions
    detection_methods text[],
    validation_criteria jsonb,
    reference_ranges jsonb,
    severity_thresholds jsonb,
    regulatory_limits jsonb,
    notes text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(category_id, name)
);

-- Add triggers for toxicity endpoints
CREATE TRIGGER update_toxicity_endpoint_categories_modtime
    BEFORE UPDATE ON toxicity_endpoint_categories
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_toxicity_endpoints_modtime
    BEFORE UPDATE ON toxicity_endpoints
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER audit_toxicity_endpoint_categories_trigger
    AFTER INSERT OR UPDATE OR DELETE ON toxicity_endpoint_categories
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_toxicity_endpoints_trigger
    AFTER INSERT OR UPDATE OR DELETE ON toxicity_endpoints
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Toxicity mechanism classification
CREATE TABLE IF NOT EXISTS toxicity_mechanism_categories (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    description text,
    level integer NOT NULL, -- Hierarchical level (1=highest)
    mechanism_type text NOT NULL, -- 'molecular', 'cellular', 'systemic', 'organ'
    parent_category_id uuid REFERENCES toxicity_mechanism_categories(id),
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    CONSTRAINT valid_level CHECK (level > 0)
);

CREATE TABLE IF NOT EXISTS toxicity_mechanisms (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    category_id uuid NOT NULL REFERENCES toxicity_mechanism_categories(id),
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
    text_search_vector tsvector GENERATED ALWAYS AS (
        setweight(to_tsvector('english', coalesce(name, '')), 'A') ||
        setweight(to_tsvector('english', coalesce(description, '')), 'B')
    ) STORED,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(category_id, name)
);

-- Add triggers for toxicity mechanisms
CREATE TRIGGER update_toxicity_mechanism_categories_modtime
    BEFORE UPDATE ON toxicity_mechanism_categories
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_toxicity_mechanisms_modtime
    BEFORE UPDATE ON toxicity_mechanisms
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER audit_toxicity_mechanism_categories_trigger
    AFTER INSERT OR UPDATE OR DELETE ON toxicity_mechanism_categories
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_toxicity_mechanisms_trigger
    AFTER INSERT OR UPDATE OR DELETE ON toxicity_mechanisms
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Organ toxicity classification
CREATE TABLE IF NOT EXISTS organ_systems (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    description text,
    major_components text[],
    key_functions text[],
    vulnerability_factors text[],
    assessment_methods text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

CREATE TABLE IF NOT EXISTS organ_toxicity_patterns (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    organ_system_id uuid NOT NULL REFERENCES organ_systems(id),
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
    text_search_vector tsvector GENERATED ALWAYS AS (
        setweight(to_tsvector('english', coalesce(name, '')), 'A') ||
        setweight(to_tsvector('english', coalesce(description, '')), 'B')
    ) STORED,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(organ_system_id, name)
);

-- Add triggers for organ toxicity
CREATE TRIGGER update_organ_systems_modtime
    BEFORE UPDATE ON organ_systems
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_organ_toxicity_patterns_modtime
    BEFORE UPDATE ON organ_toxicity_patterns
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER audit_organ_systems_trigger
    AFTER INSERT OR UPDATE OR DELETE ON organ_systems
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_organ_toxicity_patterns_trigger
    AFTER INSERT OR UPDATE OR DELETE ON organ_toxicity_patterns
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Safety thresholds and monitoring
CREATE TABLE IF NOT EXISTS safety_thresholds (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    parameter_name text NOT NULL,
    threshold_value double precision NOT NULL,
    unit text,
    severity_level text NOT NULL,
    description text,
    intervention_required boolean DEFAULT false,
    validation_method text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(parameter_name, severity_level)
);

-- Insert standard safety thresholds
INSERT INTO safety_thresholds (parameter_name, threshold_value, unit, severity_level, description) VALUES
('ALT', 40, 'U/L', 'warning', 'Upper limit for alanine aminotransferase'),
('AST', 40, 'U/L', 'warning', 'Upper limit for aspartate aminotransferase'),
('Creatinine', 1.2, 'mg/dL', 'warning', 'Upper limit for serum creatinine'),
('QTc', 450, 'ms', 'warning', 'Upper limit for corrected QT interval'),
('Neutrophils', 1500, 'cells/µL', 'warning', 'Lower limit for neutrophil count'),
('Platelets', 150000, 'cells/µL', 'warning', 'Lower limit for platelet count'),
('Heart Rate', 100, 'bpm', 'warning', 'Upper limit for resting heart rate'),
('Blood Pressure', 140, 'mmHg', 'warning', 'Upper limit for systolic blood pressure'),
('Body Temperature', 38.3, '°C', 'warning', 'Upper limit for body temperature'),
('Respiratory Rate', 20, 'breaths/min', 'warning', 'Upper limit for respiratory rate'),
('Glucose', 126, 'mg/dL', 'warning', 'Upper limit for fasting glucose'),
('Total Bilirubin', 1.2, 'mg/dL', 'warning', 'Upper limit for total bilirubin'),
('Albumin', 3.5, 'g/dL', 'warning', 'Lower limit for serum albumin'),
('eGFR', 60, 'mL/min', 'warning', 'Lower limit for estimated glomerular filtration rate'),
('Oxygen Saturation', 95, '%', 'warning', 'Lower limit for oxygen saturation');

-- Add metadata
COMMENT ON TABLE safety_thresholds IS 'Safety threshold values for various clinical and laboratory parameters';
COMMENT ON COLUMN safety_thresholds.parameter_name IS 'Name of the safety parameter being measured';
COMMENT ON COLUMN safety_thresholds.threshold_value IS 'Numerical threshold value';
COMMENT ON COLUMN safety_thresholds.unit IS 'Unit of measurement';
COMMENT ON COLUMN safety_thresholds.severity_level IS 'Severity level when threshold is exceeded';
COMMENT ON COLUMN safety_thresholds.description IS 'Description of the threshold and its significance';

-- Add triggers for safety thresholds
CREATE TRIGGER update_safety_thresholds_modtime
    BEFORE UPDATE ON safety_thresholds
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER audit_safety_thresholds_trigger
    AFTER INSERT OR UPDATE OR DELETE ON safety_thresholds
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Monitoring parameters
CREATE TABLE IF NOT EXISTS monitoring_parameters (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    category text NOT NULL,
    frequency text NOT NULL,
    monitoring_method text NOT NULL,
    alert_conditions text[],
    normal_range jsonb,
    data_type text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Insert standard monitoring parameters
INSERT INTO monitoring_parameters (name, category, frequency, monitoring_method, alert_conditions) VALUES
('Liver Function', 'biochemical', '1 week', 'Blood test', ARRAY['ALT > 3x ULN', 'AST > 3x ULN', 'ALP > 2x ULN']),
('Kidney Function', 'biochemical', '1 week', 'Blood test', ARRAY['Creatinine > 1.5x baseline', 'eGFR decrease > 25%']),
('Cardiac Function', 'physiological', '1 day', 'ECG', ARRAY['QTc > 500ms', 'QTc increase > 60ms']),
('Blood Pressure', 'physiological', '6 hours', 'Automated measurement', ARRAY['Systolic > 160', 'Diastolic > 100']),
('Blood Count', 'hematological', '1 week', 'Blood test', ARRAY['WBC < 3000/µL', 'Platelets < 100k/µL']),
('Mental Status', 'neurological', '6 hours', 'Clinical assessment', ARRAY['Confusion', 'Agitation', 'Drowsiness']),
('Body Temperature', 'physiological', '6 hours', 'Temperature measurement', ARRAY['> 38.5°C', '< 35.5°C']),
('Respiratory Rate', 'physiological', '6 hours', 'Clinical measurement', ARRAY['> 24/min', '< 8/min']),
('Oxygen Saturation', 'physiological', '6 hours', 'Pulse oximetry', ARRAY['< 92%', 'Drop > 4% from baseline']),
('Cognitive Function', 'neurological', '1 day', 'Cognitive assessment', ARRAY['MMSE decrease > 2 points', 'New onset confusion']),
('Sleep Pattern', 'behavioral', '1 day', 'Sleep diary', ARRAY['Insomnia > 2 hours', 'Excessive drowsiness']),
('Appetite', 'behavioral', '1 day', 'Food intake log', ARRAY['Decrease > 50%', 'Complete loss of appetite']),
('Mood', 'psychological', '1 day', 'Mood scale', ARRAY['Severe depression', 'Mania', 'Anxiety']),
('Movement', 'neurological', '1 day', 'Clinical observation', ARRAY['Tremor', 'Ataxia', 'Dystonia']),
('Pain Level', 'subjective', '6 hours', 'Pain scale', ARRAY['Score > 7/10', 'Acute increase > 3 points']);

-- Add metadata
COMMENT ON TABLE monitoring_parameters IS 'Parameters to monitor for safety assessment';
COMMENT ON COLUMN monitoring_parameters.name IS 'Name of the parameter to monitor';
COMMENT ON COLUMN monitoring_parameters.category IS 'Category of the monitoring parameter';
COMMENT ON COLUMN monitoring_parameters.frequency IS 'How often the parameter should be monitored';
COMMENT ON COLUMN monitoring_parameters.monitoring_method IS 'Method used to monitor this parameter';
COMMENT ON COLUMN monitoring_parameters.alert_conditions IS 'Conditions that should trigger alerts';

-- Add triggers for monitoring parameters
CREATE TRIGGER update_monitoring_params_modtime
    BEFORE UPDATE ON monitoring_parameters
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER audit_monitoring_params_trigger
    AFTER INSERT OR UPDATE OR DELETE ON monitoring_parameters
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Subjective effects categorization
-- Subjective effects classification
CREATE TABLE IF NOT EXISTS subjective_effect_categories (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    description text,
    domain text NOT NULL, -- e.g., 'cognitive', 'perceptual', 'emotional'
    level integer NOT NULL, -- Hierarchical level (1=highest)
    parent_category_id uuid REFERENCES subjective_effect_categories(id),
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    CONSTRAINT valid_level CHECK (level > 0)
);

CREATE TABLE IF NOT EXISTS subjective_effects (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    category_id uuid NOT NULL REFERENCES subjective_effect_categories(id),
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
    text_search_vector tsvector GENERATED ALWAYS AS (
        setweight(to_tsvector('english', coalesce(name, '')), 'A') ||
        setweight(to_tsvector('english', coalesce(description, '')), 'B')
    ) STORED,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(category_id, name)
);

CREATE TABLE IF NOT EXISTS effect_relationships (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    effect_id uuid NOT NULL REFERENCES subjective_effects(id),
    related_effect_id uuid NOT NULL REFERENCES subjective_effects(id),
    relationship_type text NOT NULL,
    strength double precision,
    evidence_level text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(effect_id, related_effect_id)
);

-- Research and analysis parameters
CREATE TABLE IF NOT EXISTS quantum_parameters (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    parameter_name text NOT NULL UNIQUE,
    description text,
    unit text,
    calculation_method text NOT NULL,
    typical_range jsonb,
    accuracy_metrics jsonb,
    validation_criteria jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

CREATE TABLE IF NOT EXISTS analysis_parameters (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    parameter_name text NOT NULL UNIQUE,
    description text,
    data_type text NOT NULL,
    validation_rules jsonb,
    default_value jsonb,
    allowed_range jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Web data sources and integration
CREATE TABLE IF NOT EXISTS web_source_categories (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    description text,
    source_type text NOT NULL, -- e.g., 'academic', 'community', 'regulatory'
    reliability_rating integer NOT NULL CHECK (reliability_rating BETWEEN 1 AND 5),
    validation_requirements text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

CREATE TABLE IF NOT EXISTS web_data_sources (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    category_id uuid NOT NULL REFERENCES web_source_categories(id),
    name text NOT NULL,
    base_url text NOT NULL,
    description text,
    api_endpoint text,
    access_method text NOT NULL, -- e.g., 'api', 'scraping', 'manual'
    authentication_type text, -- e.g., 'oauth', 'api_key', 'none'
    rate_limits jsonb,
    data_format text, -- e.g., 'json', 'xml', 'html'
    update_frequency text,
    last_validated timestamptz,
    validation_status text,
    data_quality_metrics jsonb,
    coverage_areas text[],
    known_limitations text[],
    usage_requirements text,
    citation_format text,
    notes text,
    active boolean NOT NULL DEFAULT true,
    text_search_vector tsvector GENERATED ALWAYS AS (
        setweight(to_tsvector('english', coalesce(name, '')), 'A') ||
        setweight(to_tsvector('english', coalesce(description, '')), 'B')
    ) STORED,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(category_id, name)
);

-- Create indexes
-- Create indexes for receptor tables
CREATE INDEX idx_receptor_family_categories_name ON receptor_family_categories(name);
CREATE INDEX idx_receptor_family_categories_type ON receptor_family_categories(receptor_type);

CREATE INDEX idx_receptor_families_category ON receptor_families(category_id);
CREATE INDEX idx_receptor_families_name ON receptor_families(name);
CREATE INDEX idx_receptor_families_type ON receptor_families(protein_type);
CREATE INDEX idx_receptor_families_signaling ON receptor_families(signaling_type);

CREATE INDEX idx_receptor_subtypes_name ON receptor_subtypes(subtype_name);
CREATE INDEX idx_receptor_subtypes_family ON receptor_subtypes(family_id);

-- Create indexes for therapeutic classes
CREATE INDEX idx_therapeutic_class_categories_name ON therapeutic_class_categories(name);
CREATE INDEX idx_therapeutic_class_categories_parent ON therapeutic_class_categories(parent_category_id);
CREATE INDEX idx_therapeutic_class_categories_level ON therapeutic_class_categories(level);

CREATE INDEX idx_therapeutic_classes_category ON therapeutic_classes(category_id);
CREATE INDEX idx_therapeutic_classes_name ON therapeutic_classes(name);
CREATE INDEX idx_therapeutic_classes_text_search ON therapeutic_classes USING gin(text_search_vector);

CREATE INDEX idx_pharmacological_classes_name ON pharmacological_classes(name);
CREATE INDEX idx_pharmacological_classes_type ON pharmacological_classes(mechanism_type);

-- Create indexes for toxicity tables
CREATE INDEX idx_toxicity_endpoint_categories_name ON toxicity_endpoint_categories(name);
CREATE INDEX idx_toxicity_endpoint_categories_type ON toxicity_endpoint_categories(measurement_type);

CREATE INDEX idx_toxicity_endpoints_category ON toxicity_endpoints(category_id);
CREATE INDEX idx_toxicity_endpoints_name ON toxicity_endpoints(name);
CREATE INDEX idx_toxicity_endpoints_unit ON toxicity_endpoints(standard_unit);

-- Create indexes for toxicity mechanisms
CREATE INDEX idx_toxicity_mechanism_categories_name ON toxicity_mechanism_categories(name);
CREATE INDEX idx_toxicity_mechanism_categories_type ON toxicity_mechanism_categories(mechanism_type);
CREATE INDEX idx_toxicity_mechanism_categories_parent ON toxicity_mechanism_categories(parent_category_id);
CREATE INDEX idx_toxicity_mechanism_categories_level ON toxicity_mechanism_categories(level);

CREATE INDEX idx_toxicity_mechanisms_category ON toxicity_mechanisms(category_id);
CREATE INDEX idx_toxicity_mechanisms_name ON toxicity_mechanisms(name);
CREATE INDEX idx_toxicity_mechanisms_text_search ON toxicity_mechanisms USING gin(text_search_vector);

-- Create indexes for organ toxicity
CREATE INDEX idx_organ_systems_name ON organ_systems(name);

CREATE INDEX idx_organ_toxicity_patterns_organ ON organ_toxicity_patterns(organ_system_id);
CREATE INDEX idx_organ_toxicity_patterns_name ON organ_toxicity_patterns(name);
CREATE INDEX idx_organ_toxicity_patterns_text_search ON organ_toxicity_patterns USING gin(text_search_vector);

-- Create indexes for safety thresholds and monitoring
CREATE INDEX idx_safety_thresholds_param ON safety_thresholds(parameter_name);
CREATE INDEX idx_safety_thresholds_severity ON safety_thresholds(severity_level);

CREATE INDEX idx_monitoring_params_name ON monitoring_parameters(name);
CREATE INDEX idx_monitoring_params_category ON monitoring_parameters(category);

-- Create indexes for subjective effects
CREATE INDEX idx_effect_categories_name ON subjective_effect_categories(name);
CREATE INDEX idx_effect_categories_domain ON subjective_effect_categories(domain);
CREATE INDEX idx_effect_categories_level ON subjective_effect_categories(level);
CREATE INDEX idx_effect_categories_parent ON subjective_effect_categories(parent_category_id);

CREATE INDEX idx_subjective_effects_category ON subjective_effects(category_id);
CREATE INDEX idx_subjective_effects_name ON subjective_effects(name);
CREATE INDEX idx_subjective_effects_text_search ON subjective_effects USING gin(text_search_vector);

CREATE INDEX idx_effect_relationships_effect ON effect_relationships(effect_id);
CREATE INDEX idx_effect_relationships_related ON effect_relationships(related_effect_id);
CREATE INDEX idx_effect_relationships_type ON effect_relationships(relationship_type);

CREATE INDEX idx_quantum_params_name ON quantum_parameters(parameter_name);
CREATE INDEX idx_quantum_params_method ON quantum_parameters(calculation_method);

CREATE INDEX idx_analysis_params_name ON analysis_parameters(parameter_name);
CREATE INDEX idx_analysis_params_type ON analysis_parameters(data_type);

-- Create indexes for web sources
CREATE INDEX idx_web_source_categories_name ON web_source_categories(name);
CREATE INDEX idx_web_source_categories_type ON web_source_categories(source_type);

CREATE INDEX idx_web_sources_category ON web_data_sources(category_id);
CREATE INDEX idx_web_sources_name ON web_data_sources(name);
CREATE INDEX idx_web_sources_active ON web_data_sources(active);
CREATE INDEX idx_web_sources_text_search ON web_data_sources USING gin(text_search_vector);
