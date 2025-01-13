-- Safety and toxicity data tables

-- Base toxicity assay data
CREATE TABLE IF NOT EXISTS toxicity_assays (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    assay_type text NOT NULL, -- e.g., 'in_vitro', 'in_vivo'
    cell_line text, -- For in vitro assays
    organism text, -- For in vivo assays
    endpoint text NOT NULL, -- e.g., 'cytotoxicity', 'genotoxicity'
    concentration double precision,
    concentration_unit text,
    exposure_time interval,
    result_value double precision,
    result_unit text,
    result_type text, -- e.g., 'IC50', 'LD50', 'EC50'
    confidence_score double precision,
    protocol_details jsonb,
    reference_doi text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Cytotoxicity data
CREATE TABLE IF NOT EXISTS cytotoxicity_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    cell_line text NOT NULL,
    tissue_type text,
    assay_method text, -- e.g., 'MTT', 'LDH', 'ATP'
    exposure_time interval,
    ic50_value double precision,
    ic50_unit text,
    cell_viability double precision, -- Percentage
    cytotoxicity_mechanism text[],
    morphological_changes text[],
    reference_doi text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Genotoxicity data
CREATE TABLE IF NOT EXISTS genotoxicity_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    test_type text NOT NULL, -- e.g., 'Ames', 'Micronucleus'
    organism text,
    metabolic_activation boolean,
    result text NOT NULL, -- e.g., 'positive', 'negative'
    mutation_type text[],
    dna_damage_type text[],
    mechanism text[],
    reference_doi text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Carcinogenicity data
CREATE TABLE IF NOT EXISTS carcinogenicity_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    study_type text NOT NULL,
    species text,
    duration interval,
    dose_levels jsonb, -- Array of dose levels and units
    tumor_types jsonb, -- Types and incidences
    mechanism text[],
    histopathology jsonb,
    reference_doi text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Metabolic toxicity
CREATE TABLE IF NOT EXISTS metabolic_toxicity (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    enzyme_affected text[], -- e.g., ['CYP3A4', 'CYP2D6']
    inhibition_type text, -- e.g., 'competitive', 'irreversible'
    ki_value double precision,
    ki_unit text,
    metabolites jsonb, -- Known toxic metabolites
    pathway_disruption text[],
    clinical_significance text,
    reference_doi text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Developmental toxicity
CREATE TABLE IF NOT EXISTS developmental_toxicity (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    study_type text NOT NULL,
    species text,
    exposure_period text,
    dose_levels jsonb,
    developmental_effects jsonb,
    mechanism text[],
    teratogenicity boolean,
    reference_doi text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Immunotoxicity data
CREATE TABLE IF NOT EXISTS immunotoxicity_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    study_type text NOT NULL,
    immune_parameters text[], -- e.g., ['T-cell function', 'cytokine levels']
    effect_type text, -- e.g., 'immunosuppression', 'immunostimulation'
    mechanism text[],
    clinical_relevance text,
    reference_doi text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Neurotoxicity data
CREATE TABLE IF NOT EXISTS neurotoxicity_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    study_type text NOT NULL,
    brain_regions_affected text[],
    behavioral_effects text[],
    cellular_effects text[],
    mechanism text[],
    reversibility text,
    long_term_effects jsonb,
    reference_doi text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Cardiotoxicity data
CREATE TABLE IF NOT EXISTS cardiotoxicity_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    study_type text NOT NULL,
    cardiac_effects text[], -- e.g., ['QT prolongation', 'arrhythmia']
    mechanism text[],
    herg_ic50 double precision,
    ecg_changes jsonb,
    hemodynamic_effects jsonb,
    reference_doi text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Hepatotoxicity data
CREATE TABLE IF NOT EXISTS hepatotoxicity_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    study_type text NOT NULL,
    liver_effects text[],
    mechanism text[],
    enzyme_elevations jsonb, -- ALT, AST, etc.
    histopathology jsonb,
    clinical_significance text,
    reference_doi text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Safety documentation and hazard classification
CREATE TABLE IF NOT EXISTS safety_documents (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    sds_url text,
    msds_url text,
    sds_last_updated date,
    sds_provider text,
    sds_version text,
    safety_data_sheet jsonb, -- Structured SDS data
    handling_precautions text[],
    storage_precautions text[],
    disposal_instructions text[],
    first_aid_measures jsonb,
    firefighting_measures jsonb,
    accidental_release_measures jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Hazard classifications (GHS, NFPA, HMIS)
CREATE TABLE IF NOT EXISTS hazard_classifications (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    ghs_classifications text[], -- e.g., ['H200', 'H301', 'H315']
    signal_word text, -- 'Danger' or 'Warning'
    pictograms text[], -- e.g., ['GHS01', 'GHS06']
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
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    CONSTRAINT valid_nfpa_ratings CHECK (
        nfpa_health BETWEEN 0 AND 4 AND
        nfpa_fire BETWEEN 0 AND 4 AND
        nfpa_reactivity BETWEEN 0 AND 4
    ),
    CONSTRAINT valid_hmis_ratings CHECK (
        hmis_health BETWEEN 0 AND 4 AND
        hmis_fire BETWEEN 0 AND 4 AND
        hmis_physical BETWEEN 0 AND 4
    )
);

-- Storage and handling safety requirements
CREATE TABLE IF NOT EXISTS storage_requirements (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
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
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Personal protective equipment requirements
CREATE TABLE IF NOT EXISTS ppe_requirements (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    eye_protection text,
    skin_protection text,
    respiratory_protection text,
    hand_protection text,
    body_protection text,
    minimum_ppe_rating text,
    special_requirements text[],
    exposure_limits jsonb, -- Various exposure limits (PEL, TLV, etc.)
    monitoring_requirements text[],
    decontamination_procedures text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Emergency response procedures
CREATE TABLE IF NOT EXISTS emergency_procedures (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    spill_response text[],
    fire_fighting_measures text[],
    first_aid_procedures jsonb,
    evacuation_criteria text[],
    emergency_contacts jsonb,
    special_hazards text[],
    cleanup_procedures text[],
    disposal_procedures text[],
    reporting_requirements text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Add audit triggers for toxicity tables
DROP TRIGGER IF EXISTS audit_toxicity_assays_trigger ON toxicity_assays;
DROP TRIGGER IF EXISTS audit_cytotoxicity_data_trigger ON cytotoxicity_data;
DROP TRIGGER IF EXISTS audit_genotoxicity_data_trigger ON genotoxicity_data;
DROP TRIGGER IF EXISTS audit_carcinogenicity_data_trigger ON carcinogenicity_data;
DROP TRIGGER IF EXISTS audit_metabolic_toxicity_trigger ON metabolic_toxicity;
DROP TRIGGER IF EXISTS audit_developmental_toxicity_trigger ON developmental_toxicity;
DROP TRIGGER IF EXISTS audit_immunotoxicity_data_trigger ON immunotoxicity_data;
DROP TRIGGER IF EXISTS audit_neurotoxicity_data_trigger ON neurotoxicity_data;
DROP TRIGGER IF EXISTS audit_cardiotoxicity_data_trigger ON cardiotoxicity_data;
DROP TRIGGER IF EXISTS audit_hepatotoxicity_data_trigger ON hepatotoxicity_data;

CREATE TRIGGER audit_toxicity_assays_trigger
    AFTER INSERT OR UPDATE OR DELETE ON toxicity_assays
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_cytotoxicity_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON cytotoxicity_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_genotoxicity_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON genotoxicity_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_carcinogenicity_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON carcinogenicity_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_metabolic_toxicity_trigger
    AFTER INSERT OR UPDATE OR DELETE ON metabolic_toxicity
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_developmental_toxicity_trigger
    AFTER INSERT OR UPDATE OR DELETE ON developmental_toxicity
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_immunotoxicity_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON immunotoxicity_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_neurotoxicity_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON neurotoxicity_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_cardiotoxicity_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON cardiotoxicity_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_hepatotoxicity_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON hepatotoxicity_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Add audit triggers for safety documentation and hazard classification
DROP TRIGGER IF EXISTS audit_safety_documents_trigger ON safety_documents;
DROP TRIGGER IF EXISTS audit_hazard_classifications_trigger ON hazard_classifications;
DROP TRIGGER IF EXISTS audit_storage_requirements_trigger ON storage_requirements;
DROP TRIGGER IF EXISTS audit_ppe_requirements_trigger ON ppe_requirements;
DROP TRIGGER IF EXISTS audit_emergency_procedures_trigger ON emergency_procedures;

CREATE TRIGGER audit_safety_documents_trigger
    AFTER INSERT OR UPDATE OR DELETE ON safety_documents
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_hazard_classifications_trigger
    AFTER INSERT OR UPDATE OR DELETE ON hazard_classifications
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_storage_requirements_trigger
    AFTER INSERT OR UPDATE OR DELETE ON storage_requirements
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_ppe_requirements_trigger
    AFTER INSERT OR UPDATE OR DELETE ON ppe_requirements
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_emergency_procedures_trigger
    AFTER INSERT OR UPDATE OR DELETE ON emergency_procedures
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Add table comments
COMMENT ON TABLE toxicity_assays IS 'General toxicity assay data for compounds';
COMMENT ON TABLE cytotoxicity_data IS 'Cell-based toxicity data';
COMMENT ON TABLE genotoxicity_data IS 'DNA and chromosome damage data';
COMMENT ON TABLE carcinogenicity_data IS 'Cancer-related toxicity data';
COMMENT ON TABLE metabolic_toxicity IS 'Metabolic enzyme interactions and toxicity';
COMMENT ON TABLE developmental_toxicity IS 'Developmental and reproductive toxicity data';
COMMENT ON TABLE immunotoxicity_data IS 'Immune system toxicity data';
COMMENT ON TABLE neurotoxicity_data IS 'Nervous system toxicity data';
COMMENT ON TABLE cardiotoxicity_data IS 'Cardiovascular system toxicity data';
COMMENT ON TABLE hepatotoxicity_data IS 'Liver toxicity data';

COMMENT ON TABLE safety_documents IS 'Safety documentation including SDS and handling instructions';
COMMENT ON TABLE hazard_classifications IS 'GHS and other hazard classification systems for compounds';
COMMENT ON TABLE storage_requirements IS 'Storage conditions and safety requirements';
COMMENT ON TABLE ppe_requirements IS 'Personal protective equipment and exposure control requirements';
COMMENT ON TABLE emergency_procedures IS 'Emergency response and spill control procedures';

-- Create views for safety analysis
CREATE OR REPLACE VIEW compound_safety_overview AS
SELECT 
    c.id as compound_id,
    c.name as compound_name,
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
    COUNT(DISTINCT ta.id) as toxicity_assay_count,
    COUNT(DISTINCT cd.id) as cytotoxicity_data_count,
    COUNT(DISTINCT gd.id) as genotoxicity_data_count
FROM compounds c
LEFT JOIN hazard_classifications hc ON c.id = hc.compound_id
LEFT JOIN storage_requirements sr ON c.id = sr.compound_id
LEFT JOIN ppe_requirements pr ON c.id = pr.compound_id
LEFT JOIN emergency_procedures ep ON c.id = ep.compound_id
LEFT JOIN toxicity_assays ta ON c.id = ta.compound_id
LEFT JOIN cytotoxicity_data cd ON c.id = cd.compound_id
LEFT JOIN genotoxicity_data gd ON c.id = gd.compound_id
GROUP BY c.id, c.name, hc.ghs_classifications, hc.signal_word, hc.pictograms,
         hc.nfpa_health, hc.nfpa_fire, hc.nfpa_reactivity, hc.nfpa_special,
         sr.storage_conditions, sr.incompatible_materials,
         pr.minimum_ppe_rating, ep.special_hazards;
