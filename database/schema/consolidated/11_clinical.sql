-- Consolidated clinical schema
-- Combines clinical trials, experience reports, and outcomes

------------------------------------------
-- Clinical Trials and Participants
------------------------------------------

-- Clinical trials
CREATE TABLE IF NOT EXISTS clinical_trials (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    trial_id text NOT NULL UNIQUE,
    title text NOT NULL,
    phase text,
    status text NOT NULL,
    start_date date,
    end_date date,
    enrollment integer,
    study_type text,
    design text[],
    primary_outcome text[],
    secondary_outcome text[],
    inclusion_criteria text[],
    exclusion_criteria text[],
    locations text[],
    results_summary jsonb,
    adverse_events jsonb,
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Patient data
CREATE TABLE IF NOT EXISTS patient_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    trial_id uuid NOT NULL REFERENCES clinical_trials(id) ON DELETE CASCADE,
    patient_id text NOT NULL,
    demographics jsonb,
    medical_history text[],
    concurrent_medications text[],
    vital_signs jsonb,
    lab_results jsonb,
    adverse_events text[],
    outcome_measures jsonb,
    follow_up_data jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE (trial_id, patient_id)
);

-- Trial demographics summary
CREATE TABLE IF NOT EXISTS trial_demographics (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    trial_id uuid NOT NULL REFERENCES clinical_trials(id) ON DELETE CASCADE,
    age_range jsonb,
    gender_distribution jsonb,
    ethnicity_distribution jsonb,
    health_status text[],
    concurrent_medications jsonb,
    comorbidities jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Protocol deviations
CREATE TABLE IF NOT EXISTS protocol_deviations (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    trial_id uuid NOT NULL REFERENCES clinical_trials(id) ON DELETE CASCADE,
    deviation_date date NOT NULL,
    deviation_type text NOT NULL,
    description text NOT NULL,
    severity_level text,
    impact_assessment text,
    affected_participants text[],
    root_cause_analysis text,
    corrective_actions text[],
    preventive_measures text[],
    reporting_status text,
    resolution_status text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Experience Reports and Outcomes
------------------------------------------

-- Clinical experience reports
CREATE TABLE IF NOT EXISTS clinical_experience_reports (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    report_type text NOT NULL,
    report_date date NOT NULL,
    setting text,
    dosage jsonb,
    route_of_administration text,
    duration interval,
    effects_timeline jsonb,
    subjective_effects text[],
    physical_effects text[],
    cognitive_effects text[],
    emotional_effects text[],
    side_effects text[],
    interactions text[],
    overall_experience text,
    therapeutic_value text,
    risk_assessment jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Therapeutic outcomes
CREATE TABLE IF NOT EXISTS therapeutic_outcomes (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    condition text NOT NULL,
    treatment_protocol jsonb,
    efficacy_rating double precision,
    response_rate double precision,
    remission_rate double precision,
    relapse_rate double precision,
    side_effect_profile jsonb,
    quality_of_life_measures jsonb,
    long_term_outcomes jsonb,
    cost_effectiveness jsonb,
    evidence_quality text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Dose response data
CREATE TABLE IF NOT EXISTS dose_response_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    administration_route text NOT NULL,
    dose_level double precision NOT NULL,
    dose_unit text NOT NULL,
    response_measure text NOT NULL,
    response_value double precision,
    response_variability jsonb,
    conditions jsonb,
    population_data jsonb,
    methodology text,
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Time course data
CREATE TABLE IF NOT EXISTS time_course_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    administration_route text NOT NULL,
    dose_amount double precision,
    dose_unit text,
    time_points jsonb,
    effect_measures jsonb,
    peak_time interval,
    peak_intensity double precision,
    duration_total interval,
    notes text,
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Treatment Protocols and Responses
------------------------------------------

-- Treatment protocols
CREATE TABLE IF NOT EXISTS treatment_protocols (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    protocol_name text NOT NULL,
    indication text NOT NULL,
    dosing_schedule jsonb,
    administration_route text,
    duration text,
    monitoring_requirements text[],
    contraindications text[],
    precautions text[],
    drug_interactions text[],
    adjustment_criteria jsonb,
    success_metrics jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Treatment responses
CREATE TABLE IF NOT EXISTS treatment_responses (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    patient_id uuid NOT NULL REFERENCES patient_data(id) ON DELETE CASCADE,
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    response_type text NOT NULL,
    response_measure double precision,
    response_duration interval,
    side_effects text[],
    tolerability_score double precision,
    adherence_rate double precision,
    discontinuation_reason text,
    predictive_factors jsonb,
    biomarkers jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Quality Metrics and Assessments
------------------------------------------

-- Clinical assessments
CREATE TABLE IF NOT EXISTS clinical_assessments (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    patient_id uuid NOT NULL REFERENCES patient_data(id) ON DELETE CASCADE,
    assessment_type text NOT NULL,
    assessment_date date NOT NULL,
    clinician text,
    symptoms text[],
    severity_scores jsonb,
    functional_measures jsonb,
    cognitive_measures jsonb,
    physical_measures jsonb,
    laboratory_results jsonb,
    imaging_results jsonb,
    recommendations text[],
    follow_up_plan text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Quality of life measures
CREATE TABLE IF NOT EXISTS quality_of_life_measures (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    patient_id uuid NOT NULL REFERENCES patient_data(id) ON DELETE CASCADE,
    assessment_date date NOT NULL,
    physical_health_score double precision,
    mental_health_score double precision,
    social_functioning_score double precision,
    daily_activities_score double precision,
    pain_interference_score double precision,
    sleep_quality_score double precision,
    overall_wellbeing_score double precision,
    detailed_measures jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Quality metrics
CREATE TABLE IF NOT EXISTS quality_metrics (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    trial_id uuid NOT NULL REFERENCES clinical_trials(id) ON DELETE CASCADE,
    metric_type text NOT NULL,
    measurement_date date NOT NULL,
    metric_value double precision,
    target_range jsonb,
    deviation_reason text,
    action_required boolean,
    action_taken text,
    follow_up_date date,
    reviewer text,
    notes text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Compliance tracking
CREATE TABLE IF NOT EXISTS compliance_tracking (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    requirement_type text NOT NULL,
    description text NOT NULL,
    due_date date,
    status text NOT NULL,
    responsible_party text,
    verification_method text,
    verification_date date,
    verification_evidence text[],
    compliance_level text,
    risk_assessment jsonb,
    mitigation_steps text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Indexes
------------------------------------------

-- Clinical trials indexes
CREATE INDEX IF NOT EXISTS idx_clinical_trials_compound ON clinical_trials(compound_id);
CREATE INDEX IF NOT EXISTS idx_clinical_trials_id ON clinical_trials(trial_id);
CREATE INDEX IF NOT EXISTS idx_clinical_trials_phase ON clinical_trials(phase);
CREATE INDEX IF NOT EXISTS idx_clinical_trials_status ON clinical_trials(status);
CREATE INDEX IF NOT EXISTS idx_clinical_trials_dates ON clinical_trials(start_date, end_date);

-- Patient data indexes
CREATE INDEX IF NOT EXISTS idx_patient_data_trial ON patient_data(trial_id);
CREATE INDEX IF NOT EXISTS idx_patient_data_patient ON patient_data(patient_id);

-- Trial demographics indexes
CREATE INDEX IF NOT EXISTS idx_trial_demographics_trial ON trial_demographics(trial_id);

-- Protocol deviation indexes
CREATE INDEX IF NOT EXISTS idx_protocol_deviations_trial ON protocol_deviations(trial_id);
CREATE INDEX IF NOT EXISTS idx_protocol_deviations_date ON protocol_deviations(deviation_date);
CREATE INDEX IF NOT EXISTS idx_protocol_deviations_type ON protocol_deviations(deviation_type);
CREATE INDEX IF NOT EXISTS idx_protocol_deviations_severity ON protocol_deviations(severity_level);

-- Experience report indexes
CREATE INDEX IF NOT EXISTS idx_experience_reports_compound ON clinical_experience_reports(compound_id);
CREATE INDEX IF NOT EXISTS idx_experience_reports_type ON clinical_experience_reports(report_type);
CREATE INDEX IF NOT EXISTS idx_experience_reports_date ON clinical_experience_reports(report_date);

-- Therapeutic outcome indexes
CREATE INDEX IF NOT EXISTS idx_therapeutic_outcomes_compound ON therapeutic_outcomes(compound_id);
CREATE INDEX IF NOT EXISTS idx_therapeutic_outcomes_condition ON therapeutic_outcomes(condition);

-- Dose response indexes
CREATE INDEX IF NOT EXISTS idx_dose_response_compound ON dose_response_data(compound_id);
CREATE INDEX IF NOT EXISTS idx_dose_response_route ON dose_response_data(administration_route);

-- Time course indexes
CREATE INDEX IF NOT EXISTS idx_time_course_compound ON time_course_data(compound_id);
CREATE INDEX IF NOT EXISTS idx_time_course_route ON time_course_data(administration_route);

-- Treatment protocol indexes
CREATE INDEX IF NOT EXISTS idx_treatment_protocols_compound ON treatment_protocols(compound_id);
CREATE INDEX IF NOT EXISTS idx_treatment_protocols_indication ON treatment_protocols(indication);

-- Treatment response indexes
CREATE INDEX IF NOT EXISTS idx_treatment_responses_patient ON treatment_responses(patient_id);
CREATE INDEX IF NOT EXISTS idx_treatment_responses_compound ON treatment_responses(compound_id);
CREATE INDEX IF NOT EXISTS idx_treatment_responses_type ON treatment_responses(response_type);

-- Clinical assessment indexes
CREATE INDEX IF NOT EXISTS idx_clinical_assessments_patient ON clinical_assessments(patient_id);
CREATE INDEX IF NOT EXISTS idx_clinical_assessments_type ON clinical_assessments(assessment_type);
CREATE INDEX IF NOT EXISTS idx_clinical_assessments_date ON clinical_assessments(assessment_date);

-- Quality of life indexes
CREATE INDEX IF NOT EXISTS idx_quality_of_life_patient ON quality_of_life_measures(patient_id);
CREATE INDEX IF NOT EXISTS idx_quality_of_life_date ON quality_of_life_measures(assessment_date);

-- Quality metric indexes
CREATE INDEX IF NOT EXISTS idx_quality_metrics_trial ON quality_metrics(trial_id);
CREATE INDEX IF NOT EXISTS idx_quality_metrics_type ON quality_metrics(metric_type);
CREATE INDEX IF NOT EXISTS idx_quality_metrics_date ON quality_metrics(measurement_date);

-- Compliance tracking indexes
CREATE INDEX IF NOT EXISTS idx_compliance_tracking_compound ON compliance_tracking(compound_id);
CREATE INDEX IF NOT EXISTS idx_compliance_tracking_type ON compliance_tracking(requirement_type);
CREATE INDEX IF NOT EXISTS idx_compliance_tracking_status ON compliance_tracking(status);
CREATE INDEX IF NOT EXISTS idx_compliance_tracking_dates ON compliance_tracking(due_date, verification_date);

------------------------------------------
-- Update Triggers
------------------------------------------

-- Clinical trials triggers
CREATE TRIGGER update_clinical_trials_modtime
    BEFORE UPDATE ON clinical_trials
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_patient_data_modtime
    BEFORE UPDATE ON patient_data
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_trial_demographics_modtime
    BEFORE UPDATE ON trial_demographics
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_protocol_deviations_modtime
    BEFORE UPDATE ON protocol_deviations
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Experience report triggers
CREATE TRIGGER update_experience_reports_modtime
    BEFORE UPDATE ON clinical_experience_reports
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_therapeutic_outcomes_modtime
    BEFORE UPDATE ON therapeutic_outcomes
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_dose_response_modtime
    BEFORE UPDATE ON dose_response_data
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_time_course_modtime
    BEFORE UPDATE ON time_course_data
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Treatment protocol triggers
CREATE TRIGGER update_treatment_protocols_modtime
    BEFORE UPDATE ON treatment_protocols
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_treatment_responses_modtime
    BEFORE UPDATE ON treatment_responses
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Assessment triggers
CREATE TRIGGER update_clinical_assessments_modtime
    BEFORE UPDATE ON clinical_assessments
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_quality_of_life_modtime
    BEFORE UPDATE ON quality_of_life_measures
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_quality_metrics_modtime
    BEFORE UPDATE ON quality_metrics
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_compliance_tracking_modtime
    BEFORE UPDATE ON compliance_tracking
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

------------------------------------------
-- Audit Triggers
------------------------------------------

-- Clinical trials audit
CREATE TRIGGER audit_clinical_trials_trigger
    AFTER INSERT OR UPDATE OR DELETE ON clinical_trials
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_patient_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON patient_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_trial_demographics_trigger
    AFTER INSERT OR UPDATE OR DELETE ON trial_demographics
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_protocol_deviations_trigger
    AFTER INSERT OR UPDATE OR DELETE ON protocol_deviations
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Experience report audit
CREATE TRIGGER audit_experience_reports_trigger
    AFTER INSERT OR UPDATE OR DELETE ON clinical_experience_reports
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_therapeutic_outcomes_trigger
    AFTER INSERT OR UPDATE OR DELETE ON therapeutic_outcomes
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_dose_response_trigger
    AFTER INSERT OR UPDATE OR DELETE ON dose_response_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_time_course_trigger
    AFTER INSERT OR UPDATE OR DELETE ON time_course_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Treatment protocol audit
CREATE TRIGGER audit_treatment_protocols_trigger
    AFTER INSERT OR UPDATE OR DELETE ON treatment_protocols
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_treatment_responses_trigger
    AFTER INSERT OR UPDATE OR DELETE ON treatment_responses
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Assessment audit
CREATE TRIGGER audit_clinical_assessments_trigger
    AFTER INSERT OR UPDATE OR DELETE ON clinical_assessments
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_quality_of_life_trigger
    AFTER INSERT OR UPDATE OR DELETE ON quality_of_life_measures
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_quality_metrics_trigger
    AFTER INSERT OR UPDATE OR DELETE ON quality_metrics
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_compliance_tracking_trigger
    AFTER INSERT OR UPDATE OR DELETE ON compliance_tracking
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

------------------------------------------
-- Views
------------------------------------------

-- Trial overview
CREATE OR REPLACE VIEW trial_overview AS
SELECT 
    t.id as trial_id,
    t.compound_id,
    t.trial_id as identifier,
    t.title,
    t.phase,
    t.status,
    t.start_date,
    t.end_date,
    t.enrollment,
    COUNT(DISTINCT p.id) as participant_count,
    COUNT(DISTINCT pd.id) as deviation_count,
    COUNT(DISTINCT tr.id) as response_count,
    COUNT(DISTINCT ca.id) as assessment_count
FROM clinical_trials t
LEFT JOIN patient_data p ON t.id = p.trial_id
LEFT JOIN protocol_deviations pd ON t.id = pd.trial_id
LEFT JOIN treatment_responses tr ON tr.patient_id IN (SELECT id FROM patient_data WHERE trial_id = t.id)
LEFT JOIN clinical_assessments ca ON ca.patient_id IN (SELECT id FROM patient_data WHERE trial_id = t.id)
GROUP BY t.id, t.compound_id, t.trial_id, t.title, t.phase, t.status, t.start_date, t.end_date, t.enrollment;

-- Patient outcome summary
CREATE OR REPLACE VIEW patient_outcome_summary AS
SELECT 
    p.id as patient_id,
    p.trial_id,
    t.compound_id,
    COUNT(DISTINCT tr.id) as response_count,
    COUNT(DISTINCT ca.id) as assessment_count,
    COUNT(DISTINCT qol.id) as qol_measurement_count,
    AVG(tr.response_measure) as avg_response_measure,
    AVG(tr.tolerability_score) as avg_tolerability,
    AVG(tr.adherence_rate) as avg_adherence,
    AVG(qol.overall_wellbeing_score) as avg_wellbeing
FROM patient_data p
JOIN clinical_trials t ON p.trial_id = t.id
LEFT JOIN treatment_responses tr ON p.id = tr.patient_id
LEFT JOIN clinical_assessments ca ON p.id = ca.patient_id
LEFT JOIN quality_of_life_measures qol ON p.id = qol.patient_id
GROUP BY p.id, p.trial_id, t.compound_id;

-- Therapeutic efficacy summary
CREATE OR REPLACE VIEW therapeutic_efficacy_summary AS
SELECT 
    compound_id,
    condition,
    COUNT(*) as outcome_count,
    AVG(efficacy_rating) as avg_efficacy,
    AVG(response_rate) as avg_response_rate,
    AVG(remission_rate) as avg_remission_rate,
    AVG(relapse_rate) as avg_relapse_rate
FROM therapeutic_outcomes
GROUP BY compound_id, condition;

-- Table comments
COMMENT ON TABLE clinical_trials IS 'Clinical trial metadata and design information';
COMMENT ON TABLE patient_data IS 'Individual patient data and outcomes';
COMMENT ON TABLE trial_demographics IS 'Demographic summary for trial participants';
COMMENT ON TABLE protocol_deviations IS 'Protocol deviation tracking and resolution';
COMMENT ON TABLE clinical_experience_reports IS 'Clinical experience and observation reports';
COMMENT ON TABLE therapeutic_outcomes IS 'Therapeutic outcome measurements and analysis';
COMMENT ON TABLE dose_response_data IS 'Dose-response relationships and data';
COMMENT ON TABLE time_course_data IS 'Time course of effects and responses';
COMMENT ON TABLE treatment_protocols IS 'Treatment protocols and guidelines';
COMMENT ON TABLE treatment_responses IS 'Individual treatment response data';
COMMENT ON TABLE clinical_assessments IS 'Clinical assessment records';
COMMENT ON TABLE quality_of_life_measures IS 'Quality of life measurements';
COMMENT ON TABLE quality_metrics IS 'Quality control metrics';
COMMENT ON TABLE compliance_tracking IS 'Compliance requirement tracking';
