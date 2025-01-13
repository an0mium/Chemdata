-- Consolidated regulatory schema
-- Combines regulatory data, compliance tracking, and enhanced monitoring

------------------------------------------
-- Core Regulatory Data
------------------------------------------

-- Regulatory submissions with enhanced tracking
CREATE TABLE IF NOT EXISTS regulatory_submissions (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    submission_id text NOT NULL UNIQUE,
    regulatory_body text NOT NULL, -- e.g., 'FDA', 'EMA', 'PMDA'
    submission_type text NOT NULL, -- e.g., 'IND', 'NDA', 'MAA'
    status text NOT NULL,
    submission_date date NOT NULL,
    decision_date date,
    review_cycle integer,
    review_phase text,
    reviewer_comments jsonb,
    proposed_indication text[],
    review_pathway text,
    priority_status text,
    orphan_status boolean,
    expedited_programs text[],
    contact_info jsonb,
    reference_numbers jsonb,
    approval_date date,
    approval_conditions text[],
    supporting_documents text[],
    -- Enhanced tracking fields
    risk_assessment jsonb,
    timeline_tracking jsonb,
    milestone_status jsonb,
    stakeholder_feedback jsonb,
    resource_allocation jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Regulatory approvals with enhanced compliance
CREATE TABLE IF NOT EXISTS regulatory_approvals (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    submission_id uuid NOT NULL REFERENCES regulatory_submissions(id) ON DELETE CASCADE,
    approval_number text NOT NULL UNIQUE,
    approval_type text NOT NULL,
    approval_date date,
    expiration_date date,
    approved_indications text[],
    marketing_status text,
    post_market_requirements text[],
    labeling_requirements jsonb,
    special_conditions text[],
    reference_documents text[],
    -- Enhanced compliance fields
    compliance_status jsonb,
    renewal_tracking jsonb,
    condition_fulfillment jsonb,
    audit_history jsonb,
    risk_mitigation_status jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Enhanced regulatory communications
CREATE TABLE IF NOT EXISTS regulatory_communications (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    submission_id uuid NOT NULL REFERENCES regulatory_submissions(id) ON DELETE CASCADE,
    communication_date date NOT NULL,
    communication_type text NOT NULL,
    direction text NOT NULL, -- 'incoming' or 'outgoing'
    sender text,
    recipient text,
    subject text NOT NULL,
    content text,
    attachments text[],
    response_required boolean,
    response_deadline date,
    response_status text,
    follow_up_actions text[],
    importance_level text,
    -- Enhanced tracking fields
    thread_tracking jsonb,
    response_metrics jsonb,
    escalation_status jsonb,
    resolution_tracking jsonb,
    stakeholder_engagement jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Compliance Management
------------------------------------------

-- Enhanced regulatory requirements
CREATE TABLE IF NOT EXISTS regulatory_requirements (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    requirement_type text NOT NULL,
    description text NOT NULL,
    regulatory_body text NOT NULL,
    region text,
    effective_date date,
    compliance_deadline date,
    status text NOT NULL,
    documentation_needed text[],
    responsible_party text,
    risk_level text,
    mitigation_plan text,
    verification_method text,
    verification_evidence text[],
    compliance_details jsonb,
    -- Enhanced tracking fields
    requirement_source jsonb,
    impact_assessment jsonb,
    dependency_mapping jsonb,
    change_history jsonb,
    stakeholder_roles jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Enhanced compliance tracking
CREATE TABLE IF NOT EXISTS compliance_tracking (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    requirement_id uuid NOT NULL REFERENCES regulatory_requirements(id) ON DELETE CASCADE,
    assessment_date date NOT NULL,
    compliance_status text NOT NULL,
    verification_method text,
    verification_evidence text[],
    reviewer text,
    findings text[],
    corrective_actions text[],
    follow_up_date date,
    risk_assessment jsonb,
    notes text,
    -- Enhanced monitoring fields
    monitoring_metrics jsonb,
    trend_analysis jsonb,
    deviation_tracking jsonb,
    remediation_progress jsonb,
    effectiveness_metrics jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Safety and Monitoring
------------------------------------------

-- Enhanced safety reports
CREATE TABLE IF NOT EXISTS safety_reports (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    report_type text NOT NULL, -- e.g., 'PSUR', 'PBRER', 'DSUR'
    report_date date NOT NULL,
    reporting_period jsonb,
    report_period_start date,
    report_period_end date,
    severity_level text,
    description text NOT NULL,
    affected_subjects text[],
    safety_concerns text[],
    benefit_risk_assessment text,
    investigation_details jsonb,
    root_cause text,
    corrective_actions text[],
    preventive_measures text[],
    regulatory_actions text[],
    report_conclusions text,
    regulatory_impact text,
    report_status text,
    reference_documents text[],
    -- Enhanced analysis fields
    trend_analysis jsonb,
    signal_detection jsonb,
    risk_categorization jsonb,
    impact_metrics jsonb,
    followup_tracking jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Enhanced safety monitoring
CREATE TABLE IF NOT EXISTS safety_monitoring (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    monitoring_type text NOT NULL,
    start_date date NOT NULL,
    end_date date,
    frequency text,
    parameters jsonb,
    thresholds jsonb,
    responsible_party text,
    reporting_requirements text[],
    alert_conditions jsonb,
    status text,
    -- Enhanced monitoring fields
    monitoring_metrics jsonb,
    alert_history jsonb,
    threshold_breaches jsonb,
    intervention_tracking jsonb,
    effectiveness_measures jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Inspections and Audits
------------------------------------------

-- Enhanced regulatory inspections
CREATE TABLE IF NOT EXISTS regulatory_inspections (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    inspection_type text NOT NULL,
    inspection_date date NOT NULL,
    regulatory_body text NOT NULL,
    inspector_names text[],
    inspector_organization text,
    inspection_scope text[],
    findings jsonb,
    classification text,
    response_required boolean,
    response_deadline date,
    response_status text,
    corrective_actions jsonb,
    follow_up_required boolean,
    inspection_report_ref text,
    -- Enhanced tracking fields
    preparation_status jsonb,
    observation_tracking jsonb,
    response_progress jsonb,
    effectiveness_verification jsonb,
    lesson_learned jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Product Management
------------------------------------------

-- Enhanced labeling changes
CREATE TABLE IF NOT EXISTS labeling_changes (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    change_type text NOT NULL,
    implementation_date date,
    sections_affected text[],
    change_description text,
    reason_for_change text,
    regulatory_requirement text,
    approval_status text,
    reference_documents text[],
    -- Enhanced tracking fields
    impact_assessment jsonb,
    implementation_status jsonb,
    verification_results jsonb,
    distribution_tracking jsonb,
    effectiveness_metrics jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Enhanced post-marketing commitments
CREATE TABLE IF NOT EXISTS post_marketing_commitments (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    approval_id uuid NOT NULL REFERENCES regulatory_approvals(id) ON DELETE CASCADE,
    commitment_type text NOT NULL,
    description text,
    due_date date,
    status text,
    completion_date date,
    submission_date date,
    regulatory_body_feedback text,
    reference_documents text[],
    -- Enhanced tracking fields
    progress_metrics jsonb,
    resource_allocation jsonb,
    milestone_tracking jsonb,
    stakeholder_updates jsonb,
    quality_metrics jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Enhanced post-market surveillance
CREATE TABLE IF NOT EXISTS post_market_surveillance (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    surveillance_type text NOT NULL,
    start_date date NOT NULL,
    end_date date,
    monitoring_parameters jsonb,
    data_sources text[],
    findings jsonb,
    signal_detection_methods text[],
    identified_risks text[],
    mitigation_strategies jsonb,
    reporting_frequency text,
    status text,
    -- Enhanced monitoring fields
    data_quality_metrics jsonb,
    trend_analysis jsonb,
    signal_validation jsonb,
    intervention_tracking jsonb,
    effectiveness_measures jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Enhanced regulatory actions
CREATE TABLE IF NOT EXISTS regulatory_actions (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    action_type text NOT NULL,
    action_date date NOT NULL,
    regulatory_body text NOT NULL,
    description text NOT NULL,
    reason text,
    impact_assessment text,
    required_response text,
    response_deadline date,
    response_status text,
    resolution_date date,
    resolution_details text,
    -- Enhanced tracking fields
    action_tracking jsonb,
    response_progress jsonb,
    effectiveness_metrics jsonb,
    stakeholder_impact jsonb,
    lesson_learned jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Indexes
------------------------------------------

-- Core regulatory indexes
CREATE INDEX IF NOT EXISTS idx_regulatory_submissions_compound ON regulatory_submissions(compound_id);
CREATE INDEX IF NOT EXISTS idx_regulatory_submissions_id ON regulatory_submissions(submission_id);
CREATE INDEX IF NOT EXISTS idx_regulatory_submissions_type ON regulatory_submissions(submission_type);
CREATE INDEX IF NOT EXISTS idx_regulatory_submissions_body ON regulatory_submissions(regulatory_body);
CREATE INDEX IF NOT EXISTS idx_regulatory_submissions_status ON regulatory_submissions(status);
CREATE INDEX IF NOT EXISTS idx_regulatory_submissions_dates ON regulatory_submissions(submission_date, decision_date);

CREATE INDEX IF NOT EXISTS idx_regulatory_approvals_submission ON regulatory_approvals(submission_id);
CREATE INDEX IF NOT EXISTS idx_regulatory_approvals_number ON regulatory_approvals(approval_number);
CREATE INDEX IF NOT EXISTS idx_regulatory_approvals_type ON regulatory_approvals(approval_type);
CREATE INDEX IF NOT EXISTS idx_regulatory_approvals_dates ON regulatory_approvals(approval_date, expiration_date);

CREATE INDEX IF NOT EXISTS idx_regulatory_communications_submission ON regulatory_communications(submission_id);
CREATE INDEX IF NOT EXISTS idx_regulatory_communications_date ON regulatory_communications(communication_date);
CREATE INDEX IF NOT EXISTS idx_regulatory_communications_type ON regulatory_communications(communication_type);
CREATE INDEX IF NOT EXISTS idx_regulatory_communications_direction ON regulatory_communications(direction);
CREATE INDEX IF NOT EXISTS idx_regulatory_communications_status ON regulatory_communications(response_status);

-- Compliance indexes
CREATE INDEX IF NOT EXISTS idx_regulatory_requirements_compound ON regulatory_requirements(compound_id);
CREATE INDEX IF NOT EXISTS idx_regulatory_requirements_type ON regulatory_requirements(requirement_type);
CREATE INDEX IF NOT EXISTS idx_regulatory_requirements_body ON regulatory_requirements(regulatory_body);
CREATE INDEX IF NOT EXISTS idx_regulatory_requirements_status ON regulatory_requirements(status);
CREATE INDEX IF NOT EXISTS idx_regulatory_requirements_dates ON regulatory_requirements(effective_date, compliance_deadline);

CREATE INDEX IF NOT EXISTS idx_compliance_tracking_requirement ON compliance_tracking(requirement_id);
CREATE INDEX IF NOT EXISTS idx_compliance_tracking_date ON compliance_tracking(assessment_date);
CREATE INDEX IF NOT EXISTS idx_compliance_tracking_status ON compliance_tracking(compliance_status);
CREATE INDEX IF NOT EXISTS idx_compliance_tracking_followup ON compliance_tracking(follow_up_date);

-- Safety indexes
CREATE INDEX IF NOT EXISTS idx_safety_reports_compound ON safety_reports(compound_id);
CREATE INDEX IF NOT EXISTS idx_safety_reports_type ON safety_reports(report_type);
CREATE INDEX IF NOT EXISTS idx_safety_reports_date ON safety_reports(report_date);
CREATE INDEX IF NOT EXISTS idx_safety_reports_period ON safety_reports(report_period_start, report_period_end);
CREATE INDEX IF NOT EXISTS idx_safety_reports_severity ON safety_reports(severity_level);
CREATE INDEX IF NOT EXISTS idx_safety_reports_status ON safety_reports(report_status);

CREATE INDEX IF NOT EXISTS idx_safety_monitoring_compound ON safety_monitoring(compound_id);
CREATE INDEX IF NOT EXISTS idx_safety_monitoring_type ON safety_monitoring(monitoring_type);
CREATE INDEX IF NOT EXISTS idx_safety_monitoring_dates ON safety_monitoring(start_date, end_date);
CREATE INDEX IF NOT EXISTS idx_safety_monitoring_status ON safety_monitoring(status);

-- Inspection indexes
CREATE INDEX IF NOT EXISTS idx_regulatory_inspections_compound ON regulatory_inspections(compound_id);
CREATE INDEX IF NOT EXISTS idx_regulatory_inspections_type ON regulatory_inspections(inspection_type);
CREATE INDEX IF NOT EXISTS idx_regulatory_inspections_date ON regulatory_inspections(inspection_date);
CREATE INDEX IF NOT EXISTS idx_regulatory_inspections_body ON regulatory_inspections(regulatory_body);
CREATE INDEX IF NOT EXISTS idx_regulatory_inspections_class ON regulatory_inspections(classification);

-- Product management indexes
CREATE INDEX IF NOT EXISTS idx_labeling_changes_compound ON labeling_changes(compound_id);
CREATE INDEX IF NOT EXISTS idx_labeling_changes_type ON labeling_changes(change_type);
CREATE INDEX IF NOT EXISTS idx_labeling_changes_date ON labeling_changes(implementation_date);
CREATE INDEX IF NOT EXISTS idx_labeling_changes_status ON labeling_changes(approval_status);

CREATE INDEX IF NOT EXISTS idx_post_marketing_approval ON post_marketing_commitments(approval_id);
CREATE INDEX IF NOT EXISTS idx_post_marketing_type ON post_marketing_commitments(commitment_type);
CREATE INDEX IF NOT EXISTS idx_post_marketing_status ON post_marketing_commitments(status);
CREATE INDEX IF NOT EXISTS idx_post_marketing_dates ON post_marketing_commitments(due_date, completion_date);

CREATE INDEX IF NOT EXISTS idx_post_market_surveillance_compound ON post_market_surveillance(compound_id);
CREATE INDEX IF NOT EXISTS idx_post_market_surveillance_type ON post_market_surveillance(surveillance_type);
CREATE INDEX IF NOT EXISTS idx_post_market_surveillance_dates ON post_market_surveillance(start_date, end_date);
CREATE INDEX IF NOT EXISTS idx_post_market_surveillance_status ON post_market_surveillance(status);

CREATE INDEX IF NOT EXISTS idx_regulatory_actions_compound ON regulatory_actions(compound_id);
CREATE INDEX IF NOT EXISTS idx_regulatory_actions_type ON regulatory_actions(action_type);
CREATE INDEX IF NOT EXISTS idx_regulatory_actions_date ON regulatory_actions(action_date);
CREATE INDEX IF NOT EXISTS idx_regulatory_actions_body ON regulatory_actions(regulatory_body);
CREATE INDEX IF NOT EXISTS idx_regulatory_actions_status ON regulatory_actions(response_status);

-- JSONB indexes for enhanced fields
CREATE INDEX IF NOT EXISTS idx_submissions_risk ON regulatory_submissions USING gin (risk_assessment);
CREATE INDEX IF NOT EXISTS idx_approvals_compliance ON regulatory_approvals USING gin (compliance_status);
CREATE INDEX IF NOT EXISTS idx_communications_metrics ON regulatory_communications USING gin (response_metrics);
CREATE INDEX IF NOT EXISTS idx_requirements_impact ON regulatory_requirements USING gin (impact_assessment);
CREATE INDEX IF NOT EXISTS idx_compliance_monitoring ON compliance_tracking USING gin (monitoring_metrics);
CREATE INDEX IF NOT EXISTS idx_safety_trends ON safety_reports USING gin (trend_analysis);
CREATE INDEX IF NOT EXISTS idx_monitoring_metrics ON safety_monitoring USING gin (monitoring_metrics);
CREATE INDEX IF NOT EXISTS idx_inspections_tracking ON regulatory_inspections USING gin (observation_tracking);
CREATE INDEX IF NOT EXISTS idx_labeling_impact ON labeling_changes USING gin (impact_assessment);
CREATE INDEX IF NOT EXISTS idx_commitments_progress ON post_marketing_commitments USING gin (progress_metrics);
CREATE INDEX IF NOT EXISTS idx_surveillance_findings ON post_market_surveillance USING gin (findings);
CREATE INDEX IF NOT EXISTS idx_actions_tracking ON regulatory_actions USING gin (action_tracking);

------------------------------------------
-- Triggers
------------------------------------------

-- Update triggers for timestamp management
DO $$ 
DECLARE
    t text;
BEGIN
    FOR t IN SELECT table_name 
             FROM information_schema.tables 
             WHERE table_schema = 'public' 
             AND table_type = 'BASE TABLE'
             AND table_name LIKE ANY (ARRAY[
                 'regulatory_%',
                 'compliance_%',
                 'safety_%',
                 'post_%',
                 'labeling_%'
             ])
    LOOP
        EXECUTE format('
            CREATE TRIGGER update_%I_modtime 
            BEFORE UPDATE ON %I 
            FOR EACH ROW
            EXECUTE FUNCTION update_updated_at_column();
        ', t, t);
    END LOOP;
END $$;

-- Audit triggers for tracking changes
DO $$ 
DECLARE
    t text;
BEGIN
    FOR t IN SELECT table_name 
             FROM information_schema.tables 
             WHERE table_schema = 'public' 
             AND table_type = 'BASE TABLE'
             AND table_name LIKE ANY (ARRAY[
                 'regulatory_%',
                 'compliance_%',
                 'safety_%',
                 'post_%',
                 'labeling_%'
             ])
    LOOP
        EXECUTE format('
            CREATE TRIGGER audit_%I_trigger
            AFTER INSERT OR UPDATE OR DELETE ON %I
            FOR EACH ROW
            EXECUTE FUNCTION audit_trigger_func();
        ', t, t);
    END LOOP;
END $$;

------------------------------------------
-- Views
------------------------------------------

-- Regulatory overview
CREATE OR REPLACE VIEW regulatory_overview AS
SELECT 
    c.id as compound_id,
    c.name as compound_name,
    COUNT(DISTINCT rs.id) as submission_count,
    COUNT(DISTINCT ra.id) as approval_count,
    COUNT(DISTINCT rc.id) as communication_count,
    COUNT(DISTINCT rr.id) as requirement_count,
    COUNT(DISTINCT ct.id) as compliance_check_count,
    COUNT(DISTINCT sr.id) as safety_report_count,
    COUNT(DISTINCT sm.id) as monitoring_count,
    COUNT(DISTINCT ri.id) as inspection_count,
    COUNT(DISTINCT lc.id) as labeling_change_count,
    COUNT(DISTINCT pmc.id) as post_marketing_count,
    COUNT(DISTINCT pms.id) as surveillance_count,
    COUNT(DISTINCT ract.id) as action_count
FROM compounds c
LEFT JOIN regulatory_submissions rs ON c.id = rs.compound_id
LEFT JOIN regulatory_approvals ra ON rs.id = ra.submission_id
LEFT JOIN regulatory_communications rc ON rs.id = rc.submission_id
LEFT JOIN regulatory_requirements rr ON c.id = rr.compound_id
LEFT JOIN compliance_tracking ct ON rr.id = ct.requirement_id
LEFT JOIN safety_reports sr ON c.id = sr.compound_id
LEFT JOIN safety_monitoring sm ON c.id = sm.compound_id
LEFT JOIN regulatory_inspections ri ON c.id = ri.compound_id
LEFT JOIN labeling_changes lc ON c.id = lc.compound_id
LEFT JOIN post_marketing_commitments pmc ON ra.id = pmc.approval_id
LEFT JOIN post_market_surveillance pms ON c.id = pms.compound_id
LEFT JOIN regulatory_actions ract ON c.id = ract.compound_id
GROUP BY c.id, c.name;

-- Compliance status
CREATE OR REPLACE VIEW compliance_status AS
SELECT 
    c.id as compound_id,
    c.name as compound_name,
    rr.requirement_type,
    rr.status as requirement_status,
    ct.compliance_status,
    ct.assessment_date,
    ct.follow_up_date,
    ct.risk_assessment,
    ct.corrective_actions,
    rr.regulatory_body,
    rr.region
FROM compounds c
JOIN regulatory_requirements rr ON c.id = rr.compound_id
LEFT JOIN compliance_tracking ct ON rr.id = ct.requirement_id;

-- Safety monitoring status
CREATE OR REPLACE VIEW safety_monitoring_status AS
SELECT 
    c.id as compound_id,
    c.name as compound_name,
    sm.monitoring_type,
    sm.status,
    sm.parameters,
    sm.thresholds,
    sm.alert_conditions,
    sm.monitoring_metrics,
    sm.threshold_breaches,
    sm.intervention_tracking,
    sm.effectiveness_measures
FROM compounds c
JOIN safety_monitoring sm ON c.id = sm.compound_id;

------------------------------------------
-- Comments
------------------------------------------

COMMENT ON TABLE regulatory_submissions IS 'Regulatory submission records with enhanced tracking';
COMMENT ON TABLE regulatory_approvals IS 'Approved regulatory submissions with enhanced compliance tracking';
COMMENT ON TABLE regulatory_communications IS 'Communications with regulatory authorities with enhanced monitoring';
COMMENT ON TABLE regulatory_requirements IS 'Regulatory requirements and obligations with enhanced tracking';
COMMENT ON TABLE compliance_tracking IS 'Compliance assessment and verification with enhanced monitoring';
COMMENT ON TABLE safety_reports IS 'Safety-related reports and findings with enhanced analysis';
COMMENT ON TABLE safety_monitoring IS 'Safety monitoring programs with enhanced metrics tracking';
COMMENT ON TABLE regulatory_inspections IS 'Regulatory inspection details with enhanced tracking';
COMMENT ON TABLE labeling_changes IS 'Product labeling changes with enhanced implementation tracking';
COMMENT ON TABLE post_marketing_commitments IS 'Post-marketing study commitments with enhanced progress tracking';
COMMENT ON TABLE post_market_surveillance IS 'Post-market safety surveillance with enhanced monitoring';
COMMENT ON TABLE regulatory_actions IS 'Regulatory authority actions with enhanced response tracking';

COMMENT ON VIEW regulatory_overview IS 'Overview of regulatory activities and metrics per compound';
COMMENT ON VIEW compliance_status IS 'Current compliance status and tracking for regulatory requirements';
COMMENT ON VIEW safety_monitoring_status IS 'Current safety monitoring status and metrics';
