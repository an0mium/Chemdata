-- Analysis schema combining SAR, research, and analysis functionality

-- Structure-activity relationship analysis
CREATE TABLE IF NOT EXISTS sar_analysis (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
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
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Activity correlations
CREATE TABLE IF NOT EXISTS activity_correlations (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    feature_id uuid NOT NULL, -- Can reference either structural_features or pharmacophore_features
    activity_type text NOT NULL,
    correlation_coefficient double precision,
    statistical_significance double precision,
    analysis_method text,
    sample_size integer,
    confidence_interval jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Structure similarity analysis
CREATE TABLE IF NOT EXISTS structure_similarity (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id_1 uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    compound_id_2 uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    similarity_metric text NOT NULL,
    similarity_score double precision,
    comparison_method text,
    fingerprint_type text,
    calculation_parameters jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- SAR patterns
CREATE TABLE IF NOT EXISTS sar_patterns (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    pattern_type text NOT NULL,
    structural_elements text[],
    activity_impact jsonb,
    confidence_score double precision,
    supporting_compounds text[],
    detection_method text,
    validation_status text,
    notes text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Activity cliffs
CREATE TABLE IF NOT EXISTS activity_cliffs (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_pair uuid[] NOT NULL, -- Array of two compound IDs
    activity_difference double precision,
    structural_similarity double precision,
    cliff_magnitude double precision,
    activity_type text,
    detection_method text,
    significance_score double precision,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Scaffold analysis
CREATE TABLE IF NOT EXISTS scaffold_analysis (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    scaffold_smiles text NOT NULL,
    compound_count integer,
    average_activity double precision,
    activity_range jsonb,
    diversity_score double precision,
    important_substitutions text[],
    analysis_method text,
    notes text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Pathway analysis
CREATE TABLE IF NOT EXISTS pathway_analysis (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    pathway_name text NOT NULL,
    pathway_type text NOT NULL,
    affected_proteins text[],
    regulation_effects jsonb,
    downstream_effects jsonb,
    feedback_mechanisms jsonb,
    pathway_crosstalk jsonb,
    temporal_dynamics jsonb,
    tissue_specificity jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Systems biology analysis
CREATE TABLE IF NOT EXISTS systems_biology_analysis (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    analysis_level text NOT NULL,
    network_effects jsonb,
    cellular_responses jsonb,
    metabolic_impact jsonb,
    signaling_cascades jsonb,
    regulatory_networks jsonb,
    adaptation_mechanisms jsonb,
    system_robustness jsonb,
    emergent_properties jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Research findings
CREATE TABLE IF NOT EXISTS research_findings (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    finding_type text NOT NULL,
    description text NOT NULL,
    methodology text[],
    experimental_data jsonb,
    statistical_analysis jsonb,
    conclusions text[],
    limitations text[],
    future_directions text[],
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Create indexes
CREATE INDEX IF NOT EXISTS idx_sar_analysis_compound ON sar_analysis(compound_id);
CREATE INDEX IF NOT EXISTS idx_sar_analysis_type ON sar_analysis(analysis_type);

CREATE INDEX IF NOT EXISTS idx_activity_correlations_compound ON activity_correlations(compound_id);
CREATE INDEX IF NOT EXISTS idx_activity_correlations_feature ON activity_correlations(feature_id);
CREATE INDEX IF NOT EXISTS idx_activity_correlations_type ON activity_correlations(activity_type);

CREATE INDEX IF NOT EXISTS idx_structure_similarity_compound1 ON structure_similarity(compound_id_1);
CREATE INDEX IF NOT EXISTS idx_structure_similarity_compound2 ON structure_similarity(compound_id_2);
CREATE INDEX IF NOT EXISTS idx_structure_similarity_metric ON structure_similarity(similarity_metric);

CREATE INDEX IF NOT EXISTS idx_sar_patterns_type ON sar_patterns(pattern_type);

CREATE INDEX IF NOT EXISTS idx_activity_cliffs_compounds ON activity_cliffs USING gin(compound_pair);
CREATE INDEX IF NOT EXISTS idx_activity_cliffs_type ON activity_cliffs(activity_type);

CREATE INDEX IF NOT EXISTS idx_scaffold_analysis_smiles ON scaffold_analysis(scaffold_smiles);

CREATE INDEX IF NOT EXISTS idx_pathway_analysis_compound ON pathway_analysis(compound_id);
CREATE INDEX IF NOT EXISTS idx_pathway_analysis_name ON pathway_analysis(pathway_name);
CREATE INDEX IF NOT EXISTS idx_pathway_analysis_type ON pathway_analysis(pathway_type);

CREATE INDEX IF NOT EXISTS idx_systems_biology_compound ON systems_biology_analysis(compound_id);
CREATE INDEX IF NOT EXISTS idx_systems_biology_level ON systems_biology_analysis(analysis_level);

CREATE INDEX IF NOT EXISTS idx_research_findings_compound ON research_findings(compound_id);
CREATE INDEX IF NOT EXISTS idx_research_findings_type ON research_findings(finding_type);

-- Add update triggers
DO $$ 
DECLARE
    t text;
BEGIN
    FOR t IN 
        SELECT table_name 
        FROM information_schema.tables 
        WHERE table_schema = 'public' 
        AND table_type = 'BASE TABLE'
        AND table_name IN (
            'sar_analysis',
            'activity_correlations', 
            'structure_similarity',
            'sar_patterns',
            'activity_cliffs',
            'scaffold_analysis',
            'pathway_analysis',
            'systems_biology_analysis',
            'research_findings'
        )
    LOOP
        EXECUTE format('
            CREATE TRIGGER update_%I_modtime 
            BEFORE UPDATE ON %I 
            FOR EACH ROW EXECUTE FUNCTION update_updated_at_column();
            
            CREATE TRIGGER audit_%I_trigger
            AFTER INSERT OR UPDATE OR DELETE ON %I
            FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();',
            t, t, t, t);
    END LOOP;
END $$;

-- Literature analysis tables
CREATE TABLE IF NOT EXISTS scientific_papers (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    doi text UNIQUE,
    title text NOT NULL,
    authors text[] NOT NULL,
    journal text,
    publication_date date,
    abstract text,
    full_text text,
    methodology text[],
    study_type text,                     -- e.g., 'Clinical Trial', 'Review'
    sample_size integer,
    study_duration interval,
    quality_metrics jsonb,               -- Study quality assessment
    evidence_level text,                 -- e.g., 'A1', 'B2'
    key_findings text[],
    limitations text[],
    compounds_studied uuid[],            -- References to compounds table
    validation_status text,              -- e.g., 'Validated', 'Pending', 'Disputed'
    peer_review_status text,             -- e.g., 'Peer-reviewed', 'Preprint'
    citation_count integer,
    impact_factor double precision,
    external_links jsonb,                -- Links to external databases
    supplementary_data jsonb,            -- Additional data and materials
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

CREATE TABLE IF NOT EXISTS literature_findings (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    paper_id uuid NOT NULL REFERENCES scientific_papers(id) ON DELETE CASCADE,
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    finding_type text NOT NULL,          -- e.g., 'Efficacy', 'Safety'
    finding_details text NOT NULL,
    statistical_significance double precision,
    confidence_interval jsonb,
    methodology_notes text,
    limitations text[],
    replication_status text,             -- e.g., 'Replicated', 'Not replicated'
    validation_method text[],
    supporting_evidence jsonb,
    contradicting_evidence jsonb,
    clinical_relevance text,
    research_implications text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

CREATE TABLE IF NOT EXISTS meta_analyses (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    title text NOT NULL,
    topic text NOT NULL,
    included_papers uuid[] NOT NULL,     -- References to scientific_papers
    methodology text NOT NULL,
    total_sample_size integer,
    pooled_effect_size double precision,
    heterogeneity_metrics jsonb,
    subgroup_analyses jsonb,
    sensitivity_analyses jsonb,
    publication_bias_assessment jsonb,
    quality_assessment_method text,
    evidence_strength text,              -- e.g., 'Strong', 'Moderate', 'Weak'
    clinical_implications text[],
    research_gaps text[],
    conclusions text,
    limitations text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Create indexes for literature analysis
CREATE INDEX IF NOT EXISTS idx_papers_doi ON scientific_papers(doi);
CREATE INDEX IF NOT EXISTS idx_papers_date ON scientific_papers(publication_date);
CREATE INDEX IF NOT EXISTS idx_papers_type ON scientific_papers(study_type);
CREATE INDEX IF NOT EXISTS idx_papers_evidence ON scientific_papers(evidence_level);
CREATE INDEX IF NOT EXISTS idx_papers_validation ON scientific_papers(validation_status);
CREATE INDEX IF NOT EXISTS idx_papers_compounds ON scientific_papers USING gin(compounds_studied);

CREATE INDEX IF NOT EXISTS idx_findings_paper ON literature_findings(paper_id);
CREATE INDEX IF NOT EXISTS idx_findings_compound ON literature_findings(compound_id);
CREATE INDEX IF NOT EXISTS idx_findings_type ON literature_findings(finding_type);
CREATE INDEX IF NOT EXISTS idx_findings_replication ON literature_findings(replication_status);

CREATE INDEX IF NOT EXISTS idx_meta_analyses_topic ON meta_analyses(topic);
CREATE INDEX IF NOT EXISTS idx_meta_analyses_evidence ON meta_analyses(evidence_strength);
CREATE INDEX IF NOT EXISTS idx_meta_analyses_papers ON meta_analyses USING gin(included_papers);

-- Add update triggers for literature analysis tables
CREATE TRIGGER update_scientific_papers_modtime
    BEFORE UPDATE ON scientific_papers
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_literature_findings_modtime
    BEFORE UPDATE ON literature_findings
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_meta_analyses_modtime
    BEFORE UPDATE ON meta_analyses
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Add audit triggers for literature analysis tables
CREATE TRIGGER audit_scientific_papers_trigger
    AFTER INSERT OR UPDATE OR DELETE ON scientific_papers
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_literature_findings_trigger
    AFTER INSERT OR UPDATE OR DELETE ON literature_findings
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_meta_analyses_trigger
    AFTER INSERT OR UPDATE OR DELETE ON meta_analyses
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Add comments
COMMENT ON TABLE sar_analysis IS 'Structure-activity relationship analysis results';
COMMENT ON TABLE activity_correlations IS 'Correlations between structural features and activities';
COMMENT ON TABLE structure_similarity IS 'Pairwise structural similarity between compounds';
COMMENT ON TABLE sar_patterns IS 'Identified structure-activity relationship patterns';
COMMENT ON TABLE activity_cliffs IS 'Activity cliff analysis between compound pairs';
COMMENT ON TABLE scaffold_analysis IS 'Analysis of molecular scaffolds and their properties';
COMMENT ON TABLE pathway_analysis IS 'Analysis of pathway effects and regulation';
COMMENT ON TABLE systems_biology_analysis IS 'Systems-level biological analysis';
COMMENT ON TABLE research_findings IS 'Research findings and conclusions';
COMMENT ON TABLE scientific_papers IS 'Scientific literature and research papers';
COMMENT ON TABLE literature_findings IS 'Specific findings from scientific literature';
COMMENT ON TABLE meta_analyses IS 'Meta-analyses of multiple research papers';
