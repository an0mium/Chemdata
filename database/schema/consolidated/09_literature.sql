-- Consolidated literature schema
-- Combines core literature tables and enhanced research functionality

-- Publications and their metadata
CREATE TABLE IF NOT EXISTS publications (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    doi text UNIQUE,
    pubmed_id text UNIQUE,
    title text NOT NULL,
    authors text[] NOT NULL,
    journal text,
    publication_date date,
    abstract text,
    full_text text,
    keywords text[],
    publication_type text,
    study_type text,
    sample_size integer,
    study_duration interval,
    methodology text[],
    evidence_level text,
    quality_metrics jsonb,
    language text,
    full_text_url text,
    pdf_path text,
    citation_count integer,
    impact_factor double precision,
    source_database text,
    compounds_studied uuid[], -- References to compounds table
    metadata jsonb,
    text_search_vector tsvector GENERATED ALWAYS AS (
        setweight(to_tsvector('english', coalesce(title, '')), 'A') ||
        setweight(to_tsvector('english', coalesce(abstract, '')), 'B') ||
        setweight(to_tsvector('english', array_to_string(coalesce(keywords, ARRAY[]::text[]), ' ')), 'C')
    ) STORED,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Compound literature associations
CREATE TABLE IF NOT EXISTS compound_publications (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    publication_id uuid NOT NULL REFERENCES publications(id) ON DELETE CASCADE,
    association_type text NOT NULL, -- e.g., 'primary_focus', 'mentioned', 'compared'
    context_summary text,
    relevance_score double precision,
    extracted_data jsonb, -- structured data extracted from the publication about this compound
    annotation_status text,
    reviewer text,
    review_date date,
    notes text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(compound_id, publication_id)
);

-- Literature citations
CREATE TABLE IF NOT EXISTS citations (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    citing_publication_id uuid NOT NULL REFERENCES publications(id) ON DELETE CASCADE,
    cited_publication_id uuid NOT NULL REFERENCES publications(id) ON DELETE CASCADE,
    citation_context text,
    citation_type text, -- e.g., 'methodology', 'results', 'discussion'
    sentiment text, -- e.g., 'supports', 'contradicts', 'neutral'
    importance_score double precision,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(citing_publication_id, cited_publication_id)
);

-- Literature reviews
CREATE TABLE IF NOT EXISTS literature_reviews (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    title text NOT NULL,
    topic text NOT NULL,
    review_type text NOT NULL, -- e.g., 'systematic', 'narrative', 'meta-analysis'
    authors text[] NOT NULL,
    review_date date NOT NULL,
    search_criteria jsonb,
    included_studies uuid[] REFERENCES publications(id),
    excluded_studies jsonb, -- includes studies and reasons for exclusion
    total_sample_size integer,
    pooled_effect_size double precision,
    heterogeneity_metrics jsonb,
    methodology text,
    findings text,
    conclusions text,
    limitations text,
    recommendations text,
    status text,
    reviewer text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Literature analysis
CREATE TABLE IF NOT EXISTS literature_analysis (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    publication_id uuid NOT NULL REFERENCES publications(id) ON DELETE CASCADE,
    analysis_type text NOT NULL,
    analysis_date date NOT NULL,
    methodology_assessment text,
    study_design_score double precision,
    methodology_score double precision,
    sample_size_adequacy double precision,
    results_assessment text,
    statistical_rigor text,
    reproducibility_assessment text,
    confounding_control text,
    replication_status text,
    reporting_quality double precision,
    bias_assessment jsonb,
    quality_score double precision,
    strength_of_evidence text,
    quality_of_evidence text,
    consistency_score double precision,
    applicability_score double precision,
    limitations_noted text[],
    reviewer_comments text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Key findings
CREATE TABLE IF NOT EXISTS publication_findings (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    publication_id uuid NOT NULL REFERENCES publications(id) ON DELETE CASCADE,
    finding_type text NOT NULL,
    finding_summary text NOT NULL,
    methodology text,
    statistical_significance double precision,
    confidence_interval jsonb,
    confidence_level text,
    limitations text[],
    implications text,
    supporting_data jsonb,
    reviewer_notes text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Literature topics
CREATE TABLE IF NOT EXISTS literature_topics (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    topic_name text NOT NULL UNIQUE,
    description text,
    parent_topic uuid REFERENCES literature_topics(id),
    keywords text[],
    related_topics uuid[] REFERENCES literature_topics(id),
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Publication topic associations
CREATE TABLE IF NOT EXISTS publication_topics (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    publication_id uuid NOT NULL REFERENCES publications(id) ON DELETE CASCADE,
    topic_id uuid NOT NULL REFERENCES literature_topics(id) ON DELETE CASCADE,
    relevance_score double precision,
    context_summary text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(publication_id, topic_id)
);

-- Literature trends
CREATE TABLE IF NOT EXISTS literature_trends (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    topic_id uuid REFERENCES literature_topics(id) ON DELETE CASCADE,
    trend_period_start date NOT NULL,
    trend_period_end date NOT NULL,
    publication_count integer,
    citation_trends jsonb,
    key_developments text[],
    emerging_themes text[],
    trend_analysis text,
    future_directions text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Author information
CREATE TABLE IF NOT EXISTS authors (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL,
    affiliations text[],
    orcid_id text UNIQUE,
    research_areas text[],
    publication_count integer,
    h_index integer,
    citation_count integer,
    contact_info jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Publication author associations
CREATE TABLE IF NOT EXISTS publication_authors (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    publication_id uuid NOT NULL REFERENCES publications(id) ON DELETE CASCADE,
    author_id uuid NOT NULL REFERENCES authors(id) ON DELETE CASCADE,
    author_position integer,
    contribution_type text[], -- e.g., ['conceptualization', 'methodology', 'writing']
    corresponding_author boolean,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(publication_id, author_id)
);

-- Enhanced research functionality
CREATE TABLE IF NOT EXISTS research_projects (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    title text NOT NULL,
    description text,
    start_date date,
    end_date date,
    status text,
    methodology text,
    objectives text[],
    key_findings text[],
    related_publications uuid[] REFERENCES publications(id),
    team_members uuid[] REFERENCES authors(id),
    funding_sources text[],
    budget_info jsonb,
    resources_used text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Create indexes
CREATE INDEX IF NOT EXISTS idx_publications_doi ON publications(doi);
CREATE INDEX IF NOT EXISTS idx_publications_pubmed ON publications(pubmed_id);
CREATE INDEX IF NOT EXISTS idx_publications_date ON publications(publication_date);
CREATE INDEX IF NOT EXISTS idx_publications_type ON publications(publication_type);
CREATE INDEX IF NOT EXISTS idx_publications_journal ON publications(journal);
CREATE INDEX IF NOT EXISTS idx_publications_study_type ON publications(study_type);
CREATE INDEX IF NOT EXISTS idx_publications_evidence ON publications(evidence_level);
CREATE INDEX IF NOT EXISTS idx_publications_compounds ON publications USING gin(compounds_studied);
CREATE INDEX IF NOT EXISTS idx_publications_search ON publications USING gin(text_search_vector);

CREATE INDEX IF NOT EXISTS idx_compound_publications_compound ON compound_publications(compound_id);
CREATE INDEX IF NOT EXISTS idx_compound_publications_publication ON compound_publications(publication_id);
CREATE INDEX IF NOT EXISTS idx_compound_publications_type ON compound_publications(association_type);
CREATE INDEX IF NOT EXISTS idx_compound_publications_status ON compound_publications(annotation_status);

CREATE INDEX IF NOT EXISTS idx_citations_citing ON citations(citing_publication_id);
CREATE INDEX IF NOT EXISTS idx_citations_cited ON citations(cited_publication_id);
CREATE INDEX IF NOT EXISTS idx_citations_type ON citations(citation_type);
CREATE INDEX IF NOT EXISTS idx_citations_sentiment ON citations(sentiment);

CREATE INDEX IF NOT EXISTS idx_literature_reviews_topic ON literature_reviews(topic);
CREATE INDEX IF NOT EXISTS idx_literature_reviews_type ON literature_reviews(review_type);
CREATE INDEX IF NOT EXISTS idx_literature_reviews_date ON literature_reviews(review_date);
CREATE INDEX IF NOT EXISTS idx_literature_reviews_status ON literature_reviews(status);

CREATE INDEX IF NOT EXISTS idx_literature_analysis_publication ON literature_analysis(publication_id);
CREATE INDEX IF NOT EXISTS idx_literature_analysis_type ON literature_analysis(analysis_type);
CREATE INDEX IF NOT EXISTS idx_literature_analysis_date ON literature_analysis(analysis_date);
CREATE INDEX IF NOT EXISTS idx_literature_analysis_quality ON literature_analysis(quality_score);
CREATE INDEX IF NOT EXISTS idx_literature_analysis_replication ON literature_analysis(replication_status);
CREATE INDEX IF NOT EXISTS idx_literature_analysis_consistency ON literature_analysis(consistency_score);
CREATE INDEX IF NOT EXISTS idx_literature_analysis_applicability ON literature_analysis(applicability_score);

CREATE INDEX IF NOT EXISTS idx_publication_findings_publication ON publication_findings(publication_id);
CREATE INDEX IF NOT EXISTS idx_publication_findings_type ON publication_findings(finding_type);

CREATE INDEX IF NOT EXISTS idx_literature_topics_name ON literature_topics(topic_name);
CREATE INDEX IF NOT EXISTS idx_literature_topics_parent ON literature_topics(parent_topic);

CREATE INDEX IF NOT EXISTS idx_publication_topics_publication ON publication_topics(publication_id);
CREATE INDEX IF NOT EXISTS idx_publication_topics_topic ON publication_topics(topic_id);

CREATE INDEX IF NOT EXISTS idx_literature_trends_topic ON literature_trends(topic_id);
CREATE INDEX IF NOT EXISTS idx_literature_trends_period ON literature_trends(trend_period_start, trend_period_end);

CREATE INDEX IF NOT EXISTS idx_authors_name ON authors(name);
CREATE INDEX IF NOT EXISTS idx_authors_orcid ON authors(orcid_id);

CREATE INDEX IF NOT EXISTS idx_publication_authors_publication ON publication_authors(publication_id);
CREATE INDEX IF NOT EXISTS idx_publication_authors_author ON publication_authors(author_id);
CREATE INDEX IF NOT EXISTS idx_publication_authors_position ON publication_authors(author_position);
CREATE INDEX IF NOT EXISTS idx_publication_authors_corresponding ON publication_authors(corresponding_author);

CREATE INDEX IF NOT EXISTS idx_research_projects_title ON research_projects(title);
CREATE INDEX IF NOT EXISTS idx_research_projects_status ON research_projects(status);
CREATE INDEX IF NOT EXISTS idx_research_projects_dates ON research_projects(start_date, end_date);

-- Add triggers for timestamp updates
DO $$ 
DECLARE
    t text;
BEGIN
    FOR t IN SELECT table_name 
             FROM information_schema.tables 
             WHERE table_schema = 'public' 
             AND table_type = 'BASE TABLE' 
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

-- Add comments
COMMENT ON TABLE publications IS 'Scientific publications and their metadata';
COMMENT ON TABLE compound_publications IS 'Associations between compounds and publications';
COMMENT ON TABLE citations IS 'Citation relationships between publications';
COMMENT ON TABLE literature_reviews IS 'Systematic and narrative literature reviews';
COMMENT ON TABLE literature_analysis IS 'Analysis and assessment of publications';
COMMENT ON TABLE publication_findings IS 'Key findings extracted from publications';
COMMENT ON TABLE literature_topics IS 'Hierarchical organization of literature topics';
COMMENT ON TABLE publication_topics IS 'Associations between publications and topics';
COMMENT ON TABLE literature_trends IS 'Analysis of trends in scientific literature';
COMMENT ON TABLE authors IS 'Author information and metrics';
COMMENT ON TABLE publication_authors IS 'Associations between publications and authors';
COMMENT ON TABLE research_projects IS 'Research project management and tracking';
