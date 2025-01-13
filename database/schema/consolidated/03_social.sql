-- Consolidated social and community data schema

-- Core social content tables
CREATE TABLE social_posts (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    platform text NOT NULL, -- 'reddit', 'twitter', 'longecity', 'psychonautwiki', 'tripsit', etc.
    platform_subdivision text, -- e.g., subreddit, forum section, hashtag group
    external_id text NOT NULL,
    url text NOT NULL,
    title text,
    content text NOT NULL,
    author_id text NOT NULL,
    author_username text NOT NULL,
    post_type text NOT NULL, -- 'text', 'link', 'image', 'video', 'article', 'experience', 'scientific'
    created_at timestamptz NOT NULL,
    engagement_metrics jsonb NOT NULL DEFAULT '{}', -- {views, likes, replies, shares, etc.}
    classification_data jsonb NOT NULL DEFAULT '{}', -- {label, score, categories}
    sentiment_data jsonb NOT NULL DEFAULT '{}', -- {score, positive_ratio, negative_ratio, neutral_ratio}
    metadata jsonb NOT NULL DEFAULT '{}', -- Platform-specific metadata
    is_scientific boolean NOT NULL DEFAULT false,
    is_experience_report boolean NOT NULL DEFAULT false,
    is_harm_reduction boolean NOT NULL DEFAULT false,
    platform_specific_data jsonb NOT NULL DEFAULT '{}', -- Platform-specific fields
    tags text[] NOT NULL DEFAULT '{}',
    research_citations text[] NOT NULL DEFAULT '{}',
    view_count integer NOT NULL DEFAULT 0,
    created_at_internal timestamptz NOT NULL DEFAULT now(),
    UNIQUE(platform, external_id),
    CONSTRAINT check_engagement_metrics CHECK (
        CASE platform
            WHEN 'reddit' THEN (engagement_metrics->>'score')::integer >= 0
            WHEN 'twitter' THEN 
                (engagement_metrics->>'retweet_count')::integer >= 0 AND
                (engagement_metrics->>'favorite_count')::integer >= 0
            ELSE true
        END
    ),
    CONSTRAINT check_classification_scores CHECK (
        (classification_data->>'spam_score')::float BETWEEN 0 AND 1 AND
        (classification_data->>'toxicity_score')::float BETWEEN 0 AND 1
    ),
    CONSTRAINT check_sentiment_scores CHECK (
        (sentiment_data->>'score')::float BETWEEN -1 AND 1 AND
        (sentiment_data->>'positive_ratio')::float BETWEEN 0 AND 1 AND
        (sentiment_data->>'negative_ratio')::float BETWEEN 0 AND 1 AND
        (sentiment_data->>'neutral_ratio')::float BETWEEN 0 AND 1
    )
);

-- Comments and threaded discussions
CREATE TABLE social_comments (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    post_id uuid NOT NULL REFERENCES social_posts(id) ON DELETE CASCADE,
    platform text NOT NULL,
    platform_subdivision text,
    external_id text NOT NULL,
    parent_id text, -- External ID of parent comment
    content text NOT NULL,
    author_id text NOT NULL,
    author_username text NOT NULL,
    created_at timestamptz NOT NULL,
    engagement_metrics jsonb NOT NULL DEFAULT '{}', -- {score, replies, etc.}
    classification_data jsonb NOT NULL DEFAULT '{}',
    sentiment_data jsonb NOT NULL DEFAULT '{}',
    is_scientific boolean NOT NULL DEFAULT false,
    has_citations boolean NOT NULL DEFAULT false,
    reported_effects text[] NOT NULL DEFAULT '{}',
    reported_side_effects text[] NOT NULL DEFAULT '{}',
    platform_specific_data jsonb NOT NULL DEFAULT '{}',
    created_at_internal timestamptz NOT NULL DEFAULT now(),
    UNIQUE(platform, external_id),
    CONSTRAINT check_comment_engagement CHECK (
        CASE platform
            WHEN 'reddit' THEN (engagement_metrics->>'score')::integer >= 0
            WHEN 'twitter' THEN (engagement_metrics->>'likes')::integer >= 0
            ELSE true
        END
    ),
    CONSTRAINT check_comment_classification CHECK (
        (classification_data->>'spam_score')::float BETWEEN 0 AND 1 AND
        (classification_data->>'toxicity_score')::float BETWEEN 0 AND 1
    ),
    CONSTRAINT check_comment_sentiment CHECK (
        (sentiment_data->>'score')::float BETWEEN -1 AND 1 AND
        (sentiment_data->>'positive_ratio')::float BETWEEN 0 AND 1 AND
        (sentiment_data->>'negative_ratio')::float BETWEEN 0 AND 1 AND
        (sentiment_data->>'neutral_ratio')::float BETWEEN 0 AND 1
    )
);

-- Platform-specific statistics
CREATE TABLE social_platform_stats (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    platform text NOT NULL,
    platform_subdivision text, -- e.g., subreddit, forum section
    compound_id uuid REFERENCES compounds(id) ON DELETE CASCADE,
    total_posts integer NOT NULL DEFAULT 0,
    total_comments integer NOT NULL DEFAULT 0,
    unique_authors integer NOT NULL DEFAULT 0,
    scientific_post_ratio numeric,
    experience_report_ratio numeric,
    harm_reduction_ratio numeric,
    top_compounds text[],
    topic_distribution jsonb,
    sentiment_distribution jsonb,
    quality_metrics jsonb,
    engagement_metrics jsonb,
    last_updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(platform, platform_subdivision, compound_id)
);

-- Stack and protocol tracking
CREATE TABLE social_stacks (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    stack_name text NOT NULL,
    creator text,
    description text,
    purpose text, -- e.g., 'Memory', 'Focus', 'Neuroprotection'
    compounds uuid[] NOT NULL, -- References to compounds table
    dosages jsonb, -- Structured dosage information
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
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Protocol tracking
CREATE TABLE social_protocols (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    protocol_name text NOT NULL,
    creator text,
    description text NOT NULL,
    target_outcome text,
    compounds uuid[] NOT NULL, -- References to compounds table
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
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Enhanced experience reports
CREATE TABLE social_experience_reports (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    post_id uuid REFERENCES social_posts(id) ON DELETE CASCADE,
    compound_id uuid REFERENCES compounds(id) ON DELETE CASCADE,
    stack_id uuid REFERENCES social_stacks(id) ON DELETE SET NULL,
    protocol_id uuid REFERENCES social_protocols(id) ON DELETE SET NULL,
    title text NOT NULL,
    content text NOT NULL,
    author text,
    platform text NOT NULL,
    platform_subdivision text,
    experience_date date,
    report_date timestamptz NOT NULL,
    substance_data jsonb[], -- Array of {substance, dose, roa, timing}
    duration interval,
    setting_data jsonb,
    intention text,
    effects_timeline jsonb[], -- Detailed timeline with intensity tracking
    timeline_data jsonb, -- Additional timeline metadata
    reported_effects text[],
    side_effects text[],
    interactions text[],
    test_kit_info jsonb, -- Test kit results and methods
    detection_time_data jsonb, -- Detection windows by test type
    after_effects_data jsonb, -- Structured after-effects tracking
    body_weight numeric,
    weight_unit text,
    gender text,
    age integer,
    experience_level text,
    harm_reduction_notes text[],
    classification_data jsonb,
    sentiment_data jsonb,
    report_version integer,
    experience_category text[], -- e.g., ['First Time', 'Difficult Experience', 'Medical Use']
    total_views integer,
    report_quality_score double precision,
    medical_conditions text[],
    medications text[],
    baseline_metrics jsonb,
    outcome_metrics jsonb,
    testing_methods text[],
    overall_rating integer,
    platform_specific_data jsonb,
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Experience categories
CREATE TABLE experience_categories (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    category_name text NOT NULL UNIQUE,
    parent_category text REFERENCES experience_categories(category_name),
    description text,
    report_count integer NOT NULL DEFAULT 0,
    created_at timestamptz NOT NULL DEFAULT now()
);


-- Scientific content
CREATE TABLE social_scientific_content (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    post_id uuid REFERENCES social_posts(id) ON DELETE CASCADE,
    compound_id uuid REFERENCES compounds(id) ON DELETE CASCADE,
    title text NOT NULL,
    content text NOT NULL,
    author text,
    platform text,
    platform_subdivision text,
    publication_date timestamptz NOT NULL,
    content_type text NOT NULL, -- 'research', 'review', 'analysis', 'protocol'
    research_topics text[],
    methodology text,
    findings text[],
    limitations text[],
    citations text[],
    peer_review_notes text[],
    quality_score numeric,
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Research reviews
CREATE TABLE social_research_reviews (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    title text NOT NULL,
    content text NOT NULL,
    author text,
    publication_date timestamptz NOT NULL,
    study_type text[], -- e.g., ['Clinical Trial', 'Meta-Analysis']
    methodology text,
    findings text[],
    limitations text[],
    research_quality_score numeric,
    citations text[],
    peer_review_notes text[],
    platform text NOT NULL,
    platform_subdivision text,
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Harm reduction content
CREATE TABLE social_harm_reduction (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    post_id uuid REFERENCES social_posts(id) ON DELETE CASCADE,
    compound_id uuid REFERENCES compounds(id) ON DELETE CASCADE,
    title text NOT NULL,
    content text NOT NULL,
    author text,
    platform text,
    platform_subdivision text,
    created_at timestamptz NOT NULL DEFAULT now(),
    category text NOT NULL, -- 'Emergency', 'General', 'Dosage', 'ROA', 'Combinations', etc.
    importance_level text, -- 'Critical', 'Important', 'Helpful'
    safety_notes text[],
    warnings text[],
    contraindications text[],
    emergency_procedures text[],
    sources text[],
    verification_status text NOT NULL DEFAULT 'unverified',
    verification_notes text,
    last_reviewed_at timestamptz,
    created_at_internal timestamptz NOT NULL DEFAULT now()
);

-- Compound combinations
CREATE TABLE social_compound_combinations (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id_1 uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    compound_id_2 uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    interaction_type text NOT NULL, -- 'dangerous', 'unsafe', 'caution', 'safe', 'synergy'
    risk_level text NOT NULL,
    description text,
    mechanism text,
    evidence_level text,
    sources text[],
    platforms text[] NOT NULL,
    reported_count integer NOT NULL DEFAULT 0,
    reports uuid[], -- References to experience_reports
    created_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(compound_id_1, compound_id_2)
);

-- Dosage and administration data
CREATE TABLE social_dosage_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
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
    created_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(compound_id, platform, route_of_administration)
);

-- Effect definitions and categorization
CREATE TABLE social_effects (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    effect_name text NOT NULL UNIQUE,
    description text,
    type text, -- e.g., 'Cognitive', 'Physical', 'Visual', 'Auditory'
    url text,
    analysis_data jsonb, -- Additional structured data about the effect
    related_effects text[],
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Effect-specific experience reports
CREATE TABLE social_effect_reports (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    effect_id uuid NOT NULL REFERENCES social_effects(id) ON DELETE CASCADE,
    report_text text NOT NULL,
    intensity integer, -- 1-5 scale
    duration interval,
    onset_time interval,
    conditions jsonb, -- Context in which effect was experienced
    platform text NOT NULL,
    platform_subdivision text,
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Enhanced substance data
CREATE TABLE social_substance_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    platform text NOT NULL,
    platform_subdivision text,
    external_id text NOT NULL,
    common_names text[],
    chemical_class text[],
    psychoactive_class text[],
    summary text,
    tolerance_data jsonb, -- {development_time, full_tolerance, baseline_return}
    roa_data jsonb, -- Routes of administration with onset, duration, dosage by route
    effect_data jsonb, -- Structured effects information
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
    dosage_data jsonb, -- Structured dosage guidelines
    duration_data jsonb, -- Structured duration information
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
    last_updated_at timestamptz,
    created_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(platform, external_id)
);

-- Cross-platform analytics
CREATE TABLE social_compound_mentions (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    platform text NOT NULL,
    platform_subdivision text,
    analysis_date date NOT NULL,
    total_mentions integer NOT NULL DEFAULT 0,
    unique_authors integer NOT NULL DEFAULT 0,
    total_engagement integer NOT NULL DEFAULT 0,
    scientific_mentions integer NOT NULL DEFAULT 0,
    experience_reports integer NOT NULL DEFAULT 0,
    harm_reduction_mentions integer NOT NULL DEFAULT 0,
    sentiment_distribution jsonb,
    topic_distribution jsonb,
    user_demographics jsonb,
    geographic_distribution jsonb,
    temporal_patterns jsonb,
    common_contexts text[],
    related_compounds text[],
    platform_specific_metrics jsonb, -- Platform-specific metrics (e.g., Reddit: {subreddit_stats, flair_stats}, Twitter: {hashtag_stats, retweet_patterns})
    cross_platform_engagement_flow jsonb, -- How content spreads across platforms
    user_influence_metrics jsonb, -- User influence scores across platforms
    content_propagation_patterns jsonb, -- How content evolves as it spreads
    created_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(compound_id, platform, platform_subdivision, analysis_date)
);

-- Trend analysis
CREATE TABLE social_compound_trends (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    trend_start_date date NOT NULL,
    trend_end_date date NOT NULL,
    platforms text[] NOT NULL,
    total_mentions integer NOT NULL DEFAULT 0,
    unique_authors integer NOT NULL DEFAULT 0,
    total_engagement integer NOT NULL DEFAULT 0,
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
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Content quality metrics
CREATE TABLE social_content_quality (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
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
    created_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(compound_id, platform, platform_subdivision, analysis_date)
);

-- Influence network analysis
CREATE TABLE social_influence_networks (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
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
    created_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(compound_id, platform, platform_subdivision, analysis_date)
);

-- Alert rules and monitoring
CREATE TABLE social_alert_rules (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL,
    description text,
    alert_type text NOT NULL, -- 'safety', 'trend', 'scientific', 'regulatory'
    severity text NOT NULL, -- 'critical', 'high', 'medium', 'low'
    platforms text[] NOT NULL,
    platform_subdivisions text[],
    compounds uuid[],
    trigger_conditions jsonb NOT NULL,
    required_metrics text[] NOT NULL,
    threshold_values jsonb NOT NULL,
    cooldown_period interval,
    is_active boolean NOT NULL DEFAULT true,
    last_triggered_at timestamptz,
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Safety incidents and alerts
CREATE TABLE social_safety_incidents (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    incident_type text NOT NULL, -- 'adverse_effect', 'interaction', 'contamination'
    severity text NOT NULL,
    description text NOT NULL,
    reported_effects text[],
    reported_causes text[],
    platforms text[] NOT NULL,
    platform_subdivisions text[],
    content_urls text[],
    verification_status text NOT NULL DEFAULT 'unverified',
    verification_notes text,
    response_actions text[],
    resolution_status text NOT NULL DEFAULT 'open',
    resolution_notes text,
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Create indexes
CREATE INDEX idx_social_posts_compound ON social_posts(compound_id);
CREATE INDEX idx_social_posts_platform ON social_posts(platform);
CREATE INDEX idx_social_posts_subdivision ON social_posts(platform_subdivision);
CREATE INDEX idx_social_posts_external_id ON social_posts(external_id);
CREATE INDEX idx_social_posts_created ON social_posts(created_at);
CREATE INDEX idx_social_posts_type ON social_posts(post_type);
CREATE INDEX idx_social_posts_scientific ON social_posts(is_scientific);

-- Content search indexes
CREATE INDEX idx_social_posts_content_search ON social_posts USING gin(to_tsvector('english', content));
CREATE INDEX idx_reddit_content_search ON social_posts USING gin(to_tsvector('english', content)) WHERE platform = 'reddit';
CREATE INDEX idx_twitter_content_search ON social_posts USING gin(to_tsvector('english', content)) WHERE platform = 'twitter';
CREATE INDEX idx_community_content_search ON social_posts USING gin(to_tsvector('english', content));

-- Engagement metric indexes
CREATE INDEX idx_reddit_score ON social_posts ((engagement_metrics->>'score')) WHERE platform = 'reddit';
CREATE INDEX idx_twitter_engagement ON social_posts ((engagement_metrics->>'retweet_count')) WHERE platform = 'twitter';

-- Classification metric indexes
CREATE INDEX idx_community_spam ON social_posts ((classification_data->>'spam_score'));
CREATE INDEX idx_community_toxicity ON social_posts ((classification_data->>'toxicity_score'));

CREATE INDEX idx_social_comments_post ON social_comments(post_id);
CREATE INDEX idx_social_comments_platform ON social_comments(platform);
CREATE INDEX idx_social_comments_subdivision ON social_comments(platform_subdivision);
CREATE INDEX idx_social_comments_external ON social_comments(external_id);
CREATE INDEX idx_social_comments_parent ON social_comments(parent_id);
CREATE INDEX idx_social_comments_created ON social_comments(created_at);
CREATE INDEX idx_social_comments_scientific ON social_comments(is_scientific);

CREATE INDEX idx_platform_stats_platform ON social_platform_stats(platform);
CREATE INDEX idx_platform_stats_subdivision ON social_platform_stats(platform_subdivision);
CREATE INDEX idx_platform_stats_compound ON social_platform_stats(compound_id);
CREATE INDEX idx_platform_stats_updated ON social_platform_stats(last_updated_at);

CREATE INDEX idx_experience_reports_post ON social_experience_reports(post_id);
CREATE INDEX idx_experience_reports_compound ON social_experience_reports(compound_id);
CREATE INDEX idx_experience_reports_platform ON social_experience_reports(platform);
CREATE INDEX idx_experience_reports_subdivision ON social_experience_reports(platform_subdivision);
CREATE INDEX idx_experience_reports_date ON social_experience_reports(experience_date);
CREATE INDEX idx_experience_reports_stack ON social_experience_reports(stack_id);
CREATE INDEX idx_experience_reports_protocol ON social_experience_reports(protocol_id);

CREATE INDEX idx_experience_categories_name ON experience_categories(category_name);
CREATE INDEX idx_experience_categories_parent ON experience_categories(parent_category);

CREATE INDEX idx_social_stacks_compounds ON social_stacks USING gin(compounds);
CREATE INDEX idx_social_stacks_purpose ON social_stacks(purpose);
CREATE INDEX idx_social_stacks_rating ON social_stacks(rating);
CREATE INDEX idx_social_stacks_platform ON social_stacks(platform);
CREATE INDEX idx_social_stacks_subdivision ON social_stacks(platform_subdivision);

CREATE INDEX idx_social_protocols_compounds ON social_protocols USING gin(compounds);
CREATE INDEX idx_social_protocols_outcome ON social_protocols(target_outcome);
CREATE INDEX idx_social_protocols_score ON social_protocols(review_score);
CREATE INDEX idx_social_protocols_platform ON social_protocols(platform);
CREATE INDEX idx_social_protocols_subdivision ON social_protocols(platform_subdivision);

CREATE INDEX idx_scientific_content_post ON social_scientific_content(post_id);
CREATE INDEX idx_scientific_content_compound ON social_scientific_content(compound_id);
CREATE INDEX idx_scientific_content_platform ON social_scientific_content(platform);
CREATE INDEX idx_scientific_content_subdivision ON social_scientific_content(platform_subdivision);
CREATE INDEX idx_scientific_content_type ON social_scientific_content(content_type);
CREATE INDEX idx_scientific_content_quality ON social_scientific_content(quality_score);

CREATE INDEX idx_social_research_compound ON social_research_reviews(compound_id);
CREATE INDEX idx_social_research_date ON social_research_reviews(publication_date);
CREATE INDEX idx_social_research_type ON social_research_reviews USING gin(study_type);
CREATE INDEX idx_social_research_quality ON social_research_reviews(research_quality_score);
CREATE INDEX idx_social_research_platform ON social_research_reviews(platform);
CREATE INDEX idx_social_research_subdivision ON social_research_reviews(platform_subdivision);

CREATE INDEX idx_harm_reduction_post ON social_harm_reduction(post_id);
CREATE INDEX idx_harm_reduction_compound ON social_harm_reduction(compound_id);
CREATE INDEX idx_harm_reduction_platform ON social_harm_reduction(platform);
CREATE INDEX idx_harm_reduction_subdivision ON social_harm_reduction(platform_subdivision);
CREATE INDEX idx_harm_reduction_category ON social_harm_reduction(category);
CREATE INDEX idx_harm_reduction_status ON social_harm_reduction(verification_status);

CREATE INDEX idx_combinations_compounds ON social_compound_combinations(compound_id_1, compound_id_2);
CREATE INDEX idx_combinations_type ON social_compound_combinations(interaction_type);
CREATE INDEX idx_combinations_risk ON social_compound_combinations(risk_level);
CREATE INDEX idx_combinations_platforms ON social_compound_combinations USING gin(platforms);

CREATE INDEX idx_dosage_compound ON social_dosage_data(compound_id);
CREATE INDEX idx_dosage_platform ON social_dosage_data(platform);
CREATE INDEX idx_dosage_subdivision ON social_dosage_data(platform_subdivision);
CREATE INDEX idx_dosage_roa ON social_dosage_data(route_of_administration);

CREATE INDEX idx_effects_name ON social_effects(effect_name);
CREATE INDEX idx_effects_type ON social_effects(type);

CREATE INDEX idx_effect_reports_compound ON social_effect_reports(compound_id);
CREATE INDEX idx_effect_reports_effect ON social_effect_reports(effect_id);
CREATE INDEX idx_effect_reports_intensity ON social_effect_reports(intensity);
CREATE INDEX idx_effect_reports_platform ON social_effect_reports(platform);
CREATE INDEX idx_effect_reports_subdivision ON social_effect_reports(platform_subdivision);

CREATE INDEX idx_substance_data_compound ON social_substance_data(compound_id);
CREATE INDEX idx_substance_data_platform ON social_substance_data(platform);
CREATE INDEX idx_substance_data_subdivision ON social_substance_data(platform_subdivision);
CREATE INDEX idx_substance_data_external ON social_substance_data(external_id);
CREATE INDEX idx_substance_data_names ON social_substance_data USING gin(common_names);
CREATE INDEX idx_substance_data_class ON social_substance_data USING gin(chemical_class, psychoactive_class);

CREATE INDEX idx_mentions_compound ON social_compound_mentions(compound_id);
CREATE INDEX idx_mentions_platform ON social_compound_mentions(platform);
CREATE INDEX idx_mentions_subdivision ON social_compound_mentions(platform_subdivision);
CREATE INDEX idx_mentions_date ON social_compound_mentions(analysis_date);

CREATE INDEX idx_trends_compound ON social_compound_trends(compound_id);
CREATE INDEX idx_trends_dates ON social_compound_trends(trend_start_date, trend_end_date);
CREATE INDEX idx_trends_platforms ON social_compound_trends USING gin(platforms);
CREATE INDEX idx_trends_correlation ON social_compound_trends(correlation_strength);
CREATE INDEX idx_trends_similarity ON social_compound_trends(temporal_similarity);

CREATE INDEX idx_quality_compound ON social_content_quality(compound_id);
CREATE INDEX idx_quality_platform ON social_content_quality(platform);
CREATE INDEX idx_quality_subdivision ON social_content_quality(platform_subdivision);
CREATE INDEX idx_quality_date ON social_content_quality(analysis_date);

CREATE INDEX idx_influence_compound ON social_influence_networks(compound_id);
CREATE INDEX idx_influence_platform ON social_influence_networks(platform);
CREATE INDEX idx_influence_subdivision ON social_influence_networks(platform_subdivision);
CREATE INDEX idx_influence_date ON social_influence_networks(analysis_date);

CREATE INDEX idx_alert_rules_type ON social_alert_rules(alert_type);
CREATE INDEX idx_alert_rules_severity ON social_alert_rules(severity);
CREATE INDEX idx_alert_rules_platforms ON social_alert_rules USING gin(platforms);
CREATE INDEX idx_alert_rules_subdivisions ON social_alert_rules USING gin(platform_subdivisions);
CREATE INDEX idx_alert_rules_compounds ON social_alert_rules USING gin(compounds);
CREATE INDEX idx_alert_rules_active ON social_alert_rules(is_active);

CREATE INDEX idx_safety_incidents_compound ON social_safety_incidents(compound_id);
CREATE INDEX idx_safety_incidents_type ON social_safety_incidents(incident_type);
CREATE INDEX idx_safety_incidents_severity ON social_safety_incidents(severity);
CREATE INDEX idx_safety_incidents_platforms ON social_safety_incidents USING gin(platforms);
CREATE INDEX idx_safety_incidents_subdivisions ON social_safety_incidents USING gin(platform_subdivisions);
CREATE INDEX idx_safety_incidents_status ON social_safety_incidents(verification_status);
CREATE INDEX idx_safety_incidents_resolution ON social_safety_incidents(resolution_status);

-- Create audit triggers
DO $$
DECLARE
    t text;
BEGIN
    FOR t IN SELECT table_name
             FROM information_schema.tables
             WHERE table_schema = 'public'
             AND table_type = 'BASE TABLE'
             AND table_name LIKE 'social_%'
    LOOP
        EXECUTE format('
            CREATE TRIGGER audit_%I_trigger
            AFTER INSERT OR UPDATE OR DELETE ON %I
            FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();',
            t, t);
    END LOOP;
END $$;

-- Add predefined experience categories
INSERT INTO experience_categories (category_name, description) VALUES
('First Time', 'First experiences with a substance'),
('Difficult Experience', 'Challenging or negative experiences'),
('Medical Use', 'Use for medical or therapeutic purposes'),
('Spiritual', 'Experiences with spiritual or mystical elements'),
('Train Wrecks & Trip Disasters', 'Particularly difficult or dangerous experiences'),
('Health Problems', 'Experiences involving health issues'),
('Retrospective / Summary', 'Looking back on past experiences'),
('Small Collection', 'Brief or collected experiences'),
('Master/Teacher Plants', 'Experiences with traditional plant medicines'),
('Combinations', 'Experiences with multiple substances'),
('Preparation / Recipes', 'Methods of preparation or consumption'),
('Bad Trips', 'Specifically difficult psychological experiences'),
('Glowing Experiences', 'Particularly positive experiences'),
('Clinical Research', 'Experiences in research settings');

-- Add predefined effects
INSERT INTO social_effects (effect_name, type, description) VALUES
('Visual Drifting', 'Visual', 'The experience of textures, surfaces, and objects appearing to move or flow'),
('Geometric Patterns', 'Visual', 'The experience of seeing various geometric patterns and forms'),
('Time Distortion', 'Cognitive', 'Alterations in the perception of time passing'),
('Euphoria', 'Physical/Cognitive', 'A state of intense happiness and well-being'),
('Enhanced Music Appreciation', 'Auditory', 'Music sounds more detailed, meaningful, or emotionally impactful'),
('Ego Dissolution', 'Cognitive', 'The experience of a decreased sense of self-identity'),
('Synesthesia', 'Cognitive', 'The mixing of sensory modalities'),
('Enhanced Tactile Sensation', 'Physical', 'Increased sensitivity to physical touch and textures'),
('Conceptual Thinking', 'Cognitive', 'Abstract thoughts become more vivid and meaningful'),
('Visual Acuity Enhancement', 'Visual', 'Improved clarity and sharpness of vision');

-- Add predefined harm reduction guidelines
INSERT INTO social_harm_reduction (compound_id, category, title, content, importance_level) VALUES
(NULL, 'Emergency', 'When to Call Emergency Services',
'Call emergency services immediately if someone experiences: severe confusion, unconsciousness, difficulty breathing, seizures, severe overheating, or chest pain.',
'Critical'),
(NULL, 'General', 'Test Your Substances',
'Always test your substances with multiple reagent tests. Never consume unidentified substances.',
'Critical'),
(NULL, 'Dosage', 'Start Low, Go Slow',
'Always start with a low dose, especially with new substances or batches. Wait sufficient time before considering redosing.',
'Critical'),
(NULL, 'ROA', 'Safe Injection Practices',
'Use clean equipment, never share needles, and practice proper hygiene. Know the proper injection techniques for harm reduction.',
'Critical'),
(NULL, 'Combinations', 'Avoid Dangerous Combinations',
'Research interactions before combining substances. Many combinations can be unexpectedly dangerous.',
'Critical');
