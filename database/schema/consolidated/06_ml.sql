-- Consolidated machine learning schema
-- Combines model definitions, predictions, and features

-- ML models
CREATE TABLE IF NOT EXISTS ml_models (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
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
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(name, version)
);

-- Model versions
CREATE TABLE IF NOT EXISTS model_versions (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    model_id uuid NOT NULL REFERENCES ml_models(id) ON DELETE CASCADE,
    version_number text NOT NULL,
    changes_description text,
    performance_metrics jsonb,
    validation_results jsonb,
    deployment_status text,
    deployed_at timestamptz,
    deprecated_at timestamptz,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(model_id, version_number)
);

-- Training history
CREATE TABLE IF NOT EXISTS training_history (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    model_id uuid NOT NULL REFERENCES ml_models(id) ON DELETE CASCADE,
    version_id uuid NOT NULL REFERENCES model_versions(id),
    training_run_id text NOT NULL,
    start_time timestamptz NOT NULL,
    end_time timestamptz,
    parameters jsonb,
    metrics jsonb,
    loss_history jsonb,
    validation_history jsonb,
    hardware_metrics jsonb,
    status text NOT NULL,
    error_logs text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(model_id, training_run_id)
);

-- Model predictions
CREATE TABLE IF NOT EXISTS model_predictions (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    model_id uuid NOT NULL REFERENCES ml_models(id),
    version_id uuid NOT NULL REFERENCES model_versions(id),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    prediction_type text NOT NULL,
    predicted_value jsonb NOT NULL,
    confidence_score double precision,
    prediction_date timestamptz NOT NULL,
    input_features jsonb,
    explanation jsonb,
    metadata jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Feature definitions
CREATE TABLE IF NOT EXISTS feature_definitions (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    description text,
    feature_type text NOT NULL,
    data_type text NOT NULL,
    calculation_method text,
    dependencies text[],
    validation_rules jsonb,
    metadata jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Feature values
CREATE TABLE IF NOT EXISTS feature_values (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    feature_id uuid NOT NULL REFERENCES feature_definitions(id),
    value jsonb NOT NULL,
    calculation_date timestamptz NOT NULL,
    confidence_score double precision,
    metadata jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(compound_id, feature_id)
);

-- Feature engineering pipelines
CREATE TABLE IF NOT EXISTS feature_pipelines (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    description text,
    steps jsonb NOT NULL,
    input_features text[],
    output_features text[],
    parameters jsonb,
    validation_rules jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Feature selection results
CREATE TABLE IF NOT EXISTS feature_selection (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    model_id uuid NOT NULL REFERENCES ml_models(id),
    selection_method text NOT NULL,
    selected_features text[] NOT NULL,
    importance_scores jsonb,
    selection_criteria jsonb,
    validation_metrics jsonb,
    selection_date timestamptz NOT NULL,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Model performance metrics
CREATE TABLE IF NOT EXISTS model_metrics (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    model_id uuid NOT NULL REFERENCES ml_models(id),
    version_id uuid NOT NULL REFERENCES model_versions(id),
    metric_type text NOT NULL,
    metric_value double precision NOT NULL,
    metric_date timestamptz NOT NULL,
    dataset_info jsonb,
    calculation_method text,
    confidence_interval jsonb,
    metadata jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Create indexes
CREATE INDEX idx_models_name ON ml_models(name);
CREATE INDEX idx_models_type ON ml_models(model_type);
CREATE INDEX idx_models_target ON ml_models(target_variable);

CREATE INDEX idx_versions_model ON model_versions(model_id);
CREATE INDEX idx_versions_status ON model_versions(deployment_status);
CREATE INDEX idx_versions_deployed ON model_versions(deployed_at);

CREATE INDEX idx_training_model ON training_history(model_id);
CREATE INDEX idx_training_version ON training_history(version_id);
CREATE INDEX idx_training_status ON training_history(status);
CREATE INDEX idx_training_date ON training_history(start_time);

CREATE INDEX idx_predictions_model ON model_predictions(model_id);
CREATE INDEX idx_predictions_version ON model_predictions(version_id);
CREATE INDEX idx_predictions_compound ON model_predictions(compound_id);
CREATE INDEX idx_predictions_type ON model_predictions(prediction_type);
CREATE INDEX idx_predictions_date ON model_predictions(prediction_date);

CREATE INDEX idx_features_name ON feature_definitions(name);
CREATE INDEX idx_features_type ON feature_definitions(feature_type);

CREATE INDEX idx_feature_values_compound ON feature_values(compound_id);
CREATE INDEX idx_feature_values_feature ON feature_values(feature_id);
CREATE INDEX idx_feature_values_date ON feature_values(calculation_date);

CREATE INDEX idx_pipelines_name ON feature_pipelines(name);
CREATE INDEX idx_pipelines_input ON feature_pipelines USING gin(input_features);
CREATE INDEX idx_pipelines_output ON feature_pipelines USING gin(output_features);

CREATE INDEX idx_selection_model ON feature_selection(model_id);
CREATE INDEX idx_selection_method ON feature_selection(selection_method);
CREATE INDEX idx_selection_date ON feature_selection(selection_date);

CREATE INDEX idx_metrics_model ON model_metrics(model_id);
CREATE INDEX idx_metrics_version ON model_metrics(version_id);
CREATE INDEX idx_metrics_type ON model_metrics(metric_type);
CREATE INDEX idx_metrics_date ON model_metrics(metric_date);

-- Training datasets
CREATE TABLE IF NOT EXISTS training_datasets (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    dataset_name text NOT NULL UNIQUE,
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
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Model validation results
CREATE TABLE IF NOT EXISTS model_validation_results (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    model_id uuid NOT NULL REFERENCES ml_models(id) ON DELETE CASCADE,
    validation_type text NOT NULL,
    validation_date timestamptz NOT NULL,
    test_dataset_info jsonb,
    validation_metrics jsonb,
    test_set_performance jsonb,
    cross_validation_results jsonb,
    error_analysis jsonb,
    validation_plots jsonb,
    validation_notes text,
    recommendations text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Model deployments
CREATE TABLE IF NOT EXISTS model_deployments (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    model_id uuid NOT NULL REFERENCES ml_models(id) ON DELETE CASCADE,
    deployment_name text NOT NULL,
    deployment_environment text NOT NULL,
    deployment_date timestamptz NOT NULL,
    status text NOT NULL,
    version_tag text NOT NULL,
    configuration jsonb,
    performance_metrics jsonb,
    monitoring_config jsonb,
    rollback_info jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Model monitoring
CREATE TABLE IF NOT EXISTS model_monitoring (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    model_id uuid NOT NULL REFERENCES ml_models(id) ON DELETE CASCADE,
    monitoring_date timestamptz NOT NULL,
    metric_name text NOT NULL,
    metric_value double precision NOT NULL,
    threshold_value double precision,
    alert_status text,
    data_drift_metrics jsonb,
    performance_metrics jsonb,
    resource_usage jsonb,
    alert_history jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Create indexes for new tables
CREATE INDEX idx_training_datasets_name ON training_datasets(dataset_name);
CREATE INDEX idx_training_datasets_source ON training_datasets(source);
CREATE INDEX idx_training_datasets_version ON training_datasets(version);

CREATE INDEX idx_validation_results_model ON model_validation_results(model_id);
CREATE INDEX idx_validation_results_type ON model_validation_results(validation_type);
CREATE INDEX idx_validation_results_date ON model_validation_results(validation_date);

CREATE INDEX idx_model_deployments_model ON model_deployments(model_id);
CREATE INDEX idx_model_deployments_name ON model_deployments(deployment_name);
CREATE INDEX idx_model_deployments_env ON model_deployments(deployment_environment);
CREATE INDEX idx_model_deployments_status ON model_deployments(status);

CREATE INDEX idx_model_monitoring_model ON model_monitoring(model_id);
CREATE INDEX idx_model_monitoring_date ON model_monitoring(monitoring_date);
CREATE INDEX idx_model_monitoring_metric ON model_monitoring(metric_name);
CREATE INDEX idx_model_monitoring_status ON model_monitoring(alert_status);

-- Add triggers for timestamp updates
CREATE TRIGGER update_models_modtime
    BEFORE UPDATE ON ml_models
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_versions_modtime
    BEFORE UPDATE ON model_versions
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_training_modtime
    BEFORE UPDATE ON training_history
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_predictions_modtime
    BEFORE UPDATE ON model_predictions
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_features_modtime
    BEFORE UPDATE ON feature_definitions
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_feature_values_modtime
    BEFORE UPDATE ON feature_values
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_pipelines_modtime
    BEFORE UPDATE ON feature_pipelines
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_selection_modtime
    BEFORE UPDATE ON feature_selection
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_metrics_modtime
    BEFORE UPDATE ON model_metrics
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();
