-- Consolidated web interface schema
-- Combines web templates, components, settings, and API configurations

-- Web templates
CREATE TABLE IF NOT EXISTS web_templates (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    description text,
    template_type text NOT NULL,
    content text NOT NULL,
    parameters jsonb,
    styling jsonb,
    scripts jsonb,
    version text NOT NULL,
    is_active boolean DEFAULT true,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Web components
CREATE TABLE IF NOT EXISTS web_components (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    description text,
    component_type text NOT NULL,
    configuration jsonb NOT NULL,
    dependencies text[],
    styling jsonb,
    client_scripts jsonb,
    server_scripts jsonb,
    version text NOT NULL,
    is_active boolean DEFAULT true,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- User interface settings
CREATE TABLE IF NOT EXISTS ui_settings (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id text NOT NULL,
    theme text,
    layout_preferences jsonb,
    display_options jsonb,
    notification_settings jsonb,
    accessibility_settings jsonb,
    custom_views jsonb[],
    dashboard_config jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(user_id)
);

-- API endpoints
CREATE TABLE IF NOT EXISTS api_endpoints (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    path text NOT NULL UNIQUE,
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
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- API keys
CREATE TABLE IF NOT EXISTS api_keys (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    key_hash text NOT NULL UNIQUE,
    user_id text NOT NULL,
    name text NOT NULL,
    permissions text[],
    rate_limit integer,
    expires_at timestamptz,
    last_used_at timestamptz,
    is_active boolean DEFAULT true,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Rate limiting
CREATE TABLE IF NOT EXISTS rate_limits (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    api_key_id uuid REFERENCES api_keys(id) ON DELETE CASCADE,
    endpoint_id uuid REFERENCES api_endpoints(id) ON DELETE CASCADE,
    requests_count integer NOT NULL DEFAULT 0,
    window_start timestamptz NOT NULL,
    window_end timestamptz NOT NULL,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(api_key_id, endpoint_id, window_start)
);

-- Web hooks
CREATE TABLE IF NOT EXISTS web_hooks (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL,
    url text NOT NULL,
    event_types text[] NOT NULL,
    headers jsonb,
    is_active boolean DEFAULT true,
    secret_key text,
    retry_config jsonb,
    timeout interval,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Web hook deliveries
CREATE TABLE IF NOT EXISTS web_hook_deliveries (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    hook_id uuid NOT NULL REFERENCES web_hooks(id) ON DELETE CASCADE,
    event_type text NOT NULL,
    payload jsonb NOT NULL,
    response_status integer,
    response_body text,
    delivery_status text NOT NULL,
    attempt_count integer DEFAULT 0,
    next_retry_at timestamptz,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- External services
CREATE TABLE IF NOT EXISTS external_services (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    service_type text NOT NULL,
    base_url text NOT NULL,
    auth_config jsonb,
    rate_limit jsonb,
    timeout interval,
    retry_config jsonb,
    is_active boolean DEFAULT true,
    health_check_path text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Service status
CREATE TABLE IF NOT EXISTS service_status (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    service_id uuid NOT NULL REFERENCES external_services(id) ON DELETE CASCADE,
    status text NOT NULL,
    response_time interval,
    error_message text,
    check_timestamp timestamptz NOT NULL,
    metrics jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Analytics and monitoring
CREATE TABLE IF NOT EXISTS trend_analysis (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
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
    next_review_date timestamptz,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

CREATE TABLE IF NOT EXISTS alert_triggers (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
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
    is_active boolean NOT NULL DEFAULT true,
    last_triggered_at timestamptz,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

CREATE TABLE IF NOT EXISTS performance_metrics (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
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
    measured_at timestamptz NOT NULL DEFAULT now(),
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

CREATE TABLE IF NOT EXISTS usage_statistics (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
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
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Mobile support
CREATE TABLE IF NOT EXISTS user_preferences (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
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
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE (user_id, device_id)
);

CREATE TABLE IF NOT EXISTS device_settings (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    device_id text NOT NULL UNIQUE,
    device_type text NOT NULL,
    os_version text,
    app_version text,
    screen_size text,
    capabilities jsonb,
    settings jsonb,
    performance_profile jsonb,
    security_settings jsonb,
    network_preferences jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

CREATE TABLE IF NOT EXISTS offline_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    device_id text NOT NULL,
    data_type text NOT NULL,
    data_id uuid NOT NULL,
    content jsonb NOT NULL,
    version integer NOT NULL,
    priority integer DEFAULT 0,
    compression_type text,
    encryption_type text,
    validation_hash text,
    expires_at timestamptz,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE (device_id, data_type, data_id)
);

CREATE TABLE IF NOT EXISTS sync_status (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    device_id text NOT NULL,
    data_type text NOT NULL,
    last_sync_at timestamptz NOT NULL,
    sync_version integer NOT NULL,
    status text NOT NULL,
    conflict_resolution jsonb,
    retry_count integer DEFAULT 0,
    next_retry_at timestamptz,
    error_details text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE (device_id, data_type)
);

-- Create indexes
CREATE INDEX idx_templates_name ON web_templates(name);

-- Analytics indexes
CREATE INDEX idx_trend_analysis_type ON trend_analysis(analysis_type);
CREATE INDEX idx_trend_analysis_target ON trend_analysis(target_type, target_id);
CREATE INDEX idx_trend_analysis_period ON trend_analysis(time_period);
CREATE INDEX idx_trend_analysis_priority ON trend_analysis(priority);
CREATE INDEX idx_trend_analysis_status ON trend_analysis(status);

CREATE INDEX idx_alert_triggers_type ON alert_triggers(trigger_type);
CREATE INDEX idx_alert_triggers_target ON alert_triggers(target_type, target_id);
CREATE INDEX idx_alert_triggers_active ON alert_triggers(is_active);
CREATE INDEX idx_alert_triggers_severity ON alert_triggers(severity);
CREATE INDEX idx_alert_triggers_priority ON alert_triggers(priority);

CREATE INDEX idx_performance_metrics_type ON performance_metrics(metric_type);
CREATE INDEX idx_performance_metrics_component ON performance_metrics(component);
CREATE INDEX idx_performance_metrics_measured ON performance_metrics(measured_at);
CREATE INDEX idx_performance_metrics_status ON performance_metrics(status);
CREATE INDEX idx_performance_metrics_tags ON performance_metrics USING gin(tags);

CREATE INDEX idx_usage_statistics_feature ON usage_statistics(feature);
CREATE INDEX idx_usage_statistics_action ON usage_statistics(action);
CREATE INDEX idx_usage_statistics_device ON usage_statistics(device_type);
CREATE INDEX idx_usage_statistics_session ON usage_statistics(session_id);
CREATE INDEX idx_usage_statistics_user ON usage_statistics(user_id);
CREATE INDEX idx_usage_statistics_created ON usage_statistics(created_at);

CREATE INDEX idx_user_preferences_user ON user_preferences(user_id);
CREATE INDEX idx_user_preferences_device ON user_preferences(device_id);

CREATE INDEX idx_device_settings_type ON device_settings(device_type);
CREATE INDEX idx_device_settings_version ON device_settings(app_version);

CREATE INDEX idx_offline_data_device ON offline_data(device_id);
CREATE INDEX idx_offline_data_type ON offline_data(data_type);
CREATE INDEX idx_offline_data_priority ON offline_data(priority);
CREATE INDEX idx_offline_data_expires ON offline_data(expires_at);

CREATE INDEX idx_sync_status_device ON sync_status(device_id);
CREATE INDEX idx_sync_status_type ON sync_status(data_type);
CREATE INDEX idx_sync_status_sync ON sync_status(last_sync_at);
CREATE INDEX idx_sync_status_status ON sync_status(status);

CREATE INDEX idx_templates_type ON web_templates(template_type);
CREATE INDEX idx_templates_active ON web_templates(is_active);

CREATE INDEX idx_components_name ON web_components(name);
CREATE INDEX idx_components_type ON web_components(component_type);
CREATE INDEX idx_components_active ON web_components(is_active);

CREATE INDEX idx_settings_user ON ui_settings(user_id);

CREATE INDEX idx_endpoints_path ON api_endpoints(path);
CREATE INDEX idx_endpoints_method ON api_endpoints(method);
CREATE INDEX idx_endpoints_version ON api_endpoints(version);
CREATE INDEX idx_endpoints_deprecated ON api_endpoints(is_deprecated);

CREATE INDEX idx_api_keys_hash ON api_keys(key_hash);
CREATE INDEX idx_api_keys_user ON api_keys(user_id);
CREATE INDEX idx_api_keys_active ON api_keys(is_active);
CREATE INDEX idx_api_keys_expires ON api_keys(expires_at);

CREATE INDEX idx_rate_limits_key ON rate_limits(api_key_id);
CREATE INDEX idx_rate_limits_endpoint ON rate_limits(endpoint_id);
CREATE INDEX idx_rate_limits_window ON rate_limits(window_start, window_end);

CREATE INDEX idx_hooks_name ON web_hooks(name);
CREATE INDEX idx_hooks_active ON web_hooks(is_active);
CREATE INDEX idx_hooks_events ON web_hooks USING gin(event_types);

CREATE INDEX idx_hook_deliveries_hook ON web_hook_deliveries(hook_id);
CREATE INDEX idx_hook_deliveries_event ON web_hook_deliveries(event_type);
CREATE INDEX idx_hook_deliveries_status ON web_hook_deliveries(delivery_status);
CREATE INDEX idx_hook_deliveries_retry ON web_hook_deliveries(next_retry_at);

CREATE INDEX idx_services_name ON external_services(name);
CREATE INDEX idx_services_type ON external_services(service_type);
CREATE INDEX idx_services_active ON external_services(is_active);

CREATE INDEX idx_service_status_service ON service_status(service_id);
CREATE INDEX idx_service_status_timestamp ON service_status(check_timestamp);
CREATE INDEX idx_service_status_status ON service_status(status);

-- Add triggers for timestamp updates
CREATE TRIGGER update_templates_modtime
    BEFORE UPDATE ON web_templates
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_components_modtime
    BEFORE UPDATE ON web_components
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_settings_modtime
    BEFORE UPDATE ON ui_settings
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_endpoints_modtime
    BEFORE UPDATE ON api_endpoints
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_api_keys_modtime
    BEFORE UPDATE ON api_keys
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_rate_limits_modtime
    BEFORE UPDATE ON rate_limits
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_hooks_modtime
    BEFORE UPDATE ON web_hooks
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_hook_deliveries_modtime
    BEFORE UPDATE ON web_hook_deliveries
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_services_modtime
    BEFORE UPDATE ON external_services
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_service_status_modtime
    BEFORE UPDATE ON service_status
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Add triggers for analytics tables
CREATE TRIGGER update_trend_analysis_modtime
    BEFORE UPDATE ON trend_analysis
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_alert_triggers_modtime
    BEFORE UPDATE ON alert_triggers
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_performance_metrics_modtime
    BEFORE UPDATE ON performance_metrics
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_usage_statistics_modtime
    BEFORE UPDATE ON usage_statistics
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_user_preferences_modtime
    BEFORE UPDATE ON user_preferences
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_device_settings_modtime
    BEFORE UPDATE ON device_settings
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_offline_data_modtime
    BEFORE UPDATE ON offline_data
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_sync_status_modtime
    BEFORE UPDATE ON sync_status
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();
