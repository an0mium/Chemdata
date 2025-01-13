-- Supplier and product information tables
CREATE TABLE IF NOT EXISTS chemical_suppliers (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL UNIQUE,
    website text,
    api_base_url text,
    catalog_url text,
    contact_info jsonb,
    quality_certifications text[],
    shipping_regions text[],
    special_requirements text[],
    notes text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Product catalog entries
CREATE TABLE IF NOT EXISTS supplier_products (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    supplier_id uuid NOT NULL REFERENCES chemical_suppliers(id),
    catalog_number text NOT NULL,
    cas_number text,
    product_name text NOT NULL,
    grade text, -- e.g., 'research', 'analytical', 'HPLC'
    purity numeric,
    purity_method text,
    form text, -- e.g., 'powder', 'solution'
    package_sizes jsonb, -- Array of available sizes and units
    prices jsonb, -- Price information by package size
    currency text,
    availability_status text,
    lead_time interval,
    shipping_conditions text[],
    regulatory_info jsonb,
    technical_info_url text,
    last_checked timestamptz,
    last_price_update timestamptz,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now(),
    UNIQUE(supplier_id, catalog_number)
);

-- Quality control data
CREATE TABLE IF NOT EXISTS supplier_qc_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    product_id uuid NOT NULL REFERENCES supplier_products(id),
    batch_number text,
    manufacture_date date,
    expiry_date date,
    analysis_date date,
    test_method text,
    test_parameter text,
    specification text,
    result text,
    units text,
    analyst text,
    equipment_used text[],
    notes text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Supplier certifications and compliance
CREATE TABLE IF NOT EXISTS supplier_certifications (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    supplier_id uuid NOT NULL REFERENCES chemical_suppliers(id),
    certification_type text NOT NULL,
    certification_number text,
    issuing_body text,
    issue_date date,
    expiry_date date,
    scope text[],
    status text,
    documentation_url text,
    audit_history jsonb,
    notes text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Create indexes for supplier tables
CREATE INDEX IF NOT EXISTS idx_supplier_products_compound ON supplier_products(compound_id);
CREATE INDEX IF NOT EXISTS idx_supplier_products_supplier ON supplier_products(supplier_id);
CREATE INDEX IF NOT EXISTS idx_supplier_products_catalog ON supplier_products(catalog_number);
CREATE INDEX IF NOT EXISTS idx_supplier_products_cas ON supplier_products(cas_number);
CREATE INDEX IF NOT EXISTS idx_supplier_qc_data_product ON supplier_qc_data(product_id);
CREATE INDEX IF NOT EXISTS idx_supplier_certifications_supplier ON supplier_certifications(supplier_id);

-- Add triggers for supplier tables
CREATE TRIGGER update_chemical_suppliers_modtime
    BEFORE UPDATE ON chemical_suppliers
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_supplier_products_modtime
    BEFORE UPDATE ON supplier_products
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_supplier_qc_data_modtime
    BEFORE UPDATE ON supplier_qc_data
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_supplier_certifications_modtime
    BEFORE UPDATE ON supplier_certifications
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Add audit triggers for supplier tables
CREATE TRIGGER audit_chemical_suppliers_trigger
    AFTER INSERT OR UPDATE OR DELETE ON chemical_suppliers
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_supplier_products_trigger
    AFTER INSERT OR UPDATE OR DELETE ON supplier_products
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_supplier_qc_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON supplier_qc_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_supplier_certifications_trigger
    AFTER INSERT OR UPDATE OR DELETE ON supplier_certifications
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Add table comments for supplier tables
COMMENT ON TABLE chemical_suppliers IS 'Information about chemical suppliers and their capabilities';
COMMENT ON TABLE supplier_products IS 'Product catalog entries from chemical suppliers';
COMMENT ON TABLE supplier_qc_data IS 'Quality control data for supplier products';
COMMENT ON TABLE supplier_certifications IS 'Supplier certifications and compliance information';

-- Create view for product availability
CREATE OR REPLACE VIEW compound_availability AS
SELECT 
    c.id as compound_id,
    c.name as compound_name,
    sp.supplier_id,
    cs.name as supplier_name,
    sp.catalog_number,
    sp.product_name,
    sp.grade,
    sp.purity,
    sp.form,
    sp.package_sizes,
    sp.prices,
    sp.currency,
    sp.availability_status,
    sp.lead_time,
    sp.shipping_conditions,
    sp.last_checked,
    sp.last_price_update
FROM compounds c
JOIN supplier_products sp ON c.id = sp.compound_id
JOIN chemical_suppliers cs ON sp.supplier_id = cs.id
WHERE sp.availability_status = 'in_stock'
ORDER BY c.name, cs.name;

-- Create view for supplier quality metrics
CREATE OR REPLACE VIEW supplier_quality_metrics AS
SELECT 
    cs.id as supplier_id,
    cs.name as supplier_name,
    COUNT(DISTINCT sp.id) as product_count,
    COUNT(DISTINCT sqd.id) as qc_records_count,
    COUNT(DISTINCT sc.id) as certification_count,
    AVG(sp.purity) as avg_product_purity,
    COUNT(DISTINCT CASE WHEN sp.availability_status = 'in_stock' THEN sp.id END) as in_stock_count,
    array_agg(DISTINCT sc.certification_type) as certifications
FROM chemical_suppliers cs
LEFT JOIN supplier_products sp ON cs.id = sp.supplier_id
LEFT JOIN supplier_qc_data sqd ON sp.id = sqd.product_id
LEFT JOIN supplier_certifications sc ON cs.id = sc.supplier_id
GROUP BY cs.id, cs.name;

-- Insert reference data for major suppliers
INSERT INTO chemical_suppliers (name, website, api_base_url, catalog_url) VALUES
('Millipore Sigma', 'https://www.sigmaaldrich.com', 'https://api.sigmaaldrich.com', 'https://www.sigmaaldrich.com/catalog'),
('Cayman Chemical', 'https://www.caymanchem.com', 'https://api.caymanchem.com', 'https://www.caymanchem.com/products'),
('Tocris', 'https://www.tocris.com', 'https://api.tocris.com', 'https://www.tocris.com/products'),
('Toronto Research Chemicals', 'https://www.trc-canada.com', 'https://api.trc-canada.com', 'https://www.trc-canada.com/products-listing');
