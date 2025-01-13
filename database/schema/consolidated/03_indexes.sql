-- Create indexes for all tables
CREATE INDEX IF NOT EXISTS idx_compounds_name ON compounds(name);
CREATE INDEX IF NOT EXISTS idx_compounds_inchi_key ON compounds(inchi_key);
CREATE INDEX IF NOT EXISTS idx_compounds_smiles ON compounds(smiles);
CREATE INDEX IF NOT EXISTS idx_compounds_cas_number ON compounds(cas_number);
CREATE INDEX IF NOT EXISTS idx_compounds_pubchem_cid ON compounds(pubchem_cid);
CREATE INDEX IF NOT EXISTS idx_compounds_chembl_id ON compounds(chembl_id);
CREATE INDEX IF NOT EXISTS idx_compounds_drugbank_id ON compounds(drugbank_id);
CREATE INDEX IF NOT EXISTS idx_compounds_unii ON compounds(unii);
CREATE INDEX IF NOT EXISTS idx_compounds_kegg_id ON compounds(kegg_id);
CREATE INDEX IF NOT EXISTS idx_compounds_chemspider_id ON compounds(chemspider_id);
CREATE INDEX IF NOT EXISTS idx_compounds_zinc_id ON compounds(zinc_id);
CREATE INDEX IF NOT EXISTS idx_compounds_chebi_id ON compounds(chebi_id);
CREATE INDEX IF NOT EXISTS idx_compounds_common_names ON compounds USING gin(common_names);

-- Indexes for enhanced properties
CREATE INDEX IF NOT EXISTS idx_compounds_drug_class ON compounds USING gin (drug_class);
CREATE INDEX IF NOT EXISTS idx_compounds_mechanism ON compounds USING gin (mechanism_categories);
CREATE INDEX IF NOT EXISTS idx_compounds_therapeutic ON compounds USING gin (therapeutic_categories);
CREATE INDEX IF NOT EXISTS idx_compounds_effects ON compounds USING gin (pharmacological_effects);
CREATE INDEX IF NOT EXISTS idx_compounds_routes ON compounds USING gin (administration_routes);

-- Indexes for descriptor tables
CREATE INDEX IF NOT EXISTS idx_descriptors_2d_compound_id ON descriptors_2d(compound_id);
CREATE INDEX IF NOT EXISTS idx_descriptors_3d_compound_id ON descriptors_3d(compound_id);
CREATE INDEX IF NOT EXISTS idx_molecular_fingerprints_compound_id ON molecular_fingerprints(compound_id);
CREATE INDEX IF NOT EXISTS idx_pharmacophore_features_compound_id ON pharmacophore_features(compound_id);

-- Indexes for quantum tables
CREATE INDEX IF NOT EXISTS idx_electronic_structure_compound_id ON electronic_structure(compound_id);
CREATE INDEX IF NOT EXISTS idx_quantum_critical_params_compound_id ON quantum_critical_params(compound_id);
CREATE INDEX IF NOT EXISTS idx_quantum_dynamics_compound_id ON quantum_dynamics(compound_id);
CREATE INDEX IF NOT EXISTS idx_phase_transitions_compound_id ON phase_transitions(compound_id);
CREATE INDEX IF NOT EXISTS idx_scaling_analysis_compound_id ON scaling_analysis(compound_id);
CREATE INDEX IF NOT EXISTS idx_quantum_observables_compound_id ON quantum_observables(compound_id);

-- Indexes for binding tables
CREATE INDEX IF NOT EXISTS idx_binding_data_compound ON binding_data(compound_id);
CREATE INDEX IF NOT EXISTS idx_binding_data_receptor ON binding_data(receptor_family_id);
CREATE INDEX IF NOT EXISTS idx_binding_data_assay ON binding_data(assay_type_id);
CREATE INDEX IF NOT EXISTS idx_binding_data_value ON binding_data(value);
CREATE INDEX IF NOT EXISTS idx_binding_data_source ON binding_data(data_source);
CREATE INDEX IF NOT EXISTS idx_binding_quality_data ON binding_data_quality(binding_data_id);
CREATE INDEX IF NOT EXISTS idx_binding_sar_compound ON binding_sar(compound_id);
CREATE INDEX IF NOT EXISTS idx_binding_sar_receptor ON binding_sar(receptor_family_id);
CREATE INDEX IF NOT EXISTS idx_binding_kinetics_binding ON binding_kinetics(binding_data_id);
CREATE INDEX IF NOT EXISTS idx_binding_site_mapping_binding ON binding_site_mapping(binding_data_id);
CREATE INDEX IF NOT EXISTS idx_binding_site_mapping_site ON binding_site_mapping(site_name);
CREATE INDEX IF NOT EXISTS idx_binding_assay_protocols_type ON binding_assay_protocols(assay_type);

-- Create GiST index for multidimensional data
CREATE INDEX IF NOT EXISTS idx_density_grid_points ON electronic_structure USING gin (density_grid_points);
CREATE INDEX IF NOT EXISTS idx_molecular_fingerprints_ecfp4 ON molecular_fingerprints USING gist (ecfp4_bits);
CREATE INDEX IF NOT EXISTS idx_molecular_fingerprints_maccs ON molecular_fingerprints USING gist (maccs_keys);

-- Add triggers for updated_at timestamps
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
            'compounds',
            'descriptors_2d',
            'descriptors_3d',
            'molecular_fingerprints',
            'pharmacophore_features',
            'electronic_structure',
            'quantum_critical_params',
            'quantum_dynamics',
            'phase_transitions',
            'scaling_analysis',
            'quantum_observables',
            'binding_assay_types',
            'binding_data',
            'binding_data_quality',
            'binding_sar',
            'binding_kinetics',
            'binding_site_mapping',
            'binding_assay_protocols'
        )
    LOOP
        EXECUTE format('
            CREATE TRIGGER update_%I_modtime 
            BEFORE UPDATE ON %I 
            FOR EACH ROW
            EXECUTE FUNCTION update_updated_at_column();
            
            CREATE TRIGGER audit_%I_trigger
            AFTER INSERT OR UPDATE OR DELETE ON %I
            FOR EACH ROW
            EXECUTE FUNCTION audit_trigger_func();',
            t, t, t, t);
    END LOOP;
END $$;
