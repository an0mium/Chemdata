-- Core compound tables and properties

-- Base compound table with enhanced fields
CREATE TABLE IF NOT EXISTS compounds (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT NOT NULL,
    smiles TEXT NOT NULL,
    inchi TEXT,
    inchi_key TEXT UNIQUE,
    cas_number TEXT UNIQUE,
    pubchem_cid TEXT,
    chembl_id TEXT,
    drugbank_id TEXT,
    unii TEXT,
    kegg_id TEXT,
    chemspider_id TEXT,
    zinc_id TEXT,
    chebi_id TEXT,
    iupac_name TEXT,
    preferred_iupac_name TEXT,
    common_names TEXT[],
    einecs_number TEXT,
    rtecs_number TEXT,
    hsdb_number TEXT,
    ccdc_number TEXT,
    reaxys_id TEXT,
    lipidmaps_id TEXT,
    nsc_number TEXT,
    qcarchive_id TEXT,
    nomad_id TEXT,
    materials_project_id TEXT,
    basis_set_id TEXT,
    method_id TEXT,
    calculation_id TEXT,
    molecular_weight FLOAT,
    molecular_formula TEXT,
    -- Enhanced properties
    stereochemistry jsonb,
    crystal_structure jsonb,
    solubility_data jsonb,
    pka_values jsonb,
    partition_coefficients jsonb,
    surface_properties jsonb,
    conformational_analysis jsonb,
    chirality_info jsonb,
    isomer_details jsonb,
    melting_point double precision,
    boiling_point double precision,
    density double precision,
    refractive_index double precision,
    optical_rotation double precision,
    drug_class text[],
    mechanism_categories text[],
    therapeutic_categories text[],
    pharmacological_effects text[],
    administration_routes text[],
    bioavailability_data jsonb,
    metabolism_data jsonb,
    distribution_data jsonb,
    legal_status jsonb,
    scheduling_info jsonb,
    approval_status jsonb,
    clinical_trial_status jsonb,
    patent_status jsonb,
    registration_numbers jsonb,
    control_status jsonb,
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT valid_smiles CHECK (validate_smiles(smiles)),
    CONSTRAINT valid_cas_number CHECK (cas_number ~ '^[0-9]{1,7}-[0-9]{2}-[0-9]$'),
    CONSTRAINT valid_pubchem_cid CHECK (pubchem_cid ~ '^[0-9]+$'),
    CONSTRAINT valid_chembl_id CHECK (chembl_id ~ '^CHEMBL[0-9]+$'),
    CONSTRAINT valid_drugbank_id CHECK (drugbank_id ~ '^DB[0-9]{5}$'),
    CONSTRAINT valid_unii CHECK (unii ~ '^[A-Z0-9]{10}$'),
    CONSTRAINT valid_kegg_id CHECK (kegg_id ~ '^C[0-9]{5}$|^D[0-9]{5}$'),
    CONSTRAINT valid_chemspider_id CHECK (chemspider_id ~ '^[0-9]+$'),
    CONSTRAINT valid_zinc_id CHECK (zinc_id ~ '^ZINC[0-9]{12}$'),
    CONSTRAINT valid_chebi_id CHECK (chebi_id ~ '^CHEBI:[0-9]+$'),
    CONSTRAINT valid_einecs CHECK (einecs_number ~ '^[0-9]{3}-[0-9]{3}-[0-9]$'),
    CONSTRAINT valid_rtecs CHECK (rtecs_number ~ '^[A-Z]{2}[0-9]{4,7}$'),
    CONSTRAINT valid_hsdb CHECK (hsdb_number ~ '^[0-9]+$'),
    CONSTRAINT valid_ccdc CHECK (ccdc_number ~ '^[0-9]{6}$'),
    CONSTRAINT valid_reaxys CHECK (reaxys_id ~ '^RX[0-9]+$'),
    CONSTRAINT valid_lipidmaps CHECK (lipidmaps_id ~ '^LM[A-Z]{2}[0-9]{8}$'),
    CONSTRAINT valid_nsc CHECK (nsc_number ~ '^NSC[0-9]+$'),
    CONSTRAINT valid_qcarchive CHECK (qcarchive_id ~ '^QCA[0-9]+$'),
    CONSTRAINT valid_nomad CHECK (nomad_id ~ '^NOMAD[0-9]+$'),
    CONSTRAINT valid_materials_project CHECK (materials_project_id ~ '^mp-[0-9]+$'),
    CONSTRAINT valid_basis_set CHECK (basis_set_id ~ '^[A-Z0-9\-]+$'),
    CONSTRAINT valid_method CHECK (method_id ~ '^[A-Z0-9\-]+$'),
    CONSTRAINT valid_calculation CHECK (calculation_id ~ '^CALC[0-9]+$')
);

-- Add comprehensive comments for all fields
COMMENT ON TABLE compounds IS 'Core compound information and properties';
COMMENT ON COLUMN compounds.id IS 'Unique identifier for each compound';
COMMENT ON COLUMN compounds.name IS 'Common or systematic name of the compound';
COMMENT ON COLUMN compounds.smiles IS 'SMILES notation representing chemical structure';
COMMENT ON COLUMN compounds.inchi IS 'InChI identifier for the compound';
COMMENT ON COLUMN compounds.inchi_key IS 'InChIKey for efficient compound lookup';
COMMENT ON COLUMN compounds.cas_number IS 'CAS Registry Number';
COMMENT ON COLUMN compounds.pubchem_cid IS 'PubChem Compound ID';
COMMENT ON COLUMN compounds.chembl_id IS 'ChEMBL database identifier';
COMMENT ON COLUMN compounds.drugbank_id IS 'DrugBank database identifier';
COMMENT ON COLUMN compounds.unii IS 'FDA Unique Ingredient Identifier';
COMMENT ON COLUMN compounds.kegg_id IS 'KEGG database identifier';
COMMENT ON COLUMN compounds.chemspider_id IS 'ChemSpider database identifier';
COMMENT ON COLUMN compounds.zinc_id IS 'ZINC database identifier';
COMMENT ON COLUMN compounds.chebi_id IS 'ChEBI identifier';
COMMENT ON COLUMN compounds.iupac_name IS 'IUPAC systematic name';
COMMENT ON COLUMN compounds.preferred_iupac_name IS 'Preferred IUPAC Name (PIN)';
COMMENT ON COLUMN compounds.common_names IS 'Array of common/trivial names';
COMMENT ON COLUMN compounds.molecular_weight IS 'Molecular weight in g/mol';
COMMENT ON COLUMN compounds.molecular_formula IS 'Molecular formula';

-- Comments for enhanced properties
COMMENT ON COLUMN compounds.stereochemistry IS 'Stereochemical configuration details';
COMMENT ON COLUMN compounds.crystal_structure IS 'Crystal structure parameters';
COMMENT ON COLUMN compounds.solubility_data IS 'Solubility in various solvents';
COMMENT ON COLUMN compounds.pka_values IS 'Acid dissociation constants';
COMMENT ON COLUMN compounds.partition_coefficients IS 'Various partition coefficients';
COMMENT ON COLUMN compounds.surface_properties IS 'Surface tension, etc.';
COMMENT ON COLUMN compounds.conformational_analysis IS 'Conformational states';
COMMENT ON COLUMN compounds.chirality_info IS 'Chirality details';
COMMENT ON COLUMN compounds.isomer_details IS 'Isomer information';
COMMENT ON COLUMN compounds.drug_class IS 'Therapeutic/pharmacological classes';
COMMENT ON COLUMN compounds.mechanism_categories IS 'Mechanism of action categories';
COMMENT ON COLUMN compounds.therapeutic_categories IS 'Therapeutic use categories';
COMMENT ON COLUMN compounds.pharmacological_effects IS 'Known pharmacological effects';
COMMENT ON COLUMN compounds.administration_routes IS 'Routes of administration';
COMMENT ON COLUMN compounds.bioavailability_data IS 'Bioavailability information';
COMMENT ON COLUMN compounds.metabolism_data IS 'Metabolic pathway information';
COMMENT ON COLUMN compounds.distribution_data IS 'Distribution information';
COMMENT ON COLUMN compounds.legal_status IS 'Legal status by jurisdiction';
COMMENT ON COLUMN compounds.scheduling_info IS 'Drug scheduling information';
COMMENT ON COLUMN compounds.approval_status IS 'Approval status by region';
COMMENT ON COLUMN compounds.clinical_trial_status IS 'Clinical trial information';
COMMENT ON COLUMN compounds.patent_status IS 'Patent status information';
COMMENT ON COLUMN compounds.registration_numbers IS 'Various registration numbers';
COMMENT ON COLUMN compounds.control_status IS 'Control status by jurisdiction';

-- 2D Descriptors (Extended)
CREATE TABLE IF NOT EXISTS descriptors_2d (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    -- Topological descriptors
    wiener_index FLOAT,
    balaban_j_index FLOAT,
    bertz_ct FLOAT,
    schultz_molecular_topological_index FLOAT,
    -- Connectivity indices
    chi_0_index FLOAT,
    chi_1_index FLOAT,
    chi_2_index FLOAT,
    chi_3_index FLOAT,
    -- Electrotopological indices
    sum_estate_indices FLOAT,
    mean_estate_indices FLOAT,
    -- Shape descriptors
    kappa_1 FLOAT,
    kappa_2 FLOAT,
    kappa_3 FLOAT,
    -- Additional physicochemical
    molar_refractivity FLOAT,
    van_der_waals_volume FLOAT,
    polarizability FLOAT,
    formal_charge INTEGER,
    -- Ring descriptors
    ring_count INTEGER,
    aromatic_ring_count INTEGER,
    aliphatic_ring_count INTEGER,
    ring_fusion_degree INTEGER,
    -- Fragment-based
    rotatable_bond_count INTEGER,
    rigid_bond_count INTEGER,
    chain_atom_count INTEGER,
    chain_bond_count INTEGER,
    -- Atom counts by type
    carbon_count INTEGER,
    nitrogen_count INTEGER,
    oxygen_count INTEGER,
    sulfur_count INTEGER,
    phosphorus_count INTEGER,
    halogen_count INTEGER,
    created_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    -- Validation
    CONSTRAINT valid_counts CHECK (
        ring_count >= 0 AND
        aromatic_ring_count >= 0 AND
        aliphatic_ring_count >= 0 AND
        rotatable_bond_count >= 0 AND
        chain_atom_count >= 0
    )
);

-- 3D Descriptors
CREATE TABLE IF NOT EXISTS descriptors_3d (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    conformer_id INTEGER,  -- For multiple conformers
    -- Geometric descriptors
    radius_of_gyration FLOAT,
    molecular_volume FLOAT,
    molecular_surface_area FLOAT,
    solvent_accessible_surface_area FLOAT,
    polar_surface_area_3d FLOAT,
    -- Shape descriptors
    spherosity FLOAT,
    asphericity FLOAT,
    eccentricity FLOAT,
    inertial_shape_factor FLOAT,
    -- Moment descriptors
    principal_moment_1 FLOAT,
    principal_moment_2 FLOAT,
    principal_moment_3 FLOAT,
    -- Distance matrices
    gravitational_index FLOAT,
    radius_of_distribution FLOAT,
    -- Surface properties
    molecular_surface_potential FLOAT,
    average_surface_charge FLOAT,
    -- Conformational energies
    total_energy FLOAT,
    strain_energy FLOAT,
    created_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    -- Validation
    CONSTRAINT valid_3d_ranges CHECK (
        molecular_volume > 0 AND
        molecular_surface_area > 0 AND
        solvent_accessible_surface_area > 0
    )
);

-- Molecular Fingerprints
CREATE TABLE IF NOT EXISTS molecular_fingerprints (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    -- Structural fingerprints
    maccs_keys BIT VARYING(166),  -- MACCS keys
    pubchem_bits BIT VARYING(881),  -- PubChem fingerprint
    -- Circular fingerprints (ECFP)
    ecfp4_bits BIT VARYING(1024),
    ecfp6_bits BIT VARYING(1024),
    -- Path-based fingerprints
    daylight_bits BIT VARYING(1024),
    -- Pharmacophore fingerprints
    pharma_bits BIT VARYING(512),
    -- Additional fingerprint types
    atom_pairs BIT VARYING(1024),
    torsion_bits BIT VARYING(1024),
    morgan_bits BIT VARYING(1024),
    created_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP
);

-- Pharmacophore Features
CREATE TABLE IF NOT EXISTS pharmacophore_features (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    -- Feature counts
    h_bond_donors_count INTEGER,
    h_bond_acceptors_count INTEGER,
    pos_charge_groups_count INTEGER,
    neg_charge_groups_count INTEGER,
    aromatic_rings_count INTEGER,
    hydrophobic_groups_count INTEGER,
    -- Feature positions (3D)
    donor_positions JSONB,  -- Array of 3D coordinates
    acceptor_positions JSONB,
    charge_positions JSONB,
    aromatic_positions JSONB,
    hydrophobic_positions JSONB,
    -- Feature properties
    donor_strengths FLOAT[],
    acceptor_strengths FLOAT[],
    charge_strengths FLOAT[],
    -- Feature type and detection
    feature_type text NOT NULL,
    coordinates jsonb,
    strength double precision,
    interaction_radius double precision,
    optional boolean,
    detection_method text,
    confidence_score double precision,
    created_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    -- Validation
    CONSTRAINT valid_feature_counts CHECK (
        h_bond_donors_count >= 0 AND
        h_bond_acceptors_count >= 0 AND
        pos_charge_groups_count >= 0 AND
        neg_charge_groups_count >= 0 AND
        aromatic_rings_count >= 0 AND
        hydrophobic_groups_count >= 0
    )
);

-- Electronic Structure Data
CREATE TABLE IF NOT EXISTS electronic_structure (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    -- Core electronic properties
    total_electronic_energy FLOAT,
    electron_correlation_energy FLOAT,
    exchange_energy FLOAT,
    kinetic_energy FLOAT,
    -- Orbital data
    homo_lumo_gap FLOAT,
    orbital_energies FLOAT[],
    orbital_occupancies INTEGER[],
    -- Electron density
    density_matrix JSONB,
    density_grid_points JSONB,
    density_values FLOAT[],
    -- Wavefunction data
    wavefunction JSONB,
    state_type VARCHAR(50),
    -- Band structure
    band_structure JSONB,
    -- Calculation metadata
    method TEXT,
    basis_set TEXT,
    convergence_criteria JSONB,
    created_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    metadata JSONB,
    CONSTRAINT valid_energies CHECK (
        total_electronic_energy < 0 AND
        homo_lumo_gap >= 0
    )
);

-- Quantum Critical Parameters
CREATE TABLE IF NOT EXISTS quantum_critical_params (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    -- Critical point parameters
    critical_temperature FLOAT,
    critical_pressure FLOAT,
    critical_field FLOAT,
    -- Order parameters
    primary_order_parameter TEXT,
    order_parameter_values JSONB,
    -- Critical exponents
    alpha FLOAT,  -- Specific heat exponent
    beta FLOAT,   -- Order parameter exponent
    gamma FLOAT,  -- Susceptibility exponent
    delta FLOAT,  -- Critical isotherm exponent
    nu FLOAT,     -- Correlation length exponent
    eta FLOAT,    -- Anomalous dimension
    -- Correlation functions
    correlation_length FLOAT,
    correlation_function JSONB,
    dynamic_exponent_z FLOAT,
    -- Phase diagram
    phase_boundaries JSONB,
    multicriticality_type TEXT,
    -- Quantum properties
    coherence_length FLOAT,
    entanglement_entropy FLOAT,
    quantum_fluctuations JSONB,
    created_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    metadata JSONB,
    CONSTRAINT valid_critical_params CHECK (
        critical_temperature >= 0 AND
        critical_pressure >= 0 AND
        correlation_length > 0
    )
);

-- Quantum Dynamics
CREATE TABLE IF NOT EXISTS quantum_dynamics (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    -- Time evolution
    time_points FLOAT[],
    wavefunction_evolution JSONB,
    density_matrix_evolution JSONB,
    -- Dynamic properties
    coherence_times FLOAT[],
    relaxation_rates FLOAT[],
    dephasing_rates FLOAT[],
    -- Transport properties
    conductivity_tensor FLOAT[],
    hall_conductance FLOAT,
    thermal_conductivity FLOAT,
    -- Spectral properties
    spectral_function JSONB,
    optical_conductivity JSONB,
    -- Quantum correlations
    entanglement_spectrum FLOAT[],
    mutual_information FLOAT,
    -- Environmental coupling
    dissipation_kernel JSONB,
    noise_spectrum JSONB,
    metadata JSONB,
    -- Note: Array validation will be handled by triggers instead of CHECK constraints
    -- since PostgreSQL doesn't support complex array validation in CHECK constraints
    created_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP
);

-- Phase Transitions
CREATE TABLE IF NOT EXISTS phase_transitions (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    transition_type VARCHAR(50) NOT NULL,
    critical_temperature FLOAT,
    critical_pressure FLOAT,
    order_parameter JSONB,
    correlation_length FLOAT,
    transition_order INTEGER,
    hysteresis_data JSONB,
    fluctuation_data JSONB,
    created_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    metadata JSONB,
    CONSTRAINT valid_transition_params CHECK (
        critical_temperature >= 0 AND
        critical_pressure >= 0 AND
        correlation_length > 0 AND
        transition_order > 0
    )
);

-- Scaling Analysis
CREATE TABLE IF NOT EXISTS scaling_analysis (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    phase_transition_id uuid REFERENCES phase_transitions(id),
    -- Scaling functions
    scaling_function_type TEXT,
    scaling_variables JSONB,
    scaling_dimensions FLOAT[],
    -- RG analysis
    rg_flow_equations JSONB,
    fixed_points JSONB,
    relevant_operators JSONB,
    -- Universal properties
    universality_class TEXT,
    central_charge FLOAT,
    operator_spectrum JSONB,
    -- Finite-size scaling
    size_scaling_exponents FLOAT[],
    correction_exponents FLOAT[],
    -- Crossover behavior
    crossover_scales JSONB,
    crossover_functions JSONB,
    -- Metadata
    analysis_method TEXT,
    confidence_metrics JSONB,
    created_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    metadata JSONB
);

-- Quantum Observables
CREATE TABLE IF NOT EXISTS quantum_observables (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    observable_type VARCHAR(50) NOT NULL,
    value FLOAT,
    uncertainty FLOAT,
    measurement_basis TEXT,
    operator_type TEXT,
    expectation_value FLOAT,
    variance FLOAT,
    created_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at timestamptz NOT NULL DEFAULT CURRENT_TIMESTAMP,
    metadata JSONB,
    CONSTRAINT valid_uncertainty CHECK (
        uncertainty >= 0
    )
);

-- Binding assay types
CREATE TABLE IF NOT EXISTS binding_assay_types (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT NOT NULL,
    description TEXT,
    method_type TEXT,
    detection_type TEXT,
    typical_unit TEXT,
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP
);

-- Binding data measurements with enhanced fields
CREATE TABLE IF NOT EXISTS binding_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    receptor_family_id uuid NOT NULL REFERENCES receptor_families(id),
    assay_type_id uuid NOT NULL REFERENCES binding_assay_types(id),
    value FLOAT NOT NULL,
    unit TEXT NOT NULL,
    confidence_score FLOAT CHECK (confidence_score BETWEEN 0 AND 1),
    data_source TEXT,
    publication_doi TEXT,
    experimental_conditions JSONB,
    measurement_date DATE,
    -- Enhanced binding fields
    binding_site text,
    binding_mode text,
    binding_kinetics jsonb,
    activity_type text,
    activity_value double precision,
    activity_unit text,
    efficacy double precision,
    potency double precision,
    assay_description text,
    assay_organism text,
    assay_type text,
    assay_conditions jsonb,
    assay_cell_line text,
    assay_tissue_type text,
    sar_data jsonb,
    pharmacophore_features jsonb,
    binding_pocket_residues text[],
    experimental_method text,
    validation_method text,
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT valid_value CHECK (value > 0)
);

-- Binding data quality metrics
CREATE TABLE IF NOT EXISTS binding_data_quality (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    binding_data_id uuid NOT NULL REFERENCES binding_data(id) ON DELETE CASCADE,
    replicate_count INTEGER,
    standard_deviation FLOAT,
    confidence_interval FLOAT,
    quality_score FLOAT CHECK (quality_score BETWEEN 0 AND 1),
    validation_notes TEXT,
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP
);

-- Structure-activity relationships
CREATE TABLE IF NOT EXISTS binding_sar (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    receptor_family_id uuid NOT NULL REFERENCES receptor_families(id),
    structural_feature TEXT,
    effect_type TEXT,
    effect_magnitude FLOAT,
    confidence_score FLOAT CHECK (confidence_score BETWEEN 0 AND 1),
    evidence_type TEXT,
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP
);

-- Binding kinetics
CREATE TABLE IF NOT EXISTS binding_kinetics (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    binding_data_id uuid NOT NULL REFERENCES binding_data(id) ON DELETE CASCADE,
    kon_value double precision,
    kon_unit text,
    koff_value double precision,
    koff_unit text,
    residence_time double precision,
    residence_time_unit text,
    temperature double precision,
    ph double precision,
    ionic_strength text,
    buffer_conditions jsonb,
    method_details text,
    equipment_used text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Binding site mapping
CREATE TABLE IF NOT EXISTS binding_site_mapping (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    binding_data_id uuid NOT NULL REFERENCES binding_data(id) ON DELETE CASCADE,
    site_name text NOT NULL,
    residues text[],
    interaction_types text[],
    binding_pocket_volume double precision,
    surface_accessibility double precision,
    conservation_score double precision,
    mutation_effects jsonb,
    structural_features jsonb,
    modeling_method text,
    confidence_score double precision,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Binding assay protocols
CREATE TABLE IF NOT EXISTS binding_assay_protocols (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name text NOT NULL,
    description text,
    assay_type text NOT NULL,
    protocol_steps jsonb,
    reagents jsonb,
    equipment text[],
    controls jsonb,
    validation_criteria jsonb,
    limitations text[],
    reference_list text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Add comments for binding tables
COMMENT ON TABLE binding_assay_types IS 'Types of binding assays and their characteristics';
COMMENT ON TABLE binding_data IS 'Core table for receptor binding measurements';
COMMENT ON TABLE binding_data_quality IS 'Quality metrics for binding measurements';
COMMENT ON TABLE binding_sar IS 'Structure-activity relationships for binding data';
COMMENT ON TABLE binding_kinetics IS 'Detailed binding kinetics measurements';
COMMENT ON TABLE binding_site_mapping IS 'Binding site characterization and mapping';
COMMENT ON TABLE binding_assay_protocols IS 'Standardized protocols for binding assays';

-- Create triggers for updating the updated_at timestamps
CREATE TRIGGER update_compounds_modtime
    BEFORE UPDATE ON compounds
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_descriptors_2d_modtime
    BEFORE UPDATE ON descriptors_2d
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_descriptors_3d_modtime
    BEFORE UPDATE ON descriptors_3d
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_molecular_fingerprints_modtime
    BEFORE UPDATE ON molecular_fingerprints
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_pharmacophore_features_modtime
    BEFORE UPDATE ON pharmacophore_features
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_electronic_structure_modtime
    BEFORE UPDATE ON electronic_structure
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_quantum_critical_params_modtime
    BEFORE UPDATE ON quantum_critical_params
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_quantum_dynamics_modtime
    BEFORE UPDATE ON quantum_dynamics
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_phase_transitions_modtime
    BEFORE UPDATE ON phase_transitions
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_scaling_analysis_modtime
    BEFORE UPDATE ON scaling_analysis
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_quantum_observables_modtime
    BEFORE UPDATE ON quantum_observables
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_binding_assay_types_modtime
    BEFORE UPDATE ON binding_assay_types
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_binding_data_modtime
    BEFORE UPDATE ON binding_data
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_binding_data_quality_modtime
    BEFORE UPDATE ON binding_data_quality
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_binding_sar_modtime
    BEFORE UPDATE ON binding_sar
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_binding_kinetics_modtime
    BEFORE UPDATE ON binding_kinetics
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_binding_site_mapping_modtime
    BEFORE UPDATE ON binding_site_mapping
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_binding_assay_protocols_modtime
    BEFORE UPDATE ON binding_assay_protocols
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Create audit triggers for all tables
CREATE TRIGGER audit_compounds_trigger
    AFTER INSERT OR UPDATE OR DELETE ON compounds
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_descriptors_2d_trigger
    AFTER INSERT OR UPDATE OR DELETE ON descriptors_2d
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_descriptors_3d_trigger
    AFTER INSERT OR UPDATE OR DELETE ON descriptors_3d
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_molecular_fingerprints_trigger
    AFTER INSERT OR UPDATE OR DELETE ON molecular_fingerprints
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_pharmacophore_features_trigger
    AFTER INSERT OR UPDATE OR DELETE ON pharmacophore_features
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_electronic_structure_trigger
    AFTER INSERT OR UPDATE OR DELETE ON electronic_structure
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_quantum_critical_params_trigger
    AFTER INSERT OR UPDATE OR DELETE ON quantum_critical_params
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_quantum_dynamics_trigger
    AFTER INSERT OR UPDATE OR DELETE ON quantum_dynamics
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_phase_transitions_trigger
    AFTER INSERT OR UPDATE OR DELETE ON phase_transitions
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_scaling_analysis_trigger
    AFTER INSERT OR UPDATE OR DELETE ON scaling_analysis
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_quantum_observables_trigger
    AFTER INSERT OR UPDATE OR DELETE ON quantum_observables
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_binding_assay_types_trigger
    AFTER INSERT OR UPDATE OR DELETE ON binding_assay_types
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_binding_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON binding_data
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_binding_data_quality_trigger
    AFTER INSERT OR UPDATE OR DELETE ON binding_data_quality
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_binding_sar_trigger
    AFTER INSERT OR UPDATE OR DELETE ON binding_sar
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_binding_kinetics_trigger
    AFTER INSERT OR UPDATE OR DELETE ON binding_kinetics
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_binding_site_mapping_trigger
    AFTER INSERT OR UPDATE OR DELETE ON binding_site_mapping
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_binding_assay_protocols_trigger
    AFTER INSERT OR UPDATE OR DELETE ON binding_assay_protocols
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

-- Create indexes for compounds table
CREATE INDEX IF NOT EXISTS idx_compounds_name ON compounds(name);
CREATE INDEX IF NOT EXISTS idx_compounds_smiles ON compounds(smiles);
CREATE INDEX IF NOT EXISTS idx_compounds_inchi ON compounds(inchi);
CREATE INDEX IF NOT EXISTS idx_compounds_cas ON compounds(cas_number);
CREATE INDEX IF NOT EXISTS idx_compounds_pubchem ON compounds(pubchem_cid);
CREATE INDEX IF NOT EXISTS idx_compounds_chembl ON compounds(chembl_id);
CREATE INDEX IF NOT EXISTS idx_compounds_drugbank ON compounds(drugbank_id);
CREATE INDEX IF NOT EXISTS idx_compounds_created ON compounds(created_at);
CREATE INDEX IF NOT EXISTS idx_compounds_updated ON compounds(updated_at);

-- Create molecular_descriptors view for compatibility
CREATE OR REPLACE VIEW molecular_descriptors AS
SELECT 
    id,
    compound_id,
    wiener_index,
    balaban_j_index,
    bertz_ct,
    schultz_molecular_topological_index,
    chi_0_index,
    chi_1_index,
    chi_2_index,
    chi_3_index,
    sum_estate_indices,
    mean_estate_indices,
    kappa_1,
    kappa_2,
    kappa_3,
    molar_refractivity,
    van_der_waals_volume,
    polarizability,
    formal_charge,
    ring_count,
    aromatic_ring_count,
    aliphatic_ring_count,
    ring_fusion_degree,
    rotatable_bond_count,
    rigid_bond_count,
    chain_atom_count,
    chain_bond_count,
    carbon_count,
    nitrogen_count,
    oxygen_count,
    sulfur_count,
    phosphorus_count,
    halogen_count,
    created_at,
    updated_at
FROM descriptors_2d;

-- Insert common assay types
INSERT INTO binding_assay_types (name, description, method_type, detection_type, typical_unit) VALUES
('Radioligand Binding', 'Direct measurement using radioactive ligands', 'Competition', 'Scintillation', 'Ki (nM)'),
('Fluorescence Binding', 'Fluorescence-based binding assays', 'Direct', 'Fluorescence', 'Kd (nM)'),
('SPR', 'Surface plasmon resonance binding', 'Direct', 'Optical', 'KD (M)'),
('FRET', 'Förster resonance energy transfer', 'Proximity', 'Fluorescence', 'EC50 (nM)'),
('BRET', 'Bioluminescence resonance energy transfer', 'Proximity', 'Luminescence', 'EC50 (nM)')
ON CONFLICT DO NOTHING;
