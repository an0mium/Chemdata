-- Quantum chemistry and computational data tables

-- Reference data tables
-- Quantum criticality tables for tracking metal-insulator transitions

CREATE TABLE IF NOT EXISTS quantum_hamiltonians (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    -- Extended Hückel matrices
    eh_hamiltonian JSONB,  -- H_rs = (K/2)(H_rr + H_ss)S_rs
    overlap_matrix JSONB,  -- S_rs overlap matrix
    transformed_hamiltonian JSONB,  -- H_L = S^(-1/2)H(EH)S^(-1/2)
    -- Critical parameters
    disorder_strength FLOAT,  -- W value, critical at W_c ≈ 16.5
    is_critical BOOLEAN,  -- True if at metal-insulator transition
    -- Metadata
    calculation_method TEXT,  -- e.g. 'Extended Hückel'
    calculation_parameters JSONB,
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS energy_level_statistics (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    -- Level spacing data
    energy_levels FLOAT[],  -- E_i eigenvalues
    level_spacings FLOAT[], -- s_i = (E_i+1 - E_i)/Δ normalized spacings
    spacing_distribution JSONB, -- P(s) histogram data
    cumulative_distribution JSONB, -- I(S) = ∫P(s)ds
    -- Classification
    distribution_type TEXT, -- 'Poisson', 'Wigner', or 'Semi-Poisson'
    gamma_parameter FLOAT, -- γ for critical distribution
    confidence_score FLOAT, -- Statistical confidence in classification
    -- Metadata
    analysis_parameters JSONB,
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT valid_gamma CHECK (gamma_parameter >= 0)
);

CREATE TABLE IF NOT EXISTS wavefunction_analysis (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    -- Multifractal analysis
    box_probabilities JSONB, -- μ_k(l) box probabilities
    scaling_exponents JSONB, -- τ(q) scaling exponents
    fractal_dimensions JSONB, -- D_q generalized dimensions
    correlation_dimension FLOAT, -- D_2 special case
    -- Localization metrics
    localization_length FLOAT, -- ξ localization length
    participation_ratio FLOAT, -- Measure of state extension
    -- Metadata
    analysis_parameters JSONB,
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT valid_dimensions CHECK (
        correlation_dimension >= 0 AND
        correlation_dimension <= 1
    )
);

CREATE TABLE IF NOT EXISTS quantum_basis_sets (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,  -- Add UNIQUE constraint for foreign key references
    description TEXT,
    basis_type TEXT,  -- minimal, double-zeta, triple-zeta, etc.
    elements TEXT[],  -- supported elements
    accuracy_level TEXT,
    computational_cost TEXT,
    reference_citation TEXT
);
COMMENT ON TABLE quantum_basis_sets IS 'Standard basis sets for quantum calculations';

CREATE TABLE IF NOT EXISTS quantum_functionals (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,  -- Add UNIQUE constraint for foreign key references
    description TEXT,
    functional_type TEXT,  -- LDA, GGA, hybrid, etc.
    properties_handled TEXT[],  -- exchange, correlation, etc.
    accuracy_metrics JSONB,
    computational_cost TEXT,
    reference_citation TEXT
);
COMMENT ON TABLE quantum_functionals IS 'Density functionals for electronic structure calculations';

CREATE TABLE IF NOT EXISTS quantum_observables_ref (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,  -- Add UNIQUE constraint for foreign key references
    description TEXT,
    operator_type TEXT,  -- energy, momentum, spin, etc.
    measurement_units TEXT,
    uncertainty_type TEXT,
    standard_deviation_typical FLOAT,
    measurement_protocol TEXT
);
COMMENT ON TABLE quantum_observables_ref IS 'Reference data for quantum mechanical observables';

CREATE TABLE IF NOT EXISTS phase_transition_types (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,  -- Add UNIQUE constraint for foreign key references
    description TEXT,
    transition_order INTEGER,  -- 1st order, 2nd order, etc.
    critical_exponents JSONB,  -- α, β, γ, δ, etc.
    universality_class TEXT,
    characteristic_properties TEXT[]
);
COMMENT ON TABLE phase_transition_types IS 'Classification of phase transitions and their properties';

-- Electronic structure data
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
    orbital_energies FLOAT[],  -- Array of orbital energy levels
    orbital_occupancies INTEGER[],  -- Electron occupancy per orbital
    -- Electron density
    density_matrix JSONB,  -- Stored as sparse matrix
    density_grid_points JSONB,  -- 3D grid points for density
    density_values FLOAT[],  -- Electron density values at grid points
    -- Wavefunction data
    wavefunction JSONB,  -- Store as complex numbers
    state_type VARCHAR(50),  -- ground, excited, etc.
    -- Band structure
    band_structure JSONB,  -- For periodic systems
    -- Calculation metadata
    method TEXT REFERENCES quantum_functionals(name),  -- Computational method used (DFT, MP2, etc.)
    basis_set TEXT REFERENCES quantum_basis_sets(name),  -- Quantum chemistry basis set
    convergence_criteria JSONB,  -- Convergence parameters
    metadata JSONB,  -- Additional metadata
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT valid_energies CHECK (
        total_electronic_energy < 0 AND  -- Total energy should be negative
        homo_lumo_gap >= 0  -- HOMO-LUMO gap must be non-negative
    )
);

-- Quantum critical parameters
CREATE TABLE IF NOT EXISTS quantum_critical_params (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    -- Critical point parameters
    critical_temperature FLOAT,
    critical_pressure FLOAT,
    critical_field FLOAT,
    -- Order parameters
    primary_order_parameter TEXT,  -- Type of order parameter
    order_parameter_values JSONB,  -- Values at different conditions
    -- Critical exponents
    alpha FLOAT,  -- Specific heat exponent
    beta FLOAT,   -- Order parameter exponent
    gamma FLOAT,  -- Susceptibility exponent
    delta FLOAT,  -- Critical isotherm exponent
    nu FLOAT,     -- Correlation length exponent
    eta FLOAT,    -- Anomalous dimension
    -- Correlation functions
    correlation_length FLOAT,
    correlation_function JSONB,  -- Spatial correlation data
    dynamic_exponent_z FLOAT,    -- Dynamic critical exponent
    -- Phase diagram
    phase_boundaries JSONB,      -- Phase transition lines
    multicriticality_type TEXT,  -- Type of multicritical behavior
    -- Quantum properties
    coherence_length FLOAT,      -- Quantum coherence length
    entanglement_entropy FLOAT,  -- Quantum entanglement measure
    quantum_fluctuations JSONB,  -- Quantum fluctuation characteristics
    metadata JSONB,  -- Additional metadata
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT valid_critical_params CHECK (
        critical_temperature >= 0 AND
        critical_pressure >= 0 AND
        correlation_length > 0
    )
);

-- Quantum dynamics
CREATE TABLE IF NOT EXISTS quantum_dynamics (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    -- Time evolution
    time_points FLOAT[],         -- Time points for dynamics
    wavefunction_evolution JSONB, -- Time-dependent wavefunction
    density_matrix_evolution JSONB, -- Time-dependent density matrix
    -- Dynamic properties
    coherence_times FLOAT[],     -- Quantum coherence timescales
    relaxation_rates FLOAT[],    -- Energy relaxation rates
    dephasing_rates FLOAT[],     -- Quantum dephasing rates
    -- Transport properties
    conductivity_tensor FLOAT[], -- Quantum conductivity
    hall_conductance FLOAT,     -- Hall effect
    thermal_conductivity FLOAT,  -- Thermal transport
    -- Spectral properties
    spectral_function JSONB,    -- Many-body spectral function
    optical_conductivity JSONB,  -- Frequency-dependent conductivity
    -- Quantum correlations
    entanglement_spectrum FLOAT[], -- Entanglement eigenvalues
    mutual_information FLOAT,     -- Quantum mutual information
    -- Environmental coupling
    dissipation_kernel JSONB,    -- System-bath coupling
    noise_spectrum JSONB,        -- Environmental noise
    metadata JSONB,  -- Additional metadata
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT valid_rates CHECK (
        NOT EXISTS (
            SELECT 1 FROM unnest(relaxation_rates) AS rate WHERE rate < 0
        ) AND
        NOT EXISTS (
            SELECT 1 FROM unnest(dephasing_rates) AS rate WHERE rate < 0
        )
    )
);

-- Phase transitions
CREATE TABLE IF NOT EXISTS phase_transitions (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    transition_type TEXT NOT NULL REFERENCES phase_transition_types(name),
    critical_temperature FLOAT,
    critical_pressure FLOAT,
    order_parameter JSONB,  -- Order parameter data
    correlation_length FLOAT,
    -- Additional fields
    transition_order INTEGER,  -- First order, second order, etc.
    hysteresis_data JSONB,    -- Hysteresis measurements
    fluctuation_data JSONB,   -- Critical fluctuations
    metadata JSONB,  -- Additional metadata
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT valid_transition_params CHECK (
        critical_temperature >= 0 AND
        critical_pressure >= 0 AND
        correlation_length > 0 AND
        transition_order > 0
    )
);

-- Scaling analysis
CREATE TABLE IF NOT EXISTS scaling_analysis (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    phase_transition_id uuid REFERENCES phase_transitions(id),
    -- Scaling functions
    scaling_function_type TEXT,   -- Type of scaling analysis
    scaling_variables JSONB,      -- Relevant scaling variables
    scaling_dimensions FLOAT[],   -- Scaling dimensions
    -- RG analysis
    rg_flow_equations JSONB,      -- RG flow characteristics
    fixed_points JSONB,           -- RG fixed points
    relevant_operators JSONB,     -- Relevant perturbations
    -- Universal properties
    universality_class TEXT,      -- Universality classification
    central_charge FLOAT,         -- Conformal central charge
    operator_spectrum JSONB,      -- Spectrum of scaling operators
    -- Finite-size scaling
    size_scaling_exponents FLOAT[], -- Finite-size scaling
    correction_exponents FLOAT[],   -- Scaling corrections
    -- Crossover behavior
    crossover_scales JSONB,       -- Characteristic scales
    crossover_functions JSONB,    -- Crossover functions
    -- Metadata
    analysis_method TEXT,         -- Method used for analysis
    confidence_metrics JSONB,     -- Quality of scaling collapse
    metadata JSONB,  -- Additional metadata
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP
);

-- Quantum observables
CREATE TABLE IF NOT EXISTS quantum_observables (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    observable_type TEXT NOT NULL REFERENCES quantum_observables_ref(name),
    value FLOAT,
    uncertainty FLOAT,
    -- Additional fields
    measurement_basis TEXT,    -- Basis for measurement
    operator_type TEXT,        -- Type of quantum operator
    expectation_value FLOAT,   -- Quantum expectation value
    variance FLOAT,            -- Quantum uncertainty
    metadata JSONB,  -- Additional metadata
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT valid_uncertainty CHECK (uncertainty >= 0)
);

-- Create indexes
CREATE INDEX IF NOT EXISTS idx_electronic_structure_compound ON electronic_structure(compound_id);
CREATE INDEX IF NOT EXISTS idx_quantum_critical_params_compound ON quantum_critical_params(compound_id);
CREATE INDEX IF NOT EXISTS idx_quantum_dynamics_compound ON quantum_dynamics(compound_id);
CREATE INDEX IF NOT EXISTS idx_phase_transitions_compound ON phase_transitions(compound_id);
CREATE INDEX IF NOT EXISTS idx_scaling_analysis_compound ON scaling_analysis(compound_id);
CREATE INDEX IF NOT EXISTS idx_scaling_analysis_transition ON scaling_analysis(phase_transition_id);
CREATE INDEX IF NOT EXISTS idx_quantum_observables_compound ON quantum_observables(compound_id);

-- Create indexes for reference tables
CREATE INDEX IF NOT EXISTS idx_quantum_basis_sets_type ON quantum_basis_sets(basis_type);
CREATE INDEX IF NOT EXISTS idx_quantum_functionals_type ON quantum_functionals(functional_type);
CREATE INDEX IF NOT EXISTS idx_quantum_observables_type ON quantum_observables_ref(operator_type);
CREATE INDEX IF NOT EXISTS idx_phase_transitions_order ON phase_transition_types(transition_order);

-- Create GiST indexes for JSONB fields
CREATE INDEX IF NOT EXISTS idx_electronic_structure_density ON electronic_structure USING gin (density_matrix);
CREATE INDEX IF NOT EXISTS idx_quantum_critical_params_boundaries ON quantum_critical_params USING gin (phase_boundaries);
CREATE INDEX IF NOT EXISTS idx_quantum_dynamics_evolution ON quantum_dynamics USING gin (wavefunction_evolution);
CREATE INDEX IF NOT EXISTS idx_phase_transitions_order ON phase_transitions USING gin (order_parameter);
CREATE INDEX IF NOT EXISTS idx_scaling_analysis_flow ON scaling_analysis USING gin (rg_flow_equations);

-- Add triggers for updated_at
DROP TRIGGER IF EXISTS update_electronic_structure_modtime ON electronic_structure;
DROP TRIGGER IF EXISTS update_quantum_critical_params_modtime ON quantum_critical_params;
DROP TRIGGER IF EXISTS update_quantum_dynamics_modtime ON quantum_dynamics;
DROP TRIGGER IF EXISTS update_phase_transitions_modtime ON phase_transitions;
DROP TRIGGER IF EXISTS update_scaling_analysis_modtime ON scaling_analysis;
DROP TRIGGER IF EXISTS update_quantum_observables_modtime ON quantum_observables;

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

-- Add audit triggers
DROP TRIGGER IF EXISTS audit_electronic_structure_trigger ON electronic_structure;
DROP TRIGGER IF EXISTS audit_quantum_critical_params_trigger ON quantum_critical_params;
DROP TRIGGER IF EXISTS audit_quantum_dynamics_trigger ON quantum_dynamics;
DROP TRIGGER IF EXISTS audit_phase_transitions_trigger ON phase_transitions;
DROP TRIGGER IF EXISTS audit_scaling_analysis_trigger ON scaling_analysis;
DROP TRIGGER IF EXISTS audit_quantum_observables_trigger ON quantum_observables;

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

-- Quantum Properties (Consolidated View)
CREATE TABLE IF NOT EXISTS quantum_properties (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    -- Electronic structure properties
    total_electronic_energy FLOAT,
    homo_lumo_gap FLOAT,
    electron_density JSONB,
    orbital_energies FLOAT[],
    -- Critical behavior
    critical_temperature FLOAT,
    critical_pressure FLOAT,
    correlation_length FLOAT,
    order_parameter JSONB,
    -- Quantum dynamics
    coherence_time FLOAT,
    relaxation_rate FLOAT,
    dephasing_rate FLOAT,
    -- Phase transition properties
    transition_type TEXT,
    transition_order INTEGER,
    universality_class TEXT,
    -- Calculation metadata
    basis_set TEXT REFERENCES quantum_basis_sets(name),
    method TEXT REFERENCES quantum_functionals(name),
    calculation_parameters JSONB,
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT valid_quantum_ranges CHECK (
        total_electronic_energy < 0 AND
        homo_lumo_gap >= 0 AND
        critical_temperature >= 0 AND
        critical_pressure >= 0 AND
        correlation_length > 0 AND
        coherence_time >= 0 AND
        relaxation_rate >= 0 AND
        dephasing_rate >= 0 AND
        transition_order > 0
    )
);

-- Create indexes for quantum_properties
CREATE INDEX IF NOT EXISTS idx_quantum_properties_compound ON quantum_properties(compound_id);
CREATE INDEX IF NOT EXISTS idx_quantum_properties_basis ON quantum_properties(basis_set);
CREATE INDEX IF NOT EXISTS idx_quantum_properties_method ON quantum_properties(method);
CREATE INDEX IF NOT EXISTS idx_quantum_properties_transition ON quantum_properties(transition_type);

-- Create triggers for quantum_properties
CREATE TRIGGER update_quantum_properties_modtime
    BEFORE UPDATE ON quantum_properties
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER audit_quantum_properties_trigger
    AFTER INSERT OR UPDATE OR DELETE ON quantum_properties
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

-- Add table comments
COMMENT ON TABLE quantum_properties IS 'Consolidated quantum mechanical properties for compounds';
COMMENT ON TABLE electronic_structure IS 'Quantum mechanical electronic structure data';
COMMENT ON TABLE quantum_critical_params IS 'Parameters characterizing quantum critical behavior';
COMMENT ON TABLE quantum_dynamics IS 'Time-dependent quantum properties and dynamics';
COMMENT ON TABLE phase_transitions IS 'Phase transition data and properties';
COMMENT ON TABLE scaling_analysis IS 'Scaling analysis and renormalization group results';
COMMENT ON TABLE quantum_observables IS 'Quantum observable measurements and uncertainties';

-- Add column comments
COMMENT ON COLUMN electronic_structure.wavefunction IS 'Wavefunction data stored as complex numbers in JSONB format';
COMMENT ON COLUMN electronic_structure.density_matrix IS 'Electron density matrix stored as sparse matrix in JSONB format';
COMMENT ON COLUMN quantum_critical_params.order_parameter_values IS 'Order parameter data at different conditions stored in JSONB format';
COMMENT ON COLUMN phase_transitions.order_parameter IS 'Order parameter data stored in JSONB format';
COMMENT ON COLUMN scaling_analysis.rg_flow_equations IS 'Renormalization group flow equations stored in JSONB format';

-- Insert standard basis sets
INSERT INTO quantum_basis_sets (name, description, basis_type, elements, accuracy_level, computational_cost) VALUES
('STO-3G', 'Minimal basis set using 3 Gaussian functions', 'minimal', ARRAY['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl'], 'low', 'very low'),
('6-31G', 'Split valence basis set', 'double-zeta', ARRAY['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl'], 'medium', 'medium'),
('cc-pVDZ', 'Correlation consistent double-zeta basis', 'double-zeta', ARRAY['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl'], 'high', 'high'),
('cc-pVTZ', 'Correlation consistent triple-zeta basis', 'triple-zeta', ARRAY['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl'], 'very high', 'very high');

-- Insert standard density functionals
INSERT INTO quantum_functionals (name, description, functional_type, properties_handled, computational_cost) VALUES
('LDA', 'Local Density Approximation', 'LDA', ARRAY['exchange', 'correlation'], 'low'),
('PBE', 'Perdew-Burke-Ernzerhof', 'GGA', ARRAY['exchange', 'correlation'], 'medium'),
('B3LYP', 'Becke 3-parameter Lee-Yang-Parr', 'hybrid', ARRAY['exchange', 'correlation'], 'high'),
('M06-2X', 'Minnesota 06 functional with 2X exchange', 'hybrid-meta-GGA', ARRAY['exchange', 'correlation', 'dispersion'], 'very high');

-- Insert standard quantum observables
INSERT INTO quantum_observables_ref (name, description, operator_type, measurement_units, uncertainty_type) VALUES
('Energy', 'Total electronic energy', 'energy', 'Hartree', 'absolute'),
('Dipole Moment', 'Electric dipole moment', 'electromagnetic', 'Debye', 'relative'),
('Spin', 'Electron spin', 'angular momentum', 'ℏ/2', 'discrete'),
('Electron Density', 'Probability density of electrons', 'density', 'e/Å³', 'statistical');

-- Insert phase transition types
INSERT INTO phase_transition_types (name, description, transition_order, critical_exponents, universality_class) VALUES
('Continuous', 'Second-order phase transition', 2, '{"alpha": 0.110, "beta": 0.325, "gamma": 1.24, "delta": 4.82, "nu": 0.63}', 'Ising'),
('Discontinuous', 'First-order phase transition', 1, NULL, NULL),
('BCS', 'Superconducting transition', 2, '{"alpha": 0, "beta": 0.5, "gamma": 1.0, "delta": 3.0, "nu": 0.5}', 'Mean Field'),
('BEC', 'Bose-Einstein condensation', 2, '{"alpha": -1, "beta": 0.5, "gamma": 1.0, "delta": 3.0, "nu": 0.5}', 'Gaussian');

-- Functions for quantum criticality analysis

-- Function to analyze level spacing statistics and classify compounds
CREATE OR REPLACE FUNCTION analyze_level_spacing(
    compound_id uuid,
    min_levels integer DEFAULT 1000  -- Minimum number of levels needed for reliable statistics
) RETURNS TABLE (
    distribution_type TEXT,  -- 'Poisson', 'Wigner', or 'Semi-Poisson'
    gamma_value FLOAT,      -- γ parameter for critical distribution
    confidence_score FLOAT, -- Statistical confidence in classification
    chi_squared FLOAT      -- Goodness of fit measure
) AS $$
DECLARE
    spacing_count INTEGER;
    mean_spacing FLOAT;
    spacing_variance FLOAT;
BEGIN
    -- Get normalized level spacings
    WITH spacings AS (
        SELECT level_spacings 
        FROM energy_level_statistics 
        WHERE compound_id = $1
    )
    SELECT 
        array_length(level_spacings, 1),
        avg(unnest),
        variance(unnest)
    INTO spacing_count, mean_spacing, spacing_variance
    FROM spacings, unnest(level_spacings);

    -- Check if we have enough statistics
    IF spacing_count < min_levels THEN
        RETURN QUERY SELECT 
            'insufficient_data'::TEXT,
            NULL::FLOAT,
            0::FLOAT,
            NULL::FLOAT;
        RETURN;
    END IF;

    -- Calculate fit scores for each distribution
    RETURN QUERY
    WITH distribution_fits AS (
        SELECT
            -- Poisson: P(s) = exp(-s)
            sum((hist_count - exp(-bin_center))^2) as poisson_chi2,
            -- Wigner: P(s) = (π/2)s*exp(-πs²/4)
            sum((hist_count - (pi()/2)*bin_center*exp(-pi()*bin_center^2/4))^2) as wigner_chi2,
            -- Semi-Poisson: P(s) = 4s*exp(-2s)
            sum((hist_count - 4*bin_center*exp(-2*bin_center))^2) as semi_poisson_chi2
        FROM (
            SELECT 
                width_bucket(unnest(level_spacings), 0, 5, 50) as bin,
                count(*)::float/spacing_count as hist_count,
                (width_bucket(unnest(level_spacings), 0, 5, 50) - 0.5)*0.1 as bin_center
            FROM energy_level_statistics
            WHERE compound_id = $1
            GROUP BY bin
        ) as histogram
    )
    SELECT
        CASE
            WHEN poisson_chi2 < LEAST(wigner_chi2, semi_poisson_chi2) THEN 'Poisson'
            WHEN wigner_chi2 < LEAST(poisson_chi2, semi_poisson_chi2) THEN 'Wigner'
            ELSE 'Semi-Poisson'
        END,
        CASE 
            WHEN semi_poisson_chi2 < LEAST(poisson_chi2, wigner_chi2) THEN 0.0  -- γ ≈ 0 for critical
            ELSE NULL
        END,
        1.0 - LEAST(poisson_chi2, wigner_chi2, semi_poisson_chi2)/
            (poisson_chi2 + wigner_chi2 + semi_poisson_chi2),
        LEAST(poisson_chi2, wigner_chi2, semi_poisson_chi2)
    FROM distribution_fits;
END;
$$ LANGUAGE plpgsql;

-- Function to perform multifractal analysis of wavefunctions
CREATE OR REPLACE FUNCTION analyze_wavefunction_multifractality(
    compound_id uuid,
    q_values integer[] DEFAULT ARRAY[-10,-5,-2,-1,0,1,2,5,10],  -- Generalized dimension orders
    min_box_size integer DEFAULT 10,                            -- Minimum box size for scaling
    max_box_size integer DEFAULT 1000                          -- Maximum box size for scaling
) RETURNS TABLE (
    q INTEGER,              -- Order of generalized dimension
    tau FLOAT,             -- Mass exponent τ(q)
    dimension FLOAT,        -- Generalized dimension D_q
    r_squared FLOAT,       -- Fit quality for scaling
    is_critical BOOLEAN    -- True if D_2 ≈ 0.5 (critical)
) AS $$
BEGIN
    RETURN QUERY
    WITH box_counting AS (
        -- Calculate box probabilities μ_k(l) at different scales
        SELECT 
            q_val,
            box_size,
            sum(power(prob, q_val)) as moment
        FROM (
            SELECT 
                q.q as q_val,
                l.l as box_size,
                sum(power(abs(coeff), 2)) as prob
            FROM unnest(q_values) as q(q),
                 generate_series(min_box_size, max_box_size, min_box_size) as l(l),
                 (SELECT wavefunction FROM electronic_structure WHERE compound_id = $1) as wf,
                 jsonb_array_elements(wf) as coeff
            GROUP BY q.q, l.l
        ) as boxes
        GROUP BY q_val, box_size
    ),
    scaling_fits AS (
        -- Fit scaling exponents τ(q)
        SELECT
            q_val,
            regr_slope(ln(moment), ln(box_size)) as tau,
            regr_r2(ln(moment), ln(box_size)) as r2
        FROM box_counting
        GROUP BY q_val
    )
    SELECT 
        q_val,
        tau,
        CASE 
            WHEN q_val = 1 THEN tau
            ELSE tau/(q_val - 1)
        END as dq,
        r2,
        CASE
            WHEN q_val = 2 THEN abs(tau/(q_val - 1) - 0.5) < 0.05  -- D_2 ≈ 0.5 criterion
            ELSE NULL
        END
    FROM scaling_fits
    ORDER BY q_val;
END;
$$ LANGUAGE plpgsql;

-- Function to calculate quantum criticality indicators
CREATE OR REPLACE FUNCTION calculate_quantum_criticality_indicators(
    compound_id uuid
) RETURNS TABLE (
    indicator_name TEXT,
    indicator_value FLOAT,
    confidence_score FLOAT
) AS $$
BEGIN
    RETURN QUERY
    WITH electronic_data AS (
        SELECT 
            homo_lumo_gap,
            total_electronic_energy,
            electron_correlation_energy
        FROM electronic_structure
        WHERE compound_id = $1
    ),
    critical_data AS (
        SELECT 
            correlation_length,
            entanglement_entropy
        FROM quantum_critical_params
        WHERE compound_id = $1
    ),
    level_stats AS (
        SELECT * FROM analyze_level_spacing($1)
    ),
    fractal_analysis AS (
        SELECT * FROM analyze_wavefunction_multifractality($1)
        WHERE q = 2  -- Get correlation dimension D_2
    )
    SELECT 'level_spacing'::TEXT,
           CASE 
               WHEN distribution_type = 'Semi-Poisson' THEN 1.0
               WHEN distribution_type = 'Poisson' THEN 0.0
               ELSE 0.5
           END,
           confidence_score
    FROM level_stats
    UNION ALL
    SELECT 'correlation_dimension',
           dimension,
           r_squared
    FROM fractal_analysis
    UNION ALL
    SELECT 'energy_gap',
           homo_lumo_gap,
           CASE 
               WHEN homo_lumo_gap > 0 THEN 1.0
               ELSE 0.5
           END
    FROM electronic_data
    UNION ALL
    SELECT 'correlation_strength',
           correlation_length,
           CASE 
               WHEN correlation_length > 10 THEN 1.0
               WHEN correlation_length > 5 THEN 0.8
               ELSE 0.6
           END
    FROM critical_data;
END;
$$ LANGUAGE plpgsql;

-- Trigger function to update quantum properties
CREATE OR REPLACE FUNCTION update_quantum_properties()
RETURNS TRIGGER AS $$
BEGIN
    -- Update electronic structure when critical parameters change
    IF TG_TABLE_NAME = 'quantum_critical_params' THEN
        UPDATE electronic_structure
        SET updated_at = CURRENT_TIMESTAMP
        WHERE compound_id = NEW.compound_id;
    END IF;
    
    -- Update scaling analysis when dynamics change
    IF TG_TABLE_NAME = 'quantum_dynamics' THEN
        UPDATE scaling_analysis
        SET updated_at = CURRENT_TIMESTAMP
        WHERE compound_id = NEW.compound_id;
    END IF;
    
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Create triggers for quantum property updates
CREATE TRIGGER trigger_quantum_critical_updates
AFTER INSERT OR UPDATE ON quantum_critical_params
FOR EACH ROW EXECUTE FUNCTION update_quantum_properties();

CREATE TRIGGER trigger_quantum_dynamics_updates
AFTER INSERT OR UPDATE ON quantum_dynamics
FOR EACH ROW EXECUTE FUNCTION update_quantum_properties();


-- Add foreign key constraints
ALTER TABLE electronic_structure
ADD CONSTRAINT fk_basis_set FOREIGN KEY (basis_set) REFERENCES quantum_basis_sets(name),
ADD CONSTRAINT fk_functional FOREIGN KEY (method) REFERENCES quantum_functionals(name);

ALTER TABLE quantum_observables
ADD CONSTRAINT fk_observable_type FOREIGN KEY (observable_type) REFERENCES quantum_observables_ref(name);

ALTER TABLE phase_transitions
ADD CONSTRAINT fk_transition_type FOREIGN KEY (transition_type) REFERENCES phase_transition_types(name);


-- Research data tables
CREATE TABLE IF NOT EXISTS quantum_research_projects (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT NOT NULL,
    description TEXT,
    start_date DATE,
    end_date DATE,
    status TEXT,
    principal_investigator TEXT,
    research_type TEXT,  -- theoretical, computational, experimental
    funding_source TEXT,
    budget DECIMAL,
    objectives TEXT[],
    metadata JSONB,
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP
);
COMMENT ON TABLE quantum_research_projects IS 'Research projects involving quantum mechanical studies';

CREATE TABLE IF NOT EXISTS quantum_calculations (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    project_id uuid REFERENCES quantum_research_projects(id),
    compound_id uuid REFERENCES compounds(id),
    calculation_type TEXT,  -- geometry optimization, frequency analysis, etc.
    basis_set TEXT REFERENCES quantum_basis_sets(name),
    functional TEXT REFERENCES quantum_functionals(name),
    calculation_status TEXT,
    start_timestamp TIMESTAMP,
    end_timestamp TIMESTAMP,
    cpu_hours FLOAT,
    memory_gb FLOAT,
    convergence_achieved BOOLEAN,
    energy_hartree FLOAT,
    energy_gradient_norm FLOAT,
    calculation_parameters JSONB,
    output_files TEXT[],
    metadata JSONB,
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT check_energy_hartree CHECK (energy_hartree < 0),
    CONSTRAINT check_energy_gradient CHECK (energy_gradient_norm >= 0),
    CONSTRAINT check_resources CHECK (cpu_hours > 0 AND memory_gb > 0)
);
COMMENT ON TABLE quantum_calculations IS 'Quantum mechanical calculations and their parameters';

CREATE TABLE IF NOT EXISTS quantum_research_findings (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    calculation_id uuid REFERENCES quantum_calculations(id),
    finding_type TEXT,  -- structural, electronic, spectroscopic, etc.
    description TEXT,
    numerical_value FLOAT,
    units TEXT,
    confidence_level FLOAT,
    methodology TEXT,
    validation_method TEXT,
    publication_reference TEXT,
    discovery_date DATE,
    significance_level TEXT,  -- high, medium, low
    metadata JSONB,
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT check_confidence CHECK (confidence_level BETWEEN 0 AND 1),
    CONSTRAINT check_significance CHECK (significance_level IN ('high', 'medium', 'low'))
);
COMMENT ON TABLE quantum_research_findings IS 'Research findings from quantum calculations';

CREATE TABLE IF NOT EXISTS quantum_structure_correlations (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid REFERENCES compounds(id),
    property_type TEXT,  -- electronic, geometric, energetic
    property_value FLOAT,
    correlation_type TEXT,  -- linear, polynomial, exponential
    correlation_coefficient FLOAT,
    statistical_significance FLOAT,
    sample_size INTEGER,
    methodology TEXT,
    validation_metrics JSONB,
    metadata JSONB,
    created_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT check_correlation CHECK (correlation_coefficient BETWEEN -1 AND 1),
    CONSTRAINT check_significance_value CHECK (statistical_significance BETWEEN 0 AND 1)
);
COMMENT ON TABLE quantum_structure_correlations IS 'Structure-property correlations from quantum studies';

-- Create additional indexes for research tables
CREATE INDEX idx_quantum_calculations_project ON quantum_calculations(project_id);
CREATE INDEX idx_quantum_calculations_compound ON quantum_calculations(compound_id);
CREATE INDEX idx_quantum_calculations_basis ON quantum_calculations(basis_set);
CREATE INDEX idx_quantum_calculations_functional ON quantum_calculations(functional);
CREATE INDEX idx_quantum_findings_type ON quantum_research_findings(finding_type);
CREATE INDEX idx_quantum_correlations_property ON quantum_structure_correlations(property_type);

-- Add triggers for updated_at
CREATE TRIGGER update_quantum_research_projects_modtime
    BEFORE UPDATE ON quantum_research_projects
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_quantum_calculations_modtime
    BEFORE UPDATE ON quantum_calculations
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_quantum_research_findings_modtime
    BEFORE UPDATE ON quantum_research_findings
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_quantum_structure_correlations_modtime
    BEFORE UPDATE ON quantum_structure_correlations
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Add audit triggers
CREATE TRIGGER audit_quantum_research_projects_trigger
    AFTER INSERT OR UPDATE OR DELETE ON quantum_research_projects
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_quantum_calculations_trigger
    AFTER INSERT OR UPDATE OR DELETE ON quantum_calculations
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_quantum_research_findings_trigger
    AFTER INSERT OR UPDATE OR DELETE ON quantum_research_findings
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_quantum_structure_correlations_trigger
    AFTER INSERT OR UPDATE OR DELETE ON quantum_structure_correlations
    FOR EACH ROW
    EXECUTE FUNCTION audit_trigger_func();

-- Function to analyze Environment-assisted Quantum Transport (ENAQT)
CREATE OR REPLACE FUNCTION analyze_enaqt_properties(
    compound_id uuid
) RETURNS TABLE (
    transport_efficiency FLOAT,      -- Overall quantum transport efficiency
    coherence_time FLOAT,           -- Quantum coherence timescale
    decoherence_rate FLOAT,         -- Environmental decoherence rate
    anti_zeno_factor FLOAT,         -- Suppression of anti-Zeno effect
    is_optimal BOOLEAN              -- True if transport is optimal at criticality
) AS $$
BEGIN
    RETURN QUERY
    WITH dynamics_data AS (
        SELECT 
            qd.coherence_times[1] as coh_time,
            qd.dephasing_rates[1] as deph_rate,
            qd.conductivity_tensor[1] as conductivity,
            els.distribution_type,
            wa.correlation_dimension as d2,
            qh.disorder_strength
        FROM quantum_dynamics qd
        LEFT JOIN energy_level_statistics els ON qd.compound_id = els.compound_id
        LEFT JOIN wavefunction_analysis wa ON qd.compound_id = wa.compound_id
        LEFT JOIN quantum_hamiltonians qh ON qd.compound_id = qh.compound_id
        WHERE qd.compound_id = $1
    )
    SELECT
        conductivity as transport_efficiency,
        coh_time as coherence_time,
        deph_rate as decoherence_rate,
        CASE 
            WHEN distribution_type = 'Semi-Poisson' 
                 AND abs(disorder_strength - 16.5) < 0.5  -- Critical disorder strength
            THEN 1.0
            ELSE deph_rate * coh_time
        END as anti_zeno_factor,
        CASE
            -- Optimal transport conditions from paper:
            -- 1. Semi-Poissonian level statistics (critical)
            -- 2. D2 ≈ 0.5 (multifractal)
            -- 3. W ≈ 16.5 (Anderson transition)
            -- 4. Algebraic coherence decay
            WHEN distribution_type = 'Semi-Poisson' 
                 AND abs(d2 - 0.5) < 0.05
                 AND abs(disorder_strength - 16.5) < 0.5
                 AND deph_rate * coh_time < 1.0  -- Coherent regime
            THEN true
            ELSE false
        END as is_optimal
    FROM dynamics_data;
END;
$$ LANGUAGE plpgsql;

-- Function to analyze quantum transport mechanism
CREATE OR REPLACE FUNCTION analyze_transport_mechanism(
    compound_id uuid
) RETURNS TABLE (
    mechanism_type TEXT,           -- 'ENAQT', 'Classical', or 'Quantum'
    transport_regime TEXT,         -- 'Coherent', 'Critical', or 'Localized'
    efficiency_score FLOAT,        -- Transport efficiency metric
    mechanism_confidence FLOAT     -- Confidence in mechanism classification
) AS $$
BEGIN
    RETURN QUERY
    WITH transport_data AS (
        SELECT 
            enaqt.*,
            els.distribution_type,
            wa.correlation_dimension as d2,
            qh.disorder_strength
        FROM analyze_enaqt_properties($1) enaqt
        LEFT JOIN energy_level_statistics els ON els.compound_id = $1
        LEFT JOIN wavefunction_analysis wa ON wa.compound_id = $1
        LEFT JOIN quantum_hamiltonians qh ON qh.compound_id = $1
    )
    SELECT
        -- Classify transport mechanism
        CASE
            WHEN is_optimal THEN 'ENAQT'
            WHEN distribution_type = 'Wigner' THEN 'Quantum'
            ELSE 'Classical'
        END as mechanism_type,
        
        -- Determine transport regime
        CASE
            WHEN distribution_type = 'Semi-Poisson' 
                 AND abs(d2 - 0.5) < 0.05 THEN 'Critical'
            WHEN decoherence_rate * coherence_time < 1.0 THEN 'Coherent'
            ELSE 'Localized'
        END as transport_regime,
        
        -- Calculate efficiency score
        transport_efficiency * 
        CASE
            WHEN is_optimal THEN 2.0  -- Boost score for optimal ENAQT
            ELSE 1.0
        END as efficiency_score,
        
        -- Calculate confidence
        CASE
            WHEN is_optimal THEN 1.0
            WHEN distribution_type = 'Semi-Poisson' 
                 OR abs(d2 - 0.5) < 0.05 THEN 0.8
            ELSE 0.6
        END as mechanism_confidence
    FROM transport_data;
END;
$$ LANGUAGE plpgsql;

-- Create views for quantum criticality analysis
CREATE OR REPLACE VIEW quantum_criticality_summary AS
SELECT 
    c.id as compound_id,
    c.name as compound_name,
    -- Level spacing analysis
    els.distribution_type,
    els.gamma_parameter,
    els.confidence_score as spacing_confidence,
    -- Multifractal analysis
    wa.correlation_dimension as d2,
    wa.participation_ratio,
    -- Critical parameters
    qh.disorder_strength,
    qh.is_critical as hamiltonian_critical,
    -- Consolidated criticality assessment
    CASE 
        WHEN els.distribution_type = 'Semi-Poisson' 
             AND abs(wa.correlation_dimension - 0.5) < 0.05
             AND qh.is_critical = true
        THEN true
        ELSE false
    END as is_critical,
    -- Confidence metrics
    (els.confidence_score + 
     CASE 
         WHEN abs(wa.correlation_dimension - 0.5) < 0.05 THEN 1.0
         ELSE 0.0
     END +
     CASE 
         WHEN qh.is_critical THEN 1.0
         ELSE 0.0
     END
    )/3.0 as criticality_confidence
FROM compounds c
LEFT JOIN energy_level_statistics els ON c.id = els.compound_id
LEFT JOIN wavefunction_analysis wa ON c.id = wa.compound_id
LEFT JOIN quantum_hamiltonians qh ON c.id = qh.compound_id;

COMMENT ON VIEW quantum_criticality_summary IS 
'Consolidated view of quantum criticality indicators for compounds based on:
1. Level spacing distribution (should be Semi-Poissonian)
2. Correlation dimension D2 ≈ 0.5 (multifractal criterion)
3. Critical disorder strength at metal-insulator transition';

-- Create views for quantum transport analysis
CREATE OR REPLACE VIEW quantum_transport_summary AS
SELECT 
    c.id as compound_id,
    c.name as compound_name,
    -- Criticality indicators
    qcs.distribution_type,
    qcs.d2 as correlation_dimension,
    qcs.is_critical,
    -- ENAQT properties
    enaqt.transport_efficiency,
    enaqt.coherence_time,
    enaqt.decoherence_rate,
    enaqt.anti_zeno_factor,
    enaqt.is_optimal as transport_optimal,
    -- Transport assessment
    CASE 
        WHEN qcs.is_critical AND enaqt.is_optimal THEN 'Optimal quantum transport'
        WHEN qcs.is_critical THEN 'Critical but suboptimal transport'
        WHEN enaqt.is_optimal THEN 'Optimal transport but not critical'
        ELSE 'Neither critical nor optimal'
    END as transport_classification,
    -- Confidence metrics
    qcs.criticality_confidence,
    (qcs.criticality_confidence + 
     CASE 
         WHEN enaqt.is_optimal THEN 1.0
         ELSE enaqt.transport_efficiency
     END)/2.0 as overall_confidence
FROM compounds c
LEFT JOIN quantum_criticality_summary qcs ON c.id = qcs.compound_id
LEFT JOIN LATERAL (SELECT * FROM analyze_enaqt_properties(c.id)) enaqt ON true;

COMMENT ON VIEW quantum_transport_summary IS 
'Consolidated view of quantum transport properties combining:
1. Quantum criticality indicators (level spacing, multifractality)
2. Environment-assisted Quantum Transport (ENAQT) metrics
3. Overall transport efficiency assessment';

-- Create views for calculation summaries
CREATE OR REPLACE VIEW quantum_calculation_summary AS
SELECT 
    c.id,
    c.compound_id,
    c.calculation_type,
    c.basis_set,
    c.functional,
    c.energy_hartree,
    c.convergence_achieved,
    COUNT(f.id) as finding_count,
    STRING_AGG(DISTINCT f.finding_type, ', ') as finding_types
FROM quantum_calculations c
LEFT JOIN quantum_research_findings f ON c.id = f.calculation_id
GROUP BY c.id;

CREATE OR REPLACE VIEW quantum_structure_insights AS
SELECT 
    c.compound_id,
    COUNT(DISTINCT c.id) as calculation_count,
    COUNT(DISTINCT f.id) as finding_count,
    COUNT(DISTINCT corr.id) as correlation_count,
    STRING_AGG(DISTINCT c.calculation_type, ', ') as calculation_types,
    STRING_AGG(DISTINCT f.finding_type, ', ') as finding_types
FROM quantum_calculations c
LEFT JOIN quantum_research_findings f ON c.id = f.calculation_id
LEFT JOIN quantum_structure_correlations corr ON c.compound_id = corr.compound_id
GROUP BY c.compound_id;

-- Create view combining transport mechanism and criticality analysis
CREATE OR REPLACE VIEW quantum_transport_mechanism_summary AS
SELECT 
    c.id as compound_id,
    c.name as compound_name,
    -- Transport mechanism analysis
    tm.mechanism_type,
    tm.transport_regime,
    tm.efficiency_score,
    tm.mechanism_confidence,
    -- Criticality indicators
    qcs.distribution_type,
    qcs.d2 as correlation_dimension,
    qcs.disorder_strength,
    qcs.is_critical,
    -- ENAQT properties
    enaqt.transport_efficiency,
    enaqt.coherence_time,
    enaqt.decoherence_rate,
    enaqt.anti_zeno_factor,
    enaqt.is_optimal as transport_optimal,
    -- Combined assessment
    CASE 
        WHEN tm.mechanism_type = 'ENAQT' AND qcs.is_critical THEN 'Optimal quantum transport at criticality'
        WHEN tm.mechanism_type = 'ENAQT' THEN 'Environment-assisted transport'
        WHEN qcs.is_critical THEN 'Critical but not environment-assisted'
        WHEN tm.mechanism_type = 'Quantum' THEN 'Quantum transport'
        ELSE 'Classical transport'
    END as transport_classification,
    -- Overall confidence
    (tm.mechanism_confidence + qcs.criticality_confidence)/2.0 as overall_confidence
FROM compounds c
LEFT JOIN LATERAL (SELECT * FROM analyze_transport_mechanism(c.id)) tm ON true
LEFT JOIN quantum_criticality_summary qcs ON c.id = qcs.compound_id
LEFT JOIN LATERAL (SELECT * FROM analyze_enaqt_properties(c.id)) enaqt ON true;

COMMENT ON VIEW quantum_transport_mechanism_summary IS 
'Comprehensive view of quantum transport properties combining:
1. Transport mechanism classification (ENAQT, Quantum, Classical)
2. Quantum criticality indicators
3. Environment-assisted transport metrics
4. Overall transport efficiency assessment';

-- Add reference data
INSERT INTO quantum_research_projects (name, description, research_type, status) VALUES
('QM Structure Analysis', 'Quantum mechanical analysis of molecular structures', 'computational', 'active'),
('Electronic Properties', 'Investigation of electronic structure properties', 'theoretical', 'active'),
('Reaction Mechanisms', 'Quantum study of reaction pathways', 'computational', 'planned');
