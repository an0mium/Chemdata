-- Consolidated pharmacology schema
-- Combines pharmacokinetics, pharmacodynamics, mechanisms, and effects

------------------------------------------
-- Core Pharmacology Data
------------------------------------------

-- Pharmacokinetic parameters
CREATE TABLE IF NOT EXISTS pharmacokinetic_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    parameter_type text NOT NULL, -- e.g., 'absorption', 'distribution', 'metabolism', 'excretion'
    species text,
    route_of_administration text,
    dose double precision,
    dose_unit text,
    value double precision,
    value_unit text,
    confidence_interval jsonb,
    study_type text,
    experimental_conditions jsonb,
    analysis_method text,
    reference_doi text,
    notes text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- ADME properties
CREATE TABLE IF NOT EXISTS adme_properties (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    property_type text NOT NULL,
    value double precision,
    unit text,
    prediction_method text,
    experimental_validation boolean,
    confidence_score double precision,
    data_source text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Pharmacodynamic data
CREATE TABLE IF NOT EXISTS pharmacodynamic_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    effect_type text NOT NULL,           -- e.g., 'Central', 'Peripheral'
    mechanism text NOT NULL,
    onset_characteristics jsonb,         -- Onset profile
    duration_profile jsonb,             -- Duration characteristics
    effect_magnitude jsonb,             -- Effect strength data
    dose_dependency jsonb,              -- Dose-response relationship
    tolerance_development jsonb,        -- Tolerance characteristics
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Molecular Biology Data
------------------------------------------

-- Receptor families table with enhanced fields
CREATE TABLE IF NOT EXISTS receptor_families (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    family_name text NOT NULL UNIQUE,
    description text,
    protein_type text,
    signaling_type text,
    receptor_class text,
    primary_effects text[],
    therapeutic_areas text[],
    known_compounds text[],
    abuse_potential text,
    safety_profile text,
    evolutionary_conservation text,
    tissue_distribution jsonb,
    regulatory_status text,
    cellular_location text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Enhanced receptor expression tracking
CREATE TABLE IF NOT EXISTS receptor_expression (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    receptor_family_id uuid REFERENCES receptor_families(id) ON DELETE CASCADE,
    tissue_type text NOT NULL,
    expression_level text CHECK (expression_level IN ('high', 'medium', 'low', 'not_detected')),
    confidence_score float CHECK (confidence_score BETWEEN 0 AND 1),
    data_source text,
    experimental_method text,
    species text,
    developmental_stage text,
    conditions jsonb,
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Enhanced receptor functions tracking
CREATE TABLE IF NOT EXISTS receptor_functions (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    receptor_family_id uuid REFERENCES receptor_families(id) ON DELETE CASCADE,
    function_name text NOT NULL,
    description text,
    pathway text,
    physiological_effect text,
    evidence_type text,
    mechanism_of_action text,
    regulatory_effects jsonb,
    cellular_responses text[],
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Enhanced receptor interactions tracking
CREATE TABLE IF NOT EXISTS receptor_interactions (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    receptor_family_id uuid REFERENCES receptor_families(id) ON DELETE CASCADE,
    interacting_protein text NOT NULL,
    interaction_type text,
    effect text,
    binding_region text,
    binding_affinity double precision,
    interaction_mechanism text,
    physiological_outcome text,
    evidence_source text,
    experimental_conditions jsonb,
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Protein-Receptor relationships
CREATE TABLE IF NOT EXISTS protein_receptor_relationships (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    protein_id uuid NOT NULL REFERENCES proteins(id),
    receptor_family_id uuid NOT NULL REFERENCES receptor_families(id),
    relationship_type text NOT NULL, -- e.g., 'subunit', 'accessory', 'modulator'
    role_description text,
    evidence_type text[],
    confidence_score double precision,
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Receptor variants and mutations
CREATE TABLE IF NOT EXISTS receptor_variants (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    subtype_id uuid REFERENCES receptor_subtypes(id) ON DELETE CASCADE,
    variant_name text NOT NULL,
    mutation_type text,
    sequence_change text,
    functional_impact text,
    population_frequency double precision,
    clinical_associations text[],
    phenotype_effects jsonb,
    structural_changes jsonb,
    pharmacological_changes jsonb,
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Receptor binding sites
CREATE TABLE IF NOT EXISTS subtype_binding_sites (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    subtype_id uuid REFERENCES receptor_subtypes(id) ON DELETE CASCADE,
    site_name text NOT NULL,
    location text,
    residues text[],
    binding_properties jsonb,
    allosteric_effects jsonb,
    structural_features jsonb,
    conservation_score double precision,
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Receptor signaling pathways
CREATE TABLE IF NOT EXISTS subtype_signaling (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    subtype_id uuid REFERENCES receptor_subtypes(id) ON DELETE CASCADE,
    pathway_name text NOT NULL,
    signaling_proteins text[],
    second_messengers text[],
    cellular_response text,
    temporal_profile jsonb,
    pathway_specificity double precision,
    regulatory_mechanisms jsonb,
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now()
);

-- Specialized receptor types
CREATE TABLE IF NOT EXISTS nuclear_receptor_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    receptor_id uuid REFERENCES receptor_families(id) ON DELETE CASCADE,
    dna_binding_domain text,
    ligand_binding_domain text,
    coregulator_interactions text[],
    response_elements text[],
    activation_mechanism text,
    tissue_specific_effects jsonb,
    post_translational_mods jsonb,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

CREATE TABLE IF NOT EXISTS peptide_receptor_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    receptor_id uuid REFERENCES receptor_families(id) ON DELETE CASCADE,
    endogenous_ligands text[],
    peptide_specificity jsonb,
    signaling_cascades text[],
    regulatory_mechanisms text[],
    physiological_roles text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Molecular Interactions
------------------------------------------

-- Receptor binding profiles (enhanced version)
CREATE TABLE IF NOT EXISTS receptor_binding_profiles (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    profile_type text NOT NULL,          -- e.g., 'Primary', 'Secondary'
    receptor_affinities jsonb,           -- Map of receptor to affinity values
    binding_ratios jsonb,               -- Relative binding ratios
    selectivity_data jsonb,             -- Receptor selectivity metrics
    functional_effects jsonb,           -- Functional responses
    methodology text,                   -- How profile was determined
    confidence_score double precision,
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Enzyme interactions
CREATE TABLE IF NOT EXISTS enzyme_interactions (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    enzyme_name text NOT NULL,
    interaction_type text NOT NULL, -- e.g., 'inhibition', 'induction', 'substrate'
    parameter_type text, -- e.g., 'Ki', 'IC50', 'Km'
    value double precision,
    unit text,
    conditions jsonb,
    experimental_method text,
    reference_doi text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Transporter interactions
CREATE TABLE IF NOT EXISTS transporter_interactions (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    transporter_name text NOT NULL,
    interaction_type text NOT NULL, -- e.g., 'substrate', 'inhibitor', 'inducer'
    parameter_type text,
    value double precision,
    unit text,
    conditions jsonb,
    experimental_method text,
    reference_doi text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Drug interactions (enhanced version)
CREATE TABLE IF NOT EXISTS drug_interactions (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    interacting_compound_id uuid NOT NULL REFERENCES compounds(id),
    interaction_type text NOT NULL,      -- e.g., 'Synergistic', 'Antagonistic'
    mechanism text,
    risk_level text NOT NULL,
    effect_description text,
    clinical_significance text,
    evidence_type text[],
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Mechanism of Action
------------------------------------------

-- Detailed mechanism data
CREATE TABLE IF NOT EXISTS mechanisms (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    mechanism_of_action text NOT NULL,
    target_name text,
    action_type text,
    binding_site text,
    pathway_involvement text[],          -- Biological pathways affected
    downstream_effects text[],           -- Downstream cellular effects
    signaling_cascades text[],          -- Affected signaling cascades
    cellular_responses text[],           -- Cellular level responses
    molecular_targets text[],            -- Additional molecular targets
    epigenetic_effects jsonb,           -- Epigenetic modifications
    gene_expression_changes jsonb,       -- Gene expression impacts
    protein_modifications jsonb,         -- Post-translational modifications
    metabolic_effects jsonb,            -- Effects on metabolism
    confidence_score double precision,
    reference_dois text[],
    evidence_type text[],               -- Types of evidence supporting mechanism
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Neurotransmitter effects
CREATE TABLE IF NOT EXISTS neurotransmitter_effects (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    neurotransmitter text NOT NULL,      -- e.g., 'Serotonin', 'Dopamine'
    effect_type text NOT NULL,           -- e.g., 'Release', 'Reuptake Inhibition'
    magnitude double precision,          -- Effect strength
    brain_regions text[],               -- Affected brain regions
    temporal_profile jsonb,             -- Time course of effects
    downstream_cascades text[],         -- Downstream signaling
    compensatory_changes text[],        -- Adaptive responses
    interaction_network jsonb,          -- Interactions with other systems
    evidence_type text[],
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Neural and Behavioral Effects
------------------------------------------

-- Synaptic effects
CREATE TABLE IF NOT EXISTS synaptic_effects (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    effect_type text NOT NULL,           -- e.g., 'LTP', 'LTD'
    synaptic_components text[],         -- Affected synaptic elements
    neurotransmitter_systems text[],    -- Involved neurotransmitters
    plasticity_changes jsonb,           -- Synaptic plasticity effects
    temporal_dynamics jsonb,            -- Time course
    spatial_extent text,                -- Local vs. distributed effects
    functional_impact text,             -- Functional consequences
    evidence_type text[],
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Neural circuit effects
CREATE TABLE IF NOT EXISTS neural_circuit_effects (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    circuit_name text NOT NULL,          -- e.g., 'Reward', 'Fear'
    effect_type text NOT NULL,
    brain_regions text[],               -- Affected regions
    connectivity_changes jsonb,         -- Network connectivity effects
    activation_patterns jsonb,          -- Activity patterns
    behavioral_correlates text[],       -- Associated behaviors
    functional_outcomes text[],         -- Functional impacts
    evidence_type text[],
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Cognitive effects
CREATE TABLE IF NOT EXISTS cognitive_effects (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    domain text NOT NULL,                -- e.g., 'Memory', 'Attention'
    effect_type text NOT NULL,
    magnitude double precision,
    onset_time interval,
    duration interval,
    dose_dependency jsonb,
    context_dependency text[],
    individual_variation jsonb,
    interaction_effects jsonb,
    evidence_type text[],
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Consciousness effects
CREATE TABLE IF NOT EXISTS consciousness_effects (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    effect_category text NOT NULL,       -- e.g., 'Psychedelic', 'Dissociative'
    subjective_effects text[],          -- Reported effects
    perceptual_changes jsonb,           -- Changes in perception
    cognitive_alterations jsonb,        -- Changes in cognition
    emotional_changes jsonb,            -- Emotional effects
    time_perception_effects jsonb,      -- Effects on time perception
    self_perception_effects jsonb,      -- Effects on self-perception
    dose_relationship jsonb,            -- Dose-dependent effects
    onset_characteristics jsonb,        -- Onset profile
    duration_profile jsonb,             -- Duration characteristics
    evidence_type text[],
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Autonomic effects
CREATE TABLE IF NOT EXISTS autonomic_effects (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    system_affected text NOT NULL,       -- e.g., 'Cardiovascular', 'Respiratory'
    effect_type text NOT NULL,
    magnitude double precision,
    temporal_profile jsonb,
    dose_dependency jsonb,
    risk_factors text[],
    compensatory_mechanisms text[],
    clinical_significance text,
    evidence_type text[],
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Neuroendocrine effects
CREATE TABLE IF NOT EXISTS neuroendocrine_effects (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    hormone_system text NOT NULL,        -- e.g., 'HPA axis', 'HPG axis'
    effect_type text NOT NULL,
    magnitude double precision,
    temporal_profile jsonb,
    feedback_mechanisms jsonb,
    systemic_impacts jsonb,
    regulatory_changes text[],
    clinical_relevance text,
    evidence_type text[],
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Behavioral effects
CREATE TABLE IF NOT EXISTS behavioral_effects (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    behavior_category text NOT NULL,     -- e.g., 'Reward', 'Anxiety'
    effect_type text NOT NULL,
    magnitude double precision,
    temporal_profile jsonb,
    context_dependency text[],
    individual_factors jsonb,
    dose_relationship jsonb,
    interaction_effects jsonb,
    evidence_type text[],
    reference_dois text[],
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Experimental Data
------------------------------------------

-- Dose-response data
CREATE TABLE IF NOT EXISTS dose_response_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    endpoint text NOT NULL,
    dose_values double precision[],
    dose_unit text,
    response_values double precision[],
    response_unit text,
    curve_parameters jsonb, -- EC50, Hill slope, etc.
    experimental_conditions jsonb,
    analysis_method text,
    reference_doi text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

-- Time course data
CREATE TABLE IF NOT EXISTS time_course_data (
    id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
    compound_id uuid NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    parameter_measured text NOT NULL,
    time_points double precision[],
    time_unit text,
    values double precision[],
    value_unit text,
    experimental_conditions jsonb,
    analysis_method text,
    reference_doi text,
    created_at timestamptz NOT NULL DEFAULT now(),
    updated_at timestamptz NOT NULL DEFAULT now()
);

------------------------------------------
-- Indexes
------------------------------------------

-- Core pharmacology indexes
CREATE INDEX IF NOT EXISTS idx_pharmacokinetic_data_compound ON pharmacokinetic_data(compound_id);
CREATE INDEX IF NOT EXISTS idx_pharmacokinetic_data_type ON pharmacokinetic_data(parameter_type);
CREATE INDEX IF NOT EXISTS idx_pharmacokinetic_data_species ON pharmacokinetic_data(species);

CREATE INDEX IF NOT EXISTS idx_adme_properties_compound ON adme_properties(compound_id);
CREATE INDEX IF NOT EXISTS idx_adme_properties_type ON adme_properties(property_type);

CREATE INDEX IF NOT EXISTS idx_pharmacodynamic_compound ON pharmacodynamic_data(compound_id);
CREATE INDEX IF NOT EXISTS idx_pharmacodynamic_type ON pharmacodynamic_data(effect_type);

-- Molecular biology indexes
CREATE INDEX IF NOT EXISTS idx_genes_symbol ON genes(symbol);
CREATE INDEX IF NOT EXISTS idx_genes_ensembl ON genes(ensembl_id);
CREATE INDEX IF NOT EXISTS idx_genes_entrez ON genes(entrez_id);
CREATE INDEX IF NOT EXISTS idx_genes_hgnc ON genes(hgnc_id);
CREATE INDEX IF NOT EXISTS idx_genes_organism ON genes(organism);

CREATE INDEX IF NOT EXISTS idx_proteins_symbol ON proteins(symbol);
CREATE INDEX IF NOT EXISTS idx_proteins_uniprot ON proteins(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_proteins_gene ON proteins(gene_id);
CREATE INDEX IF NOT EXISTS idx_proteins_organism ON proteins(organism);
CREATE INDEX IF NOT EXISTS idx_proteins_family ON proteins(protein_family);

CREATE INDEX IF NOT EXISTS idx_protein_receptor_protein ON protein_receptor_relationships(protein_id);
CREATE INDEX IF NOT EXISTS idx_protein_receptor_family ON protein_receptor_relationships(receptor_family_id);
CREATE INDEX IF NOT EXISTS idx_protein_receptor_type ON protein_receptor_relationships(relationship_type);

-- Molecular interaction indexes
CREATE INDEX IF NOT EXISTS idx_receptor_profiles_compound ON receptor_binding_profiles(compound_id);
CREATE INDEX IF NOT EXISTS idx_receptor_profiles_type ON receptor_binding_profiles(profile_type);

CREATE INDEX IF NOT EXISTS idx_enzyme_interactions_compound ON enzyme_interactions(compound_id);
CREATE INDEX IF NOT EXISTS idx_enzyme_interactions_enzyme ON enzyme_interactions(enzyme_name);
CREATE INDEX IF NOT EXISTS idx_enzyme_interactions_type ON enzyme_interactions(interaction_type);

CREATE INDEX IF NOT EXISTS idx_transporter_interactions_compound ON transporter_interactions(compound_id);
CREATE INDEX IF NOT EXISTS idx_transporter_interactions_transporter ON transporter_interactions(transporter_name);
CREATE INDEX IF NOT EXISTS idx_transporter_interactions_type ON transporter_interactions(interaction_type);

CREATE INDEX IF NOT EXISTS idx_drug_interactions_compound1 ON drug_interactions(compound_id);
CREATE INDEX IF NOT EXISTS idx_drug_interactions_compound2 ON drug_interactions(interacting_compound_id);
CREATE INDEX IF NOT EXISTS idx_drug_interactions_type ON drug_interactions(interaction_type);

-- Mechanism indexes
CREATE INDEX IF NOT EXISTS idx_mechanisms_compound ON mechanisms(compound_id);
CREATE INDEX IF NOT EXISTS idx_mechanisms_type ON mechanisms(mechanism_of_action);

CREATE INDEX IF NOT EXISTS idx_neurotransmitter_effects_compound ON neurotransmitter_effects(compound_id);
CREATE INDEX IF NOT EXISTS idx_neurotransmitter_effects_type ON neurotransmitter_effects(neurotransmitter, effect_type);

-- Neural and behavioral effect indexes
CREATE INDEX IF NOT EXISTS idx_synaptic_effects_compound ON synaptic_effects(compound_id);
CREATE INDEX IF NOT EXISTS idx_synaptic_effects_type ON synaptic_effects(effect_type);

CREATE INDEX IF NOT EXISTS idx_neural_circuit_effects_compound ON neural_circuit_effects(compound_id);
CREATE INDEX IF NOT EXISTS idx_neural_circuit_effects_circuit ON neural_circuit_effects(circuit_name);

CREATE INDEX IF NOT EXISTS idx_cognitive_effects_compound ON cognitive_effects(compound_id);
CREATE INDEX IF NOT EXISTS idx_cognitive_effects_domain ON cognitive_effects(domain);

CREATE INDEX IF NOT EXISTS idx_consciousness_effects_compound ON consciousness_effects(compound_id);
CREATE INDEX IF NOT EXISTS idx_consciousness_effects_category ON consciousness_effects(effect_category);

CREATE INDEX IF NOT EXISTS idx_autonomic_effects_compound ON autonomic_effects(compound_id);
CREATE INDEX IF NOT EXISTS idx_autonomic_effects_system ON autonomic_effects(system_affected);

CREATE INDEX IF NOT EXISTS idx_neuroendocrine_effects_compound ON neuroendocrine_effects(compound_id);
CREATE INDEX IF NOT EXISTS idx_neuroendocrine_effects_system ON neuroendocrine_effects(hormone_system);

CREATE INDEX IF NOT EXISTS idx_behavioral_effects_compound ON behavioral_effects(compound_id);
CREATE INDEX IF NOT EXISTS idx_behavioral_effects_category ON behavioral_effects(behavior_category);

-- Experimental data indexes
CREATE INDEX IF NOT EXISTS idx_dose_response_compound ON dose_response_data(compound_id);
CREATE INDEX IF NOT EXISTS idx_dose_response_endpoint ON dose_response_data(endpoint);

CREATE INDEX IF NOT EXISTS idx_time_course_compound ON time_course_data(compound_id);
CREATE INDEX IF NOT EXISTS idx_time_course_parameter ON time_course_data(parameter_measured);

-- Receptor variant and binding site indexes
CREATE INDEX IF NOT EXISTS idx_receptor_variants_subtype ON receptor_variants(subtype_id);
CREATE INDEX IF NOT EXISTS idx_receptor_variants_name ON receptor_variants(variant_name);
CREATE INDEX IF NOT EXISTS idx_receptor_variants_type ON receptor_variants(mutation_type);
CREATE INDEX IF NOT EXISTS idx_variants_phenotype ON receptor_variants USING gin (phenotype_effects);

CREATE INDEX IF NOT EXISTS idx_binding_sites_subtype ON subtype_binding_sites(subtype_id);
CREATE INDEX IF NOT EXISTS idx_binding_sites_name ON subtype_binding_sites(site_name);
CREATE INDEX IF NOT EXISTS idx_binding_sites_properties ON subtype_binding_sites USING gin (binding_properties);

CREATE INDEX IF NOT EXISTS idx_signaling_subtype ON subtype_signaling(subtype_id);
CREATE INDEX IF NOT EXISTS idx_signaling_pathway ON subtype_signaling(pathway_name);
CREATE INDEX IF NOT EXISTS idx_signaling_regulation ON subtype_signaling USING gin (regulatory_mechanisms);

CREATE INDEX IF NOT EXISTS idx_nuclear_receptor_receptor ON nuclear_receptor_data(receptor_id);
CREATE INDEX IF NOT EXISTS idx_nuclear_receptor_effects ON nuclear_receptor_data USING gin (tissue_specific_effects);

CREATE INDEX IF NOT EXISTS idx_peptide_receptor_receptor ON peptide_receptor_data(receptor_id);
CREATE INDEX IF NOT EXISTS idx_peptide_receptor_specificity ON peptide_receptor_data USING gin (peptide_specificity);

------------------------------------------
-- Triggers
------------------------------------------

-- Core pharmacology triggers
CREATE TRIGGER update_pharmacokinetic_data_modtime
    BEFORE UPDATE ON pharmacokinetic_data
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_adme_properties_modtime
    BEFORE UPDATE ON adme_properties
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_pharmacodynamic_modtime
    BEFORE UPDATE ON pharmacodynamic_data
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Molecular interaction triggers
CREATE TRIGGER update_receptor_binding_profiles_modtime
    BEFORE UPDATE ON receptor_binding_profiles
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_enzyme_interactions_modtime
    BEFORE UPDATE ON enzyme_interactions
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_transporter_interactions_modtime
    BEFORE UPDATE ON transporter_interactions
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_drug_interactions_modtime
    BEFORE UPDATE ON drug_interactions
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Mechanism triggers
CREATE TRIGGER update_mechanisms_modtime
    BEFORE UPDATE ON mechanisms
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_neurotransmitter_effects_modtime
    BEFORE UPDATE ON neurotransmitter_effects
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Neural and behavioral effect triggers
CREATE TRIGGER update_synaptic_effects_modtime
    BEFORE UPDATE ON synaptic_effects
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_neural_circuit_effects_modtime
    BEFORE UPDATE ON neural_circuit_effects
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_cognitive_effects_modtime
    BEFORE UPDATE ON cognitive_effects
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_consciousness_effects_modtime
    BEFORE UPDATE ON consciousness_effects
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_autonomic_effects_modtime
    BEFORE UPDATE ON autonomic_effects
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_neuroendocrine_effects_modtime
    BEFORE UPDATE ON neuroendocrine_effects
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_behavioral_effects_modtime
    BEFORE UPDATE ON behavioral_effects
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Experimental data triggers
CREATE TRIGGER update_dose_response_data_modtime
    BEFORE UPDATE ON dose_response_data
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_time_course_data_modtime
    BEFORE UPDATE ON time_course_data
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Receptor variant and binding site triggers
CREATE TRIGGER update_receptor_variants_modtime
    BEFORE UPDATE ON receptor_variants
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_subtype_binding_sites_modtime
    BEFORE UPDATE ON subtype_binding_sites
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_subtype_signaling_modtime
    BEFORE UPDATE ON subtype_signaling
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_nuclear_receptor_data_modtime
    BEFORE UPDATE ON nuclear_receptor_data
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_peptide_receptor_data_modtime
    BEFORE UPDATE ON peptide_receptor_data
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

------------------------------------------
-- Audit Triggers
------------------------------------------

-- Core pharmacology audit
CREATE TRIGGER audit_pharmacokinetic_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON pharmacokinetic_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_adme_properties_trigger
    AFTER INSERT OR UPDATE OR DELETE ON adme_properties
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_pharmacodynamic_trigger
    AFTER INSERT OR UPDATE OR DELETE ON pharmacodynamic_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Molecular interaction audit
CREATE TRIGGER audit_receptor_binding_profiles_trigger
    AFTER INSERT OR UPDATE OR DELETE ON receptor_binding_profiles
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_enzyme_interactions_trigger
    AFTER INSERT OR UPDATE OR DELETE ON enzyme_interactions
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_transporter_interactions_trigger
    AFTER INSERT OR UPDATE OR DELETE ON transporter_interactions
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_drug_interactions_trigger
    AFTER INSERT OR UPDATE OR DELETE ON drug_interactions
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Receptor variant and binding site audit
CREATE TRIGGER audit_receptor_variants_trigger
    AFTER INSERT OR UPDATE OR DELETE ON receptor_variants
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_subtype_binding_sites_trigger
    AFTER INSERT OR UPDATE OR DELETE ON subtype_binding_sites
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_subtype_signaling_trigger
    AFTER INSERT OR UPDATE OR DELETE ON subtype_signaling
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_nuclear_receptor_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON nuclear_receptor_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_peptide_receptor_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON peptide_receptor_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Mechanism audit
CREATE TRIGGER audit_mechanisms_trigger
    AFTER INSERT OR UPDATE OR DELETE ON mechanisms
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_neurotransmitter_effects_trigger
    AFTER INSERT OR UPDATE OR DELETE ON neurotransmitter_effects
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Neural and behavioral effect audit
CREATE TRIGGER audit_synaptic_effects_trigger
    AFTER INSERT OR UPDATE OR DELETE ON synaptic_effects
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_neural_circuit_effects_trigger
    AFTER INSERT OR UPDATE OR DELETE ON neural_circuit_effects
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_cognitive_effects_trigger
    AFTER INSERT OR UPDATE OR DELETE ON cognitive_effects
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_consciousness_effects_trigger
    AFTER INSERT OR UPDATE OR DELETE ON consciousness_effects
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_autonomic_effects_trigger
    AFTER INSERT OR UPDATE OR DELETE ON autonomic_effects
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_neuroendocrine_effects_trigger
    AFTER INSERT OR UPDATE OR DELETE ON neuroendocrine_effects
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_behavioral_effects_trigger
    AFTER INSERT OR UPDATE OR DELETE ON behavioral_effects
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

-- Experimental data audit
CREATE TRIGGER audit_dose_response_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON dose_response_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

CREATE TRIGGER audit_time_course_data_trigger
    AFTER INSERT OR UPDATE OR DELETE ON time_course_data
    FOR EACH ROW EXECUTE FUNCTION audit_trigger_func();

------------------------------------------
-- Views
------------------------------------------

-- Compound pharmacology overview
CREATE OR REPLACE VIEW compound_pharmacology_overview AS
SELECT 
    c.id as compound_id,
    c.name as compound_name,
    COUNT(DISTINCT pk.id) as pk_data_count,
    COUNT(DISTINCT pd.id) as pd_data_count,
    COUNT(DISTINCT rb.id) as receptor_binding_count,
    COUNT(DISTINCT ei.id) as enzyme_interaction_count,
    COUNT(DISTINCT ti.id) as transporter_interaction_count,
    COUNT(DISTINCT di.id) as drug_interaction_count,
    COUNT(DISTINCT m.id) as mechanism_count,
    COUNT(DISTINCT ne.id) as neurotransmitter_effect_count,
    COUNT(DISTINCT se.id) as synaptic_effect_count,
    COUNT(DISTINCT nce.id) as neural_circuit_effect_count,
    COUNT(DISTINCT ce.id) as cognitive_effect_count,
    COUNT(DISTINCT cse.id) as consciousness_effect_count,
    COUNT(DISTINCT ae.id) as autonomic_effect_count,
    COUNT(DISTINCT nee.id) as neuroendocrine_effect_count,
    COUNT(DISTINCT be.id) as behavioral_effect_count
FROM compounds c
LEFT JOIN pharmacokinetic_data pk ON c.id = pk.compound_id
LEFT JOIN pharmacodynamic_data pd ON c.id = pd.compound_id
LEFT JOIN receptor_binding_profiles rb ON c.id = rb.compound_id
LEFT JOIN enzyme_interactions ei ON c.id = ei.compound_id
LEFT JOIN transporter_interactions ti ON c.id = ti.compound_id
LEFT JOIN drug_interactions di ON c.id = di.compound_id
LEFT JOIN mechanisms m ON c.id = m.compound_id
LEFT JOIN neurotransmitter_effects ne ON c.id = ne.compound_id
LEFT JOIN synaptic_effects se ON c.id = se.compound_id
LEFT JOIN neural_circuit_effects nce ON c.id = nce.compound_id
LEFT JOIN cognitive_effects ce ON c.id = ce.compound_id
LEFT JOIN consciousness_effects cse ON c.id = cse.compound_id
LEFT JOIN autonomic_effects ae ON c.id = ae.compound_id
LEFT JOIN neuroendocrine_effects nee ON c.id = nee.compound_id
LEFT JOIN behavioral_effects be ON c.id = be.compound_id
GROUP BY c.id, c.name;

-- Mechanism summary view
CREATE OR REPLACE VIEW mechanism_summary AS
SELECT 
    c.id as compound_id,
    c.name as compound_name,
    array_agg(DISTINCT m.mechanism_of_action) as mechanisms,
    array_agg(DISTINCT m.target_name) FILTER (WHERE m.target_name IS NOT NULL) as targets,
    array_agg(DISTINCT ne.neurotransmitter) FILTER (WHERE ne.neurotransmitter IS NOT NULL) as affected_neurotransmitters,
    array_agg(DISTINCT nce.circuit_name) FILTER (WHERE nce.circuit_name IS NOT NULL) as affected_circuits,
    array_agg(DISTINCT ce.domain) FILTER (WHERE ce.domain IS NOT NULL) as cognitive_domains,
    array_agg(DISTINCT cse.effect_category) FILTER (WHERE cse.effect_category IS NOT NULL) as consciousness_effects
FROM compounds c
LEFT JOIN mechanisms m ON c.id = m.compound_id
LEFT JOIN neurotransmitter_effects ne ON c.id = ne.compound_id
LEFT JOIN neural_circuit_effects nce ON c.id = nce.compound_id
LEFT JOIN cognitive_effects ce ON c.id = ce.compound_id
LEFT JOIN consciousness_effects cse ON c.id = cse.compound_id
GROUP BY c.id, c.name;

-- Drug interaction network view
CREATE OR REPLACE VIEW drug_interaction_network AS
SELECT 
    c1.name as compound_name,
    c2.name as interacting_compound,
    di.interaction_type,
    di.mechanism,
    di.risk_level,
    di.effect_description,
    di.clinical_significance,
    di.evidence_type
FROM drug_interactions di
JOIN compounds c1 ON di.compound_id = c1.id
JOIN compounds c2 ON di.interacting_compound_id = c2.id;

-- Table comments
COMMENT ON TABLE pharmacokinetic_data IS 'Pharmacokinetic parameters and measurements';
COMMENT ON TABLE adme_properties IS 'Absorption, Distribution, Metabolism, and Excretion properties';
COMMENT ON TABLE pharmacodynamic_data IS 'Pharmacodynamic effects and characteristics';
COMMENT ON TABLE receptor_binding_profiles IS 'Detailed receptor binding and selectivity data';
COMMENT ON TABLE enzyme_interactions IS 'Interactions with metabolic enzymes';
COMMENT ON TABLE transporter_interactions IS 'Interactions with drug transporters';
COMMENT ON TABLE drug_interactions IS 'Drug-drug interaction data and mechanisms';
COMMENT ON TABLE mechanisms IS 'Detailed mechanism of action data';
COMMENT ON TABLE neurotransmitter_effects IS 'Effects on neurotransmitter systems';
COMMENT ON TABLE synaptic_effects IS 'Effects on synaptic function and plasticity';
COMMENT ON TABLE neural_circuit_effects IS 'Effects on neural circuits and networks';
COMMENT ON TABLE cognitive_effects IS 'Effects on cognitive function and domains';
COMMENT ON TABLE consciousness_effects IS 'Effects on consciousness and perception';
COMMENT ON TABLE autonomic_effects IS 'Effects on autonomic nervous system';
COMMENT ON TABLE neuroendocrine_effects IS 'Effects on neuroendocrine systems';
COMMENT ON TABLE behavioral_effects IS 'Effects on behavior and psychological function';
COMMENT ON TABLE dose_response_data IS 'Dose-response relationships and curves';
COMMENT ON TABLE time_course_data IS 'Time-dependent pharmacological measurements';

-- Receptor variant and binding site comments
COMMENT ON TABLE receptor_variants IS 'Known variants and mutations of receptor subtypes';
COMMENT ON TABLE subtype_binding_sites IS 'Binding sites specific to receptor subtypes';
COMMENT ON TABLE subtype_signaling IS 'Signaling pathways specific to receptor subtypes';
COMMENT ON TABLE nuclear_receptor_data IS 'Specialized data for nuclear receptors';
COMMENT ON TABLE peptide_receptor_data IS 'Specialized data for peptide receptors';
