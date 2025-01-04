/**
 * Structure drawing and editing functionality for ChemData application
 */

const StructureEditor = {
    ketcher: null,
    templates: {},

    /**
     * Initialize structure editor
     */
    init: function () {
        // Initialize Ketcher after modal is shown
        $('#structure-modal').on('shown.bs.modal', function () {
            if (!StructureEditor.ketcher) {
                StructureEditor.ketcher = new window.ketcher.Ketcher(
                    document.getElementById('ketcher-container')
                );
                StructureEditor._loadTemplates();
            }
        });

        // Handle structure input methods
        this._initStructureInputs();

        // Handle structure tools
        this._initStructureTools();

        // Handle template management
        this._initTemplateManagement();
    },

    /**
     * Initialize structure input methods
     */
    _initStructureInputs: function () {
        // Handle SMILES input
        $('#load-smiles').click(function () {
            const smiles = $('#smiles-input').val();
            if (smiles) {
                ProgressIndicator.show('Loading structure...');
                StructureEditor.ketcher.setMolecule(smiles)
                    .then(() => {
                        ProgressIndicator.hide();
                        StructureEditor._validateStructure();
                    })
                    .catch(error => {
                        ErrorHandler.show('Invalid SMILES string: ' + error.message);
                        ProgressIndicator.hide();
                    });
            }
        });

        // Handle MOL file input
        $('#mol-file-input').change(function (e) {
            const file = e.target.files[0];
            if (file) {
                const reader = new FileReader();
                reader.onload = function (e) {
                    ProgressIndicator.show('Loading structure...');
                    StructureEditor.ketcher.setMolecule(e.target.result)
                        .then(() => {
                            ProgressIndicator.hide();
                            StructureEditor._validateStructure();
                        })
                        .catch(error => {
                            ErrorHandler.show('Invalid MOL file: ' + error.message);
                            ProgressIndicator.hide();
                        });
                };
                reader.readAsText(file);
            }
        });

        // Handle structure use
        $('#use-structure').click(function () {
            ProgressIndicator.show('Processing structure...');
            StructureEditor.ketcher.getSmiles()
                .then(smiles => {
                    window.currentStructure = smiles;
                    $('#structure-input').val(smiles);
                    $('#structure-modal').modal('hide');
                    CompoundTable.reload();
                    ProgressIndicator.hide();
                })
                .catch(error => {
                    ErrorHandler.show('Error getting structure: ' + error.message);
                    ProgressIndicator.hide();
                });
        });
    },

    /**
     * Initialize structure manipulation tools
     */
    _initStructureTools: function () {
        // Clean structure
        $('#clean-structure').click(function () {
            ProgressIndicator.show('Cleaning structure...');
            StructureEditor.ketcher.layout()
                .then(() => {
                    ProgressIndicator.hide();
                    StructureEditor._validateStructure();
                })
                .catch(error => {
                    ErrorHandler.show('Error cleaning structure: ' + error.message);
                    ProgressIndicator.hide();
                });
        });

        // Flip horizontal
        $('#flip-horizontal').click(function () {
            StructureEditor.ketcher.runCommand('flipH');
        });

        // Flip vertical
        $('#flip-vertical').click(function () {
            StructureEditor.ketcher.runCommand('flipV');
        });

        // Clear canvas
        $('#clear-canvas').click(function () {
            StructureEditor.ketcher.clear();
            $('#structure-info').hide();
            $('#structure-warnings').hide();
        });

        // Auto-clean toggle
        $('#auto-clean').change(function () {
            StructureEditor.ketcher.settings.set('autoClean', $(this).prop('checked'));
        });

        // Show implicit H toggle
        $('#show-implicit-h').change(function () {
            StructureEditor.ketcher.settings.set('showImplicitH', $(this).prop('checked'));
        });
    },

    /**
     * Initialize template management
     */
    _initTemplateManagement: function () {
        // Load saved templates
        const savedTemplates = localStorage.getItem('structureTemplates');
        if (savedTemplates) {
            this.templates = JSON.parse(savedTemplates);
            this._updateTemplateButtons();
        }

        // Handle template saving
        $('#save-template').click(function () {
            const name = $('#template-name').val();
            const smiles = $('#template-smiles').val();
            const category = $('#template-category').val();

            if (name && smiles) {
                StructureEditor.templates[name] = {
                    smiles: smiles,
                    category: category
                };
                localStorage.setItem('structureTemplates',
                    JSON.stringify(StructureEditor.templates));

                StructureEditor._updateTemplateButtons();
                $('#template-modal').modal('hide');
            }
        });

        // Handle template use
        $('#structure-templates').on('click', '.template-btn', function () {
            const smiles = $(this).data('smiles');
            StructureEditor.ketcher.setMolecule(smiles)
                .then(() => StructureEditor._validateStructure())
                .catch(error => ErrorHandler.show('Error loading template: ' + error.message));
        });
    },

    /**
     * Load default structure templates
     */
    _loadTemplates: function () {
        const defaultTemplates = {
            'Benzene': {
                smiles: 'c1ccccc1',
                category: 'rings'
            },
            'Cyclohexane': {
                smiles: 'C1CCCCC1',
                category: 'rings'
            },
            'Pyridine': {
                smiles: 'c1ccncc1',
                category: 'rings'
            },
            'Piperidine': {
                smiles: 'C1CCNCC1',
                category: 'rings'
            },
            'Indole': {
                smiles: 'c1ccc2c(c1)cc[nH]2',
                category: 'rings'
            },
            'Phenol': {
                smiles: 'Oc1ccccc1',
                category: 'functional'
            },
            'Amine': {
                smiles: 'CN',
                category: 'functional'
            },
            'Carboxylic Acid': {
                smiles: 'CC(=O)O',
                category: 'functional'
            }
        };

        // Merge with existing templates
        this.templates = { ...defaultTemplates, ...this.templates };
        this._updateTemplateButtons();
    },

    /**
     * Update template buttons display
     */
    _updateTemplateButtons: function () {
        const container = document.getElementById('structure-templates');
        container.innerHTML = '';

        // Group templates by category
        const grouped = Object.entries(this.templates)
            .reduce((acc, [name, data]) => {
                const category = data.category || 'custom';
                if (!acc[category]) acc[category] = [];
                acc[category].push({ name, ...data });
                return acc;
            }, {});

        // Create buttons for each category
        Object.entries(grouped).forEach(([category, templates]) => {
            const div = document.createElement('div');
            div.className = 'col-12';
            div.innerHTML = `<h6 class="mt-2">${category.charAt(0).toUpperCase() + category.slice(1)}</h6>`;
            container.appendChild(div);

            templates.forEach(template => {
                const col = document.createElement('div');
                col.className = 'col';
                col.innerHTML = `
                    <button type="button" class="btn btn-outline-secondary w-100 template-btn"
                            data-smiles="${template.smiles}"
                            data-bs-toggle="tooltip" title="${template.smiles}">
                        ${template.name}
                    </button>
                `;
                container.appendChild(col);
            });
        });

        // Initialize tooltips
        $('[data-bs-toggle="tooltip"]').tooltip();
    },

    /**
     * Validate current structure and update info/warnings
     */
    _validateStructure: async function () {
        try {
            const smiles = await this.ketcher.getSmiles();
            const mol = RDKitModule.get_mol(smiles);

            if (mol) {
                // Calculate properties
                const properties = ChemDataUtils.calculateDescriptors(smiles);

                // Update info display
                const infoDiv = document.getElementById('structure-properties');
                infoDiv.innerHTML = `
                    <dl class="row mb-0">
                        <dt class="col-sm-4">Molecular Weight</dt>
                        <dd class="col-sm-8">${ChemDataUtils.formatNumber(properties.mw)}</dd>
                        <dt class="col-sm-4">LogP</dt>
                        <dd class="col-sm-8">${ChemDataUtils.formatNumber(properties.logp)}</dd>
                        <dt class="col-sm-4">TPSA</dt>
                        <dd class="col-sm-8">${ChemDataUtils.formatNumber(properties.tpsa)}</dd>
                        <dt class="col-sm-4">H-Bond Donors</dt>
                        <dd class="col-sm-8">${properties.hbd}</dd>
                        <dt class="col-sm-4">H-Bond Acceptors</dt>
                        <dd class="col-sm-8">${properties.hba}</dd>
                        <dt class="col-sm-4">Rotatable Bonds</dt>
                        <dd class="col-sm-8">${properties.rotatable_bonds}</dd>
                    </dl>
                `;
                $('#structure-info').show();

                // Check for potential issues
                const warnings = [];
                if (properties.mw > 500) warnings.push('High molecular weight (> 500)');
                if (properties.logp > 5) warnings.push('High LogP (> 5)');
                if (properties.tpsa < 40) warnings.push('Low TPSA (< 40)');
                if (properties.hbd > 5) warnings.push('Many H-bond donors (> 5)');
                if (properties.hba > 10) warnings.push('Many H-bond acceptors (> 10)');

                if (warnings.length > 0) {
                    const warningList = document.getElementById('warning-list');
                    warningList.innerHTML = warnings.map(w => `<li>${w}</li>`).join('');
                    $('#structure-warnings').show();
                } else {
                    $('#structure-warnings').hide();
                }
            } else {
                $('#structure-info').hide();
                $('#structure-warnings').hide();
            }
        } catch (error) {
            console.error('Error validating structure:', error);
            $('#structure-info').hide();
            $('#structure-warnings').hide();
        }
    }
};

// Initialize structure editor when document is ready
$(document).ready(function () {
    StructureEditor.init();
});
