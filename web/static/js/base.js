/**
 * Base JavaScript functionality for ChemData application
 */

// Global state
let RDKitModule;
let currentStructure = null;
let currentCompoundData = null;

// Initialize RDKit
window.initRDKitModule().then(function (RDKit) {
    RDKitModule = RDKit;
    console.log('RDKit initialized');
}).catch(function (error) {
    console.error('Failed to initialize RDKit:', error);
});

// Utility Functions
const ChemDataUtils = {
    /**
     * Format a number with specified precision
     */
    formatNumber: function (num, precision = 2) {
        return Number.parseFloat(num).toFixed(precision);
    },

    /**
     * Convert SMILES to molblock format
     */
    smilesToMol: function (smiles) {
        try {
            const mol = RDKitModule.get_mol(smiles);
            return mol ? mol.get_molblock() : null;
        } catch (error) {
            console.error('Error converting SMILES to mol:', error);
            return null;
        }
    },

    /**
     * Calculate molecular descriptors
     */
    calculateDescriptors: function (smiles) {
        try {
            const mol = RDKitModule.get_mol(smiles);
            if (!mol) return null;

            return {
                mw: mol.get_molecular_weight(),
                logp: mol.get_logp(),
                tpsa: mol.get_tpsa(),
                hba: mol.get_num_hba(),
                hbd: mol.get_num_hbd(),
                rotatable_bonds: mol.get_num_rotatable_bonds()
            };
        } catch (error) {
            console.error('Error calculating descriptors:', error);
            return null;
        }
    },

    /**
     * Calculate fingerprint similarity between two SMILES
     */
    calculateSimilarity: function (smiles1, smiles2) {
        try {
            const mol1 = RDKitModule.get_mol(smiles1);
            const mol2 = RDKitModule.get_mol(smiles2);
            if (!mol1 || !mol2) return 0;

            const fp1 = mol1.get_morgan_fp(2);
            const fp2 = mol2.get_morgan_fp(2);
            return RDKitModule.tanimoto_similarity(fp1, fp2);
        } catch (error) {
            console.error('Error calculating similarity:', error);
            return 0;
        }
    }
};

// Progress Indicator
const ProgressIndicator = {
    show: function (message = 'Processing...') {
        $('.progress-text').text(message);
        $('.toast').toast('show');
        $('.progress-bar').css('width', '100%');
    },

    hide: function () {
        $('.toast').toast('hide');
        $('.progress-bar').css('width', '0%');
    },

    update: function (progress, message) {
        $('.progress-bar').css('width', `${progress}%`);
        if (message) {
            $('.progress-text').text(message);
        }
    }
};

// Error Handling
const ErrorHandler = {
    show: function (message, type = 'error') {
        const alertClass = type === 'warning' ? 'alert-warning' : 'alert-danger';
        const alert = $(`
            <div class="alert ${alertClass} alert-dismissible fade show" role="alert">
                ${message}
                <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
            </div>
        `);

        // Remove existing alerts
        $('.alert').alert('close');

        // Add new alert at the top of the content
        $('.container-fluid').prepend(alert);

        // Auto-dismiss after 5 seconds
        setTimeout(() => {
            alert.alert('close');
        }, 5000);
    },

    handleAjaxError: function (error) {
        let message = 'An error occurred while processing your request.';

        if (error.responseJSON && error.responseJSON.error) {
            message = error.responseJSON.error;
        } else if (error.statusText) {
            message = error.statusText;
        }

        this.show(message);
    }
};

// Event Handlers
$(document).ready(function () {
    // Initialize tooltips
    $('[data-bs-toggle="tooltip"]').tooltip();

    // Initialize popovers
    $('[data-bs-toggle="popover"]').popover();

    // Handle AJAX errors globally
    $(document).ajaxError(function (event, jqXHR, settings, error) {
        ErrorHandler.handleAjaxError(jqXHR);
    });

    // Handle structure copy
    $('.copy-structure').click(function () {
        const smiles = $(this).data('smiles');
        navigator.clipboard.writeText(smiles).then(function () {
            // Show success message
            $(this).tooltip('hide')
                .attr('data-bs-original-title', 'Copied!')
                .tooltip('show');

            // Reset tooltip after 1 second
            setTimeout(() => {
                $(this).attr('data-bs-original-title', 'Copy SMILES');
            }, 1000);
        }).catch(function (err) {
            ErrorHandler.show('Failed to copy structure: ' + err);
        });
    });

    // Handle dark mode toggle
    $('#dark-mode-toggle').change(function () {
        $('body').toggleClass('dark-mode');
        localStorage.setItem('darkMode', $('body').hasClass('dark-mode'));
    });

    // Initialize dark mode from localStorage
    if (localStorage.getItem('darkMode') === 'true') {
        $('body').addClass('dark-mode');
        $('#dark-mode-toggle').prop('checked', true);
    }
});

// Export utilities for use in other modules
window.ChemDataUtils = ChemDataUtils;
window.ProgressIndicator = ProgressIndicator;
window.ErrorHandler = ErrorHandler;
