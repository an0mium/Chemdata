/**
 * Filter handling functionality for ChemData application
 */

const FilterManager = {
    /**
     * Initialize filter handlers
     */
    init: function () {
        this._initTextSearch();
        this._initSelectFilters();
        this._initSimilarityFilter();
        this._initMLFilters();
        this._initColumnVisibility();
        this._initFilterForm();
        this._initFilterStorage();
    },

    /**
     * Initialize text search functionality
     */
    _initTextSearch: function () {
        let searchTimeout;
        $('#search').on('input', function () {
            clearTimeout(searchTimeout);
            searchTimeout = setTimeout(() => {
                CompoundTable.reload();
            }, 500); // Debounce search input
        });
    },

    /**
     * Initialize select filters (target, activity, legal status)
     */
    _initSelectFilters: function () {
        // Update options when data changes
        this._updateFilterOptions('target-filter', '/api/targets');
        this._updateFilterOptions('activity-filter', '/api/activities');
        this._updateFilterOptions('legal-filter', '/api/legal-statuses');

        // Handle changes
        $('#target-filter, #activity-filter, #legal-filter').on('change', function () {
            CompoundTable.reload();
        });
    },

    /**
     * Initialize similarity search functionality
     */
    _initSimilarityFilter: function () {
        // Update threshold display
        $('#similarity-threshold').on('input', function () {
            $('#similarity-value').text($(this).val() + '%');
        });

        // Handle structure input changes
        $('#structure-input').on('change', function () {
            if ($(this).val()) {
                $('#similarity-controls').show();
                CompoundTable.reload();
            } else {
                $('#similarity-controls').hide();
            }
        });

        // Handle threshold changes
        $('#similarity-threshold').on('change', function () {
            if ($('#structure-input').val()) {
                CompoundTable.reload();
            }
        });
    },

    /**
     * Initialize ML-based filters
     */
    _initMLFilters: function () {
        // Handle toxicity filter changes
        $('#high-toxicity').on('change', function () {
            CompoundTable.reload();
        });

        // Handle abuse potential filter changes
        $('#high-abuse').on('change', function () {
            CompoundTable.reload();
        });

        // Initialize tooltips for ML filters
        const mlTooltips = {
            'high-toxicity': 'Show compounds with >70% predicted toxicity risk',
            'high-abuse': 'Show compounds with >70% predicted abuse potential'
        };

        Object.entries(mlTooltips).forEach(([id, text]) => {
            $(`#${id}`).parent().attr('title', text);
        });
        $('[title]').tooltip();
    },

    /**
     * Initialize column visibility controls
     */
    _initColumnVisibility: function () {
        // Group column toggles
        const groups = {};
        $('.column-toggle').each(function () {
            const group = $(this).closest('.column-group').find('label').text();
            if (!groups[group]) groups[group] = [];
            groups[group].push(this);
        });

        // Add group toggle buttons
        Object.entries(groups).forEach(([group, toggles]) => {
            const container = $(toggles[0]).closest('.column-group');
            const btn = $(`
                <button type="button" class="btn btn-sm btn-outline-secondary group-toggle mb-2">
                    Toggle ${group}
                </button>
            `);
            container.prepend(btn);

            // Handle group toggle
            btn.click(function () {
                const anyUnchecked = toggles.some(t => !t.checked);
                toggles.forEach(t => t.checked = anyUnchecked);
                CompoundTable.reload();
            });
        });

        // Handle individual column toggles
        $('.column-toggle').on('change', function () {
            CompoundTable.reload();
        });
    },

    /**
     * Initialize filter form submission
     */
    _initFilterForm: function () {
        $('#filter-form').on('submit', function (e) {
            e.preventDefault();
            CompoundTable.reload();
        });

        // Add reset button
        const resetBtn = $(`
            <button type="button" class="btn btn-outline-secondary me-2">
                Reset Filters
            </button>
        `);
        $('#filter-form button[type="submit"]').before(resetBtn);

        // Handle reset
        resetBtn.click(() => {
            this._resetFilters();
        });
    },

    /**
     * Initialize filter state storage
     */
    _initFilterStorage: function () {
        // Restore filters on page load
        const savedFilters = localStorage.getItem('compoundFilters');
        if (savedFilters) {
            const filters = JSON.parse(savedFilters);
            this._applyFilters(filters);
        }

        // Save filters when changed
        const saveFilters = () => {
            const filters = this._getCurrentFilters();
            localStorage.setItem('compoundFilters', JSON.stringify(filters));
        };

        // Add save handlers to all filter controls
        $('#search, #target-filter, #activity-filter, #legal-filter, ' +
            '#structure-input, #similarity-threshold, #high-toxicity, #high-abuse')
            .on('change', saveFilters);
    },

    /**
     * Update filter select options from API
     */
    _updateFilterOptions: async function (selectId, url) {
        try {
            const response = await fetch(url);
            if (!response.ok) throw new Error('Failed to fetch options');

            const options = await response.json();
            const select = document.getElementById(selectId);
            const currentValue = select.value;

            // Clear existing options (except first "All" option)
            while (select.options.length > 1) {
                select.remove(1);
            }

            // Add new options
            options.forEach(option => {
                const opt = document.createElement('option');
                opt.value = option.value || option;
                opt.textContent = option.label || option;
                select.appendChild(opt);
            });

            // Restore previous selection if it exists in new options
            if (currentValue && options.some(o => (o.value || o) === currentValue)) {
                select.value = currentValue;
            }
        } catch (error) {
            console.error(`Error updating ${selectId} options:`, error);
            ErrorHandler.show(`Failed to load ${selectId.replace('-', ' ')} options`);
        }
    },

    /**
     * Get current filter values
     */
    _getCurrentFilters: function () {
        return {
            search: $('#search').val(),
            target: $('#target-filter').val(),
            activity: $('#activity-filter').val(),
            legal: $('#legal-filter').val(),
            structure: $('#structure-input').val(),
            similarity: $('#similarity-threshold').val(),
            highToxicity: $('#high-toxicity').prop('checked'),
            highAbuse: $('#high-abuse').prop('checked'),
            columns: $('.column-toggle:checked').map(function () {
                return $(this).val();
            }).get()
        };
    },

    /**
     * Apply saved filters
     */
    _applyFilters: function (filters) {
        $('#search').val(filters.search);
        $('#target-filter').val(filters.target);
        $('#activity-filter').val(filters.activity);
        $('#legal-filter').val(filters.legal);
        $('#structure-input').val(filters.structure);
        $('#similarity-threshold').val(filters.similarity);
        $('#high-toxicity').prop('checked', filters.highToxicity);
        $('#high-abuse').prop('checked', filters.highAbuse);

        // Update column visibility
        $('.column-toggle').each(function () {
            $(this).prop('checked', filters.columns.includes($(this).val()));
        });

        // Update similarity threshold display
        if (filters.similarity) {
            $('#similarity-value').text(filters.similarity + '%');
        }

        // Show/hide similarity controls
        if (filters.structure) {
            $('#similarity-controls').show();
        }

        // Reload table with restored filters
        CompoundTable.reload();
    },

    /**
     * Reset all filters to default values
     */
    _resetFilters: function () {
        // Clear text inputs
        $('#search, #structure-input').val('');

        // Reset selects to first option
        $('#target-filter, #activity-filter, #legal-filter').prop('selectedIndex', 0);

        // Reset similarity threshold
        $('#similarity-threshold').val(70);
        $('#similarity-value').text('70%');
        $('#similarity-controls').hide();

        // Uncheck ML filters
        $('#high-toxicity, #high-abuse').prop('checked', false);

        // Reset column visibility
        $('.column-toggle').prop('checked', true);

        // Clear saved filters
        localStorage.removeItem('compoundFilters');

        // Reload table
        CompoundTable.reload();
    }
};

// Initialize filters when document is ready
$(document).ready(function () {
    FilterManager.init();
});
