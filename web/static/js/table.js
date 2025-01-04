/**
 * Table initialization and handling for ChemData application
 */

const CompoundTable = {
    table: null,

    /**
     * Initialize the compounds DataTable
     */
    init: function () {
        this.table = $('#compounds-table').DataTable({
            serverSide: true,
            ajax: {
                url: '/api/compounds',
                data: function (d) {
                    return {
                        search: $('#search').val(),
                        target: $('#target-filter').val(),
                        activity: $('#activity-filter').val(),
                        legal: $('#legal-filter').val(),
                        structure: window.currentStructure || '',
                        similarity: $('#similarity-threshold').val() / 100,
                        high_toxicity: $('#high-toxicity').prop('checked'),
                        high_abuse: $('#high-abuse').prop('checked'),
                        columns: $('.column-toggle:checked').map(function () {
                            return $(this).val();
                        }).get(),
                        page: Math.floor(d.start / d.length) + 1,
                        per_page: d.length,
                        draw: d.draw
                    };
                }
            },
            columns: [
                {
                    data: 'smiles',
                    render: function (data) {
                        return data ?
                            `<img src="/api/structure/${encodeURIComponent(data)}" class="structure-img">` :
                            '';
                    }
                },
                {
                    data: 'name',
                    render: function (data, type, row) {
                        if (type === 'display') {
                            return `<a href="#" class="compound-detail" data-compound='${JSON.stringify(row)}'>${data}</a>`;
                        }
                        return data;
                    }
                },
                {
                    data: 'target_type',
                    render: function (data) {
                        return data ? `<span class="badge bg-info">${data}</span>` : '';
                    }
                },
                {
                    data: 'binding_summary',
                    render: function (data) {
                        if (data) {
                            let html = '<div class="activity-data">';
                            if (data.strongest_binding) {
                                html += `<strong>${data.strongest_binding.target}</strong>: `;
                                html += `${data.strongest_binding.value} ${data.strongest_binding.unit}<br>`;
                            }
                            if (data.primary_targets && data.primary_targets.length > 0) {
                                html += `<small>Primary targets: ${data.primary_targets.join(', ')}</small>`;
                            }
                            html += '</div>';
                            return html;
                        }
                        return '';
                    }
                },
                {
                    data: null,
                    render: function (data, type, row) {
                        if (type === 'display') {
                            const props = [];
                            if (row.molecular_weight) props.push(`MW: ${ChemDataUtils.formatNumber(row.molecular_weight)}`);
                            if (row.logp) props.push(`LogP: ${ChemDataUtils.formatNumber(row.logp)}`);
                            if (row.tpsa) props.push(`TPSA: ${ChemDataUtils.formatNumber(row.tpsa)}`);
                            return props.join('<br>');
                        }
                        return '';
                    }
                },
                {
                    data: null,
                    render: function (data, type, row) {
                        if (type === 'display') {
                            let html = '<div class="predictions-summary">';

                            // Add toxicity predictions
                            if (row.toxicity_predictions) {
                                const highRisks = Object.entries(row.toxicity_predictions)
                                    .filter(([_, pred]) => pred.probability > 0.7);
                                if (highRisks.length > 0) {
                                    html += `<div class="alert alert-danger py-1 px-2 mb-1">`;
                                    html += `High toxicity risk (${highRisks.length})`;
                                    html += `</div>`;
                                }
                            }

                            // Add abuse potential
                            if (row.abuse_potential) {
                                const highRisks = Object.entries(row.abuse_potential)
                                    .filter(([_, pred]) => pred.probability > 0.7);
                                if (highRisks.length > 0) {
                                    html += `<div class="alert alert-warning py-1 px-2 mb-1">`;
                                    html += `High abuse potential (${highRisks.length})`;
                                    html += `</div>`;
                                }
                            }

                            html += '</div>';
                            return html;
                        }
                        return '';
                    }
                },
                {
                    data: 'regulatory_status',
                    render: function (data) {
                        if (data) {
                            let html = '';
                            if (data.scheduling) {
                                html += `<span class="badge bg-danger">${data.scheduling}</span><br>`;
                            }
                            if (data.warnings && data.warnings.length > 0) {
                                html += `<small class="text-danger">${data.warnings.join(', ')}</small>`;
                            }
                            return html;
                        }
                        return '';
                    }
                },
                {
                    data: 'urls',
                    render: function (data) {
                        if (data) {
                            var links = [];
                            for (var key in data) {
                                var name = key.replace('_url', '').toUpperCase();
                                links.push(
                                    `<li><a href="${data[key]}" target="_blank">${name}</a></li>`
                                );
                            }
                            return `<ul class="url-list">${links.join('')}</ul>`;
                        }
                        return '';
                    }
                }
            ],
            pageLength: 25,
            order: [[1, 'asc']],
            drawCallback: function () {
                // Initialize tooltips on newly loaded rows
                $('[data-bs-toggle="tooltip"]').tooltip();
            }
        });

        // Handle compound detail click
        $('#compounds-table').on('click', '.compound-detail', function (e) {
            e.preventDefault();
            const data = JSON.parse($(this).attr('data-compound'));
            CompoundDetail.show(data);
        });

        // Handle filter form submission
        $('#filter-form').on('submit', function (e) {
            e.preventDefault();
            CompoundTable.reload();
        });

        // Handle column visibility changes
        $('.column-toggle').on('change', function () {
            CompoundTable.reload();
        });

        // Handle similarity threshold changes
        $('#similarity-threshold').on('input', function () {
            $('#similarity-value').text($(this).val() + '%');
        });

        // Progress handling
        this.table.on('processing.dt', function (e, settings, processing) {
            if (processing) {
                ProgressIndicator.show('Loading compounds...');
            } else {
                ProgressIndicator.hide();
            }
        });

        // Error handling
        this.table.on('error.dt', function (e, settings, techNote, message) {
            ErrorHandler.show('Error loading compounds: ' + message);
        });
    },

    /**
     * Reload the table with current filters
     */
    reload: function () {
        this.table.ajax.reload();
    },

    /**
     * Export current view to various formats
     */
    export: {
        toCSV: function () {
            window.location.href = '/api/compounds/export/csv?' + $.param(this._getExportParams());
        },

        toJSON: function () {
            window.location.href = '/api/compounds/export/json?' + $.param(this._getExportParams());
        },

        toSDF: function () {
            window.location.href = '/api/compounds/export/sdf?' + $.param(this._getExportParams());
        },

        _getExportParams: function () {
            return {
                search: $('#search').val(),
                target: $('#target-filter').val(),
                activity: $('#activity-filter').val(),
                legal: $('#legal-filter').val(),
                structure: window.currentStructure || '',
                similarity: $('#similarity-threshold').val() / 100,
                high_toxicity: $('#high-toxicity').prop('checked'),
                high_abuse: $('#high-abuse').prop('checked')
            };
        }
    }
};

// Initialize table when document is ready
$(document).ready(function () {
    CompoundTable.init();
});
