/**
 * Compound detail view handling for ChemData application
 */

const CompoundDetail = {
    currentData: null,
    charts: {},

    /**
     * Show compound detail modal with data
     */
    show: function (data) {
        this.currentData = data;

        // Update basic info
        this._updateBasicInfo();

        // Update structure display
        this._updateStructureDisplay();

        // Update predictions
        this._updatePredictions();

        // Update activity data
        this._updateActivityData();

        // Update community data
        this._updateCommunityData();

        // Update literature data
        this._updateLiteratureData();

        // Update similar compounds
        this._updateSimilarCompounds();

        // Show modal
        $('#compound-detail-modal').modal('show');
    },

    /**
     * Update basic compound information
     */
    _updateBasicInfo: function () {
        const data = this.currentData;

        $('#detail-name').text(data.name);
        $('#detail-smiles').text(data.smiles);
        $('#detail-mw').text(ChemDataUtils.formatNumber(data.molecular_weight));
        $('#detail-logp').text(ChemDataUtils.formatNumber(data.logp));
        $('#detail-tpsa').text(ChemDataUtils.formatNumber(data.tpsa));
        $('#detail-hbd').text(data.hbd);
        $('#detail-hba').text(data.hba);
        $('#detail-rotatable').text(data.rotatable_bonds);

        // Update copy button data
        $('.copy-structure').data('smiles', data.smiles);
    },

    /**
     * Update structure display
     */
    _updateStructureDisplay: function () {
        const img = new Image();
        img.src = `/api/structure/${encodeURIComponent(this.currentData.smiles)}`;
        img.onload = function () {
            const container = document.getElementById('structure-display');
            container.innerHTML = '';
            container.appendChild(img);
        };
    },

    /**
     * Update ML predictions
     */
    _updatePredictions: function () {
        this._updateToxicityPredictions();
        this._updateAbusePredictions();
        this._updateBindingPredictions();
    },

    _updateToxicityPredictions: function () {
        const container = document.getElementById('toxicity-predictions');
        container.innerHTML = '';

        const predictions = this.currentData.toxicity_predictions;
        if (!predictions) {
            container.innerHTML = '<p>No toxicity predictions available</p>';
            return;
        }

        Object.entries(predictions)
            .sort((a, b) => b[1].probability - a[1].probability)
            .forEach(([type, pred]) => {
                const div = document.createElement('div');
                div.className = 'mb-2';

                const probability = pred.probability * 100;
                const badgeClass = probability > 70 ? 'bg-danger' :
                    probability > 30 ? 'bg-warning' : 'bg-success';

                div.innerHTML = `
                    <span class="prediction-badge ${badgeClass}">
                        ${type}: ${probability.toFixed(1)}%
                    </span>
                    <div class="confidence-bar">
                        <div class="bar" style="width: ${pred.confidence * 100}%"></div>
                    </div>
                    ${pred.explanation ? `<small class="text-muted">${pred.explanation}</small>` : ''}
                `;

                container.appendChild(div);
            });
    },

    _updateAbusePredictions: function () {
        const container = document.getElementById('abuse-predictions');
        container.innerHTML = '';

        const predictions = this.currentData.abuse_potential;
        if (!predictions) {
            container.innerHTML = '<p>No abuse potential predictions available</p>';
            return;
        }

        Object.entries(predictions)
            .sort((a, b) => b[1].probability - a[1].probability)
            .forEach(([type, pred]) => {
                const div = document.createElement('div');
                div.className = 'mb-2';

                const probability = pred.probability * 100;
                const badgeClass = probability > 70 ? 'bg-danger' :
                    probability > 30 ? 'bg-warning' : 'bg-success';

                div.innerHTML = `
                    <span class="prediction-badge ${badgeClass}">
                        ${type}: ${probability.toFixed(1)}%
                    </span>
                    <div class="confidence-bar">
                        <div class="bar" style="width: ${pred.confidence * 100}%"></div>
                    </div>
                `;

                if (pred.mechanisms && pred.mechanisms.length > 0) {
                    div.innerHTML += `
                        <ul class="mechanism-list">
                            ${pred.mechanisms.map(m => `<li>${m}</li>`).join('')}
                        </ul>
                    `;
                }

                container.appendChild(div);
            });
    },

    _updateBindingPredictions: function () {
        const container = document.getElementById('binding-predictions');
        container.innerHTML = '';

        const predictions = this.currentData.binding_predictions;
        if (!predictions || predictions.length === 0) {
            container.innerHTML = '<p>No binding predictions available</p>';
            return;
        }

        predictions
            .sort((a, b) => b.probability - a.probability)
            .slice(0, 5) // Top 5 predictions
            .forEach(pred => {
                const div = document.createElement('div');
                div.className = 'mb-2';

                const probability = pred.probability * 100;
                const badgeClass = probability > 70 ? 'bg-primary' :
                    probability > 30 ? 'bg-info' : 'bg-secondary';

                div.innerHTML = `
                    <span class="prediction-badge ${badgeClass}">
                        ${pred.target}: ${probability.toFixed(1)}%
                    </span>
                    <div class="confidence-bar">
                        <div class="bar" style="width: ${pred.confidence * 100}%"></div>
                    </div>
                `;

                container.appendChild(div);
            });
    },

    /**
     * Update activity data
     */
    _updateActivityData: function () {
        this._updateBindingSummary();
        this._updateActivityChart();
        this._updateActivityTable();
    },

    _updateBindingSummary: function () {
        const container = document.getElementById('binding-summary');
        const summary = this.currentData.binding_summary;

        if (!summary) {
            container.innerHTML = '<p>No binding data available</p>';
            return;
        }

        container.innerHTML = `
            <div class="mb-3">
                <h6>Primary Targets:</h6>
                ${summary.primary_targets.map(target =>
            `<span class="badge bg-primary me-1">${target}</span>`
        ).join('')}
            </div>
            ${summary.strongest_binding ? `
                <div class="mb-3">
                    <h6>Strongest Binding:</h6>
                    <p class="mb-0">
                        ${summary.strongest_binding.target}: 
                        ${summary.strongest_binding.value} ${summary.strongest_binding.unit}
                    </p>
                </div>
            ` : ''}
        `;
    },

    _updateActivityChart: function () {
        const data = this.currentData.binding_summary;
        if (!data || !data.target_families) return;

        const chartData = [{
            type: 'bar',
            x: data.target_families,
            y: data.target_families.map(f =>
                data.binding_data.filter(b => b.target_family === f).length
            ),
            marker: {
                color: 'rgb(55, 83, 109)'
            }
        }];

        const layout = {
            title: 'Target Family Distribution',
            xaxis: { title: 'Target Family' },
            yaxis: { title: 'Number of Interactions' },
            margin: { l: 50, r: 50, t: 50, b: 100 },
            template: 'plotly_white'
        };

        Plotly.newPlot('activity-chart', chartData, layout, {
            displayModeBar: false,
            responsive: true
        });
    },

    _updateActivityTable: function () {
        const tbody = document.getElementById('activity-data');
        tbody.innerHTML = '';

        const activities = this.currentData.binding_summary?.binding_data || [];

        activities.forEach(activity => {
            const row = document.createElement('tr');
            row.innerHTML = `
                <td>${activity.target}</td>
                <td>${activity.type}</td>
                <td>${ChemDataUtils.formatNumber(activity.value)}</td>
                <td>${activity.unit}</td>
                <td>
                    ${activity.reference ?
                    `<a href="${activity.reference}" target="_blank">
                            <i class="bi bi-link-45deg"></i>
                        </a>` :
                    ''}
                </td>
            `;
            tbody.appendChild(row);
        });
    },

    /**
     * Update community data
     */
    _updateCommunityData: function () {
        const container = document.getElementById('community-data');
        const summary = this.currentData.community_summary;

        if (!summary) {
            container.innerHTML = '<p>No community data available</p>';
            return;
        }

        container.innerHTML = `
            <p><strong>Total Reports:</strong> ${summary.total_reports}</p>
            <h6>Reported Effects:</h6>
            <div class="mb-3">
                ${summary.reported_effects.map(effect =>
            `<span class="badge bg-info me-1">${effect}</span>`
        ).join('')}
            </div>
            ${summary.safety_concerns.length > 0 ? `
                <div class="alert alert-warning">
                    <h6>Safety Concerns:</h6>
                    <ul class="mb-0">
                        ${summary.safety_concerns.map(concern =>
            `<li>${concern}</li>`
        ).join('')}
                    </ul>
                </div>
            ` : ''}
        `;

        // Update experience reports
        const reportsContainer = document.getElementById('experience-reports');
        if (summary.experience_reports && summary.experience_reports.length > 0) {
            reportsContainer.innerHTML = `
                <h6>Recent Reports:</h6>
                <div class="list-group">
                    ${summary.experience_reports.slice(0, 3).map(report => `
                        <div class="list-group-item">
                            <div class="d-flex justify-content-between">
                                <small class="text-muted">${report.date}</small>
                                <small class="text-muted">${report.source}</small>
                            </div>
                            <p class="mb-1">${report.summary}</p>
                        </div>
                    `).join('')}
                </div>
            `;
        } else {
            reportsContainer.innerHTML = '';
        }
    },

    /**
     * Update literature data
     */
    _updateLiteratureData: function () {
        const container = document.getElementById('literature-data');
        const data = this.currentData;

        container.innerHTML = `
            <div class="mb-3">
                <h6>References:</h6>
                <ul class="list-unstyled">
                    ${data.reference_dois.map(doi =>
            `<li><a href="https://doi.org/${doi}" target="_blank">${doi}</a></li>`
        ).join('')}
                </ul>
            </div>
        `;

        // Update patent data
        const patentContainer = document.getElementById('patent-data');
        if (data.patent_count > 0) {
            patentContainer.innerHTML = `
                <h6>Patents (${data.patent_count}):</h6>
                <ul class="list-unstyled">
                    ${data.patent_numbers.map((num, i) =>
                `<li>${num}: ${data.patent_titles[i] || ''}</li>`
            ).join('')}
                </ul>
            `;
        } else {
            patentContainer.innerHTML = '';
        }
    },

    /**
     * Update similar compounds
     */
    _updateSimilarCompounds: function () {
        const container = document.getElementById('similar-compounds');
        const similar = this.currentData.similar_compounds;

        if (!similar || similar.length === 0) {
            container.innerHTML = '<p>No similar compounds found</p>';
            return;
        }

        container.innerHTML = similar.map(compound => `
            <div class="col-md-4 col-lg-3">
                <div class="card mb-3">
                    <img src="/api/structure/${encodeURIComponent(compound.smiles)}" 
                         class="card-img-top p-2" alt="${compound.name}">
                    <div class="card-body">
                        <h6 class="card-title">${compound.name}</h6>
                        <p class="card-text">
                            Similarity: ${(compound.similarity * 100).toFixed(1)}%
                        </p>
                        <button class="btn btn-sm btn-outline-primary compound-detail"
                                data-compound='${JSON.stringify(compound)}'>
                            View Details
                        </button>
                    </div>
                </div>
            </div>
        `).join('');
    }
};

// Initialize event handlers
$(document).ready(function () {
    // Handle export buttons
    $('#export-json').click(function (e) {
        e.preventDefault();
        if (CompoundDetail.currentData) {
            const blob = new Blob(
                [JSON.stringify(CompoundDetail.currentData, null, 2)],
                { type: 'application/json' }
            );
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = `${CompoundDetail.currentData.name}.json`;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
        }
    });

    $('#export-sdf').click(function (e) {
        e.preventDefault();
        if (CompoundDetail.currentData) {
            window.location.href = `/api/compounds/${CompoundDetail.currentData.id}/sdf`;
        }
    });

    $('#export-csv').click(function (e) {
        e.preventDefault();
        if (CompoundDetail.currentData) {
            window.location.href = `/api/compounds/${CompoundDetail.currentData.id}/csv`;
        }
    });

    $('#export-report').click(function (e) {
        e.preventDefault();
        if (CompoundDetail.currentData) {
            window.location.href = `/api/compounds/${CompoundDetail.currentData.id}/report`;
        }
    });

    // Handle similar compound clicks
    $('#similar-compounds').on('click', '.compound-detail', function (e) {
        e.preventDefault();
        const data = JSON.parse($(this).data('compound'));
        CompoundDetail.show(data);
    });

    // Clean up charts when modal is hidden
    $('#compound-detail-modal').on('hidden.bs.modal', function () {
        Object.values(CompoundDetail.charts).forEach(chart => {
            if (chart && typeof chart.destroy === 'function') {
                chart.destroy();
            }
        });
        CompoundDetail.charts = {};
    });
});
