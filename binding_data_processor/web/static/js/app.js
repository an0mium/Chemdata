// Initialize components when document is ready
document.addEventListener('DOMContentLoaded', () => {
    initializeStructureViewer();
    initializePlots();
    initializeSearch();
    initializeFilters();
    initializeExport();
    initializeTooltips();
    initializeEventHandlers();
});

// Structure Viewer
function initializeStructureViewer() {
    const viewer = document.getElementById('structure-viewer');
    if (!viewer) return;

    // Initialize RDKit viewer
    window.viewer = new RDKit.MolDraw2D(viewer);

    // Load initial structure if available
    const smiles = viewer.dataset.smiles;
    if (smiles) {
        window.viewer.drawMolecule(smiles);
    }

    // Add interaction handlers
    viewer.addEventListener('mousewheel', handleStructureZoom);
    viewer.addEventListener('mousedown', handleStructureRotation);
}

function handleStructureZoom(event) {
    event.preventDefault();
    const delta = event.deltaY * -0.01;
    window.viewer.zoom(delta);
}

function handleStructureRotation(event) {
    if (event.buttons !== 1) return;

    const startX = event.clientX;
    const startAngle = window.viewer.getRotation();

    function onMouseMove(e) {
        const dx = e.clientX - startX;
        const newAngle = startAngle + (dx * 0.5);
        window.viewer.setRotation(newAngle);
    }

    function onMouseUp() {
        document.removeEventListener('mousemove', onMouseMove);
        document.removeEventListener('mouseup', onMouseUp);
    }

    document.addEventListener('mousemove', onMouseMove);
    document.addEventListener('mouseup', onMouseUp);
}

// Plots
function initializePlots() {
    // Initialize binding profile plot
    const bindingPlot = document.getElementById('binding-plot');
    if (bindingPlot && bindingPlot.dataset.plot) {
        const plotData = JSON.parse(bindingPlot.dataset.plot);
        Plotly.newPlot('binding-plot', plotData.data, plotData.layout);
    }

    // Initialize activity profile plot
    const activityPlot = document.getElementById('activity-plot');
    if (activityPlot && activityPlot.dataset.plot) {
        const plotData = JSON.parse(activityPlot.dataset.plot);
        Plotly.newPlot('activity-plot', plotData.data, plotData.layout);
    }

    // Initialize safety profile plot
    const safetyPlot = document.getElementById('safety-plot');
    if (safetyPlot && safetyPlot.dataset.plot) {
        const plotData = JSON.parse(safetyPlot.dataset.plot);
        Plotly.newPlot('safety-plot', plotData.data, plotData.layout);
    }
}

// Search
function initializeSearch() {
    const searchForm = document.getElementById('searchForm');
    if (!searchForm) return;

    // Handle search submission
    searchForm.addEventListener('submit', async (event) => {
        event.preventDefault();
        const formData = new FormData(searchForm);

        try {
            showLoading('search-results');
            const response = await fetch('/api/compounds/search?' + new URLSearchParams(formData));
            const data = await response.json();

            if (data.success) {
                updateCompoundList(data.compounds);
            } else {
                showError('search-error', data.error);
            }
        } catch (error) {
            showError('search-error', 'Search failed. Please try again.');
        } finally {
            hideLoading('search-results');
        }
    });

    // Handle structure search toggle
    const structureToggle = document.getElementById('searchStructure');
    const similarityInput = document.getElementById('similarityThreshold');
    if (structureToggle && similarityInput) {
        structureToggle.addEventListener('change', () => {
            similarityInput.parentElement.style.display =
                structureToggle.checked ? 'block' : 'none';
        });
    }
}

// Filters
function initializeFilters() {
    const filterForm = document.getElementById('filterForm');
    if (!filterForm) return;

    // Handle filter submission
    filterForm.addEventListener('submit', async (event) => {
        event.preventDefault();
        const formData = new FormData(filterForm);

        try {
            showLoading('compound-list');
            const response = await fetch('/api/compounds/filter?' + new URLSearchParams(formData));
            const data = await response.json();

            if (data.success) {
                updateCompoundList(data.compounds);
            } else {
                showError('filter-error', data.error);
            }
        } catch (error) {
            showError('filter-error', 'Filtering failed. Please try again.');
        } finally {
            hideLoading('compound-list');
        }
    });

    // Handle filter reset
    filterForm.addEventListener('reset', () => {
        setTimeout(() => filterForm.dispatchEvent(new Event('submit')), 0);
    });
}

// Export
function initializeExport() {
    const exportForm = document.getElementById('exportForm');
    if (!exportForm) return;

    // Handle export submission
    exportForm.addEventListener('submit', async (event) => {
        event.preventDefault();
        const formData = new FormData(exportForm);

        try {
            showLoading('export-status');
            const response = await fetch('/api/compounds/export', {
                method: 'POST',
                body: formData
            });

            if (response.ok) {
                const blob = await response.blob();
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = `compounds.${formData.get('format')}`;
                document.body.appendChild(a);
                a.click();
                window.URL.revokeObjectURL(url);
                document.body.removeChild(a);
            } else {
                const data = await response.json();
                showError('export-error', data.error);
            }
        } catch (error) {
            showError('export-error', 'Export failed. Please try again.');
        } finally {
            hideLoading('export-status');
        }
    });
}

// Tooltips
function initializeTooltips() {
    const tooltips = document.querySelectorAll('[data-bs-toggle="tooltip"]');
    tooltips.forEach(tooltip => {
        new bootstrap.Tooltip(tooltip);
    });
}

// Event Handlers
function initializeEventHandlers() {
    // Compound list row click handler
    const compoundList = document.querySelector('.compound-list tbody');
    if (compoundList) {
        compoundList.addEventListener('click', (event) => {
            const row = event.target.closest('tr');
            if (row) {
                const casNumber = row.dataset.cas;
                window.location.href = `/compounds/${casNumber}`;
            }
        });
    }

    // Analysis section collapse handlers
    const collapseTriggers = document.querySelectorAll('[data-bs-toggle="collapse"]');
    collapseTriggers.forEach(trigger => {
        trigger.addEventListener('click', () => {
            const icon = trigger.querySelector('i');
            if (icon) {
                icon.classList.toggle('fa-chevron-down');
                icon.classList.toggle('fa-chevron-up');
            }
        });
    });
}

// Utility Functions
function showLoading(elementId) {
    const element = document.getElementById(elementId);
    if (element) {
        element.classList.add('loading');
    }
}

function hideLoading(elementId) {
    const element = document.getElementById(elementId);
    if (element) {
        element.classList.remove('loading');
    }
}

function showError(elementId, message) {
    const element = document.getElementById(elementId);
    if (element) {
        element.textContent = message;
        element.style.display = 'block';
        setTimeout(() => {
            element.style.display = 'none';
        }, 5000);
    }
}

function updateCompoundList(compounds) {
    const tbody = document.querySelector('.compound-list tbody');
    if (!tbody) return;

    tbody.innerHTML = compounds.map(compound => `
        <tr data-cas="${compound.cas_number}">
            <td class="structure-cell">
                <div class="structure-viewer" data-smiles="${compound.smiles}"></div>
            </td>
            <td>${compound.name}</td>
            <td>${compound.cas_number}</td>
            <td>${compound.psychoactive_class || ''}</td>
            <td>
                <div class="badge-container">
                    ${compound.binding_data ? `
                        <span class="badge bg-primary">
                            <i class="fas fa-dna"></i> Binding Data
                        </span>
                    ` : ''}
                    ${compound.activity_predictions ? `
                        <span class="badge bg-success">
                            <i class="fas fa-chart-line"></i> Activity
                        </span>
                    ` : ''}
                    ${compound.safety_predictions ? `
                        <span class="badge bg-warning">
                            <i class="fas fa-shield-alt"></i> Safety
                        </span>
                    ` : ''}
                </div>
            </td>
        </tr>
    `).join('');

    // Reinitialize structure viewers
    tbody.querySelectorAll('.structure-viewer').forEach(viewer => {
        const rdkitViewer = new RDKit.MolDraw2D(viewer);
        rdkitViewer.drawMolecule(viewer.dataset.smiles);
    });
}
