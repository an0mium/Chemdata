/**
 * @jest-environment jsdom
 */

import '@testing-library/jest-dom';
import { fireEvent, waitFor } from '@testing-library/dom';

// Mock external dependencies
global.RDKit = {
    MolDraw2D: jest.fn().mockImplementation(() => ({
        drawMolecule: jest.fn(),
        zoom: jest.fn(),
        getRotation: jest.fn().mockReturnValue(0),
        setRotation: jest.fn()
    }))
};

global.Plotly = {
    newPlot: jest.fn()
};

global.bootstrap = {
    Tooltip: jest.fn()
};

// Mock fetch
global.fetch = jest.fn();

// Helper to setup document body
function setupDOM(html) {
    document.body.innerHTML = html;
}

// Test data
const testCompounds = [
    {
        name: 'Caffeine',
        smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
        cas_number: '58-08-2',
        psychoactive_class: 'STIMULANT',
        binding_data: true,
        activity_predictions: true,
        safety_predictions: true
    },
    {
        name: 'Amphetamine',
        smiles: 'CC(N)CC1=CC=CC=C1',
        cas_number: '300-62-9',
        psychoactive_class: 'STIMULANT',
        binding_data: true,
        activity_predictions: true,
        safety_predictions: true
    }
];

// Tests
describe('Structure Viewer', () => {
    beforeEach(() => {
        setupDOM(`
            <div id="structure-viewer" data-smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C"></div>
        `);
        require('../app.js');
    });

    test('initializes RDKit viewer', () => {
        expect(RDKit.MolDraw2D).toHaveBeenCalled();
    });

    test('handles zoom', () => {
        const viewer = document.getElementById('structure-viewer');
        fireEvent.wheel(viewer, { deltaY: -100 });
        expect(global.RDKit.MolDraw2D.mock.results[0].value.zoom).toHaveBeenCalled();
    });

    test('handles rotation', async () => {
        const viewer = document.getElementById('structure-viewer');
        fireEvent.mouseDown(viewer, { buttons: 1, clientX: 0 });
        fireEvent.mouseMove(document, { clientX: 50 });
        fireEvent.mouseUp(document);
        expect(global.RDKit.MolDraw2D.mock.results[0].value.setRotation).toHaveBeenCalled();
    });
});

describe('Plots', () => {
    beforeEach(() => {
        setupDOM(`
            <div id="binding-plot" data-plot='{"data":[],"layout":{}}'></div>
            <div id="activity-plot" data-plot='{"data":[],"layout":{}}'></div>
            <div id="safety-plot" data-plot='{"data":[],"layout":{}}'></div>
        `);
        require('../app.js');
    });

    test('initializes all plots', () => {
        expect(Plotly.newPlot).toHaveBeenCalledTimes(3);
    });
});

describe('Search', () => {
    beforeEach(() => {
        setupDOM(`
            <form id="searchForm">
                <input type="text" name="query" value="test">
                <button type="submit">Search</button>
            </form>
            <div id="search-results"></div>
            <div id="search-error"></div>
        `);
        require('../app.js');
    });

    test('handles successful search', async () => {
        fetch.mockImplementationOnce(() =>
            Promise.resolve({
                ok: true,
                json: () => Promise.resolve({
                    success: true,
                    compounds: testCompounds
                })
            })
        );

        const form = document.getElementById('searchForm');
        fireEvent.submit(form);

        await waitFor(() => {
            expect(fetch).toHaveBeenCalled();
        });
    });

    test('handles search error', async () => {
        fetch.mockImplementationOnce(() =>
            Promise.resolve({
                ok: true,
                json: () => Promise.resolve({
                    success: false,
                    error: 'Search failed'
                })
            })
        );

        const form = document.getElementById('searchForm');
        fireEvent.submit(form);

        await waitFor(() => {
            const error = document.getElementById('search-error');
            expect(error.textContent).toBe('Search failed');
        });
    });
});

describe('Filters', () => {
    beforeEach(() => {
        setupDOM(`
            <form id="filterForm">
                <select name="target">
                    <option value="5-HT2A">5-HT2A</option>
                </select>
                <button type="submit">Apply</button>
                <button type="reset">Reset</button>
            </form>
            <div id="compound-list"></div>
            <div id="filter-error"></div>
        `);
        require('../app.js');
    });

    test('handles successful filtering', async () => {
        fetch.mockImplementationOnce(() =>
            Promise.resolve({
                ok: true,
                json: () => Promise.resolve({
                    success: true,
                    compounds: testCompounds
                })
            })
        );

        const form = document.getElementById('filterForm');
        fireEvent.submit(form);

        await waitFor(() => {
            expect(fetch).toHaveBeenCalled();
        });
    });

    test('handles filter reset', async () => {
        const form = document.getElementById('filterForm');
        const resetButton = form.querySelector('button[type="reset"]');
        fireEvent.click(resetButton);

        await waitFor(() => {
            expect(fetch).toHaveBeenCalled();
        });
    });
});

describe('Export', () => {
    beforeEach(() => {
        setupDOM(`
            <form id="exportForm">
                <select name="format">
                    <option value="tsv">TSV</option>
                </select>
                <button type="submit">Export</button>
            </form>
            <div id="export-status"></div>
            <div id="export-error"></div>
        `);
        require('../app.js');
    });

    test('handles successful export', async () => {
        const blob = new Blob(['test'], { type: 'text/tab-separated-values' });
        fetch.mockImplementationOnce(() =>
            Promise.resolve({
                ok: true,
                blob: () => Promise.resolve(blob)
            })
        );

        const form = document.getElementById('exportForm');
        fireEvent.submit(form);

        await waitFor(() => {
            expect(fetch).toHaveBeenCalled();
        });
    });

    test('handles export error', async () => {
        fetch.mockImplementationOnce(() =>
            Promise.resolve({
                ok: false,
                json: () => Promise.resolve({
                    error: 'Export failed'
                })
            })
        );

        const form = document.getElementById('exportForm');
        fireEvent.submit(form);

        await waitFor(() => {
            const error = document.getElementById('export-error');
            expect(error.textContent).toBe('Export failed');
        });
    });
});

describe('Event Handlers', () => {
    beforeEach(() => {
        setupDOM(`
            <table class="compound-list">
                <tbody>
                    <tr data-cas="58-08-2">
                        <td>Caffeine</td>
                    </tr>
                </tbody>
            </table>
            <div data-bs-toggle="collapse">
                <i class="fas fa-chevron-down"></i>
            </div>
        `);
        require('../app.js');
    });

    test('handles compound row click', () => {
        const row = document.querySelector('tr');
        fireEvent.click(row);
        expect(window.location.href).toContain('/compounds/58-08-2');
    });

    test('handles collapse toggle', () => {
        const trigger = document.querySelector('[data-bs-toggle="collapse"]');
        const icon = trigger.querySelector('i');
        fireEvent.click(trigger);
        expect(icon.classList.contains('fa-chevron-up')).toBe(true);
    });
});

describe('Utility Functions', () => {
    beforeEach(() => {
        setupDOM(`
            <div id="test-loading"></div>
            <div id="test-error"></div>
            <table class="compound-list">
                <tbody></tbody>
            </table>
        `);
        require('../app.js');
    });

    test('shows and hides loading state', () => {
        showLoading('test-loading');
        expect(document.getElementById('test-loading').classList.contains('loading')).toBe(true);
        hideLoading('test-loading');
        expect(document.getElementById('test-loading').classList.contains('loading')).toBe(false);
    });

    test('shows error message', () => {
        showError('test-error', 'Test error');
        expect(document.getElementById('test-error').textContent).toBe('Test error');
    });

    test('updates compound list', () => {
        updateCompoundList(testCompounds);
        const tbody = document.querySelector('.compound-list tbody');
        expect(tbody.children.length).toBe(2);
        expect(tbody.querySelector('tr').dataset.cas).toBe('58-08-2');
    });
});
