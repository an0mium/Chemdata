// Import testing libraries and utilities
import '@testing-library/jest-dom';
import { configure } from '@testing-library/dom';
import { TextEncoder, TextDecoder } from 'util';

// Configure testing library
configure({
    testIdAttribute: 'data-testid',
    asyncUtilTimeout: 2000,
});

// Mock browser APIs
global.TextEncoder = TextEncoder;
global.TextDecoder = TextDecoder;

// Mock window location
Object.defineProperty(window, 'location', {
    value: {
        href: '',
        pathname: '',
        search: '',
        hash: '',
        assign: jest.fn(),
        replace: jest.fn()
    },
    writable: true
});

// Mock storage APIs
const storageMock = {
    getItem: jest.fn(),
    setItem: jest.fn(),
    removeItem: jest.fn(),
    clear: jest.fn()
};
Object.defineProperty(window, 'localStorage', { value: storageMock });
Object.defineProperty(window, 'sessionStorage', { value: storageMock });

// Mock fetch API with binding_data_processor specific defaults
global.fetch = jest.fn(() =>
    Promise.resolve({
        ok: true,
        json: () => Promise.resolve({
            compounds: [],
            predictions: {},
            analysis: {},
            web_data: {}
        }),
        text: () => Promise.resolve(''),
        blob: () => Promise.resolve(new Blob()),
    })
);

// Mock observers
global.ResizeObserver = class ResizeObserver {
    observe() { }
    unobserve() { }
    disconnect() { }
};

global.IntersectionObserver = class IntersectionObserver {
    observe() { }
    unobserve() { }
    disconnect() { }
};

// Mock window.matchMedia
Object.defineProperty(window, 'matchMedia', {
    writable: true,
    value: jest.fn().mockImplementation(query => ({
        matches: false,
        media: query,
        onchange: null,
        addListener: jest.fn(),
        removeListener: jest.fn(),
        addEventListener: jest.fn(),
        removeEventListener: jest.fn(),
        dispatchEvent: jest.fn(),
    })),
});

// Mock console with binding_data_processor specific formatting
global.console = {
    ...console,
    error: jest.fn((message, ...args) => {
        if (message.includes('binding_data_processor')) {
            console._error(`[BDP Error] ${message}`, ...args);
        }
    }),
    warn: jest.fn(),
    info: jest.fn(),
    debug: jest.fn()
};

// Custom matchers for binding_data_processor
expect.extend({
    toHaveBeenCalledWithCompound(received, compound) {
        const pass = received.mock.calls.some(call =>
            call.some(arg => arg?.cas_number === compound.cas_number)
        );
        return {
            pass,
            message: () =>
                `expected ${received.getMockName()} to have been called with compound ${compound.cas_number}`,
        };
    },
    toHaveBeenCalledWithStructure(received, smiles) {
        const pass = received.mock.calls.some(call =>
            call.some(arg => arg?.smiles === smiles)
        );
        return {
            pass,
            message: () =>
                `expected ${received.getMockName()} to have been called with structure ${smiles}`,
        };
    },
    toHaveValidBindingData(received) {
        const hasValidData = received?.binding_data?.every(data =>
            typeof data.affinity === 'number' &&
            typeof data.confidence === 'number' &&
            typeof data.target === 'string'
        );
        return {
            pass: hasValidData,
            message: () => 'expected compound to have valid binding data',
        };
    },
});

// Mock RDKit with binding_data_processor specific functionality
jest.mock('rdkit', () => ({
    __esModule: true,
    default: {
        Molecule: jest.fn().mockImplementation(() => ({
            toSmiles: jest.fn().mockReturnValue('C'),
            toMolBlock: jest.fn().mockReturnValue(''),
            getNumAtoms: jest.fn().mockReturnValue(1),
            calculateProperties: jest.fn().mockReturnValue({
                mw: 100,
                logp: 1,
                psa: 50,
                hba: 2,
                hbd: 1,
                rotatable_bonds: 0
            }),
        })),
        MoleculeDetails: jest.fn().mockImplementation(() => ({
            getMolecularWeight: jest.fn().mockReturnValue(100),
            getLogP: jest.fn().mockReturnValue(1),
            getPSA: jest.fn().mockReturnValue(50),
            getHBondDonorCount: jest.fn().mockReturnValue(1),
            getHBondAcceptorCount: jest.fn().mockReturnValue(2),
            getRotatableBondCount: jest.fn().mockReturnValue(0),
        })),
    },
}));

// Mock Plotly with binding_data_processor specific defaults
jest.mock('plotly.js', () => ({
    __esModule: true,
    default: {
        newPlot: jest.fn(),
        react: jest.fn(),
        purge: jest.fn(),
        relayout: jest.fn(),
        restyle: jest.fn(),
    },
}));

// Cleanup
afterEach(() => {
    jest.clearAllMocks();
    fetch.mockClear();
    localStorage.clear();
    sessionStorage.clear();
});
