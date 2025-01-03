"""Chemical structure utilities and manipulation."""

from typing import Dict, Any, List, Optional, Tuple, Set
import re

from rdkit import Chem, DataStructs
from rdkit.Chem import (
    AllChem, Draw, rdDepictor, rdFMCS,
    rdMolDescriptors, rdMolTransforms, rdMolAlign,
    Crippen
)
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem.rdDepictor import Compute2DCoords

from logger import LogManager


class StructureUtils:
    """Utilities for chemical structure manipulation and analysis."""
    
    # Common substructure SMARTS patterns
    FUNCTIONAL_GROUPS = {
        'alcohol': '[OH]',
        'amine': '[NH2,NH1,NH0]',
        'carboxyl': '[CX3](=O)[OX2H1]',
        'ether': '[OX2]([#6])[#6]',
        'ester': '[#6][CX3](=O)[OX2][#6]',
        'ketone': '[#6][CX3](=O)[#6]',
        'aldehyde': '[CX3H1](=O)[#6]',
        'amide': '[NX3][CX3](=[OX1])[#6]',
        'nitro': '[NX3](=O)=O',
        'sulfonamide': '[SX4](=[OX1])(=[OX1])([NX3])',
        'phosphate': '[PX4](=[OX1])([OX2])',
        'halogen': '[F,Cl,Br,I]'
    }
    
    # Receptor ligand SMARTS patterns
    RECEPTOR_PATTERNS = {
        # NMDA receptor patterns
        'nmda_core': 'C1CCNCC1',  # Basic piperidine core
        'ketamine': 'O=C1(CCCCC1)NC1CCCCC1',  # Ketamine core
        'pcp': 'C1CCN(CC1)C1CCCCC1',  # PCP core
        'memantine': 'CC12CC3CC(CC(C3)(C1)C2)N',  # Memantine core
        'mk801': 'CC12CC3CC(CC(C3)(C1)C2)NC1=NCCN1',  # MK-801 core
        'dxm': 'COc1ccc2C3CCNCC3COc2c1',  # DXM core
        'nitrous': '[N-]=[N+]=O',  # Nitrous oxide
        'xenon': '[Xe]',  # Xenon
        
        # 5-HT2 patterns
        # Core structures
        'indole': 'c12ccccc1[nH]cc2',  # Basic indole structure
        'tryptamine': 'c12ccccc1[nH]cc2CCN',  # Tryptamine core
        'phenethylamine': 'c1ccccc1CCN',  # Phenethylamine core
        'ergoline': 'C1CN2CCc3c([nH]c4ccccc34)C2C1',  # Ergoline core
        'quinazoline': 'c1ccc2ncnc(c2c1)',  # Quinazoline core
        'benzofuran': 'c1ccc2occc2c1',  # Benzofuran core
        'indazole': 'c1ccc2[nH]ncc2c1',  # Indazole core
        'azepine': 'C1CCCNCc1',  # Azepine core
        'tetrahydropyridine': 'C1=CCNCC1',  # Tetrahydropyridine core
        'pyrrolopyridine': 'c1ccnc2[nH]ccc12',  # Pyrrolopyridine core
        'piperazine': 'C1CNCCN1',  # Added piperazine core
        'piperidine': 'C1CCNCC1',  # Added piperidine core
        'morpholine': 'C1COCCN1',  # Added morpholine core
        'thiophene': 'c1ccsc1',  # Added thiophene core
        'pyrrole': 'c1cc[nH]c1',  # Added pyrrole core
        'imidazole': 'c1c[nH]cn1',  # Added imidazole core
        'oxazole': 'c1cocn1',  # Added oxazole core
        'thiazole': 'c1cscn1',  # Added thiazole core
        'pyrazole': 'c1cn[nH]c1',  # Added pyrazole core
        'triazole': 'c1nc[nH]n1',  # Added triazole core
        
        # Lysergamide patterns
        'lysergamide': 'C1CN2CCc3c([nH]c4ccccc34)[C@H]2CC1C(=O)N',  # Basic lysergamide
        'ald52': 'C1CN2CCc3c([nH]c4ccccc34)[C@H]2CC1C(=O)NC(=O)',  # ALD-52 core
        'eth_lad': 'C1CN2CCc3c([nH]c4ccccc34)[C@H]2CC1C(=O)NCC',  # ETH-LAD core
        '1p_lsd': 'C1CN2CCc3c([nH]c4ccccc34)[C@H]2CC1C(=O)N(C(=O)CC)C',  # Added 1P-LSD
        '1cp_lsd': 'C1CN2CCc3c([nH]c4ccccc34)[C@H]2CC1C(=O)N(C(=O)CCC)C',  # Added 1cP-LSD
        '1v_lsd': 'C1CN2CCc3c([nH]c4ccccc34)[C@H]2CC1C(=O)N(C(=O)CCCC)C',  # Added 1V-LSD
        
        # Phenethylamine patterns
        '2c': 'c1cc(OC)c(cc1CCN)OC',  # 2C-x core
        '2cb': 'c1cc(OC)c(cc1CCN)Br',  # 2C-B specific
        'nbome': 'c1cc(OC)c(cc1CCN)OCc2ccccc2OC',  # NBOMe core
        'nboh': 'c1cc(OC)c(cc1CCN)OCc2ccccc2O',  # NBOH core
        'nbf': 'c1cc(OC)c(cc1CCN)OCc2ccccc2F',  # NBF core
        'nbmd': 'c1cc(OC)c(cc1CCN)OCc2cc3OCOc3cc2',  # NBMD core
        'n1nap': 'c1cc(OC)c(cc1CCNCc2cccc3ccccc23)[N+](=O)[O-]',  # 25N-N1-Nap
        'nbnpome': 'c1cc(OC)c(cc1CCN)OCc2cccc3ccccc23',  # Added NBNPOMe
        'nbnoh': 'c1cc(OC)c(cc1CCN)OCc2cccc3ccccc23O',  # Added NBNOH
        
        # DOx patterns
        'dox': 'c1cc(OC)c(cc1C(C)CN)OC',  # DOx core
        'dom': 'c1cc(OC)c(cc1C(C)CN)OC',  # DOM specific
        'doi': 'c1cc(OC)c(cc1C(C)CN)I',  # DOI specific
        'dob': 'c1cc(OC)c(cc1C(C)CN)Br',  # DOB specific
        'doc': 'c1cc(OC)c(cc1C(C)CN)Cl',  # Added DOC
        'dof': 'c1cc(OC)c(cc1C(C)CN)F',  # Added DOF
        'doet': 'c1cc(OC)c(cc1C(CC)CN)OC',  # Added DOET
        'dopr': 'c1cc(OC)c(cc1C(CCC)CN)OC',  # Added DOPR
        
        # Benzofuran patterns
        'apb': 'c1ccc2c(c1)CC(CN)O2',  # APB core
        'mapb': 'c1ccc2c(c1)CC(CNC)O2',  # MAPB core
        'dhp_dmt': 'c1ccc2c(c1)CC(CCN(C)C)O2',  # 4,5-DHP-DMT
        'eapb': 'c1ccc2c(c1)CC(CNCC)O2',  # Added EAPB
        'bk_2cb': 'c1cc(OC)c(cc1CC(=O)CN)Br',  # Added bk-2C-B
        'bk_ebdp': 'c1ccc2c(c1)CC(C(=O)CNCC)O2',  # Added bk-EBDP
        
        # Novel structures
        'pha57378': 'C1(OCC2)=C(N2C3=C4CCNCC3)C4=CC=C1',  # PHA-57378
        'pnu22394': 'C3CNCCc2c3c1ccccc1n2C',  # PNU-22394
        'tbg': 'CN1CCC2=C(CC1)NC3=C2C=CC(=C3)OC',  # Tabernanthalog
        'dm506': 'CN1CCC2=C(CC1)NC3=CC=CC=C23',  # DM-506
        'iti1549': 'c1ccc2c(c1)n3c(n2)CCNCC3',  # ITI-1549 core
        'al34662': 'CC(N)Cn1ncc2ccc(O)cc12',  # AL-34662
        'al38022a': 'c13CCCOc3ccc(cn2)c1n2CC(C)N',  # AL-38022A
        'ihch7113': 'CN1CCN2[C@H]3CCNC[C@H]3C4=C2C1=CC=C4',  # IHCH-7113
        'r69': 'C[C@@H]1C=C(CNC1)c1c[nH]c2ncccc12',  # (R)-69
        'sn22': 'CN1CCC(CC1)c1c[nH]c2c1cccc2',  # SN-22
        'vu6067416': 'Brc1cc2c(n[nH]c2cc1)C=1CNCCC=1',  # VU6067416
        'rs13449': 'CC1=C2C(=CC=C1)NC=C2C3=CCCNC3',  # RS134-49
        'z3517967757': 'CC(C1=NC=CC=N1)N2CCCC(C2)C3=CC=C(C=C3)O',  # Z3517967757
        'cp132484': 'Cn1cc(c2c1ccc3c2CCCO3)CCN',  # CP-132,484
        'al37350a': 'O2c1ccc3c(c1CCC2)c(c[nH]3)C[C@@H](N)C',  # AL-37350A
        'cp809101': 'CN1CCC[C@H]1Cc2c[nH]c3ccc(F)cc23',  # Added CP-809,101
        'cp135807': 'CN1CCC[C@H]1Cc2c[nH]c3ccc(OC)cc23',  # Added CP-135,807
        'cp93129': 'Cc1ccc2[nH]cc(CCN3CCC[C@H]3C)c2c1',  # Added CP-93,129
        
        # Other important patterns
        'fly': 'C1OC2c3ccccc3OC2C1CN',  # DragonFLY core
        'mescaline': 'c1c(OC)c(OC)c(cc1CCN)OC',  # Mescaline core
        'amphetamine': 'c1ccccc1CC(N)C',  # Amphetamine core
        'mdma': 'c1ccc2c(c1)OCO2CC(NC)C',  # MDMA core
        'bromo_dragonfly': 'c1cc(Br)c2OC3CNCC(O2)C3c1',  # Added Bromo-DragonFLY
        'tcb_2': 'c1cc(OC)c(cc1CCN)C#C',  # Added TCB-2
        'al_lad': 'C1CN2CCc3c([nH]c4ccccc34)[C@H]2CC1C(=O)NCC=C',  # Added AL-LAD
        'pro_lad': 'C1CN2CCc3c([nH]c4ccccc34)[C@H]2CC1C(=O)NC(C)CC',  # Added PRO-LAD
        
        # MPMI and related compounds
        'mpmi': 'c1ccc2c(c1)[nH]cc2CC3CCCN3C',  # MPMI core
        '5_meo_mpmi': 'COc1ccc2c(c1)[nH]cc2CC3CCCN3C',  # 5-MeO-MPMI
        '4_ho_mpmi': 'Oc1ccc2c(c1)[nH]cc2CC3CCCN3C',  # 4-HO-MPMI
        '5f_mpmi': 'Fc1ccc2c(c1)[nH]cc2CC3CCCN3C',  # 5F-MPMI
        '4_aco_mpmi': 'CC(=O)Oc1ccc2c(c1)[nH]cc2CC3CCCN3C',  # Added 4-AcO-MPMI
        '5_meo_dmt': 'COc1ccc2c(c1)[nH]cc2CCN(C)C',  # Added 5-MeO-DMT
        '4_ho_dmt': 'Oc1ccc2c(c1)[nH]cc2CCN(C)C',  # Added 4-HO-DMT
        
        # Pyrrolidine tryptamines
        'pyr_t': 'c1ccc2c(c1)[nH]cc2CCN3CCCC3',  # pyr-T core
        '5_meo_pyr_t': 'COc1ccc2c(c1)[nH]cc2CCN3CCCC3',  # 5-MeO-pyr-T
        '4_ho_pyr_t': 'Oc1ccc2c(c1)[nH]cc2CCN3CCCC3',  # 4-HO-pyr-T
        '5_meo_dpt': 'COc1ccc2c(c1)[nH]cc2CCN(CC)CC',  # Added 5-MeO-DPT
        '4_aco_dmt': 'CC(=O)Oc1ccc2c(c1)[nH]cc2CCN(C)C',  # Added 4-AcO-DMT
        
        # CP series
        'cp135807': 'CN1CCC[C@H]1Cc2c[nH]c3ccc(OC)cc23',  # CP-135,807
        'cp132484': 'Cn1cc(c2c1ccc3c2CCCO3)CCN',  # Added CP-132,484
        'cp93129': 'Cc1ccc2[nH]cc(CCN3CCC[C@H]3C)c2c1',  # Added CP-93,129
        'cp94253': 'Cc1ccc2[nH]cc(CCN3CCC[C@H]3CC)c2c1'  # Added CP-94,253
    }

    
    # Drug class patterns
    DRUG_CLASSES = {
        'nootropics': {
            'racetams': [
                'Piracetam', 'Aniracetam', 'Oxiracetam', 'Pramiracetam',
                'Phenylpiracetam', 'Nefiracetam', 'Coluracetam', 'Fasoracetam'
            ],
            'cholinergics': [
                'Alpha-GPC', 'CDP-Choline', 'Centrophenoxine', 'DMAE',
                'Huperzine A', 'Galantamine'
            ],
            'ampakines': [
                'Sunifiram', 'Unifiram', 'CX-717', 'IDRA-21'
            ],
            'peptides': [
                'Noopept', 'Semax', 'Selank', 'P21', 'Cerebrolysin'
            ]
        },
        'analgesics': {
            'opioids': [
                'Morphine', 'Codeine', 'Oxycodone', 'Hydrocodone',
                'Fentanyl', 'Tramadol', 'Buprenorphine'
            ],
            'nsaids': [
                'Ibuprofen', 'Naproxen', 'Aspirin', 'Diclofenac',
                'Celecoxib', 'Meloxicam'
            ]
        },
        'anesthetics': {
            'general': [
                'Propofol', 'Etomidate', 'Ketamine', 'Sevoflurane',
                'Desflurane', 'Isoflurane'
            ],
            'local': [
                'Lidocaine', 'Bupivacaine', 'Ropivacaine', 'Mepivacaine',
                'Tetracaine', 'Procaine'
            ]
        },
        'stimulants': {
            'amphetamines': [
                'Amphetamine', 'Methamphetamine', 'Lisdexamfetamine',
                'MDMA', 'MDA', 'MDEA'
            ],
            'phenidates': [
                'Methylphenidate', 'Ethylphenidate', 'Isopropylphenidate',
                '4F-MPH', '4-Me-TMP'
            ],
            'eugeroics': [
                'Modafinil', 'Armodafinil', 'Adrafinil', 'Hydrafinil',
                'Flmodafinil'
            ]
        },
        'antidepressants': {
            'ssris': [
                'Fluoxetine', 'Sertraline', 'Paroxetine', 'Citalopram',
                'Escitalopram', 'Fluvoxamine'
            ],
            'snris': [
                'Venlafaxine', 'Desvenlafaxine', 'Duloxetine',
                'Levomilnacipran', 'Milnacipran'
            ],
            'tricyclics': [
                'Amitriptyline', 'Nortriptyline', 'Imipramine',
                'Desipramine', 'Clomipramine'
            ]
        },
        'anxiolytics': {
            'benzodiazepines': [
                'Diazepam', 'Alprazolam', 'Clonazepam', 'Lorazepam',
                'Temazepam', 'Midazolam'
            ],
            'azapirones': [
                'Buspirone', 'Tandospirone', 'Gepirone', 'Ipsapirone'
            ]
        },
        'antipsychotics': {
            'typical': [
                'Haloperidol', 'Chlorpromazine', 'Fluphenazine',
                'Perphenazine', 'Thioridazine'
            ],
            'atypical': [
                'Clozapine', 'Olanzapine', 'Quetiapine', 'Risperidone',
                'Aripiprazole'
            ]
        },
        'mood_stabilizers': {
            'anticonvulsants': [
                'Lamotrigine', 'Valproate', 'Carbamazepine',
                'Oxcarbazepine', 'Gabapentin'
            ],
            'lithium': ['Lithium']
        },
        'cognitive_enhancers': {
            'acetylcholinesterase_inhibitors': [
                'Donepezil', 'Rivastigmine', 'Galantamine',
                'Huperzine A'
            ],
            'nmda_modulators': [
                'Memantine', 'Dextromethorphan', 'Ketamine',
                'Nitrous Oxide'
            ]
        }
    }

    # Additional patterns for specific compound classes
    COMPOUND_CLASSES = {
        'phenethylamines': {
            '2c_series': [
                '2C-B', '2C-C', '2C-D', '2C-E', '2C-I', '2C-N', '2C-P',
                '2C-T-2', '2C-T-7', '2C-T-21', '2C-G', '2C-H', '2C-iP',
                '2C-O', '2C-TFM', '2C-YN', '2C-V', '2C-EF', '2C-G-1',
                '2C-G-2', '2C-G-3', '2C-G-4', '2C-G-5', '2C-G-6', '2C-G-N',
                '2C-T-4', '2C-T-8', '2C-T-9', '2C-T-13', '2C-T-15', '2C-T-17',
                '2C-T-19', '2C-T-20', '2C-T-25'
            ],
            'nbome_series': [
                '25B-NBOMe', '25C-NBOMe', '25I-NBOMe', '25N-NBOMe',
                '25D-NBOMe', '25E-NBOMe', '25G-NBOMe', '25H-NBOMe',
                '25P-NBOMe', '25T-NBOMe', '25TFM-NBOMe', '25CN-NBOMe',
                '25iP-NBOMe'
            ],
            'nboh_series': [
                '25B-NBOH', '25C-NBOH', '25I-NBOH', '25N-NBOH',
                '25D-NBOH', '25E-NBOH', '25H-NBOH', '25P-NBOH',
                '25CN-NBOH', '2C-B-DragonFLY-NBOH'
            ],
            'nbf_series': [
                '25B-NBF', '25C-NBF', '25I-NBF', '25D-NBF',
                '25E-NBF', '25H-NBF', '25P-NBF', '25T2-NBF',
                '25T7-NBF', '25TFM-NBF'
            ],
            'nbmd_series': [
                '25B-NBMD', '25C-NBMD', '25I-NBMD', '25D-NBMD',
                '25E-NBMD', '25H-NBMD', '25P-NBMD', '25T2-NBMD',
                '25T7-NBMD', '25TFM-NBMD'
            ],
            'other_derivatives': [
                '25B-N1POMe', '25B-NAcPip', '25B-NB23DM', '25B-NB25DM',
                '25C-NBCl', '25C-NBOEt', '25C-NBOiPr', '25I-N2Nap1OH',
                '25I-N3MT2M', '25I-N4MT3M', '25I-NB34MD', '25I-NBAm',
                '25I-NBBr', '25I-NBMeOH', '25I-NBTFM'
            ]
        },
        'tryptamines': {
            'pyrrolidine_tryptamines': [
                'MPMI', '5-MeO-MPMI', '4-HO-MPMI', '5F-MPMI',
                'pyr-T', '5-MeO-pyr-T', '4-HO-pyr-T',
                'CP-135,807'
            ],
            'base_tryptamines': [
                'DMT', 'DET', 'DPT', 'DiPT', 'MiPT', 'EiPT',
                'DALT', '4-HO-DALT', '5-MeO-DALT', 'DBT', 'DCPT',
                'EPT', 'MPT', 'PiPT'
            ],
            'substituted_tryptamines': [
                '4-HO-DMT', '4-AcO-DMT', '5-MeO-DMT', '4-HO-MET',
                '4-AcO-MET', '5-MeO-MET', '4-HO-DET', '4-AcO-DET',
                '5-MeO-DET', '4-HO-DiPT', '4-AcO-DiPT', '5-MeO-DiPT',
                '4-HO-MiPT', '4-AcO-MiPT', '5-MeO-MiPT', '4-HO-EPT',
                '4-HO-McPT', '4-HO-MPT', '5-MeO-EiPT', '5-MeO-MALT',
                '5-MeO-MPMI', '4-HO-MPMI', '5F-MPMI'  # Added MPMI series
            ],
            'alpha_alkyltryptamines': [
                'αMT', 'α-ET', '5-MeO-αMT', '5-MeO-α-ET',
                '4,5-DHP-α-MT', '4-HO-αMT', '5-F-αMT'
            ]
        },
        'lysergamides': {
            'lsd_analogs': [
                'LSD', '1P-LSD', 'ALD-52', 'ETH-LAD', 'AL-LAD',
                '1cP-LSD', '1V-LSD', '1B-LSD', 'LSZ', 'LSA',
                'PRO-LAD', 'PARGY-LAD', 'MIPLA', 'LAMPA', 'LSH',
                'LSM-775', 'LSD-Pip', 'MLD-41', 'ECPLA', 'BU-LAD',
                'LAE-32', 'LSP', 'LPD-824'
            ]
        },
        'benzofurans': {
            'apb_series': [
                '5-APB', '6-APB', '5-MAPB', '6-MAPB',
                '5-EAPB', '6-EAPB', '5-APDB', '6-APDB',
                '5-MeO-BFE', '5-MeO-DiBF'
            ],
            'benzofuran_others': [
                '2C-B-FLY', '2C-B-DragonFLY', 'Bromo-DragonFLY',
                'TFMFly', '2C-E-FLY', '2CBFly-NBOMe', 'F-2',
                'F-22', 'TFMFly'
            ]
        },
        'novel_compounds': {
            'non_hallucinogenic': [
                'PHA-57378', 'PNU-22394', 'Tabernanthalog',
                'DM-506', 'ITI-1549', 'AL-34662', 'AAZ-A-154',
                'IHCH-7086', 'Z3517967757', 'CP-135,807'  # Added CP-135,807
            ],
            'research_compounds': [
                'IHCH-7079', 'IHCH-7086', 'IHCH-7113',
                'AL-38022A', 'Z3517967757', 'RS134-49',
                'VU6067416', '25N-N1-Nap', '(R)-69', 'SN-22',
                'CP-132,484', 'AL-37350A', '4,5-DHP-DMT'
            ]
        }
    }
    
    
    def __init__(self):
        """Initialize structure utilities."""
        self.logger = LogManager().get_logger("structure_utils")

    def standardize_smiles(self, smiles: str) -> Optional[str]:
        """
        Standardize SMILES string to canonical form.
        
        Args:
            smiles: Input SMILES string
            
        Returns:
            Standardized SMILES string or None if invalid
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
                
            # Remove hydrogens
            mol = Chem.RemoveHs(mol)
            
            # Kekulize
            Chem.Kekulize(mol)
            
            # Generate 2D coordinates
            rdDepictor.Compute2DCoords(mol)
            
            # Return canonical SMILES
            return Chem.MolToSmiles(
                mol,
                isomericSmiles=True,
                canonical=True,
                kekuleSmiles=True
            )
            
        except Exception as e:
            self.logger.error(f"Error standardizing SMILES: {str(e)}")
            return None

    def get_substructures(self, mol: Chem.Mol) -> Dict[str, int]:
        """
        Identify common substructures in molecule.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Dictionary of substructure counts
        """
        counts = {}
        
        try:
            # Check common functional groups
            for name, smarts in self.FUNCTIONAL_GROUPS.items():
                pattern = Chem.MolFromSmarts(smarts)
                if pattern:
                    matches = mol.GetSubstructMatches(pattern)
                    counts[name] = len(matches)
            
            # Check 5-HT2 receptor ligand patterns
            for name, smarts in self.SEROTONIN_PATTERNS.items():
                pattern = Chem.MolFromSmarts(smarts)
                if pattern:
                    matches = mol.GetSubstructMatches(pattern)
                    counts[f"serotonin_{name}"] = len(matches)
                    
        except Exception as e:
            self.logger.error(f"Error getting substructures: {str(e)}")
            
        return counts


    def is_potential_ligand(self, mol: Chem.Mol) -> Tuple[bool, List[str]]:
        """
        Check if molecule has structural features common to receptor ligands.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Tuple of (is_potential_ligand, matching_patterns)
        """
        matching_patterns = []
        
        try:
            # Check each receptor pattern
            for name, smarts in self.RECEPTOR_PATTERNS.items():
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    matching_patterns.append(name)
            
            # Additional checks for common features
            descriptors = {
                'MW': rdMolDescriptors.CalcExactMolWt(mol),
                'LogP': Crippen.MolLogP(mol),
                'TPSA': Chem.rdMolDescriptors.CalcTPSA(mol),
                'RotBonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
                'HBA': rdMolDescriptors.CalcNumHBA(mol),  # Added Hydrogen Bond Acceptors
                'HBD': rdMolDescriptors.CalcNumHBD(mol),  # Added Hydrogen Bond Donors
                'Rings': rdMolDescriptors.CalcNumRings(mol),  # Added Ring Count
                'ArRings': rdMolDescriptors.CalcNumAromaticRings(mol),  # Added Aromatic Ring Count
                'SP3': rdMolDescriptors.CalcFractionCSP3(mol),  # Added SP3 Character
                'MR': Crippen.MolMR(mol)  # Added Molar Refractivity
            }
            
            # Typical ranges for 5-HT2 ligands based on known actives
            if (
                200 <= descriptors['MW'] <= 600 and      # Molecular weight range
                1 <= descriptors['LogP'] <= 6 and        # LogP range
                20 <= descriptors['TPSA'] <= 90 and      # TPSA range
                descriptors['RotBonds'] <= 7 and         # Rotatable bonds limit
                1 <= descriptors['HBA'] <= 7 and         # H-bond acceptors range
                0 <= descriptors['HBD'] <= 3 and         # H-bond donors range
                1 <= descriptors['Rings'] <= 5 and       # Ring count range
                1 <= descriptors['ArRings'] <= 3 and     # Aromatic rings range
                0.2 <= descriptors['SP3'] <= 0.8 and     # SP3 character range
                40 <= descriptors['MR'] <= 150           # Molar refractivity range
            ):
                matching_patterns.append('physicochemical_properties')
            
            # Check for key functional groups
            for name, smarts in self.FUNCTIONAL_GROUPS.items():
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    matching_patterns.append(f'functional_group_{name}')
            
            # Check for specific compound classes
            for class_name, compounds in self.COMPOUND_CLASSES.items():
                for subclass, patterns in compounds.items():
                    for pattern in patterns:
                        if pattern.lower() in Chem.MolToSmiles(mol).lower():
                            matching_patterns.append(f'{class_name}_{subclass}')
            
            return bool(matching_patterns), matching_patterns
            
        except Exception as e:
            self.logger.error(f"Error checking 5-HT2 ligand potential: {str(e)}")
            return False, []

    def get_scaffold(self, mol: Chem.Mol) -> Optional[str]:
        """
        Get Murcko scaffold SMILES.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Scaffold SMILES or None if failed
        """
        try:
            scaffold = AllChem.MurckoDecompose(mol)
            return Chem.MolToSmiles(scaffold) if scaffold else None
            
        except Exception as e:
            self.logger.error(f"Error getting scaffold: {str(e)}")
            return None

    def get_fragments(self, mol: Chem.Mol) -> List[str]:
        """
        Break molecule into fragments at rotatable bonds.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            List of fragment SMILES
        """
        fragments = []
        
        try:
            # Find rotatable bonds
            rot_bonds = Chem.rdMolDescriptors.FindAllRotatableBonds(mol)
            
            # Break at each rotatable bond
            for bond_idx in rot_bonds:
                bond = mol.GetBondWithIdx(bond_idx)
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # Get atoms connected by bond
                    atom1 = bond.GetBeginAtom()
                    atom2 = bond.GetEndAtom()
                    
                    # Break bond and get fragments
                    fragments.extend(
                        Chem.MolToSmiles(frag)
                        for frag in Chem.rdmolops.GetMolFrags(
                            mol,
                            asMols=True,
                            sanitizeFrags=False
                        )
                    )
                    
        except Exception as e:
            self.logger.error(f"Error getting fragments: {str(e)}")
            
        return list(set(fragments))  # Remove duplicates

    def align_structures(
        self,
        ref_mol: Chem.Mol,
        probe_mol: Chem.Mol
    ) -> Tuple[float, Chem.Mol]:
        """
        Align probe molecule to reference using MCS.
        
        Args:
            ref_mol: Reference molecule
            probe_mol: Probe molecule to align
            
        Returns:
            Tuple of (RMSD, aligned molecule)
        """
        try:
            # Find maximum common substructure
            mcs = rdFMCS.FindMCS([ref_mol, probe_mol])
            if mcs and mcs.numAtoms > 0:
                # Get atom mappings
                mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
                ref_match = ref_mol.GetSubstructMatch(mcs_mol)
                probe_match = probe_mol.GetSubstructMatch(mcs_mol)
                
                if ref_match and probe_match:
                    # Create conformers if needed
                    if not ref_mol.GetNumConformers():
                        AllChem.EmbedMolecule(ref_mol)
                    if not probe_mol.GetNumConformers():
                        AllChem.EmbedMolecule(probe_mol)
                    
                    # Align using MCS atoms
                    rmsd = rdMolAlign.AlignMol(
                        probe_mol,
                        ref_mol,
                        atomMap=list(zip(probe_match, ref_match))
                    )
                    return rmsd, probe_mol
                    
        except Exception as e:
            self.logger.error(f"Error aligning structures: {str(e)}")
            
        return float('inf'), probe_mol

    def generate_conformers(
        self,
        mol: Chem.Mol,
        n_conf: int = 10,
        optimize: bool = True
    ) -> List[float]:
        """
        Generate multiple conformers and return their energies.
        
        Args:
            mol: RDKit molecule
            n_conf: Number of conformers to generate
            optimize: Whether to optimize conformers
            
        Returns:
            List of conformer energies
        """
        energies = []
        
        try:
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Generate conformers
            AllChem.EmbedMultipleConfs(
                mol,
                numConfs=n_conf,
                randomSeed=42,
                pruneRmsThresh=0.5
            )
            
            if optimize:
                # Optimize each conformer
                for conf_id in range(mol.GetNumConformers()):
                    # MMFF optimization
                    AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
                    
                    # Calculate energy
                    props = AllChem.MMFFGetMoleculeProperties(mol)
                    ff = AllChem.MMFFGetMoleculeForceField(
                        mol,
                        props,
                        confId=conf_id
                    )
                    if ff:
                        energy = ff.CalcEnergy()
                        energies.append(energy)
                        
        except Exception as e:
            self.logger.error(f"Error generating conformers: {str(e)}")
            
        return energies

    def get_fingerprint_similarity(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol,
        fp_type: str = 'morgan'
    ) -> float:
        """
        Calculate fingerprint-based similarity between molecules.
        
        Args:
            mol1: First molecule
            mol2: Second molecule
            fp_type: Fingerprint type ('morgan', 'maccs', 'topological')
            
        Returns:
            Tanimoto similarity score
        """
        try:
            if fp_type == 'morgan':
                fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
                fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
            elif fp_type == 'maccs':
                fp1 = AllChem.GetMACCSKeysFingerprint(mol1)
                fp2 = AllChem.GetMACCSKeysFingerprint(mol2)
            elif fp_type == 'topological':
                fp1 = Chem.RDKFingerprint(mol1)
                fp2 = Chem.RDKFingerprint(mol2)
            else:
                raise ValueError(f"Unknown fingerprint type: {fp_type}")
                
            return DataStructs.TanimotoSimilarity(fp1, fp2)
            
        except Exception as e:
            self.logger.error(f"Error calculating similarity: {str(e)}")
            return 0.0

    def enumerate_tautomers(self, mol: Chem.Mol) -> List[str]:
        """
        Enumerate possible tautomers of molecule.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            List of tautomer SMILES
        """
        tautomers = set()
        
        try:
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Define SMARTS patterns for common tautomeric transformations
            patterns = [
                # Keto-enol
                ('[CX3](=[OX1])[CX4]', '[CX3]([OX2H1])=[CX3]'),
                # Imine-enamine
                ('[CX3]=[NX2]-[CX4]', '[CX4]-[NX3]-[CX3]'),
                # Amide-imidic acid
                ('[CX3](=[OX1])[NX3]', '[CX3]([OX2H1])=[NX2]')
            ]
            
            # Apply each transformation
            for pattern in patterns:
                reactant_smarts, product_smarts = pattern
                
                # Create reaction
                rxn = AllChem.ReactionFromSmarts(
                    f'{reactant_smarts}>>{product_smarts}'
                )
                
                # Apply reaction
                products = rxn.RunReactants((mol,))
                for product_tuple in products:
                    for product in product_tuple:
                        try:
                            Chem.SanitizeMol(product)
                            tautomers.add(Chem.MolToSmiles(product))
                        except Exception:
                            continue
                            
        except Exception as e:
            self.logger.error(f"Error enumerating tautomers: {str(e)}")
            
        return list(tautomers)

    def standardize_stereochemistry(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Standardize stereochemistry representation.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Molecule with standardized stereochemistry
        """
        try:
            # Find all stereocenters
            stereo_centers = Chem.FindMolChiralCenters(
                mol,
                includeUnassigned=True
            )
            
            # Assign random stereochemistry to unspecified centers
            Chem.AssignStereochemistry(mol, force=True)
            
            # Clean up
            Chem.AssignStereochemistry(mol, cleanIt=True)
            
            return mol
            
        except Exception as e:
            self.logger.error(f"Error standardizing stereochemistry: {str(e)}")
            return mol

    def get_largest_fragment(self, mol: Chem.Mol) -> Optional[str]:
        """
        Get SMILES of largest fragment in molecule.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            SMILES of largest fragment or None if failed
        """
        try:
            # Split into fragments
            fragments = Chem.GetMolFrags(mol, asMols=True)
            if not fragments:
                return None
                
            # Find largest fragment by number of atoms
            largest = max(fragments, key=lambda m: m.GetNumAtoms())
            
            return Chem.MolToSmiles(largest)
            
        except Exception as e:
            self.logger.error(f"Error getting largest fragment: {str(e)}")
            return None

    def neutralize_charges(self, mol: Chem.Mol) -> Optional[str]:
        """
        Neutralize formal charges where possible.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            SMILES of neutralized molecule or None if failed
        """
        try:
            # Patterns for neutralization
            patterns = [
                # Carboxylate to carboxylic acid
                ('[CX3](=O)[O-]>>[CX3](=O)[OH]', ),
                # Ammonium to amine
                ('[NX4+]>>[NX3]', ),
                # Sulfonate to sulfonic acid
                ('[SX4](=O)(=O)[O-]>>[SX4](=O)(=O)[OH]', ),
                # Phosphonate to phosphonic acid
                ('[PX4](=O)([O-])[O-]>>[PX4](=O)([OH])[OH]', )
            ]
            
            # Apply each neutralization pattern
            for pattern in patterns:
                rxn = AllChem.ReactionFromSmarts(pattern[0])
                products = rxn.RunReactants((mol,))
                if products:
                    # Take first product
                    product = products[0][0]
                    try:
                        Chem.SanitizeMol(product)
                        mol = product
                    except Exception:
                        continue
                        
            return Chem.MolToSmiles(mol)
            
        except Exception as e:
            self.logger.error(f"Error neutralizing charges: {str(e)}")
            return None
