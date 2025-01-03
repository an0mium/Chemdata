#!/usr/bin/env python3
"""Flask application for browsing compound data."""

from flask import Flask, render_template, jsonify, request, send_file
import pandas as pd
import os
from typing import List, Dict, Any
import json
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from io import BytesIO

app = Flask(__name__)

# Configure pandas display options
pd.set_option('display.max_colwidth', None)


class DataBrowser:
    def __init__(self, data_dir: str = "../data"):
        """Initialize with data directory path."""
        self.data_dir = data_dir
        self.current_file = self._get_latest_tsv()
        self.df = pd.read_csv(self.current_file, sep='\t')
        
        # Convert SMILES and InChI to strings
        self.df['smiles'] = self.df['smiles'].fillna('').astype(str)
        self.df['inchi'] = self.df['inchi'].fillna('').astype(str)
        
        # Get column info
        self.columns = self.df.columns.tolist()
        self.url_columns = [col for col in self.columns if col.endswith('_url')]
        
        # Define column groups
        self.column_groups = {
            'identifiers': ['name', 'cas_number', 'smiles', 'inchi', 'molecular_formula'],
            'properties': ['molecular_weight', 'logp', 'hbd', 'hba', 'tpsa', 'rotatable_bonds'],
            'activity': ['target_type', 'activity_type', 'activity_value', 'activity_unit'],
            'pharmacology': ['mechanism', 'primary_target', 'secondary_targets'],
            'legal': ['legal_status', 'scheduling', 'controlled_status'],
            'references': self.url_columns
        }
        
    def _get_latest_tsv(self) -> str:
        """Get path of most recent TSV file."""
        tsv_files = [
            f for f in os.listdir(self.data_dir) 
            if f.startswith('receptor_compounds_') and f.endswith('.tsv')
        ]
        if not tsv_files:
            raise FileNotFoundError("No compound data files found")
            
        latest = max(tsv_files)
        return os.path.join(self.data_dir, latest)
        
    def get_compounds(
        self,
        search: str = '',
        target_type: str = '',
        structure_search: str = '',
        similarity_threshold: float = 0.7,
        activity_type: str = '',
        legal_status: str = '',
        selected_columns: List[str] = None,
        page: int = 1,
        per_page: int = 50
    ) -> Dict[str, Any]:
        """Get paginated compound data with optional filtering."""
        df = self.df
        
        # Apply filters
        if search:
            df = df[
                df['name'].str.contains(search, case=False) |
                df['cas_number'].str.contains(search, case=False) |
                df['molecular_formula'].str.contains(search, case=False)
            ]
            
        if target_type:
            df = df[df['target_type'] == target_type]
            
        if structure_search:
            query_mol = Chem.MolFromSmiles(structure_search)
            if query_mol:
                query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2)
                
                def calc_similarity(smiles):
                    try:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
                            return DataStructs.TanimotoSimilarity(query_fp, fp)
                    except Exception as e:
                        print(f"Error calculating similarity: {str(e)}")
                        return 0
                    return 0
                
                df['similarity'] = df['smiles'].apply(calc_similarity)
                df = df[df['similarity'] >= similarity_threshold]
                
        if activity_type:
            df = df[df['activity_type'] == activity_type]
            
        if legal_status:
            df = df[df['legal_status'] == legal_status]
            
        # Select columns
        if selected_columns:
            df = df[selected_columns]
            
        # Calculate pagination
        total = len(df)
        start = (page - 1) * per_page
        end = start + per_page
        
        # Get page of data
        page_df = df.iloc[start:end]
        
        # Convert to list of dicts
        compounds = []
        for _, row in page_df.iterrows():
            compound = row.to_dict()
            
            # Add URL links
            compound['urls'] = {
                col: url for col, url in row.items()
                if col in self.url_columns and url
            }
            
            compounds.append(compound)
            
        return {
            'compounds': compounds,
            'total': total,
            'page': page,
            'per_page': per_page,
            'pages': (total + per_page - 1) // per_page
        }
        
    def get_target_types(self) -> List[str]:
        """Get list of unique target types."""
        return sorted(self.df['target_type'].unique().tolist())
        
    def get_activity_types(self) -> List[str]:
        """Get list of unique activity types."""
        return sorted(self.df['activity_type'].unique().tolist())
        
    def get_legal_statuses(self) -> List[str]:
        """Get list of unique legal statuses."""
        return sorted(self.df['legal_status'].unique().tolist())
        
    def get_column_groups(self) -> Dict[str, List[str]]:
        """Get column groups."""
        return self.column_groups
        
    def get_stats(self) -> Dict[str, Any]:
        """Get summary statistics."""
        return {
            'total_compounds': len(self.df),
            'with_activity': len(self.df[self.df['activity_data'].str.len() > 0]),
            'by_target': self.df['target_type'].value_counts().to_dict(),
            'by_activity': self.df['activity_type'].value_counts().to_dict(),
            'by_legal_status': self.df['legal_status'].value_counts().to_dict()
        }
        
    def structure_to_image(self, smiles: str) -> bytes:
        """Convert SMILES to PNG image."""
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol)
            img_io = BytesIO()
            img.save(img_io, 'PNG')
            img_io.seek(0)
            return img_io.getvalue()
        return None


# Initialize data browser
browser = DataBrowser()


@app.route('/')
def index():
    """Render main page."""
    return render_template(
        'index.html',
        target_types=browser.get_target_types(),
        activity_types=browser.get_activity_types(),
        legal_statuses=browser.get_legal_statuses(),
        column_groups=browser.get_column_groups(),
        stats=browser.get_stats()
    )


@app.route('/api/compounds')
def get_compounds():
    """API endpoint for compound data."""
    search = request.args.get('search', '')
    target = request.args.get('target', '')
    structure = request.args.get('structure', '')
    similarity = float(request.args.get('similarity', 0.7))
    activity = request.args.get('activity', '')
    legal = request.args.get('legal', '')
    columns = request.args.getlist('columns')
    page = int(request.args.get('page', 1))
    per_page = int(request.args.get('per_page', 50))
    
    data = browser.get_compounds(
        search=search,
        target_type=target,
        structure_search=structure,
        similarity_threshold=similarity,
        activity_type=activity,
        legal_status=legal,
        selected_columns=columns,
        page=page,
        per_page=per_page
    )
    return jsonify(data)


@app.route('/api/structure/<smiles>')
def get_structure_image(smiles: str):
    """API endpoint for structure images."""
    img_data = browser.structure_to_image(smiles)
    if img_data:
        return send_file(
            BytesIO(img_data),
            mimetype='image/png'
        )
    return '', 404


if __name__ == '__main__':
    app.run(debug=True, port=5000)