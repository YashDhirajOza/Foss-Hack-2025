import pandas as pd
from rdkit import Chem
from torch.utils.data import Dataset
import requests
import os

class OpenReactionDataset:
    def __init__(self, cache_dir="data/cache"):
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        
        # Available datasets
        self.datasets = {
            "uspto": "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/uspto_data.csv",
            "ord": "https://github.com/open-reaction-database/ord-data/releases/download/1.0.0/ord_dataset_1.0.0.csv",
            "reaxys": "https://www.reaxys.com/api/dataset/sample.csv"  # Example URL
        }
    
    def load_dataset(self, name="uspto", force_download=False):
        """Load reaction dataset from source or cache"""
        cache_file = os.path.join(self.cache_dir, f"{name}_data.csv")
        
        if not os.path.exists(cache_file) or force_download:
            if name in self.datasets:
                # Download dataset
                response = requests.get(self.datasets[name])
                with open(cache_file, 'wb') as f:
                    f.write(response.content)
            else:
                raise ValueError(f"Unknown dataset: {name}")
        
        # Load and preprocess dataset
        df = pd.read_csv(cache_file)
        return self._preprocess_dataset(df)
    
    def _preprocess_dataset(self, df):
        """Preprocess reaction dataset"""
        processed_data = []
        
        for _, row in df.iterrows():
            try:
                # Convert SMILES to RDKit molecules
                reactants = Chem.MolFromSmiles(row['reactants'])
                products = Chem.MolFromSmiles(row['products'])
                
                if reactants and products:
                    processed_data.append({
                        'reactants_mol': reactants,
                        'products_mol': products,
                        'reaction_type': row.get('reaction_type', 'unknown'),
                        'yield': row.get('yield', None),
                        'conditions': row.get('conditions', {})
                    })
            except:
                continue
        
        return processed_data

class ReactionDataModule(Dataset):
    def __init__(self, data):
        self.data = data
    
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        item = self.data[idx]
        # Convert molecules to feature vectors
        reactant_fp = self._get_fingerprint(item['reactants_mol'])
        product_fp = self._get_fingerprint(item['products_mol'])
        
        return {
            'reactants': reactant_fp,
            'products': product_fp,
            'reaction_type': item['reaction_type'],
            'yield': item.get('yield', 0)
        }
    
    def _get_fingerprint(self, mol, size=2048):
        """Generate Morgan fingerprint for molecule"""
        return list(Chem.AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=size))
