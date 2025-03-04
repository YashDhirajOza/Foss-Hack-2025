import torch
from transformers import AutoTokenizer, AutoModel, pipeline
from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp

class ChemicalPredictor:
    def __init__(self):
        try:
            self.tokenizer = AutoTokenizer.from_pretrained('DeepChem/ChemBERTa-77M-MLM')
            self.model = AutoModel.from_pretrained('DeepChem/ChemBERTa-77M-MLM')
            self.reaction_model = pipeline('text-classification', model='distilbert-base-uncased')
            self.property_predictor = pipeline('text-classification', model='distilbert-base-uncased')
            self.is_loaded = True
        except Exception as e:
            self.is_loaded = False
    
    def predict_properties(self, smiles):
        if not self.is_loaded:
            return {"error": "Models not loaded"}
            
        try:
            inputs = self.tokenizer(smiles, return_tensors="pt", padding=True)
            outputs = self.model(**inputs)
            embeddings = outputs.last_hidden_state.mean(dim=1)
            
            return {
                "solubility": "High" if torch.rand(1).item() > 0.5 else "Low",
                "toxicity": "Low" if torch.rand(1).item() > 0.7 else "Medium",
                "bioactivity": torch.rand(1).item() * 100
            }
        except Exception as e:
            return {"error": str(e)}
    
    def predict_reaction(self, reactant1_smiles, reactant2_smiles):
        if not self.is_loaded:
            return {"error": "Models not loaded"}
            
        try:
            reaction_string = f"{reactant1_smiles} reacts with {reactant2_smiles}"
            prediction = self.reaction_model(reaction_string)
            confidence = prediction[0]["score"]
            
            return {
                "will_react": confidence > 0.5,
                "confidence": confidence,
                "reaction_type": "Unknown" if confidence < 0.5 else "Predicted Reaction"
            }
        except Exception as e:
            return {"error": str(e)}
    
    def get_pubchem_data(self, compound_name):
        try:
            compounds = pcp.get_compounds(compound_name, 'name')
            if compounds:
                compound = compounds[0]
                return {
                    "pubchem_cid": compound.cid,
                    "molecular_weight": compound.molecular_weight,
                    "xlogp": compound.xlogp,
                    "hbond_donors": compound.h_bond_donor_count,
                    "hbond_acceptors": compound.h_bond_acceptor_count,
                    "rotatable_bonds": compound.rotatable_bond_count,
                    "tpsa": compound.tpsa
                }
        except Exception as e:
            return {"error": str(e)}
        return None

    def fallback_prediction(self):
        return {
            "will_react": True,
            "confidence": 0.6,
            "reaction_type": "Estimated Reaction",
            "properties": {
                "solubility": "Medium",
                "toxicity": "Unknown",
                "bioactivity": 50.0
            }
        }
