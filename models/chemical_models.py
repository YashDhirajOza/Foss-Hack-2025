from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import torch
from transformers import AutoModelForSequenceClassification, AutoTokenizer
from data.chemical_data import REACTION_TEMPLATES, CHEMICALS_DATABASE  # Changed from SAFETY_DATABASE

class ChemicalPredictor:
    def __init__(self):
        # Use expanded reaction templates
        self.reaction_template_db = REACTION_TEMPLATES
        
        # Load Hugging Face's ChemBERTa for reaction classification
        try:
            model_name = "seyonec/ChemBERTa-zinc-base-v1"
            self.tokenizer = AutoTokenizer.from_pretrained(model_name)
            self.model = AutoModelForSequenceClassification.from_pretrained(model_name)
        except:
            self.tokenizer = None
            self.model = None
    
    def predict_reaction(self, reactants_smiles):
        """Predict products using RDKit and reaction templates"""
        try:
            # Convert SMILES to RDKit molecules
            mols = [Chem.MolFromSmiles(s.strip()) for s in reactants_smiles.split('.')]
            
            # Simple template matching
            combined_smiles = '.'.join(sorted(reactants_smiles.split('.')))
            if combined_smiles in self.reaction_template_db:
                template = self.reaction_template_db[combined_smiles]
                product_mol = Chem.MolFromSmiles(template['product'])
                return {
                    'success': True,
                    'product_smiles': template['product'],
                    'product_img': Draw.MolToImage(product_mol),
                    'reaction_name': template['name'],
                    'conditions': template['conditions']
                }
            
            return {'success': False, 'error': 'No matching reaction template found'}
            
        except Exception as e:
            return {'success': False, 'error': str(e)}
    
    def analyze_reaction_feasibility(self, smiles):
        """Use ChemBERTa to analyze reaction feasibility"""
        if not self.tokenizer or not self.model:
            return {'feasibility': 'Model not loaded'}
            
        try:
            inputs = self.tokenizer(smiles, return_tensors="pt", padding=True, truncation=True)
            outputs = self.model(**inputs)
            prediction = torch.nn.functional.softmax(outputs.logits, dim=1)
            feasibility_score = prediction[0][1].item()
            
            return {
                'feasibility_score': round(feasibility_score * 100, 2),
                'confidence': 'High' if feasibility_score > 0.8 else 'Medium' if feasibility_score > 0.5 else 'Low'
            }
        except:
            return {'feasibility': 'Error in analysis'}
