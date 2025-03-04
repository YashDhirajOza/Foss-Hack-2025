import torch
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

class MLReactionPredictor:
    def __init__(self):
        self.tokenizer = AutoTokenizer.from_pretrained("seyonec/ChemBERTa-zinc-base-v1")
        self.model = AutoModelForSequenceClassification.from_pretrained("seyonec/ChemBERTa-zinc-base-v1")
        self.fingerprint_size = 2048
        
    def predict_unknown_reaction(self, reactant1, reactant2):
        """Predict reaction outcome for unknown combinations"""
        # Generate molecular fingerprints
        fp1 = self._get_fingerprint(reactant1)
        fp2 = self._get_fingerprint(reactant2)
        
        # Combine fingerprints
        combined = np.concatenate([fp1, fp2])
        
        # Make prediction
        inputs = self.tokenizer(str(combined.tolist()), return_tensors="pt")
        outputs = self.model(**inputs)
        
        # Process prediction
        probabilities = torch.nn.functional.softmax(outputs.logits, dim=1)
        prediction = torch.argmax(probabilities, dim=1)
        
        return {
            'probability': probabilities[0][prediction].item(),
            'will_react': bool(prediction.item()),
            'confidence': 'High' if probabilities[0][prediction].item() > 0.8 else 'Medium'
        }
    
    def _get_fingerprint(self, smiles):
        """Generate Morgan fingerprint for molecule"""
        mol = Chem.MolFromSmiles(smiles)
        return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, self.fingerprint_size))
