from rdkit import Chem
from rdkit.Chem import Draw
import streamlit as st

def draw_2d_structure(smiles, size=(300, 300), title=""):
    """Create 2D chemical structure visualization"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol, size=size, legend=title)
            return img
    except:
        return None

def draw_reaction_equation(reactant1_smiles, reactant2_smiles, product_smiles):
    """Draw complete reaction equation with structures"""
    try:
        # Create individual molecule images
        r1 = Chem.MolFromSmiles(reactant1_smiles)
        r2 = Chem.MolFromSmiles(reactant2_smiles)
        p = Chem.MolFromSmiles(product_smiles)
        
        # Combine into reaction scheme
        img = Draw.ReactionToImage([r1, r2], [p])
        return img
    except:
        return None
