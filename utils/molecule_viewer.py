from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol
import streamlit.components.v1 as components

def convert_chemical_to_smiles(identifier):
    smiles_map = {
        'H2SO4': 'OS(=O)(=O)O',
        'CH3OH': 'CO',
        'CH3COOH': 'CC(=O)O',
        'H2O': 'O',
        'HCl': '[H+].[Cl-]',
        'NaOH': '[Na+].[OH-]',
        'NH3': 'N',
        'O2': 'O=O',
        'N2': 'N#N',
        'CO2': 'O=C=O',
        'CaCO3': '[Ca+2].[C+0]([O-])([O-])=O',
        'AlCl3': '[Al+3].[Cl-].[Cl-].[Cl-]',
        'KCl': '[K+].[Cl-]',
        'H2O2': 'OO',
        'AgNO3': '[Ag+].[N+](=O)([O-])[O-]',
        'K2Cr2O7': '[K+].[K+].[O-]Cr(=O)(=O)OCr(=O)(=O)[O-]',
        'ZnCl2': '[Zn+2].[Cl-].[Cl-]',
        'Na2CO3': '[Na+].[Na+].[O-]C(=O)[O-]'
    }
    return smiles_map.get(identifier, identifier)

def generate_3d_structure(identifier):
    try:
        smiles = convert_chemical_to_smiles(identifier)
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            mol = Chem.MolFromSmiles(identifier)
            
        if mol is None:
            return None
            
        mol = Chem.AddHs(mol)
        success = AllChem.EmbedMolecule(mol, randomSeed=42)
        
        if success == -1:
            AllChem.Compute2DCoords(mol)
            success = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
        
        if success != -1:
            AllChem.MMFFOptimizeMolecule(mol)
            
        return mol
    except Exception as e:
        return None

def show_molecule_3d(mol, size=(300, 300), style="stick"):
    try:
        if mol is None:
            return None
            
        viewer = py3Dmol.view(width=size[0], height=size[1])
        
        try:
            pdb = Chem.MolToPDBBlock(mol)
            viewer.addModel(pdb, "pdb")
        except:
            mol_block = Chem.MolToMolBlock(mol)
            viewer.addModel(mol_block, "mol")
        
        try:
            viewer.setStyle({'stick':{}})
            viewer.zoomTo()
        except:
            viewer.setStyle({'sphere':{}})
        
        return viewer
    except Exception as e:
        return None

def get_2d_depiction(mol):
    if mol is None:
        return None
    return Draw.MolToImage(mol)
