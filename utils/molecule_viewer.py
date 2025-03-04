"""Module for handling molecular structure visualization using RDKit and Py3DMol."""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol

SMILES_MAP = {
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

def convert_chemical_to_smiles(identifier: str) -> str:
    """Convert chemical formula or name to SMILES notation.
    
    Args:
        identifier: Chemical formula or name
        
    Returns:
        SMILES string representation
    """
    return SMILES_MAP.get(identifier, identifier)

def generate_3d_structure(identifier: str) -> Chem.Mol:
    """Generate 3D molecular structure from identifier.
    
    Args:
        identifier: Chemical formula or SMILES string
        
    Returns:
        RDKit molecule object with 3D coordinates
    """
    try:
        smiles = convert_chemical_to_smiles(identifier)
        mol = Chem.MolFromSmiles(smiles)
        
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
    except Exception:
        return None

def show_molecule_3d(mol: Chem.Mol, size: tuple = (300, 300)) -> py3Dmol.view:
    """Display molecule in 3D using py3Dmol.
    
    Args:
        mol: RDKit molecule object
        size: Tuple of (width, height) for viewer
        
    Returns:
        py3Dmol viewer object
    """
    if mol is None:
        return None
        
    try:
        viewer = py3Dmol.view(width=size[0], height=size[1])
        
        try:
            mol_block = Chem.MolToMolBlock(mol)
            viewer.addModel(mol_block, "mol")
            viewer.setStyle({'stick':{}})
            viewer.zoomTo()
        except Exception:
            return None
        
        return viewer
    except Exception:
        return None

def get_2d_depiction(mol: Chem.Mol) -> AllChem.Image:
    """Get 2D depiction of molecule.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        PIL Image object
    """
    if mol is None:
        return None
    return Draw.MolToImage(mol)
