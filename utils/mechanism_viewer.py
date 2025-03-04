from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import numpy as np

SPECIFIC_MECHANISMS = {
    ("HCl", "NaOH"): [
        {
            "name": "Initial State",
            "smiles": "[Na+].[OH-].[Cl-].[H+]",
            "description": "Dissociated ions in solution",
            "arrow_text": "Dissociation in water"
        },
        {
            "name": "Ion Combination",
            "smiles": "[Na+].[Cl-].O",
            "description": "H+ and OH- combine to form H2O",
            "arrow_text": "H+ + OH- → H2O"
        }
    ],
    ("H2SO4", "NaOH"): [
        {
            "name": "Initial State",
            "smiles": "[Na+].[OH-].O=S(=O)(O)O",
            "description": "Reactants in solution",
            "arrow_text": "Dissociation"
        },
        {
            "name": "Final State",
            "smiles": "[Na+].[Na+].[O-]S(=O)(=O)[O-]",
            "description": "Complete neutralization",
            "arrow_text": "Salt formation"
        }
    ],
    ("AgNO3", "NaCl"): [
        {
            "name": "Dissociation",
            "smiles": "[Ag+].[N+](=O)[O-].[Na+].[Cl-]",
            "description": "Ions in solution",
            "arrow_text": "Ion exchange"
        },
        {
            "name": "Precipitation",
            "smiles": "[Ag]Cl.[Na+].[N+](=O)[O-]",
            "description": "Silver chloride precipitate forms",
            "arrow_text": "Precipitation"
        }
    ],
    ("CaCO3", "HCl"): [
        {
            "name": "Initial Attack",
            "smiles": "O=C(O[Ca])O.[H+].[Cl-]",
            "description": "Acid attacks carbonate",
            "arrow_text": "Proton transfer"
        },
        {
            "name": "Decomposition",
            "smiles": "[Ca+2].[Cl-].[Cl-].O=C=O.O",
            "description": "CO2 evolution",
            "arrow_text": "Gas formation"
        }
    ],
    ("Fe", "HCl"): [
        {
            "name": "Initial State",
            "smiles": "[Fe].[H+].[Cl-]",
            "description": "Metal in acid",
            "arrow_text": "Oxidation"
        },
        {
            "name": "Product Formation",
            "smiles": "[Fe+2].[Cl-].[Cl-].[H][H]",
            "description": "H2 gas evolution",
            "arrow_text": "Reduction"
        }
    ]
}

def create_mechanism_image(step_data):
    """Create a detailed mechanism step image"""
    try:
        mol = Chem.MolFromSmiles(step_data["smiles"])
        if mol is None:
            return None
            
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Create smaller image with details
        legend = f"{step_data['name']}\n{step_data['description']}"
        img = Draw.MolToImage(mol, size=(250, 250), legend=legend)
        
        return img, step_data["description"], step_data["arrow_text"]
    except:
        return None

def draw_reaction_mechanism(reaction_type, reactants, products):
    """Create visual representation of reaction mechanism"""
    # Check for specific mechanisms first
    if len(reactants) == 2:
        reaction_key = tuple(sorted(reactants))
        if reaction_key in SPECIFIC_MECHANISMS:
            mechanism_steps = []
            for step in SPECIFIC_MECHANISMS[reaction_key]:
                result = create_mechanism_image(step)
                if result:
                    img, desc, arrow = result
                    mechanism_steps.append((step["name"], img, desc, arrow))
            return mechanism_steps
    
    # Fall back to generic mechanisms
    generic_mechanisms = {
        "Neutralization": [
            {
                "name": "Dissociation",
                "smiles": "[H+].[OH-]",
                "description": "Acid and base dissociation",
                "arrow_text": "Ionization"
            },
            {
                "name": "Water Formation",
                "smiles": "O",
                "description": "H+ and OH- combine",
                "arrow_text": "Neutralization"
            }
        ],
        # Add more generic mechanisms...
    }
    
    if reaction_type in generic_mechanisms:
        mechanism_steps = []
        for step in generic_mechanisms[reaction_type]:
            result = create_mechanism_image(step)
            if result:
                img, desc, arrow = result
                mechanism_steps.append((step["name"], img, desc, arrow))
        return mechanism_steps
    
    return None

def get_available_mechanisms():
    """Return list of reaction types with available mechanisms"""
    return [
        "Neutralization",
        "Oxidation",
        "Esterification",
        "Double Displacement",
        "Addition"
    ]

def get_reaction_type(chemical1, chemical2):
    """Determine reaction type from reactants"""
    if chemical1 in ["HCl", "H2SO4"] and chemical2 in ["NaOH", "KOH"]:
        return "Neutralization"
    elif "NO3" in chemical1 + chemical2 and "Cl" in chemical1 + chemical2:
        return "Double Displacement"
    elif chemical1 in ["Fe", "Zn"] and chemical2 in ["HCl", "H2SO4"]:
        return "Single Displacement"
    elif chemical1 == "CaCO3" and chemical2 in ["HCl", "H2SO4"]:
        return "Decomposition"
    return "Unknown"
