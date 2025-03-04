SMILES_MAPPING = {
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

CHEMICALS_DATABASE = {
    "H2SO4": {
        "name": "Sulfuric Acid",
        "formula": "H2SO4",
        "molar_mass": "98.079 g/mol",
        "density": "1.83 g/cm³",
        "melting_point": "10°C",
        "boiling_point": "337°C",
        "description": "Strong mineral acid, highly corrosive, dehydrating agent",
        "uses": ["Battery acid", "Fertilizer production", "Mineral processing"],
        "hazards": "Corrosive! Causes severe skin burns and eye damage.",
        "precautions": "Wear PPE (gloves, goggles, lab coat). Use in fume hood.",
        "first_aid": "Flush with water for 15 mins. Neutralize with sodium bicarbonate.",
        "storage": "Store in acid cabinet. Keep away from metals.",
        "incompatible": "Bases, metals, organic materials",
        "smiles": "OS(=O)(=O)O"
    },
    "CH3OH": {
        "name": "Methanol",
        "formula": "CH3OH",
        "molar_mass": "32.04 g/mol",
        "density": "0.792 g/cm³",
        "melting_point": "-97.6°C",
        "boiling_point": "64.7°C",
        "description": "Simplest alcohol, colorless liquid, used as solvent",
        "uses": ["Solvent", "Fuel", "Antifreeze"],
        "hazards": "Highly flammable! Toxic if swallowed.",
        "precautions": "No flames/sparks. Use ventilation.",
        "first_aid": "Medical attention if swallowed.",
        "storage": "Flammables cabinet.",
        "incompatible": "Oxidizers",
        "smiles": "CO"
    },
    "NaOH": {
        "name": "Sodium Hydroxide",
        "formula": "NaOH",
        "molar_mass": "40.00 g/mol",
        "density": "2.13 g/cm³",
        "melting_point": "318°C",
        "boiling_point": "1388°C",
        "description": "Strong base, commonly known as caustic soda",
        "uses": ["Soap making", "pH control", "Chemical manufacturing"],
        "hazards": "Highly corrosive to skin and eyes",
        "precautions": "Use proper PPE, avoid contact with acids",
        "first_aid": "Flush with water for 15 minutes, seek medical attention",
        "storage": "Keep in dry place, sealed container",
        "incompatible": "Acids, metals, organic materials",
        "smiles": "[Na+].[OH-]"
    },
    "HCl": {
        "name": "Hydrochloric Acid",
        "formula": "HCl",
        "molar_mass": "36.46 g/mol",
        "density": "1.19 g/cm³",
        "melting_point": "-27°C",
        "boiling_point": "48°C",
        "description": "Strong mineral acid",
        "uses": ["pH adjustment", "Metal cleaning", "Chemical synthesis"],
        "hazards": "Corrosive, respiratory irritant",
        "precautions": "Use in fume hood, wear acid-resistant PPE",
        "first_aid": "Flush with water, seek immediate medical attention",
        "storage": "Acid cabinet, ventilated area",
        "incompatible": "Bases, oxidizers, metals",
        "smiles": "[H+].[Cl-]"
    },
    "AgNO3": {
        "name": "Silver Nitrate",
        "formula": "AgNO3",
        "molar_mass": "169.87 g/mol",
        "density": "4.35 g/cm³",
        "melting_point": "210°C",
        "boiling_point": "440°C",
        "description": "Inorganic compound used in silver reactions",
        "uses": ["Photography", "Mirror making", "Chemical testing"],
        "hazards": "Oxidizer, corrosive, stains skin",
        "precautions": "Avoid organic materials, use PPE",
        "first_aid": "Flush with water, remove contaminated clothing",
        "storage": "Dark bottle, away from light",
        "incompatible": "Organic materials, reducing agents"
    },
    "CaCO3": {
        "name": "Calcium Carbonate",
        "formula": "CaCO3",
        "molar_mass": "100.09 g/mol",
        "density": "2.71 g/cm³",
        "melting_point": "825°C",
        "boiling_point": "Decomposes",
        "description": "Common mineral, found in limestone",
        "uses": ["Construction", "Antacid", "Calcium supplement"],
        "hazards": "Low toxicity",
        "precautions": "Dust mask when handling powder",
        "first_aid": "Remove from exposure, fresh air",
        "storage": "Keep dry",
        "incompatible": "Strong acids"
    },
    "NH3": {
        "name": "Ammonia",
        "formula": "NH3",
        "molar_mass": "17.03 g/mol",
        "density": "0.73 g/cm³ (liquid)",
        "melting_point": "-77.7°C",
        "boiling_point": "-33.3°C",
        "description": "Colorless gas with pungent odor",
        "uses": ["Fertilizer production", "Cleaning agent", "Refrigerant"],
        "hazards": "Toxic gas, corrosive to respiratory system",
        "precautions": "Use in well-ventilated area, gas mask",
        "first_aid": "Fresh air, respiratory support if needed",
        "storage": "Pressurized container, cool area",
        "incompatible": "Acids, oxidizers, halogens"
    },
    "O2": {
        "name": "Oxygen",
        "formula": "O2",
        "molar_mass": "32.00 g/mol",
        "density": "1.429 g/L (gas)",
        "melting_point": "-218.79°C",
        "boiling_point": "-182.96°C",
        "description": "Colorless gas essential for life",
        "uses": ["Medical treatment", "Combustion", "Industrial processes"],
        "hazards": "Strong oxidizer, supports combustion",
        "precautions": "Keep away from oils and grease",
        "first_aid": "No direct hazard",
        "storage": "Pressurized cylinder, secure",
        "incompatible": "Flammable materials, oils"
    },
    "Fe": {
        "name": "Iron",
        "formula": "Fe",
        "molar_mass": "55.845 g/mol",
        "density": "7.874 g/cm³",
        "melting_point": "1538°C",
        "boiling_point": "2862°C",
        "description": "Metallic element, magnetic",
        "uses": ["Steel production", "Construction", "Manufacturing"],
        "hazards": "Combustible as powder",
        "precautions": "Avoid strong oxidizers",
        "first_aid": "Remove from exposure",
        "storage": "Keep dry to prevent rust",
        "incompatible": "Strong acids, oxidizers"
    },
    "KMnO4": {
        "name": "Potassium Permanganate",
        "formula": "KMnO4",
        "molar_mass": "158.034 g/mol",
        "density": "2.703 g/cm³",
        "melting_point": "240°C",
        "boiling_point": "Decomposes",
        "description": "Purple crystalline solid, strong oxidizer",
        "uses": ["Water treatment", "Disinfectant", "Chemical synthesis"],
        "hazards": "Strong oxidizer, environmental hazard",
        "precautions": "Avoid organic materials, wear PPE",
        "first_aid": "Flush with water, seek medical attention",
        "storage": "Cool, dry place, separate from reducers",
        "incompatible": "Organic materials, reducing agents"
    },
    "H2O2": {
        "name": "Hydrogen Peroxide",
        "formula": "H2O2",
        "molar_mass": "34.01 g/mol",
        "density": "1.45 g/cm³",
        "melting_point": "-0.43°C",
        "boiling_point": "150.2°C",
        "description": "Clear liquid, strong oxidizer",
        "uses": ["Disinfectant", "Bleaching agent", "Rocket propellant"],
        "hazards": "Strong oxidizer, corrosive to skin",
        "precautions": "Use PPE, avoid organic materials",
        "first_aid": "Flush with water, seek medical attention",
        "storage": "Cool, dark place, vented cap",
        "incompatible": "Organic materials, metals",
        "smiles": "OO"
    },
    "CH3COOH": {
        "name": "Acetic Acid",
        "formula": "CH3COOH",
        "molar_mass": "60.05 g/mol",
        "density": "1.049 g/cm³",
        "melting_point": "16.6°C",
        "boiling_point": "118.1°C",
        "description": "Clear liquid with pungent odor",
        "uses": ["Food additive", "Chemical synthesis", "Cleaning agent"],
        "hazards": "Corrosive, flammable liquid",
        "precautions": "Use in ventilated area, wear PPE",
        "first_aid": "Flush with water, neutralize with base",
        "storage": "Sealed container, away from bases",
        "incompatible": "Bases, oxidizers",
        "smiles": "CC(=O)O"
    },
    "K2Cr2O7": {
        "name": "Potassium Dichromate",
        "formula": "K2Cr2O7",
        "molar_mass": "294.18 g/mol",
        "density": "2.676 g/cm³",
        "melting_point": "398°C",
        "boiling_point": "500°C (decomposes)",
        "description": "Orange-red crystals, strong oxidizer",
        "uses": ["Oxidizing agent", "Chrome plating", "Laboratory reagent"],
        "hazards": "Toxic, carcinogenic, strong oxidizer",
        "precautions": "Full PPE required, avoid skin contact",
        "first_aid": "Immediate medical attention required",
        "storage": "Sealed container, away from reducers",
        "incompatible": "Reducing agents, organic materials",
        "smiles": "[K+].[K+].[O-]Cr(=O)(=O)OCr(=O)(=O)[O-]"
    },
    "ZnCl2": {
        "name": "Zinc Chloride",
        "formula": "ZnCl2",
        "molar_mass": "136.29 g/mol",
        "density": "2.91 g/cm³",
        "melting_point": "290°C",
        "boiling_point": "732°C",
        "description": "White crystalline solid, hygroscopic",
        "uses": ["Metal treatment", "Flux agent", "Dehydrating agent"],
        "hazards": "Corrosive, hygroscopic",
        "precautions": "Avoid moisture, use PPE",
        "first_aid": "Flush with water, seek medical attention",
        "storage": "Airtight container, dry area",
        "incompatible": "Water sensitive",
        "smiles": "[Zn+2].([Cl-]).[Cl-]"
    },
    "Na2CO3": {
        "name": "Sodium Carbonate",
        "formula": "Na2CO3",
        "molar_mass": "105.99 g/mol",
        "density": "2.54 g/cm³",
        "melting_point": "851°C",
        "boiling_point": "1600°C",
        "description": "White powder, basic salt",
        "uses": ["Water treatment", "Glass making", "Detergent"],
        "hazards": "Eye and respiratory irritant",
        "precautions": "Avoid dust formation",
        "first_aid": "Flush with water",
        "storage": "Sealed container, dry area",
        "incompatible": "Acids",
        "smiles": "[Na+].[Na+].[O-]C(=O)[O-]"
    },
    "KCl": {
        "name": "Potassium Chloride",
        "formula": "KCl",
        "molar_mass": "74.55 g/mol",
        "density": "1.98 g/cm³",
        "melting_point": "770°C",
        "boiling_point": "1500°C",
        "description": "White crystalline solid",
        "uses": ["Fertilizer", "Medicine", "Food additive"],
        "hazards": "Low toxicity",
        "precautions": "Standard lab practices",
        "first_aid": "Flush with water",
        "storage": "Standard storage",
        "incompatible": "Strong oxidizers",
        "smiles": "[K+].[Cl-]"
    },
    "AlCl3": {
        "name": "Aluminum Chloride",
        "formula": "AlCl3",
        "molar_mass": "133.34 g/mol",
        "density": "2.48 g/cm³",
        "melting_point": "190°C",
        "boiling_point": "182.7°C (sublimes)",
        "description": "White or pale yellow solid",
        "uses": ["Catalyst", "Chemical synthesis", "Dye production"],
        "hazards": "Corrosive, water-reactive",
        "precautions": "Keep dry, use in fume hood",
        "first_aid": "Flush with water, seek medical attention",
        "storage": "Airtight container, dry area",
        "incompatible": "Water, alcohols",
        "smiles": "[Al+3].[Cl-].[Cl-].[Cl-]"
    }
}

for chemical, data in CHEMICALS_DATABASE.items():
    if chemical in SMILES_MAPPING:
        data['smiles'] = SMILES_MAPPING[chemical]

REACTION_TEMPLATES = {
    "Hydrogenation": {
        "reactants": ["C=C", "H2"],
        "products": ["CH3-CH3"],
        "equation": "C=C + H2 -> CH3-CH3",
        "conditions": {
            "temperature": "25-50°C",
            "pressure": "1-10 atm",
            "catalyst": "Pd/C",
            "solvent": "Ethanol"
        },
        "type": "Addition",
        "mechanism": "Syn addition of H2 across double bond",
        "yield": "90-95%"
    },
    "Esterification": {
        "reactants": ["CH3COOH", "CH3OH"],
        "products": ["CH3COOCH3", "H2O"],
        "equation": "CH3COOH + CH3OH -> CH3COOCH3 + H2O",
        "conditions": {
            "temperature": "60-70°C",
            "catalyst": "H2SO4",
            "time": "2-4 hours"
        },
        "type": "Condensation",
        "mechanism": "Nucleophilic acyl substitution",
        "yield": "75-85%"
    },
    "Grignard Reaction": {
        "reactants": ["R-X", "Mg", "Ether"],
        "products": ["R-MgX"],
        "equation": "R-X + Mg -> R-MgX",
        "conditions": {
            "temperature": "-10 to 25°C",
            "solvent": "Dry ether",
            "workup": "H3O+"
        },
        "type": "Addition",
        "mechanism": "Nucleophilic addition of alkyl group",
        "yield": "80-95%"
    },
}

def format_reaction(template_name):
    """Format reaction for display"""
    rxn = REACTION_TEMPLATES.get(template_name, {})
    return {
        "display": f"{' + '.join(rxn.get('reactants', []))} → {' + '.join(rxn.get('products', []))}",
        "details": rxn
    }

def get_smiles(chemical_formula):
    """Get correct SMILES notation for a chemical formula"""
    return SMILES_MAPPING.get(chemical_formula, '')
