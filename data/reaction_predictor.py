import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

class ReactionPredictor:
    def predict_product(self, chemical1, chemical2, conditions):
        """Predict reaction products and outcomes"""
        try:
            # Try both combinations of the pair key
            pair_key1 = f"{chemical1['formula']}-{chemical2['formula']}"
            pair_key2 = f"{chemical2['formula']}-{chemical1['formula']}"
            
            reaction = None
            if pair_key1 in KNOWN_REACTIONS:
                reaction = KNOWN_REACTIONS[pair_key1]
            elif pair_key2 in KNOWN_REACTIONS:
                reaction = KNOWN_REACTIONS[pair_key2]
            
            if reaction:
                # Calculate yield based on conditions
                base_yield = reaction['base_yield']
                temp_factor = self._calculate_temperature_effect(conditions['temperature'], reaction['optimal_temp'])
                pressure_factor = self._calculate_pressure_effect(conditions['pressure'], reaction['optimal_pressure'])
                catalyst_factor = 1.2 if conditions['catalyst'] in reaction['suitable_catalysts'] else 0.8
                
                predicted_yield = min(99, base_yield * temp_factor * pressure_factor * catalyst_factor)
                
                return {
                    'success': True,
                    'products': reaction['products'],
                    'yield': predicted_yield,
                    'reaction_type': reaction['type'],
                    'mechanism': reaction['mechanism'],
                    'energy': reaction['energy_required'],
                    'time': reaction['estimated_time'],
                    'hazards': reaction['potential_hazards']
                }
            
            return {
                'success': False,
                'error': f"No known reaction between {chemical1['formula']} and {chemical2['formula']}"
            }
            
        except Exception as e:
            return {
                'success': False,
                'error': f"Error during prediction: {str(e)}"
            }
    
    def _calculate_temperature_effect(self, actual_temp, optimal_temp):
        """Calculate temperature effect on yield"""
        diff = abs(actual_temp - optimal_temp)
        return np.exp(-0.001 * diff)
    
    def _calculate_pressure_effect(self, actual_pressure, optimal_pressure):
        """Calculate pressure effect on yield"""
        ratio = actual_pressure / optimal_pressure
        return 2 / (1 + np.exp(-ratio)) - 0.5

KNOWN_REACTIONS = {
    "H2SO4-NaOH": {
        "products": ["Na2SO4", "H2O"],
        "base_yield": 95,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Neutralization",
        "mechanism": "Acid-base reaction",
        "energy_required": -57.7,
        "estimated_time": 0.1,
        "potential_hazards": ["Exothermic", "Corrosive"]
    },
    "NaOH-HCl": {
        "products": ["NaCl", "H2O"],
        "base_yield": 98,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Neutralization",
        "mechanism": "Acid-base reaction",
        "energy_required": -57.1,
        "estimated_time": 0.1,
        "potential_hazards": ["Exothermic"]
    },
    "CH3OH-H2SO4": {
        "products": ["CH3HSO4", "H2O"],
        "base_yield": 75,
        "optimal_temp": 20,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Esterification",
        "mechanism": "Nucleophilic substitution",
        "energy_required": -25.0,
        "estimated_time": 0.5,
        "potential_hazards": ["Corrosive", "Flammable"]
    },
    "NaCl-AgNO3": {
        "products": ["AgCl", "NaNO3"],
        "base_yield": 95,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Double Displacement",
        "mechanism": "Precipitation reaction",
        "energy_required": -12.3,
        "estimated_time": 0.1,
        "potential_hazards": ["Silver compounds", "Mild irritant"]
    },
    "CaCO3-HCl": {
        "products": ["CaCl2", "H2O", "CO2"],
        "base_yield": 92,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Decomposition",
        "mechanism": "Acid-base reaction",
        "energy_required": -28.4,
        "estimated_time": 0.2,
        "potential_hazards": ["CO2 evolution", "Corrosive acid"]
    },
    "H2-O2": {
        "products": ["H2O"],
        "base_yield": 99,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["Pt", "Pd"],
        "type": "Synthesis",
        "mechanism": "Oxidation-reduction",
        "energy_required": -285.8,
        "estimated_time": 0.1,
        "potential_hazards": ["Explosive", "Heat generation"]
    },
    "CH3COOH-CH3OH": {
        "products": ["CH3COOCH3", "H2O"],
        "base_yield": 67,
        "optimal_temp": 65,
        "optimal_pressure": 1,
        "suitable_catalysts": ["H2SO4", "H+"],
        "type": "Esterification",
        "mechanism": "Nucleophilic acyl substitution",
        "energy_required": 45.5,
        "estimated_time": 2.0,
        "potential_hazards": ["Flammable", "Irritant"]
    },
    "Fe-O2": {
        "products": ["Fe2O3"],
        "base_yield": 85,
        "optimal_temp": 500,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Oxidation",
        "mechanism": "Direct combination",
        "energy_required": -824.2,
        "estimated_time": 1.0,
        "potential_hazards": ["High temperature", "Metal oxide formation"]
    },
    "NH3-HCl": {
        "products": ["NH4Cl"],
        "base_yield": 98,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Acid-Base",
        "mechanism": "Direct combination",
        "energy_required": -176.2,
        "estimated_time": 0.1,
        "potential_hazards": ["Gas evolution", "Irritant"]
    },
    "Na2CO3-CaCl2": {
        "products": ["CaCO3", "2NaCl"],
        "base_yield": 90,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Double Displacement",
        "mechanism": "Precipitation",
        "energy_required": -32.1,
        "estimated_time": 0.2,
        "potential_hazards": ["Mild irritant"]
    },
    "CH4-Cl2": {
        "products": ["CH3Cl", "HCl"],
        "base_yield": 75,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["UV light"],
        "type": "Substitution",
        "mechanism": "Free radical halogenation",
        "energy_required": -103.4,
        "estimated_time": 0.5,
        "potential_hazards": ["Flammable", "Toxic gas"]
    },
    "KMnO4-H2O2": {
        "products": ["MnO2", "O2", "KOH", "H2O"],
        "base_yield": 88,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Redox",
        "mechanism": "Decomposition",
        "energy_required": -198.6,
        "estimated_time": 0.3,
        "potential_hazards": ["Strong oxidizer", "O2 evolution"]
    },
    "Zn-HCl": {
        "products": ["ZnCl2", "H2"],
        "base_yield": 95,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Single Displacement",
        "mechanism": "Redox reaction",
        "energy_required": -153.9,
        "estimated_time": 0.2,
        "potential_hazards": ["H2 evolution", "Corrosive"]
    },
    "KOH-HNO3": {
        "products": ["KNO3", "H2O"],
        "base_yield": 96,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Neutralization",
        "mechanism": "Acid-base reaction",
        "energy_required": -56.2,
        "estimated_time": 0.1,
        "potential_hazards": ["Exothermic", "Corrosive"]
    },
    "Fe2O3-CO": {
        "products": ["Fe", "CO2"],
        "base_yield": 88,
        "optimal_temp": 800,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Reduction",
        "mechanism": "Metal oxide reduction",
        "energy_required": 24.5,
        "estimated_time": 1.0,
        "potential_hazards": ["High temperature", "CO gas"]
    },
    "Cu-HNO3": {
        "products": ["Cu(NO3)2", "NO2", "H2O"],
        "base_yield": 92,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Oxidation",
        "mechanism": "Metal oxidation",
        "energy_required": -147.3,
        "estimated_time": 0.5,
        "potential_hazards": ["Toxic gas", "Corrosive"]
    },
    "Na2CO3-HCl": {
        "products": ["NaCl", "H2O", "CO2"],
        "base_yield": 95,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Neutralization",
        "mechanism": "Acid-base decomposition",
        "energy_required": -63.4,
        "estimated_time": 0.2,
        "potential_hazards": ["CO2 evolution"]
    },
    "BaCl2-Na2SO4": {
        "products": ["BaSO4", "NaCl"],
        "base_yield": 98,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Double Displacement",
        "mechanism": "Precipitation",
        "energy_required": -21.8,
        "estimated_time": 0.1,
        "potential_hazards": ["Barium compounds"]
    },
    "Al-Fe2O3": {
        "products": ["Al2O3", "Fe"],
        "base_yield": 92,
        "optimal_temp": 950,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Single Displacement",
        "mechanism": "Thermite reaction",
        "energy_required": -852.4,
        "estimated_time": 0.1,
        "potential_hazards": ["Extremely exothermic", "High temperature"]
    },
    "H2O2-KI": {
        "products": ["H2O", "O2", "KOH", "I2"],
        "base_yield": 85,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Decomposition",
        "mechanism": "Catalytic decomposition",
        "energy_required": -196.1,
        "estimated_time": 0.3,
        "potential_hazards": ["O2 evolution", "Oxidizer"]
    },
    "NaHCO3-CH3COOH": {
        "products": ["CH3COONa", "H2O", "CO2"],
        "base_yield": 94,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Acid-Base",
        "mechanism": "Neutralization with decomposition",
        "energy_required": -48.2,
        "estimated_time": 0.2,
        "potential_hazards": ["CO2 evolution"]
    },
    "K2Cr2O7-HCl": {
        "products": ["CrCl3", "Cl2", "KCl", "H2O"],
        "base_yield": 89,
        "optimal_temp": 60,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Redox",
        "mechanism": "Oxidation-reduction",
        "energy_required": -128.5,
        "estimated_time": 0.5,
        "potential_hazards": ["Toxic chemicals", "Oxidizer", "Chlorine gas"]
    },
    "Ca(OH)2-CO2": {
        "products": ["CaCO3", "H2O"],
        "base_yield": 96,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Precipitation",
        "mechanism": "Gas absorption",
        "energy_required": -113.0,
        "estimated_time": 0.3,
        "potential_hazards": ["Basic solution"]
    },
    "NH4NO3-Heat": {
        "products": ["N2O", "H2O"],
        "base_yield": 99,
        "optimal_temp": 200,
        "optimal_pressure": 1,
        "suitable_catalysts": ["None"],
        "type": "Decomposition",
        "mechanism": "Thermal decomposition",
        "energy_required": 36.4,
        "estimated_time": 0.2,
        "potential_hazards": ["Explosive potential", "Toxic gas"]
    },
    "C2H5OH-O2": {
        "products": ["CO2", "H2O"],
        "base_yield": 99,
        "optimal_temp": 25,
        "optimal_pressure": 1,
        "suitable_catalysts": ["Pt"],
        "type": "Combustion",
        "mechanism": "Complete oxidation",
        "energy_required": -1367.0,
        "estimated_time": 0.1,
        "potential_hazards": ["Highly flammable", "Fire hazard"]
    }
}
