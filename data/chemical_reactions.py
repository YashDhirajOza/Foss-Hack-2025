INDUSTRIAL_REACTIONS = {
    "Haber Process": {
        "simple_equation": "N₂ + 3H₂ → 2NH₃",
        "reactant_a": "N₂",
        "reactant_b": "3H₂",
        "product": "2NH₃",
        "reactants": ["N2", "H2"],
        "products": ["NH3"],
        "equation": "N2 + 3H2 ⇌ 2NH3",
        "industrial_scale": {
            "annual_production": "170 million tonnes",
            "major_producers": ["BASF", "CF Industries", "Yara"],
            "costs": {"capital": "High", "operational": "Medium-High"}
        },
        "conditions": {
            "temperature": "450-500°C",
            "pressure": "150-300 atm",
            "catalyst": "Iron with K2O, Al2O3 promoters",
            "reaction_time": "Continuous flow",
            "yield": "15% per pass, >98% with recycling"
        },
        "economics": {
            "raw_material_cost": "$200-400/ton NH3",
            "energy_consumption": "30-35 GJ/ton NH3",
            "market_price": "$600-1000/ton NH3"
        },
        "process_controls": {
            "critical_parameters": ["Pressure", "Temperature", "H2:N2 ratio"],
            "monitoring": ["Gas composition", "Catalyst activity", "Heat exchange"],
            "safety_measures": ["Pressure relief", "Gas detection", "Emergency shutdown"]
        },
        "applications": ["Fertilizers", "Explosives", "Plastics"],
        "environmental_impact": {
            "CO2_emissions": "1.5-3.0 tons/ton NH3",
            "energy_efficiency": "Medium",
            "sustainable_alternatives": "Green hydrogen-based production"
        },
        "optimization_metrics": {
            "conversion_vs_temperature": {
                "temps": [350, 400, 450, 500, 550],
                "conversion": [10, 12, 15, 13, 11]
            },
            "pressure_effect": {
                "pressures": [100, 150, 200, 250, 300],
                "yield": [8, 12, 15, 17, 18]
            }
        }
    },
    "Sulfuric Acid": {
        "simple_equation": "SO₃ + H₂O → H₂SO₄",
        "reactant_a": "SO₃",
        "reactant_b": "H₂O",
        "product": "H₂SO₄",
        "conditions": {
            "temperature": "80-100°C",
            "pressure": "1-2 atm",
            "catalyst": "None required",
            "yield": "98-99%"
        },

    },

}

REACTION_CATEGORIES = {
    "Petrochemical": ["Cracking", "Reforming", "Alkylation"],
    "Polymer": ["Polymerization", "Condensation", "Addition"],
    "Pharmaceutical": ["Coupling", "Hydrogenation", "Protection"],
    "Agricultural": ["Nitration", "Amination", "Oxidation"],
    "Fine Chemicals": ["Esterification", "Acylation", "Sulfonation"]
}

def get_reaction_metrics(reaction_name):
    """Calculate economic and process metrics"""
    reaction = INDUSTRIAL_REACTIONS[reaction_name]
    return {
        "roi": calculate_roi(reaction),
        "efficiency": calculate_efficiency(reaction),
        "sustainability": calculate_sustainability_score(reaction)
    }

def calculate_roi(reaction):
    """Calculate return on investment"""
    try:
        # Extract costs and revenue data
        raw_material_cost = float(reaction['economics']['raw_material_cost'].split('/')[0].replace('$','').split('-')[0])
        energy_cost = float(reaction['economics']['energy_consumption'].split()[0]) * 0.1  # Assuming $0.1 per GJ
        revenue = float(reaction['economics']['market_price'].split('/')[0].replace('$','').split('-')[0])
        
        # Calculate production metrics
        annual_production = float(reaction['industrial_scale']['annual_production'].split()[0])
        yield_percentage = float(reaction['conditions']['yield'].split('%')[0].split(',')[0])
        
        # Calculate costs and revenue
        total_cost = (raw_material_cost + energy_cost) * annual_production
        total_revenue = revenue * annual_production * (yield_percentage/100)
        
        # Calculate ROI
        roi = ((total_revenue - total_cost) / total_cost) * 100
        return round(roi, 2)
    except:
        return 0

def calculate_efficiency(reaction):
    """Calculate process efficiency"""
    try:

        energy_input = float(reaction['economics']['energy_consumption'].split()[0])
        theoretical_energy = energy_input * 0.6  # Theoretical minimum energy requirement
        energy_efficiency = (theoretical_energy / energy_input) * 100
        
        actual_yield = float(reaction['conditions']['yield'].split('%')[0].split(',')[0])
        yield_efficiency = actual_yield
        

        resource_efficiency = 90 if reaction['environmental_impact']['energy_efficiency'] == "High" else \
                            70 if reaction['environmental_impact']['energy_efficiency'] == "Medium" else 50
        
        # Calculate overall efficiency
        overall_efficiency = (energy_efficiency + yield_efficiency + resource_efficiency) / 3
        
        return {
            'overall': round(overall_efficiency, 2),
            'energy': round(energy_efficiency, 2),
            'yield': round(yield_efficiency, 2),
            'resource': round(resource_efficiency, 2)
        }
    except:
        return {
            'overall': 0,
            'energy': 0,
            'yield': 0,
            'resource': 0
        }

def calculate_sustainability_score(reaction):
    """Calculate sustainability metrics"""
    try:
        # Initialize base score
        score = 100
        
        # CO2 emissions impact (-10 points per ton of CO2 per ton of product)
        co2_emissions = float(reaction['environmental_impact']['CO2_emissions'].split()[0].split('-')[0])
        score -= co2_emissions * 10
        
        # Energy efficiency impact
        if reaction['environmental_impact']['energy_efficiency'] == "Low":
            score -= 30
        elif reaction['environmental_impact']['energy_efficiency'] == "Medium":
            score -= 15
        
        # Resource utilization
        if 'recycling' in reaction['conditions']['yield'].lower():
            score += 20
        
        # Process safety impact
        hazard_count = len(reaction['process_controls']['safety_measures'])
        score -= hazard_count * 5
        
        # Sustainable alternatives bonus
        if 'sustainable_alternatives' in reaction['environmental_impact']:
            score += 10
            
        # Normalize score between 0 and 100
        final_score = max(0, min(100, score))
        
        return {
            'score': round(final_score, 2),
            'metrics': {
                'emissions_impact': round(co2_emissions * 10, 2),
                'energy_rating': reaction['environmental_impact']['energy_efficiency'],
                'safety_rating': 100 - (hazard_count * 5),
                'sustainability_potential': 'High' if final_score > 75 else 'Medium' if final_score > 50 else 'Low'
            }
        }
    except:
        return {
            'score': 0,
            'metrics': {
                'emissions_impact': 0,
                'energy_rating': 'Unknown',
                'safety_rating': 0,
                'sustainability_potential': 'Unknown'
            }
        }

def format_industrial_equation(reaction_data):
    """Format reaction equation in A + B → C format"""
    return f"{reaction_data['reactant_a']} + {reaction_data['reactant_b']} → {reaction_data['product']}"
