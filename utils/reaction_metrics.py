def calculate_atom_economy(reactant1, reactant2, products):
    """Calculate atom economy and other green chemistry metrics"""
    try:
        # Calculate molecular weights
        reactant1_mw = float(reactant1['molar_mass'].split()[0])
        reactant2_mw = float(reactant2['molar_mass'].split()[0])
        products_mw = sum(float(p.split()[0]) for p in [prod['molar_mass'] for prod in products])
        
        # Calculate metrics
        atom_economy = (products_mw / (reactant1_mw + reactant2_mw)) * 100
        
        return {
            'atom_economy': round(atom_economy, 2),
            'waste_factor': round(100 - atom_economy, 2),
            'efficiency_class': 'High' if atom_economy > 80 else 'Medium' if atom_economy > 60 else 'Low'
        }
    except:
        return None
