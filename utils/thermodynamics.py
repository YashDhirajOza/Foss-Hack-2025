import numpy as np

class ThermodynamicsCalculator:
    def __init__(self):
        self.R = 8.314  # Gas constant J/(mol·K)
        
    def calculate_gibbs_energy(self, delta_h, delta_s, temperature):
        """Calculate Gibbs free energy"""
        return delta_h - temperature * delta_s
    
    def calculate_equilibrium_constant(self, delta_g, temperature):
        """Calculate equilibrium constant"""
        return np.exp(-delta_g / (self.R * temperature))
    
    def calculate_reaction_quotient(self, concentrations, stoichiometry):
        """Calculate reaction quotient"""
        Q = 1
        for species, coeff in stoichiometry.items():
            if species in concentrations:
                Q *= concentrations[species] ** coeff
        return Q
    
    def predict_spontaneity(self, delta_g):
        """Determine reaction spontaneity"""
        if delta_g < 0:
            return "Spontaneous"
        elif delta_g > 0:
            return "Non-spontaneous"
        return "At equilibrium"
