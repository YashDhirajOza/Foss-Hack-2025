"""Module for thermodynamic calculations in chemical reactions."""

import numpy as np
from typing import Dict, Union, List

class ThermodynamicsCalculator:
    """Calculator for thermodynamic properties of chemical reactions."""
    
    def __init__(self):
        self.gas_constant = 8.314  # J/(mol·K)
    
    def calculate_gibbs_energy(self, delta_h: float, delta_s: float, 
                             temperature: float) -> float:
        """Calculate Gibbs free energy.
        
        Args:
            delta_h: Enthalpy change (J/mol)
            delta_s: Entropy change (J/mol·K)
            temperature: Temperature (K)
            
        Returns:
            Gibbs free energy change (J/mol)
        """
        return delta_h - temperature * delta_s
    
    def calculate_equilibrium_constant(self, delta_g: float, 
                                    temperature: float) -> float:
        """Calculate equilibrium constant.
        
        Args:
            delta_g: Gibbs free energy change (J/mol)
            temperature: Temperature (K)
            
        Returns:
            Equilibrium constant
        """
        return np.exp(-delta_g / (self.gas_constant * temperature))
    
    def calculate_reaction_quotient(self, concentrations: Dict[str, float], 
                                  stoichiometry: Dict[str, int]) -> float:
        """Calculate reaction quotient.
        
        Args:
            concentrations: Dict of species concentrations
            stoichiometry: Dict of stoichiometric coefficients
            
        Returns:
            Reaction quotient
        """
        quotient = 1.0
        for species, coeff in stoichiometry.items():
            if species in concentrations:
                quotient *= concentrations[species] ** coeff
        return quotient
    
    def predict_spontaneity(self, delta_g: float) -> str:
        """Determine reaction spontaneity.
        
        Args:
            delta_g: Gibbs free energy change (J/mol)
            
        Returns:
            Description of spontaneity
        """
        if delta_g < 0:
            return "Spontaneous"
        if delta_g > 0:
            return "Non-spontaneous"
        return "At equilibrium"
