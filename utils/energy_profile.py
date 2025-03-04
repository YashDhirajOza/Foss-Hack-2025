import plotly.graph_objects as go
import numpy as np

def create_energy_profile(reaction_data, conditions):
    """Generate interactive energy profile diagram"""
    fig = go.Figure()
    
    # Calculate energy levels
    initial_energy = 0
    transition_energy = reaction_data['energy'] + 20 if reaction_data['energy'] > 0 else 20
    final_energy = reaction_data['energy']
    
    # Temperature effect
    temp_factor = conditions['temperature'] / 298.15
    transition_energy *= temp_factor
    
    # Create smooth curve
    x = np.linspace(0, 100, 200)
    y = np.concatenate([
        50 * np.exp(-(x[:66] - 33)**2 / 400),
        50 * np.exp(-(x[66:] - 66)**2 / 400)
    ])
    
    # Add traces
    fig.add_trace(go.Scatter(x=x, y=y * temp_factor + initial_energy,
                            mode='lines', name='Energy Path'))
    
    # Add key points
    fig.add_trace(go.Scatter(
        x=[10, 50, 90],
        y=[initial_energy, transition_energy, final_energy],
        mode='markers+text',
        name='Energy States',
        text=['Reactants', 'Transition State', 'Products'],
        textposition='top center'
    ))
    
    fig.update_layout(
        title='Reaction Energy Profile',
        xaxis_title='Reaction Progress',
        yaxis_title='Energy (kJ/mol)'
    )
    
    return fig
