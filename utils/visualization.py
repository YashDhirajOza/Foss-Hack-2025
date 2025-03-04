import plotly.graph_objects as go
import numpy as np

def create_energy_diagram(reaction_data):
    """Create reaction energy diagram"""
    fig = go.Figure()
    
    # Energy levels
    energy_levels = {
        'reactants': 0,
        'transition_state': reaction_data['energy'] + 20 if reaction_data['energy'] > 0 else 20,
        'products': reaction_data['energy']
    }
    
    # Add energy curve
    x = np.linspace(0, 100, 100)
    y = np.concatenate([
        50 * np.exp(-(x[:33] - 20)**2 / 400),
        50 * np.exp(-(x[33:66] - 50)**2 / 400),
        50 * np.exp(-(x[66:] - 80)**2 / 400)
    ])
    
    fig.add_trace(go.Scatter(x=x, y=y + energy_levels['reactants'],
                            mode='lines', name='Energy Path'))
    
    # Add annotations
    fig.add_annotation(x=10, y=energy_levels['reactants'],
                      text='Reactants', showarrow=False)
    fig.add_annotation(x=50, y=energy_levels['transition_state'],
                      text='Transition State', showarrow=False)
    fig.add_annotation(x=90, y=energy_levels['products'],
                      text='Products', showarrow=False)
    
    fig.update_layout(
        title='Reaction Energy Diagram',
        xaxis_title='Reaction Progress',
        yaxis_title='Energy (kJ/mol)'
    )
    
    return fig

def create_kinetics_plot(time, conversion, temperature):
    """Create reaction kinetics visualization"""
    fig = go.Figure()
    
    # Rate constant (Arrhenius equation)
    k = 0.01 * np.exp(-50/(0.008314 * (temperature + 273.15)))
    
    # Theoretical curve
    t = np.linspace(0, max(time), 100)
    theoretical = 100 * (1 - np.exp(-k * t))
    
    fig.add_trace(go.Scatter(x=time, y=conversion,
                            mode='markers', name='Experimental'))
    fig.add_trace(go.Scatter(x=t, y=theoretical,
                            mode='lines', name='Theoretical'))
    
    fig.update_layout(
        title='Reaction Kinetics',
        xaxis_title='Time (min)',
        yaxis_title='Conversion (%)'
    )
    
    return fig
