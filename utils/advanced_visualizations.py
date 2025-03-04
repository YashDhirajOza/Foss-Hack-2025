import plotly.graph_objects as go
import plotly.express as px
import numpy as np

def create_3d_energy_surface(temperature_range, pressure_range, reaction_data):
    """Create 3D surface plot of energy landscape"""
    temp_points = np.linspace(temperature_range[0], temperature_range[1], 50)
    press_points = np.linspace(pressure_range[0], pressure_range[1], 50)
    
    T, P = np.meshgrid(temp_points, press_points)
    # Calculate Gibbs energy surface
    G = -8.314 * T * np.log(P) + reaction_data['energy']
    
    fig = go.Figure(data=[go.Surface(z=G, x=T, y=P)])
    fig.update_layout(
        title='Reaction Energy Landscape',
        scene = dict(
            xaxis_title='Temperature (K)',
            yaxis_title='Pressure (atm)',
            zaxis_title='Gibbs Energy (kJ/mol)'
        )
    )
    return fig

def create_reaction_network(reactants, intermediates, products):
    """Create reaction network visualization"""
    fig = go.Figure()
    
    # Create nodes for species
    nodes_x = [0, 1, 2]  # Positions for reactants, intermediates, products
    nodes_y = [0, 0, 0]
    
    # Add nodes
    fig.add_trace(go.Scatter(
        x=[0], y=[0],
        mode='markers+text',
        name='Reactants',
        text=reactants,
        marker=dict(size=20, color='blue')
    ))
    
    fig.add_trace(go.Scatter(
        x=[1], y=[0],
        mode='markers+text',
        name='Intermediates',
        text=intermediates,
        marker=dict(size=20, color='green')
    ))
    
    fig.add_trace(go.Scatter(
        x=[2], y=[0],
        mode='markers+text',
        name='Products',
        text=products,
        marker=dict(size=20, color='red')
    ))
    
    # Add arrows
    for i in range(len(nodes_x)-1):
        fig.add_annotation(
            x=nodes_x[i], y=nodes_y[i],
            ax=nodes_x[i+1], ay=nodes_y[i+1],
            xref='x', yref='y',
            axref='x', ayref='y',
            showarrow=True,
            arrowhead=2,
            arrowsize=2,
            arrowwidth=2
        )
    
    fig.update_layout(
        title='Reaction Network',
        showlegend=True,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
    )
    return fig

def create_kinetics_visualization(time_data, concentration_data, temperature):
    """Create interactive kinetics visualization"""
    # Calculate rate constants using Arrhenius equation
    k = 0.001 * np.exp(-50000/(8.314 * (temperature + 273.15)))
    theoretical = 100 * (1 - np.exp(-k * time_data))
    
    fig = go.Figure()
    
    # Add experimental data
    fig.add_trace(go.Scatter(
        x=time_data,
        y=concentration_data,
        mode='markers',
        name='Experimental'
    ))
    
    # Add theoretical curve
    fig.add_trace(go.Scatter(
        x=time_data,
        y=theoretical,
        mode='lines',
        name='Theoretical'
    ))
    
    # Add rate law equation
    fig.add_annotation(
        x=0.5,
        y=0.95,
        xref='paper',
        yref='paper',
        text=f'Rate Constant (k) = {k:.2e}',
        showarrow=False
    )
    
    fig.update_layout(
        title='Reaction Kinetics',
        xaxis_title='Time (s)',
        yaxis_title='Concentration (%)',
        hovermode='x unified'
    )
    return fig

def create_catalyst_performance_chart(catalyst_data):
    """Create catalyst performance comparison"""
    fig = px.bar(
        x=list(catalyst_data.keys()),
        y=[data['activity'] for data in catalyst_data.values()],
        color=[data['selectivity'] for data in catalyst_data.values()],
        labels={
            'x': 'Catalyst',
            'y': 'Activity (%)',
            'color': 'Selectivity (%)'
        },
        title='Catalyst Performance Comparison'
    )
    
    # Add turnover frequency annotations
    for i, cat in enumerate(catalyst_data.keys()):
        fig.add_annotation(
            x=i,
            y=catalyst_data[cat]['activity'],
            text=f"TOF: {catalyst_data[cat]['tof']:.1f} s⁻¹",
            showarrow=False,
            yshift=10
        )
    
    return fig

def create_process_flow_diagram(unit_operations):
    """Create interactive process flow diagram"""
    fig = go.Figure()
    
    # Position nodes in a flow diagram layout
    positions = {
        'Feed': (0, 0),
        'Reactor': (1, 0),
        'Separator': (2, 0),
        'Recycle': (1.5, 1),
        'Product': (3, 0)
    }
    
    # Add nodes
    for unit, pos in positions.items():
        fig.add_trace(go.Scatter(
            x=[pos[0]], y=[pos[1]],
            mode='markers+text',
            name=unit,
            text=unit,
            marker=dict(size=30, symbol='square'),
            hoverinfo='text',
            hovertext=f"{unit}\nTemp: {unit_operations[unit]['temperature']}°C\nPress: {unit_operations[unit]['pressure']} atm"
        ))
    
    # Add connecting lines
    connections = [
        ('Feed', 'Reactor'),
        ('Reactor', 'Separator'),
        ('Separator', 'Product'),
        ('Separator', 'Recycle'),
        ('Recycle', 'Reactor')
    ]
    
    for start, end in connections:
        fig.add_trace(go.Scatter(
            x=[positions[start][0], positions[end][0]],
            y=[positions[start][1], positions[end][1]],
            mode='lines+text',
            line=dict(width=2, color='black'),
            showlegend=False
        ))
    
    fig.update_layout(
        title='Process Flow Diagram',
        showlegend=False,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
    )
    return fig
