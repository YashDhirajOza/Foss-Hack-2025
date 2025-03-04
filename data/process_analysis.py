import plotly.graph_objects as go
import numpy as np

def create_process_visualization(reaction_data):
    """Create industrial process visualizations"""
    figs = []
    
    # Temperature vs Conversion
    temp_data = reaction_data['optimization_metrics']['conversion_vs_temperature']
    fig1 = go.Figure()
    fig1.add_trace(go.Scatter(
        x=temp_data['temps'],
        y=temp_data['conversion'],
        mode='lines+markers',
        name='Conversion'
    ))
    fig1.update_layout(
        title="Temperature vs Conversion",
        xaxis_title="Temperature (°C)",
        yaxis_title="Conversion (%)"
    )
    figs.append(fig1)
    
    # Pressure vs Yield
    pressure_data = reaction_data['optimization_metrics']['pressure_effect']
    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(
        x=pressure_data['pressures'],
        y=pressure_data['yield'],
        mode='lines+markers',
        name='Yield'
    ))
    fig2.update_layout(
        title="Pressure vs Yield",
        xaxis_title="Pressure (atm)",
        yaxis_title="Yield (%)"
    )
    figs.append(fig2)
    
    # Economics Radar Chart
    fig3 = go.Figure()
    fig3.add_trace(go.Scatterpolar(
        r=[80, 60, 90, 70, 85],
        theta=['ROI', 'Efficiency', 'Safety', 'Sustainability', 'Quality'],
        fill='toself'
    ))
    fig3.update_layout(
        title="Process Performance Metrics",
        polar=dict(radialaxis=dict(visible=True, range=[0, 100]))
    )
    figs.append(fig3)
    
    return figs
