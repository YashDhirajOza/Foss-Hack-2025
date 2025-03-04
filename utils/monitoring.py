import numpy as np
import plotly.graph_objects as go

class ReactionMonitor:
    def __init__(self):
        self.data = {
            'temperature': [],
            'pressure': [],
            'conversion': [],
            'time': []
        }
    
    def update(self, time, temp, pressure, conversion):
        """Update reaction parameters"""
        self.data['time'].append(time)
        self.data['temperature'].append(temp)
        self.data['pressure'].append(pressure)
        self.data['conversion'].append(conversion)
    
    def create_dashboard(self):
        """Create real-time monitoring dashboard"""
        fig = go.Figure()
        
        # Temperature trace
        fig.add_trace(go.Scatter(
            x=self.data['time'],
            y=self.data['temperature'],
            name='Temperature (°C)',
            yaxis='y1'
        ))
        
        # Conversion trace
        fig.add_trace(go.Scatter(
            x=self.data['time'],
            y=self.data['conversion'],
            name='Conversion (%)',
            yaxis='y2'
        ))
        
        fig.update_layout(
            title='Real-time Reaction Monitoring',
            yaxis1=dict(title='Temperature (°C)'),
            yaxis2=dict(title='Conversion (%)', overlaying='y', side='right')
        )
        
        return fig
