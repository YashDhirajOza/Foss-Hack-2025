import time
import numpy as np
from datetime import datetime
import plotly.graph_objects as go
from collections import deque

class ReactionMonitor:
    def __init__(self, buffer_size=100):
        self.temperature_buffer = deque(maxlen=buffer_size)
        self.pressure_buffer = deque(maxlen=buffer_size)
        self.conversion_buffer = deque(maxlen=buffer_size)
        self.timestamps = deque(maxlen=buffer_size)
        
    def add_datapoint(self, temp, pressure, conversion):
        """Add new monitoring datapoint"""
        self.temperature_buffer.append(temp)
        self.pressure_buffer.append(pressure)
        self.conversion_buffer.append(conversion)
        self.timestamps.append(datetime.now())
        
    def get_real_time_plot(self):
        """Generate real-time monitoring plot"""
        fig = go.Figure()
        
        # Temperature trace
        fig.add_trace(go.Scatter(
            x=list(self.timestamps), 
            y=list(self.temperature_buffer),
            name='Temperature (°C)',
            yaxis='y1'
        ))
        
        # Pressure trace
        fig.add_trace(go.Scatter(
            x=list(self.timestamps), 
            y=list(self.pressure_buffer),
            name='Pressure (atm)',
            yaxis='y2'
        ))
        
        # Conversion trace
        fig.add_trace(go.Scatter(
            x=list(self.timestamps), 
            y=list(self.conversion_buffer),
            name='Conversion (%)',
            yaxis='y3'
        ))
        
        fig.update_layout(
            title='Real-time Reaction Monitoring',
            height=600,
            yaxis=dict(title='Temperature (°C)', side='left'),
            yaxis2=dict(title='Pressure (atm)', overlaying='y', side='right'),
            yaxis3=dict(title='Conversion (%)', overlaying='y', side='right'),
            showlegend=True
        )
        
        return fig
    
    def get_safety_alerts(self):
        """Check for safety conditions"""
        alerts = []
        
        # Temperature alerts
        if len(self.temperature_buffer) > 0:
            current_temp = self.temperature_buffer[-1]
            if current_temp > 500:
                alerts.append({"level": "Critical", "message": "Temperature exceeds safety limit!"})
            elif current_temp > 400:
                alerts.append({"level": "Warning", "message": "Temperature approaching critical level"})
        
        # Pressure alerts
        if len(self.pressure_buffer) > 0:
            current_pressure = self.pressure_buffer[-1]
            if current_pressure > 300:
                alerts.append({"level": "Critical", "message": "Pressure exceeds safety limit!"})
            elif current_pressure > 250:
                alerts.append({"level": "Warning", "message": "Pressure approaching critical level"})
        
        return alerts
