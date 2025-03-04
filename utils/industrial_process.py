import plotly.graph_objects as go
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional

@dataclass
class ProcessParameters:
    temp_setpoint: float
    pressure: float
    flow_rate: float
    catalyst_loading: float
    residence_time: float
    recycle_ratio: float

class IndustrialProcessSimulator:
    def __init__(self, process_name, parameters):
        self.process = process_name
        self.params = parameters
        self.time_steps = np.linspace(0, 24, 100)
        self.equipment_states = {
            'reactor': 'running',
            'heat_exchanger': 'normal',
            'pump': 'running',
            'separator': 'normal'
        }
    
    def simulate_batch_reactor(self):
        """Simulate batch reactor process with detailed parameters"""
        # Temperature profile with disturbances
        random_disturbance = np.random.normal(0, 1, len(self.time_steps))
        temp_profile = (self.params['temp_setpoint'] + 
                       np.sin(self.time_steps * np.pi/12) * 2 + 
                       random_disturbance)
        
        # Pressure profile
        pressure_profile = self.params['pressure'] * (1 + 0.05 * np.sin(self.time_steps * np.pi/6))
        
        # Conversion profile with catalyst deactivation
        catalyst_activity = np.exp(-0.05 * self.time_steps)
        conversion = 100 * (1 - np.exp(-self.time_steps/5)) * catalyst_activity
        
        # Mass balance
        feed_rate = self.params.get('feed_rate', 100)  # kg/hr
        product_rate = feed_rate * (conversion/100)
        accumulated_product = np.cumsum(product_rate) * (self.time_steps[1] - self.time_steps[0])
        
        # Energy consumption
        heating_power = 50 * temp_profile/100  # kW
        pump_power = 20 * np.ones_like(self.time_steps)  # kW
        total_power = heating_power + pump_power
        
        # Cost calculations
        raw_material_cost = feed_rate * self.params['raw_material_cost']
        energy_cost = total_power * self.params['energy_cost']
        labor_cost = self.params['labor_cost'] * np.ones_like(self.time_steps)
        maintenance_cost = 0.1 * raw_material_cost  # 10% of raw material cost
        
        total_cost = (raw_material_cost + energy_cost + labor_cost + maintenance_cost)
        
        return {
            'temperature': temp_profile,
            'pressure': pressure_profile,
            'conversion': conversion,
            'catalyst_activity': catalyst_activity,
            'product_rate': product_rate,
            'accumulated_product': accumulated_product,
            'power_consumption': total_power,
            'costs': {
                'raw_material': raw_material_cost,
                'energy': energy_cost,
                'labor': labor_cost,
                'maintenance': maintenance_cost,
                'total': total_cost
            },
            'total_cost': np.sum(total_cost),
            'final_yield': conversion[-1],
            'equipment_status': self.monitor_equipment()
        }
    
    def monitor_equipment(self):
        """Monitor equipment status"""
        return {
            'reactor': {
                'status': self.equipment_states['reactor'],
                'temperature': self.params['temp_setpoint'],
                'pressure': self.params['pressure'],
                'agitation': 'normal'
            },
            'heat_exchanger': {
                'status': self.equipment_states['heat_exchanger'],
                'efficiency': 85,
                'fouling_factor': 0.02
            },
            'pump': {
                'status': self.equipment_states['pump'],
                'flow_rate': self.params.get('feed_rate', 100),
                'head': 'normal'
            },
            'separator': {
                'status': self.equipment_states['separator'],
                'efficiency': 90
            }
        }
    
    def create_dashboard(self, simulation_results):
        """Create comprehensive process dashboard"""
        fig = go.Figure()
        
        # Process variables
        fig.add_trace(go.Scatter(
            x=self.time_steps,
            y=simulation_results['temperature'],
            name='Temperature (°C)',
            yaxis='y1'
        ))
        
        fig.add_trace(go.Scatter(
            x=self.time_steps,
            y=simulation_results['pressure'],
            name='Pressure (atm)',
            yaxis='y2'
        ))
        
        fig.add_trace(go.Scatter(
            x=self.time_steps,
            y=simulation_results['conversion'],
            name='Conversion (%)',
            yaxis='y3'
        ))
        
        # Production rate
        fig.add_trace(go.Scatter(
            x=self.time_steps,
            y=simulation_results['product_rate'],
            name='Production Rate (kg/hr)',
            yaxis='y4'
        ))
        
        fig.update_layout(
            title='Industrial Process Dashboard',
            height=800,
            grid={'rows': 2, 'columns': 2, 'pattern': 'independent'},
            yaxis1=dict(title='Temperature (°C)', domain=[0.5, 1]),
            yaxis2=dict(title='Pressure (atm)', domain=[0.5, 1]),
            yaxis3=dict(title='Conversion (%)', domain=[0, 0.4]),
            yaxis4=dict(title='Production (kg/hr)', domain=[0, 0.4])
        )
        
        return fig
    
    def optimize_conditions(self, objective='profit'):
        """Optimize process conditions"""
        results = []
        
        # Grid search over temperature and pressure
        temperatures = np.linspace(self.params['temp_setpoint']*0.8, 
                                 self.params['temp_setpoint']*1.2, 5)
        pressures = np.linspace(self.params['pressure']*0.8,
                              self.params['pressure']*1.2, 5)
        
        for temp in temperatures:
            for press in pressures:
                test_params = self.params.copy()
                test_params['temp_setpoint'] = temp
                test_params['pressure'] = press
                
                sim = IndustrialProcessSimulator(self.process, test_params)
                result = sim.simulate_batch_reactor()
                
                if objective == 'profit':
                    metric = result['accumulated_product'][-1] - result['total_cost']
                elif objective == 'conversion':
                    metric = result['conversion'][-1]
                
                results.append({
                    'temperature': temp,
                    'pressure': press,
                    'metric': metric
                })
        
        best_result = max(results, key=lambda x: x['metric'])
        return best_result
