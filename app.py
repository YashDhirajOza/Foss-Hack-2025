import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from models.chemical_models import ChemicalPredictor
from utils.molecule_viewer import generate_3d_structure, show_molecule_3d
from utils.advanced_visualizations import *
from models.pretrained_models import ChemicalPredictor
from data.reaction_datasets import OpenReactionDataset, ReactionDataModule

# Must be the first Streamlit command
st.set_page_config(page_title="ChemAI FOSS", layout="wide")

# Mock data for demo
def mock_reaction_prediction(reactants, temperature):
    return {
        "products": ["CCO", "H2O"],
        "yield": min(95, 30 + temperature//2),
        "energy_efficiency": "High" if temperature < 150 else "Medium"
    }

def mock_safety_query(chemical):
    return CHEMICALS_DATABASE.get(chemical, {
        "name": "Unknown Chemical",
        "hazards": "No data available",
        "precautions": "Exercise caution",
        "first_aid": "Seek professional advice",
        "storage": "Follow general chemical storage guidelines",
        "incompatible": "Unknown"
    })

def mock_process_optimization(reactants, parameters):
    # Simulate process optimization
    temps = np.linspace(parameters['temp_min'], parameters['temp_max'], 20)
    yields = [30 + t/2 + np.random.normal(0, 5) for t in temps]
    optimal_temp = temps[np.argmax(yields)]
    return {
        "optimal_temp": optimal_temp,
        "max_yield": max(yields),
        "temp_curve": {"temps": temps, "yields": yields}
    }

# Initialize predictors and resources
@st.cache_resource
def load_resources():
    return {
        'chemical': ChemicalPredictor(),
        'reaction': ReactionPredictor(),
        'datasets': OpenReactionDataset(),
        'predictor': ChemicalPredictor(),
        'uspto_data': OpenReactionDataset().load_dataset('uspto')
    }

resources = load_resources()

# Add new constant at the top of file
REACTION_TYPES = {
    "Acid-Base": ["Neutralization", "Salt Formation"],
    "Redox": ["Oxidation", "Reduction", "Single Displacement"],
    "Organic": ["Esterification", "Addition", "Substitution"],
    "Precipitation": ["Double Displacement", "Complex Formation"],
    "Gas Evolution": ["Decomposition", "Metal-Acid"]
}

# Streamlit UI
st.title("🧪 Open-Source Chemical Assistant")

# Sidebar
with st.sidebar:
    st.header("Settings")
    st.selectbox("Theme", ["Light", "Dark"])
    st.slider("Font Size", 12, 24, 16)

# Main Tabs
tab1, tab2, tab3 = st.tabs(["Reaction Lab", "Process Optimizer", "Safety Guide"])

# Enhanced Reaction Lab Tab
with tab1:
    st.header("Chemical Reaction Simulator")
    
    # Add reaction type selector
    reaction_mode = st.radio("Simulation Mode", ["Simple Reaction", "Industrial Process"], horizontal=True)
    
    if reaction_mode == "Simple Reaction":
        # Add reaction category selector
        col1, col2 = st.columns([1, 2])
        with col1:
            reaction_category = st.selectbox("Reaction Category", list(REACTION_TYPES.keys()))
        with col2:
            reaction_subtype = st.selectbox("Reaction Type", REACTION_TYPES[reaction_category])
        
        st.markdown("### Select Reactants")
        col1, col2, col3 = st.columns([2, 1, 2])
        
        with col1:
            chemical_1 = st.selectbox("Chemical 1", list(CHEMICALS_DATABASE.keys()))
            amount_1 = st.number_input("Amount (g)", 1.0, 1000.0, 100.0, key="amount1")
            
            # Add concentration for solutions
            concentration_1 = st.number_input("Concentration (M)", 0.1, 18.0, 1.0, key="conc1")
            
        with col2:
            st.markdown("<div style='text-align: center; font-size: 40px; padding: 30px'>+</div>", unsafe_allow_html=True)
            st.markdown("<div style='text-align: center;'>Reaction Type:<br/>" + reaction_subtype + "</div>", unsafe_allow_html=True)
        
        with col3:
            chemical_2 = st.selectbox("Chemical 2", list(CHEMICALS_DATABASE.keys()))
            amount_2 = st.number_input("Amount (g)", 1.0, 1000.0, 100.0, key="amount2")
            concentration_2 = st.number_input("Concentration (M)", 0.1, 18.0, 1.0, key="conc2")
        

        # In the 3D structure section
        if chemical_1 and chemical_2:
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("### 3D Structure - Reactant 1")
                try:
                    # Try to get structure from formula if SMILES not available
                    formula_1 = CHEMICALS_DATABASE[chemical_1]['formula']
                    mol_3d_1 = generate_3d_structure(formula_1)
                    if mol_3d_1:
                        show_molecule_3d(mol_3d_1, size=(300, 300))
                    else:
                        st.info(f"Cannot generate 3D structure for {formula_1}")
                except Exception as e:
                    st.info("3D visualization not available for this molecule")
            
            with col2:
                st.markdown("### 3D Structure - Reactant 2")
                try:
                    formula_2 = CHEMICALS_DATABASE[chemical_2]['formula']
                    mol_3d_2 = generate_3d_structure(formula_2)
                    if mol_3d_2:
                        show_molecule_3d(mol_3d_2, size=(300, 300))
                    else:
                        st.info(f"Cannot generate 3D structure for {formula_2}")
                except Exception as e:
                    st.info("3D visualization not available for this molecule")

        
        # Reaction Conditions with more options
        st.markdown("### Reaction Conditions")
        col1, col2, col3 = st.columns(3)
        with col1:
            temperature = st.slider("Temperature (°C)", -50, 500, 25)
            pressure = st.slider("Pressure (atm)", 1, 300, 1)
        with col2:
            catalyst = st.selectbox("Catalyst", ["None", "H+", "OH-", "Pd/C", "Pt", "Fe", "Ni", "V2O5", "MnO2"])
            solvent = st.selectbox("Solvent", ["None", "Water", "Ethanol", "Methanol", "THF", "DMF", "DMSO", "Acetone"])
        with col3:
            st.markdown("### Additional Parameters")
            stirring = st.checkbox("Stirring")
            if stirring:
                st.slider("Stirring Speed (rpm)", 100, 1000, 400)
            heat_rate = st.slider("Heating Rate (°C/min)", 1, 50, 10)

        # Enhanced Simulation Button
        if st.button("🧪 Simulate Reaction", type="primary", use_container_width=True):
            try:
                with st.spinner("Simulating reaction..."):
                    conditions = {
                        'temperature': temperature,
                        'pressure': pressure,
                        'catalyst': catalyst,
                        'solvent': solvent,
                        'stirring': stirring,
                        'heat_rate': heat_rate,
                        'concentration_1': concentration_1,
                        'concentration_2': concentration_2
                    }
                    
                    # Check if chemicals exist in database
                    if chemical_1 not in CHEMICALS_DATABASE or chemical_2 not in CHEMICALS_DATABASE:
                        st.error("One or both chemicals not found in database!")
                    
                    
                    prediction = predictors['reaction'].predict_product(
                        CHEMICALS_DATABASE[chemical_1],
                        CHEMICALS_DATABASE[chemical_2],
                        conditions
                    )
                    
                    if prediction and prediction.get('success', False):
                        # Display reaction equation
                        st.markdown(f"""
                        ## Reaction Result
                        <div style='text-align: center; font-size: 24px; padding: 20px; 
                             background-color: #262730; color: #ffffff; 
                             border-radius: 10px; border: 2px solid #4CAF50; margin: 10px 0;
                             box-shadow: 0 2px 4px rgba(0,0,0,0.1);'>
                            <span style='color: #4CAF50'>{CHEMICALS_DATABASE[chemical_1]['formula']}</span>
                            <span style='color: #ffffff'> + </span>
                            <span style='color: #4CAF50'>{CHEMICALS_DATABASE[chemical_2]['formula']}</span>
                            <span style='color: #ffffff'> → </span>
                            <span style='color: #FFA500'>{' + '.join(prediction['products'])}</span>
                        </div>
                        """, unsafe_allow_html=True)
                        
                        # Update mechanism visualization section
                        st.subheader("Reaction Mechanism")
                        # Get reaction type and predict mechanism
                        current_reaction_type = get_reaction_type(chemical_1, chemical_2)
                        st.info(f"Reaction Type: {current_reaction_type}")
                        
                        mechanism_steps = draw_reaction_mechanism(
                            current_reaction_type,
                            [chemical_1, chemical_2],
                            prediction['products']
                        )
                        
                        if mechanism_steps:
                            st.markdown("""
                            <style>
                                .mechanism-step { 
                                    background-color: #1e1e1e;
                                    padding: 15px;
                                    border-radius: 10px;
                                    margin: 5px 0;
                                }
                                .step-description {
                                    color: #4CAF50;
                                    font-size: 14px;
                                    margin: 5px 0;
                                }
                                .arrow-text {
                                    color: #FFA500;
                                    text-align: center;
                                    font-size: 16px;
                                    margin: 5px 0;
                                }
                                .mechanism-grid {
                                    display: grid;
                                    grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
                                    gap: 10px;
                                }
                            </style>
                            """, unsafe_allow_html=True)
                            
                            st.markdown("<div class='mechanism-grid'>", unsafe_allow_html=True)
                            for i, (step_name, img, description, arrow) in enumerate(mechanism_steps, 1):
                                st.markdown(f"""
                                <div class='mechanism-step'>
                                    <h4>Step {i}: {step_name}</h4>
                                    <div class='step-description'>{description}</div>
                                """, unsafe_allow_html=True)
                                st.image(img, use_column_width=True)
                                if i < len(mechanism_steps):
                                    st.markdown(f"<div class='arrow-text'>{arrow} ↓</div>", 
                                              unsafe_allow_html=True)
                                st.markdown("</div>", unsafe_allow_html=True)
                            st.markdown("</div>", unsafe_allow_html=True)
                        
                        else:
                            st.warning("No detailed mechanism available for this reaction type.")

                
                        col1, col2 = st.columns(2)
                        with col1:
                            st.subheader("Reaction Details")
                            st.metric("Predicted Yield", f"{prediction['yield']:.1f}%")
                            st.metric("Reaction Time", f"{prediction['time']} hours")
                            st.metric("Energy Change", f"{prediction['energy']} kJ/mol")
                            st.markdown(f"**Mechanism**: {prediction['mechanism']}")
                            st.markdown(f"**Type**: {prediction['reaction_type']}")
                        
                        with col2:
                            st.subheader("Safety & Conditions")
                            for hazard in prediction['hazards']:
                                st.warning(f"⚠️ {hazard}")
                            
                            st.info(f"""
                            **Optimal Conditions:**
                            - Temperature: {conditions['temperature']}°C
                            - Pressure: {conditions['pressure']} atm
                            - Catalyst: {conditions['catalyst']}
                            - Solvent: {conditions['solvent']}
                            """)
                        
                        st.subheader("Reaction Progress")
                        progress_fig = go.Figure()
                        times = np.linspace(0, prediction['time']*60, 20)
                        progress = 100 * (1 - np.exp(-times/30))
                        progress_fig.add_trace(go.Scatter(
                            x=times, y=progress,
                            mode='lines+markers',
                            name='Conversion'
                        ))
                        progress_fig.update_layout(
                            title="Reaction Progress",
                            xaxis_title="Time (minutes)",
                            yaxis_title="Conversion (%)"
                        )
                        st.plotly_chart(progress_fig, use_container_width=True)

                        st.subheader("Reaction Visualization")
                        col1, col2 = st.columns(2)
                        with col1:
                            # Energy profile
                            energy_profile = create_energy_profile(prediction, conditions)
                            st.plotly_chart(energy_profile, use_container_width=True)
                        
                        with col2:
                            # Green chemistry metrics
                            metrics = calculate_atom_economy(
                                CHEMICALS_DATABASE[chemical_1],
                                CHEMICALS_DATABASE[chemical_2],
                                [CHEMICALS_DATABASE.get(p, {}) for p in prediction['products']]
                            )
                            if metrics:
                                st.subheader("Green Chemistry Metrics")
                                st.metric("Atom Economy", f"{metrics['atom_economy']}%")
                                st.metric("Waste Factor", f"{metrics['waste_factor']}%")
                                st.info(f"Efficiency Class: {metrics['efficiency_class']}")
                    
                        # After displaying the basic simulation results, add new visualizations
                        st.subheader("Advanced Visualizations")
                        
                        # Create tabs for different visualizations
                        viz_tab1, viz_tab2, viz_tab3 = st.tabs(["Energy Analysis", "Kinetics", "Process"])
                        
                        with viz_tab1:
                            # 3D Energy Surface
                            energy_surface = create_3d_energy_surface(
                                [temperature-100, temperature+100],
                                [pressure-50, pressure+50],
                                prediction
                            )
                            st.plotly_chart(energy_surface, use_container_width=True)
                            
                            # Reaction Network
                            network = create_reaction_network(
                                [CHEMICALS_DATABASE[chemical_1]['formula'], 
                                 CHEMICALS_DATABASE[chemical_2]['formula']],
                                ["Intermediate Complex"],
                                prediction['products']
                            )
                            st.plotly_chart(network, use_container_width=True)
                        
                        with viz_tab2:
                            # Kinetics Visualization
                            times = np.linspace(0, prediction['time']*3600, 50)  # Convert hours to seconds
                            concentrations = 100 * (1 - np.exp(-times/1000))  # Mock concentration data
                            kinetics = create_kinetics_visualization(times, concentrations, temperature)
                            st.plotly_chart(kinetics, use_container_width=True)
                            
                            # Catalyst Performance
                            if catalyst != "None":
                                catalyst_data = {
                                    "None": {"activity": 40, "selectivity": 60, "tof": 0.5},
                                    catalyst: {"activity": 85, "selectivity": 90, "tof": 2.5},
                                }
                                catalyst_chart = create_catalyst_performance_chart(catalyst_data)
                                st.plotly_chart(catalyst_chart, use_container_width=True)
                        
                        with viz_tab3:
                            # Process Flow Diagram
                            unit_operations = {
                                "Feed": {"temperature": 25, "pressure": 1},
                                "Reactor": {"temperature": temperature, "pressure": pressure},
                                "Separator": {"temperature": temperature-20, "pressure": pressure-5},
                                "Recycle": {"temperature": temperature-10, "pressure": pressure},
                                "Product": {"temperature": 30, "pressure": 1}
                            }
                            flow_diagram = create_process_flow_diagram(unit_operations)
                            st.plotly_chart(flow_diagram, use_container_width=True)
                    
                    else:
                        error_msg = prediction.get('error', 'Unknown error occurred') if prediction else 'No prediction returned'
                        st.error(f"Simulation Error: {error_msg}")
                        st.info("Try different chemicals or conditions")
                        
            except Exception as e:
                st.error(f"An error occurred: {str(e)}")
                st.info("Please try again with different parameters")

    else:
        # Industrial Process Interface
        st.markdown("### Industrial Process Interface")
        
        # Process Selection
        process_type = st.selectbox("Select Process", 
                                  ["Haber Process", "Sulfuric Acid Production", 
                                   "Methanol Synthesis", "Ethylene Oxide Production"])
        
        # Parameter Configuration
        st.subheader("Process Parameters")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            temp_setpoint = st.number_input("Temperature Setpoint (°C)", 0, 1000, 450)
            pressure = st.number_input("Pressure (atm)", 1, 500, 200)
            feed_rate = st.number_input("Feed Rate (kg/hr)", 10, 1000, 100)
            
        with col2:
            catalyst_loading = st.number_input("Catalyst Loading (%)", 0.0, 100.0, 5.0)
            energy_cost = st.number_input("Energy Cost ($/kWh)", 0.01, 1.0, 0.12)
            raw_material_cost = st.number_input("Raw Material Cost ($/kg)", 0.1, 100.0, 2.0)
            
        with col3:
            labor_cost = st.number_input("Labor Cost ($/hr)", 10, 200, 50)
            maintenance_factor = st.slider("Maintenance Factor", 0.0, 1.0, 0.1)
            recycle_ratio = st.slider("Recycle Ratio", 0.0, 1.0, 0.8)
        
        # Process Parameters
        params = {
            'temp_setpoint': temp_setpoint,
            'pressure': pressure,
            'feed_rate': feed_rate,
            'catalyst_loading': catalyst_loading,
            'energy_cost': energy_cost,
            'raw_material_cost': raw_material_cost,
            'labor_cost': labor_cost
        }
        
        # Simulation Control
        col1, col2 = st.columns(2)
        with col1:
            if st.button("🏭 Start Simulation", type="primary", use_container_width=True):
                with st.spinner("Running industrial process simulation..."):
                    simulator = IndustrialProcessSimulator(process_type, params)
                    results = simulator.simulate_batch_reactor()
                    
                    # Display Dashboard
                    st.plotly_chart(simulator.create_dashboard(results), use_container_width=True)
                    
                    # Display Key Metrics
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Final Conversion", f"{results['final_yield']:.1f}%")
                        st.metric("Total Product", f"{results['accumulated_product'][-1]:.0f} kg")
                    
                    with col2:
                        st.metric("Energy Usage", f"{np.mean(results['power_consumption']):.1f} kW")
                        st.metric("Operating Cost", f"${results['total_cost']:.2f}/hr")
                    
                    with col3:
                        st.metric("Catalyst Activity", f"{results['catalyst_activity'][-1]*100:.1f}%")
                        roi = ((results['accumulated_product'][-1] * 100) - 
                               results['total_cost']) / results['total_cost'] * 100
                        st.metric("ROI", f"{roi:.1f}%")
                    
                    # Equipment Status
                    st.subheader("Equipment Status")
                    status = results['equipment_status']
                    for equipment, details in status.items():
                        with st.expander(f"{equipment.title()} Status"):
                            for key, value in details.items():
                                st.text(f"{key.title()}: {value}")
        
        with col2:
            if st.button("🎯 Optimize Process", type="primary", use_container_width=True):
                with st.spinner("Optimizing process conditions..."):
                    simulator = IndustrialProcessSimulator(process_type, params)
                    optimal = simulator.optimize_conditions(objective='profit')
                    
                    st.success("Optimization Complete!")
                    st.json({
                        "Optimal Temperature": f"{optimal['temperature']:.1f}°C",
                        "Optimal Pressure": f"{optimal['pressure']:.1f} atm",
                        "Expected Improvement": f"{optimal['metric']:.1f}%"
                    })

# Process Optimization Tab
with tab2:
    st.header("Process Optimization")
    
    col1, col2 = st.columns(2)
    with col1:
        reactants = st.text_input("Reactants (SMILES)", "C=O, [H]", key="opt_reactants")
        params = {
            'temp_min': st.number_input("Min Temperature (°C)", 0, 200, 25),
            'temp_max': st.number_input("Max Temperature (°C)", 0, 300, 150)
        }
    
    if st.button("Optimize Process", type="primary"):
        opt_result = mock_process_optimization(reactants.split(","), params)
        
        with col2:
            st.subheader("Optimization Results")
            st.metric("Optimal Temperature", f"{opt_result['optimal_temp']:.1f}°C")
            st.metric("Maximum Yield", f"{opt_result['max_yield']:.1f}%")
            
            # Temperature vs Yield curve
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=opt_result['temp_curve']['temps'],
                y=opt_result['temp_curve']['yields'],
                mode='lines+markers',
                name='Yield Curve'
            ))
            fig.update_layout(
                title="Temperature vs. Yield Optimization",
                xaxis_title="Temperature (°C)",
                yaxis_title="Yield (%)"
            )
            st.plotly_chart(fig, use_container_width=True)

# Safety Assistant Tab
with tab3:
    st.header("Chemical Safety Advisor")
    chemical = st.selectbox("Select Chemical", list(CHEMICALS_DATABASE.keys()))
    
    if chemical:
        info = CHEMICALS_DATABASE[chemical]
        
        st.subheader(f"{info['name']} ({info['formula']})")
        
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("### ⚛️ Properties")
            st.info(f"""
            - Molar Mass: {info['molar_mass']}
            - Density: {info['density']}
            - Melting Point: {info['melting_point']}
            - Boiling Point: {info['boiling_point']}
            """)
            
            st.markdown("### 🔬 Uses")
            st.write(", ".join(info['uses']))
        
        with col2:
            st.markdown("### ⚠️ Safety")
            st.warning(info['hazards'])
            st.error(f"First Aid: {info['first_aid']}")
            
        st.markdown("### 📝 Description")
        st.write(info['description'])

def calculate_reaction_costs(chem1, amt1, chem2, amt2, conditions):
    """Calculate reaction costs including materials, energy, and operation"""
    # Implementation of cost calculation
    return {
        'materials': 100,
        'energy': 50,
        'operation': 30,
        'total': 180
    }

