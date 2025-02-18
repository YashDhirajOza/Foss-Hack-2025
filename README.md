🧪 AI-Powered Chemical Reaction Optimizer 
📢 Hackathon Project for FOSS-Hack-2025

The AI-Powered Chemical Reaction Optimizer is an open-source platform designed to revolutionize chemical R&D, industrial process optimization, and laboratory safety using machine learning and AI-driven automation. It combines chemical reaction prediction, real-time process optimization, and an AI-powered safety assistant, making it a powerful tool for chemists, researchers, and industrial engineers. This solution helps pharmaceutical companies, material scientists, and chemical industries optimize reaction conditions, reduce waste, improve safety compliance, and accelerate innovation.

🔗 Objective: Build an open-source AI-powered platform to:
✅ Predict chemical reactions using ML models.
✅ Optimize industrial chemical processes using simulation & data-driven tuning.
✅ Ensure chemical safety with an NLP-powered chatbot.
✅ Visualize key insights in an interactive Streamlit dashboard.

🛠️ Technology Stack (100% Open Source)
Component	Open-Source Technology (License)
|ML Model	      |PyTorch (MIT), RDKit (BSD), Scikit-Learn
------------------------------------------------------
|Chatbot NLP	  |Hugging Face Transformers (Apache 2.0), SpaCy (MIT)
------------------------------------------------------
|Backend	      |Django (BSD), FastAPI (MIT)
------------------------------------------------------
|Frontend	      |Streamlit (Apache 2.0)
------------------------------------------------------
|Database	      |PostgreSQL (PostgreSQL License)
------------------------------------------------------
|Visualization  |Plotly (MIT), Matplotlib (PSF License)
------------------------------------------------------
|Simulation	    |ChemPy (BSD), SciPy (BSD)
------------------------------------------------------
|Deployment	    |Docker (Apache 2.0), NGINX (BSD)

🚀 Key Features
🔬 1. AI-Powered Reaction Prediction
✅ Data Source: Uses PubChem Reaction Dataset (open data).
✅ ML Model: Trained on USPTO public reaction datasets with RDKit-based featurization.
✅ Output: Predicts reaction products and yield.
✅ Visualization: Displays reaction pathways using Plotly.

🏭 2. Process Optimization & Simulation
✅ Simulation: Uses ChemPy for modeling reaction kinetics.
✅ Optimization: Utilizes SciPy to optimize temperature, pressure, catalysts for efficiency.

🛑 3. AI-Powered Safety Assistant
✅ NLP Model: Fine-tuned DistilBERT model on open MSDS (Material Safety Data Sheets) databases.
✅ Data Source: PubChem, OSHA MSDS sheets (scraped with licensing compliance).
✅ Functionality: Provides safety protocols, handling guidelines, emergency actions.

📊 4. Real-Time Dashboard
✅ Streamlit UI: Displays ML predictions, safety alerts, and reaction optimization results.
✅ Graphical Insights: Cost/waste reduction metrics using Plotly.
✅ License Alerts: Detects chemicals with restrictive patents via PubChem data.

📂 Project Structure

Edit
├── LICENSE               # MIT or Apache 2.0 License
├── README.md             # Setup guide with FOSS ethos
├── requirements.txt      # Dependencies
├── app
│   ├── backend           # FastAPI routes for ML & NLP
│   ├── frontend          # Streamlit UI
│   ├── ml                # PyTorch/RDKit models
│   ├── safety_bot        # Hugging Face NLP pipeline
├── data                  # Sample PubChem datasets
└── tests                 # pytest scripts
