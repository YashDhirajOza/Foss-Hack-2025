import pandas as pd
from pathlib import Path
import json

class ExportManager:
    def __init__(self, output_dir="exports"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
    
    def export_reaction_results(self, reaction_data, format="excel"):
        if format == "excel":
            df = pd.DataFrame(reaction_data)
            df.to_excel(self.output_dir / "reaction_results.xlsx")
        elif format == "csv":
            df = pd.DataFrame(reaction_data)
            df.to_csv(self.output_dir / "reaction_results.csv")
        elif format == "json":
            with open(self.output_dir / "reaction_results.json", 'w') as f:
                json.dump(reaction_data, f, indent=2)
