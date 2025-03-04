import yaml
from pathlib import Path

def load_config():
    config_path = Path(__file__).parent / "config.yaml"
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

class Settings:
    def __init__(self):
        config = load_config()
        self.debug_mode = config.get('debug', False)
        self.api_key = config.get('api_key', '')
        self.max_reactions = config.get('max_reactions', 1000)
        self.model_path = config.get('model_path', 'models/')
        self.database_url = config.get('database_url', 'sqlite:///./chemistry.db')
