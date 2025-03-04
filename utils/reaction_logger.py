import logging
import json
from datetime import datetime
from pathlib import Path

class ReactionLogger:
    def __init__(self, log_dir="logs"):
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(exist_ok=True)
        self._setup_logger()
    
    def _setup_logger(self):
        """Setup logging configuration"""
        self.logger = logging.getLogger('ReactionLogger')
        self.logger.setLevel(logging.INFO)
        
        # File handler for detailed logs
        fh = logging.FileHandler(
            self.log_dir / f'reactions_{datetime.now().strftime("%Y%m%d")}.log'
        )
        fh.setFormatter(logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s'
        ))
        self.logger.addHandler(fh)
    
    def log_reaction(self, reaction_data):
        """Log reaction details"""
        # Log to file
        self.logger.info(json.dumps(reaction_data))
        
        # Save detailed data
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        with open(self.log_dir / f'reaction_{timestamp}.json', 'w') as f:
            json.dump(reaction_data, f, indent=2)
    
    def log_error(self, error_data):
        """Log error details"""
        self.logger.error(json.dumps(error_data))
    
    def get_reaction_history(self, days=7):
        """Retrieve reaction history"""
        history = []
        for log_file in self.log_dir.glob('reaction_*.json'):
            if (datetime.now() - datetime.fromtimestamp(log_file.stat().st_mtime)).days <= days:
                with open(log_file) as f:
                    history.append(json.load(f))
        return history
