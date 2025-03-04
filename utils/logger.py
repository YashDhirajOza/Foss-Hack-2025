import logging
from datetime import datetime

def setup_logger():
    logger = logging.getLogger('ChemAI')
    logger.setLevel(logging.INFO)
    
    # File handler
    fh = logging.FileHandler(f'logs/chem_ai_{datetime.now().strftime("%Y%m%d")}.log')
    fh.setLevel(logging.INFO)
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)
    
    # Formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger
