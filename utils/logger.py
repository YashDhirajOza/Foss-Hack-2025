"""Module for handling application logging."""

import logging
from datetime import datetime
from pathlib import Path

def setup_logger(name: str, log_dir: str = "logs") -> logging.Logger:
    """Configure and return a logger instance.
    
    Args:
        name: Logger name
        log_dir: Directory for log files
        
    Returns:
        Configured logger instance
    """
    log_path = Path(log_dir)
    log_path.mkdir(exist_ok=True)
    
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    file_handler = logging.FileHandler(
        log_path / f"{name}_{datetime.now():%Y%m%d}.log",
        encoding='utf-8'
    )
    file_handler.setFormatter(formatter)
    
    logger.addHandler(file_handler)
    return logger
