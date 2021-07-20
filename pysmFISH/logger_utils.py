"""
Set utils to control logging in htcondor and dask
"""
from typing import *
import logging
import time
from pathlib import Path
from pythonjsonlogger import jsonlogger
from logging import Logger


def selected_logger()->Logger:
    """Logger function used inside all the modules, classes and
    functions. If you want to change the logging procedure in the code 
    just replace the content of this function

    Returns:
        logger (logger): selected type of logger
    """
        
    # Redirect warnings to logger
    logging.captureWarnings(True)
    # Create logger
    logger = logging.getLogger("distributed.worker")
    logger.setLevel(logging.INFO)
    return logger

def json_logger(log_path:str,name:str)->logging.logger:
    """Save the logs in a JSON file that can be easily
    parsed.

    Args:
        log_path (str): path where to save the log file
        name (str): name of the logging file

    Returns:
        logger (logger): json logger
    """

    log_path = Path(log_path)
    date_tag = time.strftime("%y%m%d_%H_%M_%S")
    fname = log_path / (name + '_' + date_tag + '.log')
    handler = logging.FileHandler(str(fname))  # Or FileHandler or anything else
    # Configure the fields to include in the JSON output. message is the main log string itself
    format_str = '%(message)%(levelname)%(name)%(asctime)'
    formatter = jsonlogger.JsonFormatter(format_str)
    handler.setFormatter(formatter)
    logger = logging.getLogger()
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    # Normally we would attach the handler to the root logger, and this would be unnecessary
    logger.propagate = True
    logger.info('Created logger')
    return logger