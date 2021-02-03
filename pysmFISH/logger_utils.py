"""
Set utils to control logging 
"""
from typing import *
import logging
import time
from pathlib import Path
from pythonjsonlogger import jsonlogger

from logging.handlers import RotatingFileHandler



import prefect
from prefect import task
from prefect.utilities.logging import get_logger


def selected_logger():
    """
    Logger function. If you want to change the logging
    procedure replace the content of this function
    """
    # Redirect warnings to logger
    logging.captureWarnings(True)
    # Create logger
    logger = logging.getLogger("distributed.worker")
    logger.setLevel(logging.INFO)
    return logger

def json_logger(log_path,name):
    log_path = Path(log_path)
    date_tag = time.strftime("%y%m%d_%H_%M_%S")
    fname = log_path / (name + '_' + date_tag + '_pipeline_run.log')
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