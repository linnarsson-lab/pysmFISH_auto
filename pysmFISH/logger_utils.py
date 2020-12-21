"""
Set utils to control logging 
"""
from typing import *
import logging
import time
from pathlib import Path
from pythonjsonlogger import jsonlogger

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
