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



def setup_logger(log_fpath:str):
    """
    Logger setup for writing logs output to files on monod
    This logs works when activated under @task decorator.
    In order to get the output from all mapped task make sure
    to set 'a' in the file handler.

    Args:
    -----
        log_fpath:str
            fpath to the log file
    
    Returns:
    --------
        logger: prefect logger
            logger to use in the specific task set

    """
    logger = prefect.context.get("logger")
    file_handler = logging.FileHandler(str(log_fpath),'a')
    file_handler.setLevel(logging.DEBUG)
    format_str = '%(message)%(levelname)%(name)%(asctime)'
    formatter = jsonlogger.JsonFormatter(format_str)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return logger







@task(name='setup_extra_loggers')
def setup_extra_loggers():
    """
    Function to format and add extra logs for  monitoring
    the different steps of the pipeline
    """
    extra_logs_list = ['parsing','preprocessing','test']
    prefect.config.logging.extra_loggers = extra_logs_list



def prefect_logging_setup(extra_log_name:str):
    """
    Function to select a preset extra log to collect the
    output on prefect local or cloud

    Args:
        extra_log_name: str
            logger to activate, the name must be included in 
            prefect.config.logging.extra_loggers
    """

    # general_logger = prefect.context.get("logger")
    # general_logger.debug(f'{prefect.config.logging.extra_loggers}')
    # #assert extra_log_name in prefect.config.logging.extra_loggers, general_logger.error('logger name not in prefect.config.logging.extra_loggers')
    try:
        logger = get_logger(extra_log_name)
        logger.setLevel(logging.DEBUG)
        return logger
    except:
        raise ValueError('logger name not in prefect.config.logging.extra_loggers')



@task()
def test_logger(log_stdout=True):
    logger = prefect.utilities.logging.get_logger("preprocessing")
    logger.setLevel(logging.DEBUG)
    logger.debug("An info message.")
    logger.warning("A warning message.")


# Make written output as extra util function
# @task()

def test_write_logs_to_file(logs_folder_path:str, log_name:str, log_stdout=True):
    logs_folder_path = Path(logs_folder_path)
    logger = prefect.utilities.logging.get_logger(log_name)
    date_tag = time.strftime("%y%m%d_%H_%M_%S")
    # fname = logs_folder_path / (log_name + '_' + date_tag + '_pysmFISH_analysis.log')
    fname = logs_folder_path / (log_name + '.log')
    # logger = logging.getLogger(log_name)
    handler = logging.FileHandler(str(fname))
    format_str = '%(message)%(levelname)%(name)%(asctime)'
    formatter = jsonlogger.JsonFormatter(format_str)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.info("logger started.")
    
    return logger


def simple_writing_logger():
    logger = prefect.context.get("logger")
    file_handler = logging.FileHandler('/wsfish/smfish_ssd/LBEXP20201014_EEL_Mouse_2420um_auto/prefect_logs/my.log')
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)
    return logger
