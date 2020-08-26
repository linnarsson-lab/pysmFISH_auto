from typing import *
import os
import shutil
from prefect import task
from prefect.engine import signals
from pathlib import Path

from pysmFISH.logger_utils import prefect_logging_setup


@task(name='move_data')
def move_data(path_source_location:str,path_destination:str):
    """
    Function used to transfer the files to another location

    Args:
        path_source_location: str
            path to the data to be moved
        path_destination: str
            path to the destination of the transfer

    Returns:
        experiment_path: Posix
            path to the data after transferring
    """

    logger = prefect_logging_setup('transfer_data')
    path_source_location = Path(path_source_location)
    path_destination = Path(path_destination)

    try:
        os.stat(path_source_location)
    except:
        logger.error(f' The {path_source_location} directory is missing')
        fail_signal = signals.FAIL('The source directory is missing')
        fail_signal.flag = True
        fail_signal.value = None
        raise fail_signal
    
    try:
        os.stat(path_destination)
        shutil.move(path_source_location.as_posix(),path_destination.as_posix())
        logger.info(f'data moved from {path_source_location} to {path_destination}')
        return path_destination / path_source_location.stem
    except:
        logger.info(f' The {path_destination} directory is missing')
        fail_signal = signals.FAIL('The destination directory is missing')
        fail_signal.flag = True
        fail_signal.value = None
        raise fail_signal
