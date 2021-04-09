#!/home/simone/mini/envs/test_d/bin/python
from dask.distributed import Client
from pysmFISH.logger_utils import selected_logger
import time

logger = selected_logger()

client= Client()
logger.debug(f'dask client: {client}')
time.sleep(500)