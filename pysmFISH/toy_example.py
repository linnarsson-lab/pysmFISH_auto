from dask import distributed, Client
from pysmFISH.logger_utils import selected_logger
import time

logger = selected_logger()

client= Client()
logger.debug(f'dask client: {client}')
time.sleep(500)