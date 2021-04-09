from dask import distributed, Client
from pysmFISH import selected_logger

logger = selected_logger()

client= Client()
logger.debug(f'dask client: {client}')

