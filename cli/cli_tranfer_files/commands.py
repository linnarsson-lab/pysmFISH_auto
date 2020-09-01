import click
import time
from pathlib import Path

from pysmFISH.flows.transfer_data_flow import transfer_data_flow

@click.group('data-transfer')
def data_transfer():
    """
    group of commands used to transfer data between HD
    """
    pass

@click.command('transfer-data-flow')
@click.option('--config_db_path', type=str, help='Path to the directory with the config files with the \
                 parameters for data transfering')
def transfer_data_to_processing_hd(config_db_path:str):
    """
    This prefect flow is used to transfer data from a temporary storage to the 
    processing HD.
    The flow scan for new data once a day

    """
    click.echo('Start looking for data completely tranferred from the machines')
    transfer_data_flow(config_db_path)
    click.echo('--------------------------------------------------------------')
    click.echo('scanning of the folder terminated'      )
    click.echo('----------------------------------------------------------------')
