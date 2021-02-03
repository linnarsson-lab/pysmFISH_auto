import click
import time
from pathlib import Path

from flows.data_organization import experiment_transfer_to_processing_hd

from pysmFISH.logger_utils import selected_logger

@click.group('data-organization')
def data_organization_utilities():
    """
    group of commands used to organize the data and moved them
    in different HD for processing or storage
    """
    pass

@data_organization_utilities.command('exp-transfer-processing-hd')

@click.option('--path_source', type=str, help='Path to the folder to scan')

@click.option('--path_destination', type=str, help='Path to the folder/hd where the \
                the processing is run')
@click.option('--flag_file_key', type=str, default= '_transfer_to_monod_completed.txt' ,
                help='String used to identify the  \
                the files used  indicate that the experiment have been fully \
                transfer from the machine \
                ex. _transfer_to_monod_completed.txt')

@click.option('--completion_pattern', type=str, default= 'processing_completed.txt' ,
                help='String used to identify the files  \
                used to confirm that the processing is completed and a new experiment \
                can be transfered \
                ex. processing_completed.txt')


@click.option('--min_free_space', type=int, default=1000, 
                help='Minimum storage space in the  \
                destination required to trigger transfer expressed in Gb \
                ex. 1000')

@click.option('--monitor_frequency', type=int, default=43200, help='Scanning frequency of the source  \
                folder in seconds  \
                ex. 43200')

                                        


def exp_proc_transf(path_source,path_destination,
                    flag_file_key,completion_pattern,min_free_space,monitor_frequency):
    
    logger = selected_logger()
    click.echo('    ')
    click.echo('--------------------------------------------------------------')
    click.echo('Transfer the experiments to the processing HD')
    click.echo('--------------------------------------------------------------')
    experiment_transfer_to_processing_hd(path_source,path_destination,
                    flag_file_key,completion_pattern,min_free_space,monitor_frequency)
