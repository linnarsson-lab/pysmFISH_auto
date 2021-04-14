import click
import time
from pathlib import Path

from pysmFISH.directory_scanner import experiment_processing_runner

from pysmFISH.logger_utils import selected_logger

@click.group('experiment-folder-monitoring')
def folder_monitoring_utilities():
    """
    group of commands used to monitor the folder to process
    experiments
    """
    pass

@folder_monitoring_utilities.command('exp-directory-monitor')


@click.option('--path_source', type=str, help='Path to the folder to scan')

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

                                        


def dir_monitoring(path_source,
                    flag_file_key,completion_pattern,min_free_space,monitor_frequency):
    
    logger = selected_logger()
    click.echo('    ')
    click.echo('--------------------------------------------------------------')
    click.echo('monitor directory')
    click.echo('--------------------------------------------------------------')
    experiment_processing_runner(path_source,path_destination,
                    flag_file_key,completion_pattern,min_free_space,monitor_frequency)
