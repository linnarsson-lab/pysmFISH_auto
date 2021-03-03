import click
import time
from pathlib import Path

from pysmFISH.logger_utils import selected_logger


from flows.flow_human_embryo_single_fov import flow_human_embryo



@click.group('flows_runner')
def flows_runner():
    """
    group of commands used to setup a processing env
    """
    pass

@flows_runner.command('eel-human-embryo')

@click.option('--experiment_fpath', type=str, help='Path to the folder \
                where the experiments will be processed')
@click.option('--run_type', type=str, help='Type of run \
                select between new/re-run')
@click.option('--parsing_type', type=str, help='key to select the type of data parsing to run \
                - original: parse the data out from robofish system \
                - reparsing_from_processing_folder: parse the raw data stored in the \
                                experiment folder in the processing HD \
                - reparsing_from_storage: parse the raw data stored in the \
                                experiment folder in the storage HD \
                - no_parsing: skip parsing step')


def run_eel_human_embryo(experiment_fpath:str,run_type:str,parsing_type:str):
    logger = selected_logger()
    start_time = time.time()
    click.echo('    ')
    click.echo('--------------------------------------------------------------')
    click.echo('Run eel human embryo processing')
    click.echo('--------------------------------------------------------------')
    
    flow_human_embryo(experiment_fpath,run_type,parsing_type)

    click.echo('    ')
    click.echo('--------------------------------------------------------------')
    click.echo(f'processing wall time: {(time.time()-start_time)/60} min')
    click.echo('--------------------------------------------------------------')