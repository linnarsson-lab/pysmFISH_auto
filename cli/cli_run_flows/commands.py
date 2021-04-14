import click
import time
from pathlib import Path

from pysmFISH.logger_utils import selected_logger


from flows.flow_human_embryo_single_fov import flow_human_embryo
from flows.flow_human_adult_single_fov import flow_human_adult


@click.group('flows_runner')
def flows_runner():
    """
    group of commands used to setup a processing env
    """
    pass

@flows_runner.command('general-flow')

@click.option('--experiment_fpath', type=str, help='Path to the folder \
                where the experiments will be processed')
@click.option('--run_type', type=str,default='new', help='Type of run \
                select between new/re-run')
@click.option('--parsing_type', type=str, default='original', help='key to select the type of data parsing to run \
                - original: parse the data out from robofish system \
                - reparsing_from_processing_folder: parse the raw data stored in the \
                                experiment folder in the processing HD \
                - reparsing_from_storage: parse the raw data stored in the \
                                experiment folder in the storage HD \
                - no_parsing: skip parsing step')

@click.option('--fresh_nuclei_processing', type=bool,default=False, help='Preprocess \
                the fresh nuclei staining used in eel for cell segmentation')
@click.option('--raw_data_folder_storage_path', type=str,default='/fish/rawdata', 
                help='Path to the folder where the output will be stored')
@click.option('--dataset_folder_storage_path', type=str,default='/fish/fish_datasets', 
                help='Path to the folder where the stored the parsed raw data')
@click.option('--save_intermediate_steps', type=bool,default=False, 
                help='Set to True if you want to save the intermediate steps of the \
                    processing')
@click.option('--store_dataset', type=bool,default=False, 
                help='Set to True if you want to save the parsed raw data')



def run_general_flow(experiment_fpath:str,
                    run_type:str,
                    parsing_type:str,
                    fresh_nuclei_processing:bool = False,
                    raw_data_folder_storage_path:str = '/fish/rawdata',
                    dataset_folder_storage_path:str = '/fish/fish_datasets',
                    save_intermediate_steps:bool = False,
                    store_dataset:bool = True):
    
    logger = selected_logger()
    start_time = time.time()
    click.echo('    ')
    click.echo('--------------------------------------------------------------')
    click.echo('Start processing')
    click.echo('--------------------------------------------------------------')
    
    flow_human_embryo(experiment_fpath,run_type,parsing_type)

    click.echo('    ')
    click.echo('--------------------------------------------------------------')
    click.echo(f'processing wall time: {(time.time()-start_time)/60} min')
    click.echo('--------------------------------------------------------------')