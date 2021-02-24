import click
import time
from pathlib import Path

from pysmFISH.logger_utils import selected_logger
from pysmFISH.utils import create_dir

from pysmFISH.configuration_files import create_processing_env_config_file
from pysmFISH.configuration_files import create_general_analysis_config_file



@click.group('processing_env')
def setup_processing_env():
    """
    group of commands used to setup a processing env
    """
    pass

@setup_processing_env.command('setup_de_novo_env')

@click.option('--processing_folder_path', type=str, help='Path to the folder \
                where the experiments will be processed')

def de_novo_env(processing_folder_path:str):
    logger = selected_logger()
    click.echo('    ')
    click.echo('--------------------------------------------------------------')
    click.echo('Setup the requirements for the processing environment')
    click.echo('--------------------------------------------------------------')
    
    processing_folder_path = Path(processing_folder_path)
    config_path = processing_folder_path / 'config_db'
    codebooks_path = processing_folder_path / 'codebooks'
    probes_path = processing_folder_path / 'probes_sets'
    
    # Create the directories that will contain the files required for the processing
    create_dir(config_path)
    create_dir(codebooks_path)
    create_dir(probes_path)

    # Create the analysis master files
    create_processing_env_config_file(config_path)
    create_general_analysis_config_file(config_path)


def recreate_general_analysis_config(processing_folder_path:str):
    logger = selected_logger()
    click.echo('    ')
    click.echo('--------------------------------------------------------------')
    click.echo('Re-create general analysis config file')
    click.echo('--------------------------------------------------------------')
    
    processing_folder_path = Path(processing_folder_path)

    config_path = processing_folder_path / 'config_db'

    create_general_analysis_config_file(config_path)