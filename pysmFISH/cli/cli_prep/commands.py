import click
import time
from pathlib import Path

from pysmFISH.logger_utils import selected_logger
from pysmFISH.utils import create_dir

from pysmFISH.configuration_files import create_general_analysis_config_file
from pysmFISH.barcodes_analysis import simplify_barcodes_reference



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
    create_general_analysis_config_file(config_path)


@setup_processing_env.command('recreate_general_analysis_config')

@click.option('--processing_folder_path', type=str, help='Path to the folder \
                where the experiments will be processed')

def recreate_general_analysis_config(processing_folder_path:str):
    logger = selected_logger()
    click.echo('    ')
    click.echo('--------------------------------------------------------------')
    click.echo('Re-create general analysis config file')
    click.echo('--------------------------------------------------------------')
    
    processing_folder_path = Path(processing_folder_path)

    config_path = processing_folder_path / 'config_db'

    create_general_analysis_config_file(config_path)



@setup_processing_env.command('convert_codebook')

@click.option('--codebook_path', type=str, help='Path to the folder \
                with the codebook to be parsed. The codebook excel file \
                should be in the config_db folder')

def convert_codebook(codebook_path:str):
    logger = selected_logger()
    click.echo('    ')
    click.echo('--------------------------------------------------------------')
    click.echo('Convert the codebook in a lighter version')
    click.echo('--------------------------------------------------------------')
    parser = simplify_barcodes_reference(codebook_path)
    parser.convert_barcode()
