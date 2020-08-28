"""
create configuration files needed to run the pipeline

"""
from typing import *
import yaml
from pathlib import Path
from collections import OrderedDict

from prefect import task
from prefect.engine import signals

from pysmFISH.logger_utils import prefect_logging_setup


# to avoid reference for nested structures
# https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml (comment)
yaml.SafeDumper.ignore_aliases = lambda *args : True

def create_processing_env_config_file(experimnent_fpath:str):
    """
    This function creates a configuration files with the parameters requested for the 
    setup of the processing cluster. 
    The current file is defaulted to htcondor

    Args:
        experiment_fpath: str
            path to the experiment that will be processed
    """

    processing_env_config_fpath = experimnent_fpath / 'pipeline_config' / 'processing_env_config.yaml'
    processing_env_config = OrderedDict()

    processing_env_config['htcondor'] = {}

    processing_env_config['htcondor']['smfish-serial'] = {}
    processing_env_config['htcondor']['smfish-serial']['cores'] = 1
    processing_env_config['htcondor']['smfish-serial']['memory'] = '10GB'
    processing_env_config['htcondor']['smfish-serial']['disk'] = '0.1GB'
    processing_env_config['htcondor']['smfish-serial']['local_directory'] = '/tmp'

    processing_env_config['htcondor']['smfish-barcoded'] = {}
    processing_env_config['htcondor']['smfish-barcoded']['cores'] = 1
    processing_env_config['htcondor']['smfish-barcoded']['memory'] = '10GB'
    processing_env_config['htcondor']['smfish-barcoded']['disk'] = '0.1GB'
    processing_env_config['htcondor']['smfish-barcoded']['local_directory'] = '/tmp'

    processing_env_config['htcondor']['eel-barcoded'] = {}
    processing_env_config['htcondor']['eel-barcoded']['cores'] = 1
    processing_env_config['htcondor']['eel-barcoded']['memory'] = '10GB'
    processing_env_config['htcondor']['eel-barcoded']['disk'] = '0.1GB'
    processing_env_config['htcondor']['eel-barcoded']['local_directory'] = '/tmp'

    with open(processing_env_config_fpath, 'w') as new_config:
            yaml.safe_dump(dict(processing_env_config), new_config,default_flow_style=False,sort_keys=False)

@task(name='create_analysis_config_file')
def create_analysis_config_file(experiment_fpath:str, experiment_info:Dict):
    """
    This function creates a basic standard configuration files with the parameters used for running
    analysis
    The selection of the paramenters is dependent of the machine used for image acquisition. If the
    machine is not defined in the experiment_info dictionary a generic processing set up is defined

    Args:
        experiment_fpath: str
            path to the experiment that will be processed

        experiment_info: dict
            dictionary containing all the info generated by robofish
    """
    logger = prefect_logging_setup('create-analasys-config-file')
    analysis_config_fpath = experiment_fpath / 'pipeline_config' / 'analysis_config.yaml'
    analysis_config = OrderedDict()

    try:
        machine = experiment_info['Machine']
    except NameError:
        machine = 'NOT_DEFINED'
    
    analysis_config['fish_signal'] = {}
    analysis_config['fish_signal']['preprocessing'] = {}
    analysis_config['fish_signal']['preprocessing']['fish'] = {}
    analysis_config['fish_signal']['preprocessing']['small-beads-registration'] = {}
    analysis_config['fish_signal']['preprocessing']['large-beads-registration'] = {}

    analysis_config['fish_signal']['rounds_registration'] = {}
    analysis_config['fish_signal']['rounds_registration']['reference_hybridization'] = 1

    analysis_config['fish_signal']['counting'] = {}
    analysis_config['fish_signal']['counting']['fish'] = {}
    analysis_config['fish_signal']['counting']['small-beads-registration'] = {}
    analysis_config['fish_signal']['counting']['large-beads-registration'] = {}

    analysis_config['fish_signal']['barcodes_extraction'] = {}

    analysis_config['staining'] = {}
    analysis_config['staining']['preprocessing'] = {}

    if machine == 'ROBOFISH1':

        analysis_config['fish_signal']['preprocessing']['fish']['flat_field_kernel'] = (2,100,100)
        analysis_config['fish_signal']['preprocessing']['fish']['filtering_small_kernel'] = (1,8,8)
        analysis_config['fish_signal']['preprocessing']['fish']['filtering_laplacian_kernel'] = (0.2,0.5,0.5)
        analysis_config['fish_signal']['preprocessing']['fish']['large_obj_removal_percentile'] = 99
        analysis_config['fish_signal']['preprocessing']['fish']['large_obj_removal_min_obj_size'] = 50
        analysis_config['fish_signal']['preprocessing']['fish']['large_obj_removal_selem'] = 3 

    
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['flat_field_kernel'] = (2,100,100)
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['filtering_small_kernel'] = (1,8,8)
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['filtering_laplacian_kernel'] = (0.2,0.5,0.5)
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['large_obj_removal_percentile'] = 99
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['large_obj_removal_min_obj_size'] = 50
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['large_obj_removal_selem'] = 3 
        
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['flat_field_kernel'] = (2,100,100)
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['filtering_small_kernel'] = (1,8,8)
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['filtering_laplacian_kernel'] = (0.2,0.5,0.5)
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['large_obj_removal_percentile'] = 99
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['large_obj_removal_min_obj_size'] = 50
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['large_obj_removal_selem'] = 3

        analysis_config['fish_signal']['counting']['fish']['min_distance'] = 2
        analysis_config['fish_signal']['counting']['fish']['min_obj_size'] = 2
        analysis_config['fish_signal']['counting']['fish']['max_obj_size'] = 200
        analysis_config['fish_signal']['counting']['fish']['num_peaks_per_label'] = 1

        analysis_config['fish_signal']['counting']['small-beads-registration']['min_distance'] = 2
        analysis_config['fish_signal']['counting']['small-beads-registration']['min_obj_size'] = 2
        analysis_config['fish_signal']['counting']['small-beads-registration']['max_obj_size'] = 200
        analysis_config['fish_signal']['counting']['small-beads-registration']['num_peaks_per_label'] = 1

        analysis_config['fish_signal']['counting']['large-beads-registration']['min_distance'] = 2
        analysis_config['fish_signal']['counting']['large-beads-registration']['min_obj_size'] = 2
        analysis_config['fish_signal']['counting']['large-beads-registration']['max_obj_size'] = 200
        analysis_config['fish_signal']['counting']['large-beads-registration']['num_peaks_per_label'] = 1

        analysis_config['fish_signal']['barcodes_extraction']['pxl_hood_size'] = 3

        analysis_config['staining']['preprocessing']['flat_field_kernel'] = (2,100,100)

        if 'barcoded' in experiment_info['experiment_type']:
            analysis_config['fresh_nuclei'] = {}
            analysis_config['fresh_nuclei']['preprocessing'] = {}
            analysis_config['fresh_nuclei']['preprocessing']['large_kernel_size'] = (5,50,50)

    elif machine == 'ROBOFISH2':

        analysis_config['fish_signal']['preprocessing']['fish']['flat_field_kernel'] = (2,100,100)
        analysis_config['fish_signal']['preprocessing']['fish']['filtering_small_kernel'] = (1,8,8)
        analysis_config['fish_signal']['preprocessing']['fish']['filtering_laplacian_kernel'] = (0.2,0.5,0.5)
        analysis_config['fish_signal']['preprocessing']['fish']['large_obj_removal_percentile'] = 99
        analysis_config['fish_signal']['preprocessing']['fish']['large_obj_removal_min_obj_size'] = 50
        analysis_config['fish_signal']['preprocessing']['fish']['large_obj_removal_selem'] = 3 

    
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['flat_field_kernel'] = (2,100,100)
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['filtering_small_kernel'] = (1,8,8)
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['filtering_laplacian_kernel'] = (0.2,0.5,0.5)
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['large_obj_removal_percentile'] = 99
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['large_obj_removal_min_obj_size'] = 50
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['large_obj_removal_selem'] = 3 
        
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['flat_field_kernel'] = (2,100,100)
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['filtering_small_kernel'] = (1,8,8)
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['filtering_laplacian_kernel'] = (0.2,0.5,0.5)
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['large_obj_removal_percentile'] = 99
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['large_obj_removal_min_obj_size'] = 50
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['large_obj_removal_selem'] = 3

        analysis_config['fish_signal']['counting']['fish']['min_distance'] = 2
        analysis_config['fish_signal']['counting']['fish']['min_obj_size'] = 2
        analysis_config['fish_signal']['counting']['fish']['max_obj_size'] = 200
        analysis_config['fish_signal']['counting']['fish']['num_peaks_per_label'] = 1

        analysis_config['fish_signal']['counting']['small-beads-registration']['min_distance'] = 2
        analysis_config['fish_signal']['counting']['small-beads-registration']['min_obj_size'] = 2
        analysis_config['fish_signal']['counting']['small-beads-registration']['max_obj_size'] = 200
        analysis_config['fish_signal']['counting']['small-beads-registration']['num_peaks_per_label'] = 1

        analysis_config['fish_signal']['counting']['large-beads-registration']['min_distance'] = 2
        analysis_config['fish_signal']['counting']['large-beads-registration']['min_obj_size'] = 2
        analysis_config['fish_signal']['counting']['large-beads-registration']['max_obj_size'] = 200
        analysis_config['fish_signal']['counting']['large-beads-registration']['num_peaks_per_label'] = 1

        analysis_config['fish_signal']['barcodes_extraction']['pxl_hood_size'] = 3

        analysis_config['staining']['preprocessing']['flat_field_kernel'] = (2,100,100)

        if 'barcoded' in experiment_info['experiment_type']:
            analysis_config['fresh_nuclei'] = {}
            analysis_config['fresh_nuclei']['preprocessing'] = {}
            analysis_config['fresh_nuclei']['preprocessing']['large_kernel_size'] = (5,50,50)

    elif machine == 'NOT_DEFINED':

        analysis_config['fish_signal']['preprocessing']['fish']['flat_field_kernel'] = (2,100,100)
        analysis_config['fish_signal']['preprocessing']['fish']['filtering_small_kernel'] = (1,8,8)
        analysis_config['fish_signal']['preprocessing']['fish']['filtering_laplacian_kernel'] = (0.2,0.5,0.5)
        analysis_config['fish_signal']['preprocessing']['fish']['large_obj_removal_percentile'] = 99
        analysis_config['fish_signal']['preprocessing']['fish']['large_obj_removal_min_obj_size'] = 50
        analysis_config['fish_signal']['preprocessing']['fish']['large_obj_removal_selem'] = 3 

    
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['flat_field_kernel'] = (2,100,100)
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['filtering_small_kernel'] = (1,8,8)
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['filtering_laplacian_kernel'] = (0.2,0.5,0.5)
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['large_obj_removal_percentile'] = 99
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['large_obj_removal_min_obj_size'] = 50
        analysis_config['fish_signal']['preprocessing']['small-beads-registration']['large_obj_removal_selem'] = 3 
        
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['flat_field_kernel'] = (2,100,100)
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['filtering_small_kernel'] = (1,8,8)
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['filtering_laplacian_kernel'] = (0.2,0.5,0.5)
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['large_obj_removal_percentile'] = 99
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['large_obj_removal_min_obj_size'] = 50
        analysis_config['fish_signal']['preprocessing']['large-beads-registration']['large_obj_removal_selem'] = 3

        analysis_config['fish_signal']['counting']['fish']['min_distance'] = 2
        analysis_config['fish_signal']['counting']['fish']['min_obj_size'] = 2
        analysis_config['fish_signal']['counting']['fish']['max_obj_size'] = 200
        analysis_config['fish_signal']['counting']['fish']['num_peaks_per_label'] = 1

        analysis_config['fish_signal']['counting']['small-beads-registration']['min_distance'] = 2
        analysis_config['fish_signal']['counting']['small-beads-registration']['min_obj_size'] = 2
        analysis_config['fish_signal']['counting']['small-beads-registration']['max_obj_size'] = 200
        analysis_config['fish_signal']['counting']['small-beads-registration']['num_peaks_per_label'] = 1

        analysis_config['fish_signal']['counting']['large-beads-registration']['min_distance'] = 2
        analysis_config['fish_signal']['counting']['large-beads-registration']['min_obj_size'] = 2
        analysis_config['fish_signal']['counting']['large-beads-registration']['max_obj_size'] = 200
        analysis_config['fish_signal']['counting']['large-beads-registration']['num_peaks_per_label'] = 1

        analysis_config['fish_signal']['barcodes_extraction']['pxl_hood_size'] = 3

        analysis_config['staining']['preprocessing']['flat_field_kernel'] = (2,100,100)

        if 'barcoded' in experiment_info['experiment_type']:
            analysis_config['fresh_nuclei'] = {}
            analysis_config['fresh_nuclei']['preprocessing'] = {}
            analysis_config['fresh_nuclei']['preprocessing']['large_kernel_size'] = (5,50,50)
    try:
        with open(analysis_config_fpath, 'w') as new_config:
                yaml.safe_dump(dict(analysis_config), new_config,default_flow_style=False,sort_keys=False)
    except:
        logger.error(f'cannot save the analysis_config_file')
        signals.FAIL(f'cannot save the analysis_config_file')


# @task(name='load-processing-env-config')
def load_processing_env_config_file(experiment_fpath:str):
    logger = prefect_logging_setup('load-processing-env-config')
    experiment_fpath = Path(experiment_fpath)
    processing_env_config_fpath = experiment_fpath / 'pipeline_config' / 'processing_env_config.yaml'
    try:
        processing_env_config = OrderedDict(yaml.safe_load(open(processing_env_config_fpath, 'rb')))
    except (FileExistsError,NameError,FileNotFoundError) as e:
        logger.debug(f'{processing_env_config_fpath} missing')
        try:
            create_processing_env_config_file(experiment_fpath)
        except:
            logger.error('cannot create the processing config file')
            signals.FAIL('cannot create the processing config file')
        else:
            processing_env_config = OrderedDict(yaml.safe_load(open(processing_env_config_fpath, 'rb')))
            return processing_env_config
    else:
        return processing_env_config


@task(name='load-analysis-config')
def load_analysis_config_file(experiment_fpath:str):
    logger = prefect_logging_setup('load-analysis-config')
    experiment_fpath = Path(experiment_fpath)
    analysis_config_fpath = experiment_fpath / 'pipeline_config' / 'analysis_config.yaml'
    try:
        analysis_config = OrderedDict(yaml.safe_load(open(analysis_config_fpath, 'rb')))
    except (FileExistsError,NameError,FileNotFoundError) as e:
        logger.debug(f'{analysis_config_fpath} missing')
        signals.FAIL(f'{analysis_config_fpath} missing')
    else:
        return analysis_config


def load_transferring_config(config_db_fpath:str)->Dict:
    """
    Function used to load the parameters used for transfering
    the data from the storage to the processing HD

    Args:
        config_db_fpath; str
            path to the folder containing the data_transfer_config.yaml
    """
    
    logger = prefect_logging_setup(f'load-transfer-config')
    config_db_fpath = Path(config_db_fpath)

    transfer_config_fpath = config_db_fpath / 'data_transfer_config.yaml'
    try:
         transfer_config = OrderedDict(yaml.safe_load(open(transfer_config_fpath, 'rb')))
    except (FileExistsError,NameError,FileNotFoundError) as e:
        logger.debug(f'{transfer_config_fpath} missing')
        signals.FAIL(f'{transfer_config_fpath} missing')
    else:
        return  transfer_config

# @task(name='load_experiments_info_file')
def load_experiment_config_file(experiment_fpath:str):
    """
    Function that load the experiment general information generated by the
    machines.

    Args:
        experiment_fpath: str
            location of the folder to be processed
    Return:
        experiment_info: ordered dict
            ordered dict with the parsed info generated by the instrument

    """
    
    logger = prefect_logging_setup('load_experiments_info_file')

    experiment_fpath = Path(experiment_fpath)
    experiment_name = experiment_fpath.stem
    search_key = experiment_name + '_config.yaml'
    
    try:
        experiment_info_fpath = list(experiment_fpath.glob(search_key))[0]
    except:
        logger.error(f'No experiment info file in {experiment_fpath}')
        fail_signal = signals.FAIL("No experiment info file in the folder")
        fail_signal.flag = True
        fail_signal.value = None
        raise fail_signal
    
    try:
        experiment_info = OrderedDict(yaml.safe_load(open(experiment_info_fpath, 'rb')))
        return experiment_info
    except:
        logger.error(f'Experiment info file has the wrong name in {experiment_fpath}')
        fail_signal = signals.FAIL("Experiment info file has the wrong name in the folder")
        fail_signal.flag = True
        fail_signal.value = None
        raise fail_signal

