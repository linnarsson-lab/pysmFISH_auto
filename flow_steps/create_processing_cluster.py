import sys

from pysmFISH.logger_utils import selected_logger

from pysmFISH.configuration_files import load_experiment_config_file
from pysmFISH.configuration_files import load_processing_env_config_file

from pysmFISH.processing_cluster_setup import start_processing_env

def create_processing_cluster(processing_env_config_fpath:str,experiment_fpath:str):
    """
    Function to create the cluster that will be used to process
    the entire pipeline

    Args:
         processing_env_config_fpath: str
           str to the location of the file with the config parameters
           for setting up the cluster
        experiment_fpath: str
            path to the experiment to process

    Return:
        cluster: dask-cluster-obj
                cluster responsible for the data processing
    """
    logger = selected_logger()

    try:
        processing_env_config = load_processing_env_config_file(processing_env_config_fpath)
    except:
        logger.error(f'error loading {processing_env_config_fpath}')
        sys.exit(f'error loading {processing_env_config_fpath}')
    else:
        try:
            experiment_info = load_experiment_config_file(experiment_fpath)
        except:
            logger.error(f'error loading {experiment_fpath}')
            sys.exit(f'error loading {experiment_fpath}')
        else:
            # Cluster setup
            cluster = start_processing_env(processing_env_config=processing_env_config, 
                                    experiment_info=experiment_info,
                                    experiment_fpath=experiment_fpath)

            return cluster