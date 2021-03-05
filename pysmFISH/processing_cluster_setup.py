from typing import *
import yaml
import dask
import sys
from dask_jobqueue.htcondor import HTCondorCluster
from dask.distributed import Client, LocalCluster
from psutil import virtual_memory
from pathlib import Path


from pysmFISH.logger_utils import selected_logger

def htcondor_cluster_setup(htcondor_cluster_setup: dict):

    """
    Utility class used to create a dask cluster in HTCondor.

    Parameters:
    -----------
        cluster_setup_dictionary
        dictionary with the info for the cluster setup
        {
            cores: number of cores
            memory: RAM to use
            disk: space on disk to use
            local_directory: Dask worker local directory for file spilling.        
        }

    """

    logger = selected_logger()
    cores = htcondor_cluster_setup['cores']
    memory = htcondor_cluster_setup['memory']
    disk = htcondor_cluster_setup['disk']
    local_directory = htcondor_cluster_setup['local_directory']
    log_directory = htcondor_cluster_setup['experiment_fpath'] + '/logs/'
    cluster = HTCondorCluster(cores=cores, memory=memory, 
                        disk=disk,local_directory=local_directory,
                        log_directory=log_directory,
                        death_timeout=5000)
    logger.info(f'created cluster with {cores} cores and {memory} memory')
    return cluster

def local_cluster_setup():
    
    """
    Utility to set up a dask cluster on a local computer. I will use
    all the cpus-1 and scatter the memory
    """

    total_ram = virtual_memory()
    total_ram = total_ram.available
    cores = dask.multiprocessing.multiprocessing.cpu_count()-2

    # Calculate the total ram to use for each worker
    # worker_memory_limit = 0.9
    # worker_memory = (total_ram*worker_memory_limit)/cores

    cores = 5
    worker_memory = 10000000000
    cluster = LocalCluster(n_workers=cores, threads_per_worker=1, memory_limit=worker_memory)


    return cluster


def start_processing_env(processing_env_config:Dict,experiment_info:Dict,experiment_fpath:str):
    """
    Function to start the processing env. In the current setup
    is set up to run on the local computer or in a HPC cluster 
    manged by HTCondor

    Args;
        processing_env_config: Dict
            Dict with the parameters for starting the cluster
        experiment_info: Dict
            dictionary with all the info describing the experiment
        experiment_fpath: str
            path to the experiment to process
    Return:
        cluster: dask-cluster-obj
                cluster responsible for the data processing

    """
    logger = selected_logger()
    processing_engine = processing_env_config['processing_engine']
    if processing_engine == 'htcondor':
        experiment_type = experiment_info['Experiment_type']
        cluster_config_parameters = processing_env_config[processing_engine][experiment_type]
        cluster_config_parameters['experiment_fpath'] = experiment_fpath
        cluster = htcondor_cluster_setup(cluster_config_parameters)
        # cluster.scale(jobs=50)
        minimum_jobs = 5
        maximum_jobs = 800
        cluster.adapt(minimum_jobs=minimum_jobs,maximum_jobs=maximum_jobs)
        logger.info(f'adaptive dask cluster with {minimum_jobs} minimum jobs')
        return cluster
    elif processing_engine == 'local':
        cluster = local_cluster_setup()
        return cluster
    else:
        logger.error(f'the processing engine is not defined check the name')
        sys.exit(f'the processing engine is not defined check the name')


def start_transfering_env(processing_hd_location:str):
    """
    Function to start the environment used for transfering the files that
    will be processed in a cluster. 

    It requires a data_transfer_config.yaml in a folder called config_db in
    the folder where the experiments will be processed

    Args;
        processing_hd_location: str
            path to the folder where the experiments are processed

    """

    logger = selected_logger()

    processing_hd_location = Path(processing_hd_location)
    transfer_env_config_fpath = processing_hd_location / 'config_db' / 'data_transfer_config.yaml'

    try:
        transfer_env_config = yaml.safe_load(open(transfer_env_config_fpath,'rb'))
    except (FileExistsError,NameError) as e:
        logger.error(f'missing transfer env configuration file')
    else:
        processing_engine = transfer_env_config['processing_engine']
        if processing_engine == 'htcondor':
            cluster_config_parameters = transfer_env_config[processing_engine]
            cluster = htcondor_cluster_setup(cluster_config_parameters)
            return cluster
        elif processing_engine == 'local':
            cluster = local_cluster_setup()
            return cluster
        else:
            logger.error(f'the processing engine is not defined check the name')
            sys.exit(f'the processing engine is not defined check the name')
