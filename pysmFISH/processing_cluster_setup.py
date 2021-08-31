from typing import *
import dask
import sys
from dask_jobqueue.htcondor import HTCondorCluster
from dask.distributed import LocalCluster
from psutil import virtual_memory

from pysmFISH.logger_utils import selected_logger


def htcondor_cluster_setup(htcondor_cluster_setup: dict):
    """Utility function to start a HTCondor cluster

    Args:
        htcondor_cluster_setup (dict): cluster_setup_dictionary
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
    log_directory = htcondor_cluster_setup['logs_directory']
    cluster = HTCondorCluster(cores=cores, memory=memory, 
                        disk=disk,local_directory=local_directory,
                        log_directory=log_directory,
                        processes=1)
    logger.info(f'created cluster with {cores} cores and {memory} memory')
    # make cluster more resilient
    cluster.scheduler.allowed_failures = 1000
    return cluster

# extra = ["--lifetime", "500m",
#                         "--death_timeout","5000" ]

#  "--lifetime-stagger", "10m",
# death_timeout=5000,


def local_cluster_setup():
    """Utility to set up a dask cluster on a local computer. I will use
    all the cpus-1 and scatter the memory. In thi

    Returns:
       cluster: dask cluster
    """

    total_ram = virtual_memory()
    total_ram = total_ram.available
    cores = dask.multiprocessing.multiprocessing.cpu_count()-1

    # Calculate the total ram to use for each worker
    worker_memory_limit = 0.9
    worker_memory = (total_ram*worker_memory_limit)/cores

    #cores = 5
    # worker_memory = 10000000000
    # cluster = LocalCluster(n_workers=cores, threads_per_worker=1, memory_limit=worker_memory)
    cluster = LocalCluster(n_workers=cores, memory_limit=worker_memory)


    return cluster


# TODO Run experiment with fixed size cluster (No adaptive)

def start_processing_env(processing_env_config:Dict):
    """Function to start the processing env. In the current setup
    is set up to run on the local computer or in a HPC cluster 
    manged by HTCondor. The max number of jobs in htcondor is hardcoded to 15
    and the cluster is adaptive.

    Args:
        processing_env_config (Dict): Parameters that define the 
                    cluster characteristics.

    """
    
    logger = selected_logger()
    processing_engine = processing_env_config['processing_engine']
    if processing_engine == 'htcondor':
        cluster = htcondor_cluster_setup(processing_env_config)
        cluster.scale(jobs=1)
        # Always put a minimum to avoid the cluster to shut down
        minimum_jobs = 1
        maximum_jobs = 12 #15
        # cluster.adapt(minimum_jobs=minimum_jobs,maximum_jobs=maximum_jobs)
        cluster.scale(maximum_jobs)
        logger.info(f'adaptive dask cluster with {minimum_jobs} minimum jobs')
        return cluster
    elif processing_engine == 'local':
        cluster = local_cluster_setup()
        return cluster
    else:
        logger.error(f'the processing engine is not defined check the name')
        sys.exit(f'the processing engine is not defined check the name')
