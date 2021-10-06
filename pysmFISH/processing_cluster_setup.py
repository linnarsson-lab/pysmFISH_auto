from typing import *
import dask
import os
import signal
import sys
from dask_jobqueue.htcondor import HTCondorCluster
from dask.distributed import LocalCluster, SSHCluster
from psutil import virtual_memory

from pysmFISH.logger_utils import selected_logger

def htcondor_cluster_setup(htcondor_cluster_setup: dict):
    """Utility function to start a HTCondor cluster

    Args:
        htcondor_cluster_setup (dict): cluster_setup_dictionary
        dictionary with the info for the cluster setup

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


def local_cluster_setup(cores:int, memory:str):
    """Utility to set up a dask cluster on a local computer. I will use
    all the cpus-1 and scatter the memory. In thi

    Args:
        cores (int): number of cores of the computer to use for processing
        memory (str): memory for each core (ex. 5GB) 
    Returns:
       cluster: dask cluster
    """

    # total_ram = virtual_memory()
    # total_ram = total_ram.available
    
    # cores = dask.multiprocessing.multiprocessing.cpu_count()-1

    # Calculate the total ram to use for each worker
    # worker_memory_limit = 0.9
    # worker_memory = (total_ram*worker_memory_limit)/cores

    #cores = 5
    # worker_memory = 10000000000
    # cluster = LocalCluster(n_workers=cores, threads_per_worker=1, memory_limit=worker_memory)
    # cluster = LocalCluster(n_workers=cores, memory_limit=worker_memory)
    cluster = LocalCluster(n_workers=cores, memory_limit=memory, processes=True,threads_per_worker=1)



    return cluster


def unmanaged_cluster_setup(htcondor_cluster_setup:Dict):
    """Create and start a unmanaged cluster. The cluster is
    created by ssh into the machines. Because there is bug in
    OpenSSH in not enough to kill the cluster by using
    client.close()
    cluster.close()
    it is necessary to kill the distributed process using
    kill_distributed_process()

    Args:
        htcondor_cluster_setup (Dict): dictionary with the info for the cluster setup

    Returns:
       cluster: dask cluster
    """
    
    logger = selected_logger()
    cores = htcondor_cluster_setup['cores']
    memory = htcondor_cluster_setup['memory']
    disk = htcondor_cluster_setup['disk']
    local_directory = htcondor_cluster_setup['local_directory']
    log_directory = htcondor_cluster_setup['logs_directory']
    scheduler_port = htcondor_cluster_setup['scheduler_port']
    dashboard_port = htcondor_cluster_setup['dashboard_port']
    nprocs = htcondor_cluster_setup['nprocs']
    scheduler_address = htcondor_cluster_setup['scheduler_address']
    workers_addresses_list = htcondor_cluster_setup['workers_addresses_list']
    nthreads = htcondor_cluster_setup['nthreads']
    
    
    all_addresses = [scheduler_address] + workers_addresses_list
    
    worker_options = {"nprocs":nprocs,
                     "cores":cores,
                     "memory_limit":memory,
                     "nthreads":nthreads,
                     "local_directory":local_directory}
    
    
    cluster = SSHCluster(
        all_addresses,
        connect_options={"known_hosts": None},
        worker_options=worker_options,
        scheduler_options={"port": 0, 
                       "dashboard_address":25399}
        )

    return cluster



def kill_process(process_name:str='distributed.cli.dask.scheduler'):
    """General function used to kill a process by name directly from a
    script (https://www.geeksforgeeks.org/kill-a-process-by-name-using-python/)
    
    I needed the function to compensate for the bug in OpenSSH that
    doesn't allow the destruction of the cluster using the standard
    dask commands ( https://github.com/dask/distributed/issues/3420 ). 
    By default the distributed process will be killed.

    Args:
        process_name (str, optional): [description]. Defaults to 'distributed.cli.dask.scheduler'.
    """
   # https://github.com/dask/distributed/issues/3420  
    # Ask user for the name of process
    
    logger = selected_logger()
    try:
        # iterating through each instance of the process
        for line in os.popen("ps ax | grep " + process_name + " | grep -v grep"):
            fields = line.split()
             
            # extracting Process ID from the output
            pid = fields[0]
             
            # terminating process
            os.kill(int(pid), signal.SIGKILL)
        logger.info("Process Successfully terminated")
         
    except:
        logger.error("Error Encountered while running script")




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
        if processing_env_config['adaptive']:
            cluster.adapt(minimum_jobs=minimum_jobs,maximum_jobs=processing_env_config['maximum_jobs'])
            logger.info(f"Started adaptive cluster")
        else:
            cluster.scale(jobs=processing_env_config['maximum_jobs'])
            logger.info(f"Started non adaptive cluster")
        return cluster
    elif processing_engine == 'local':
        cluster = local_cluster_setup(processing_env_config['cores'],processing_env_config['memory'])
        return cluster
    elif processing_engine == 'unmanaged cluster':
        cluster = unmanaged_cluster_setup(processing_env_config)
        logger.info(f"Started unmanaged cluster")
        return cluster
    else:
        logger.error(f'the processing engine is not defined check the name')
        sys.exit(f'the processing engine is not defined check the name')
