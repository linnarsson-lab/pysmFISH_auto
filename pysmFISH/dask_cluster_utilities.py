from typing import *
import logging
import dask
from dask_jobqueue.htcondor import HTCondorCluster
from dask.distributed import Client, LocalCluster
from psutil import virtual_memory
from typing import *

class htcondor_cluster_setup():

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

    def __init__(self,htcondor_cluster_setup: dict):

        self.logger = logging.getLogger(__name__)

        self.cores = htcondor_cluster_setup['cores']
        self.memory = htcondor_cluster_setup['memory']
        self.disk = htcondor_cluster_setup['disk']
        self.local_directory = htcondor_cluster_setup['local_directory']

        self.cluster = HTCondorCluster(cores=self.cores, memory=self.memory, disk=self.disk,local_directory=self.local_directory,death_timeout=2400)
        self.logger.info(f'created cluster with {self.cores} cores and {self.memory} memory')
        self.logger.info(f'URL to dask dashboard {self.cluster.dashboard_link}')


class local_cluster_setup():
    """
    Utility to set up a dask cluster on a local computer. I will use
    all the cpus-1 and scatter the memory
    """

    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.total_ram = virtual_memory()
        self.total_ram = self.total_ram.total
        self.cores = dask.multiprocessing.multiprocessing.cpu_count()-2

        # Calculate the total ram to use for each worker
        self.worker_memory_limit = 0.9
        self.worker_memory = (self.total_ram*self.worker_memory_limit)/self.cores
        self.cluster = LocalCluster(n_workers=self.cores, memory_limit=self.worker_memory)
