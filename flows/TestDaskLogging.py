import dask
import logging
import time

from dask_jobqueue.htcondor import HTCondorCluster
from dask.distributed import Client, LocalCluster


def logger_test():
    # Redirect warnings to logger
    logging.captureWarnings(True)
    # Create logger
    logger = logging.getLogger("distributed.worker")
    logger.setLevel(logging.INFO)
    return logger

def test_parallel(a):
    logger = logger_test()
    logger.info(f'rattolomeo baisizzi')
    time.sleep(10)
    logger.info(f'processing ended')


def main():
    cores = 1
    memory = '10GB'
    disk = '1GB'
    local_directory = '/tmp'
    death_timeout = 2000
    cluster_size = 2
    logs = 'debug'
    log_directory = '/fish/work_std/test_dask_htcondor_logging'

    # logs=logs,

    list_all_process = range(10)
    cluster = HTCondorCluster(cores=cores, memory=memory, disk=disk,local_directory=local_directory,
                death_timeout=death_timeout,log_directory=log_directory,processes=1)
    cluster.adapt(minimum_jobs=cluster_size) 
    
    client = Client(cluster,asynchronous=True)
    
    L = client.map(test_parallel,list_all_process)
    # all_futures = []
    # for idx in list_all_process:
    #     future = dask.delayed(test_parallel)(idx)
    #     all_futures.append(future)
    
    # logger.info(f'all futures have been collected')
    # _ = dask.compute(*all_futures)

    client.gather(L)

    client.close()
    cluster.close()


if __name__ == "__main__":
    main()