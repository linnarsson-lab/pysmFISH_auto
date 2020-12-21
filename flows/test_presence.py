import dask
import time
import pickle

from pysmFISH.dask_cluster_utilities_tasks import start_processing_env
from pysmFISH.configuration_files_tasks import load_processing_env_config_file, load_experiment_config_file

from dask_jobqueue.htcondor import HTCondorCluster
from dask.distributed import Client, LocalCluster


def test_parallel(a):
    print('test_parallel')
    time.sleep(10)

if __name__ == "__main__":
    
    cores = 1
    memory = '1GB'
    disk = '1GB'
    local_directory = '/tmp'
    death_timeout_list = [2400,90,1000,60,30,5000,None]
    # death_timeout = None

    cluster_size_list = [50,100,150]
    total_time = 0
    step_length = 30
    max_length = 300
    for cluster_size in cluster_size_list:
        for death_timeout in death_timeout_list:
            print(f'running cluster {cluster_size}, death_timeout {death_timeout}')
            cluster = HTCondorCluster(cores=cores, memory=memory, disk=disk,local_directory=local_directory,death_timeout=death_timeout)
            all_steps_increase_cluster_size = []

            # print(f'cluster observed: {len(cluster.observed)}')
            # print(f'cluster plan: {len(cluster.plan)}')

            all_steps_increase_cluster_size.append(len(cluster.observed))
            cluster.scale(jobs=cluster_size) 

            # print(f'cluster observed: {len(cluster.observed)}')
            # print(f'cluster plan: {len(cluster.plan)}')

            all_steps_increase_cluster_size.append(len(cluster.observed))

            while True:
                # print(f'cluster observed: {len(cluster.observed)}')
                # print(f'cluster plan: {len(cluster.plan)}')
                time.sleep(step_length)
                total_time += step_length
                all_steps_increase_cluster_size.append(len(cluster.observed))
                if total_time > max_length:
                    break

            cluster.close()

            fpath = '/home/simone/test_cluster_loading/cluster_scale_'+ str(cluster_size) + '_jobs_step_length_' + str(step_length) + '_death_time_out_' + str(death_timeout) +'.pkl'
            pickle.dump(all_steps_increase_cluster_size,open(fpath,'wb'))
            print(f'cluster size {cluster_size}, death_timeout {death_timeout} all step increase {all_steps_increase_cluster_size}')

            time.sleep(200)
#     all_futures = []
#     for idx in np.range(100):
#         future = dask.delayed(test_parallel)(idx)
#         all_futures.append(future)


# dask.compute(*all_futures)
config_db_fpath ='/wsfish/smfish_ssd/config_db'
experiment_fpath = '/wsfish/smfish_ssd/LBEXP20201207_EEL_HE_test2'

processing_env_config = load_processing_env_config_file(config_db_fpath)

experiment_info_loader = load_experiment_config_file()
experiment_info = experiment_info_loader.run(experiment_fpath)

# # QC for experiment info if contains all the info

# # Activate processing cluster
cluster = start_processing_env(processing_env_config,experiment_info)

address = cluster.scheduler_address

