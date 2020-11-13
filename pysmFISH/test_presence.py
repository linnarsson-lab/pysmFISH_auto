from pysmFISH.dask_cluster_utilities_tasks import start_processing_env
from pysmFISH.configuration_files_tasks import load_processing_env_config_file, load_experiment_config_file

config_db_fpath ='/wsfish/smfish_ssd/config_db'
experiment_fpath = '/wsfish/smfish_ssd/LBEXP20201014_EEL_Mouse_2420um_auto'

processing_env_config = load_processing_env_config_file(config_db_fpath)

experiment_info = load_experiment_config_file(experiment_fpath)

# QC for experiment info if contains all the info


# Activate processing cluster
cluster = start_processing_env(processing_env_config,experiment_info)

address = cluster.scheduler_address