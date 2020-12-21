import prefect
from prefect import task, Flow, Parameter, flatten, unmapped
from prefect.engine.executors import DaskExecutor
from prefect.environments import LocalEnvironment
from prefect import Task
from prefect.environments.storage import Local

from pysmFISH.utilities_tasks import check_ready_experiments
from pysmFISH.utilities_tasks import free_space
from pysmFISH.utilities_tasks import transfer_data


from pysmFISH.configuration_files_tasks import load_transferring_config

# transferring_config_fpath = '/wsfish/smfish_ssd/config_db'
# load_tconfig = load_transferring_config()
# transfer_config = load_tconfig.run(transferring_config_fpath)


# Load parameters that regulates data transfer 
# transferring_config_fpath = Parameter('transferring_config_fpath',default='/wsfish/smfish_ssd/config_db')
# load_tconfig = load_transferring_config()
# transfer_config = load_tconfig(transferring_config_fpath.default)

with Flow("transfering-flow",environment=LocalEnvironment(DaskExecutor(address='tcp://193.10.16.58:32833')),
            storage=Local(directory='/home/simone/tmp_code/flows')) as flow:
   
    # Required processing parameters
    path_tmp_storage_server = Parameter('path_tmp_storage_server',default='/fish/work_std')
    flag_file_key = Parameter('flag_file_key',default='transfer_to_monod_completed.txt')
    processing_hd_location = Parameter('processing_hd_location',default='/wsfish/smfish_ssd')
    minimum_hd_free_space = Parameter('minimum_hd_free_space',default= 1000)
    
    # Determine if there is enough space in the destination HD
    check_space = free_space()
    allowed_space = check_space(processing_hd_location,minimum_hd_free_space)
    
    # Determine if there is a new experiment to transfer
    ready_experiment = check_ready_experiments()
    experiment_fpath = ready_experiment(path_tmp_storage_server,flag_file_key)
    experiment_fpath.set_upstream(allowed_space)

    # Tranfer the experiment to the ssd drive
    move_experiment = transfer_data()
    move_experiment(experiment_fpath, processing_hd_location, flag_file_key)

flow.register(project_name="test")
# flow.run_agent()