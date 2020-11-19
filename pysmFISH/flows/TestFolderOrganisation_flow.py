import prefect
from prefect import task, Flow, Parameter, flatten, unmapped
from prefect.engine.executors import DaskExecutor
from prefect.environments import LocalEnvironment
from prefect import Task
from prefect.environments.storage import Local

from pysmFISH.utilities_tasks import check_ready_experiments
from pysmFISH.utilities_tasks import create_folder_structure
from pysmFISH.utilities_tasks import organise_files
from pysmFISH.configuration_files_tasks import load_experiment_config_file


with Flow("folder-organization",environment=LocalEnvironment(DaskExecutor(address='tcp://193.10.16.58:32833')),
            storage=Local(directory='/home/simone/tmp_code/flows')) as flow:
   
    # Required processing parameters
    processing_hd_location = Parameter('processing_hd_location',default='/wsfish/smfish_ssd')
    flag_file_key = Parameter('flag_file_key',default='transfer_to_monod_completed.txt')

    # Determine if the experiment has been completely transferred to the ssd drive
    ready_experiment = check_ready_experiments()
    experiment_fpath = ready_experiment(processing_hd_location,flag_file_key)

    # Create the processing folders structure
    create_folders = create_folder_structure()
    folders = create_folders(experiment_fpath)

    folders = create_folder_structure(experiment_fpath)

    # Load experiment configuration file generated by robofish machines
    # load_exp_cfg = load_experiment_config_file()
    # experiment_info = load_exp_cfg(experiment_fpath)

  

    # Organize the data in the different folders
    # org_files = organise_files()
    # organisation = org_files(experiment_fpath,experiment_info)
    # organisation.set_upstream(folders)

flow.register(project_name="test")
# flow.run_agent()