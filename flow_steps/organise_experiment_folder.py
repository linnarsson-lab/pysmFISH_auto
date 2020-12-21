from pysmFISH.utils import check_ready_experiments
from pysmFISH.configuration_files import load_experiment_config_file
from pysmFISH.utils import create_folder_structure
from pysmFISH.utils import organise_files_tmp
from pysmFISH.utils import sort_data_into_folders

def organise_experiment_folder_standard(processing_hd_location:str, flag_file_key:str):
    """
    Function to organise the experimental folder and the files to process.

    Args:
        processing_hd_location: str
            path of the HD with the experiments to process
        flag_file_key: str
            key used to determine if an experiment has been completely
            transferred

    """
    # Determine if the experiment has been completely transferred to the ssd drive
    experiment_fpath = check_ready_experiments(processing_hd_location,flag_file_key)

    # Load experiment configuration file generated by robofish machines
    experiment_info = load_experiment_config_file(experiment_fpath)
    
    # Create the processing folders structure
    folders = create_folder_structure(experiment_fpath)

    # Organize the dark_images, codebooks and probes
    # Temporary approach, will be replaced by shoji
    # org_files = organise_files_tmp()
    # organisation = org_files(experiment_fpath,experiment_info)
    # organisation.set_upstream(folders)

    sorted_data = sort_data_into_folders(experiment_fpath,experiment_info)

    return experiment_fpath