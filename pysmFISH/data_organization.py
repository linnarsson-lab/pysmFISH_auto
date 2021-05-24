from pathlib import Path
import os
import sys
import time
import shutil

from pysmFISH.utils import free_space, create_dir, copytree
from pysmFISH.logger_utils import selected_logger


def transfer_experiment(flagged_file_path:str, path_destination:str, flag_file_key:str, completion_fpath:str):
    """
    Function to scan the folder where the data are transferred from the machines.
    It looks for files that are generated upon transfer completion (flag_file).
    To keep thing simple and not overload the system only one of the flag_files
    identified in the scanning of the folder will be processed.
    
    NB: In order to process only the experiments designed for the pipeline 
    the folder name finish with auto ExpName_auto

    Args: 
    path_tmp_storage_server: str
        path where the data are transferred from the microscopes
    path_destination: str
            path to the destination of the transfer
    flag_file_key: str
        string that define the flag_files. The flag key should not have _auto_
        in the name.

    Returns:
        experiment_path: str
            experiment path in the tmp folder
    """
    
    logger = selected_logger()
    flagged_file_path = Path(flagged_file_path)
    path_destination = Path(path_destination)

    try:
        experiment_name = (flagged_file_path.name).split(flag_file_key)[0]
        logger.info(f'{experiment_name} ready to be processed')
    except:
        logger.error(f'No experiment folder matching the flagged file')
        logger.error(f'The flagged file is removed')
        os.remove(flagged_file_path.as_posix())
    else:
        experiment_fpath = flagged_file_path.parent / experiment_name
        
        try:
            os.stat(path_destination)
        except:
            logger.info(f' The {path_destination} directory is missing')
        else:
            if 'auto' in experiment_fpath.stem:
                shutil.move(experiment_fpath.as_posix(),path_destination.as_posix())
                tag_file_name = path_destination / (experiment_fpath.stem + flag_file_key)
                open(tag_file_name,'w').close()
                os.remove(flagged_file_path)
                os.remove(completion_fpath)
                logger.info(f'data moved from {experiment_fpath.parent} to {path_destination}')
            else:
                logger.info(f' auto tag not in experiment name {experiment_fpath.stem}')
                logger.error(f'The flagged file is removed')
                os.remove(flagged_file_path)


def experiment_transfer_to_processing_hd(path_source:str,path_destination:str,
                                         flag_file_key:str='_transfer_to_monod_completed.txt',
                                         completion_pattern:str='processing_completed.txt',
                                         min_free_space:int=1000,
                                         monitor_frequency:int=43200):

    logger = selected_logger()
    path_source = Path(path_source)
    path_destination = Path(path_destination)
    try:

        while True:
            files = list(path_source.glob('*' + flag_file_key))
            files.sort(key=lambda x: os.path.getmtime(x))
            logger.info(f'processing files {files}')
            logger.info(f'experiments to process {len(files)}')
            # Check that the previous file has been processed
            completion_file = list(path_source.glob(completion_pattern))

            if completion_file:
                logger.info(f'completion file present')
                # Check if there is enough space to run transfer and processing
                allowed_space = free_space(path_destination, min_free_space)
                for flagged_file_path in files:
                    if allowed_space:
                        transfer_experiment(flagged_file_path, path_destination, flag_file_key,completion_file[0])
                        allowed_space = free_space(path_destination, min_free_space)
            else:
                logger.info(f'completion file is missing')

            time.sleep(monitor_frequency)

    except KeyboardInterrupt:
        logger.info(f'scanning interrupted')
        pass
        # log output



def reorganize_processing_dir(experiment_fpath:str,
                            storage_fpath:str,
                            store_dataset=False,
                            dataset_storage_fpath=None,
                            results_storage_fpath=None):
    """
    Function to transfer and copy the data in the data storage HD

    Args:
    ----
    experiment_fpath: str
        path to the folder where the processing will be run

    storage_fpath: str
        path of the experiment in the storage HD that contains the raw data
    """
    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)
    storage_fpath = Path(storage_fpath)

    experiment_name = experiment_fpath.name



    # Remove the tmp data
    folders_to_remove = ['tmp','logs']
    for folder_name in folders_to_remove:
        shutil.rmtree((experiment_fpath / folder_name).as_posix(),ignore_errors=False)
    
    # removed or store parsed data
    parsed_data_fpath = experiment_fpath / (experiment_name + '_img_data.zarr')
    filtered_data_fpath = experiment_fpath / (experiment_name + '_preprocessed_img_data.zarr')
    if store_dataset:
        if dataset_storage_fpath:
            try:
                dataset_file = list(experiment_fpath.glob('*_img_data_dataset.parquet'))[0]
            except IndexError:
                logger.error(f"The dataset file is missing")
            else:
                dataset_file_new_name = Path(dataset_storage_fpath) / dataset_file.name
                shutil.copy2(dataset_file,dataset_file_new_name)
                dataset_storage_fpath = Path(dataset_storage_fpath) / parsed_data_fpath.name
                try:
                    shutil.move(parsed_data_fpath,dataset_storage_fpath)
                except FileNotFoundError:
                    logger.error(f"The parsed raw data zarr file is missing")
                else:
                    filtered_storage_fpath = Path(filtered_data_fpath) / filtered_data_fpath.name
                    try:
                        shutil.move(filtered_storage_fpath,filtered_storage_fpath)
                    except FileNotFoundError:
                        logger.error(f"The filtered data zarr file is missing")
        else:
            logger.error(f'the dataset_storage_fpath is missing')

    else:
        if results_storage_fpath:
            results_fpath = results_storage_fpath / ('results_' + experiment_fpath.stem)
            create_dir(results_fpath)
            results_original = experiment_fpath / 'results'
            output_fig_original = experiment_fpath / 'output_figures'

            results_new = results_fpath / 'results'
            output_fig_new = results_fpath / 'output_figures'

            shutil.move(results_original,results_new)
            shutil.move(output_fig_original,output_fig_new)
        
        shutil.rmtree(parsed_data_fpath.as_posix())


    # Transfer the remaining folder to the raw data folder
    new_folder_location = storage_fpath / experiment_fpath.stem
    try:
        _ = shutil.move(experiment_fpath.as_posix(),new_folder_location.as_posix())
    except OSError:
        logger.error(f'the experiment folder cannot be moved')
   

def transfer_data_to_storage(experiment_fpath:str,storage_fpath:str,):
    """
    Function to transfer and copy the data in the data storage HD

    Args:
    ----
    experiment_fpath: str
        path to the folder where the processing will be run

    storage_fpath: str
        path of the experiment in the storage HD that contains the raw data
    """
    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)
    storage_fpath = Path(storage_fpath)

    folders_to_transfer = ['raw_data','extra_files','original_robofish_logs']
    folders_to_copy = ['extra_processing_data',
                        'pipeline_config',
                        'codebook',
                        'probes',
                        'fresh_nuclei']

    try:
        config_file_fpath = list(experiment_fpath.glob('*.yaml'))[0]
    except:
        logger.error(f'The experiment .yaml config file is missing')
        sys.exit(f'The experiment .yaml config file is missing')
    else:
        experiment_name = experiment_fpath.stem
        stored_experiment_fpath = storage_fpath / experiment_name
        create_dir(stored_experiment_fpath)  

        for folder_name in folders_to_transfer:
            src = experiment_fpath / folder_name
            dst = stored_experiment_fpath / folder_name
            shutil.move(src.as_posix(),dst.as_posix())

        for folder_name in folders_to_copy:
            src = (experiment_fpath / folder_name).as_posix()
            dst = (stored_experiment_fpath / folder_name).as_posix()
            copytree(src,dst)
        
        # Copy configuration file
        dst_config = (stored_experiment_fpath / config_file_fpath.name).as_posix()
        _ = shutil.copy2(config_file_fpath,dst_config)

        # move the pkl files
        pkl_fpaths = experiment_fpath.glob('*.pkl')
        for pkl_path in pkl_fpaths:
            dst = stored_experiment_fpath / pkl_path.name
            shutil.move(pkl_path.as_posix(),dst.as_posix())


def transfer_files_from_storage(storage_experiment_fpath:str,experiment_fpath:str):
    """
    Function used to transfer the files required for data processing when
    the raw images are collected from the storage directory

    Args:
    ----
    storage_experiment_fpath: str
        path of the experiment in the storage HD that contains the raw data
    
    experiment_fpath: str
        path to the folder where the processing will be run
    """
    logger = selected_logger()

    storage_experiment_fpath = Path(storage_experiment_fpath)
    logger.info(f'storage directory {storage_experiment_fpath.as_posix()}')
    experiment_fpath = Path(experiment_fpath)
    logger.info(f'target directory {experiment_fpath.as_posix()}')

    folders_to_copy=[
        'extra_processing_data',
        'microscope_tiles_coords',
        'codebook',
        'probes',
        'pipeline_config'
    ]
    
    list_dir = [x[0] for x in os.walk(storage_experiment_fpath)]
    try:
        config_file_fpath = list(storage_experiment_fpath.glob('*.yaml'))[0]
    except:
        logger.error(f'The experiment .yaml config file is missing')
        sys.exit(f'The experiment .yaml config file is missing')
    else:
        for folder_name in folders_to_copy:
            src = (storage_experiment_fpath / folder_name).as_posix()
            if src in list_dir:
                dst = (experiment_fpath / folder_name).as_posix()
                copytree(src,dst)
            else:
                logger.error(f'missing {folder_name}')
        dst_config = (experiment_fpath / config_file_fpath.name).as_posix()
        _ = shutil.copy2(config_file_fpath,dst_config)
