from pathlib import Path
import os
import time
import shutil

from pysmFISH.utils import free_space
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