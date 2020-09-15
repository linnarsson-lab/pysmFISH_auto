from typing import *
import os
import zarr
import shutil
import yaml
import datetime
import xarray as xr
from pathlib import Path
from collections import OrderedDict

from prefect import task
from prefect.engine import signals

from pysmFISH.logger_utils import prefect_logging_setup


# to avoid reference for nested structures
# https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml (comment)
yaml.SafeDumper.ignore_aliases = lambda *args : True


# @task(name='check_completed_transfer_to_monod')
def check_completed_transfer_to_monod(path_tmp_storage_server:str, flag_file_key:str):
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
    flag_file_key: str
        string that define the flag_files. The flag key should not have _auto_
        in the name.


    Returns:
        experiment_path: Posix
            experiment path in the tmp folder
    """

    logger = prefect_logging_setup('check_completed_transfer_to_monod')
    path_tmp_storage_server = Path(path_tmp_storage_server)
    flag_file_key_general = '*_auto_' + flag_file_key
    flag_file_key = '_' + flag_file_key
    flagged_files_list = list(path_tmp_storage_server.glob(flag_file_key_general))
    

    if flagged_files_list:
        flagged_file_path = flagged_files_list[0]
        experiment_name = (flagged_file_path.name).split(flag_file_key)[0]
        os.remove(flagged_file_path)
        logger.info(f'{experiment_name} ready to be processed')
        return flagged_file_path.parent / experiment_name
    else:
        logger.info(f'No new experiments to be processed')
        skip_signal = signals.SKIP("No new experiments to be processed")
        skip_signal.flag = True
        skip_signal.value = None
        raise skip_signal



@task(name='move_data')
def move_data(path_source_location:str,path_destination:str, flag_file_key:str):
    """
    Function used to transfer the files to another location

    Args:
        path_source_location: str
            path to the data to be moved
        path_destination: str
            path to the destination of the transfer
        flag_file_key: str
        string that define the flag_files. The flag key should not have _auto_
        in the name.


    Returns:
        experiment_path: Posix
            path to the data after transferring
    """

    logger = prefect_logging_setup('transfer_data')
    path_source_location = Path(path_source_location)
    path_destination = Path(path_destination)

    try:
        os.stat(path_source_location)
    except:
        logger.error(f' The {path_source_location} directory is missing')
        fail_signal = signals.FAIL('The source directory is missing')
        fail_signal.flag = True
        fail_signal.value = None
        raise fail_signal
    
    if os.stat(path_destination):
        os.stat(path_destination)
        shutil.move(path_source_location.as_posix(),path_destination.as_posix())
        tag_file_name = path_destination / (path_source_location.stem + '_' + flag_file_key)
        open(tag_file_name,'w').close()
        logger.info(f'data moved from {path_source_location} to {path_destination}')
    else:
        logger.info(f' The {path_destination} directory is missing')
        fail_signal = signals.FAIL('The destination directory is missing')
        fail_signal.flag = True
        fail_signal.value = None
        raise fail_signal






@task(name='create_folder_structure')
def create_folder_structure(experiment_fpath:str):
    """
    Function used to create the folder structure where to sort the files
    generated by the machines and the saving the data created during the
    processing. It creates the backbone structure common to all analysis

    original_robofish_logs: contains all the original robofish logs.
	extra_files: contains the extra files acquired during imaging.
	extra_processing_data: contains extra files used in the analysis 
												like the dark images for flat field correction.
    pipeline_config: contains all the configuration files.
    raw_data: contains the renamed .nd2 files and the corresponding 
						pickle configuration files. It is the directory that is 
						backed up on the server.
	output_figures: contains the reports and visualizations
    notebooks: will contain potential notebooks used for processing the data
    probes: will contains the fasta file with the probes used in the experiment
    tmp: save temporary data
    

    Args:
        experiment_fpath: str
            folder path of the experiment
    """
    logger = prefect_logging_setup('created_raw_tmp_dir')
    experiment_fpath = Path(experiment_fpath)
    folders_list = ['raw_data',
                    'original_robofish_logs',
                    'extra_processing_data',
                    'extra_files',
                    'pipeline_config',
                    'output_figures',
                    'notebooks',
                    'probes',
                    'tmp']
    for folder_name in folders_list:
        try:
            os.stat(experiment_fpath / folder_name )
            logger.info(f'{folder_name} already exist')
        except FileNotFoundError:
            os.mkdir(experiment_fpath / folder_name)
            os.chmod(experiment_fpath / folder_name,0o777)


@task(name='upload-extra-files')
def collect_extra_files(experiment_fpath:str, experiment_info:Dict):
    """
    Function used to collect extra files required for the processing
    - instrument_name_dark_img for field correction
    - codebook if the experiment is barcoded

    """
    logger = prefect_logging_setup('collect_extra_files')
    experiment_fpath = Path(experiment_fpath)

    try:
        machine = experiment_info['Machine']
    except NameError:
        machine = 'NOT_DEFINED'

    dark_img_fpath = experiment_fpath.parent / 'dark_imgs' / (experiment_info['Machine'] + '_dark_img.npy')
    try:
        shutil.copy(dark_img_fpath, (experiment_fpath / 'extra_processing_data'))
    except FileNotFoundError:
        logger.error('missing dark image')
        raise signals.FAIL('missing dark image')
    else:
        processing_env_config_fpath = experiment_fpath.parent / 'config_db' / 'processing_env_config.yaml'
        try:
            shutil.copy(processing_env_config_fpath, (experiment_fpath / 'pipeline_config'))
        except FileNotFoundError:
            logger.error('missing pipeline env config file')
            raise signals.FAIL('missing pipeline env config file')
        else:
            probes_fpath = experiment_fpath.parent / 'probes_sets' / (experiment_info['Probes'] + '.fasta')
            try:
                shutil.copy(probes_fpath, (experiment_fpath / 'probes'))
            except FileNotFoundError:
                logger.error('missing probes set file')
                raise signals.FAIL('missing probes set file')
            else:
                if 'barcoded' in experiment_info['Experiment_type']:
                    codebooks_folder = experiment_fpath.parent / 'codebooks'
                    codebook_code = experiment_info['Codebook']
                    
                    codebook_fpath = codebooks_folder / (codebook_code + '_codebook.parquet')
                    
                    # Create codebook folder in the experiment folder
                    try:
                        os.stat(experiment_fpath / 'codebook' )
                        logger.info(f'codebook folder already exist')
                    except FileNotFoundError:
                        os.mkdir(experiment_fpath / 'codebook')
                        os.chmod(experiment_fpath / 'codebook',0o777)
                    try:
                        shutil.copy(codebook_fpath, (experiment_fpath / 'codebook'))
                    except FileNotFoundError:
                        logger.error('codebook is missing')
                        raise signals.FAIL('codebook is missing')


@task(name = 'load_data_array')
def load_data_array(input_tuple):
    """
    Function used to load the images out memory as xarray data array
    
    Args:
        input_tuple: tuple
        contains the path to the zarr file to process (str) and the
        number of the fov to analyse

    Return:
        img_data_array: xarray data array
        data array with the image and all the required metadata
    """
    logger = prefect_logging_setup('load-data-array')
    zarr_store, fov_num = input_tuple

    try:
        dataset = xr.open_zarr(zarr_store)
    except:
        logger.error(f'cannot load the dataset, check the zarr store {zarr_store}')
        signals.FAIL(f'cannot load the dataset, check the zarr store {zarr_store}')
    else:
        try:
            img_data_array = dataset[str(fov_num)]
        except:
            logger.error(f'cannot load the data array, check if fov {fov_num} is missing')
            signals.FAIL(f'cannot load the data array, check if fov {fov_num} is missing')
        else:
            return img_data_array


@task(name='create_empty_zarr_file')
def create_empty_zarr_file(experiment_fpath:str,tag:str)-> str:
    """
    Function that create and empty zarr file 

    Args:
        experiment_fpath: str
            location of the folder to be processed
        tag: str
            tag to add to the file name
    Return:
        empty_fpath: str
            path of the created file

    """
    
    logger = prefect_logging_setup(f'create empty {tag} file')
    experiment_fpath = Path(experiment_fpath)
    experiment_name = experiment_fpath.stem
    zarr_fpath = experiment_fpath / (experiment_name + '_' + tag + '.zarr')
    
    dataset = xr.Dataset()
    try:
        dataset.to_zarr(zarr_fpath, mode='w', consolidated=True)
    except:
        logger.error(f'cannot create {tag} file')
        signals.FAIL(f'cannot create {tag} file')
    else:
        return zarr_fpath



@task(name = 'sort-data-folder')
def sort_data_folder(experiment_fpath:str,experiment_info:Dict):
    """
    Function used to sort the data in the experiment folder in the
    subfolders
    
    Args:
        experiment_fpath: str
            location of the folder to be processed
        experiment_info: ordered dict
            ordered dict with the parsed info generated by the instrument
    
    """
    logger = prefect_logging_setup('sort-data-folder')
    experiment_fpath = Path(experiment_fpath)

    # Move the robofish generated logs
    robofish_logs = experiment_fpath.glob('*.log')
    for log in robofish_logs:
        shutil.move(log.as_posix(), (experiment_fpath / 'original_robofish_logs').as_posix())

    # Move the crosses images
    crosses_files = experiment_fpath.glob('*_cross.nd2')
    for cross in crosses_files:
        shutil.move(cross.as_posix(), (experiment_fpath / 'extra_files').as_posix())

    # Move the coords files
    coords_files = experiment_fpath.glob('*.xml')
    for coord in coords_files:
        shutil.move(coord.as_posix(), (experiment_fpath / 'extra_files').as_posix())

    # Move the fresh nuclei file for eel
    if 'eel' in experiment_info['Experiment_type']:
        nuclei_files = list(experiment_fpath.glob('*ChannelDAPI*'))
        if len(nuclei_files):
            try:
                os.stat(experiment_fpath / 'fresh_nuclei' )
                logger.info(f'fresh_nuclei already exist')
            except FileNotFoundError:
                os.mkdir(experiment_fpath / 'fresh_nuclei')
                os.chmod(experiment_fpath / 'fresh_nuclei',0o777)
            
            for nuclei in nuclei_files:
                shutil.move(nuclei.as_posix(), (experiment_fpath / 'fresh_nuclei').as_posix())
        else:
            logger.info(f'The experiment do not have images of fresh nuclei')


@task(name = 'consolidate-metadata')
def consolidate_zarr_metadata(parsed_raw_data_fpath:str):
    """
    Function to consolidate all the zarr metadata in one unique
    json file for eady indexing and searching

    Args:
        parsed_raw_data_fpath: str
            path to the file with all the parsed images
    
    Returns:
        consolidated_grp: zarr group
            zarr groups instance with the consolidated metadata
    """

    logger = prefect_logging_setup(f'consolidate-metadata')
    
    try:
        store = zarr.DirectoryStore(parsed_raw_data_fpath)
        consolidated_grp = zarr.consolidate_metadata(store)
    except:
        logger.error(f'cannot consolidate metadata of the parsed zarr file')
        signals.FAIL(f'cannot consolidate metadata of the parsed zarr file')
    
    else:
        return consolidated_grp