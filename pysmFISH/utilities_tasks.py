from typing import *
import os
import zarr
import shutil
import yaml
import datetime
import xarray as xr
import numpy as np
from pathlib import Path
from collections import OrderedDict

import prefect
from prefect import Task
from prefect import task
from prefect.engine import signals

from pysmFISH.logger_utils import prefect_logging_setup


# to avoid reference for nested structures
# https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml (comment)
yaml.SafeDumper.ignore_aliases = lambda *args : True

class free_space(Task):
    """
    A task used to determine if there is enough space in the
    HD where the experiment will be processed

    Args:
    -----
        hd_path:str
            pathway of the target HD where the processing will be run
        min_free_space:int
            minimum space required for processing in Gb (ex: 1000)
    
    """
    def run(self, hd_path:str, min_free_space:int):
        """
        Function used to determine if there is enough free space in the
        HD where the experiment will be processed

        The hard coded minimum current space is 

        Args:
        -----
            hd_path:str
                pathway of the target HD where the processing will be run
            min_free_space:int
                minimum space required for processing in Gb (ex: 1000)
        
        Returns:
        --------
            True/False: bool
                True if there is enough free space for running the experiment
        """
        total, used, free = shutil.disk_usage(hd_path)
        free_space_giga = free // (2**30)
        if free_space_giga <= min_free_space:
            self.logger.info(f'Free space in the HD: {free_space_giga} Gb data cannot be transferred,\
                            not enough space on the HD')
            skip_signal = signals.FAIL("Not enogh space in the processing hd")
            skip_signal.flag = True
            skip_signal.value = None
            raise skip_signal
        else:
            self.logger.info(f'Free space in the HD: {free_space_giga} Gb data can be transferred')
            return True

class check_ready_experiments(Task):
    """
    class to scan the folder where the data are transferred from the machines.
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
        experiment_path: str
            experiment path in the tmp folder
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self, path_tmp_storage_server:str, flag_file_key:str):
        path_tmp_storage_server = Path(path_tmp_storage_server)
        flag_file_key_general = '*_auto_' + flag_file_key
        flag_file_key = '_' + flag_file_key
        flagged_files_list = list(path_tmp_storage_server.glob(flag_file_key_general))
        
        if flagged_files_list:
            flagged_file_path = flagged_files_list[0]
            experiment_name = (flagged_file_path.name).split(flag_file_key)[0]
            os.remove(flagged_file_path)
            self.logger.info(f'{experiment_name} ready to be processed')
            return str(flagged_file_path.parent / experiment_name)
        else:
            self.logger.error(f'No new experiments to be processed')
            skip_signal = signals.FAIL(f"No new experiments to be processed")
            skip_signal.flag = True
            skip_signal.value = None
            raise skip_signal



class transfer_data(Task):
    """
    Class used to move files to another location

    Args:
        path_source_location: str
            path to the data to be moved
        path_destination: str
            path to the destination of the transfer
        flag_file_key: str
        string that define the flag_files. The flag key should not have _auto_
        in the name.
    """

    def run(self, path_source_location:str,path_destination:str, flag_file_key:str):
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
        """

        path_source_location = Path(path_source_location)
        path_destination = Path(path_destination)

        try:
            os.stat(path_source_location)
        except:
            self.logger.error(f' The {path_source_location} directory is missing')
            fail_signal = signals.FAIL('The source directory is missing')
            fail_signal.flag = True
            fail_signal.value = None
            raise fail_signal
        else:
            try:
                os.stat(path_destination)
            except:
                self.logger.info(f' The {path_destination} directory is missing')
                fail_signal = signals.FAIL('The destination directory is missing')
                fail_signal.flag = True
                fail_signal.value = None
                raise fail_signal
            else:
                shutil.move(path_source_location.as_posix(),path_destination.as_posix())
                tag_file_name = path_destination / (path_source_location.stem + '_' + flag_file_key)
                open(tag_file_name,'w').close()
                self.logger.info(f'data moved from {path_source_location} to {path_destination}')


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



# @task(name='create_folder_structure')
# def create_folder_structure(experiment_fpath:str):
#     """
#     Function used to create the folder structure where to sort the files
#     generated by the machines and the saving the data created during the
#     processing. It creates the backbone structure common to all analysis

#     original_robofish_logs: contains all the original robofish logs.
# 	extra_files: contains the extra files acquired during imaging.
# 	extra_processing_data: contains extra files used in the analysis 
# 												like the dark images for flat field correction.
#     pipeline_config: contains all the configuration files.
#     raw_data: contains the renamed .nd2 files and the corresponding 
# 						pickle configuration files. It is the directory that is 
# 						backed up on the server.
# 	output_figures: contains the reports and visualizations
#     notebooks: will contain potential notebooks used for processing the data
#     probes: will contains the fasta file with the probes used in the experiment
#     tmp: save temporary data
    

#     Args:
#         experiment_fpath: str
#             folder path of the experiment
#     """
#     logger = prefect_logging_setup('created_raw_tmp_dir')
#     experiment_fpath = Path(experiment_fpath)
#     folders_list = ['raw_data',
#                     'original_robofish_logs',
#                     'extra_processing_data',
#                     'extra_files',
#                     'pipeline_config',
#                     'output_figures',
#                     'notebooks',
#                     'probes',
#                     'tmp']
#     for folder_name in folders_list:
#         try:
#             os.stat(experiment_fpath / folder_name )
#             logger.info(f'{folder_name} already exist')
#         except FileNotFoundError:
#             os.mkdir(experiment_fpath / folder_name)
#             os.chmod(experiment_fpath / folder_name,0o777)


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


@task(name = 'open-consolidated-metadata')
def open_consolidated_metadata(parsed_raw_data_fpath:str):
    logger = prefect_logging_setup(f'consolidate-metadata')
    
    try:
        store = zarr.DirectoryStore(parsed_raw_data_fpath)
    except:
        logger.error(f'the metadata are not consolidated')
        signals.FAIL(f'the metadata are not consolidated')
    else:
        consolidated_grp = zarr.open_consolidated(store)
        return consolidated_grp



# @task(name = 'load-raw-images-and-filtering-attrs')
def load_raw_images(zarr_grp_name:str,parsed_raw_data_fpath:str)->np.ndarray:
    """
    Function used to load a raw image and metadata from the 
    parsed raw file and the attrs for the filtering
        parsed_raw_data_fpath: str
            fpath to zarr store containing the parsed raw images
        zarr_grp_name: str
            fpath to the group to process. The group contain the raw images and the 
            corresponding metadata

            grp = experiment_name_channel_fov_X
                dataset = raw_data_fov_X

    """
    logger = prefect_logging_setup(f'consolidate-metadata')
    st = zarr.DirectoryStore(parsed_raw_data_fpath)
    root = zarr.group(store=st,overwrite=False)

    metadata = root[zarr_grp_name].attrs
    img = root[zarr_grp_name][metadata['fov_name']][...]
    return (img, dict(metadata))

@task(name='sort-images-according-processing-type')
def sorting_grps(grps, experiment_info, analysis_parameters):
    """
    Function used to separate the group names according to 
    the processing type

        grps = zarr.hierarchy.Group
            consolidate hierarchy zarr groups for fast access to metadata
        experiment_info: ordered dict
            ordered dict with the parsed info generated by the instrument
        analysis_parameters: dict
            dict with all the parameters to use for analysis  
    """

    fish_grp = []
    beads_grp = []
    staining_grp = []
    for name, grp in grps.items():
        if grp.attrs['stitching_channel'] == grp.attrs['channel']:
            beads_grp.append(name)
        elif '_ST' in grp.attrs['target_name']:
            staining_grp.append(name)
        else:
            fish_grp.append(name)

    logger = prefect_logging_setup(f'sort-groups')
    if experiment_info['Stitching_type'] == 'small-beads':
        beads_selected_parameters  = analysis_parameters['small-beads']
        beads_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']
    elif experiment_info['Stitching_type'] == 'large-beads':
        beads_selected_parameters  = analysis_parameters['small-beads']
        beads_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']
    elif experiment_info['Stitching_type'] == 'both-beads':
        beads_selected_parameters  = analysis_parameters['small-beads']
        beads_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']
    else:
        logger.error(f'Wrong stitching type in the experiment file')
        err = signals.FAIL(f'Wrong stitching type in the experiment file')
        raise err
    
    fish_selected_parameters = analysis_parameters['fish']
    fish_selected_parameters['BarcodesExtractionResolution'] = analysis_parameters['BarcodesExtractionResolution']
    fish_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']

    staining_selected_parameters = analysis_parameters['staining']
    staining_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']

    return fish_grp, fish_selected_parameters, beads_grp, beads_selected_parameters, staining_grp, staining_selected_parameters


@task(name='sort-images-according-processing-type')
def sorting_grps_fov(grps, experiment_info, analysis_parameters):
    """
    Function used to separate the group names according to 
    the processing type and by fov to reduce the number of calls
    to the writing.

    Args:
    -----
        grps = zarr.hierarchy.Group
            consolidate hierarchy zarr groups for fast access to metadata
        experiment_info: ordered dict
            ordered dict with the parsed info generated by the instrument
        analysis_parameters: dict
            dict with all the parameters to use for analysis  
    """

    fish_grp = {}
    beads_grp = {}
    staining_grp = []
    for name, grp in grps.items():
        if grp.attrs['stitching_channel'] == grp.attrs['channel']:
            if grp.attrs['channel'] not in beads_grp.keys():
                beads_grp[grp.attrs['channel']] = {}
            if grp.attrs['fov_num'] not in beads_grp[grp.attrs['channel']].keys():
                beads_grp[grp.attrs['channel']][grp.attrs['fov_num']] = []
                beads_grp[grp.attrs['channel']][grp.attrs['fov_num']].append(name)
            else:
                beads_grp[grp.attrs['channel']][grp.attrs['fov_num']].append(name)
        elif '_ST' in grp.attrs['target_name']:
            staining_grp.append(name)
        else:
            if grp.attrs['channel'] not in fish_grp.keys():
                fish_grp[grp.attrs['channel']] = {}
            if  grp.attrs['fov_num'] not in fish_grp[grp.attrs['channel']].keys():
                fish_grp[grp.attrs['channel']][grp.attrs['fov_num']] = []
                fish_grp[grp.attrs['channel']][grp.attrs['fov_num']].append(name)
            else:
                fish_grp[grp.attrs['channel']][grp.attrs['fov_num']].append(name)

    fish_list = []
    beads_list = []
    for channel in fish_grp.keys():
        for fov, names in fish_grp[channel].items():
            fish_list.append(names)
    
    if len(fish_grp.keys()) > 1:
        fish_list = [el for sgr in fish_list for el in sgr]
    
    for channel in beads_grp.keys():
        for fov, names in beads_grp[channel].items():
            beads_list.append(names)
    
    if len(beads_grp.keys()) > 1:
        beads_list = [el for sgr in beads_list for el in sgr]


    logger = prefect_logging_setup(f'sort-groups')
    if experiment_info['Stitching_type'] == 'small-beads':
        beads_selected_parameters  = analysis_parameters['small-beads']
        beads_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']
    elif experiment_info['Stitching_type'] == 'large-beads':
        beads_selected_parameters  = analysis_parameters['small-beads']
        beads_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']
    else:
        logger.error(f'Wrong stitching type in the experiment file')
        err = signals.FAIL(f'Wrong stitching type in the experiment file')
        raise err
    
    fish_selected_parameters = analysis_parameters['fish']
    fish_selected_parameters['BarcodesExtractionResolution'] = analysis_parameters['BarcodesExtractionResolution']
    fish_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']

    staining_selected_parameters = analysis_parameters['staining']
    staining_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']

    logger.info(f'fish group {fish_list}')
    logger.info(f'beads group {beads_list}')

    return fish_list, fish_selected_parameters, beads_list, beads_selected_parameters, staining_grp, staining_selected_parameters