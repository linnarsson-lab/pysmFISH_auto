from typing import *
import logging
import shutil
import yaml
import sys
import os
import click
import pickle
import datetime
import nd2reader
import numpy as np
from collections import OrderedDict
from skimage import img_as_float64
from pathlib import Path
from scipy.ndimage import fourier_shift

from pysmFISH.logger_utils import selected_logger

def free_space(hd_path:str, min_free_space:int):
    """
    Function used to determine if there is enough space in the
    HD where the experiment will be processed

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
    logger = selected_logger()
    total, used, free = shutil.disk_usage(hd_path)
    free_space_giga = free // (2**30)
    if free_space_giga <= min_free_space:
        logger.info(f'Free space in the HD: {free_space_giga} Gb data cannot be transferred,\
                        not enough space on the HD')
        return False
    else:
        logger.info(f'Free space in the HD: {free_space_giga} Gb data can be transferred')
        return True



def create_folder_structure(experiment_fpath:str,run_type:str):
    """
    Class used to create the folder structure where to sort the files
    generated by the machines and the saving the data created during the
    processing. It creates the backbone structure common to all analysis

    FOLDER STRUCTURE
    - original_robofish_logs: contains all the original robofish logs.
    - extra_files: contains the extra files acquired during imaging.
    - extra_processing_data: contains extra files used in the analysis 
                                                like the dark images for flat field correction.
    - pipeline_config: contains all the configuration files.
    - raw_data: contains the renamed .nd2 files and the corresponding 
                        pickle configuration files. It is the directory that is 
                        backed up on the server.
    - output_figures: contains the reports and visualizations
    - notebooks: will contain potential notebooks used for processing the data
    - probes: will contains the fasta file with the probes used in the experiment
    - tmp: save temporary data
    
    Args:
        experiment_fpath: str
            folder path of the experiment
    """

    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)

    if run_type == 'new':

        folders_list = ['raw_data',
                        'original_robofish_logs',
                        'extra_processing_data',
                        'extra_files',
                        'pipeline_config',
                        'output_figures',
                        'notebooks',
                        'probes',
                        'tmp',
                        'logs',
                        'results',
                        'microscope_tiles_coords']
        subfolders_tmp = ['raw_counts','filtered_images','registered_counts',
                        'combined_rounds_images']
    
    else:
        folders_list = ['extra_processing_data',
                    'pipeline_config',
                    'output_figures',
                    'probes',
                    'tmp',
                    'logs',
                    'results',
                    'microscope_tiles_coords']
        subfolders_tmp = ['raw_counts','filtered_images','registered_counts',
                    'combined_rounds_images']
    
    for folder_name in folders_list:
        try:
            os.stat(experiment_fpath / folder_name )
            logger.info(f'{folder_name} already exist')
        except FileNotFoundError:
            os.mkdir(experiment_fpath / folder_name)
            os.chmod(experiment_fpath / folder_name,0o777)
    
    for folder_name in subfolders_tmp:
        try:
            os.stat(experiment_fpath / 'tmp' / folder_name )
            logger.info(f'{folder_name} already exist')
        except FileNotFoundError:
            os.mkdir(experiment_fpath / 'tmp' / folder_name)
            os.chmod(experiment_fpath / 'tmp' / folder_name,0o777)



def collect_processing_files(experiment_fpath:str, experiment_info:Dict):

    """
    Task used to sort the files in the experiment folder and
    gather the extra files required for the processing
    
    - instrument_name_dark_img for field correction
    - codebook if the experiment is barcoded

    """
    logger = selected_logger()

    experiment_fpath = Path(experiment_fpath)

    try:
        machine = experiment_info['Machine']
    except NameError:
        machine = 'NOT_DEFINED'

    # This step can be also be removed in case we won't use the processing env config 
    processing_env_config_fpath = experiment_fpath.parent / 'config_db' / 'processing_env_config.yaml'
    try:
        shutil.copy(processing_env_config_fpath, (experiment_fpath / 'pipeline_config'))
    except FileNotFoundError:
        logger.error('missing pipeline env config file')
        sys.exit('missing pipeline env config file')
    else:
        # This step will be modified when the rpobes will be stored in shoji
        probes_fpath = experiment_fpath.parent / 'probes_sets' / (experiment_info['Probe_FASTA_name'])
        try:
            shutil.copy(probes_fpath, (experiment_fpath / 'probes'))
        except FileNotFoundError:
            logger.error('missing probes set file')
            sys.exit('missing probes set file')
        else:
            # This step will be modified when the codebook will be stored in shoji
            if 'barcoded' in experiment_info['Experiment_type']:
                codebooks_folder = experiment_fpath.parent / 'codebooks'
                codebook_code = experiment_info['Codebook']
                
                codebook_fpath = codebooks_folder / codebook_code
                
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
                    sys.exit('codebook is missing')


def sort_data_into_folders(experiment_fpath:str,experiment_info:Dict):
    """
    Function used to sort the files created by the machine and contained 
    in the experiment folder in the subfolders
    
    Args:
        experiment_fpath: str
            location of the folder to be processed
        experiment_info: ordered dict
            ordered dict with the parsed info generated by the instrument
    
    """
    
    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)

    # Move the robofish generated logs
    robofish_logs = list(experiment_fpath.glob('*.log'))
    if len(robofish_logs):
        for log in robofish_logs:
            shutil.move(log.as_posix(), (experiment_fpath / 'original_robofish_logs').as_posix())
            logger.debug(f'moved {log.stem} into original_robofish_logs')
    else:
        logger.debug(f'there are no original robofish logs in the experiment folder')


    # Move the crosses images
    crosses_files = list(experiment_fpath.glob('*_cross.nd2'))
    if len(crosses_files):
        for cross in crosses_files:
            shutil.move(cross.as_posix(), (experiment_fpath / 'extra_files').as_posix())
            logger.debug(f'moved {cross.stem} into extra_files')
    else:
        logger.debug(f'there are no crosses images in the experiment folders')
    
    # Move the coords files
    coords_files = list(experiment_fpath.glob('*.xml'))
    if len(coords_files):
        for coord in coords_files:
            shutil.move(coord.as_posix(), (experiment_fpath / 'extra_files').as_posix())
            logger.debug(f'moved {coord.stem} into extra_files')
    else:
        logger.debug(f'there are no coords files in the experiment folders')
    
    # Move the fresh nuclei file for eel
    if 'eel' in experiment_info['Experiment_type']:
        nuclei_files = list(experiment_fpath.glob('*ChannelEuropium_Cy3*'))
        if len(nuclei_files):
            try:
                os.stat(experiment_fpath / 'fresh_nuclei' )
                logger.debug(f'fresh_nuclei already exist')
            except FileNotFoundError:
                os.mkdir(experiment_fpath / 'fresh_nuclei')
                os.chmod(experiment_fpath / 'fresh_nuclei',0o777)
            
            for nuclei in nuclei_files:
                shutil.move(nuclei.as_posix(), (experiment_fpath / 'fresh_nuclei').as_posix())
                logger.debug(f'moved {nuclei.stem} into fresh nuclei')
        else:
            logger.debug(f'The experiment does not have images of fresh nuclei')
        
        large_views = list(experiment_fpath.glob('*ChannelDAPI*'))
        if len(large_views):
            for large_view in large_views:
                shutil.move(large_view.as_posix(), (experiment_fpath / 'extra_files').as_posix())
                logger.debug(f'moved {large_view.stem} into extra_files')
        else:
            logger.debug(f'The experiment does not have large view images of fresh nuclei')

    # move remaining .nd2 files
    all_nd2 = list(experiment_fpath.glob('*.nd2'))
    remaining = [el for el in all_nd2 if 'Count' not in el.stem]
    if len(remaining):
        for fremaining in remaining:
            shutil.move(fremaining.as_posix(), (experiment_fpath / 'extra_files').as_posix())
            logger.debug(f'moved {fremaining.stem} into extra_files')
    else:
        logger.debug(f'There are no extra files to be moved')



def create_dark_img(experiment_fpath,experiment_info):
    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)
    experiment_name = experiment_info['EXP_name']
    machine = experiment_info['Machine']
    # Check if the dark image is already present
    try:
        pres = list((experiment_fpath / 'extra_processing_data').glob('*_dark_img.npy'))[0]
    except:
        nd2_list = list((experiment_fpath / 'extra_files').glob('*.nd2'))
        nd2_blank = [el for el in nd2_list if 'Blank' in el.stem]
        if nd2_blank:
            nd2fh = nd2reader.ND2Reader(nd2_blank[0])
            nd2fh.bundle_axes = 'zyx'
            nd2fh.iter_axes = 'z'
            # Image with single channel
            channel = nd2fh.metadata['channels'][0]
            # I can collect just one z because error of the reader
            # that created both z and t planes
            dark_img = np.median(nd2fh[0],axis=0)
            dark_img = dark_img.astype(np.uint16)
            fname = experiment_fpath / 'extra_processing_data' / (experiment_name + '_' + machine + '_dark_img.npy')
            np.save(fname, dark_img)
            logger.debug(f'Created dark image')
        else:
            logger.error(f'the Blank .nd2 for the dark image is missing')
    else:
        logger.debug(f'the dark image is already present')


def sorting_grps(grps, experiment_info, analysis_parameters):
    """
    Function used to separate the group names according to 
    the processing type
    Args:
    -----
        grps = zarr.hierarchy.Group
            consolidate hierarchy zarr groups for fast access to metadata
        experiment_info: ordered dict
            ordered dict with the parsed info generated by the instrument
        analysis_parameters: dict
            dict with all the parameters to use for analysis  
    Returns:
    --------
        fish_grp: list
            list of the groups containing fish images
        fish_selected_parameters: list
            matching parameters used for fish analysis
        beads_grp: list
            list of the groups containing beads images
        beads_selected_parameters: list
            matching parameters used for beads analysis
        staining_grp: list
            list of the groups containing IF/staining images
        staining_selected_parameters:list
            matching parameters used for IF/staining analysis
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

    if experiment_info['Stitching_type'] == 'small-beads':
        beads_selected_parameters  = analysis_parameters['small-beads']
        beads_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']
    elif experiment_info['Stitching_type'] == 'large-beads':
        beads_selected_parameters  = analysis_parameters['large-beads']
        beads_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']
    elif experiment_info['Stitching_type'] == 'both-beads':
        beads_selected_parameters  = analysis_parameters['both-beads']
        beads_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']
    elif experiment_info['Stitching_type'] == 'nuclei' :
        beads_selected_parameters  = analysis_parameters['nuclei']
        beads_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']
    else:
        logger.error(f'Wrong stitching type in the experiment file')
        sys.exit(f'Wrong stitching type in the experiment file')
    
    fish_selected_parameters = analysis_parameters['fish']
    fish_selected_parameters['BarcodesExtractionResolution'] = analysis_parameters['BarcodesExtractionResolution']
    fish_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']

    staining_selected_parameters = analysis_parameters['staining']
    staining_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']

    combined_grps ={'fish':(fish_grp, fish_selected_parameters),
                    'beads':(beads_grp, beads_selected_parameters),
                    'staining': (staining_grp, staining_selected_parameters)
                    }

    return combined_grps


def sorting_grps_for_fov_processing(consolidated_grp, experiment_info, analysis_parameters):

    registration_channel = experiment_info['StitchingChannel']
    experiment_name =  experiment_info['EXP_name']
    key = experiment_name + '_Hybridization01_' + registration_channel + '_fov_0'
    fovs = consolidated_grp[key].attrs['fields_of_view']

    if experiment_info['Stitching_type'] == 'small-beads':
            registration_selected_parameters  = analysis_parameters['small-beads']
            registration_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']
    elif experiment_info['Stitching_type'] == 'large-beads':
        registration_selected_parameters  = analysis_parameters['large-beads']
        registration_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']
    elif experiment_info['Stitching_type'] == 'both-beads':
        registration_selected_parameters  = analysis_parameters['both-beads']
        registration_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']
    elif experiment_info['Stitching_type'] == 'nuclei' :
        registration_selected_parameters  = analysis_parameters['nuclei']
        registration_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']
    else:
        logger.error(f'Wrong stitching type in the experiment file')
        sys.exit(f'Wrong stitching type in the experiment file')

    fish_selected_parameters = analysis_parameters['fish']
    fish_selected_parameters['BarcodesExtractionResolution'] = analysis_parameters['BarcodesExtractionResolution']
    fish_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']

    staining_selected_parameters = analysis_parameters['staining']
    staining_selected_parameters['RegistrationReferenceHybridization'] = analysis_parameters['RegistrationReferenceHybridization']


    # Group the names by fov
    grpd_by_fov = {}
    for fov in fovs:
        grpd_by_fov[fov] = {}
        grpd_by_fov[fov]['split'] = {}
        grpd_by_fov[fov]['all_names'] = []
    for name, grp in consolidated_grp.items():
        fov = int(name.split('_fov_')[-1])
        grpd_by_fov[fov]['all_names'].append(name)


    # Regroup by processing type
    for fov, data_dict in grpd_by_fov.items():
        grpd_by_fov[fov]['split']['fish'] = {}
        grpd_by_fov[fov]['split']['registration'] = {}
        grpd_by_fov[fov]['split']['staining'] = {}
        
        for name in data_dict['all_names']:
            if consolidated_grp[name].attrs['stitching_channel'] == consolidated_grp[name].attrs['channel']:
                if consolidated_grp[name].attrs['channel'] in grpd_by_fov[fov]['split']['registration'].keys():
                    grpd_by_fov[fov]['split']['registration'][consolidated_grp[name].attrs['channel']][0].append(name)
                else:
                    grpd_by_fov[fov]['split']['registration'][consolidated_grp[name].attrs['channel']] = ([],registration_selected_parameters)
                    grpd_by_fov[fov]['split']['registration'][consolidated_grp[name].attrs['channel']][0].append(name)
            elif '_ST' in consolidated_grp[name].attrs['target_name']:
                if consolidated_grp[name].attrs['channel'] in grpd_by_fov[fov]['split']['staining'].keys():
                    grpd_by_fov[fov]['split']['staining'][consolidated_grp[name].attrs['channel']][0].append(name)
                else:
                    grpd_by_fov[fov]['split']['staining'][consolidated_grp[name].attrs['channel']] = ([],staining_selected_parameters)
                    grpd_by_fov[fov]['split']['staining'][consolidated_grp[name].attrs['channel']][0].append(name)
            else:
                if consolidated_grp[name].attrs['channel'] in grpd_by_fov[fov]['split']['fish'].keys():
                    grpd_by_fov[fov]['split']['fish'][consolidated_grp[name].attrs['channel']][0].append(name)
                else:
                    grpd_by_fov[fov]['split']['fish'][consolidated_grp[name].attrs['channel']] = ([],fish_selected_parameters)
                    grpd_by_fov[fov]['split']['fish'][consolidated_grp[name].attrs['channel']][0].append(name)

    return grpd_by_fov

def create_dir(dir_path):
    try:
        os.stat(dir_path)
    except:
        os.mkdir(dir_path)
        os.chmod(dir_path,0o777)



def register_images(img, shift):
        offset_image = fourier_shift(np.fft.fftn(img), shift)
        offset_image = np.fft.ifftn(offset_image).real
        return offset_image


def combine_rounds_images(images_path_list,experiment_fpath, experiment_info,all_rounds_shifts, save=True):
    experiment_fpath = Path(experiment_fpath)
    
    rounds_num = experiment_info['Barcode_length']
    fov = images_path_list[0].stem.split('_')[-2]
    channel = images_path_list[0].stem.split('_')[-4]
    experiment_name = experiment_info['EXP_name']

    fpath_dots = images_path_list[0]
    fname = (fpath_dots.stem).split('_dots')[0] + '_filtered.pkl'
    path_img = fpath_dots.parent.parent / 'filtered_images' / fname
    img, meta = pickle.load(open(path_img,'rb'))
    img_stack = np.zeros([rounds_num,img.shape[0],img.shape[1]])

    for fpath_dots in images_path_list:
        fname = (fpath_dots.stem).split('_dots')[0] + '_filtered.pkl'
        path_img = fpath_dots.parent.parent / 'filtered_images'
        fpath = path_img / fname
        round_num = int(fname.split('_')[-5].split('Hybridization')[-1])
        img, meta = pickle.load(open(fpath,'rb'))
        img[img<0] = 0
        shift = all_rounds_shifts[round_num]
        shifted_img = register_images(img,shift)
        img_stack[round_num-1,:,:] = shifted_img

    fpath = experiment_fpath / 'tmp' / 'combined_rounds_images' / (experiment_name + '_' + channel + '_combined_img_fov_' + fov + '.npy')
    
    if save:
        # Add conversion to more compress ftype
        np.save(fpath, img_stack)
    return img_stack



def combine_filtered_images(output_list,experiment_fpath,metadata, save=True):
    experiment_fpath = Path(experiment_fpath)
     
    img_stack = np.zeros([metadata['total_rounds'],metadata['img_width'],metadata['img_height']])

    for img, img_meta in output_list:
        round_num = img_meta.round_num
        img_stack[round_num-1,:,:] = shifted_img

    if save:
        # Add conversion to more compress ftype
        channel = img_meta.channel
        fov = img_meta.fov_num
        fpath = experiment_fpath / 'tmp' / 'combined_rounds_images' / (experiment_name + '_' + channel + '_combined_img_fov_' + fov + '.npy')
        np.save(fpath, img_stack)
    
    return img_stack






def register_combined_rounds_images(combined_round_images,all_rounds_shifts):
    
    img_stack = np.zeros_like(combined_round_images)

    for round_num, shift in all_rounds_shifts.items():
        img = combined_round_images[round_num-1,:,:]
        img[img<0] = 0
        shifted_img = register_images(img,shift)
        img_stack[round_num-1,:,:] = shifted_img

    return img_stack

# https://stackoverflow.com/questions/1868714/how-do-i-copy-an-entire-directory-of-files-into-an-existing-directory-using-pyth
def copytree(src, dst, symlinks=False, ignore=None):
    if not os.path.exists(dst):
        os.makedirs(dst)
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            copytree(s, d, symlinks, ignore)
        else:
            if not os.path.exists(d) or os.stat(s).st_mtime - os.stat(d).st_mtime > 1:
                shutil.copy2(s, d)


def not_run_counting_sorted_grps(experiment_fpath, sorted_grps):
    """
    Helper function to identify the files not processed for counting
    because a failure of dask/htcondor. Identify the files that
    remain to be processed.
    """
    experiment_fpath = Path(experiment_fpath)
    counts_fpath = experiment_fpath / 'tmp' / 'raw_counts'
    
    non_processed_sorted_grps = {}
    
    files_processed = list(counts_fpath.glob('*_dots.pkl'))
    files_processed = [el.stem for el in files_processed]
    for key, (grp, param) in sorted_grps.items():
        non_processed_grp = []
        for zarr_grp_name in grp:
            if zarr_grp_name in files_processed:
                files_processed.remove(zarr_grp_name)
            else:
                non_processed_grp.append(zarr_grp_name)
        non_processed_sorted_grps[key] = (non_processed_grp, param)
    return non_processed_sorted_grps






# # to avoid reference for nested structures
# # https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml (comment)
# yaml.SafeDumper.ignore_aliases = lambda *args : True


# # to avoid reference for nested structures
# # https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml (comment)
# yaml.SafeDumper.ignore_aliases = lambda *args : True

def convert_to_uint16_full_float64_range(image):
    # Clip the values above 1
    image = image / np.finfo(np.float64).max
    # image[image > 1] = 1
    # in the current processing step the images are float
    # but not anymore between (0,1)
    # Scale to the max of the uint16
    image *= np.iinfo(np.uint16).max
    # Round and convert to integer
    image = np.uint16(np.rint(image))
    return image

def convert_to_uint16(image):
    """
    Consider that the original images were uint16
    that is why i normalize for np.iinfo(np.uint16).max*1.0
    otherwise the values will be to small
    """
 
    image = image / (np.iinfo(np.uint16).max*1.0)
    image *= np.iinfo(np.uint16).max
    # Round and convert to integer
    image = image.astype(np.uint16)
    return image

def convert_from_uint16_to_float64(image):
    image = image / np.iinfo(np.uint16).max
    image = img_as_float64(image)
    return image


def load_pipeline_config_file(pipeline_config_fpath):

    try:
        pipeline_config_dict = OrderedDict(yaml.safe_load(open(pipeline_config_fpath, 'rb')))
        return pipeline_config_dict
    except FileExistsError:
        logging.exception(f'{pipeline_config_fpath} missing')
        sys.exit(f'{pipeline_config_fpath} missing')
    except NameError:
        logging.exception(f'{pipeline_config_fpath} wrong pathway name')
        sys.exit(f'{pipeline_config_fpath} wrong pathway name')



def load_running_analysis_config_file(experiment_fpath):
    experiment_fpath = Path(experiment_fpath)
    config_fpath = experiment_fpath /  'pipeline_config'
    running_analysis_config_fpath = list(config_fpath.glob('*_analysis_run.yaml'))
    if len(running_analysis_config_fpath) == 0:
        logging.exception(f'missing running analysis config file')
        sys.exit(f'missing running analysis config file')
    elif len(running_analysis_config_fpath) > 1:
        logging.exception(f'too many running analysis config file')
        sys.exit(f'too many running analysis config file')
    else:
        running_analysis_config_fpath = running_analysis_config_fpath[0]
        running_analysis_config = load_pipeline_config_file(running_analysis_config_fpath)
        return running_analysis_config



def modify_images_config(experiment_fpath:str,keyword:str):
    """
    Function used to batch modify the values of the
    processing parameters in the images_config.yaml
    files

    Parameters:
    -----------
    experiment_fpath: str
        path of the experiment to process
    keyword: str
        selection criteria used to process a subset of
        configuration files

    """
    experiment_fpath = Path(experiment_fpath)
    pipeline_config_dpath = experiment_fpath / 'pipeline_config'

    for config_fpath in pipeline_config_dpath.glob('*_images_config.yaml'):
        if keyword in config_fpath.name:
            img_config_dict = load_pipeline_config_file(config_fpath)
            for fov_name, fov_data in img_config_dict.items():
                for round_name, round_data in fov_data.items():
                    # round_data['analysis_parameters'] = {}
                    # round_data['analysis_parameters']['preprocessing'] = {}
                    # round_data['analysis_parameters']['preprocessing']['flat_field_kernel'] = (2,100,100)
                    # round_data['analysis_parameters']['preprocessing']['filtering_small_kernel'] = (1,8,8)
                    # round_data['analysis_parameters']['preprocessing']['filtering_laplacian_kernel'] = (0.2,0.5,0.5)
                    # round_data['analysis_parameters']['preprocessing']['large_obj_removal_percentile'] = 99
                    # round_data['analysis_parameters']['preprocessing']['large_obj_removal_min_obj_size']= 50
                    # round_data['analysis_parameters']['preprocessing']['large_obj_removal_selem'] = 3
                    round_data['analysis_parameters']['counting'] = {}
                    round_data['analysis_parameters']['counting']['min_distance'] = 4
                    round_data['analysis_parameters']['counting']['min_obj_size'] = 4
                    round_data['analysis_parameters']['counting']['max_obj_size'] = 200
                    round_data['analysis_parameters']['beads_registration'] = {}

            with open(config_fpath, 'w') as new_config:
                    yaml.safe_dump(dict(img_config_dict), new_config,default_flow_style=False,sort_keys=False)



def create_stringency_selection_fov_yaml(pipeline_config_fpath):
    
    pipeline_config_dict = load_pipeline_config_file(pipeline_config_fpath)
    fovs = (pipeline_config_dict['channel_coords_dict'].keys())
    stringency = pipeline_config_dict['pipeline_required_parameters']['counting_settings']['stringency']
    pipeline_config_dict['pipeline_required_parameters']['counting_settings']['stringency'] = {}
    pipeline_config_dict['pipeline_required_parameters']['counting_settings']['stringency'] = {k:int(stringency) for k in fovs}
    with open(pipeline_config_fpath, 'w') as new_config:
            yaml.safe_dump(dict(pipeline_config_dict), new_config,default_flow_style=False,sort_keys=False)

def create_specific_fov_counting_thr_yaml(pipeline_config_fpath):
    pipeline_config_dict = load_pipeline_config_file(pipeline_config_fpath)
    fovs = (pipeline_config_dict['channel_coords_dict'].keys())
    stringency = pipeline_config_dict['pipeline_required_parameters']['counting_settings']['thr_fov'] = {}
    pipeline_config_dict['pipeline_required_parameters']['counting_settings']['thr_fov'] = {k:None for k in fovs}
    with open(pipeline_config_fpath, 'w') as new_config:
            yaml.safe_dump(dict(pipeline_config_dict), new_config,default_flow_style=False,sort_keys=False)

def select_dir(processing_directory,keyword):
    list_hyb_fpaths = []
    processed_dirs = next(os.walk(processing_directory))[1]
    
    if keyword:
        for dir_name in processed_dirs:
            if keyword in dir_name:
                dir_path = processing_directory / dir_name
                list_hyb_fpaths.append(dir_path)
    else:
        for dir_name in processed_dirs:
            dir_path = processing_directory / dir_name
            list_hyb_fpaths.append(dir_path)

    return list_hyb_fpaths



class OptionEatAll(click.Option):
    '''
    https://stackoverflow.com/questions/48391777/nargs-equivalent-for-options-in-click

    Option class to ingest undefined number of options.
    '''

    def __init__(self, *args, **kwargs):
        self.save_other_options = kwargs.pop('save_other_options', True)
        nargs = kwargs.pop('nargs', -1)
        assert nargs == -1, 'nargs, if set, must be -1 not {}'.format(nargs)
        super(OptionEatAll, self).__init__(*args, **kwargs)
        self._previous_parser_process = None
        self._eat_all_parser = None

    def add_to_parser(self, parser, ctx):

        def parser_process(value, state):
            # method to hook to the parser.process
            done = False
            value = [value]
            if self.save_other_options:
                # grab everything up to the next option
                while state.rargs and not done:
                    for prefix in self._eat_all_parser.prefixes:
                        if state.rargs[0].startswith(prefix):
                            done = True
                    if not done:
                        value.append(state.rargs.pop(0))
            else:
                # grab everything remaining
                value += state.rargs
                state.rargs[:] = []
            value = tuple(value)

            # call the actual process
            self._previous_parser_process(value, state)

        retval = super(OptionEatAll, self).add_to_parser(parser, ctx)
        for name in self.opts:
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break
        return retval