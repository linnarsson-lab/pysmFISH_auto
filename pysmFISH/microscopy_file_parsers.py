
from typing import *
import pickle
import yaml
import sys
import re
import shutil
import numpy as np
import zarr
import nd2reader
import xarray as xr
from pathlib import Path


import prefect
from prefect import Task
from prefect import task
from prefect.engine import signals

from pysmFISH.logger_utils import selected_logger


def nd2_raw_files_selector(experiment_fpath: str) -> list:
    """
    Identify the nd2 raw microscopy files generated by
    the robofish machine. The files must contain CountXXXXX in the name. 

    Args:
        experiment_fpath: str 
            Path to the folder to process. It need to contain the '_auto'
            suffix in order to be process with the automated pipeline

    Returns:
        all_files_to_process: list
            List of PosixPath of the microscopy files to process
        
    """
    logger = selected_logger()

    experiment_fpath = Path(experiment_fpath)
    assert '_auto' in experiment_fpath.stem, sys.exit('no _auto in the experiment name')

    searching_key = '*Count*.nd2'
    all_files_to_process = list(experiment_fpath.glob(searching_key))

    assert all_files_to_process, sys.exit('no .nd2 raw files to process')
    
    logger.debug(f'Number of files to process {len(all_files_to_process)}.')
    return all_files_to_process




def nd2_raw_files_selector_general(folder_fpath: str) -> list:

    """
    Function used to identify the .nd2 files in a folder. The files
    do not need to have the CountXXX or the _auto suffix. This
    class can be used to identify the .nd2 files that need to be
    reparsed.

    Args:
        folder_fpath: str 
            Path to the folder to process. 
    Returns:
        all_files_to_process: list
            List of PosixPath of the microscopy files to process
        
    """
    logger = selected_logger()
    folder_fpath = Path(folder_fpath)
    
    searching_key = '*.nd2'
    all_files_to_process = list(folder_fpath.glob(searching_key))

    assert all_files_to_process, sys.exit('no .nd2 raw files to process')
    
    logger.debug(f'Number of files to process {len(all_files_to_process)}.')
    return all_files_to_process



def nikon_nd2_autoparser_zarr(nd2_file_path, parsed_raw_data_fpath, experiment_info):

    """
    Function to parse nikon .nd2 files generated by robofish.
    This parser not consider the possibility to have multiple experiment running at
    the same time with the raw imaging data present in the same folder. 
    Once the data are parsed by hybridization are saved in
    the same folder.

    NB: The file and the corresponding metadata file generated by the microscope
        must contain 'CountXXXXX' in order to be processed

    The parsed_tmp and raw_data directories are created in a previous step of the 
    pipeline during data sorting

    QC step run before parsing take care of checking if the nd2 files are present
    and in the corresponding pickle configuration files are in the folder

    This nd2 parser process one .nd2 file

    Args:
        nd2_file_path: str
            Path to the .nd2 file to be parsed
        parsed_raw_data_fpath: str
            Path to the zarr file that will store the parsed data
        experiment_info: dict
            Dictionary with overall experiment info

    Returns:
        processing_info: list of tuples
        each tuple contain the path of the zarr storage and the number of fov
    """

    logger = selected_logger()

    nd2_file_path = Path(nd2_file_path)
    nd2_fname = nd2_file_path.stem

    logger.debug(f'processing file {nd2_fname}')

    experiment_fpath = nd2_file_path.parent
    experiment_name = nd2_file_path.parent.stem
    experiment_name = experiment_name.split('_auto')[0]
    
    raw_files_dir = experiment_fpath / 'raw_data'
    
    # Extract the Count code from the file name
    count_code = re.search(r'(Count)\d{5}', nd2_file_path.stem).group()
    
    all_info_files = experiment_fpath.glob('*.pkl')
    
    info_file = [info_file for info_file in all_info_files if count_code  in info_file.stem][0]
    
    info_data = pickle.load(open(info_file, 'rb'))
    logger.debug(f'loaded info data file {info_file.stem}')
    
    hybridization_name = info_data['channels']['Hybridization']
    logger.debug(f'processing hybridization {hybridization_name}')
    
    try:
        nd2fh = nd2reader.ND2Reader(nd2_file_path)
    except:
        logger.error('Cannot load the nd2 file')
        err = signals.FAIL('Cannot load the nd2 file')
        raise err
    else:
        # Collect metadata
        all_metadata = nd2fh.parser._raw_metadata
        parsed_metadata = nd2fh.parser._raw_metadata.get_parsed_metadata()
        
        channel = parsed_metadata['channels'][0] # works because there is only one channel for file
        logger.debug(f'processing channel {channel}')

        img_shape = np.array([parsed_metadata['height'],parsed_metadata['width']])
        pixel_microns = parsed_metadata['pixel_microns']
        z_levels = parsed_metadata['z_levels']
        fields_of_view = parsed_metadata['fields_of_view']
        
        # Collect FOV coords
        x_data = np.array(all_metadata.x_data)
        x_data = x_data[:,np.newaxis]
        y_data = np.array(all_metadata.y_data)
        y_data = y_data[:,np.newaxis]
        z_data = np.array(all_metadata.z_data)
        z_data = z_data[:,np.newaxis]
        all_coords = np.hstack((z_data,x_data,y_data))
        fov_coords = all_coords[0::len(z_levels),:]
        
        tag_name = experiment_name + '_' + hybridization_name + '_' + channel
        
        
        # Save the file as zarr
        store = zarr.DirectoryStore(parsed_raw_data_fpath)
        root = zarr.group(store=store,overwrite=False)


        # Save the fov_coords
        fname = experiment_fpath / 'tmp' / 'microscope_tiles_coords' / (tag_name + '_fovs_coords.npy')
        np.save(fname, fov_coords)

        nd2fh.bundle_axes = 'zyx'
        # set iteration over the fields of view
        nd2fh.iter_axes = 'v'
        
        # Save coords of the FOV
        rows = np.arange(parsed_metadata['width'])
        cols = np.arange(parsed_metadata['height'])
        hybridization_num = int(hybridization_name.split('Hybridization')[-1])
        for fov in fields_of_view:
            img = np.array(nd2fh[fov],dtype=np.uint16)      
            array_name = tag_name + '_fov_' + str(fov)
            dgrp = root.create_group(array_name)
            fov_name = 'raw_data_fov_' + str(fov)
            # Remember that attrs must be JSON-serializable to be stored
            # in zarr
            dgrp.attrs['grp_name'] = array_name
            dgrp.attrs['fov_name'] = fov_name
            dgrp.attrs['channel'] = channel
            dgrp.attrs['target_name'] = info_data['channels'][channel]
            dgrp.attrs['img_width'] = parsed_metadata['width']
            dgrp.attrs['img_height'] = parsed_metadata['height']
            dgrp.attrs['pixel_microns'] = pixel_microns
            dgrp.attrs['z_levels'] = list(z_levels)
            dgrp.attrs['fields_of_view'] = list(fields_of_view)
            dgrp.attrs['fov_num'] = fov
            dgrp.attrs['stitching_channel'] = info_data['StitchingChannel']
            dgrp.attrs['stitching_type'] = experiment_info['Stitching_type']
            dgrp.attrs['experiment_type'] = experiment_info['Experiment_type']
            dgrp.attrs['hybridization_num'] = hybridization_num
            dgrp.attrs['experiment_name'] = experiment_name
            dgrp.attrs['fov_acquisition_coords_x'] = fov_coords[fov,1]
            dgrp.attrs['fov_acquisition_coords_y'] = fov_coords[fov,2]
            dgrp.attrs['fov_acquisition_coords_z'] = fov_coords[fov,0]


            if info_data['StitchingChannel'] == channel:
                dgrp.attrs['processing_type'] = dgrp.attrs['stitching_type']
            elif '_ST' in dgrp.attrs['target_name']:
                dgrp.attrs['processing_type'] = 'staining'
            else:
                dgrp.attrs['processing_type'] = 'fish'

            dset = dgrp.create_dataset(fov_name, data=img, shape=img.shape, chunks=(1,None,None),overwrite=True)

        # Rename the nd2 files
        # new_file_name = tag_name + '.nd2'
        # logger.debug(f'.nd2 file renamed {new_file_name}')
        # new_file_path = raw_files_dir / new_file_name
        # nd2_file_path.rename(new_file_path)
        # nd2_file_path = new_file_path
        
        # # Copy the pkl files
        # new_file_name = tag_name + '_info.pkl'
        # logger.debug(f'info data file renamed {new_file_name}')
        # new_file_path = raw_files_dir / new_file_name
        # # Must copy the pkl file in order to be able to use the file for the other channels
        # shutil.copy(str(info_file), str(new_file_path))
        
        # Save the fov_coords
        # fname = experiment_fpath / 'tmp' / (tag_name + '_fovs_coords.npy')
        # np.save(fname, fov_coords)




def nikon_nd2_reparser_zarr(nd2_file_path,parsed_raw_data_fpath,experiment_info):
    """
    This function is used to reparse the raw data stored in the raw_data folder during
    the processing. 

    NB: The file and the corresponding metadata file generated by the microscope
        have the same starting name:
        LBEXP20201014_EEL_Mouse_2420um_Hybridization01_Cy5.nd2
        LBEXP20201014_EEL_Mouse_2420um_Hybridization01_Cy5_info.pkl

    The parsed_tmp directory is created in a previous step of the 
    pipeline during data sorting


    This nd2 parser process one .nd2 file

    Args:
        nd2_file_path: str
            Path to the .nd2 file to be parsed
        parsed_raw_data_fpath: str
            Path to the zarr file that will store the parsed data
        experiment_info: dict
            Dictionary with overall experiment info

    Returns:
        processing_info: list of tuples
        each tuple contain the path of the zarr storage and the number of fov
    """

    logger = selected_logger()
    nd2_file_path = Path(nd2_file_path)
    nd2_fname = nd2_file_path.stem

    logger.debug(f'processing file {nd2_fname}')

    experiment_fpath = nd2_file_path.parent.parent
    experiment_name = nd2_file_path.parent.parent.stem
    experiment_name = experiment_name.split('_auto')[0]
    
    info_file = nd2_file_path.parent / (nd2_fname + '_info.pkl')
    info_data = pickle.load(open(info_file, 'rb'))
    hybridization_name = info_data['channels']['Hybridization']
    
    try:
        nd2fh = nd2reader.ND2Reader(nd2_file_path)
    except:
        logger.error('Cannot load the {nd2_fname} nd2 file')
        sys.exit('Cannot load the {nd2_fname} nd2 file')
    else:
        # Collect metadata
        all_metadata = nd2fh.parser._raw_metadata
        parsed_metadata = nd2fh.parser._raw_metadata.get_parsed_metadata()
        
        channel = parsed_metadata['channels'][0] # works because there is only one channel for file
        img_shape = np.array([parsed_metadata['height'],parsed_metadata['width']])
        pixel_microns = parsed_metadata['pixel_microns']
        z_levels = parsed_metadata['z_levels']
        fields_of_view = parsed_metadata['fields_of_view']
        
        # Collect FOV coords
        x_data = np.array(all_metadata.x_data)
        x_data = x_data[:,np.newaxis]
        y_data = np.array(all_metadata.y_data)
        y_data = y_data[:,np.newaxis]
        z_data = np.array(all_metadata.z_data)
        z_data = z_data[:,np.newaxis]
        all_coords = np.hstack((z_data,x_data,y_data))
        fov_coords = all_coords[0::len(z_levels),:]
        
        tag_name = experiment_name + '_' + hybridization_name + '_' + channel
    
        # Save the fov_coords
        fname = experiment_fpath / 'tmp' / 'microscope_tiles_coords' / (tag_name + '_fovs_coords.npy')
        np.save(fname, fov_coords)


        # Save the file as zarr
        store = zarr.DirectoryStore(parsed_raw_data_fpath)
        root = zarr.group(store=store,overwrite=False)

        
        nd2fh.bundle_axes = 'zyx'
        # set iteration over the fields of view
        nd2fh.iter_axes = 'v'
        
        # Save coords of the FOV
        rows = np.arange(parsed_metadata['width'])
        cols = np.arange(parsed_metadata['height'])
        hybridization_num = int(hybridization_name.split('Hybridization')[-1])
        for fov in fields_of_view:
            img = np.array(nd2fh[fov],dtype=np.uint16)      
            array_name = tag_name + '_fov_' + str(fov)
            dgrp = root.create_group(array_name)
            fov_name = 'raw_data_fov_' + str(fov)
            # Remember that attrs must be JSON-serializable to be stored
            # in zarr
            dgrp.attrs['grp_name'] = array_name
            dgrp.attrs['fov_name'] = fov_name
            dgrp.attrs['channel'] = channel
            dgrp.attrs['target_name'] = info_data['channels'][channel]
            dgrp.attrs['img_width'] = parsed_metadata['width']
            dgrp.attrs['img_height'] = parsed_metadata['height']
            dgrp.attrs['pixel_microns'] = pixel_microns
            dgrp.attrs['z_levels'] = list(z_levels)
            dgrp.attrs['fields_of_view'] = list(fields_of_view)
            dgrp.attrs['fov_num'] = fov
            dgrp.attrs['stitching_channel'] = info_data['StitchingChannel']
            dgrp.attrs['stitching_type'] = experiment_info['Stitching_type']
            dgrp.attrs['experiment_type'] = experiment_info['Experiment_type']
            dgrp.attrs['hybridization_num'] = hybridization_num
            dgrp.attrs['experiment_name'] = experiment_name
            dgrp.attrs['fov_acquisition_coords_x'] = fov_coords[fov,1]
            dgrp.attrs['fov_acquisition_coords_y'] = fov_coords[fov,2]
            dgrp.attrs['fov_acquisition_coords_z'] = fov_coords[fov,0]


            if info_data['StitchingChannel'] == channel:
                dgrp.attrs['processing_type'] = dgrp.attrs['stitching_type']
            elif '_ST' in dgrp.attrs['target_name']:
                dgrp.attrs['processing_type'] = 'staining'
            else:
                dgrp.attrs['processing_type'] = 'fish'

            dset = dgrp.create_dataset(fov_name, data=img, shape=img.shape, chunks=(1,None,None),overwrite=True)
