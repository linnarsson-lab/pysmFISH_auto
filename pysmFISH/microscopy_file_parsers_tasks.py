
from typing import *
import pickle
import yaml
import re
import shutil
import numpy as np
import zarr
import nd2reader
import xarray as xr
from pathlib import Path


import prefect
from prefect import task
from prefect.engine import signals

from pysmFISH.logger_utils import prefect_logging_setup


# from pysmFISH.utils import load_pipeline_config_file, create_dir, load_running_analysis_config_file


""" 
The parsing of the nikon files require the bftools. We use only the inf
command and the location is inferred using os.system('inf')
bftool 
"""


@task(name='nd2_files_selection')
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

    logger = prefect.utilities.logging.get_logger("parsing")
    
    assert '_auto' in experiment_fpath.stem, signals.FAIL('no _auto in the experiment name')

    experiment_fpath = Path(experiment_fpath)
    searching_key = '*Count*.nd2'
    all_files_to_process = list(experiment_fpath.glob(searching_key))

    assert all_files_to_process, signals.FAIL('no .nd2 raw files to process')
    
    logger.debug(f'Number of files to process {len(all_files_to_process)}.')
    return all_files_to_process



@task(name='nd2_autoparser')
def nikon_nd2_autoparser(nd2_file_path,parsed_raw_data_fpath):
    """
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

    Returns:
        processing_info: list of tuples
        each tuple contain the path of the zarr storage and the number of fov
    """

    logger = prefect_logging_setup("auto_parser")

    nd2_file_path = Path(nd2_file_path)
    nd2_fname = nd2_file_path.stem

    logger.debug(f'processing file {nd2_fname}')

    experiment_fpath = nd2_file_path.parent
    experiment_name = nd2_file_path.parent.stem
    experiment_name = experiment_name.split('_auto')[0]
       
    raw_files_dir = experiment_fpath / 'raw_data'
    parsed_tmp = experiment_fpath / 'parsed_tmp'
    

    # Extract the Count code from the file name
    count_code = re.search(r'(Count)\d{5}', nd2_file_path.stem).group()
    
    all_info_files = experiment_fpath.glob('*.pkl')
    
    info_file = [info_file for info_file in all_info_files if count_code  in info_file.stem][0]
    
    info_data = pickle.load(open(info_file, 'rb'))
    hybridization_name = info_data['channels']['Hybridization']
    
    try:
        nd2fh = nd2reader.ND2Reader(nd2_file_path)
    except:
        logger.error('Cannot load the nd2 file')
        signals.FAIL('Cannot load the nd2 file')
    else:
        # Collect metadata
        all_metadata = nd2fh.parser._raw_metadata
        parsed_metadata = nd2fh.parser._raw_metadata.get_parsed_metadata()
        
        channel = parsed_metadata['channels'][0] # works because there is only one channel for file
        img_height = parsed_metadata['height']
        img_width = parsed_metadata['width']
        pixel_microns = parsed_metadata['pixel_microns']
        z_levels = parsed_metadata['z_levels']
        fields_of_view = parsed_metadata['fields_of_view']
        
        # Collect FOV coords
        x_data = np.array(all_metadata.x_data)
        x_data = x_data[:,np.newaxis]
        y_data = np.array(all_metadata.x_data)
        y_data = y_data[:,np.newaxis]
        z_data = np.array(all_metadata.z_data)
        z_data = z_data[:,np.newaxis]
        all_coords = np.hstack((z_data,x_data,y_data))
        fov_coords = all_coords[0::len(z_levels),:]
        
        tag_name = experiment_name + '_' + hybridization_name + '_' + channel
        
        
        # Save the file as zarr
        # zarr_store = parsed_tmp / (experiment_name + '_' + hybridization_name + '_' + channel + '_raw_images_tmp.zarr')
        nd2fh.bundle_axes = 'zyx'
        # set iteration over the fields of view
        nd2fh.iter_axes = 'v'
        
        # Save coords of the FOV
        rows = np.arange(img_width)
        cols = np.arange(img_height)
        hybridization_num = int(hybridization_name.split('Hybridization')[-1])
        for fov in fields_of_view:
            img = np.array(nd2fh[fov],dtype=np.uint16)
            fov_attrs = {'channel': channel,
                'target_name': info_data['channels'][channel],
                'img_height': img_height,
                'img_width': img_width,
                'pixel_microns': pixel_microns,
                'z_levels':list(z_levels),
                'fov_num': fov,
                'StitchingChannel': info_data['StitchingChannel'],
                'hybridization_num': hybridization_num,
                'experiment_name' : experiment_name} 
            datarray_name = tag_name + '_fov_' + str(fov) 
            img_xarray = xr.DataArray(img, coords={'z_levels':z_levels,'rows':rows, 'cols':cols}, 
                                    dims=['z_levels','rows','cols'],attrs=fov_attrs, name=datarray_name)
            img_xarray = img_xarray.chunk(chunks=(1,img_width,img_height))
            ds = img_xarray.to_dataset(name = datarray_name)
            ds.to_zarr(parsed_raw_data_fpath, mode='a', consolidated=True)
                
        # Rename the nd2 files
        # new_file_name = tag_name + '.nd2'
        # new_file_path = raw_files_dir / new_file_name
        # nd2_file_path.rename(new_file_path)
        # nd2_file_path = new_file_path
        
        # # Copy the pkl files
        # new_file_name = tag_name + '_info.pkl'
        # new_file_path = raw_files_dir / new_file_name
        # # Must copy the pkl file in order to be able to use the file for the other channels
        # shutil.copy(str(info_file), str(new_file_path))
        
        # # Save the fov_coords
        # fname = experiment_fpath / 'tmp' / (tag_name + '_fovs_coords.npy')
        # np.save(fname, fov_coords)

        # fovs = list(fields_of_view)
        # list_store = [zarr_store] * len(fovs)
        # processing_info = list(zip(list_store,fovs))

        # return processing_info


@task(name='nd2_autoparser_single_files')
def nikon_nd2_autoparser_single_files(nd2_file_path,parsed_raw_data_fpath):
    """
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

    Returns:
        processing_info: list of tuples
        each tuple contain the path of the zarr storage and the number of fov
    """

    logger = prefect_logging_setup("auto_parser_single_files")

    nd2_file_path = Path(nd2_file_path)
    nd2_fname = nd2_file_path.stem

    logger.debug(f'processing file {nd2_fname}')

    experiment_fpath = nd2_file_path.parent
    experiment_name = nd2_file_path.parent.stem
    experiment_name = experiment_name.split('_auto')[0]
       
    raw_files_dir = experiment_fpath / 'raw_data'
    parsed_tmp = experiment_fpath / 'parsed_tmp'
    

    # Extract the Count code from the file name
    count_code = re.search(r'(Count)\d{5}', nd2_file_path.stem).group()
    
    all_info_files = experiment_fpath.glob('*.pkl')
    
    info_file = [info_file for info_file in all_info_files if count_code  in info_file.stem][0]
    
    info_data = pickle.load(open(info_file, 'rb'))
    hybridization_name = info_data['channels']['Hybridization']
    
    try:
        nd2fh = nd2reader.ND2Reader(nd2_file_path)
    except:
        logger.error('Cannot load the nd2 file')
        signals.FAIL('Cannot load the nd2 file')
    else:
        # Collect metadata
        all_metadata = nd2fh.parser._raw_metadata
        parsed_metadata = nd2fh.parser._raw_metadata.get_parsed_metadata()
        
        channel = parsed_metadata['channels'][0] # works because there is only one channel for file
        img_height = parsed_metadata['height']
        img_width = parsed_metadata['width']
        pixel_microns = parsed_metadata['pixel_microns']
        z_levels = parsed_metadata['z_levels']
        fields_of_view = parsed_metadata['fields_of_view']
        
        # Collect FOV coords
        x_data = np.array(all_metadata.x_data)
        x_data = x_data[:,np.newaxis]
        y_data = np.array(all_metadata.x_data)
        y_data = y_data[:,np.newaxis]
        z_data = np.array(all_metadata.z_data)
        z_data = z_data[:,np.newaxis]
        all_coords = np.hstack((z_data,x_data,y_data))
        fov_coords = all_coords[0::len(z_levels),:]
        
        tag_name = experiment_name + '_' + hybridization_name + '_' + channel
        
        
        # Save the file as zarr
        zarr_store = parsed_tmp / (experiment_name + '_' + hybridization_name + '_' + channel + '_raw_images_tmp.zarr')
        
        # Create empty zarr file
        dataset = xr.Dataset()
        dataset.to_zarr(zarr_store, mode='w', consolidated=True)
        
        nd2fh.bundle_axes = 'zyx'
        # set iteration over the fields of view
        nd2fh.iter_axes = 'v'
        
        # Save coords of the FOV
        rows = np.arange(img_width)
        cols = np.arange(img_height)
        hybridization_num = int(hybridization_name.split('Hybridization')[-1])
        for fov in fields_of_view:
            img = np.array(nd2fh[fov],dtype=np.uint16)
            fov_attrs = {'channel': channel,
                'target_name': info_data['channels'][channel],
                'img_height': img_height,
                'img_width': img_width,
                'pixel_microns': pixel_microns,
                'z_levels':list(z_levels),
                'fov_num': fov,
                'StitchingChannel': info_data['StitchingChannel'],
                'hybridization_num': hybridization_num,
                'experiment_name' : experiment_name} 
            datarray_name = tag_name + '_fov_' + str(fov) 
            img_xarray = xr.DataArray(img, coords={'z_levels':z_levels,'rows':rows, 'cols':cols}, 
                                    dims=['z_levels','rows','cols'],attrs=fov_attrs, name=datarray_name)
            img_xarray = img_xarray.chunk(chunks=(1,img_width,img_height))
            ds = img_xarray.to_dataset(name = datarray_name)
            ds.to_zarr(zarr_store, mode='a', consolidated=True)
                
        # # Rename the nd2 files
        # new_file_name = tag_name + '.nd2'
        # new_file_path = raw_files_dir / new_file_name
        # nd2_file_path.rename(new_file_path)
        # nd2_file_path = new_file_path
        
        # # Copy the pkl files
        # new_file_name = tag_name + '_info.pkl'
        # new_file_path = raw_files_dir / new_file_name
        # # Must copy the pkl file in order to be able to use the file for the other channels
        # shutil.copy(str(info_file), str(new_file_path))
        
        # Save the fov_coords
        # fname = experiment_fpath / 'tmp' / (tag_name + '_fovs_coords.npy')
        # np.save(fname, fov_coords)

        # fovs = list(fields_of_view)
        # list_store = [zarr_store] * len(fovs)
        # processing_info = list(zip(list_store,fovs))

        # return processing_info

@task(name='nd2_autoparser')
def nikon_nd2_autoparser_zarr_single_files(nd2_file_path,parsed_raw_data_fpath):
    """
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

    Returns:
        processing_info: list of tuples
        each tuple contain the path of the zarr storage and the number of fov
    """

    logger = prefect_logging_setup("auto_parser")

    nd2_file_path = Path(nd2_file_path)
    nd2_fname = nd2_file_path.stem

    logger.debug(f'processing file {nd2_fname}')

    experiment_fpath = nd2_file_path.parent
    experiment_name = nd2_file_path.parent.stem
    experiment_name = experiment_name.split('_auto')[0]
       
    raw_files_dir = experiment_fpath / 'raw_data'
    parsed_tmp = experiment_fpath / 'parsed_tmp'
    

    # Extract the Count code from the file name
    count_code = re.search(r'(Count)\d{5}', nd2_file_path.stem).group()
    
    all_info_files = experiment_fpath.glob('*.pkl')
    
    info_file = [info_file for info_file in all_info_files if count_code  in info_file.stem][0]
    
    info_data = pickle.load(open(info_file, 'rb'))
    hybridization_name = info_data['channels']['Hybridization']
    
    try:
        nd2fh = nd2reader.ND2Reader(nd2_file_path)
    except:
        logger.error('Cannot load the nd2 file')
        signals.FAIL('Cannot load the nd2 file')
    else:
        # Collect metadata
        all_metadata = nd2fh.parser._raw_metadata
        parsed_metadata = nd2fh.parser._raw_metadata.get_parsed_metadata()
        
        channel = parsed_metadata['channels'][0] # works because there is only one channel for file
        img_height = parsed_metadata['height']
        img_width = parsed_metadata['width']
        pixel_microns = parsed_metadata['pixel_microns']
        z_levels = parsed_metadata['z_levels']
        fields_of_view = parsed_metadata['fields_of_view']
        
        # Collect FOV coords
        x_data = np.array(all_metadata.x_data)
        x_data = x_data[:,np.newaxis]
        y_data = np.array(all_metadata.x_data)
        y_data = y_data[:,np.newaxis]
        z_data = np.array(all_metadata.z_data)
        z_data = z_data[:,np.newaxis]
        all_coords = np.hstack((z_data,x_data,y_data))
        fov_coords = all_coords[0::len(z_levels),:]
        
        tag_name = experiment_name + '_' + hybridization_name + '_' + channel
        
        
        # Save the file as zarr
        zarr_store = parsed_tmp / (experiment_name + '_' + hybridization_name + '_' + channel + '_raw_images_tmp.zarr')
        store = zarr.DirectoryStore(zarr_store)
        root = zarr.group(store=store,overwrite=True)

        
        nd2fh.bundle_axes = 'zyx'
        # set iteration over the fields of view
        nd2fh.iter_axes = 'v'
        
        # Save coords of the FOV
        rows = np.arange(img_width)
        cols = np.arange(img_height)
        hybridization_num = int(hybridization_name.split('Hybridization')[-1])
        for fov in fields_of_view:
            img = np.array(nd2fh[fov],dtype=np.uint16)        
            dset = root.create_dataset(fov, data=img, shape=img.shape, chunks=(1,None,None),overwrite=True)
            dset.attrs['channel'] = channel
            dset.attrs['target_name'] = info_data['channels'][channel]
            dset.attrs['img_height'] = img_height
            dset.attrs['img_width'] = img_width
            dset.attrs['pixel_microns'] = pixel_microns
            dset.attrs['z_levels'] = list(z_levels)
            dset.attrs['fov_num'] = fov
            dset.attrs['StitchingChannel'] = info_data['StitchingChannel']
            dset.attrs['hybridization_num'] = hybridization_num
            dset.attrs['experiment_name'] = experiment_name

        zarr.convenience.consolidate_metadata(store)                
        # Rename the nd2 files
        # new_file_name = tag_name + '.nd2'
        # new_file_path = raw_files_dir / new_file_name
        # nd2_file_path.rename(new_file_path)
        # nd2_file_path = new_file_path
        
        # # Copy the pkl files
        # new_file_name = tag_name + '_info.pkl'
        # new_file_path = raw_files_dir / new_file_name
        # # Must copy the pkl file in order to be able to use the file for the other channels
        # shutil.copy(str(info_file), str(new_file_path))
        
        # # Save the fov_coords
        # fname = experiment_fpath / 'tmp' / (tag_name + '_fovs_coords.npy')
        # np.save(fname, fov_coords)

@task(name='nd2_autoparser')
def nikon_nd2_autoparser_zarr(nd2_file_path,parsed_raw_data_fpath):
    """
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

    Returns:
        processing_info: list of tuples
        each tuple contain the path of the zarr storage and the number of fov
    """

    logger = prefect_logging_setup("auto_parser")

    nd2_file_path = Path(nd2_file_path)
    nd2_fname = nd2_file_path.stem

    logger.debug(f'processing file {nd2_fname}')

    experiment_fpath = nd2_file_path.parent
    experiment_name = nd2_file_path.parent.stem
    experiment_name = experiment_name.split('_auto')[0]
       
    raw_files_dir = experiment_fpath / 'raw_data'
    parsed_tmp = experiment_fpath / 'parsed_tmp'
    

    # Extract the Count code from the file name
    count_code = re.search(r'(Count)\d{5}', nd2_file_path.stem).group()
    
    all_info_files = experiment_fpath.glob('*.pkl')
    
    info_file = [info_file for info_file in all_info_files if count_code  in info_file.stem][0]
    
    info_data = pickle.load(open(info_file, 'rb'))
    hybridization_name = info_data['channels']['Hybridization']
    
    try:
        nd2fh = nd2reader.ND2Reader(nd2_file_path)
    except:
        logger.error('Cannot load the nd2 file')
        signals.FAIL('Cannot load the nd2 file')
    else:
        # Collect metadata
        all_metadata = nd2fh.parser._raw_metadata
        parsed_metadata = nd2fh.parser._raw_metadata.get_parsed_metadata()
        
        channel = parsed_metadata['channels'][0] # works because there is only one channel for file
        img_height = parsed_metadata['height']
        img_width = parsed_metadata['width']
        pixel_microns = parsed_metadata['pixel_microns']
        z_levels = parsed_metadata['z_levels']
        fields_of_view = parsed_metadata['fields_of_view']
        
        # Collect FOV coords
        x_data = np.array(all_metadata.x_data)
        x_data = x_data[:,np.newaxis]
        y_data = np.array(all_metadata.x_data)
        y_data = y_data[:,np.newaxis]
        z_data = np.array(all_metadata.z_data)
        z_data = z_data[:,np.newaxis]
        all_coords = np.hstack((z_data,x_data,y_data))
        fov_coords = all_coords[0::len(z_levels),:]
        
        tag_name = experiment_name + '_' + hybridization_name + '_' + channel
        
        
        # Save the file as zarr
        store = zarr.DirectoryStore(parsed_raw_data_fpath)
        root = zarr.group(store=store,overwrite=False)

        
        nd2fh.bundle_axes = 'zyx'
        # set iteration over the fields of view
        nd2fh.iter_axes = 'v'
        
        # Save coords of the FOV
        rows = np.arange(img_width)
        cols = np.arange(img_height)
        hybridization_num = int(hybridization_name.split('Hybridization')[-1])
        for fov in fields_of_view:
            img = np.array(nd2fh[fov],dtype=np.uint16)      
            array_name = tag_name + '_fov_' + str(fov)
            dgrp = root.create_group(array_name)
            fov_name = 'raw_data_fov_' + str(fov)
            dset = dgrp.create_dataset(fov_name, data=img, shape=img.shape, chunks=(1,None,None),overwrite=True)
            dset.attrs['grp_name'] = array_name
            dset.attrs['fov_name'] = fov_name
            dset.attrs['channel'] = channel
            dset.attrs['target_name'] = info_data['channels'][channel]
            dset.attrs['img_height'] = img_height
            dset.attrs['img_width'] = img_width
            dset.attrs['pixel_microns'] = pixel_microns
            dset.attrs['z_levels'] = list(z_levels)
            dset.attrs['fov_num'] = fov
            dset.attrs['StitchingChannel'] = info_data['StitchingChannel']
            dset.attrs['hybridization_num'] = hybridization_num
            dset.attrs['experiment_name'] = experiment_name
               
        # Rename the nd2 files
        # new_file_name = tag_name + '.nd2'
        # new_file_path = raw_files_dir / new_file_name
        # nd2_file_path.rename(new_file_path)
        # nd2_file_path = new_file_path
        
        # # Copy the pkl files
        # new_file_name = tag_name + '_info.pkl'
        # new_file_path = raw_files_dir / new_file_name
        # # Must copy the pkl file in order to be able to use the file for the other channels
        # shutil.copy(str(info_file), str(new_file_path))
        
        # Save the fov_coords
        fname = experiment_fpath / 'tmp' / (tag_name + '_fovs_coords.npy')
        np.save(fname, fov_coords)
