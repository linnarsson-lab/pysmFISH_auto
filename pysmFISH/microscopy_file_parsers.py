
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
from dask import array as da
from pathlib import Path


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
    # assert '_auto' in experiment_fpath.stem, sys.exit('no _auto in the experiment name')

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
        new_file_name = tag_name + '.nd2'
        logger.debug(f'.nd2 file renamed {new_file_name}')
        new_file_path = raw_files_dir / new_file_name
        nd2_file_path.rename(new_file_path)
        nd2_file_path = new_file_path
        
        # Copy the pkl files
        new_file_name = tag_name + '_info.pkl'
        logger.debug(f'info data file renamed {new_file_name}')
        new_file_path = raw_files_dir / new_file_name
        # Must copy the pkl file in order to be able to use the file for the other channels
        shutil.copy(str(info_file), str(new_file_path))
        


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
    # experiment_name = experiment_name.split('_auto')[0]
    
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
        new_exp_location = parsed_raw_data_fpath.parent
        fname = new_exp_location / 'tmp' / 'microscope_tiles_coords' / (tag_name + '_fovs_coords.npy')
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






def nikon_nd2_autoparser_xarray_zarr(nd2_file_path, parsed_raw_data_fpath, experiment_info):

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

    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
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
        
    

        # Save the fov_coords
        fname = experiment_fpath / 'microscope_tiles_coords' / (tag_name + '_fovs_coords.npy')
        np.save(fname, fov_coords)

        nd2fh.bundle_axes = 'zyx'
        # set iteration over the fields of view
        nd2fh.iter_axes = 'v'
        
        round_num = int(nd2_file_path.stem.split('_')[-2].split('Hybridization')[-1])
        z = np.arange(nd2fh.sizes['z'])
        r = np.arange(nd2fh.sizes['y'])
        c = np.arange(nd2fh.sizes['x'])
        channel_arr = np.array([channel])
        round_num = np.array([round_num])
        experiment_name_arr = np.array([experiment_name])
        chunks_dict = {'experiment_name':1,'channel':1,'fov':1,'round_num':1,
                       'zstack':nd2fh.sizes['z'],'rows':nd2fh.sizes['y'],'columns':nd2fh.sizes['x']}
        chunks_tuple = (1,1,1,1,nd2fh.sizes['z'],nd2fh.sizes['y'],nd2fh.sizes['x'])
        
        hybridization_num = int(hybridization_name.split('Hybridization')[-1])
        for fov in fields_of_view:
    
            array_name = tag_name + '_fov_' + str(fov)
            attrs = {}
            attrs['grp_name'] = array_name
            attrs['fov_name'] = 'fov_' + str(fov) 
            attrs['channel'] = channel
            attrs['target_name'] = info_data['channels'][channel]
            attrs['img_width'] = parsed_metadata['width']
            attrs['img_height'] = parsed_metadata['height']
            attrs['pixel_microns'] = pixel_microns
            attrs['z_levels'] = list(z_levels)
            attrs['fields_of_view'] = list(fields_of_view)
            attrs['fov_num'] = fov
            attrs['stitching_channel'] = experiment_info['StitchingChannel']
            attrs['stitching_type'] = experiment_info['Stitching_type']
            attrs['experiment_type'] = experiment_info['Experiment_type']
            attrs['hybridization_num'] = hybridization_num
            attrs['experiment_name'] = experiment_name
            attrs['fov_acquisition_coords_x'] = fov_coords[fov,1]
            attrs['fov_acquisition_coords_y'] = fov_coords[fov,2]
            attrs['fov_acquisition_coords_z'] = fov_coords[fov,0]


            if info_data['StitchingChannel'] == channel:
                attrs['processing_type'] = attrs['stitching_type']
            elif '_ST' in attrs['target_name']:
                attrs['processing_type'] = 'staining'
            else:
                attrs['processing_type'] = 'fish'


            fov_arr = np.array([fov])
            arr = np.array(nd2fh[fov],dtype=np.uint16)
            arr = arr[np.newaxis,np.newaxis,np.newaxis,np.newaxis, :,:,:]
            arr = da.from_array(arr, chunks=(None,None,None,None,None,None,None))
            data_xr =  xr.DataArray(arr,
                                   coords = {
                                       'experiment_name':experiment_name_arr,
                                       'channel':channel_arr, 
                                       'fov': fov_arr, 
                                       'round_num':round_num,
                                       'zstack':z,
                                       'rows':r,
                                       'columns':c},
                                    dims=['experiment_name','channel','fov','round_num',
                                          'zstack','rows','columns'],
                                    attrs= attrs
                                  )
            data_xr = data_xr.chunk(chunks_dict)
            data_dset_name = tag_name + '_fov_' + str(fov)
            data_dset_path = (parsed_raw_data_fpath / (data_dset_name + '.zarr')).as_posix()
#             data_dset = xr.Dataset({data_dset_name: data_xr})
            data_dset = xr.Dataset({'raw_data': data_xr},attrs=attrs)
            data_dset.to_zarr(data_dset_path)
#             data_dset.to_netcdf(parsed_raw_data_fpath + '/' + data_dset_name + '.nc')

        # Rename the nd2 files
        new_file_name = tag_name + '.nd2'
        logger.debug(f'.nd2 file renamed {new_file_name}')
        new_file_path = raw_files_dir / new_file_name
        nd2_file_path.rename(new_file_path)
        nd2_file_path = new_file_path
        
        # Copy the pkl files
        new_file_name = tag_name + '_info.pkl'
        logger.debug(f'info data file renamed {new_file_name}')
        new_file_path = raw_files_dir / new_file_name
        # Must copy the pkl file in order to be able to use the file for the other channels
        shutil.copy(str(info_file), str(new_file_path))
        



def nikon_nd2_reparser_xarray_zarr(nd2_file_path,parsed_raw_data_fpath,experiment_info):
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

    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)

    logger.debug(f'processing file {nd2_fname}')

    experiment_fpath = nd2_file_path.parent.parent
    experiment_name = nd2_file_path.parent.parent.stem
    # experiment_name = experiment_name.split('_auto')[0]
    
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
        new_exp_location = parsed_raw_data_fpath.parent
        fname = new_exp_location /  'microscope_tiles_coords' / (tag_name + '_fovs_coords.npy')
        np.save(fname, fov_coords)


        # Save the file as zarr
        store = zarr.DirectoryStore(parsed_raw_data_fpath)
        root = zarr.group(store=store,overwrite=False)

        
        nd2fh.bundle_axes = 'zyx'
        # set iteration over the fields of view
        nd2fh.iter_axes = 'v'
        
        round_num = int(nd2_file_path.stem.split('_')[-2].split('Hybridization')[-1])
        z = np.arange(nd2fh.sizes['z'])
        r = np.arange(nd2fh.sizes['y'])
        c = np.arange(nd2fh.sizes['x'])
        channel_arr = np.array([channel])
        round_num = np.array([round_num])
        experiment_name_arr = np.array([experiment_name])
        chunks_dict = {'experiment_name':1,'channel':1,'fov':1,'round_num':1,
                       'zstack':nd2fh.sizes['z'],'rows':nd2fh.sizes['y'],'columns':nd2fh.sizes['x']}
        chunks_tuple = (1,1,1,1,nd2fh.sizes['z'],nd2fh.sizes['y'],nd2fh.sizes['x'])
        
        hybridization_num = int(hybridization_name.split('Hybridization')[-1])
        for fov in fields_of_view:
    
            array_name = tag_name + '_fov_' + str(fov)
            attrs = {}
            attrs['grp_name'] = array_name
            attrs['fov_name'] = 'fov_' + str(fov) 
            attrs['channel'] = channel
            attrs['target_name'] = info_data['channels'][channel]
            attrs['img_width'] = parsed_metadata['width']
            attrs['img_height'] = parsed_metadata['height']
            attrs['pixel_microns'] = pixel_microns
            attrs['z_levels'] = list(z_levels)
            attrs['fields_of_view'] = list(fields_of_view)
            attrs['fov_num'] = fov
            attrs['stitching_channel'] = experiment_info['StitchingChannel']
            attrs['stitching_type'] = experiment_info['Stitching_type']
            attrs['experiment_type'] = experiment_info['Experiment_type']
            attrs['hybridization_num'] = hybridization_num
            attrs['experiment_name'] = experiment_name
            attrs['fov_acquisition_coords_x'] = fov_coords[fov,1]
            attrs['fov_acquisition_coords_y'] = fov_coords[fov,2]
            attrs['fov_acquisition_coords_z'] = fov_coords[fov,0]


            if info_data['StitchingChannel'] == channel:
                attrs['processing_type'] = attrs['stitching_type']
            elif '_ST' in attrs['target_name']:
                attrs['processing_type'] = 'staining'
            else:
                attrs['processing_type'] = 'fish'


            fov_arr = np.array([fov])
            arr = np.array(nd2fh[fov],dtype=np.uint16)
            arr = arr[np.newaxis,np.newaxis,np.newaxis,np.newaxis, :,:,:]
            arr = da.from_array(arr, chunks=(None,None,None,None,None,None,None))
            data_xr =  xr.DataArray(arr,
                                   coords = {
                                       'experiment_name':experiment_name_arr,
                                       'channel':channel_arr, 
                                       'fov': fov_arr, 
                                       'round_num':round_num,
                                       'zstack':z,
                                       'rows':r,
                                       'columns':c},
                                    dims=['experiment_name','channel','fov','round_num',
                                          'zstack','rows','columns'],
                                    attrs= attrs
                                  )
            data_xr = data_xr.chunk(chunks_dict)
            data_dset_name = tag_name + '_fov_' + str(fov)
            data_dset_path = (parsed_raw_data_fpath / (data_dset_name + '.zarr')).as_posix()
#             data_dset = xr.Dataset({data_dset_name: data_xr})
            data_dset = xr.Dataset({'raw_data': data_xr},attrs=attrs)
            data_dset.to_zarr(data_dset_path)











def single_nikon_nd2_parser_simple(nd2_file_path):
    """
    Utility class used to parse a single .nd2 file
    The output will be a .zarr file with
    root/frame/channel/fov and the image chunked on the z-level
    This is for small files, i do not make used of the ram
    """

    nd2_file_path = Path(nd2_file_path)
    fresh_nuclei_dpath = nd2_file_path.parent
    logger = selected_logger()
    parsed_fpath = fresh_nuclei_dpath / (nd2_file_path.stem + '.zarr')
    parse_st = zarr.DirectoryStore(parsed_fpath.as_posix())
    parsed_root = zarr.group(store=parse_st,overwrite=True)

    with nd2reader.ND2Reader(nd2_file_path.as_posix()) as nd2fh:
        
        # Collect metadata
        all_metadata = nd2fh.parser._raw_metadata
        parsed_metadata = nd2fh.parser._raw_metadata.get_parsed_metadata()
                
        if 'z' in nd2fh.axes:
            nd2fh.bundle_axes = 'zyx'

            # Collect FOV coords
            z_levels = parsed_metadata['z_levels']
            x_data = np.array(all_metadata.x_data)
            x_data = x_data[:,np.newaxis]
            y_data = np.array(all_metadata.y_data)
            y_data = y_data[:,np.newaxis]
            z_data = np.array(all_metadata.z_data)
            z_data = z_data[:,np.newaxis]
            all_coords = np.hstack((z_data,x_data,y_data))
            fov_coords = all_coords[0::len(z_levels),:]

        else:
            nd2fh.bundle_axes = 'yx'
            
            # Collect FOV coords
            x_data = np.array(all_metadata.x_data)
            x_data = x_data[:,np.newaxis]
            y_data = np.array(all_metadata.y_data)
            y_data = y_data[:,np.newaxis]
            fov_coords = np.hstack((x_data,y_data))

        # Save the fov_coords
        fname = nd2_file_path.parent / (nd2_file_path.stem + '_fovs_coords.npy')
        np.save(fname, fov_coords)

        # Save metadata
        fname = nd2_file_path.parent / (nd2_file_path.stem + '_parsed_metadata.pkl')
        pickle.dump(parsed_metadata, open(fname,'wb'))

        axis_iter_order = []
        # Different type of conversion depending on the axis:
        if 'c' in nd2fh.axes:
            axis_iter_order.append('c')
        if 'v' in nd2fh.axes:
            axis_iter_order.append('v')

        nd2fh.iter_axes = axis_iter_order
        
        counter = 0            

        # loop through the frames (time dimension)
        for frame in nd2fh.metadata['frames']:
            for channel in nd2fh.metadata['channels']:
                for fov in nd2fh.metadata['fields_of_view']:
                    img = np.array(nd2fh.get_frame(counter),dtype=np.uint16)
                    if len(img.shape) == 3:
                        parsed_root.create_dataset(fov, data=img, shape=img.shape, chunks=(1,None,None))
                    elif len(img.shape) == 2:
                        parsed_root.create_dataset(fov, data=img, shape=img.shape, chunks=(None,None))
                    counter +=1