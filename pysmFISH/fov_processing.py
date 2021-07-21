import pickle
import gc
import dask
import sys
import pandas as pd
import numpy as np
import zarr
from typing import *
from pathlib import Path
from dask import delayed
from dask import dataframe as dd
from dask.base import tokenize


import pysmFISH
from pysmFISH import io
from pysmFISH import fovs_registration
from pysmFISH import barcodes_analysis
from pysmFISH import stitching
from pysmFISH import preprocessing
from pysmFISH import configuration_files
from pysmFISH import utils
from pysmFISH import microscopy_file_parsers
from pysmFISH import data_models
from pysmFISH.logger_utils import selected_logger


def combine_filtered_images(output_list: list,experiment_fpath: str,
                            metadata: pd.DataFrame, save:bool=False):
    """Function used to combine all the filtered images for a fov/channel in a single
        image stack

    Args:
        output_list (list): list containing the output of preprocessing 
        experiment_fpath (str): path to the experiment to process
        metadata (pd.DataFrame): dataframe containing the metadata
        save (bool, optional): Determine if the filtered images should be stored Defaults to False.

    Returns:
        img_stack (np.ndarray): image stack of all the images for a fov. The position in the
                stack correspond to round_num-1
    """
    experiment_fpath = Path(experiment_fpath)
     
    img_stack = np.zeros([metadata['total_rounds'],metadata['img_width'],metadata['img_height']])

    for img, img_meta in output_list:
        round_num = img_meta.round_num
        img_stack[round_num-1,:,:] = img

    if save:
        # Add conversion to more compress ftype
        img_meta = output_list[0][1]
        channel = img_meta.channel
        fov = img_meta.fov_num
        fpath = experiment_fpath / 'results' / (experiment_fpath.stem + '_' + channel + '_combined_img_fov_' + fov + '.npy')
        np.save(fpath, img_stack)
    
    return img_stack



def combine_steps(*args):
    pass


def single_fov_round_processing_eel(fov_subdataset,
                                   analysis_parameters,
                                   running_functions,
                                   dark_img,
                                   experiment_fpath,
                                   preprocessed_zarr_fpath,
                                   save_steps_output=False):

    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)

    # Path of directory where to save the intermediate results
    filtered_img_path = experiment_fpath / 'tmp' / 'filtered_images'
    raw_counts_path = experiment_fpath / 'tmp' / 'raw_counts'
    
    experiment_name = fov_subdataset.experiment_name
    pipeline = fov_subdataset.pipeline
    processing_type = fov_subdataset.processing_type
    zarr_grp_name = fov_subdataset.grp_name
    
    raw_data_location = Path(fov_subdataset.raw_data_location)
    parsed_raw_data_fpath = raw_data_location.parent

    if processing_type == 'fish':
        processing_parameters = analysis_parameters['fish']
        filtering_fun = running_functions['fish_channels_preprocessing']
        counting_fun = running_functions['fish_channels_dots_calling']
    
    elif 'beads' in processing_type:
        processing_parameters = analysis_parameters[processing_type]
        filtering_fun = running_functions['reference_channels_preprocessing']
        counting_fun = running_functions['reference_channels_dots_calling']

    elif processing_type != 'staining':
        processing_parameters = analysis_parameters[processing_type]
        filtering_fun = running_functions['reference_channels_preprocessing']
        counting_fun = running_functions['reference_channels_dots_calling']


    filt_out = getattr(pysmFISH.preprocessing,filtering_fun)(
                                                    zarr_grp_name,
                                                    parsed_raw_data_fpath,
                                                    processing_parameters,
                                                    dark_img)

    counts = getattr(pysmFISH.dots_calling,counting_fun)(
                                                        filt_out[0][0],
                                                        fov_subdataset,
                                                        processing_parameters)                                              

    if save_steps_output:

        # Save the file as zarr
        store = zarr.DirectoryStore(preprocessed_zarr_fpath)
        root = zarr.group(store=store,overwrite=False)
        tag_name = experiment_name + '_' + fov_subdataset.channel + '_round_' + str(fov_subdataset.round_num) + '_fov_' + str(fov_subdataset.fov_num)
        dgrp = root.create_group(tag_name,overwrite=True)
        for k, v in filt_out[1].items():
            dgrp.attrs[k] = v
        fov_name = 'preprocessed_data_fov_' + str(fov_subdataset.fov_num)
        dgrp.attrs['fov_name'] = fov_name
        img = utils.convert_to_uint16(filt_out[0][-1])
        dset = dgrp.create_dataset(fov_name, data=img, shape=img.shape, chunks=None,overwrite=True)

        # counts.to_parquet(raw_counts_path / (fname + '.parquet'),index=False)

    # return counts, (fov_subdataset.channel,fov_subdataset.round_num,img)
    return counts, filt_out


def single_fov_round_processing_serial_nuclei(fov_subdataset,
                                   analysis_parameters,
                                   running_functions,
                                   dark_img,
                                   experiment_fpath,
                                   preprocessed_image_tag,
                                   preprocessed_zarr_fpath,
                                   save_steps_output=False):

    """
    Process the nuclei for registration
    """
    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)

    # Path of directory where to save the intermediate results
    filtered_img_path = experiment_fpath / 'tmp' / 'filtered_images'
    raw_counts_path = experiment_fpath / 'tmp' / 'raw_counts'
    
    experiment_name = fov_subdataset.experiment_name
    pipeline = fov_subdataset.pipeline
    processing_type = fov_subdataset.processing_type
    zarr_grp_name = fov_subdataset.grp_name
    
    raw_data_location = Path(fov_subdataset.raw_data_location)
    parsed_raw_data_fpath = raw_data_location.parent

    processing_parameters = analysis_parameters[processing_type]

    filt_out = getattr(pysmFISH.preprocessing,running_functions['reference_channels_preprocessing'])(
                                                                        zarr_grp_name,
                                                                        parsed_raw_data_fpath,
                                                                        processing_parameters)

    if save_steps_output:

        # Save the file as zarr
        store = zarr.DirectoryStore(preprocessed_zarr_fpath)
        root = zarr.group(store=store,overwrite=False)
        tag_name = experiment_name + '_' + fov_subdataset.channel + '_round_' + str(fov_subdataset.round_num) + '_fov_' + str(fov_subdataset.fov_num)
        dgrp = root.create_group(tag_name,overwrite=True)
        for k, v in filt_out[1].items():
            dgrp.attrs[k] = v
        fov_name = 'preprocessed_data_fov_' + str(fov_subdataset.fov_num)
        dgrp.attrs['fov_name'] = fov_name
        img = utils.convert_to_uint16(filt_out[0][-1])
        dset = dgrp.create_dataset(fov_name, data=img, shape=img.shape, chunks=None,overwrite=True)

    return (img,fov_subdataset)


def processing_barcoded_eel_fov_graph(experiment_fpath,analysis_parameters,
                                    running_functions, tiles_org,metadata,
                                    grpd_fovs,save_intermediate_steps, 
                                    preprocessed_image_tag, client, chunks_size, save_bits_int):
        """ 
        This method create a processing graph XXXXXXXXX
        Need to chunks the processing of the data because getting the intensity dots
        is slow and if i process everything in parallel the workers can get lost and
        the processing will crush
        """
        experiment_fpath = Path(experiment_fpath)
        io.create_empty_zarr_file(experiment_fpath.as_posix(), preprocessed_image_tag)
        preprocessed_zarr_fpath = experiment_fpath / (experiment_fpath.stem + '_' + preprocessed_image_tag + '.zarr')

        microscopy_file_parsers.create_dark_img(experiment_fpath,metadata)

        dark_img = preprocessing.load_dark_image(experiment_fpath)
        
        # did this conversion to avoid to pass self to dask
        # analysis_parameters = analysis_parameters
        # running_functions = running_functions
        # tile_corners_coords_pxl = tile_corners_coords_pxl
        
        
        dark_img = delayed(dark_img)
        analysis_parameters = delayed(analysis_parameters)
        running_functions = delayed(running_functions)
        tile_corners_coords_pxl = delayed(tiles_org.tile_corners_coords_pxl)

        list_all_channels = metadata['list_all_channels']
        stitching_channel = metadata['stitching_channel']
        fish_channels = set(list_all_channels).difference(stitching_channel)

        codebook_dict = configuration_files.load_codebook(experiment_fpath,metadata)
        codebook_df = delayed(codebook_dict)
        
        all_processing = []
        all_filtered_images = []
        all_fovs_channels = list(grpd_fovs.groups.keys())
        chunks = [all_fovs_channels[x:x+chunks_size] for x in range(0, len(all_fovs_channels), chunks_size)]
        for chunk in chunks:
            all_processing = []
            all_filtered_images = []
            all_counts_fov = {}
            for channel in list_all_channels:
                all_counts_fov[channel] = []
            for fov_channel in chunk:
                fov_num = fov_channel[0]
                group = grpd_fovs.get_group(fov_channel)

                for index_value, fov_subdataset in group.iterrows():
                    round_num = fov_subdataset.round_num
                    channel = fov_subdataset.channel
                    fov = fov_subdataset.fov_num
                    experiment_name = fov_subdataset.experiment_name
                    dask_delayed_name = 'filt_count_' +experiment_name + '_' + channel + \
                                    '_round_' + str(round_num) + '_fov_' +str(fov) + '-' + tokenize()
                    fov_out = delayed(single_fov_round_processing_eel, name=dask_delayed_name,nout=2)(fov_subdataset,
                                                analysis_parameters,
                                                running_functions,
                                                dark_img,
                                                experiment_fpath,
                                                preprocessed_zarr_fpath,
                                                save_steps_output=save_intermediate_steps,
                                                dask_key_name=dask_delayed_name)
                    counts, filt_out = fov_out[0], fov_out[1]
                    
                    all_counts_fov[channel].append(counts)
                    
                    if save_bits_int:
                        if channel != fov_subdataset.stitching_channel:
                            all_filtered_images.append(filt_out)
                
                all_counts_fov_concat = {}
                for processing_channel in list_all_channels:
                    name = 'concat_' +experiment_name + '_' + channel + '_' \
                                        + '_fov_' +str(fov) + '-' + tokenize()
                    all_counts_fov_concat[channel] = delayed(pd.concat,name=name)(all_counts_fov[channel],axis=0,ignore_index=True)
                
                if save_intermediate_steps:
                    
                        name = 'save_raw_counts_' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize()
                        saved_raw_counts = delayed(pickle.dump,name=name)(all_counts_fov_concat, Path(experiment_fpath) / 'results'/ (experiment_name + \
                                '_raw_fov_' + str(fov) + '.pkl'))

                        all_processing.append(saved_raw_counts)


                name = 'register_' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize()
                registered_counts = delayed(fovs_registration.beads_based_registration,name=name)(all_counts_fov,
                                                    analysis_parameters)

                
                registration_stitching_channel_output = delayed(fovs_registration.beads_based_registration_stitching_channel,name=name)(all_counts_fov[stitching_channel],
                                                        analysis_parameters)

                stitching_channel_df, all_rounds_shifts, all_rounds_matching_dots = registration_stitching_channel_output[0], \
                                                                                    registration_stitching_channel_output[1], \
                                                                                    registration_stitching_channel_output[2]


                all_stitched_coords = []

                for processing_channel in fish_channels:
                    
                    # Register fish
                    name = 'register_fish_channels_' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize()

                    registered_counts = delayed(fovs_registration.beads_based_registration_fish,name=name)(all_counts_fov[processing_channel],
                                                        all_rounds_shifts, all_rounds_matching_dots, analysis_parameters)

                    # Decoded fish
                    name = 'decode_' +experiment_name + '_' + channel + '_' \
                                        + '_fov_' +str(fov) + '-' + tokenize()
                    decoded = delayed(barcodes_analysis.extract_barcodes_NN_fast_multicolor,name=name)(registered_counts, 
                                                                              analysis_parameters,codebook_df)                                                        
                
                    # Stitch to the microscope reference coords
                    name = 'stitch_to_mic_coords_' +experiment_name + '_' + channel + '_' \
                                        + '_fov_' +str(fov) + '-' + tokenize()  
                    stitched_coords = delayed(stitching.stitch_using_microscope_fov_coords,name=name)(decoded[1],tile_corners_coords_pxl)
                
                    all_stitched_coords.append(stitched_coords)

                
                all_stitched_coords.append(stitching_channel_df)

                
                name = 'concat_' +experiment_name + '_' + channel + '_' \
                                        + '_fov_' +str(fov) + '-' + tokenize()
                all_stitched_coords = delayed(pd.concat,name=name)(all_stitched_coords,axis=0,ignore_index=True) 
                    
                    
                name = 'save_df_' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize() 
                saved_file = delayed(stitched_coords.to_parquet,name=name)(Path(experiment_fpath) / 'results'/ (experiment_name + \
                                '_decoded_fov_' + str(fov) + '.parquet'),index=False)
            
            
                all_processing.append(saved_file)
            
            _ = dask.compute(*all_processing)

        io.consolidate_zarr_metadata(preprocessed_zarr_fpath)
 
        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        selected_Hdistance = 3 / metadata['barcode_length']
        stitching_selected = 'microscope_stitched'
        io.simple_output_plotting(experiment_fpath, stitching_selected, selected_Hdistance, client,file_tag='decoded')
        # ----------------------------------------------------------------  


def processing_barcoded_eel_fov_after_dots_graph(experiment_fpath,analysis_parameters,
                                    running_functions, tiles_org,metadata,
                                    grpd_fovs, 
                                    preprocessed_image_tag, client, chunks_size):
        """ 
        This method create a processing graph to run the registration and t
        
        """
        experiment_fpath = Path(experiment_fpath)
        experiment_name = experiment_fpath.stem
    
        preprocessed_zarr_fpath = experiment_fpath / (experiment_fpath.stem + '_' + preprocessed_image_tag + '.zarr')

        analysis_parameters = delayed(analysis_parameters)
        running_functions = delayed(running_functions)
        tile_corners_coords_pxl = delayed(tiles_org.tile_corners_coords_pxl)

        codebook = configuration_files.load_codebook(experiment_fpath,metadata)
        codebook_df = delayed(codebook)
        

        all_processing = []

        all_fovs = list(grpd_fovs.groups.keys())
        chunks = [all_fovs[x:x+chunks_size] for x in range(0, len(all_fovs), chunks_size)]
        for chunk in chunks:
            all_processing = []
            for fov_num in chunk:
                
                # Modify for channels name

                counts_fpath = list((experiment_fpath / 'results').glob('*_raw_fov_'+str(fov_num)+'.parquet'))[0]

                name = 'load_counts_' +experiment_name + '_' \
                                    + '_fov_' +str(fov_num) + '-' + tokenize()
                all_counts_fov = delayed(pd.read_parquet,name=name)(counts_fpath)
                
                name = 'register_' +experiment_name + '_' \
                                    + '_fov_' +str(fov_num) + '-' + tokenize()
                registered_counts = delayed(fovs_registration.beads_based_registration,name=name)(all_counts_fov,
                                                    analysis_parameters)

                name = 'decode_' +experiment_name + '_' \
                                    + '_fov_' +str(fov_num) + '-' + tokenize()
                decoded = delayed(barcodes_analysis.extract_barcodes_NN_fast,name=name)(registered_counts, 
                                                                            analysis_parameters,codebook_df)                                                        
                
                name = 'stitch_to_mic_coords_' +experiment_name + '_' \
                                    + '_fov_' +str(fov_num) + '-' + tokenize()  


                stitched_coords = delayed(stitching.stitch_using_coords_general_df,name=name)(decoded[1],tile_corners_coords_pxl,
                                                            tiles_org.reference_corner_fov_position,
                                                            metadata, tag='microscope_stitched')
                
                name = 'save_df_' +experiment_name + '_' \
                                    + '_fov_' +str(fov_num) + '-' + tokenize() 
                saved_file = delayed(stitched_coords.to_parquet,name=name)(Path(experiment_fpath) / 'results'/ (experiment_name + \
                                '_decoded_fov_' + str(fov_num) + '.parquet'),index=False)
            
                all_processing.append(saved_file)
            
            _ = dask.compute(*all_processing)

        io.consolidate_zarr_metadata(preprocessed_zarr_fpath)
 
        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        selected_Hdistance = 3 / metadata['barcode_length']
        stitching_selected = 'microscope_stitched'
        io.simple_output_plotting(experiment_fpath, stitching_selected, selected_Hdistance, client,file_tag='decoded')
        # ----------------------------------------------------------------  




def processing_serial_fish_fov_graph(experiment_fpath,analysis_parameters,
                                    running_functions, tiles_org,metadata,
                                    grpd_fovs,save_intermediate_steps, 
                                    preprocessed_image_tag, client,chunks_size):
        """ 
        This method create a processing graph XXXXXXXXX

        """

        experiment_fpath = Path(experiment_fpath)
        io.create_empty_zarr_file(experiment_fpath, preprocessed_image_tag)
        preprocessed_zarr_fpath = experiment_fpath / (experiment_fpath.stem + '_' + preprocessed_image_tag + '.zarr')

        microscopy_file_parsers.create_dark_img(experiment_fpath,metadata)


        dark_img = preprocessing.load_dark_image(experiment_fpath)
        
        # did this conversion to avoid to pass self to dask
        analysis_parameters = analysis_parameters
        running_functions = running_functions
        tile_corners_coords_pxl = tiles_org.tile_corners_coords_pxl
        
        dark_img = delayed(dark_img)
        analysis_parameters = delayed(analysis_parameters)
        running_functions = delayed(running_functions)
        tile_corners_coords_pxl = delayed(tile_corners_coords_pxl)

        all_processing = []
        all_filtered_images = []
        all_fovs = list(grpd_fovs.groups.keys())
        chunks = [all_fovs[x:x+chunks_size] for x in range(0, len(all_fovs), chunks_size)]
        for chunk in chunks:
            all_processing = []
            all_filtered_images = []
            for fov_num in chunk:
                group = grpd_fovs.get_group(fov_num)
        # for fov_num, group in grpd_fovs:
                all_counts_fov = []
                all_nuclei_fov = []
                for index_value, fov_subdataset in group.iterrows():
                    round_num = fov_subdataset.round_num
                    channel = fov_subdataset.channel
                    fov = fov_subdataset.fov_num
                    stitching_type = fov_subdataset.stitching_type
                    experiment_name = fov_subdataset.experiment_name
                    processing_type = fov_subdataset.processing_type

                    if processing_type == 'nuclei':
                        dask_delayed_name = 'filt_' +experiment_name + '_' + channel + \
                                        '_round_' + str(round_num) + '_fov_' +str(fov) + '-' + tokenize()

                        out_nuclei = delayed(single_fov_round_processing_serial_nuclei,name=dask_delayed_name)(fov_subdataset,
                                                analysis_parameters,
                                                running_functions,
                                                dark_img,
                                                experiment_fpath,
                                                preprocessed_image_tag,
                                                preprocessed_zarr_fpath,
                                                save_steps_output=save_intermediate_steps,
                                                dask_key_name=dask_delayed_name)
                        all_nuclei_fov.append(out_nuclei)


                    else:
                        dask_delayed_name = 'filt_count_' +experiment_name + '_' + channel + \
                                        '_round_' + str(round_num) + '_fov_' +str(fov) + '-' + tokenize()
                        fov_out = delayed(single_fov_round_processing_eel,name=dask_delayed_name)(fov_subdataset,
                                                    analysis_parameters,
                                                    running_functions,
                                                    dark_img,
                                                    experiment_fpath,
                                                    preprocessed_zarr_fpath,
                                                    save_steps_output=save_intermediate_steps,
                                                    dask_key_name=dask_delayed_name)
                        
                        counts, filt_out = fov_out[0], fov_out[1]
                        all_counts_fov.append(counts)
                        # if channel != fov_subdataset.stitching_channel:
                        #     all_filtered_images.append(filt_out)
                        

                name = 'concat_fish_' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize()
                all_counts_fov = delayed(pd.concat,name=name)(all_counts_fov,axis=0,ignore_index=True)
                
                if stitching_type == 'nuclei':

                    name = 'create_nuclei_stack' +experiment_name + '_' + channel + '_' \
                                        + '_fov_' +str(fov) + '-' + tokenize()
                    filtered_nuclei_stack = delayed(combine_filtered_images,name=name)(all_nuclei_fov,experiment_fpath,metadata)

                    name = 'register_' +experiment_name + '_' + channel + '_' \
                                        + '_fov_' +str(fov) + '-' + tokenize()
                    registered_counts = delayed(fovs_registration.nuclei_based_registration,name=name)(all_counts_fov,
                                                        filtered_nuclei_stack,
                                                        analysis_parameters)

                else:

                    name = 'register_' +experiment_name + '_' + channel + '_' \
                                        + '_fov_' +str(fov) + '-' + tokenize()
                    registered_counts = delayed(fovs_registration.beads_based_registration,name=name)(all_counts_fov,
                                                        analysis_parameters)
                                                                                                    
                name = 'stitch_to_mic_coords_' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize()  

                stitched_coords = delayed(stitching.stitch_using_microscope_fov_coords,name=name)(registered_counts,tile_corners_coords_pxl,
                                                            tiles_org.reference_corner_fov_position,
                                                            metadata, tag='microscope_stitched')
                
                name = 'register_and_combine_filt_imgs' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize() 
                
                # combined_images = delayed(fovs_registration.combine_register_filtered_images,name=name)(all_filtered_images,stitched_coords,
                #                                                                                 fov_subdataset.stitching_channel)

                name = 'save_df_' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize() 
                saved_file = delayed(stitched_coords.to_parquet,name=name)(Path(experiment_fpath) / 'results'/ (experiment_name + \
                                '_decoded_fov_' + str(fov) + '.parquet'),index=False)
            
                all_processing.append(saved_file) 
            
            _ = dask.compute(*all_processing)

        io.consolidate_zarr_metadata(preprocessed_zarr_fpath)

        # # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        stitching_selected = 'microscope_stitched'
        io.simple_output_plotting_serial(experiment_fpath, stitching_selected, client, file_tag='decoded')
        # ----------------------------------------------------------------  
  



def single_fov_fresh_tissue_beads(processing_tag,
                                   fov_subdataset,
                                   analysis_parameters,
                                   running_functions,
                                   dark_img,
                                   experiment_fpath,
                                   preprocessed_zarr_fpath,
                                   save_steps_output=True):

    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)

    experiment_name = fov_subdataset.experiment_name
    zarr_grp_name = fov_subdataset.grp_name

    parsed_raw_data_fpath = experiment_fpath / 'fresh_tissue'/ (experiment_name +'_img_data.zarr')
    
    processing_parameters = analysis_parameters['fresh-tissue'][processing_tag]

    if processing_tag == 'beads':
        filtering_fun = running_functions['fresh_sample_reference_preprocessing']
        counting_fun = running_functions['fresh_sample_reference_dots_calling']


        filt_out = getattr(pysmFISH.preprocessing,filtering_fun)(
                                                        zarr_grp_name,
                                                        parsed_raw_data_fpath,
                                                        processing_parameters,
                                                        dark_img)

        counts = getattr(pysmFISH.dots_calling,counting_fun)(
                                                            filt_out[0][0],
                                                            fov_subdataset,
                                                            processing_parameters)                                              

    elif processing_tag == 'nuclei':
        filtering_fun = running_functions['fresh_sample_nuclei_preprocessing']
        filt_out = getattr(pysmFISH.preprocessing,filtering_fun)(
                                                        zarr_grp_name,
                                                        parsed_raw_data_fpath,
                                                        processing_parameters)


    if save_steps_output:

        # Save the file as zarr
        store = zarr.DirectoryStore(preprocessed_zarr_fpath)
        root = zarr.group(store=store,overwrite=False)
        tag_name = experiment_name + '_fresh_tissue_' + processing_tag + '_fov_' + str(fov_subdataset.fov_num)
        dgrp = root.create_group(tag_name,overwrite=True)
        for k, v in filt_out[1].items():
            dgrp.attrs[k] = v
        fov_name = 'preprocessed_data_fov_' + str(fov_subdataset.fov_num)
        dgrp.attrs['fov_name'] = fov_name
        img = utils.convert_to_uint16(filt_out[0][-1])
        dset = dgrp.create_dataset(fov_name, data=img, shape=img.shape, chunks=None,overwrite=True)


    if processing_tag == 'beads':
        return counts, filt_out

    elif processing_tag == 'nuclei':
        return filt_out


def process_fresh_sample_graph(experiment_fpath, running_functions, 
                                analysis_parameters, client, tag_ref_beads, 
                                tag_nuclei,parsing=True):
    logger = selected_logger()
    all_parsing = []
    
    if parsing:
        presence_nuclei = 0
        presence_beads = 0
        all_fresh_tissue_fpath = list((Path(experiment_fpath) / 'fresh_tissue').glob('*.nd2'))
        if all_fresh_tissue_fpath:
            for fpath in all_fresh_tissue_fpath:
                if tag_ref_beads in fpath.stem:
                    parsed_beads_fpath = experiment_fpath / 'fresh_tissue'/ (fpath.stem +'_img_data.zarr')
                    parsing_future = client.submit(microscopy_file_parsers.nikon_nd2_parser_simple_mfov,fpath)
                    all_parsing.append(parsing_future)
                    presence_beads = 1
                elif tag_nuclei in fpath.stem:
                    parsed_nuclei_fpath = experiment_fpath / 'fresh_tissue'/ (fpath.stem +'_img_data.zarr')
                    parsing_future = client.submit(microscopy_file_parsers.nikon_nd2_parser_simple_mfov,fpath)
                    all_parsing.append(parsing_future)
                    presence_nuclei = 1

            if presence_beads and presence_nuclei:
                _ = client.gather(all_parsing)
                io.consolidate_zarr_metadata(parsed_beads_fpath)
                io.consolidate_zarr_metadata(parsed_nuclei_fpath)
            else:
                if not presence_beads:
                    logger.error(f'missing fresh-tissue beads file')
                    sys.exit(f'missing fresh-tissue beads file')
                elif not presence_nuclei:
                    logger.error(f'missing fresh-tissue nuclei file')
                    sys.exit(f'missing fresh-tissue nuclei file')
                else:
                    logger.error(f'missing fresh-tissue beads and nuclei files')
                    sys.exit(f'missing fresh-tissue beads and nuclei files')

    else:
        all_fresh_tissue_fpath = list((Path(experiment_fpath) / 'fresh_tissue').glob('*.zarr'))
        if all_fresh_tissue_fpath:
            for fpath in all_fresh_tissue_fpath:
                if tag_ref_beads in fpath.stem:
                    parsed_beads_fpath = fpath
                    presence_beads = 1
                elif tag_nuclei in fpath.stem:
                    parsed_nuclei_fpath = fpath
                    presence_nuclei = 1
              
        if not presence_beads:
            logger.error(f'missing fresh-tissue beads parsed file')
            sys.exit(f'missing fresh-tissue beads parsed file')
        elif not presence_nuclei:
            logger.error(f'missing fresh-tissue nuclei parsed file')
            sys.exit(f'missing fresh-tissue nuclei parsed file')
        elif (not presence_nuclei) & (not presence_beads    ):
            logger.error(f'missing fresh-tissue beads and nuclei parsed files')
            sys.exit(f'missing fresh-tissue beads and nuclei parsed files')  

  
    # Create dataset
    
    ds_beads = data_models.Dataset()
    ds_nuclei = data_models.Dataset()
    ds_beads.create_full_dataset_from_zmetadata(parsed_beads_fpath)
    ds_nuclei.create_full_dataset_from_zmetadata(parsed_nuclei_fpath)

    beads_grpd_fovs = ds_beads.dataset.groupby('fov_num')
    nuclei_grpd_fovs = ds_nuclei.dataset.groupby('fov_num')
    
    # In this case I fake a dark image. It must be collected from the 
    # robofish system
    img_width = ds_nuclei.dataset.iloc[0].img_width
    img_height = ds_nuclei.dataset.iloc[0].img_height
    dark_img = np.zeros([img_width,img_height])

    base_path = experiment_fpath / 'fresh_tissue'
    # nuclei_base = parsed_nuclei_fpath.stem.split('_img_data.zarr')[0]
    # beads_base = parsed_beads_fpath.stem.split('_img_data.zarr')[0]
    nuclei_filtered_fpath = base_path /  (base_path.stem + '_nuclei_preprocessed_img_data.zarr')
    io.create_empty_zarr_file(base_path.as_posix(), tag='nuclei_preprocessed_img_data')
    beads_filtered_fpath = base_path /  (base_path.stem + '_beads_preprocessed_img_data.zarr')
    io.create_empty_zarr_file(base_path.as_posix(), tag='beads_preprocessed_img_data')
    

    all_counts_beads = []
    all_processing_nuclei = []
    processing_tag = 'beads'
    for fov_num, group in beads_grpd_fovs:
        for index_value, fov_subdataset in group.iterrows():
            round_num = fov_subdataset.round_num
            channel = fov_subdataset.channel
            fov = fov_subdataset.fov_num
            experiment_name = fov_subdataset.experiment_name
            dask_delayed_name = 'filt_count_beads_fov'+ str(fov) + '_' + tokenize()
            fov_out = delayed(single_fov_fresh_tissue_beads, name=dask_delayed_name)(
                                            processing_tag,
                                            fov_subdataset,
                                            analysis_parameters,
                                            running_functions,
                                            dark_img,
                                            experiment_fpath,
                                            preprocessed_zarr_fpath=beads_filtered_fpath,
                                            save_steps_output=True,
                                            dask_key_name=dask_delayed_name)
            counts, filt_out = fov_out[0], fov_out[1]
            all_counts_beads.append(counts)
    

    name = 'concat_all_counts_beads_fresh_tissue'+ '-' + tokenize()
    all_counts_fov = delayed(pd.concat,name=name)(all_counts_beads,axis=0,ignore_index=True)

    # Add registration and recalculation of all the coords

    name = 'save_df_beads_fresh_tissue' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize() 
    saved_file = delayed(all_counts_fov.to_parquet,name=name)(Path(experiment_fpath) / 'results'/ (experiment_name + \
                    '_counts_beads_fresh_tissue.parquet'),index=False)

    processing_tag='nuclei'
    for fov_num, group in nuclei_grpd_fovs:
        for index_value, fov_subdataset in group.iterrows():
            round_num = fov_subdataset.round_num
            channel = fov_subdataset.channel
            fov = fov_subdataset.fov_num
            experiment_name = fov_subdataset.experiment_name
            dask_delayed_name = 'filt_nuclei_fov'+ str(fov) + '_' + tokenize()
            fov_out = delayed(single_fov_fresh_tissue_beads, name=dask_delayed_name)(
                                            processing_tag,
                                            fov_subdataset,
                                            analysis_parameters,
                                            running_functions,
                                            dark_img,
                                            experiment_fpath,
                                            preprocessed_zarr_fpath=nuclei_filtered_fpath,
                                            save_steps_output=True,
                                            dask_key_name=dask_delayed_name)
            
            all_processing_nuclei.append(fov_out)

 
    end = delayed(combine_steps)(saved_file,all_processing_nuclei)
      
    _ = dask.compute(end)
    


# def collect_bits_intensity_graph(dataset:pd.Dataframe, experiment_fpath:str, 
#                             grpd_fovs,metadata:Dict,chunks_size,codebooks,client):

#     # Write to add stitching reference to the dataset to make in easy to run the
#     # bits analysis
    

#     # Need to run for fov

#     bit_channels = list(set(metadata['list_all_channels']).difference(metadata['stitching_channel'])))
#     experiment_fpath = Path(experiment_fpath)
#     experiment_name = metadata['experiment_name']
#     filtered_images_path =  Path(experiment_fpath) / (metadata['experiment_name'] + 'preprocessed_img_data.zarr')

#     all_fovs = list(grpd_fovs.groups.keys())
#     chunks = [all_fovs[x:x+chunks_size] for x in range(0, len(all_fovs), chunks_size)]
#     for chunk in chunks:
#         all_processing = []
#         all_filtered_images = []
#         for fov_num in chunk:
#             group = grpd_fovs.get_group(fov_num)
#             grpd_channel = group.groupby('channel')

#             # Load counts
#             name = 'Load_counts' +experiment_name + '_' + \
#                                         + '_fov_' +str(fov_num) + '-' + tokenize() 
                    
#             counts_fpath = experiment_fpath / 'results' / (experiment_name + '_decoded_fov_' +str(fov_num) + '.parquet')
#             counts_df = delayed(pd.read_parquet,name=name)(counts_fpath)


#             for channel, channel_grp in grpd_channel:
#                 all_filtered_images = []
#                 for index_value, fov_subdataset in channel_grp.iterrows():
#                     # Create ((img,), metadata) list to match the one used in the eel graph in order
#                     # to used the same set of functions

#                     name =  'load_processed_images_' + fov_subdataset.grp_name+ '-' + tokenize()

#                     img = delayed(io.load_general_zarr,name=name)(fov_subdataset,filtered_images_path)
#                     filt_out = ((img,), metadata)
#                     all_filtered_images.append(filt_out)

                
#                 name = 'register_and_combine_filt_imgs' +experiment_name + '_' + channel + '_' \
#                                         + '_fov_' +str(fov_num) + '-' + tokenize() 
                    
#                 combined_images = delayed(fovs_registration.combine_register_filtered_images,name=name)(all_filtered_images,counts_df,
#                                                                                                 fov_subdataset.stitching_channel)


#                 name = 'extract_barcodes_int' +experiment_name + '_' + channel + '_' \
#                                     + '_fov_' +str(fov_num) + '-' + tokenize()
#                 barcodes_max = delayed(barcodes_analysis.extract_dots_images,name=name)(counts_df,combined_images,
#                                                                                         experiment_fpath)


#                 name = 'extract_bit_flip_direction' +experiment_name + '_' + channel + '_' \
#                                     + '_fov_' +str(fov) + '-' + tokenize()
#                 flip_direction = delayed(barcodes_analysis.define_flip_direction,name=name)(codebook,
#                                                                                         experiment_fpath,
#                                                                                         stitched_coords)



    # """Collected the intensity of the 1 and 0 bits and the flipping direction.
    #     This is an extra function to run after the data are collected to avoid
    #     to slow down the initial dots counting.

    # Args:
    #     dataset ([type]): [description]
    # """