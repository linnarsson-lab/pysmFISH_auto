import pickle
import gc
import dask
import pandas as pd
import numpy as np
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

from pysmFISH.logger_utils import selected_logger


def single_fov_round_processing_eel(fov_subdataset,
                                   analysis_parameters,
                                   running_functions,
                                   dark_img,
                                   experiment_fpath,
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
                                                        filt_out[0],
                                                        fov_subdataset,
                                                        processing_parameters)                                              

    if save_steps_output:
        fname = experiment_name + '_' + fov_subdataset.channel + '_round_' + str(fov_subdataset.round_num) + '_fov_' + str(fov_subdataset.fov_num)
        np.save(filtered_img_path / (fname + '.npy'),filt_out[-1] )
        counts.to_parquet(raw_counts_path / (fname + '.parquet'),index=False)

    # return counts, (fov_subdataset.channel,fov_subdataset.round_num,img)
    return counts


def single_fov_round_processing_serial_nuclei(fov_subdataset,
                                   analysis_parameters,
                                   running_functions,
                                   dark_img,
                                   experiment_fpath,
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

    img = getattr(pysmFISH.preprocessing,running_functions['reference_channels_preprocessing'])(
                                                                        zarr_grp_name,
                                                                        parsed_raw_data_fpath,
                                                                        processing_parameters)

    if save_steps_output:
        fname = experiment_name + '_' + fov_subdataset.channel + '_round_' + str(fov_subdataset.round_num) + '_fov_' + str(fov_subdataset.fov_num)
        np.save(filtered_img_path / (fname + '.npy'),img )

    return (img,fov_subdataset)


def processing_barcoded_eel_fov_graph(experiment_fpath,analysis_parameters,
                                    running_functions, tile_corners_coords_pxl,metadata,
                                    grpd_fovs,save_intermediate_steps, client):
        """ 
        This method create a processing graph XXXXXXXXX

        """
        dark_img = preprocessing.load_dark_image(experiment_fpath)
        
        # did this conversion to avoid to pass self to dask
        analysis_parameters = analysis_parameters
        running_functions = running_functions
        tile_corners_coords_pxl = tile_corners_coords_pxl
        
        
        dark_img = delayed(dark_img)
        analysis_parameters = delayed(analysis_parameters)
        running_functions = delayed(running_functions)
        tile_corners_coords_pxl = delayed(tile_corners_coords_pxl)

        codebook = configuration_files.load_codebook(experiment_fpath,metadata)
        codebook_df = delayed(codebook)
        
        all_processing = []
        
        for fov_num, group in grpd_fovs:
            all_counts_fov = []
            for index_value, fov_subdataset in group.iterrows():
                round_num = fov_subdataset.round_num
                channel = fov_subdataset.channel
                fov = fov_subdataset.fov_num
                experiment_name = fov_subdataset.experiment_name
                dask_delayed_name = 'filt_count_' +experiment_name + '_' + channel + \
                                '_round_' + str(round_num) + '_fov_' +str(fov) + '-' + tokenize()
                counts = delayed(single_fov_round_processing_eel, name=dask_delayed_name)(fov_subdataset,
                                            analysis_parameters,
                                            running_functions,
                                            dark_img,
                                            experiment_fpath,
                                            save_steps_output=save_intermediate_steps)
                all_counts_fov.append(counts)

            name = 'concat_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            all_counts_fov = delayed(pd.concat,name=name)(all_counts_fov,axis=0,ignore_index=True)
            
            name = 'register_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            registered_counts = delayed(fovs_registration.beads_based_registration,name=name)(all_counts_fov,
                                                analysis_parameters)

            name = 'decode_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            decoded = delayed(barcodes_analysis.extract_barcodes_NN_fast,name=name)(registered_counts, 
                                                                        analysis_parameters,codebook_df)                                                        
            
            name = 'stitch_to_mic_coords_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()  
            stitched_coords = delayed(stitching.stitch_using_microscope_fov_coords_new,name=name)(decoded[1],tile_corners_coords_pxl)
            
            name = 'save_df_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize() 
            saved_file = delayed(stitched_coords.to_parquet,name=name)(Path(experiment_fpath) / 'results'/ (experiment_name + \
                            '_decoded_fov_' + str(fov) + '.parquet'),index=False)
        
            all_processing.append(saved_file) 
        _ = dask.compute(*all_processing)
 
        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        selected_Hdistance = 3 / metadata['barcode_length']
        stitching_selected = 'microscope_stitched'
        io.simple_output_plotting(experiment_fpath, stitching_selected, selected_Hdistance, client)
        # ----------------------------------------------------------------  


def processing_serial_fish_fov_graph(experiment_fpath,analysis_parameters,
                                    running_functions, tile_corners_coords_pxl,metadata,
                                    grpd_fovs,save_intermediate_steps, client):
        """ 
        This method create a processing graph XXXXXXXXX

        """

        dark_img = preprocessing.load_dark_image(experiment_fpath)
        
        # did this conversion to avoid to pass self to dask
        analysis_parameters = analysis_parameters
        running_functions = running_functions
        tile_corners_coords_pxl = tile_corners_coords_pxl
        
        dark_img = delayed(dark_img)
        analysis_parameters = delayed(analysis_parameters)
        running_functions = delayed(running_functions)
        tile_corners_coords_pxl = delayed(tile_corners_coords_pxl)

        all_processing = []

        for fov_num, group in grpd_fovs:
            all_counts_fov = []
            all_nuclei_fov = []
            for index_value, fov_subdataset in group.iterrows():
                round_num = fov_subdataset.round_num
                channel = fov_subdataset.channel
                fov = fov_subdataset.fov_num
                experiment_name = fov_subdataset.experiment_name
                processing_type = fov_subdataset.processing_type

                if processing_type == 'nuclei':
                    dask_delayed_name = 'filt_' +experiment_name + '_' + channel + \
                                    '_round_' + str(round_num) + '_fov_' +str(fov) + '-' + tokenize()

                    out_nuclei = delayed(fov_processing.single_fov_round_processing_serial_nuclei,name=dask_delayed_name)(fov_subdataset,
                                            analysis_parameters,
                                            running_functions,
                                            dark_img,
                                            experiment_fpath,
                                            save_steps_output=save_intermediate_steps)
                    all_nuclei_fov.append(out_nuclei)

                else:
                    dask_delayed_name = 'filt_count_' +experiment_name + '_' + channel + \
                                    '_round_' + str(round_num) + '_fov_' +str(fov) + '-' + tokenize()
                    counts = delayed(fov_processing.single_fov_round_processing_eel,name=dask_delayed_name)(fov_subdataset,
                                                analysis_parameters,
                                                running_functions,
                                                dark_img,
                                                experiment_fpath,
                                                save_steps_output=save_intermediate_steps)
                    all_counts_fov.append(counts)
            
            name = 'concat_fish_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            all_counts_fov = delayed(pd.concat,name=name)(all_counts_fov,axis=0,ignore_index=True)
            
            name = 'create_nuclei_stack' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            filtered_nuclei_stack = delayed(utils.combine_filtered_images,name=name)(all_nuclei_fov,experiment_fpath,metadata)

            name = 'register_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            registered_counts = delayed(fovs_registration.nuclei_based_registration,name=name)(all_counts_fov,
                                                filtered_nuclei_stack,
                                                analysis_parameters)
                                                                                                
            name = 'stitch_to_mic_coords_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()  
            stitched_coords = delayed(stitching.stitch_using_microscope_fov_coords_new,name=name)(registered_counts,tile_corners_coords_pxl)
            
            name = 'save_df_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize() 
            saved_file = delayed(stitched_coords.to_parquet,name=name)(Path(experiment_fpath) / 'results'/ (experiment_name + \
                            '_decoded_fov_' + str(fov) + '.parquet'),index=False)
        
            all_processing.append(saved_file) 
        _ = dask.compute(*all_processing)

        # # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        stitching_selected = 'microscope_stitched'
        io.simple_output_plotting_serial(experiment_fpath, stitching_selected, client)
        # ----------------------------------------------------------------  
  
