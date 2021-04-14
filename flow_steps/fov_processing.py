import pickle
import pandas as pd
import numpy as np
from pathlib import Path

import gc

import pysmFISH
import flow_steps

from pysmFISH.utils import register_combined_rounds_images
from pysmFISH.barcodes_analysis import extract_dots_images
from pysmFISH.barcodes_analysis import define_flip_direction

from pysmFISH.logger_utils import selected_logger

from pysmFISH.preprocessing import standard_not_norm_preprocessing
from pysmFISH.preprocessing import filter_remove_large_objs
from pysmFISH.preprocessing import nuclei_registration_filtering

from pysmFISH.dots_calling import osmFISH_peak_based_detection_test
from pysmFISH.dots_calling import osmFISH_barcoded_peak_based_detection_masked_thr_test
from pysmFISH.dots_calling import osmFISH_peak_based_detection_fast

from pysmFISH.fovs_registration import calculate_shift_hybridization_fov_test
from pysmFISH.fovs_registration import register_fish_test
from pysmFISH.fovs_registration import calculate_shift_hybridization_fov_nuclei
from pysmFISH.fovs_registration import register_fish_on_nuclei

from pysmFISH.barcodes_analysis import extract_barcodes_NN_test

from flow_steps.filtering_counting import single_fish_filter_count_standard_not_norm
from flow_steps.filtering_counting import single_fish_filter_count_standard_not_norm_test
from flow_steps.filtering_counting import single_fish_filter_count_avoid_large_obj_test
from flow_steps.filtering_counting import filtering_counting_both_beads
from flow_steps.filtering_counting import filtering_counting_both_beads_test
from flow_steps.filtering_counting import filtering_counting_large_beads_test

from pysmFISH.stitching import stitch_using_microscope_fov_coords_test



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
        if running_functions['fish_channels_preprocessing'] == 'filter_remove_large_objs':

            img, masked_img = getattr(pysmFISH.preprocessing,running_functions['fish_channels_preprocessing'])(
                                                            zarr_grp_name,
                                                            parsed_raw_data_fpath,
                                                            processing_parameters,
                                                            dark_img)

            counts = getattr(pysmFISH.dots_calling,running_functions['fish_channels_dots_calling'])(
                                                                masked_img,
                                                                fov_subdataset,
                                                                processing_parameters)                                              

        else:

            img = getattr(pysmFISH.preprocessing,running_functions['fish_channels_preprocessing'])(
                                                            zarr_grp_name,
                                                            parsed_raw_data_fpath,
                                                            processing_parameters,
                                                            dark_img)


            counts = getattr(pysmFISH.dots_calling,running_functions['fish_channels_dots_calling'])(
                                                                            img,
                                                                            fov_subdataset,
                                                                            processing_parameters)


    elif processing_type == 'staining':
            pass
    

    # process all type of registration
    else:
        processing_parameters = analysis_parameters[processing_type]
        
        img = getattr(pysmFISH.preprocessing,running_functions['reference_channels_preprocessing'])(
                                                                        zarr_grp_name,
                                                                        parsed_raw_data_fpath,
                                                                        processing_parameters,
                                                                        dark_img)



        counts = getattr(pysmFISH.dots_calling,running_functions['reference_channels_dots_calling'])(
                                                                            img,
                                                                            fov_subdataset,
                                                                            processing_parameters)

    
    if save_steps_output:
        fname = experiment_name + '_' + fov_subdataset.channel + '_round_' + str(fov_subdataset.round_num) + '_fov_' + str(fov_subdataset.fov_num)
        np.save(filtered_img_path / (fname + '.npy'),img )
        counts.to_parquet((fname + '.parquet'),index=False)

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

    img = getattr(pysmFISH.preprocessing,running_functions['reference_channels_preprocessing'])(
                                                                        zarr_grp_name,
                                                                        parsed_raw_data_fpath,
                                                                        processing_parameters,
                                                                        dark_img)

    if save_steps_output:
        fname = experiment_name + '_' + fov_subdataset.channel + '_round_' + str(fov_subdataset.round_num) + '_fov_' + str(fov_subdataset.fov_num)
        np.save(filtered_img_path / (fname + '.npy'),img )

    return (img,fov_subdataset)
