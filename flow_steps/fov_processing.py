import pickle
import pandas as pd
import numpy as np
from pathlib import Path


import pysmFISH
import flow_steps

from pysmFISH.utils import register_combined_rounds_images
from pysmFISH.barcodes_analysis import extract_dots_images
from pysmFISH.barcodes_analysis import define_flip_direction

from pysmFISH.logger_utils import selected_logger

from pysmFISH.fovs_registration import calculate_shift_hybridization_fov_test
from pysmFISH.fovs_registration import register_fish_test

from pysmFISH.barcodes_analysis import extract_barcodes_NN_test

from flow_steps.filtering_counting import single_fish_filter_count_standard_not_norm
from flow_steps.filtering_counting import single_fish_filter_count_standard_not_norm_test
from flow_steps.filtering_counting import filtering_counting_both_beads
from flow_steps.filtering_counting import filtering_counting_both_beads_test

from pysmFISH.stitching import stitch_using_microscope_fov_coords_test


def fov_processing_eel_barcoded(fov,
                    sorted_grp,
                    experiment_info,
                    analysis_parameters,
                    experiment_fpath,
                    parsed_raw_data_fpath,
                    running_functions,
                    img_width,
                    img_height,
                    tile_corners_coords_pxl,
                    codebook,
                    selected_genes,
                    correct_hamming_distance,
                    dark_img,
                    save_steps_output=False):
    
    
    processing_grp_split = sorted_grp['split']
    fish_img_stacks = {}

    logger = selected_logger()

    experiment_fpath = Path(experiment_fpath)
    experiment_name = experiment_fpath.stem
    
    # Path of directory where to save the intermediate results
    filtered_img_path = experiment_fpath / 'tmp' / 'filtered_images'
    raw_counts_path = experiment_fpath / 'tmp' / 'raw_counts'
    registered_counts_path = experiment_fpath / 'tmp' / 'registered_counts'
    combined_images_path = experiment_fpath / 'tmp' / 'combined_rounds_images'

    total_rounds = experiment_info['Barcode_length']
    
    counts_output = {}
    

    # # Filter and count
    # for processing_type, processing_input in processing_grp_split.items():
    #     if processing_type == 'fish':
            
    #         counts_output['fish'] = {}
    #         for channel, data_info in processing_input.items():
    #             img_stack = np.zeros([total_rounds,img_height,img_width])
    #             counts_output['fish'][channel] = {}
    #             names = data_info[0]
    #             processing_parameters = data_info[1]
    #             for zarr_grp_name in names:
    #                 round_num = int(zarr_grp_name.split('_')[-4].split('Hybridization')[-1])
    #                 counts_output['fish'][channel][zarr_grp_name], img = single_fish_filter_count_standard_not_norm_test(
    #                                                                 zarr_grp_name,
    #                                                                 parsed_raw_data_fpath,
    #                                                                 processing_parameters,
    #                                                                 dark_img)
    #                 img_stack[round_num-1,:,:] = img
    #                 if save_steps_output:
    #                      np.save(filtered_img_path / (zarr_grp_name + '.npy'),img )
                
    #             fish_img_stacks[channel] = img_stack
        
        
    #     elif processing_type == 'registration':
    #         counts_output['registration'] = {}
    #         for channel, data_info in processing_input.items():
    #             counts_output['registration'][channel] = {}
    #             names = data_info[0]
    #             processing_parameters = data_info[1]
    #             for zarr_grp_name in names:
    #                 counts_output['registration'][channel][zarr_grp_name], img = filtering_counting_both_beads_test(
    #                             zarr_grp_name,
    #                             parsed_raw_data_fpath,
    #                             processing_parameters,
    #                             dark_img)

    #                 if save_steps_output:
    #                      np.save(filtered_img_path / (zarr_grp_name + '.npy'),img)

    #     elif processing_type == 'staining':
    #         pass

    #     if save_steps_output:
    #         fname = raw_counts_path / (experiment_info['EXP_name'] + '_fov_' +str(fov) +'.pkl')
    #         pickle.dump(counts_output,open(fname,'wb'))
    
    
    # # Register the reference channel
    # registered_reference_channel_df, all_rounds_shifts, status = calculate_shift_hybridization_fov_test(
    #                                         fov,
    #                                         counts_output,
    #                                         analysis_parameters, 
    #                                         experiment_info)
    
    # if save_steps_output:
    #         fname = registered_counts_path / (experiment_info['EXP_name'] + '_fov_' +str(fov) +'.pkl')
    #         pickle.dump((registered_reference_channel_df, all_rounds_shifts),open(fname,'wb'))
    
    
    
    # # Register and decode the fish channels
    # output_channel = {}
    # for channel, counts in counts_output['fish'].items():
        
    #     output_channel[channel] = {}
    #     registered_fish_df, status = register_fish_test(
    #                         fov,
    #                         channel,
    #                         counts_output,
    #                         registered_reference_channel_df,
    #                         all_rounds_shifts,
    #                         analysis_parameters,
    #                         status)
    
    #     process_barcodes = extract_barcodes_NN_test(
    #                                 fov,
    #                                 channel,
    #                                 registered_fish_df,
    #                                 analysis_parameters,
    #                                 experiment_info,
    #                                 codebook,
    #                                 status)

    #     process_barcodes.run_extraction()

    #     registered_mic_df = stitch_using_microscope_fov_coords_test(process_barcodes.barcoded_fov_df,
    #                                                                 fov,
    #                                                                 tile_corners_coords_pxl)
         
    #     # Save the decoded data
    #     fname =  registered_counts_path / (experiment_name + '_' + channel + '_decoded_fov_' +str(fov) + '.parquet')
    #     registered_mic_df.to_parquet(fname,index=False) 
    #     # fname =  registered_counts_path / (experiment_name + '_' + channel + '_decoded_fov_' +str(fov) + '.pkl')
    #     # registered_mic_df.to_pickle(fname) 
 
    #     if process_barcodes.status == 'SUCCESS':
    #         registered_image = register_combined_rounds_images(fish_img_stacks[channel],all_rounds_shifts)
            
    #         all_regions = extract_dots_images(process_barcodes.barcoded_fov_df,
    #                                 registered_image,experiment_fpath,save=False)
            
    #         flipped_df = define_flip_direction(codebook,experiment_fpath,process_barcodes.barcoded_fov_df, 
    #                                     selected_genes, correct_hamming_distance,save=False)

    #         fname = combined_images_path / (experiment_name + '_' + channel + '_max_array_dict_' +str(fov) + '.npy')
    #         pickle.dump(all_regions, open(fname,'wb' ))
    #         fname = combined_images_path / (experiment_name + '_' + channel + '_flip_direction_' +str(fov) + '.parquet')
    #         flipped_df.to_parquet(fname,index=False)
    #         # fname = combined_images_path / (experiment_name + '_' + channel + '_flip_direction_' +str(fov) + '.pkl')
    #         # flipped_df.to_pickle(fname)

    #         if save_steps_output:
    #             fname = combined_images_path / (experiment_name + '_' + channel + '_combined_img_fov_' +str(fov) + '.npy')
    #             np.save(fname,registered_image)