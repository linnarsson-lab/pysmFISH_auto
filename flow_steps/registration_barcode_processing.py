from pathlib import Path
import numpy as np

from pysmFISH.fovs_registration import calculate_shift_hybridization_fov
from pysmFISH.fovs_registration import register_fish
from pysmFISH.barcodes_analysis import extract_barcodes_NN
from pysmFISH.barcodes_analysis import extract_dots_images
from pysmFISH.barcodes_analysis import define_flip_direction

from pysmFISH.utils import combine_rounds_images



def registration_barcode_detection_basic(processing_grps,
                                        analysis_parameters,
                                        experiment_info,
                                        experiment_fpath,
                                        codebook,
                                        selected_genes,
                                        correct_hamming_distance):

    # Register fov using the registration channel
    registered_counts_df, all_rounds_shifts, file_tags, status = calculate_shift_hybridization_fov(processing_grps[0],analysis_parameters,
                                                                save=True)
    

    # Register the fish dots coords
    registered_fish_df, file_tags, status = register_fish(processing_grps[1],analysis_parameters,registered_counts_df,
                                            all_rounds_shifts,file_tags,status,save=True)

    # Extract and decode the barcodes
    process_barcodes = extract_barcodes_NN(registered_fish_df,
                                                analysis_parameters,
                                                experiment_info,
                                                codebook,
                                                file_tags,
                                                status)
    process_barcodes.run_extraction()


    # Create combined images
    # Remember to consider the status input in order to create dictionary
    if process_barcodes.status == 'SUCCESS':
        channels = np.unique(process_barcodes.barcoded_fov_df['dot_channel'])
        for channel in channels:
            grp = [el for el in processing_grps[1] if channel in el.stem]
            img_stack = combine_rounds_images(grp,experiment_fpath, 
                        experiment_info,all_rounds_shifts, save=True)

            # Isolate the dots_subimages
            all_regions = extract_dots_images(process_barcodes.barcoded_fov_df,
                                            img_stack,experiment_fpath,save=True)

            # Define flip position and direction
            define_flip_direction(experiment_fpath,output_df, selected_genes, correct_hamming_distance,save=True)
