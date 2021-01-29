from pathlib import Path
from pysmFISH.fovs_registration import calculate_shift_hybridization_fov
from pysmFISH.fovs_registration import register_fish
from pysmFISH.barcodes_analysis import extract_barcodes_NN
from pysmFISH.stitching import stitch_using_microscope_fov_coords

def registration_barcode_detection_basic(processing_grps,
                                        analysis_parameters,
                                        experiment_info,
                                        experiment_fpath,
                                        codebook):

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
                                                status,
                                                save=True)
    process_barcodes.run_extraction()


    # Create combined images
    # Remember to consider the status input in order to create dictionary



    # Isolate the dots_subimage

    # Stitched the data
    stitched_df = stitch_using_microscope_fov_coords(process_barcodes.barcoded_fov_df,experiment_info)
    
