from pathlib import Path
from pysmFISH.fovs_registration import calculate_shift_hybridization_fov
from pysmFISH.fovs_registration import register_fish
from pysmFISH.barcodes_analysis import extract_barcodes_NN

def registration_barcode_detection_basic(processing_grps,
                                        analysis_parameters,
                                        experiment_info,
                                        experiment_fpath,
                                        codebook):

    # Register fov using the registration channel
    registered_counts_df, all_rounds_shifts = calculate_shift_hybridization_fov(processing_grps[0],analysis_parameters)
    
    # Register the fish dots coords
    registered_fish_df = register_fish(processing_grps[1],analysis_parameters,registered_counts_df,all_rounds_shifts)

    # Extract and decode the barcodes
    process_barcodes = extract_barcodes_NN(registered_fish_df,
                                                analysis_parameters,
                                                experiment_info,
                                                codebook)
    process_barcodes.run_extraction()

    channel = process_barcodes.barcoded_fov_df.loc[0,'dot_channel']
    experiment_name = process_barcodes.barcoded_fov_df.loc[0,'experiment_name']
    fov_num = process_barcodes.barcoded_fov_df.loc[0,'fov_num']
    fname = Path(experiment_fpath) / 'tmp' / 'registered_counts' / (experiment_name + '_' + channel + '_decoded_fov_' + str(fov_num) + '.parquet')
    process_barcodes.barcoded_fov_df.to_parquet(fname,index=False)