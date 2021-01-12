from pysmFISH.fovs_registration import calculate_shift_hybridization_fov
from pysmFISH.barcodes_analysis import extract_barcodes_NN

def registration_barcode_detection_basic(processing_grps,
                                        analysis_parameters,
                                        experiment_config,
                                        experiment_fpath):

    registered_counts_df = calculate_shift_hybridization_fov(processing_grps[0],analysis_parameters)
    # extracted_barcodes_df = extract_barcodes_NN(registered_counts_df,
    #                                             analysis_parameters,
    #                                             experiment_config)


                                                