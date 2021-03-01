import time

import pandas as pd
from pathlib import Path
from dask.distributed import Client

from pysmFISH.logger_utils import json_logger
from pysmFISH.logger_utils import selected_logger

from pysmFISH.configuration_files import load_experiment_config_file
from pysmFISH.configuration_files import load_processing_env_config_file
from pysmFISH.configuration_files import load_analysis_config_file
from pysmFISH.configuration_files import create_specific_analysis_config_file

from pysmFISH.utils import create_folder_structure
from pysmFISH.utils import transfer_files_from_storage
from pysmFISH.utils import collect_processing_files
from pysmFISH.utils import sort_data_into_folders
from pysmFISH.utils import create_dark_img


from pysmFISH.io import create_empty_zarr_file
from pysmFISH.io import consolidate_zarr_metadata
from pysmFISH.io import open_consolidated_metadata # to remove

from pysmFISH.microscopy_file_parsers import nd2_raw_files_selector_general
from pysmFISH.microscopy_file_parsers import nd2_raw_files_selector
from pysmFISH.microscopy_file_parsers import nikon_nd2_autoparser_zarr
from pysmFISH.microscopy_file_parsers import nikon_nd2_reparser_zarr

from pysmFISH.utils import sorting_grps
from pysmFISH.utils import not_run_counting_sorted_grps
from pysmFISH.utils import sorting_grps_for_fov_processing

from pysmFISH.fovs_registration import create_registration_grps

from flow_steps.create_processing_cluster import create_processing_cluster
from flow_steps.filtering_counting import single_fish_filter_count_standard
from flow_steps.filtering_counting import single_fish_filter_count_standard_not_norm
from flow_steps.filtering_counting import filtering_counting_both_beads
from flow_steps.registration_barcode_processing import registration_barcode_detection_basic


from pysmFISH.fovs_registration import calculate_shift_hybridization_fov_test
from pysmFISH.fovs_registration import register_fish_test

from pysmFISH.barcodes_analysis import extract_barcodes_NN_test

from flow_steps.fov_processing import fov_processing_eel_barcoded

from pysmFISH.stitching import organize_square_tiles
from pysmFISH.stitching import stitch_using_microscope_fov_coords

from pysmFISH.qc_utils import QC_registration_error
from pysmFISH.qc_utils import check_experiment_yaml_file

pipeline_start = time.time()

# ----------------------------------------------------------------
# PARAMETERS DEFINITION
# Experiment fpath will be loaded from the scanning function

experiment_fpath = '/wsfish/smfish_ssd/LBEXP20210209_EEL_HE_3680um'

raw_data_folder_storage_path = '/fish/rawdata'
results_data_folder_storage_path = '/fish/results'
parsed_image_tag = 'img_data'

# parsing type (str) can be:
# original
# reparsing_from_processing_folder
# reparsing_from_storage 
# None if parsing not to be performed

parsing_type = None

running_functions ={
                    'fish_channels_filtering_counting':single_fish_filter_count_standard_not_norm,
                   'registration_channel_filtering_counting':filtering_counting_both_beads,
                   'registration_reference':calculate_shift_hybridization_fov_test,
                   'registration_fish': register_fish_test,
                   'barcode_extraction': extract_barcodes_NN_test}


storage_experiment_fpath = (Path(raw_data_folder_storage_path) / Path(experiment_fpath).stem).as_posix()

# ----------------------------------------------------------------

# # ----------------------------------------------------------------
# # CREATE FOLDERS STRUCTURE
# create_folder_structure(experiment_fpath)
# # ----------------------------------------------------------------

# # ----------------------------------------------------------------
# # QC Experiment info file
# check_experiment_yaml_file(experiment_fpath)
# # ----------------------------------------------------------------


# # ----------------------------------------------------------------
# TRANSFER REQUIRED FILES FOR THE PROCESSING IF THE ANALYSIS START
# FROM RAW DATA IN THE STORAGE HD
if parsing_type == 'reparsing_from_storage':
    transfer_files_from_storage(storage_experiment_fpath, experiment_fpath)
# # ----------------------------------------------------------------

# ----------------------------------------------------------------
# LOAD CONFIGURATION FILES
processing_env_config = load_processing_env_config_file(experiment_fpath)
experiment_info = load_experiment_config_file(experiment_fpath)
# Add check if an analysis file is already present
create_specific_analysis_config_file(experiment_fpath, experiment_info)
analysis_parameters = load_analysis_config_file(experiment_fpath)

codebook = pd.read_parquet(Path(experiment_fpath) / 'codebook' / experiment_info['Codebook'])
# ----------------------------------------------------------------

# # ----------------------------------------------------------------
# # ORGANIZE DATA IN FOLDERS
# collect_processing_files(experiment_fpath, experiment_info)
# sort_data_into_folders(experiment_fpath, experiment_info)
# # ----------------------------------------------------------------

# ----------------------------------------------------------------
# START PIPELINE LOGGER
logger = json_logger((experiment_fpath + '/logs'),'pipeline_run')
# logger = selected_logger()
# ----------------------------------------------------------------


# # ----------------------------------------------------------------
# # CREATE DARK IMAGES
# create_dark_img(experiment_fpath,experiment_info)
# # ----------------------------------------------------------------


# ----------------------------------------------------------------
# CREATE PROCESSING CLUSTER
start = time.time()
logger.info(f'start cluster creation')
cluster = create_processing_cluster(experiment_fpath)
client = Client(cluster)
logger.info(f'cluster creation completed in {(time.time()-start)/60} min')
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# PARSING THE MICROSCOPY DATA
start = time.time()
logger.info(f'start parsing raw data')

if parsing_type == None:
    logger.info(f'data parsing is not required')
    exp_path = Path(experiment_fpath)
    parsed_raw_data_fpath = exp_path / (exp_path.stem + '_' + parsed_image_tag + '.zarr') 
    consolidated_grp = open_consolidated_metadata(parsed_raw_data_fpath.as_posix())
else:
    # Create empty zarr file for the parse data
    parsed_raw_data_fpath = create_empty_zarr_file(experiment_fpath=experiment_fpath,
                                        tag=parsed_image_tag)
    if parsing_type == 'original':
        all_raw_nd2 = nd2_raw_files_selector(experiment_fpath)

        parsing_futures = client.map(nikon_nd2_autoparser_zarr,
                                all_raw_nd2,
                                parsed_raw_data_fpath=parsed_raw_data_fpath,
                                experiment_info=experiment_info)

        _ = client.gather(parsing_futures)

        consolidated_grp = consolidate_zarr_metadata(parsed_raw_data_fpath)
    
    else:
        if parsing_type == 'reparsing_from_processing_folder':
            raw_files_fpath = experiment_fpath + '/raw_data'
        elif parsing_type == 'reparsing_from_storage':
            raw_files_fpath = storage_experiment_fpath + '/raw_data'

        all_raw_nd2 = nd2_raw_files_selector_general(folder_fpath=raw_files_fpath)
        parsing_futures = client.map(nikon_nd2_reparser_zarr,
                                all_raw_nd2,
                                parsed_raw_data_fpath=parsed_raw_data_fpath,
                                experiment_info=experiment_info)

        _ = client.gather(parsing_futures)

        consolidated_grp = consolidate_zarr_metadata(parsed_raw_data_fpath)

logger.info(f'reparsing completed in {(time.time()-start)/60} min')
 # ----------------------------------------------------------------


# ----------------------------------------------------------------
# IMAGE PREPROCESSING AND DOTS COUNTING
start = time.time()
logger.info(f'start preprocessing and dots counting')
sorted_grps = sorting_grps_for_fov_processing(consolidated_grp, experiment_info, analysis_parameters)

# PROCESSING PARAMETERS
registration_channel = experiment_info['StitchingChannel']
key = Path(experiment_fpath).stem + '_Hybridization01_' + registration_channel + '_fov_0'
fovs = consolidated_grp[key].attrs['fields_of_view']
img_width = consolidated_grp[key].attrs['img_width']
img_height = consolidated_grp[key].attrs['img_height']
registration_reference_hybridization = analysis_parameters['RegistrationReferenceHybridization']
selected_genes = 'below3Hdistance_genes'
correct_hamming_distance = 'zeroHdistance_genes' 

all_futures = []
for fov,sorted_grp in sorted_grps.items():
    future = client.submit(fov_processing_eel_barcoded,
                                        fov,
                                        sorted_grp,
                                        experiment_info=experiment_info,
                                        analysis_parameters=analysis_parameters,
                                        experiment_fpath=experiment_fpath,
                                        parsed_raw_data_fpath=parsed_raw_data_fpath,
                                        running_functions=running_functions,
                                        img_width=img_width,
                                        img_height=img_height,
                                        codebook=codebook,
                                        selected_genes=selected_genes,
                                        correct_hamming_distance=correct_hamming_distance,
                                        save_steps_output=False)
    
    all_futures.append(future)

start = time.time()
_ = client.gather(all_futures)
logger.info(f'preprocessing and dots counting completed in {(time.time()-start)/60} min')
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# STITCHING
start = time.time()
logger.info(f'start stitching using microscope coords')
round_num = analysis_parameters['RegistrationReferenceHybridization']
tiles_org = organize_square_tiles(experiment_fpath,experiment_info,consolidated_grp,round_num)
tiles_org.run_tiles_organization()

decoded_files = list((Path(experiment_fpath) / 'tmp' / 'registered_counts').glob('*_decoded_*'))

all_futures = client.map(stitch_using_microscope_fov_coords,decoded_files,
                        tile_corners_coords_pxl = tiles_org.tile_corners_coords_pxl)       

_ = client.gather(all_futures)  

logger.info(f'stitching using microscope coords completed in {(time.time()-start)/60} min')
# ----------------------------------------------------------------

# # ----------------------------------------------------------------
# # QC REGISTRATION ERROR
# start = time.time()
# logger.info(f'plot registration error')

# registration_error = QC_registration_error(client, experiment_fpath, analysis_parameters, 
#                                             tiles_org.tile_corners_coords_pxl, 
#                                             tiles_org.img_width, tiles_org.img_height)

# registration_error.run_qc()

# logger.info(f'plotting of the registration error completed in {(time.time()-start)/60} min')
# # ----------------------------------------------------------------


logger.info(f'pipeline run completed in {(time.time()-pipeline_start)/60} min')

client.close()
cluster.close()