import time
import traceback
import pickle

import pandas as pd
from pathlib import Path
from dask.distributed import Client
from dask.distributed import as_completed, wait
from dask import dataframe as dd

from pysmFISH.logger_utils import json_logger
from pysmFISH.logger_utils import selected_logger

from pysmFISH.configuration_files import load_experiment_config_file
from pysmFISH.configuration_files import load_analysis_config_file
from pysmFISH.configuration_files import create_specific_analysis_config_file
from pysmFISH.configuration_files import create_function_runner

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
from flow_steps.filtering_counting import load_dark_image

from flow_steps.fov_processing import fov_processing_eel_barcoded
from flow_steps.fov_processing import fov_processing_eel_barcoded_dev

from pysmFISH.stitching import organize_square_tiles
from pysmFISH.stitching import stitch_using_microscope_fov_coords
# from pysmFISH.stitching import get_dots_in_overlapping_regions
# from pysmFISH.stitching import removed_duplicated_dots

from pysmFISH.qc_utils import QC_registration_error
from pysmFISH.qc_utils import check_experiment_yaml_file



# def general_flow(experiment_fpath:str, run_type:str='new', parsing_type:str='original'):

#     """
#     Flows for running human embryo eel experiment

#     Args:
#     -----
#         experiment_fpath: str
#             path to the folder with the experiment to process
#         run_type: str
#             type of flow run
#             - new
#             - re-run
#         parsing_type: str
#             key to select the type of data parsing to run
#             - original: parse the data out from robofish system
#             - reparsing_from_processing_folder: parse the raw data stored in the
#                                 experiment folder in the processing HD
#             - reparsing_from_storage: parse the raw data stored in the
#                                 experiment folder in the storage HD
#             - no_parsing: skip parsing step
#     """

pipeline_start = time.time()

# ----------------------------------------------------------------
# PARAMETERS DEFINITION
# Experiment fpath will be loaded from the scanning function

experiment_fpath = '/fish/work_std/LBEXP20210211_EEL_HE_3340um'

raw_data_folder_storage_path = '/fish/rawdata'
results_data_folder_storage_path = '/fish/results'
parsed_image_tag = 'img_data'

# run type can be:
# new
# re-run
run_type = 'new'

# parsing type (str) can be:
# original
# reparsing_from_processing_folder
# reparsing_from_storage 
# None if parsing not to be performed

parsing_type = 'original'

storage_experiment_fpath = (Path(raw_data_folder_storage_path) / Path(experiment_fpath).stem).as_posix()

# ----------------------------------------------------------------

if run_type == 'new': 
    # ----------------------------------------------------------------
    # CREATE FOLDERS STRUCTURE
    create_folder_structure(experiment_fpath)
    # ----------------------------------------------------------------

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
experiment_info = load_experiment_config_file(experiment_fpath)
# Add check if an analysis file is already present
create_specific_analysis_config_file(experiment_fpath, experiment_info)
analysis_parameters = load_analysis_config_file(experiment_fpath)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# CREATE RUNNING FUNCTIONS
# used to select the sets of functions to run preprocessing and dots
# calling according to the type of processing and sample
running_functions = create_function_runner(experiment_fpath,experiment_info)
# ----------------------------------------------------------------

if run_type == 'new': 
    # ----------------------------------------------------------------
    # ORGANIZE DATA IN FOLDERS
    collect_processing_files(experiment_fpath, experiment_info)
    sort_data_into_folders(experiment_fpath, experiment_info)
    # ----------------------------------------------------------------

# ----------------------------------------------------------------
# START PIPELINE LOGGER
logger = json_logger((experiment_fpath + '/logs'),'pipeline_run')
logger_print = selected_logger()
# logger = selected_logger()
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# CREATE DARK IMAGES
create_dark_img(experiment_fpath,experiment_info)
# ----------------------------------------------------------------


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

if parsing_type == 'no_parsing':
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

        # wait(parsing_futures)
        # _ = client.gather(parsing_futures)
    
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

        # _ = client.gather(parsing_futures)
    wait(parsing_futures)
    consolidated_grp = consolidate_zarr_metadata(parsed_raw_data_fpath)
    del parsing_futures

logger.info(f'reparsing completed in {(time.time()-start)/60} min')
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# IMAGE PREPROCESSING, DOTS COUNTING, REGISTRATION TO MICROSCOPE COORDS
start = time.time()
logger.info(f'start preprocessing and dots counting')

codebook = pd.read_parquet(Path(experiment_fpath) / 'codebook' / experiment_info['Codebook'])
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

tiles_org = organize_square_tiles(experiment_fpath,experiment_info,
                                    consolidated_grp,
                                    registration_reference_hybridization)
tiles_org.run_tiles_organization()
tile_corners_coords_pxl = tiles_org.tile_corners_coords_pxl

dark_img = load_dark_image(experiment_fpath)

# Scattering will be beneficial but causes error on HTCondor
# scatter the data to different workers to save timr
# remote_tile_corners_coords_pxl = client.scatter(tile_corners_coords_pxl)
# remote_codebook = client.scatter(codebook)
# remote_dark_img = client.scatter(dark_img)

logger_print.info(f'check if the logger is printing')

all_futures = []
start = time.time()

fname = Path(experiment_fpath) / 'tmp' / 'sorted_groups.pkl'
pickle.dump(sorted_grps, open(fname,'wb'))

for fov,sorted_grp in sorted_grps.items():
    future = client.submit(fov_processing_eel_barcoded_dev,
                                        fov=fov,
                                        sorted_grp=sorted_grp,
                                        experiment_info=experiment_info,
                                        analysis_parameters=analysis_parameters,
                                        experiment_fpath=experiment_fpath,
                                        parsed_raw_data_fpath=parsed_raw_data_fpath,
                                        running_functions=running_functions,
                                        img_width=img_width,
                                        img_height=img_height,
                                        tile_corners_coords_pxl= tile_corners_coords_pxl,
                                        codebook=codebook,
                                        selected_genes=selected_genes,
                                        correct_hamming_distance=correct_hamming_distance,
                                        dark_img = dark_img,
                                        save_steps_output=False,
                                        key= ('processing-fov-'+str(fov)))
        

    all_futures.append(future)

_ = client.gather(all_futures)
# tracebacks = {}
# for future in as_completed(all_futures):
#     logger_print.info(f'processed {future.key} in {time.time()-start} sec')
#     tracebacks[future.key] = traceback.format_tb(future.traceback())
#     del future

# wait(all_futures)
# fname = Path(experiment_fpath) / 'tmp' / 'tracebacks_processing_decoding.pkl'
# pickle.dump(tracebacks, open(fname,'wb'))
# logger.info(f'preprocessing and dots counting completed in {(time.time()-start)/60} min')

# del all_futures
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

# ----------------------------------------------------------------
# REMOVE DUPLICATED DOTS FROM THE OVERLAPPING REGIONS
# start = time.time()
# logger.info(f'plot registration error')

# unfolded_overlapping_regions_dict = {key:value for (k,v) in tiles_org.overlapping_regions.items() for (key,value) in v.items()}
# corrected_overlapping_regions_dict = {}
# for key, value in unfolded_overlapping_regions_dict.items():
#     corrected_overlapping_regions_dict[key] = np.array(value)-img_width

# # Prepare the dataframe
# select_genes = 'below3Hdistance_genes'
# stitching_selected = 'microscope_stitched'
# r_tag = 'r_px_' + stitching_selected
# c_tag = 'c_px_' + stitching_selected

# counts_dd = dd.read_parquet(experiment_fpath / 'tmp' / 'registered_counts' / '*decoded*.parquet')
# counts_dd = counts_dd.loc[counts_dd.dot_id == counts_dd.barcode_reference_dot_id,:]
# counts_df = counts_dd.dropna(subset=[select_genes]).compute()
# grpd = counts_df.groupby(select_genes)
# for gene, count_df in grpd:
    




logger.info(f'pipeline run completed in {(time.time()-pipeline_start)/60} min')

client.close()
cluster.close()