import time
import traceback
import pickle
import zarr
import dask

import numpy as np
import pandas as pd
from pathlib import Path
from itertools import groupby
from dask.distributed import Client
from dask.distributed import as_completed, wait
from dask import dataframe as dd
from dask.base import tokenize

from pysmFISH.logger_utils import json_logger
from pysmFISH.logger_utils import selected_logger

from pysmFISH.data_organization import reorganize_processing_dir
from pysmFISH.data_organization import transfer_files_from_storage

from pysmFISH.data_models import Dataset

from pysmFISH.configuration_files import load_experiment_config_file
from pysmFISH.configuration_files import load_analysis_config_file
from pysmFISH.configuration_files import create_specific_analysis_config_file
from pysmFISH.configuration_files import create_function_runner

from pysmFISH.utils import create_folder_structure
from pysmFISH.utils import collect_processing_files
from pysmFISH.utils import sort_data_into_folders
from pysmFISH.utils import create_dark_img
from pysmFISH.utils import combine_filtered_images


from pysmFISH.io import create_empty_zarr_file
from pysmFISH.io import consolidate_zarr_metadata
from pysmFISH.io import open_consolidated_metadata
from pysmFISH.io import simple_output_plotting
from pysmFISH.io import simple_output_plotting_serial

from pysmFISH.microscopy_file_parsers import nd2_raw_files_selector_general
from pysmFISH.microscopy_file_parsers import nd2_raw_files_selector
from pysmFISH.microscopy_file_parsers import nikon_nd2_autoparser_zarr
from pysmFISH.microscopy_file_parsers import nikon_nd2_reparser_zarr
from pysmFISH.microscopy_file_parsers import single_nikon_nd2_parser_simple


from pysmFISH.utils import sorting_grps
from pysmFISH.utils import not_run_counting_sorted_grps
from pysmFISH.utils import sorting_grps_for_fov_processing

from pysmFISH.fovs_registration import beads_based_registration
from pysmFISH.fovs_registration import nuclei_based_registration
from pysmFISH.barcodes_analysis import decoder_fun
from pysmFISH.fovs_registration import create_registration_grps

from flow_steps.create_processing_cluster import create_processing_cluster
from flow_steps.filtering_counting import load_dark_image


from flow_steps.fov_processing import single_fov_round_processing_eel
from flow_steps.fov_processing import single_fov_round_processing_serial_nuclei

from pysmFISH.stitching import organize_square_tiles
from pysmFISH.stitching import stitch_using_microscope_fov_coords
from pysmFISH.stitching import remove_overlapping_dots_fov
from pysmFISH.stitching import clean_from_duplicated_dots
from pysmFISH.stitching import stitch_using_microscope_fov_coords_new

from pysmFISH.barcodes_analysis import extract_barcodes_NN_fast
from pysmFISH.preprocessing import fresh_nuclei_filtering

from pysmFISH.qc_utils import QC_registration_error
from pysmFISH.qc_utils import check_experiment_yaml_file


pipeline_start = time.time()

# ----------------------------------------------------------------
# PARAMETERS DEFINITION
# Experiment fpath will be loaded from the scanning function

experiment_fpath = '/datb/sl/camiel/Simone/CMEXP20210311'

raw_data_folder_storage_path = '/fish/rawdata'
dataset_folder_storage_path = '/fish/fish_datasets'
parsed_image_tag = 'img_data'

# run type can be:
# new
# re-run
run_type = 're-run'

# parsing type (str) can be:
# original
# reparsing_from_processing_folder
# reparsing_from_storage 
# no_parsing if parsing not to be performed

parsing_type = 'no_parsing'


fresh_nuclei_processing = False

save_intermediate_steps = True


storage_experiment_fpath = (Path(raw_data_folder_storage_path) / Path(experiment_fpath).stem).as_posix()

# ----------------------------------------------------------------


# ----------------------------------------------------------------
# CREATE FOLDERS STRUCTURE
create_folder_structure(experiment_fpath,run_type)
# ----------------------------------------------------------------


# ----------------------------------------------------------------
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
        _ = client.gather(parsing_futures)
    
    else:
        # add error if not correct parsing type
        if parsing_type == 'reparsing_from_processing_folder':
            raw_files_fpath = experiment_fpath + '/raw_data'
            logger.info(f'raw_files_fpath {raw_files_fpath}')
        elif parsing_type == 'reparsing_from_storage':
            raw_files_fpath = storage_experiment_fpath + '/raw_data'
        
        all_raw_nd2 = nd2_raw_files_selector_general(folder_fpath=raw_files_fpath)
        parsing_futures = client.map(nikon_nd2_reparser_zarr,
                                all_raw_nd2,
                                parsed_raw_data_fpath=parsed_raw_data_fpath,
                                experiment_info=experiment_info)

    _ = client.gather(parsing_futures)
    # wait(parsing_futures)
    consolidated_grp = consolidate_zarr_metadata(parsed_raw_data_fpath)
    # del parsing_futures

logger.info(f'reparsing completed in {(time.time()-start)/60} min')
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# CREATE DATASET
start = time.time()
logger.info(f'start dataset creation')
ds = Dataset()
ds.create_full_dataset_from_zmetadata(parsed_raw_data_fpath)

# ds.load_dataset('/datb/sl/camiel/Simone/CMEXP20210311/210503_06_58_44_CMEXP20210311_img_data_dataset_corrected.parquet')

metadata = ds.collect_metadata(ds.dataset)
# # ds.dataset.loc[:,'stitching_channel'] = 'Europium'
# # ds.dataset.loc[ds.dataset.channel == 'Europium','processing_type'] = 'large-beads'

logger.info(f'dataset creation completed in {(time.time()-start)/60} min')



# ----------------------------------------------------------------
# CREATE RUNNING FUNCTIONS
# used to select the sets of functions to run preprocessing and dots
# calling according to the type of processing and sample
running_functions = create_function_runner(experiment_fpath,metadata)
# ----------------------------------------------------------------


# # ----------------------------------------------------------------
# DETERMINE TILES ORGANIZATION
start = time.time()
logger.info(f'start calculation of tiles organization')
reference_round = analysis_parameters['RegistrationReferenceHybridization']
tiles_org = organize_square_tiles(experiment_fpath,metadata,
                                reference_round)
tiles_org.run_tiles_organization()
tile_corners_coords_pxl = tiles_org.tile_corners_coords_pxl

logger.info(f'calculation of tiles organization completed in {(time.time()-start)/60} min')

dark_img = load_dark_image(experiment_fpath)
dark_img = dask.delayed(dark_img)
analysis_parameters = dask.delayed(analysis_parameters)
running_functions = dask.delayed(running_functions)
tile_corners_coords_pxl = dask.delayed(tile_corners_coords_pxl)


grpd_fovs = ds.dataset.groupby('fov_num')

# ----------------------------------------------------------------
# PROCESSING SERIAL smFISH
# ----------------------------------------------------------------

start = time.time()
logger.info(f'start preprocessing and dots counting')

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

            out_nuclei = dask.delayed(single_fov_round_processing_serial_nuclei)(fov_subdataset,
                                    analysis_parameters,
                                    running_functions,
                                    dark_img,
                                    experiment_fpath,
                                    save_steps_output=save_intermediate_steps,
                                                dask_key_name = dask_delayed_name )
            all_nuclei_fov.append(out_nuclei)

        else:
            dask_delayed_name = 'filt_count_' +experiment_name + '_' + channel + \
                            '_round_' + str(round_num) + '_fov_' +str(fov) + '-' + tokenize()
            counts = dask.delayed(single_fov_round_processing_eel)(fov_subdataset,
                                        analysis_parameters,
                                        running_functions,
                                        dark_img,
                                        experiment_fpath,
                                        save_steps_output=save_intermediate_steps,
                                                    dask_key_name = dask_delayed_name )
            all_counts_fov.append(counts)
    
    name = 'concat_fish_' +experiment_name + '_' + channel + '_' \
                        + '_fov_' +str(fov) + '-' + tokenize()
    all_counts_fov = dask.delayed(pd.concat)(all_counts_fov,axis=0,ignore_index=True)
    
    name = 'create_nuclei_stack' +experiment_name + '_' + channel + '_' \
                        + '_fov_' +str(fov) + '-' + tokenize()
    
    filtered_nuclei_stack = dask.delayed(combine_filtered_images)(all_nuclei_fov,experiment_fpath,metadata)


    name = 'register_' +experiment_name + '_' + channel + '_' \
                        + '_fov_' +str(fov) + '-' + tokenize()
    
    registered_counts = dask.delayed(nuclei_based_registration)(all_counts_fov,
                                        filtered_nuclei_stack,
                                        analysis_parameters)
                                                                                        
    
    name = 'stitch_to_mic_coords_' +experiment_name + '_' + channel + '_' \
                        + '_fov_' +str(fov) + '-' + tokenize()  
    stitched_coords = dask.delayed(stitch_using_microscope_fov_coords_new)(registered_counts,tile_corners_coords_pxl)
    
    name = 'save_file_' +experiment_name + '_' + channel + '_' \
                        + '_fov_' +str(fov) + '-' + tokenize() 
    saved_file = dask.delayed(stitched_coords.to_parquet)(Path(experiment_fpath) / 'results'/ (experiment_name + \
                    '_decoded_fov_' + str(fov) + '.parquet'),index=False)

    all_processing.append(saved_file) 

    
_ = dask.compute(*all_processing)

# # ----------------------------------------------------------------
# GENERATE OUTPUT FOR PLOTTING
stitching_selected = 'microscope_stitched'
simple_output_plotting_serial(experiment_fpath, stitching_selected, client)
# ----------------------------------------------------------------  


logger.info(f'preprocessing and dots counting completed in {(time.time()-start)/60} min')

logger.info(f'pipeline run completed in {(time.time()-pipeline_start)/60} min')

client.close()
cluster.close()