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

from pysmFISH.data_organization import transfer_data_to_storage
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


from pysmFISH.io import create_empty_zarr_file
from pysmFISH.io import consolidate_zarr_metadata
from pysmFISH.io import open_consolidated_metadata
from pysmFISH.io import simple_output_plotting

from pysmFISH.microscopy_file_parsers import nd2_raw_files_selector_general
from pysmFISH.microscopy_file_parsers import nd2_raw_files_selector
from pysmFISH.microscopy_file_parsers import nikon_nd2_autoparser_zarr
from pysmFISH.microscopy_file_parsers import nikon_nd2_reparser_zarr
from pysmFISH.microscopy_file_parsers import single_nikon_nd2_parser_simple


from pysmFISH.utils import sorting_grps
from pysmFISH.utils import not_run_counting_sorted_grps
from pysmFISH.utils import sorting_grps_for_fov_processing

from pysmFISH.fovs_registration import create_registration_grps

from flow_steps.create_processing_cluster import create_processing_cluster
from flow_steps.filtering_counting import load_dark_image

from flow_steps.fov_processing import fov_processing_eel_barcoded
from flow_steps.fov_processing import fov_processing_eel_barcoded_dev
from flow_steps.fov_processing import single_fov_round_processing_eel

from pysmFISH.stitching import organize_square_tiles
from pysmFISH.stitching import stitch_using_microscope_fov_coords
from pysmFISH.stitching import remove_overlapping_dots_fov
from pysmFISH.stitching import clean_from_duplicated_dots


from pysmFISH.preprocessing import fresh_nuclei_filtering

from pysmFISH.qc_utils import QC_registration_error
from pysmFISH.qc_utils import check_experiment_yaml_file

def cane(x):
    return x
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

experiment_fpath = '/wsfish/smfish_ssd/test_new_dataset/JJEXP20201123_hGBM_Amine_test'

raw_data_folder_storage_path = '/fish/rawdata'
results_data_folder_storage_path = '/fish/results'
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


storage_experiment_fpath = (Path(raw_data_folder_storage_path) / Path(experiment_fpath).stem).as_posix()

# ----------------------------------------------------------------

if (run_type == 'new') or (parsing_type == 'reparsing_from_storage'): 
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
# # CREATE DATASET
# start = time.time()
# logger.info(f'start dataset creation')
# ds = Dataset()
# ds.create_full_dataset_from_zmetadata(experiment_fpath, 
#              experiment_info,
#              parsed_raw_data_fpath)

# logger.info(f'dataset creation completed in {(time.time()-start)/60} min')

# ----------------------------------------------------------------
# IMAGE PREPROCESSING, DOTS COUNTING,

ds = Dataset()
ds.load_dataset('/wsfish/smfish_ssd/test_new_dataset/JJEXP20201123_hGBM_Amine_test/210405_14_41_49_JJEXP20201123_hGBM_Amine_test_dataset.parquet')

start = time.time()
logger.info(f'start preprocessing and dots counting')
dark_img = load_dark_image(experiment_fpath)
dark_img = dask.delayed(dark_img)
all_futures_filtering_counting = []
all_imgs_fov = ds.select_all_imgs_fov(ds.dataset,100)
for index_value, fov_subdataset in all_imgs_fov.iterrows():
    round_num = fov_subdataset.round_num
    channel = fov_subdataset.channel
    fov = fov_subdataset.fov_num
    experiment_name = fov_subdataset.experiment_name
    dask_delayed_name = experiment_name + '_' + channel + \
                    '_round_' + str(round_num) + '_fov_' +str(fov) + '-' + tokenize()
    future = dask.delayed(single_fov_round_processing_eel)(fov_subdataset,
                                    analysis_parameters,
                                    running_functions,
                                    dark_img,
                                    experiment_fpath,
                                    save_steps_output=False,
                                            dask_key_name = dask_delayed_name )
    
    # future = client.submit(single_fov_round_processing_eel,fov_subdataset,
    #                                 analysis_parameters,
    #                                 running_functions,
    #                                 dark_img,
    #                                 experiment_fpath,
    #                                 save_steps_output=False)
    
    all_futures_filtering_counting.append(future)


d = dask.delayed(cane)(all_futures_filtering_counting)
z = client.compute(d)
# _ = client.gather(all_futures_filtering_counting)

logger.info(f'preprocessing and dots counting completed in {(time.time()-start)/60} min')



# # ----------------------------------------------------------------
# # IMAGE PREPROCESSING, DOTS COUNTING, REGISTRATION TO MICROSCOPE COORDS
# start = time.time()
# logger.info(f'start preprocessing and dots counting')

# codebook = pd.read_parquet(Path(experiment_fpath) / 'codebook' / experiment_info['Codebook'])
# sorted_grps = sorting_grps_for_fov_processing(consolidated_grp, experiment_info, analysis_parameters)

# # PROCESSING PARAMETERS
# registration_channel = experiment_info['StitchingChannel']
# key = Path(experiment_fpath).stem + '_Hybridization01_' + registration_channel + '_fov_0'
# fovs = consolidated_grp[key].attrs['fields_of_view']
# img_width = consolidated_grp[key].attrs['img_width']
# img_height = consolidated_grp[key].attrs['img_height']
# registration_reference_hybridization = analysis_parameters['RegistrationReferenceHybridization']
# selected_genes = 'below3Hdistance_genes'
# correct_hamming_distance = 'zeroHdistance_genes' 

# tiles_org = organize_square_tiles(experiment_fpath,experiment_info,
#                                     consolidated_grp,
#                                     registration_reference_hybridization)
# tiles_org.run_tiles_organization()
# tile_corners_coords_pxl = tiles_org.tile_corners_coords_pxl

# dark_img = load_dark_image(experiment_fpath)

# # Scattering will be beneficial but causes error on HTCondor
# # scatter the data to different workers to save timr
# # remote_tile_corners_coords_pxl = client.scatter(tile_corners_coords_pxl)
# # remote_codebook = client.scatter(codebook)
# # remote_dark_img = client.scatter(dark_img)

# logger_print.info(f'check if the logger is printing')

# all_futures = []
# start = time.time()

# fname = Path(experiment_fpath) / 'tmp' / 'sorted_groups.pkl'
# pickle.dump(sorted_grps, open(fname,'wb'))

# for fov,sorted_grp in sorted_grps.items():
#     future = client.submit(fov_processing_eel_barcoded_dev,
#                                         fov=fov,
#                                         sorted_grp=sorted_grp,
#                                         experiment_info=experiment_info,
#                                         analysis_parameters=analysis_parameters,
#                                         experiment_fpath=experiment_fpath,
#                                         parsed_raw_data_fpath=parsed_raw_data_fpath,
#                                         running_functions=running_functions,
#                                         img_width=img_width,
#                                         img_height=img_height,
#                                         tile_corners_coords_pxl= tile_corners_coords_pxl,
#                                         codebook=codebook,
#                                         selected_genes=selected_genes,
#                                         correct_hamming_distance=correct_hamming_distance,
#                                         dark_img = dark_img,
#                                         save_steps_output=False,
#                                         key= ('processing-fov-'+str(fov)))
        

#     all_futures.append(future)
# # wait(all_futures)
# _ = client.gather(all_futures)
# # tracebacks = {}
# # for future in as_completed(all_futures):
# #     logger_print.info(f'processed {future.key} in {time.time()-start} sec')
# #     tracebacks[future.key] = traceback.format_tb(future.traceback())
# #     del future

# # wait(all_futures)
# # fname = Path(experiment_fpath) / 'tmp' / 'tracebacks_processing_decoding.pkl'
# # pickle.dump(tracebacks, open(fname,'wb'))
# logger.info(f'preprocessing and dots counting completed in {(time.time()-start)/60} min')

# # del all_futures
# ----------------------------------------------------------------

# ----------------------------------------------------------------
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
# # REMOVE DUPLICATED DOTS FROM THE OVERLAPPING REGIONS
# start = time.time()
# logger.info(f'start removal of duplicated dots')

# unfolded_overlapping_regions_dict = {key:value for (k,v) in tiles_org.overlapping_regions.items() for (key,value) in v.items()}
# corrected_overlapping_regions_dict = {}
# for key, value in unfolded_overlapping_regions_dict.items():
#     corrected_overlapping_regions_dict[key] = np.array(value)-img_width

# # Prepare the dataframe
# select_genes = 'below3Hdistance_genes'
# stitching_selected = 'microscope_stitched'
# same_dot_radius = 100
# r_tag = 'r_px_' + stitching_selected
# c_tag = 'c_px_' + stitching_selected


# all_futures = []

# for cpl,chunk_coords in corrected_overlapping_regions_dict.items():
#     future = client.submit(remove_overlapping_dots_fov,
#                             cpl = cpl,
#                             chunk_coords=chunk_coords,
#                             experiment_fpath=experiment_fpath,
#                             stitching_selected=stitching_selected,
#                             select_genes=select_genes,
#                             same_dot_radius = same_dot_radius)

#     all_futures.append(future)

# to_remove = client.gather(all_futures)  
# to_remove = [el for tg in to_remove for el in tg]
# removed_dot_dict = {}
# for k, g in groupby(to_remove, key=lambda x: int(x.split('_')[0])):
#     removed_dot_dict.update({k:list(g)})

# for fov in fovs:
#     if fov not in removed_dot_dict.keys():
#         removed_dot_dict.update({fov:[]})

# logger_print.info(f'{removed_dot_dict.keys()}')

# for fov,dots_id_to_remove in removed_dot_dict.items():
#     future = client.submit(clean_from_duplicated_dots,
#                             fov = fov,
#                             dots_id_to_remove=dots_id_to_remove,
#                             experiment_fpath=experiment_fpath)

#     all_futures.append(future)

# _ = client.gather(all_futures)


# logger.info(f'removal of duplicated dots completed in {(time.time()-start)/60} min')
# # # ----------------------------------------------------------------


# # # ----------------------------------------------------------------
# # GENERATE OUTPUT FOR PLOTTING
# simple_output_plotting(experiment_fpath, stitching_selected, select_genes, client)

# ----------------------------------------------------------------
# # PROCESS FRESH NUCLEI
# start = time.time()
# logger.info(f'start processing of the fresh nuclei')
# if fresh_nuclei_processing:
#     if parsing_type == 'reparsing_from_storage':
#         pass
#     else:
#         try:
#             nuclei_fpath = list((Path(experiment_fpath) / 'fresh_nuclei').glob('*.nd2'))[0]
#         except:
#             logger.error(f'missing images of the fresh nuclei')
#         else:
#             parsing_future = client.submit(single_nikon_nd2_parser_simple,nuclei_fpath)
#             _ = client.gather(parsing_future)

#             # create zarr file
#             filtered_fpath = nuclei_fpath.parent / (nuclei_fpath.stem + '_filtered.zarr')
#             # create_empty_zarr_file(nuclei_fpath.parent.as_posix(), tag='filtered')

#             # filtering all the fovs
#             zarr_fpath = nuclei_fpath.parent / (nuclei_fpath.stem + '.zarr')
#             parsed_store = zarr.DirectoryStore(zarr_fpath)
#             parsed_root = zarr.group(store=parsed_store,overwrite=False)
#             fovs = list(parsed_root.keys())

#             all_futures = []
#             for fov in fovs:
#                 parsing_future = client.submit(fresh_nuclei_filtering,
#                                 parsed_raw_data_fpath=zarr_fpath,
#                                 filtered_raw_data_fpath=filtered_fpath,
#                                 fov=fov,
#                                 processing_parameters=analysis_parameters)
#                 all_futures.append(parsing_future)
#             _ = client.gather(all_futures)

# logger.info(f'processing of the fresh nuclei completed in {(time.time()-start)/60} min')
# ----------------------------------------------------------------



# # ----------------------------------------------------------------
# # TRANSFER THE RAW DATA TO STORAGE FOLDER
# start = time.time()
# logger.info(f'start data transfer to storage folder')

# if parsing_type != 'reparsing_from_storage':
#     transfer_data_to_storage(experiment_fpath,raw_data_folder_storage_path)

# logger.info(f'data transfer to storage folder completed in {(time.time()-start)/60} min')
# # ----------------------------------------------------------------





logger.info(f'pipeline run completed in {(time.time()-pipeline_start)/60} min')

client.close()
cluster.close()