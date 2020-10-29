"""
Test processing flow using dask delayed
"""
import prefect
from prefect import task, Flow, Parameter, flatten, unmapped
from prefect.engine.executors import DaskExecutor
from prefect.utilities.debug import raise_on_exception
from datetime import timedelta, datetime
from prefect.schedules import IntervalSchedule



from pysmFISH.dask_cluster_utilities_tasks import start_processing_env, local_cluster_setup

from pysmFISH.configuration_files_tasks import load_processing_env_config_file, load_experiment_config_file
from pysmFISH.data_model import create_shoji_db

from pysmFISH.microscopy_file_parsers_tasks import nd2_raw_files_selector, nikon_nd2_autoparser, nikon_nd2_autoparser_single_files, nikon_nd2_autoparser_zarr, nikon_nd2_autoparser_zarr_single_files
from pysmFISH.qc_tasks import check_matching_metadata_robofish
from pysmFISH.utilities_tasks import check_completed_transfer_to_monod, sort_data_folder, create_empty_zarr_file
from pysmFISH.utilities_tasks import create_folder_structure, collect_extra_files,load_data_array,consolidate_zarr_metadata,sorting_grps,load_raw_images, sorting_grps_fov
from pysmFISH.dots_calling import osmFISH_peak_based_detection

from pysmFISH.io import load_analysis_parameters, save_images_metadata, save_dots_data

from pysmFISH.fovs_registration import hybridizations_registration_grps, calculate_shift_hybridization_fov


from pysmFISH.notifications_tasks import report_input_files_errors
from pysmFISH.preprocessing_tasks import preprocessing_dot_raw_image, load_dark_image,test_preprocessing_large_scale
from pysmFISH.logger_utils import setup_extra_loggers, prefect_logging_setup

# testing import
from pathlib import Path
import zarr
import dask

if __name__ == '__main__':

    flag_file_key = 'transfer_to_monod_completed.txt'
    processing_hd_location = '/wsfish/smfish_ssd'

    config_db_fpath = Path(processing_hd_location.default) / 'config_db'
    processing_env_config = load_processing_env_config_file(config_db_fpath)

    experiment_fpath = check_completed_transfer_to_monod(processing_hd_location,flag_file_key)
    experiment_info = load_experiment_config_file(experiment_fpath)

    cluster = start_processing_env(processing_env_config,experiment_info)

    ref = create_shoji_db(experiment_info)
    analysis_parameters = load_analysis_parameters(experiment_name=experiment_info['EXP_number'])

    parsed_raw_data_fpath = '/wsfish/smfish_ssd/LBEXP20200708_EEL_Mouse_oPool5_auto/LBEXP20200708_EEL_Mouse_oPool5_auto_img_data.zarr'
    consolidated_zarr_grp = consolidate_zarr_metadata(parsed_raw_data_fpath)
    sorted_grps = sorting_grps(consolidated_zarr_grp,experiment_info,analysis_parameters)

    
    raw_fish_images_meta = []
    for grp in sorted_grps[0]:
        raw_fish_images_meta.append(dask.delayed(load_raw_images)(grp,parsed_raw_data_fpath=parsed_raw_data_fpath))

    filtered_data = []
    for img_meta in raw_fish_images_meta:
        filtered_data.append(dask.delayed(test_preprocessing_large_scale)(img_meta,experiment_fpath=experiment_fpath,
                            FlatFieldKernel=sorted_grps[1]['PreprocessingFishFlatFieldKernel'],
                            FilteringSmallKernel=sorted_grps[1]['PreprocessingFishFilteringSmallKernel'],
                            LaplacianKernel=sorted_grps[1]['PreprocessingFishFilteringLaplacianKernel']))

    dask.compute(*filtered_data)
        
    