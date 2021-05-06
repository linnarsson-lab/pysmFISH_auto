"""
Collection of the steps use to build the analysis pipeline
"""

from pysmFISH.logger_utils import selected_logger

from pysmFISH.configuration_files import load_experiment_config_file
from pysmFISH.configuration_files import create_specific_analysis_config_file

from pysmFISH.io import create_empty_zarr_file
from pysmFISH.io import consolidate_zarr_metadata

from pysmFISH.utils import collect_processing_files
from pysmFISH.utils import sort_data_into_folders

from pysmFISH.data_organization import transfer_files_from_storage

from pysmFISH.microscopy_file_parsers import nd2_raw_files_selector, nd2_raw_files_selector_general
from pysmFISH.microscopy_file_parsers import nikon_nd2_autoparser_zarr, nikon_nd2_reparser_zarr


def nikon_nd2_parsing_step(experiment_fpath,parsed_image_tag,storage_experiment_fpath,client):
    """
    explain all the steps that are run by this function
    """
    experiment_info = load_experiment_config_file(experiment_fpath)
    create_specific_analysis_config_file(experiment_fpath, experiment_info)
    
    # Create empty zarr file for the parse data
    parsed_raw_data_fpath = create_empty_zarr_file(experiment_fpath=experiment_fpath,
                                        tag=parsed_image_tag)
    if parsing_type == 'original':
        collect_processing_files(experiment_fpath, experiment_info)
        sort_data_into_folders(experiment_fpath, experiment_info)
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
            transfer_files_from_storage(storage_experiment_fpath, experiment_fpath)
            raw_files_fpath = storage_experiment_fpath + '/raw_data'
        
        all_raw_nd2 = nd2_raw_files_selector_general(folder_fpath=raw_files_fpath)
        parsing_futures = client.map(nikon_nd2_reparser_zarr,
                                all_raw_nd2,
                                parsed_raw_data_fpath=parsed_raw_data_fpath,
                                experiment_info=experiment_info)

    _ = client.gather(parsing_futures)
    consolidated_grp = consolidate_zarr_metadata(parsed_raw_data_fpath)



