
import dask

from pysmFISH.configuration_files import load_experiment_config_file
from pysmFISH.configuration_files import load_processing_env_config_file

from pysmFISH.processing_cluster_setup import start_processing_env

from pysmFISH.io import create_empty_zarr_file

from pysmFISH.microscopy_file_parsers import nd2_raw_files_selector_general
from pysmFISH.microscopy_file_parsers import nikon_nd2_reparser_zarr

experiment_fpath = '/wsfish/smfish_ssd/LBEXP20201207_EEL_HE_test2'
processing_env_config_fpath = '/wsfish/smfish_ssd/config_db'
raw_files_fpath = experiment_fpath + '/raw_data'
parsed_image_tag = 'img_data'


processing_env_config = load_processing_env_config_file(processing_env_config_fpath)
experiment_info = load_experiment_config_file(experiment_fpath)

# Cluster setup
cluster = start_processing_env(processing_env_config=processing_env_config, 
                        experiment_info=experiment_info,
                        experiment_fpath=experiment_fpath)


# Create empty zarr file for the parse data
parsed_raw_data_fpath = create_empty_zarr_file(experiment_fpath=experiment_fpath,
                                    tag=parsed_image_tag)

# Reparse the data
all_raw_nd2 = nd2_raw_files_selector_general(folder_fpath=raw_files_fpath)

all_futures = []
for nd2_file in all_raw_nd2:
    future = dask.delayed(nikon_nd2_reparser_zarr)(nd2_file,
                                            parsed_raw_data_fpath=parsed_raw_data_fpath,
                                            experiment_info=experiment_info)
    all_futures.append(future)

dask.compute(*all_futures)


print('cane')

cluster.close()