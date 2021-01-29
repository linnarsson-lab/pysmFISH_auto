import time

import pandas as pd
from pathlib import Path
from dask.distributed import Client


from pysmFISH.configuration_files import load_experiment_config_file
from pysmFISH.configuration_files import load_processing_env_config_file
from pysmFISH.configuration_files import load_analysis_config_file
from pysmFISH.configuration_files import create_specific_analysis_config_file

from pysmFISH.io import create_empty_zarr_file
from pysmFISH.io import consolidate_zarr_metadata
from pysmFISH.io import open_consolidated_metadata # to remove

from pysmFISH.microscopy_file_parsers import nd2_raw_files_selector_general
from pysmFISH.microscopy_file_parsers import nikon_nd2_reparser_zarr

from pysmFISH.utils import sorting_grps

from pysmFISH.fovs_registration import create_registration_grps

from flow_steps.create_processing_cluster import create_processing_cluster
from flow_steps.filtering_counting import single_fish_filter_count_standard
from flow_steps.filtering_counting import single_fish_filter_count_standard_not_norm
from flow_steps.registration_barcode_processing import registration_barcode_detection_basic

from pysmFISH.stitching import organize_square_tiles
from pysmFISH.stitching import stitch_using_microscope_fov_coords

pipeline_start = time.time()

# ----------------------------------------------------------------
# PARAMETERS DEFINITION
experiment_fpath = '/wsfish/smfish_ssd/LBEXP20201207_EEL_HE_test2'
processing_env_config_fpath = '/wsfish/smfish_ssd/LBEXP20201207_EEL_HE_test2/config_db'
raw_files_fpath = experiment_fpath + '/raw_data'
parsed_image_tag = 'img_data'
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# LOAD CONFIGURATION FILES
start = time.time()
print(f'start loading configuration files')
processing_env_config = load_processing_env_config_file(experiment_fpath)
experiment_info = load_experiment_config_file(experiment_fpath)

# Add check if an analysis file is already present
create_specific_analysis_config_file(experiment_fpath, experiment_info)
analysis_parameters = load_analysis_config_file(experiment_fpath)
print(f'config files loading completed in {(time.time()-start)/60} min')
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# CREATE PROCESSING CLUSTER
start = time.time()
print(f'start cluster creation')
cluster = create_processing_cluster(processing_env_config_fpath,experiment_fpath)
client = Client(cluster)
print(f'cluster creation completed in {(time.time()-start)/60} min')
# ----------------------------------------------------------------


# # ----------------------------------------------------------------
# # REPARSING THE MICROSCOPY DATA
# start = time.time()
# print(f'start reparsing raw data')
# # Create empty zarr file for the parse data
# parsed_raw_data_fpath = create_empty_zarr_file(experiment_fpath=experiment_fpath,
#                                     tag=parsed_image_tag)

# # Reparse the data
# all_raw_nd2 = nd2_raw_files_selector_general(folder_fpath=raw_files_fpath)

# parsing_futures = client.map(nikon_nd2_reparser_zarr,
#                             all_raw_nd2,
#                             parsed_raw_data_fpath=parsed_raw_data_fpath,
#                             experiment_info=experiment_info)

# _ = client.gather(parsing_futures)

# print(f'reparsing completed in {(time.time()-start)/60} min')
# # ----------------------------------------------------------------


# ----------------------------------------------------------------
# IMAGE PREPROCESSING AND DOTS COUNTING
# start = time.time()
# print(f'start preprocessing and dots counting')
# # consolidated_grp = consolidate_zarr_metadata(parsed_raw_data_fpath)
parsed_raw_data_fpath = '/wsfish/smfish_ssd/LBEXP20201207_EEL_HE_test2/LBEXP20201207_EEL_HE_test2_img_data.zarr'
consolidated_grp = open_consolidated_metadata(parsed_raw_data_fpath)
sorted_grps = sorting_grps(consolidated_grp, experiment_info, analysis_parameters)


# Staining has different processing fun
all_futures = []
for grp, grp_data in sorted_grps.items():
    if grp in ['fish','beads']:
        for el in grp_data[0]:
            future = client.submit(single_fish_filter_count_standard,
                            el,
                            parsed_raw_data_fpath = parsed_raw_data_fpath,
                            processing_parameters=sorted_grps['fish'][1])
            all_futures.append(future)

start = time.time()
_ = client.gather(all_futures)
print(f'preprocessing and dots counting completed in {(time.time()-start)/60} min')
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# REGISTRATION AND BARCODE PROCESSING
start = time.time()
print(f'start registration and barcode processing')
registration_channel = 'Europium' # must be corrected in the config file
key = Path(experiment_fpath).stem + '_Hybridization01_' + registration_channel + '_fov_0'
fovs = consolidated_grp[key].attrs['fields_of_view']
codebook = pd.read_parquet(Path(experiment_fpath) / 'codebook' / 'gene_HE_V5_extended_EELV2_codebook_16_6_5Alex647N_positive_bits.parquet')
all_grps = create_registration_grps(experiment_fpath,registration_channel, fovs)


all_futures = client.map(registration_barcode_detection_basic, all_grps,
                        analysis_parameters = analysis_parameters,
                        experiment_info = experiment_info,
                        experiment_fpath = experiment_fpath,
                        codebook = codebook)
_ = client.gather(all_futures)

print(f'registration and barcode processing completed in {(time.time()-start)/60} min')
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# STITCHING
start = time.time()
print(f'start stitching using microscope coords')
round_num =1
tiles_org = organize_square_tiles(experiment_fpath,experiment_info,consolidated_grp,round_num)
tiles_org.run_tiles_organization()

decoded_files = list((Path(experiment_fpath) / 'tmp' / 'registered_counts').glob('*_decoded_*'))

all_futures = client.map(stitch_using_microscope_fov_coords,decoded_files,
                        tile_corners_coords_pxl = tiles_org.tile_corners_coords_pxl)       

_ = client.gather(all_futures)  

print(f'stitching using microscope coords completed in {(time.time()-start)/60} min')
# ----------------------------------------------------------------

print(f'pipeline run completed in {(time.time()-pipeline_start)/60} min')

client.close()
cluster.close()