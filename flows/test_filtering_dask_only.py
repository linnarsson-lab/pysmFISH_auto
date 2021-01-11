import time

from pathlib import Path
from dask.distributed import Client

from flow_steps.filtering_counting import filtering_counting_runner,single_fish_filter_count_standard
from flow_steps.filtering_counting import single_fish_filter_count_standard
from flow_steps.create_processing_cluster import create_processing_cluster

from pysmFISH.configuration_files import create_specific_analysis_config_file
from pysmFISH.configuration_files import load_experiment_config_file
from pysmFISH.configuration_files import load_processing_env_config_file
from pysmFISH.configuration_files import load_analysis_config_file
from pysmFISH.io import consolidate_zarr_metadata, open_consolidated_metadata
from pysmFISH.utils import sorting_grps

from pysmFISH.fov_registration import create_hybridizations_registration_grps
from pysmFISH.fovs_registration import calculate_shift_hybridization_fov


experiment_fpath = '/wsfish/smfish_ssd/LBEXP20201207_EEL_HE_test2'
processing_env_config_fpath = '/wsfish/smfish_ssd/config_db'
parsed_raw_data_fpath = '/wsfish/smfish_ssd/LBEXP20201207_EEL_HE_test2/LBEXP20201207_EEL_HE_test2_img_data.zarr'

processing_env_config = load_processing_env_config_file(processing_env_config_fpath)
experiment_info = load_experiment_config_file(experiment_fpath)

create_specific_analysis_config_file(experiment_fpath, experiment_info)
analysis_parameters = load_analysis_config_file(experiment_fpath)

#consolidated_grp = consolidate_zarr_metadata(parsed_raw_data_fpath)
consolidated_grp = open_consolidated_metadata(parsed_raw_data_fpath)
sorted_grps = sorting_grps(consolidated_grp, experiment_info, analysis_parameters)

start = time.time()
cluster = create_processing_cluster(processing_env_config_fpath,experiment_fpath)
client = Client(cluster)

print(f'cluster starting {time.time()-start}')


print(f'start filtering')
# start = time.time()
# all_futures = []
# # Filtering smFISH
# fish_futures = client.map(single_fish_filter_count_standard,
#                             sorted_grps['fish'][0],
#                             parsed_raw_data_fpath = parsed_raw_data_fpath,
#                             processing_parameters=sorted_grps['fish'][1])
# all_futures.append(fish_futures) 
# # _ = client.gather(all_futures)

# # Filtering beads
# beads_futures = client.map(single_fish_filter_count_standard,
#                             sorted_grps['beads'][0],
#                             parsed_raw_data_fpath = parsed_raw_data_fpath,
#                             processing_parameters=sorted_grps['fish'][1])
# all_futures.append(beads_futures) 
# # _ = client.gather(all_futures)

# all_futures = [ft for grp_ft in all_futures for ft in grp_ft]

try:

    start = time.time()
    # Staing has different processing fun
    all_futures = []
    for grp, grp_data in sorted_grps.items():
        if grp in ['fish','beads']:
            for el in grp_data[0]:
                future = client.submit(single_fish_filter_count_standard,
                                el,
                                parsed_raw_data_fpath = parsed_raw_data_fpath,
                                processing_parameters=sorted_grps['fish'][1])
                all_futures.append(future)

    print(f'future created {time.time()-start}')

    print(f'total number of futures to process {len(all_futures)}')

    start = time.time()
    _ = client.gather(all_futures)

    print(f'future gathered {time.time()-start}')

    # start = time.time()
    # print(f'starting registration')
    # # Registration fovs
    # # registration_channel = experiment_config['StitchingChannel']
    # registration_channel = 'Europium' # must be corrected in the config file
    # fovs = consolidated_grp[[list(consolidated_grp.keys()][2]]['fields_of_view']
    # all_grps = create_hybridizations_registration_grps(experiment_fpath,registration_channel, fovs)
    
    # all_futures = client.map(calculate_shift_hybridization_fov, all_grps,
    #                         analysis_parameters = analysis_parameters)

    # _ = client.gather(all_futures)

    # print(f'future for registration gathered {time.time()-start}')

    cluster.close()

except OSError:
    print(f' ERROR IN RECONNECTING TO THE CLUSTER')
    cluster.close()

print(f'processing completed')

# Registration
# registration_channel =
# barcode_length = 
# stitching_type = 




# grp_name = 'fish'

# sorted_images_list = list(Path(parsed_raw_data_fpath).glob('*Europium*'))
# sorted_images_list = [el.stem for el in sorted_images_list]

# # sorted_images_list = ['LBEXP20201207_EEL_HE_test2_Hybridization08_Cy5_fov_88',
# #                     'LBEXP20201207_EEL_HE_test2_Hybridization08_Cy5_fov_10',
# #                     'LBEXP20201207_EEL_HE_test2_Hybridization08_Cy5_fov_21'
# # ]


# processing_parameters =  {}
# processing_parameters['PreprocessingFishFlatFieldKernel'] = (100,100) 
# processing_parameters['PreprocessingFishFilteringSmallKernel'] = (8,8)
# processing_parameters['PreprocessingFishFilteringLaplacianKernel'] = (0.5,0.5)
# processing_parameters['CountingFishMinObjDistance'] = 2
# processing_parameters['CountingFishMinObjSize'] = 2
# processing_parameters['CountingFishMaxObjSize'] = 200
# processing_parameters['CountingFishNumPeaksPerLabel'] = 1

# filtering_counting_runner(cluster,
#                             single_fish_filter_count_standard,
#                             parsed_raw_data_fpath,
#                             grp_name,
#                             sorted_images_list,
#                             processing_parameters)


# cluster.close()