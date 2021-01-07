from flow_steps.filtering_counting import filtering_counting_runner
from flow_steps.filtering_counting import single_fish_filter_count_standard
from flow_steps.create_processing_cluster import create_processing_cluster

from pysmFISH.configuration_files import load_experiment_config_file
from pysmFISH.configuration_files import load_processing_env_config_file
from pysmFISH.io import consolidate_zarr_metadata
from pysmFISH.utils import sorting_grps


from pathlib import Path

experiment_fpath = '/wsfish/smfish_ssd/LBEXP20201207_EEL_HE_test2'
processing_env_config_fpath = '/wsfish/smfish_ssd/config_db'
parsed_raw_data_fpath = '/wsfish/smfish_ssd/LBEXP20201207_EEL_HE_test2/LBEXP20201207_EEL_HE_test2_img_data.zarr'


processing_env_config = load_processing_env_config_file(processing_env_config_fpath)
experiment_info = load_experiment_config_file(experiment_fpath)

cluster = create_processing_cluster(processing_env_config_fpath,experiment_fpath)

consolidated_grp = consolidate_zarr_metadata(parsed_raw_data_fpath)

sorted_grps = sorting_grps(consolidated_grp, experiment_info, analysis_parameters)


grp_name = 'fish'

sorted_images_list = list(Path(parsed_raw_data_fpath).glob('*Europium*'))
sorted_images_list = [el.stem for el in sorted_images_list]

# sorted_images_list = ['LBEXP20201207_EEL_HE_test2_Hybridization08_Cy5_fov_88',
#                     'LBEXP20201207_EEL_HE_test2_Hybridization08_Cy5_fov_10',
#                     'LBEXP20201207_EEL_HE_test2_Hybridization08_Cy5_fov_21'
# ]


processing_parameters =  {}
processing_parameters['PreprocessingFishFlatFieldKernel'] = (100,100) 
processing_parameters['PreprocessingFishFilteringSmallKernel'] = (8,8)
processing_parameters['PreprocessingFishFilteringLaplacianKernel'] = (0.5,0.5)
processing_parameters['CountingFishMinObjDistance'] = 2
processing_parameters['CountingFishMinObjSize'] = 2
processing_parameters['CountingFishMaxObjSize'] = 200
processing_parameters['CountingFishNumPeaksPerLabel'] = 1

filtering_counting_runner(cluster,
                            single_fish_filter_count_standard,
                            parsed_raw_data_fpath,
                            grp_name,
                            sorted_images_list,
                            processing_parameters)


cluster.close()