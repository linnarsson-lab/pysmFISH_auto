import prefect
from prefect import task, Flow, Parameter, flatten, unmapped
from prefect.engine.executors import DaskExecutor
from prefect.environments import LocalEnvironment
from prefect import Task
from prefect.environments.storage import Local

from pysmFISH.configuration_files_tasks import load_experiment_config_file


from prefect.utilities.debug import raise_on_exception


with Flow("parsing",environment=LocalEnvironment(DaskExecutor(address='tcp://193.10.16.58:12270')),
            storage=Local(directory='/home/simone/tmp_code/flows')) as flow:
   
    experiment_fpath = Parameter('experiment_fpath', default = '/wsfish/smfish_ssd/AMEXP20201110_EEL_HumanH1930001V1C_auto')
    parsed_raw_data_fpath = Parameter('parsed_raw_data_fpath',default='/wsfish/smfish_ssd/wsfish/smfish_ssd/AMEXP20201110_EEL_HumanH1930001V1C_auto/AMEXP20201110_EEL_HumanH1930001V1C_auto_img_data.zarr')
     
    # Load experiment configuration file generated by robofish machines
    load_exp_cfg = load_experiment_config_file()
    experiment_info = load_exp_cfg(experiment_fpath)
    
    
    # PORT 
    # Create the shoji database that will contain the data
    ref = create_shoji_db(experiment_info)
    
    #PORT
    # Get the list of raw image groups to preprocess
    analysis_parameters = load_analysis_parameters(experiment_name=experiment_info['EXP_name'],upstream_tasks=[ref])

    #PORT
    consolidated_zarr_grp = consolidate_zarr_metadata(parsed_raw_data_fpath,upstream_tasks=[analysis_parameters])
    consolidated_zarr_grp = consolidate_zarr_metadata(parsed_raw_data_fpath,upstream_tasks=[autoparser])        
    
    # PORT
    # Sort the type of images according to processing

    # Order of output from the sorting_grps:
    # fish_grp, fish_selected_parameters, beads_grp, beads_selected_parameters,\
    # staining_grp, staining_selected_parameters
    sorted_grps = sorting_grps(consolidated_zarr_grp,experiment_info,analysis_parameters,upstream_tasks=[consolidated_zarr_grp])



    # PORT
    filtered_fish_images_metadata = single_fish.map(zarr_grp_name=sorted_grps[0][0:10],
                parsed_raw_data_fpath=unmapped(parsed_raw_data_fpath),
                experiment_fpath=unmapped(experiment_fpath),
                FlatFieldKernel=unmapped(sorted_grps[1]['PreprocessingFishFlatFieldKernel']),
                FilteringSmallKernel=unmapped(sorted_grps[1]['PreprocessingFishFilteringSmallKernel']),
                LaplacianKernel=unmapped(sorted_grps[1]['PreprocessingFishFilteringLaplacianKernel']),
                min_distance=unmapped(sorted_grps[1]['CountingFishMinObjDistance']),
                min_obj_size=unmapped(sorted_grps[1]['CountingFishMinObjSize']),
                max_obj_size=unmapped(sorted_grps[1]['CountingFishMaxObjSize']),
                num_peaks_per_label=unmapped(sorted_grps[1]['CountingFishNumPeaksPerLabel']))

    # PORT
    save_images_metadata.map(filtered_fish_images_metadata)

    
    # Load experiment configuration file generated by robofish machines
    load_exp_cfg = load_experiment_config_file()
    experiment_info = load_exp_cfg(experiment_fpath)

    # Selected the raw .nd2 files
    nd2_selector = nd2_raw_files_selector()
    all_raw_nd2 = nd2_selector(experiment_fpath)

    # Create the zarr file that will contain the parsed data
    create_zarr = create_empty_zarr_file()
    parsed_raw_data_fpath = create_zarr(experiment_fpath,parsed_image_tag)


    # Parse all the nd2 files
    autoparser = nikon_nd2_autoparser_zarr(task_run_name=lambda **kwargs: f"parsing-{kwargs['nd2_file_path']}")
    parsed_data = autoparser.map(nd2_file_path=all_raw_nd2,parsed_raw_data_fpath=unmapped(parsed_raw_data_fpath),
                                    experiment_info=unmapped(experiment_info))
    parsed_data.set_upstream([parsed_raw_data_fpath,all_raw_nd2])

# with raise_on_exception():
#     flow.run()

flow.register(project_name="test")
# flow.run_agent()