import prefect
from prefect import task, Flow, Parameter, flatten, unmapped
from prefect.engine.executors import DaskExecutor
from prefect.utilities.debug import raise_on_exception
from datetime import timedelta, datetime
from prefect.schedules import IntervalSchedule



from pysmFISH.dask_cluster_utilities_tasks import start_processing_env, local_cluster_setup

from pysmFISH.configuration_files_tasks import load_processing_env_config_file, load_experiment_config_file
from pysmFISH.data_handling import create_shoji_db

from pysmFISH.microscopy_file_parsers_tasks import nd2_raw_files_selector, nikon_nd2_autoparser, nikon_nd2_autoparser_single_files, nikon_nd2_autoparser_zarr, nikon_nd2_autoparser_zarr_single_files
from pysmFISH.qc_tasks import check_matching_metadata_robofish
from pysmFISH.utilities_tasks import check_completed_transfer_to_monod, sort_data_folder, create_empty_zarr_file, load_analysis_parameters
from pysmFISH.utilities_tasks import create_folder_structure, collect_extra_files,load_data_array,consolidate_zarr_metadata,sorting_grps,load_raw_images
from pysmFISH.notifications_tasks import report_input_files_errors
from pysmFISH.preprocessing_tasks import preprocessing_dot_raw_image, load_dark_image
from pysmFISH.logger_utils import setup_extra_loggers, prefect_logging_setup



# testing import
from pathlib import Path
import zarr

@task
def test_parallel():
    print('test_parallel')
    return 22

@task
def selection(parsed_raw_data_fpath):
    st = zarr.DirectoryStore(parsed_raw_data_fpath)
    consolidated_zarr_grp = zarr.consolidate_metadata(st)
    all_grps = list(consolidated_zarr_grp.keys())
    fish_grp = all_grps[0:3]
    beads_grp = all_grps[3:]
    return (fish_grp, beads_grp)


if __name__ == '__main__':

    # Add all the components to check for _auto files in the folder
    # or add it in the flow below

    # -------------------------------------------------------
    # Parameters for folder scanning for determine if an experiment
    # has been transferred to monod
    # -------------------------------------------------------
    logger = prefect_logging_setup("logger testing")
    
    flag_file_key = Parameter('flag_file_key', default='transfer_to_monod_completed.txt')
    # processing_hd_location = Parameter('processing_hd_location',default='/wsfish/smfish_ssd')
    processing_hd_location = Parameter('processing_hd_location',default='/Users/simone/Documents/local_data_storage/prefect_test/whd')

     # get info for submitting the error notification to github
    config_db_fpath = Path(processing_hd_location.default) / 'config_db'
    processing_env_config = load_processing_env_config_file(config_db_fpath)
    git_repo = processing_env_config['git_repo']
    git_token = processing_env_config['git_token']

    experiment_fpath = check_completed_transfer_to_monod(processing_hd_location.default,flag_file_key.default)
    experiment_info = load_experiment_config_file(experiment_fpath)

    # QC for experiment info if contains all the info


    # Activate processing cluster
    cluster = start_processing_env(processing_env_config,experiment_info)

    # schedule = IntervalSchedule(
    # start_date=datetime.utcnow() + timedelta(seconds=1),
    # interval=timedelta(minutes=1),)

    # with Flow("test_running",schedule=schedule) as flow:
    with Flow("test_running") as flow:

        # --------------------------------------------------
        #                 HOUSEKEEPING
        # --------------------------------------------------

        # Adjust folder structure and data
        # create_folder_structure(experiment_fpath)
        # collect_extra_files(experiment_fpath=experiment_fpath,experiment_info=experiment_info)
        # sort_data_folder(experiment_fpath,experiment_info)

        experiment_fpath = Parameter('experiment_fpath',default=experiment_fpath)
        experiment_info = Parameter('experiment_info',default=experiment_info)

        # Create the shoji database that will contain the data
        ref = create_shoji_db(experiment_info)
        # Get the list of raw image groups to preprocess
        #analysis_parameters = load_analysis_parameters(experiment_name=experiment_info['EXP_number'],upstream_tasks=[ref])
        
        PreprocessingFishFlatFieldKernel = Parameter('PreprocessingFishFlatFieldKernel',default=(100,100))
        FilteringSmallKernel = Parameter('FilteringSmallKernel',default=(8,8))
        PreprocessingFishFilteringLaplacianKernel = Parameter('PreprocessingFishFilteringLaplacianKernel',default=(1,1))

        PreprocessingBeadsRegistrationFlatFieldKernel = Parameter('PreprocessingBeadsRegistrationFlatFieldKernel',default=(100,100))
        PreprocessingBeadsRegistrationFilteringSmallKernel = Parameter('PreprocessingBeadsRegistrationFilteringSmallKernel',default=(8,8))
        PreprocessingBeadsRegistrationFilteringLaplacianKernel = Parameter('PreprocessingBeadsRegistrationFilteringLaplacianKernel',default=(1,1))


        # --------------------------------------------------
        #                     PARSING
        # --------------------------------------------------
        # Get all the .nd2 files to process
        # all_raw_files = nd2_raw_files_selector(experiment_fpath=experiment_fpath)
        
        # Run the crosscheck for all the pkl files
        # check_matching_metadata_robofish(all_raw_files)
        # report_input_files_errors(git_repo,experiment_fpath,git_token)
        # # Parse .nd2 files
        # tag = 'img_data'
        # parsed_raw_data_fpath = create_empty_zarr_file(experiment_fpath,tag)
        # nikon_nd2_autoparser_zarr.map(nd2_file_path=all_raw_files,parsed_raw_data_fpath=unmapped(parsed_raw_data_fpath),
        #                             experiment_info=unmapped(experiment_info))
        
        #parsed_raw_data_fpath = Parameter('parsed_raw_data_fpath',default='/wsfish/smfish_ssd/LBEXP20200708_EEL_Mouse_oPool5_auto/LBEXP20200708_EEL_Mouse_oPool5_auto_img_data.zarr')
        parsed_raw_data_fpath = Parameter('parsed_raw_data_fpath',default='/Users/simone/Documents/local_data_storage/prefect_test/whd/LBEXP20200708_EEL_Mouse_oPool5_auto/LBEXP20200708_EEL_Mouse_oPool5_auto_img_data.zarr')
        
        # consolidated_zarr_grp = consolidate_zarr_metadata(parsed_raw_data_fpath,upstream_tasks=[ref])        
        #sotre = zarr.DirectoryStore(parsed_raw_data_fpath.default)
        #consolidated_zarr_grp =zarr.open_consolidated(sotre)
        # Sort the type of images according to processing

        

        # fish_grp, fish_selected_parameters, beads_grp, beads_selected_parameters,\
        # staining_grp, staining_selected_parameters = sorting_grps(consolidated_zarr_grp,experiment_info,analysis_parameters)
        # # --------------------------------------------------
        out = selection(parsed_raw_data_fpath)
        # --------------------------------------------------
        #         PREPROCESSING AND DOTS CALLING                        
        # --------------------------------------------------
        dark_img = load_dark_image(experiment_fpath,upstream_tasks=[out])
        raw_fish_images_meta = load_raw_images.map(zarr_grp_name=out[0],
                                experiment_name=unmapped(experiment_info['EXP_number']),
                                parsed_raw_data_fpath=unmapped(parsed_raw_data_fpath))
        
        # filtered_fish_images = preprocessing_dot_raw_image.map(raw_fish_images_meta,
        #                     dark_img=unmapped(dark_img),
        #                     FlatFieldKernel=unmapped(fish_selected_parameters['PreprocessingFishFlatFieldKernel']),
        #                     FilteringSmallKernel=unmapped(fish_selected_parameters['FilteringSmallKernel']),
        #                     LaplacianKernel=unmapped(fish_selected_parameters['PreprocessingFishFilteringLaplacianKernel']))

        filtered_fish_images = preprocessing_dot_raw_image.map(raw_fish_images_meta,
                            dark_img=unmapped(dark_img),
                            FlatFieldKernel=unmapped(PreprocessingFishFlatFieldKernel),
                            FilteringSmallKernel=unmapped(FilteringSmallKernel),
                            LaplacianKernel=unmapped(PreprocessingFishFilteringLaplacianKernel))

        raw_beads_images_meta = load_raw_images.map(zarr_grp_name=out[1],
                                experiment_name=unmapped(experiment_info['EXP_number']),
                                parsed_raw_data_fpath=unmapped(parsed_raw_data_fpath))
        
        # filtered_beads_images = preprocessing_dot_raw_image.map(raw_beads_images_meta,
        #                     dark_img=unmapped(dark_img),
        #                     FlatFieldKernel=unmapped(beads_selected_parameters['PreprocessingBeadsRegistrationFlatFieldKernel']),
        #                     FilteringSmallKernel=unmapped(beads_selected_parameters['PreprocessingBeadsRegistrationFilteringSmallKernel']),
        #                     LaplacianKernel=unmapped(beads_selected_parameters['PreprocessingBeadsRegistrationFilteringLaplacianKernel']))

        filtered_beads_images = preprocessing_dot_raw_image.map(raw_beads_images_meta,
                            dark_img=unmapped(dark_img),
                            FlatFieldKernel=unmapped(PreprocessingBeadsRegistrationFlatFieldKernel),
                            FilteringSmallKernel=unmapped(PreprocessingBeadsRegistrationFilteringSmallKernel),
                            LaplacianKernel=unmapped(PreprocessingBeadsRegistrationFilteringLaplacianKernel))


        # experiment_fpath = Path('/Users/simone/Documents/local_data_storage/prefect_test/whd/exp_pre_auto')
        

        # LOAD DATA ARRAY OF RAW IMAGES
        # img_data_arrays = load_data_array.map(input_tuple=flatten(output_parsing))

        # FILTERING
        # NB when it works we can do different filtering depending on the type of tissue to run
        # use the case function in prefect

        # test = test_parallel(upstream_tasks=[check_matching_metadata_robofish(all_raw_files)])
        # create_raw_tmp_folders(experiment_fpath=experiment_fpath,upstream_tasks=[check_matching_metadata_robofish(all_raw_files)])

        # test2 = test_parallel(upstream_tasks=[test])    
        # test3 = test_parallel(upstream_tasks=[test])

        # map the parsing function on all the files

        # flow.register(project_name="project_test")
    
    executor = DaskExecutor(address=cluster.scheduler_address)
    # with raise_on_exception():


    flow_state = flow.run(executor=executor)
    # flow_state = flow.run(executor=executor)

    flow.visualize(flow_state=flow_state)
    cluster.close()