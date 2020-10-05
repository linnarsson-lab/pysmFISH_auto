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

from pysmFISH.notifications_tasks import report_input_files_errors
from pysmFISH.preprocessing_tasks import preprocessing_dot_raw_image, load_dark_image
from pysmFISH.logger_utils import setup_extra_loggers, prefect_logging_setup



# testing import
from pathlib import Path
import zarr

@task
def test_parallel(x):
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
    processing_hd_location = Parameter('processing_hd_location',default='/wsfish/smfish_ssd')
    # processing_hd_location = Parameter('processing_hd_location',default='/Users/simone/Documents/local_data_storage/prefect_test/whd')

     # get info for submitting the error notification to github
    config_db_fpath = Path(processing_hd_location.default) / 'config_db'
    processing_env_config = load_processing_env_config_file(config_db_fpath)
    # git_repo = processing_env_config['git_repo']
    # git_token = processing_env_config['git_token']

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
        analysis_parameters = load_analysis_parameters(experiment_name=experiment_info['EXP_number'],upstream_tasks=[ref])


        # SOMEWHERE COLLECT THE INFO OF NUMBER OF HYBRIDIZATION AND FOVS


        # --------------------------------------------------
        #                     PARSING
        # --------------------------------------------------
        # Get all the .nd2 files to process
        all_raw_files = nd2_raw_files_selector(experiment_fpath=experiment_fpath,upstream_tasks=[analysis_parameters])
        
        # Run the crosscheck for all the pkl files
        # check_matching_metadata_robofish(all_raw_files)
        # report_input_files_errors(git_repo,experiment_fpath,git_token)
        # # Parse .nd2 files
        tag = 'img_data'
        parsed_raw_data_fpath = create_empty_zarr_file(experiment_fpath,tag,upstream_tasks=[all_raw_files])
        autoparser = nikon_nd2_autoparser_zarr.map(nd2_file_path=all_raw_files,parsed_raw_data_fpath=unmapped(parsed_raw_data_fpath),
                                    experiment_info=unmapped(experiment_info),upstream_tasks=[parsed_raw_data_fpath])
        
       # parsed_raw_data_fpath = Parameter('parsed_raw_data_fpath',default='/wsfish/smfish_ssd/LBEXP20200708_EEL_Mouse_oPool5_auto/LBEXP20200708_EEL_Mouse_oPool5_auto_img_data.zarr')
       # parsed_raw_data_fpath = Parameter('parsed_raw_data_fpath',default='/Users/simone/Documents/local_data_storage/prefect_test/whd/LBEXP20200708_EEL_Mouse_oPool5_auto/LBEXP20200708_EEL_Mouse_oPool5_auto_img_data.zarr')
        

        consolidated_zarr_grp = consolidate_zarr_metadata(parsed_raw_data_fpath,upstream_tasks=[autoparser])        
        
    
        # Sort the type of images according to processing

        # Order of output from the sorting_grps:
        # fish_grp, fish_selected_parameters, beads_grp, beads_selected_parameters,\
        # staining_grp, staining_selected_parameters
        sorted_grps = sorting_grps(consolidated_zarr_grp,experiment_info,analysis_parameters,upstream_tasks=[consolidated_zarr_grp])
        # # --------------------------------------------------
    
        # --------------------------------------------------
        #         PREPROCESSING AND DOTS CALLING                        
        # --------------------------------------------------
        dark_img = load_dark_image(experiment_fpath,upstream_tasks=[sorted_grps[0]])
        raw_fish_images_meta = load_raw_images.map(zarr_grp_name=sorted_grps[0],
                                parsed_raw_data_fpath=unmapped(parsed_raw_data_fpath))
        

        filtered_fish_images_metadata = preprocessing_dot_raw_image.map(raw_fish_images_meta,
                            dark_img=unmapped(dark_img),
                            FlatFieldKernel=unmapped(sorted_grps[1]['PreprocessingFishFlatFieldKernel']),
                            FilteringSmallKernel=unmapped(sorted_grps[1]['PreprocessingFishFilteringSmallKernel']),
                            LaplacianKernel=unmapped(sorted_grps[1]['PreprocessingFishFilteringLaplacianKernel']))

        save_images_metadata.map(filtered_fish_images_metadata)
        

        fish_counts = osmFISH_peak_based_detection.map(filtered_fish_images_metadata,
                    min_distance=unmapped(sorted_grps[1]['CountingFishMinObjDistance']),
                    min_obj_size=unmapped(sorted_grps[1]['CountingFishMinObjSize']),
                    max_obj_size=unmapped(sorted_grps[1]['CountingFishMaxObjSize']),
                    num_peaks_per_label=unmapped(sorted_grps[1]['CountingFishNumPeaksPerLabel']))

        save_dots_data.map(fish_counts)


        raw_beads_images_meta = load_raw_images.map(zarr_grp_name=sorted_grps[2],
                                parsed_raw_data_fpath=unmapped(parsed_raw_data_fpath))
        

        filtered_beads_images_metadata = preprocessing_dot_raw_image.map(raw_beads_images_meta,
                            dark_img=unmapped(dark_img),
                            FlatFieldKernel=unmapped(sorted_grps[3]['PreprocessingBeadsRegistrationFlatFieldKernel']),
                            FilteringSmallKernel=unmapped(sorted_grps[3]['PreprocessingBeadsRegistrationFilteringSmallKernel']),
                            LaplacianKernel=unmapped(sorted_grps[3]['PreprocessingBeadsRegistrationFilteringLaplacianKernel']))

        save_images_metadata.map(filtered_beads_images_metadata)
        

        beads_counts = osmFISH_peak_based_detection.map(filtered_beads_images_metadata,
                    min_distance=unmapped(sorted_grps[3]['CountingBeadsRegistrationMinObjDistance']),
                    min_obj_size=unmapped(sorted_grps[3]['CountingBeadsRegistratiohMinObjSize']),
                    max_obj_size=unmapped(sorted_grps[3]['CountingBeadsRegistrationMaxObjSize']),
                    num_peaks_per_label=unmapped(sorted_grps[3]['CountingBeadsRegistrationNumPeaksPerLabel']))

        save_dots_data.map(beads_counts)
        

        # experiment_fpath = Path('/Users/simone/Documents/local_data_storage/prefect_test/whd/exp_pre_auto')
        

        # REGISTRATION RFERENCE CHANNEL
        # Load the counts
        # Run registration
        # Save the output

        # REGISTRATION FISH
        # For each FOV/ hybridization load the corresponding shift
        # Calculate image shift
        
        # CREATE DASK ARRARY OF REGISTERED IMAGES

        # IDENTIFICATION BARCODES AND EXTRACT IMAGES REGION OF BACODE

        # LABEL BARCODES

        # FOR EACH ROUNDS MAP THE PROBABILITY OF 1 OR ZERO

        # RUN BITS CORRECTION AND GENE RECALLING



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

    flow.visualize(flow_state=flow_state)
    cluster.close()