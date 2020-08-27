import prefect
from prefect import task, Flow, Parameter, flatten
from prefect.engine.executors import DaskExecutor
from prefect.utilities.debug import raise_on_exception
from datetime import timedelta, datetime
from prefect.schedules import IntervalSchedule



from pysmFISH.dask_cluster_utilities import local_cluster_setup

from pysmFISH.microscopy_file_parsers_tasks import nd2_raw_files_selector, nikon_nd2_autoparser
from pysmFISH.qc_tasks import check_matching_metadata_robofish
from pysmFISH.utilities_tasks import check_completed_transfer_to_monod, load_experiment_config_file
from pysmFISH.utilities_tasks import create_folder_structure, collect_extra_files,load_data_array
from pysmFISH.data_transfer_tasks import move_data

from pysmFISH.logger_utils import setup_extra_loggers, prefect_logging_setup



# testing import
from pathlib import Path


@task
def test_parallel():
    print('test_parallel')
    return 22


if __name__ == '__main__':

    # Add all the components to check for _auto files in the folder
    # or add it in the flow below

    # schedule = IntervalSchedule(
    # start_date=datetime.utcnow() + timedelta(seconds=1),
    # interval=timedelta(minutes=1),)


    # with Flow("test_running",schedule=schedule) as flow:
    with Flow("test_running") as flow:

        # Run function to select cluster locally or htcondor
        cluster = local_cluster_setup()

        # Define the parameters
        
        
        # -------------------------------------------------------
        # Parameters for folder scanning for determine if an experiment
        # has been transferred to monod
        # -------------------------------------------------------
        # path_tmp_storage_server = Parameter('path_tmp_storage_server')
        # flag_file_key = Parameter('flag_file_key')

        # Path to the location of the HD where the data are processed
        # In out case we run all the preprocessing in an SSD HD
 
        # processing_hd_location = Parameter('processing_hd_location')
        

        # GOOD
        # experiment_fpath = check_completed_transfer_to_monod(path_tmp_storage_server,flag_file_key)
        # experiment_fpath = move_data(experiment_fpath, processing_hd_location)
        # experiment_info = load_experiment_config_file(experiment_fpath)
        # create_folder_structure(experiment_fpath)
        # collect_extra_files(experiment_fpath=experiment_fpath,experiment_info=experiment_info)
        #  Get all the .nd2 files to process
        # all_raw_files = nd2_raw_files_selector(experiment_fpath=experiment_fpath)


        # # Run the crosscheck for all the pkl files
        # check_matching_metadata_robofish(all_raw_files)

        experiment_fpath = Path('/Users/simone/Documents/local_data_storage/prefect_test/whd/exp_pre_auto')
        all_raw_files = list(experiment_fpath.glob('*.nd2'))
        
        # PARSING
        output_parsing = nikon_nd2_autoparser.map(all_raw_files)

        # LOAD DATA ARRAY OF RAW IMAGES
        img_data_arrays = load_data_array.map(input_tuple=flatten(output_parsing))

        # FILTERING
        # NB when it works we can do different filtering depending on the type of tissue to run
        # use the case function in prefect

        # test = test_parallel(upstream_tasks=[check_matching_metadata_robofish(all_raw_files)])
        # create_raw_tmp_folders(experiment_fpath=experiment_fpath,upstream_tasks=[check_matching_metadata_robofish(all_raw_files)])

        # test2 = test_parallel(upstream_tasks=[test])    
        # test3 = test_parallel(upstream_tasks=[test])

        # map the parsing function on all the files


    executor = DaskExecutor(address=cluster.cluster.scheduler_address)
    # with raise_on_exception():

    # parameters=dict(path_tmp_storage_server='/Users/simone/Documents/local_data_storage/prefect_test/',
    #             flag_file_key= 'transfer_to_monod_completed.txt',
    #             processing_hd_location = '/Users/simone/Documents/local_data_storage/prefect_test/whd',
    #             )


    # flow_state = flow.run(executor=executor,parameters=parameters)
    flow_state = flow.run(executor=executor)
    
    flow.visualize(flow_state=flow_state)
    









# Path to the location of the HD where the data are processed
# In out case we run all the preprocessing in an SSD HD
 

 
# -------------------------------------------------------
# Parameter for loading extra processing files 
# -------------------------------------------------------

# codebooks location
codebooks_folder: '/fish/work_sd/codebooks'

# dark images for correction location
dark_imgs_folder: '/fish/work_sd/dark_imgs'

# folder with probes set and probes scheme
probes_info_folder: '/fish/work_sd/probes_info'


# -------------------------------------------------------
# Parameter for output and data storage 
# -------------------------------------------------------

# raw data long term storage
rawdata_storage_location: '/fish/rawdata'

# results / output long term storage
results_storage_location: '/fish/results'
