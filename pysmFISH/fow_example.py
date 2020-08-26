import prefect
from prefect import task, Flow, Parameter
from prefect.engine.executors import DaskExecutor
from prefect.utilities.debug import raise_on_exception
from datetime import timedelta, datetime
from prefect.schedules import IntervalSchedule



from pysmFISH.dask_cluster_utilities import local_cluster_setup

from pysmFISH.microscopy_file_parsers_tasks import nd2_raw_files_selector, nikon_nd2_autoparser
from pysmFISH.qc_tasks import check_matching_metadata_robofish
from pysmFISH.utilities_tasks import check_completed_transfer_to_monod, load_experiment_config_file
from pysmFISH.utilities_tasks import create_raw_tmp_folders
from pysmFISH.data_transfer_tasks import move_data

from pysmFISH.logger_utils import setup_extra_loggers, prefect_logging_setup



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

        
        path_tmp_storage_server = Parameter('path_tmp_storage_server')
        flag_file_key = Parameter('flag_file_key')
        processing_hd_location = Parameter('processing_hd_location')
        

        experiment_fpath = check_completed_transfer_to_monod(path_tmp_storage_server,flag_file_key)
        experiment_fpath = move_data(experiment_fpath, processing_hd_location)
        experiment_info = load_experiment_config_file(experiment_fpath)
        # Get all the .nd2 files to process
        all_raw_files = nd2_raw_files_selector(experiment_fpath=experiment_fpath,upstream_tasks=[setup_extra_loggers])


        # Run the crosscheck for all the pkl files
        check_matching_metadata_robofish(all_raw_files)

        # Create the folders for storing raw data and tmp parsed
        # This folders are independent from the type of experiment
        # that need to be run

        # test = test_parallel(upstream_tasks=[check_matching_metadata_robofish(all_raw_files)])
        # create_raw_tmp_folders(experiment_fpath=experiment_fpath,upstream_tasks=[check_matching_metadata_robofish(all_raw_files)])

        # test2 = test_parallel(upstream_tasks=[test])    
        # test3 = test_parallel(upstream_tasks=[test])

        # map the parsing function on all the files


    executor = DaskExecutor(address=cluster.cluster.scheduler_address)
    # with raise_on_exception():

    parameters=dict(path_tmp_storage_server='/Users/simone/Documents/local_data_storage/prefect_test/',
                flag_file_key= 'transfer_to_monod_completed.txt',
                processing_hd_location = '/Users/simone/Documents/local_data_storage/prefect_test/whd'
    )


    flow_state = flow.run(executor=executor,parameters=parameters)
    flow.visualize(flow_state=flow_state)
    