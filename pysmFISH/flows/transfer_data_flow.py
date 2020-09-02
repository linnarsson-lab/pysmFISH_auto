"""
This prefect flow is used to transfer data from a temporary storage to the 
processing HD.
The flow scan for new data once a day

"""

import prefect
from prefect import task, Flow, Parameter, flatten, unmapped
from prefect.engine.executors import DaskExecutor
from prefect.utilities.debug import raise_on_exception
from datetime import timedelta, datetime
from prefect.schedules import IntervalSchedule


from pysmFISH.dask_cluster_utilities_tasks import start_transfering_env
from pysmFISH.utilities_tasks import check_completed_transfer_to_monod, move_data
from pysmFISH.configuration_files_tasks import load_transferring_config

from pysmFISH.logger_utils import setup_extra_loggers, prefect_logging_setup



def transfer_data_flow(transfer_config_fpath):
    
    transfer_config = load_transferring_config(transfer_config_fpath)


    schedule = IntervalSchedule(
    start_date=datetime.utcnow() + timedelta(seconds=1),
    interval=timedelta(days=1),)


    # with Flow("test_running",schedule=schedule) as flow:
    with Flow("transfering_data_to_processing_directory",schedule=schedule) as flow:

        # -------------------------------------------------------
        # Parameters for folder scanning for determine if an experiment
        # has been transferred to monod
        # -------------------------------------------------------
        path_tmp_storage_server = Parameter('path_tmp_storage_server',default=transfer_config['path_tmp_storage_server'])
        flag_file_key = Parameter('flag_file_key',default=transfer_config['flag_file_key'])

        # Path to the location of the HD where the data are processed
        # In out case we run all the preprocessing in an SSD HD

        processing_hd_location = Parameter('processing_hd_location',default=transfer_config['processing_hd_location'])

        
        cluster = start_transfering_env(processing_hd_location=transfer_config['processing_hd_location'])
        cluster.adapt(minimum_jobs=2)

        
        experiment_fpath = check_completed_transfer_to_monod(path_tmp_storage_server,flag_file_key)
        move_data(experiment_fpath, processing_hd_location, flag_file_key)
       
    executor = DaskExecutor(address=cluster.scheduler_address)

    flow_state = flow.run(executor=executor)



# if __name__ == "__main__":
#     config_db_fpath = '/Users/simone/Documents/local_data_storage/prefect_test/whd/config_db'
#     main(config_db_fpath)
