"""
This prefect flow is used to transfer data from a temporary storage to the 
processing HD.
The flow scan for new data once a day

"""

import prefect
import shutil

from prefect import task, Flow, Parameter, flatten, unmapped
from prefect.engine.executors import DaskExecutor
from prefect.utilities.debug import raise_on_exception
from datetime import timedelta, datetime
from prefect.schedules import IntervalSchedule
from prefect.tasks.control_flow.case import case


from pysmFISH.dask_cluster_utilities_tasks import start_transfering_env
from pysmFISH.utilities_tasks import check_completed_transfer_to_monod, move_data
from pysmFISH.configuration_files_tasks import load_transferring_config

from pysmFISH.logger_utils import setup_logger

@task('determine_free_space')
def free_space(hd_path):
    logger = setup_logger()
    total, used, free = shutil.disk_usage(hd_path)
    free_space_giga = free // (2**30)
    logger.info(f'Free space in the HD: {free_space_giga} Gb')
    if free_space_giga <= 1000:
        return False
    else:
        return True

def transfer_data_flow(transfer_config_fpath):
    
    transfer_config = load_transferring_config(transfer_config_fpath)

    schedule = IntervalSchedule(
    start_date=datetime.utcnow() + timedelta(seconds=1),
    interval=timedelta(days=1),)

    allowed_transfer = free_space(transfer_config['processing_hd_location'])


    cluster = start_transfering_env(processing_hd_location=transfer_config['processing_hd_location'])
    cluster.adapt(minimum_jobs=2)
    executor = DaskExecutor(address=cluster.scheduler_address)

    with case(allowed_transfer, 'True'):
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

            experiment_fpath = check_completed_transfer_to_monod(path_tmp_storage_server,flag_file_key)
            move_data(experiment_fpath, processing_hd_location, flag_file_key)
       
        
        flow.register(project_name="test")
        flow.run_agent()
        flow_state = flow.run(executor=executor)
