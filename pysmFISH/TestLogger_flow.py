import prefect
from prefect import task, Flow, Parameter, flatten, unmapped
from prefect.engine.executors import DaskExecutor
from prefect.utilities.debug import raise_on_exception
from prefect.utilities.logging import get_logger
from datetime import timedelta, datetime
from prefect.schedules import IntervalSchedule

from pythonjsonlogger import jsonlogger

from pysmFISH.logger_utils import setup_extra_loggers, prefect_logging_setup, test_write_logs_to_file,simple_writing_logger
from pysmFISH.dask_cluster_utilities_tasks import start_processing_env, local_cluster_setup
from pysmFISH.configuration_files_tasks import load_processing_env_config_file, load_experiment_config_file

import logging
import time
from pathlib import Path

@task(name='testing-logger',log_stdout=True)
def slowmo(x):
    logger = logging.get_logger()
    logger.info(f'start sleep')
    time.sleep(5)
    logger.info(f'done sleep')

@task(name='testing-logger-context',log_stdout=True)
def tarca(x):
    logger = prefect.context.get("logger")
    logger.info(f'start sleep')
    time.sleep(5)
    logger.info(f'done sleep')


def setup_logger():
    fname = '/wsfish/smfish_ssd/LBEXP20201014_EEL_Mouse_2420um_auto/prefect_logs/my.log'
    logger = prefect.context.get("logger")
    file_handler = logging.FileHandler(str(fname),'w+')
    file_handler.setLevel(logging.DEBUG)
    format_str = '%(message)%(levelname)%(name)%(asctime)'
    # formatter = jsonlogger.JsonFormatter(format_str)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s: %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return logger




@task(name='testing-logger-writing-logs')
def wlog(x):
    logger = setup_logger()
    logger.info(f'start sleep')
    time.sleep(5)
    logger.info(f'done sleep')



experiment_fpath = Path('/wsfish/smfish_ssd/LBEXP20201014_EEL_Mouse_2420um_auto')
a = range(10)

if __name__ == '__main__':

    flag_file_key = Parameter('flag_file_key', default='transfer_to_monod_completed.txt')
    processing_hd_location = Parameter('processing_hd_location',default='/wsfish/smfish_ssd')
    # processing_hd_location = Parameter('processing_hd_location',default='/Users/simone/Documents/local_data_storage/prefect_test/whd')

     # get info for submitting the error notification to github
    config_db_fpath = Path(processing_hd_location.default) / 'config_db'
    processing_env_config = load_processing_env_config_file(config_db_fpath)

    experiment_info = load_experiment_config_file(experiment_fpath)

    # Activate processing cluster
    cluster = start_processing_env(processing_env_config,experiment_info)

    # with Flow("test_running",schedule=schedule) as flow:
    with Flow("test_logging") as flow:

        out_task = wlog.map(a)
        out_task = wlog.map(a)

    executor = DaskExecutor(address=cluster.scheduler_address)
    

    with raise_on_exception():
        flow_state = flow.run(executor=executor)
    cluster.close()