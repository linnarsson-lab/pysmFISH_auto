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
from pysmFISH.io import load_analysis_parameters
from pysmFISH.microscopy_file_parsers_tasks import nd2_raw_files_selector_general
from pysmFISH.utilities_tasks import create_empty_zarr_file
from pysmFISH.data_model import create_shoji_db
from pysmFISH.microscopy_file_parsers_tasks import nikon_nd2_reparser_zarr

import logging
import time
from pathlib import Path
from prefect import Client
from prefect.utilities.debug import is_serializable
from prefect.engine import signals

def setup_logger():
    fname = '/wsfish/smfish_ssd/LBEXP20201014_EEL_Mouse_2420um_auto/prefect_logs/my.log'
    logger = prefect.context.get("logger")
    file_handler = logging.FileHandler(str(fname),'a')
    file_handler.setLevel(logging.DEBUG)
    format_str = '%(message)%(levelname)%(name)%(asctime)'
    formatter = jsonlogger.JsonFormatter(format_str)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return logger


@task(task_run_name=lambda **kwargs: f"testing-logger-writing-logs-{kwargs['x']}-suiname")
def wlog(x):
    logger = setup_logger()
    logger.info(f'start sleep')
    time.sleep(5)
    logger.info(f'done sleep')

@task
def hello_task():
    logger = prefect.context.get("logger")
    logger.info("Hello, Cloud!")

a = range(10)
experiment_fpath = '/wsfish/smfish_ssd/LBEXP20201014_EEL_Mouse_2420um_auto'

flag_file_key = Parameter('flag_file_key', default='transfer_to_monod_completed.txt')
processing_hd_location = Parameter('processing_hd_location',default='/wsfish/smfish_ssd')
# processing_hd_location = Parameter('processing_hd_location',default='/Users/simone/Documents/local_data_storage/prefect_test/whd')

# get info for submitting the error notification to github
config_db_fpath = processing_hd_location.default + '/config_db'
processing_env_config = load_processing_env_config_file(config_db_fpath)
experiment_info = load_experiment_config_file(experiment_fpath)
cluster = start_processing_env(processing_env_config,experiment_info)

experiment_info = load_experiment_config_file(experiment_fpath)
executor = DaskExecutor(address=cluster.scheduler_address)

with Flow("thello-flow") as flow:
    
    experiment_fpath = Parameter('experiment_fpath',default=experiment_fpath)
    experiment_info = Parameter('experiment_info',default=experiment_info)

    # Create the shoji database that will contain the data
    ref = create_shoji_db(experiment_info)
    # Get the list of raw image groups to preprocess
    analysis_parameters = load_analysis_parameters(experiment_name=experiment_info['EXP_name'],upstream_tasks=[ref])
    
    # Reparsing .nd2 files stored in the raw_data subfolder
    raw_files_fpath = Parameter('raw_files_fpath',default=(experiment_fpath.default + '/raw_data'))
    all_raw_files = nd2_raw_files_selector_general(folder_fpath=raw_files_fpath,upstream_tasks=[raw_files_fpath])
    
    # # Run the crosscheck for all the pkl files
    # check_matching_metadata_robofish(all_raw_files)
    # report_input_files_errors(git_repo,experiment_fpath,git_token)
    # # Parse .nd2 files
    tag = 'img_data'
    parsed_raw_data_fpath = create_empty_zarr_file(experiment_fpath,tag,upstream_tasks=[all_raw_files])
    autoparser = nikon_nd2_reparser_zarr.map(nd2_file_path=all_raw_files,parsed_raw_data_fpath=unmapped(parsed_raw_data_fpath),
                                experiment_info=unmapped(experiment_info))
    

    # out_task = wlog.map(a)

flow.register(project_name="test")
flow.run_agent()
flow_state = flow.run(executor=executor)