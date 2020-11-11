import prefect
from prefect import task, Flow, Parameter, flatten, unmapped
from prefect.engine.executors import DaskExecutor
from prefect.utilities.debug import raise_on_exception
from datetime import timedelta, datetime
from prefect.schedules import IntervalSchedule
from prefect import context

from pythonjsonlogger import jsonlogger

from pysmFISH.logger_utils import setup_extra_loggers, prefect_logging_setup, test_write_logs_to_file,simple_writing_logger
from pysmFISH.dask_cluster_utilities_tasks import start_processing_env, local_cluster_setup
from pysmFISH.configuration_files_tasks import load_processing_env_config_file, load_experiment_config_file

import logging
import time
from pathlib import Path
from prefect import Client
from prefect.utilities.debug import is_serializable
from prefect.engine import signals
from prefect.environments import RemoteDaskEnvironment,LocalEnvironment

from pysmFISH.logger_utils import setup_extra_loggers,prefect_logging_setup
import logging

@task(task_run_name=lambda **kwargs: f"testing-logger-writing-logs-{kwargs['x']}-suiname")
def wlog(x):
    
    logger = context.get("logger")
    logger.debug('la culara')
    # logger = prefect_logging_setup('test')
    logger.info(f'start sleep')
    time.sleep(20)
    logger.info(f'done sleep')


a = list(range(10))
# with Flow("test_running",schedule=schedule) as flow:
with Flow("logging-flow",environment=LocalEnvironment(DaskExecutor(address='tcp://193.10.16.58:18938'))) as flow:
    # logger = prefect.utilities.logging.get_logger()
    # logger.info('this log is generated in the flow')
    out_task = wlog.map(a)
    # logger.info('done')

flow.register(project_name="test")