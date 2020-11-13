import prefect
from prefect import task, Flow, Parameter, flatten, unmapped
from prefect.engine.executors import DaskExecutor
from prefect.environments import LocalEnvironment

from logger_utils import function_logger
import time

def inner():
    logger = function_logger()
    logger.info('i am the inner function--new logger')
    time.sleep(10)

@task(task_run_name=lambda **kwargs: f"testing-logger-writing-logs-{kwargs['x']}-suiname")
def wlog(x):
    from prefect import context
    logger = context.get("logger")
    logger.debug('i am the task')
    inner()
    # logger = prefect_logging_setup('test')
    logger.info(f'start sleep')
    time.sleep(20)
    logger.info(f'done sleep')

a = list(range(10))
with Flow("logging-flow",environment=LocalEnvironment(DaskExecutor(address='tcp://193.10.16.58:14211'))) as flow:
    # logger = prefect.utilities.logging.get_logger()
    # logger.info('this log is generated in the flow')
    out_task = wlog.map(a)
    # logger.info('done')

flow.register(project_name="test")
flow.run_agent()