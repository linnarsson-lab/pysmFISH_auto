import prefect
from prefect import task, Flow, Parameter, flatten, unmapped
from prefect.engine.executors import DaskExecutor
from prefect.environments import LocalEnvironment
from prefect import Task
from prefect.environments.storage import Local

from logger_utils import function_logger
import time

def inner_log():
    logger = logger = prefect.utilities.logging.get_logger()
    logger.info('i am the inner function--new logger')
    time.sleep(10)

@task(task_run_name=lambda **kwargs: f"testing-logger-writing-logs-{kwargs['x']}-suiname")
def wlog(x):
    from prefect import context
    logger = context.get("logger")
    logger.debug('i am the task')
    inner_log()
    # logger = prefect_logging_setup('test')
    logger.info(f'start sleep')
    time.sleep(20)
    logger.info(f'done sleep')


class AddTask(Task):
    def run(self, x):
        self.x = x
        time.sleep(20)
        

a = list(range(10))
with Flow("logging-flow",environment=LocalEnvironment(DaskExecutor(address='tcp://193.10.16.58:32833')),
            storage=Local(directory='/home/simone/tmp_code/flows')) as flow:
    # logger = prefect.utilities.logging.get_logger()
    # logger.info('this log is generated in the flow')
    # out_task = wlog.map(a)
    # logger.info('done')
    add_task = AddTask(task_run_name=lambda **kwargs: f"testing-logger-writing-logs-{kwargs['x']}-suiname")
    add_task.map(a)

flow.register(project_name="test")
# flow.run_agent()