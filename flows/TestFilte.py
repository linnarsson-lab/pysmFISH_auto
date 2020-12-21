import prefect
from prefect import task, Flow, Parameter, flatten, unmapped
from prefect.engine.executors import DaskExecutor
from prefect import Task
from prefect.storage.local import Local
from prefect.run_configs import LocalRun

from prefect.utilities.debug import raise_on_exception

import time

@task
def test_parallel(a):
    print('test_parallel')
    time.sleep(20)




# adapt_kwargs={"minimum": 50}
# address='tcp://193.10.16.58:18049',
with Flow("filtering-counting",run_config=LocalRun(), executor= DaskExecutor(address='tcp://193.10.16.58:3111'),storage=Local()) as flow:
   

    all_data = Parameter('all_data',default = list(range(10)))
    test_parallel.map(all_data)

    
# with raise_on_exception():
#     flow.run()

flow.register(project_name="test")
