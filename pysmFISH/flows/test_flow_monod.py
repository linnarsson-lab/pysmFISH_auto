import prefect
import time
from prefect import task, Flow, Parameter, unmapped
from prefect.engine.executors import DaskExecutor

from pysmFISH.dask_cluster_utilities_tasks import htcondor_cluster_setup


# testing import
from pathlib import Path


@task
def test_parallel(a):
    print('test_parallel')
    time.sleep(200)


if __name__ == '__main__':

    monod_cluster = {
            'cores': 1,
            'memory': '2GB',
            'disk': '0.1GB',
            'local_directory': '/tmp'      
        }

    cluster = htcondor_cluster_setup(monod_cluster)
    cluster.adapt(minimum_jobs=2)
    print(cluster)
    # with Flow("test_running",schedule=schedule) as flow:
    with Flow("test_running") as flow:

        a = list(range(20))
        test_parallel.map(a)


    executor = DaskExecutor(address=cluster.scheduler_address)
    # with raise_on_exception():

    flow_state = flow.run(executor=executor)
    print('done')
    



