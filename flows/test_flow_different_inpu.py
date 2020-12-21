import prefect
import time
from prefect import task, Flow, Parameter, unmapped
from prefect.engine.executors import DaskExecutor, LocalDaskExecutor


# testing import
from pathlib import Path


@task
def test_parallel(a):
    print('test_parallel')
    time.sleep(1)
    return 1


if __name__ == '__main__':

    # with Flow("test_running",schedule=schedule) as flow:
    with Flow("test_running") as flow:

        a = list(range(2))
        b = list(range(2))
        c = list(range(2))
        d = list(range(2))
        va = test_parallel.map(a)
        vb = test_parallel.map(b)
        vc = test_parallel.map(c)
        vd = test_parallel.map(d)

        z = test_parallel.map(vd)

    executor = DaskExecutor()
    # with raise_on_exception():

    flow_state = flow.run(executor=executor)
    flow.visualize(flow_state=flow_state)

