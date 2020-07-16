import prefect
from prefect import task, Flow, Parameter
from prefect.engine.executors import DaskExecutor


from pysmFISH.dask_cluster_utilities import local_cluster_setup
from pysmFISH.prefect_tasks import setup_extra_loggers, nd2_raw_files_selector


cluster = local_cluster_setup()

with Flow("test_running") as flow:
    experiment_fpath = Parameter('experiment_fpath')
    #setup_extra_loggers()
    all_raw_files = nd2_raw_files_selector(experiment_fpath=experiment_fpath,upstream_tasks=[setup_extra_loggers])
    #select_files_to_parse.set_upstream(setup_extra_loggers)
    


# flow = Flow("test_running")
# experiment_fpath = Parameter('experiment_fpath')
# flow.add_task(setup_extra_loggers)
# flow.add_task(select_files_to_parse)
# select_files_to_parse.set_upstream(setup_extra_loggers, flow=flow)
# select_files_to_parse.bind(experiment_fpath=experiment_fpath, flow=flow)

executor = DaskExecutor(address=cluster.cluster.scheduler_address)
# flow.register(project_name= 'test prefect flows')
# flow.run_agent(token=oQMUEhG3q-jy2Z6rG63Kpw)


flow.run(executor=executor,parameters=dict(experiment_fpath='/Users/simone/Documents/local_data_storage/prefect_test/exp_auto'))