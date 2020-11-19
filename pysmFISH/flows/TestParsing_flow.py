import prefect
from prefect import task, Flow, Parameter, flatten, unmapped
from prefect.engine.executors import DaskExecutor
from prefect.environments import LocalEnvironment
from prefect import Task
from prefect.environments.storage import Local

from pysmFISH.microscopy_file_parsers_tasks import nd2_raw_files_selector
from pysmFISH.utilities_tasks import create_empty_zarr_file



with Flow("parsing",environment=LocalEnvironment(DaskExecutor(address='tcp://193.10.16.58:32833')),
            storage=Local(directory='/home/simone/tmp_code/flows')) as flow:
   
    experiment_fpath = Parameter('experiment_fpath', default = '/wsfish/smfish_ssd/AMEXP20201110_EEL_HumanH1930001V1C_auto')
    parsed_image_tag = Parameter('parsed_image_tag', default = 'img_data')

    # Selected the raw .nd2 files
    nd2_selector = nd2_raw_files_selector()
    all_raw_nd2 = nd2_selector(experiment_fpath)

    # Create the zarr file that will contain the parsed data
    create_zarr = create_empty_zarr_file()
    parsed_raw_data_fpath = create_zarr(experiment_fpath,parsed_image_tag)

flow.register(project_name="test")
# flow.run_agent()