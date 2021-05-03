from typing import *

from datetime import datetime


from pysmFISH.logger_utils import selected_logger
from pysmFISH.processing_cluster_setup import local_cluster_setup, htcondor_cluster_setup

# From Sten cytograph shoji
# https://github.com/linnarsson-lab/cytograph-shoji/blob/6389e8864c755f056ab7c9b51892650e5ed4f040/cytograph/pipeline/workflow.py#L12
def nice_deltastring(delta):
	result = []
	s = delta.total_seconds()
	h = s // 60 // 60
	if h >= 1:
		result.append(f"{int(h)}h")
		s -= h * 60 * 60
	m = s // 60
	if m >= 1:
		result.append(f"{int(m)}m")
		s -= m * 60
	if s >= 1:
		result.append(f"{int(s)}s")
	if len(result) > 0:
		return " ".join(result)
	else:
		return f"{delta.microseconds // 1000} ms"



class Pipleline():

    def __init__(self, name, experiment_fpath):
        self.logger = selected_logger()        
        self.parsed_image_tag = 'img_data'

    def add_step(self):
        pass

    def start_dask_processing_env(self, engine:str = 'htcondor', **kwargs):
        if engine == 'local':
            self.cluster = local_cluster_setup()
            logger = logger.debug(f'started local cluster at ')
        elif engine == 'htcondor':
            cluster = htcondor_cluster_setup(cluster_config_parameters)
            cluster.scale(jobs=1)
            minimum_jobs = 1
            maximum_jobs = 15
            cluster.adapt(minimum_jobs=minimum_jobs,maximum_jobs=maximum_jobs)

    # Add qc to check if the configuration data are available in the 