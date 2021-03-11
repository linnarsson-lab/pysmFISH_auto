"""
Set of function and graph for preprocessing the
images of the fresh nuclei collected during the 
eel experiment

"""
from typing import *
import numpy as np
import logging
import sys
from pathlib import Path

from pysmFISH.microscopy_file_parsers import single_nikon_nd2_parser_simple
from pysmFISH import preprocessing
from pysmFISH.analysis_clusters import dask_cluster
from pysmFISH.zarr_utilities import load_parsed_img, save_img, load_zarr_file, create_empty_zarr
from pysmFISH.utils import load_running_analysis_config_file, create_dir, convert_to_uint16

class eel_fresh_nuclei_processing():

    def __init__(self,experiment_fpath:str):
        self.experiment_fpath = Path(experiment_fpath)
        self.logger = logging.getLogger(__name__)

    @staticmethod
    def nd2_nuclei_parser(img_fpath:str):
        parser = single_nikon_nd2_parser_simple(img_fpath)
        parser.parse_file()
    
    @staticmethod   
    def nuclei_filtering(processing_tuple):
        img_fpath, filtered_img_fpath,fov,experiment_fpath = processing_tuple
        analysis_config = load_running_analysis_config_file(experiment_fpath)
        img_stack = load_parsed_img(img_fpath,fov)
        large_kernel_size_sigma = analysis_config['analysis_parameters']['fresh_nuclei']['preprocessing']['large_kernel_size_sigma']
        filtered_img_stack = preprocessing.large_structures_filtering(img_stack,large_kernel_size_sigma)
        filtered_img = np.amax(filtered_img_stack,axis=0)
        # save_img(filtered_img_fpath,fov,filtered_img)
        filtered_img_fpath = Path(filtered_img_fpath)
        fname = filtered_img_fpath.parent / filtered_img_fpath.stem / ('fov_' + str(fov) + '.npy')
        # filtered_img = convert_to_uint16(filtered_img)
        np.save(fname,filtered_img)

    def run_computation(self):

        # Load the analysis
        self.analysis_config = load_running_analysis_config_file(self.experiment_fpath)

        try:
            self.raw_img_fpath = list((self.experiment_fpath / 'fresh_nuclei').glob('*.nd2'))[0]
        except:
            self.logger.error(f'missing .nd2 fresh nuclei file')
            sys.exit(f'missing .nd2 fresh nuclei file')

        # Start the filtering cluster
        self.processing_cluster = dask_cluster(self.experiment_fpath,'eel_processing_cluster_setup')
        self.processing_cluster.start()
        
        # Parsing the raw data on a single worker (wait until is done before the filtering)
        parsing_future = self.processing_cluster.client.submit(self.nd2_nuclei_parser,self.raw_img_fpath)
        parsing_future.result()

        try:
            self.parsed_img_fpath = list((self.experiment_fpath / 'fresh_nuclei').glob('*.zarr'))[0]
        except:
            self.logger.error(f'missing pased fresh nuclei file')
            sys.exit(f'missing parsed fresh nuclei file')

        # create zarr file
        self.filtered_fpath = self.raw_img_fpath.parent / (self.raw_img_fpath.stem + '_filtered.zarr')
        create_empty_zarr(str(self.filtered_fpath))
        
        # np storage
        self.npy_filtered_fpath = self.raw_img_fpath.parent / (self.raw_img_fpath.stem + '_filtered')
        create_dir(self.npy_filtered_fpath)

        self.root = load_zarr_file(self.parsed_img_fpath)
        fovs = list(self.root.keys())
        tot_fovs = len(fovs)
        all_experiments_fpath = [self.experiment_fpath]*tot_fovs
        all_filtered_fpath = [self.filtered_fpath]*tot_fovs
        all_img_fpath = [self.parsed_img_fpath]*tot_fovs
        all_processing_tuples = list(zip(all_img_fpath,all_filtered_fpath,fovs,all_experiments_fpath))

        futures_processes=self.processing_cluster.client.map(self.nuclei_filtering,all_processing_tuples)
        self.processing_cluster.client.gather(futures_processes)
        self.processing_cluster.client.close()
        self.processing_cluster.cluster.close()