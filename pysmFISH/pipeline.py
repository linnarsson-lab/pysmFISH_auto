from typing import *

import os
import dask
import pandas as pd

from pathlib import Path
from datetime import datetime

from dask.distributed import Client
from dask.base import tokenize
from dask import dataframe as dd
from dask import delayed


# Import from pysmFISH package

from pysmFISH import logger_utils
from pysmFISH import utils
from pysmFISH import configuration_files
from pysmFISH import io
from pysmFISH import data_organization
from pysmFISH import microscopy_file_parsers
from pysmFISH import processing_cluster_setup
from pysmFISH import data_models
from pysmFISH import preprocessing
from pysmFISH import fovs_registration
from pysmFISH import barcodes_analysis
from pysmFISH import fov_processing
from pysmFISH import stitching



# utils.nice_deltastring for calculate the time
class Pipeline(object):
    """
    Generale pipeline for running barcoded eel or serial smFISH experiments



    """

    def __init__(self, pipeline_run_name:str, experiment_fpath:str,
                run_type:str = 'new', parsing_type:str = 'original', **kwarg):

        self.pipeline_run_name = pipeline_run_name
        self.experiment_fpath = Path(experiment_fpath)
        self.run_type = run_type
        self.parsing_type = parsing_type
        self.logger = logger_utils.json_logger(self.experiment_fpath / 'logs', pipeline_run_name +'.log')   

        # Collect some of the parameters. If missing a predefined value is assigned
        self.raw_data_folder_storage_path = kwarg.pop('raw_data_folder_storage_path', '/fish/rawdata')
        self.parsed_image_tag = kwarg.pop('parsed_image_tag','img_data')
        self.filtered_image_tag = kwarg.pop('filtered_image_tag','img_data')
        self.dataset_folder_storage_path = kwarg.pop('dataset_folder_storage_path','/fish/fish_datasets')
        self.save_intermediate_steps = kwarg.pop('save_intermediate_steps',False)
        self.store_dataset = kwarg.pop('store_dataset',True)

        # Parameters for processing in htcondor
        self.processing_env_config = {}
        self.processing_env_config['processing_engine'] = kwarg.pop('processing_engine','htcondor')
        self.processing_env_config['cores'] = kwarg.pop('cores',20)
        self.processing_env_config['memory'] = kwarg.pop('memory','200GB')
        self.processing_env_config['disk'] = kwarg.pop('disk','0.1GB')
        self.processing_env_config['local_directory'] = kwarg.pop('local_directory','/tmp')
        self.processing_env_config['logs_directory'] = (self.experiment_fpath / 'logs').as_posix()

        self.dataset_path = kwarg.pop('dataset_path','')

        # Define the experiment folder location in the storage HD
        self.storage_experiment_fpath = (Path(self.raw_data_folder_storage_path) / self.experiment_fpath.stem).as_posix()



    # -----------------------------------
    # PIPELINE STEPS
    # ------------------------------------

    def create_folders_step(self):
        utils.create_folder_structure(self.experiment_fpath, self.run_type)

    
    def load_config_files_step(self):
        self.analysis_parameters = configuration_files.load_analysis_config_file(self.experiment_fpath)
        self.parsed_raw_data_fpath = self.experiment_fpath / (self.experiment_fpath.stem + '_' + self.parsed_image_tag + '.zarr') 
    

    def create_analysis_config_file_from_dataset_step(self):
        configuration_files.create_analysis_config_file_from_dataset(self.experiment_fpath, self.metadata)


    def processing_cluster_init_step(self):
        # Start processing environment
        self.logger.info(f'Starting processing environment')
        self.cluster = processing_cluster_setup.start_processing_env(self.processing_env_config)
        self.client = Client(self.cluster)
        self.logger.debug(f'Dask dashboard info {self.client.scheduler_info()}')

    def nikon_nd2_parsing_graph_step(self):
        microscopy_file_parsers.nikon_nd2_parsing_graph(self.experiment_fpath,
                                    self.parsing_type,parsed_image_tag, 
                                    self.storage_experiment_fpath, self.raw_files_fpath,
                                    self.client)
    

    def prepare_processing_dataset_step(self):
        self.data = data_models.Dataset()
        if self.dataset_path:
            self.data.load_dataset(self.dataset_path)
        else:
            self.data.create_full_dataset_from_zmetadata(self.parsed_raw_data_fpath)
        
        self.metadata = self.data.collect_metadata(self.data.dataset)
        self.grpd_fovs = self.data.dataset.groupby('fov_num')


    def determine_tiles_organization(self):
        self.reference_round = self.analysis_parameters['RegistrationReferenceHybridization']
        self.tiles_org = stitching.organize_square_tiles(self.experiment_fpath,
                                    self.data.dataset,self.metadata,
                                    self.reference_round)
        self.tiles_org.run_tiles_organization()
        self.tile_corners_coords_pxl = self.tiles_org.tile_corners_coords_pxl

    def create_dark_img_step(self):
        utils.create_dark_img(self.experiment_fpath,self.metadata)

    def create_running_functions_step(self):
        self.running_functions = configuration_files.create_function_runner(self.experiment_fpath,self.metadata)

    def processing_barcoded_eel_step(self):
        fov_processing.processing_barcoded_eel_fov_graph(self.experiment_fpath,self.analysis_parameters,
                                    self.running_functions, self.tile_corners_coords_pxl,self.metadata,
                                    self.grpd_fovs,self.save_intermediate_steps, self.client)

    def processing_serial_fish_step(self):
        fov_processing.processing_serial_fish_fov_graph(self.experiment_fpath,self.analysis_parameters,
                                    self.running_functions, self.tile_corners_coords_pxl,self.metadata,
                                    self.grpd_fovs,self.save_intermediate_steps, self.client)

    
    # --------------------------------
    # RUNNING OPTIONS
    # --------------------------------
    
    def run_parsing_only(self):
        """
        Full run from raw images from nikon or parsed images
        """
        self.create_folders_step(self.experiment_fpath)
        self.logger.info(f'Folder structure completed')

        self.load_config_files_step()
        self.processing_cluster_init_step.start_processing_env()

        # Run parsing only if required
        if self.parsing_type != 'no_parsing':
            nikon_nd2_parsing_step()



    def run_full(self):
        """
        Full run from raw images from nikon or parsed images
        """
        self.run_parsing_only()

        # Create or load processing dataset
        self.prepare_processing_dataset_step()

        self.create_analysis_config_file_from_dataset_step()


        # CREATE DARK IMAGES
        self.create_dark_img_step()

        # CREATE RUNNING FUNCTIONS
        # used to select the sets of functions to run preprocessing and dots
        # calling according to the type of processing and sample
        self.create_running_functions_step()

        # CREATE DARK IMAGES
        self.determine_tiles_organization()

        if self.metadata['experiment_type'] == 'eel-barcoded':
            self.processing_barcoded_eel_step()
        elif self.metadata['experiment_type'] == 'smfish-serial':
            self.processing_serial_fish_step()
        else:
            self.logger.error(f"the experiment type {self.metadata['experiment_type']} is unknown")
            sys.exit(f"the experiment type {self.metadata['experiment_type']} is unknown")

    