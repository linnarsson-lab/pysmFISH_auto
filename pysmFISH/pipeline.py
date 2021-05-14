"""Module containing the classes used to create pipelines.

The master class Pipeline contains all the option to run
standard eel and smfish experiments. Use the methods as step for
building the pipeline. The methods do not require input  
however the underline functions do and must be included as 
attributes. This may not be the best solution but allows
flexibility and easy to use. Make sure to look at the docstrings
when chaining different modules.

Assert will be used to make sure that the attributes created by
a step and required from another are present

If additional functionality are
required subclass Pipeline and add/replace the different functionality
"""
from typing import *

import os
import dask
import pandas as pd
import numpy as np

from pathlib import Path
from datetime import datetime

from dask.distributed import Client
from dask.base import tokenize
from dask import dataframe as dd
from dask import delayed


# Import from pysmFISH package
import pysmFISH

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
from pysmFISH import qc_utils



# utils.nice_deltastring for calculate the time
class Pipeline(object):
    """
    General pipeline class used for running barcoded eel or serial smFISH experiments.
    The modules are used as steps to build a pipeline

    The exposed attributes allow to modify the properties of the runnning pipline
    on the fly


    Args:
        pipeline_run_name (str): Name of the running pipeline made of current datetime and
                experiment name (ex. '210507_13_25_21_LBEXP20210226_EEL_HE_2100um')
        experiment_fpath (str): Path to the experiment folder
        run_type (str): Define if it is a 'new' or 're-run' (default: new)
        parsing_type (str): Define the parsing type. Can be: 
                        original/no_parsing/reparsing_from_processing_folder/reparsing_from_storage
                        (default: original)

        Optional KWargs:
            raw_data_folder_storage_path (str): Path to the cold storage hard drive (default: /fish/rawdata)
            parsed_image_tag (str): Tag to identify the zarr file with parsed images (default: img_data)
            preprocessed_image_tag (str): Tag to identify the zarr file with preprocessed images 
                                        (default: preprocessed_img_data)
            dataset_folder_storage_path (str): Path to the location where the dataset are stored (default: /fish/fish_datasets)
            save_intermediate_steps (bool): Determine if the processed images will be saved (default: True)
            dataset_path (str): Path to an existing dataset that will be used in the processing

            processing_engine (str): Define the name of the system that will run the processing. Can be local/htcondor
                                    (default htcondor). If engine == local the parameters that define the cluster
                                    will be ignored
            cores (int): Number of cores to use in htcondor (default 20)
            memory (str): Total memory for all the cores (default 200GB)
            disk (str): Size of the spillover disk for dask (default 0.1GB)
            local_directory (str): Directory where to spill over on the node (default /tmp)
            logs_directory: (str): Directory where to store dask and htcondor logs

    Attributes:
        storage_experiment_fpath: Path to folder in the storage HD where to store (or are stored) the raw data for
                                    the current experiment
        parsed_raw_data_fpath: Path to the zarr file containing the raw data
        analysis_parameters: Parameters used for running the analysis
        metadata: Dictionary with the parameters that characterize the current experiment
        preprocessed_zarr_fpath: Path to the zarr file containing the preprocessed images
        cluster: dask cluster running locally or on htcondor
        client: dask client used to manage the cluster
        data: Dataset that define the experiment
        grpd_fovs: Dataset grouped-by field of view number
        reference_round: Round used as starting point for the decoding of the codebooks
        tiles_org: Describe the organization of the tiles in the dataset
        tile_corners_coords_pxl: Coords of the tiles acquired by the microscope in the reference round
        running_functions: Dictionary of functions used to define the type of filtering used for the 
                            fish and reference channels

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
        self.preprocessed_image_tag = kwarg.pop('filtered_image_tag','preprocessed_img_data')
        self.dataset_folder_storage_path = kwarg.pop('dataset_folder_storage_path','/fish/fish_datasets')
        self.save_intermediate_steps = kwarg.pop('save_intermediate_steps',True)
        self.dataset_path = kwarg.pop('dataset_path','')
        # self.store_dataset = kwarg.pop('store_dataset',True)

        # Parameters for processing in htcondor
        self.processing_env_config = {}
        self.processing_env_config['processing_engine'] = kwarg.pop('processing_engine','htcondor')
        self.processing_env_config['cores'] = kwarg.pop('cores',20)
        self.processing_env_config['memory'] = kwarg.pop('memory','200GB')
        self.processing_env_config['disk'] = kwarg.pop('disk','0.1GB')
        self.processing_env_config['local_directory'] = kwarg.pop('local_directory','/tmp')
        self.processing_env_config['logs_directory'] = (self.experiment_fpath / 'logs').as_posix()


        # Define the experiment folder location in the storage HD
        self.storage_experiment_fpath = (Path(self.raw_data_folder_storage_path) / self.experiment_fpath.stem).as_posix()
        self.parsed_raw_data_fpath = self.experiment_fpath / (self.experiment_fpath.stem + '_' + self.parsed_image_tag + '.zarr') 
    

    # -----------------------------------
    # PROCESSING STEPS
    # ------------------------------------

    def create_folders_step(self):
        """
        Create the folder structure used for the data processing. If the folder
        is already present it won't overwrite it

        Folder structure
            - original_robofish_logs: contains all the original robofish logs.
            - extra_files: contains the extra files acquired during imaging.
            - extra_processing_data: contains extra files used in the analysis 
                                                like the dark images for flat field correction.
            - pipeline_config: contains all the configuration files.
            - raw_data: contains the renamed .nd2 files and the corresponding 
                        pickle configuration files. It is the directory that is 
                        backed up on the server.
            - output_figures: contains the reports and visualizations
            - notebooks: will contain potential notebooks used for processing the data
            - probes: will contains the fasta file with the probes used in the experiment
            - tmp: save temporary data
        """

        utils.create_folder_structure(self.experiment_fpath, self.run_type)

    
    def prepare_processing_dataset_step(self):
        """
        If a path to an existing dataset is entered it will be loaded otherwise
        it will create a new dataset that will have all the info that characterize the
        experiment.

        Args:
            zarr_file_path (Path): Path of the file used to build the dataset

        """
        self.data = data_models.Dataset()
        if self.dataset_path:
            try:
                self.data.load_dataset(self.dataset_path)
            except:
                self.logger.error(f"can't load the dataset from {self.dataset_path}")
                sys.exit(f"can't load the dataset from {self.dataset_path}")
        else:
            try:
                self.data.create_full_dataset_from_zmetadata(self.parsed_raw_data_fpath)
            except:
                self.logger.error(f"can't create dataset from {self.parsed_raw_data_fpath}")
                sys.exit(f"can't create dataset from {self.parsed_raw_data_fpath}")
        
        self.metadata = self.data.collect_metadata(self.data.dataset)
        self.grpd_fovs = self.data.dataset.groupby('fov_num')
    
    
    def create_analysis_config_file_from_dataset_step(self):
        """
        Load or create the yaml file with all the parameters for running the analysis. It will first
        load the analysis_config.yaml file present in the pipeline_config folder. If not it 
        will create one using the master template stored in the config_db directory

        The following attributes created by another step must be accessible:
        - metadata

        """
        assert self.metadata, self.logger.error(f'cannot load/create the analysis config because missing metadata attr')
        self.analysis_parameters = configuration_files.create_analysis_config_file_from_dataset(self.experiment_fpath, self.metadata)
        

    def processing_cluster_init_step(self):
        """
        Start the processing dask cluster (dask-jobqueue for htcondor) and client

        """
        # Start processing environment
       
        self.cluster = processing_cluster_setup.start_processing_env(self.processing_env_config)
        self.client = Client(self.cluster)
        self.logger.debug(f'Dask dashboard info {self.client.scheduler_info()}')


    def nikon_nd2_parsing_graph_step(self):
        """
        Run the parsing according to what is specified by the parsing_type
        argument.

        """
        microscopy_file_parsers.nikon_nd2_parsing_graph(self.experiment_fpath,
                                    self.parsing_type,self.parsed_image_tag, 
                                    self.storage_experiment_fpath,
                                    self.client)
    


    def determine_tiles_organization(self):
        """
        Determine the organization of the field of views in the dataset

        The following attributes created by another step must be accessible:
        - analysis_parameters
        - dataset
        - metadata

        """
        assert self.metadata, self.logger.error(f'cannot determine tiles organization because missing metadata attr')
        assert isinstance(self.data.dataset, pd.DataFrame), self.logger.error(f'cannot determine tiles organization because missing dataset attr')
        assert self.analysis_parameters, self.logger.error(f'cannot determine tiles organization because missing analysis_parameters attr')
        self.reference_round = self.analysis_parameters['RegistrationReferenceHybridization']
        self.tiles_org = stitching.organize_square_tiles(self.experiment_fpath,
                                    self.data.dataset,self.metadata,
                                    self.reference_round)
        self.tiles_org.run_tiles_organization()
        self.tile_corners_coords_pxl = self.tiles_org.tile_corners_coords_pxl


    def create_running_functions_step(self):
        """
        Create the dictionary with the function names used to run a specific pipeline as defined
        in metadata['pipeline']

        The following attributes created by another step must be accessible:
        - metadata

        """
        assert self.metadata, self.logger.error(f'cannot create running functions because missing metadata attr')
        self.running_functions = configuration_files.create_function_runner(self.experiment_fpath,self.metadata)



    def processing_barcoded_eel_step(self):
        """
        Create and run a dask delayed task graph used to process barcoded eel experiments
        It runs:
        (1) Image preprocessing
        (2) Dot calling
        (3) Field of view registration
        (4) Barcode decoding
        (5) Registration to the microscope coords
        (6) Consolidate the processed images zarr file metadata
        (7) Create a simple output for quick visualization
        
        The following attributes created by another step must be accessible:
        - metadata
        - analysis_parameters
        - running_functions
        - grpd_fovs
        - client

        """
        assert self.metadata, self.logger.error(f'cannot process eel fovs because missing metadata attr')
        assert self.analysis_parameters, self.logger.error(f'cannot process eel fovs because missing analysis_parameters attr')
        assert self.running_functions, self.logger.error(f'cannot process eel fovs because missing running_functions attr')
        assert self.grpd_fovs, self.logger.error(f'cannot process eel fovs because missing grpd_fovs attr')
        assert self.client, self.logger.error(f'cannot process eel fovs because missing client attr')
    
        fov_processing.processing_barcoded_eel_fov_graph(self.experiment_fpath,self.analysis_parameters,
                                    self.running_functions, self.tile_corners_coords_pxl,self.metadata,
                                    self.grpd_fovs,self.save_intermediate_steps, 
                                    self.preprocessed_image_tag,self.client)



    def processing_serial_fish_step(self):
        """
        Create and run a dask delayed task graph used to process serial smFISH experiments
        It runs:
        (1) Image preprocessing
        (2) Dot calling
        (3) Field of view registration
        (4) Registration to the microscope coords
        (5) Consolidate the processed images zarr file metadata
        (6) Create a simple output for quick visualization
        
        The following attributes created by another step must be accessible:
        - metadata
        - analysis_parameters
        - running_functions
        - grpd_fovs
        - client

        """
        assert self.metadata, self.logger.error(f'cannot process smFISH fovs because missing metadata attr')
        assert self.analysis_parameters, self.logger.error(f'cannot process smFISH fovs because missing analysis_parameters attr')
        assert self.running_functions, self.logger.error(f'cannot process smFISH fovs because missing running_functions attr')
        assert self.grpd_fovs, self.logger.error(f'cannot process smFISH fovs because missing grpd_fovs attr')
        assert self.client, self.logger.error(f'cannot process smFISH fovs because missing client attr')
    
        
        fov_processing.processing_serial_fish_fov_graph(self.experiment_fpath,self.analysis_parameters,
                                    self.running_functions, self.tile_corners_coords_pxl,self.metadata,
                                    self.grpd_fovs,self.save_intermediate_steps, 
                                    self.preprocessed_image_tag,self.client)

    

    def remove_duplicated_dots_graph_step(self, hamming_distance:int=3,same_dot_radius:float=10,
                                            stitching_selected:str= 'microscope_stitched'):

        """
        Function to remove the duplicated barcodes present in the overlapping regions of the
        tiles

        Args:
        ----
        hamming_distance (int): Value to select the barcodes that are passing the 
            screening (< hamming_distance). Default = 3
        same_dot_radius (int): Searching distance that define two dots as identical
            Default = 10
        stitching_selected (str): barcodes coords set where the duplicated dots will be
            removed

        The following attributes created by another step must be accessible:
        - dataset
        - tiles_org
        - client

        """
        assert self.client, self.logger.error(f'cannot remove duplicated dots because missing client attr')
        assert isinstance(self.data.dataset, pd.DataFrame), self.logger.error(f'cannot remove duplicated dots because missing dataset attr')
        assert isinstance(self.tiles_org,pysmFISH.stitching.organize_square_tiles), \
                            self.logger.error(f'cannot remove duplicated dots because tiles_org is missing attr')
        
        self.hamming_distance = hamming_distance / self.metadata['barcode_length']
        self.same_dot_radius = same_dot_radius
        self.stitching_selected = stitching_selected
        stitching.remove_duplicated_dots_graph(self.experiment_fpath,self.data.dataset,self.tiles_org,
                                    self.hamming_distance,self.same_dot_radius, 
                                    self.stitching_selected, self.client)



    def processing_fresh_tissue_step(self,parsing=True,
                                    tag_ref_beads='_ChannelEuropium_Cy3_',
                                    tag_nuclei='_ChannelCy3_'):
        """
        This function create and run a processing graph that parse and filter the nuclei staining in fresh tissue
        and parse, filter and counts the Europium beads used for the registration with the smFISH images at high
        power.

        Args:
        ----
        tag_ref_beads (str): The tag reference of the .nd2 file containing the raw beads images. Default: '_ChannelEuropium_Cy3_'
        tag_ref_nuclei (str): The tag reference of the .nd2 file containing the raw images of nuclei. Default: '_ChannelCy3_'

        """
        assert self.client, self.logger.error(f'cannot process fresh tissue because missing client attr')
        assert self.analysis_parameters, self.logger.error(f'cannot process fresh tissue because missing analysis_parameters attr')
        assert self.running_functions, self.logger.error(f'cannot process fresh tissue because missing running_functions attr')
        fov_processing.process_fresh_sample_graph(self.experiment_fpath,self.running_functions,
                                                self.analysis_parameters, self.client,
                                                tag_ref_beads= tag_ref_beads,
                                                tag_nuclei= tag_nuclei,
                                                parsing=parsing)


    # --------------------------------
    # QC STEPS
    # --------------------------------


    def QC_matching_nd2_metadata_robofish_step(self):
        """
        This function is used to check that each of the nd2 files
        generated by the microscope has a matching pkl metadata
        file generated by robofish
        """
        all_raw_files = microscopy_file_parsers.nd2_raw_files_selector(self.experiment_fpath)
        qc_utils.QC_matching_nd2_metadata_robofish(all_raw_files)

    def QC_check_experiment_yaml_file_step(self):
        """
        This QC function check in the fields required in the ExperimentName_config.yaml file
        are present and have the expected values. This file is important for the parsing of the
        data.
        The following keys are required:
        'Stitching_type', 'Experiment_type', 'Barcode_length', 'Barcode', 'Codebook', 'Machine', 'Operator',
        'Overlapping_percentage','Probes_FASTA_name','Species','Start_date','Strain','Tissue','Pipeline'
        """

        qc_utils.QC_check_experiment_yaml_file(self.experiment_fpath)



    def QC_registration_error_step(self):
        """
        Visualise the error in the registration. It plots the error for each fov and the number of matching
        beads identified after the registration. The FOV number, the round number with the lowest number
        of matching beads and the lowest number of matching beads.

        The following attributes created by another step must be accessible:
        - analysis_parameters
        - tile_corners_coords_pxl
        - client

        """
        assert self.client, self.logger.error(f'cannot run QC on registration because missing client attr')
        assert self.analysis_parameters, self.logger.error(f'cannot run QC on registration because missing analysis_parameters attr')
        assert isinstance(self.tile_corners_coords_pxl, np.ndarray), self.logger.error(f'cannot run QC on registration because missing tile_corners_coords_pxl attr')
        
        qc_reg = qc_utils.QC_registration_error(self.client, self.experiment_fpath, 
                    self.analysis_parameters, self.tile_corners_coords_pxl)

        qc_reg.run_qc()
    
    # --------------------------------
    # RUNNING OPTIONS
    # --------------------------------
    

    def run_parsing_only(self):
        """
        Pipeline running the data organization and the parsing
        of the .nd2 files of the entire experiment
        """
        self.logger.info(f"Start parsing")
        start = datetime.now()

        self.create_folders_step(self.experiment_fpath)
        self.logger.info(f'Folder structure completed')

        self.create_analysis_config_file_from_dataset_step()
        self.logger.info(f'Created analysis_config.yaml file')

        self.processing_cluster_init_step.start_processing_env()
        self.logger.info(f'Started dask processing cluster')

        # Run parsing only if required
        self.logger.info(f'Parsing started')
        if self.parsing_type != 'no_parsing':
            nikon_nd2_parsing_step()
        end = datetime.now()
        self.logger.info(f"{fname}: Parsing completed in {nice_deltastring(end - start)}.")
        self.logger.info(f"")

        self.logger.info(f'Started creation of the dataset')
        self.prepare_processing_dataset_step()
        end = datetime.now()
        self.logger.info(f"{fname}: Dataset creation completed in {nice_deltastring(end - start)}.")
        self.logger.info(f"")

    
    def run_required_steps(self):
        """
        Short pipeline used to make sure that the basic required
        step are run and will be included in more complex pipelines

        """
        start = datetime.now()
        self.create_analysis_config_file_from_dataset_step()
        self.logger.info(f'Created analysis_config.yaml file')
        self.determine_tiles_organization()
        self.logger.info(f'Determined the tiles organization')
        self.create_running_functions_step()
        self.logger.info(f'Created the running function dictionary')
        end = datetime.now()
        self.logger.info(f"{fname}: Required steps completed in {nice_deltastring(end - start)}.")
        self.logger.info(f"")
    
    def run_full(self):
        """
        Full run from raw images from nikon or parsed images
        """
        
        self.run_parsing_only()
        self.run_required_steps()    
        
        self.create_preprocessed_zarr_step()
        self.create_dark_img_step()

        if self.metadata['experiment_type'] == 'eel-barcoded':
            self.processing_barcoded_eel_step()
        elif self.metadata['experiment_type'] == 'smfish-serial':
            self.processing_serial_fish_step()
        else:
            self.logger.error(f"the experiment type {self.metadata['experiment_type']} is unknown")
            sys.exit(f"the experiment type {self.metadata['experiment_type']} is unknown")

    