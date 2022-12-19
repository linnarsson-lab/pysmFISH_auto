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
import pickle
import os
import dask
import sys
import gc
import yaml
import time
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
from pysmFISH import data_organization
from pysmFISH import processing_cluster_setup
from pysmFISH import segmentation
from pysmFISH import FOV_alignment

# utils.nice_deltastring for calculate the time
class Pipeline:
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

    __Optional KWargs__:
        raw_data_folder_storage_path (str): Path to the cold storage hard drive (default: /fish/rawdata)
        parsed_image_tag (str): Tag to identify the zarr file with parsed images (default: img_data)
        preprocessed_image_tag (str): Tag to identify the zarr file with preprocessed images
                                    (default: preprocessed_img_data)
        dataset_folder_storage_path (str): Path to the location where the dataset are stored (default: /fish/fish_datasets)
        results_folder_storage_path (str): Path to the location where the dataset are stored (default: /fish/fish_results)
        save_intermediate_steps (bool): Determine if the processed images will be saved (default: True)
        dataset_path (str): Path to an existing dataset that will be used in the processing
        chunk_size (int): Number of FOV to process in parallel
        same_dot_radius_duplicate_dots (float): Searching distance that define two dots as identical
                                (default: 5)
        stitching_selected (str): Define the stitched counts on which the overlapping dotes will be removed
                                (default: microscope_stitched)
        hamming_distance (int): Value to select the barcodes that are passing the
                                screening (< hamming_distance). (default: 3)


        processing_engine (str): Define the name of the system that will run the processing. Can be local/htcondor
                                (default htcondor). If engine == local the parameters that define the cluster
                                will be ignored
        cores (int): Number of cores/job to use in htcondor or in the local processing (default 20). In the
                        the unmanaged cluster correspond to the nummber of core for each process (nprocs)
        memory (str): Total memory for all the cores in condor (default 200GB) or per core in local setup
                    or per process (nprocs) in the unmanaged cluster
        disk (str): Size of the spillover disk for dask (default 0.1GB)
        local_directory (str): Directory where to spill over on the node (default /tmp)
        logs_directory: (str): Directory where to store dask and htcondor logs
        adaptive: (bool): Decide if the cluster can increase/decrease the number of worker accroding to
                            the processig required. (default True)
        maximum_jobs (int): Max number of jobs to run in htcondor
        scheduler_port (int): define the dask scheduler port. Used for the unmanaged cluster (default 23875)
        dashboard_port (int): define the dask dashboard port: Used for the unmanaged cluser (default 8787)
        scheduler_address (str): Address of the dask scheduler. Used for the unmanaged cluser.
                            'localhost' if running of the main node (default 'localhost)

        workers_addresses_list (list[str]): Addresses of the workers (default [monod10,monod11,monod12,monod33])

        nprocs (int): number of processes for each workers (unmanaged cluster) (default 40 for single node monod)
        nthreads (int): number threads/process (default 1)
        save_bits_int: (bool): Save the intensity of the bits and the flipping direction
        start_from_preprocessed_imgs (bool): Run the processing starting from the counting
            using preprocessed images. default: False
        resume: (bool): Restart the processsing. Determine automatically which files are already processed by checking
                        the *_*decoded_* files in the results folder
        reuse_cluster (str): Connect the pipeline to a previously created cluster (default False). Can be: 'connect_to_client' ,'connect_to_scheduler'
        active_cluster (dask_cluster): Already active cluster to reconnect to when you want to reuse a cluster
                                        (default None)
        active_client (dask_client): Already active client to reconnect to when you want to reuse a cluster
                                        (default None)
        active_scheduler_address (str): Running cluster to connect when you want reuse a cluster

        fresh_tissue_segmentation_engine (str): Stardist is the default because faster and doesn't require diameter
        diameter_size (int): Size of the diameter of the cells to segment using cellpose
        min_overlapping_pixels_segmentation (int): Size of the overlapping label objects
        fov_alignment_mode (str): 'clip' or 'merge' (default clipped).
        clip_size (str): only if fov_alignment_mode is merge it will clip the clip_size length of borders.
        remove_distinct_genes (bool): when stitching it will also remove overlapping dots of different genes if set to
            true. Defaults to true
        bead_alignment_centering_mode (str): Mode to use for centering the bead alignment. Choose from: 'scan', 'middle'
            or 'mean'.
        bead_alignment_radius (float): Search radius to look for corresponding beads. Defaults to 2000.


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
    print('version cpu')
    def __init__(
        self,
        pipeline_run_name: str,
        experiment_fpath: str,
        run_type: str = "new",
        parsing_type: str = "original",
        **kwarg,
    ):
        self.pipeline_run_name = pipeline_run_name
        self.experiment_fpath = Path(experiment_fpath)
        self.run_type = run_type
        self.parsing_type = parsing_type

        date_tag = time.strftime("%y%m%d_%H_%M_%S")
        self.pipeline_run_name = date_tag + "_" + self.experiment_fpath.stem

        # Collect some of the parameters. If missing a predefined value is assigned
        self.raw_data_folder_storage_path = kwarg.pop(
            "raw_data_folder_storage_path", "/fish/raw_data"
        )
        self.parsed_image_tag = kwarg.pop("parsed_image_tag", "img_data")
        self.preprocessed_image_tag = kwarg.pop(
            "filtered_image_tag", "preprocessed_img_data"
        )
        self.dataset_folder_storage_path = kwarg.pop(
            "dataset_folder_storage_path", "/fish/fish_datasets"
        )
        self.results_folder_storage_path = kwarg.pop(
            "results_folder_storage_path", "/fish/fish_results"
        )
        self.save_intermediate_steps = kwarg.pop("save_intermediate_steps", True)
        self.dataset_path = kwarg.pop("dataset_path", "")
        self.store_dataset = kwarg.pop("store_dataset", True)
        self.chunk_size = kwarg.pop("chunk_size", 20)
        self.same_dot_radius_duplicate_dots = kwarg.pop(
            "same_dot_radius_duplicate_dots", 5
        )
        self.stitching_selected = kwarg.pop("stitching_selected", "microscope_stitched")

        self.fresh_tissue_segmentation_engine = kwarg.pop(
            "fresh_tissue_segmentation_engine", "stardist"
        )
        self.diameter_size = kwarg.pop("diameter_size", 25)

        self.hamming_distance = kwarg.pop("hamming_distance", 3)
        self.save_bits_int = kwarg.pop("save_bits_int", True)
        self.adaptive = kwarg.pop("adaptive", True)
        self.maximum_jobs = kwarg.pop("maximum_jobs", 15)
        self.scheduler_port = kwarg.pop("scheduler_port", 8786)
        self.dashboard_port = kwarg.pop("dashboard_port", 8787)
        self.nprocs = kwarg.pop("nprocs", 40)
        self.nthreads = kwarg.pop("nthreads", 1)
        self.scheduler_address = kwarg.pop("scheduler_address", "localhost")
        self.workers_addresses_list = kwarg.pop(
            "workers_addresses_list", ["monod10", "monod11", "monod12", "monod33"]
        )
        self.reuse_cluster = kwarg.pop("reuse_cluster", False)
        self.active_client = kwarg.pop("active_client", None)
        self.active_cluster = kwarg.pop("active_cluster", None)
        self.active_scheduler_address = kwarg.pop("active_scheduler_address", None)

        self.start_from_preprocessed_imgs = kwarg.pop(
            "start_from_preprocessed_imgs", False
        )
        self.resume = kwarg.pop("resume", False)
        self.processing_engine = kwarg.pop("processing_engine", "htcondor")

        # Parameters for processing in htcondor
        self.processing_env_config = {}
        self.processing_env_config["processing_engine"] = self.processing_engine
        self.processing_env_config["cores"] = kwarg.pop("cores", 20)
        self.processing_env_config["memory"] = kwarg.pop("memory", "200GB")
        self.processing_env_config["disk"] = kwarg.pop("disk", "0.1GB")
        self.processing_env_config["local_directory"] = kwarg.pop(
            "local_directory", "/tmp"
        )
        self.processing_env_config["logs_directory"] = (
            self.experiment_fpath / "logs"
        ).as_posix()
        self.processing_env_config["adaptive"] = self.adaptive
        self.processing_env_config["maximum_jobs"] = self.maximum_jobs
        self.processing_env_config["scheduler_port"] = self.scheduler_port
        self.processing_env_config["dashboard_port"] = self.dashboard_port
        self.processing_env_config["scheduler_address"] = self.scheduler_address
        self.processing_env_config[
            "workers_addresses_list"
        ] = self.workers_addresses_list
        self.processing_env_config["nprocs"] = self.nprocs
        self.processing_env_config["nthreads"] = self.nthreads

        # Define the experiment folder location in the storage HD
        self.storage_experiment_fpath = (
            Path(self.raw_data_folder_storage_path) / self.experiment_fpath.stem
        ).as_posix()
        self.parsed_raw_data_fpath = self.experiment_fpath / (
            self.experiment_fpath.stem + "_" + self.parsed_image_tag + ".zarr"
        )

        self.min_overlapping_pixels_segmentation = kwarg.pop(
            "min_overlapping_pixels_segmentation", 20
        )

        self.max_expansion_radius = kwarg.pop("max_expansion_radius", 18)
        self.fov_alignment_mode = kwarg.pop("fov_alignment_mode", "clip")
        self.clip_size = kwarg.pop("clip_size", 0)
        self.remove_distinct_genes = kwarg.pop("remove_distinct_genes", False)
        self.bead_alignment_centering_mode = kwarg.pop("bead_alignment_centering_mode", 'scan')
        self.bead_alignment_radius = kwarg.pop("bead_alignment_radius", 5000)
        self.keep_intensities = False

    # -----------------------------------
    # PROCESSING STEPS
    # ------------------------------------

    def create_folders_step(self):
        """
        Create the folder structure used for the data processing. If the folder
        is already present it won't overwrite it

        Folder structure
            - codebook: contain the codebooks used for the analysis
            - original_robofish_logs: contains all the original robofish logs.
            - extra_files: contains the extra files acquired during imaging.
            - extra_processing_data: contains extra files used in the analysis
                like the dark images for flat field correction.
            - pipeline_config: contains all the configuration files.
            - raw_data: contains the renamed .nd2 files and the corresponding .pkl metadata
                                    files
            - output_figures: contains the reports and visualizations
            - notebooks: contain potential notebooks used for processing the data
            - probes: contains the fasta file with the probes used in the experiment
            - fresh_tissue: contain the images and the process data obtained from
                                                    imaging the fresh tissue before eel processing
            - logs: contains the dask and htcondor logs
            - microscope_tiles_coords: contain the coords of the FOVs according to the
                                                                                            microscope stage.
            - results: contains all the processing results
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
                self.logger.error(
                    f"can't create dataset from {self.parsed_raw_data_fpath}"
                )
                sys.exit(f"can't create dataset from {self.parsed_raw_data_fpath}")

        self.metadata = self.data.collect_metadata(self.data.dataset)
        self.grpd_fovs = self.data.dataset.groupby("fov_num")

    def create_analysis_config_file_from_dataset_step(self):
        """
        Load or create the yaml file with all the parameters for running the analysis. It will first
        load the analysis_config.yaml file present in the pipeline_config folder. If not it
        will create one using the master template stored in the config_db directory

        The following attributes created by another step must be accessible:
        - metadata

        """
        assert self.metadata, self.logger.error(
            f"cannot load/create the analysis config because missing metadata attr"
        )
        self.analysis_parameters = (
            configuration_files.create_analysis_config_file_from_dataset(
                self.experiment_fpath, self.metadata
            )
        )

    def processing_cluster_init_step(self):
        """Create new cluster and client or reuse a cluster and client previously created
        Can connect direclty to the client/cluster using the 'connect_to_client' flag or
        to the scheduler with the 'connect_to_scheduler'.
        """

        # Start processing environment
        if self.reuse_cluster == "connect_to_client":
            self.cluster = self.active_cluster
            self.client = self.active_client
            self.client.run(gc.collect)
            self.client.run(utils.trim_memory())
        elif self.reuse_cluster == "connect_to_scheduler":
            self.client = Client(self.active_scheduler_address)
            self.client.run(gc.collect)
            self.client.run(utils.trim_memory())
        else:
            self.cluster = processing_cluster_setup.start_processing_env(
                self.processing_env_config
            )
            self.client = Client(self.cluster, asynchronous=True)
            # self.logger.debug(f'Dask dashboard info {self.client.scheduler_info()}')

    def nikon_nd2_parsing_graph_step(self):
        """
        Run the parsing according to what is specified by the parsing_type
        argument.

        """
        microscopy_file_parsers.nikon_nd2_parsing_graph(
            self.experiment_fpath,
            self.parsing_type,
            self.parsed_image_tag,
            self.storage_experiment_fpath,
            self.client,
        )

    def determine_tiles_organization(self):
        """
        Determine the organization of the field of views in the dataset

        The following attributes created by another step must be accessible:
        - analysis_parameters
        - dataset
        - metadata

        """
        assert self.metadata, self.logger.error(
            f"cannot determine tiles organization because missing metadata attr"
        )
        assert isinstance(self.data.dataset, pd.DataFrame), self.logger.error(
            f"cannot determine tiles organization because missing dataset attr"
        )
        assert self.analysis_parameters, self.logger.error(
            f"cannot determine tiles organization because missing analysis_parameters attr"
        )
        self.reference_round = self.analysis_parameters[
            "RegistrationReferenceHybridization"
        ]
        self.tiles_org = stitching.organize_square_tiles(
            self.experiment_fpath,
            self.data.dataset,
            self.metadata,
            self.reference_round,
        )
        self.tiles_org.run_tiles_organization()
        self.tile_corners_coords_pxl = self.tiles_org.tile_corners_coords_pxl

    def determine_tiles_organization_before_room_reorganisation(self):
        """
        Determine the organization of the field of views in the dataset

        The following attributes created by another step must be accessible:
        - analysis_parameters
        - dataset
        - metadata

        """
        assert self.metadata, self.logger.error(
            f"cannot determine tiles organization because missing metadata attr"
        )
        assert isinstance(self.data.dataset, pd.DataFrame), self.logger.error(
            f"cannot determine tiles organization because missing dataset attr"
        )
        assert self.analysis_parameters, self.logger.error(
            f"cannot determine tiles organization because missing analysis_parameters attr"
        )
        self.reference_round = self.analysis_parameters[
            "RegistrationReferenceHybridization"
        ]
        self.tiles_org = stitching.organize_square_tiles_old_room(
            self.experiment_fpath,
            self.data.dataset,
            self.metadata,
            self.reference_round,
        )
        self.tiles_org.run_tiles_organization()
        self.tile_corners_coords_pxl = self.tiles_org.tile_corners_coords_pxl

    def create_running_functions_step(self):
        """
        Create the dictionary with the function names used to run a specific pipeline as defined
        in metadata['pipeline']

        The following attributes created by another step must be accessible:
        - metadata

        """
        assert self.metadata, self.logger.error(
            f"cannot create running functions because missing metadata attr"
        )
        self.running_functions = configuration_files.create_function_runner(
            self.experiment_fpath, self.metadata
        )

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
        assert self.metadata, self.logger.error(
            f"cannot process eel fovs because missing metadata attr"
        )
        assert self.analysis_parameters, self.logger.error(
            f"cannot process eel fovs because missing analysis_parameters attr"
        )
        assert self.running_functions, self.logger.error(
            f"cannot process eel fovs because missing running_functions attr"
        )
        assert self.grpd_fovs, self.logger.error(
            f"cannot process eel fovs because missing grpd_fovs attr"
        )
        assert self.client, self.logger.error(
            f"cannot process eel fovs because missing client attr"
        )
        assert self.tiles_org, self.logger.error(
            f"cannot process eel fovs because missing tiles organization attr"
        )

        fov_processing.processing_barcoded_eel_fov_graph(
            self.experiment_fpath,
            self.analysis_parameters,
            self.running_functions,
            self.tiles_org,
            self.metadata,
            self.grpd_fovs,
            self.save_intermediate_steps,
            self.preprocessed_image_tag,
            self.client,
            self.chunk_size,
            self.save_bits_int,
            self.start_from_preprocessed_imgs,
        )

        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        selected_Hdistance = 3 / self.metadata["barcode_length"]
        stitching_selected = "microscope_stitched"
        io.simple_output_plotting(
            self.experiment_fpath,
            stitching_selected,
            selected_Hdistance,
            self.client,
            input_file_tag="decoded_fov",
            file_tag="microscope_stitched",
        )
        # ----------------------------------------------------------------

    def rerun_decoding_step(self):
        """
        Create and run a dask delayed task graph used to process barcoded eel experiments
        It runs:
        (2) Barcode decoding
        (3) Registration to the microscope coords

        The following attributes created by another step must be accessible:
        - metadata
        - analysis_parameters
        - grpd_fovs
        - client

        """
        assert self.metadata, self.logger.error(
            f"cannot process eel fovs because missing metadata attr"
        )
        assert self.analysis_parameters, self.logger.error(
            f"cannot process eel fovs because missing analysis_parameters attr"
        )
        assert self.grpd_fovs, self.logger.error(
            f"cannot process eel fovs because missing grpd_fovs attr"
        )
        assert self.client, self.logger.error(
            f"cannot process eel fovs because missing client attr"
        )
        assert self.tiles_org, self.logger.error(
            f"cannot process eel fovs because missing tiles organization attr"
        )

        fov_processing.processing_barcoded_eel_fov_graph_from_decoding(
            self.experiment_fpath,
            self.analysis_parameters,
            self.tiles_org,
            self.metadata,
            self.grpd_fovs,
            self.client,
            self.chunk_size,
        )

    def rerun_from_registration_step(self):
        """
        Create and run a dask delayed task graph from the registration step.
        Requires the raw_counts files

        It runs:
        (1) Field of view registration
        (2) Barcode decoding
        (3) Registration to the microscope coords

        The following attributes created by another step must be accessible:
        - metadata
        - analysis_parameters
        - running_functions
        - grpd_fovs
        - client
        - tiles_org

        """

        assert self.metadata, self.logger.error(
            f"cannot process eel fovs because missing metadata attr"
        )
        assert self.analysis_parameters, self.logger.error(
            f"cannot process eel fovs because missing analysis_parameters attr"
        )
        assert self.running_functions, self.logger.error(
            f"cannot process eel fovs because missing running_functions attr"
        )
        assert self.grpd_fovs, self.logger.error(
            f"cannot process eel fovs because missing grpd_fovs attr"
        )
        assert self.client, self.logger.error(
            f"cannot process eel fovs because missing client attr"
        )
        assert self.tiles_org, self.logger.error(
            f"cannot process eel fovs because missing tiles organization attr"
        )

        fov_processing.processing_barcoded_eel_fov_starting_from_registration_graph(
            self.experiment_fpath,
            self.analysis_parameters,
            self.running_functions,
            self.tiles_org,
            self.metadata,
            self.grpd_fovs,
            self.preprocessed_image_tag,
            self.client,
            self.chunk_size,
            self.save_bits_int,
        )

        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        selected_Hdistance = 3 / self.metadata["barcode_length"]
        stitching_selected = "microscope_stitched"
        io.simple_output_plotting(
            self.experiment_fpath,
            stitching_selected,
            selected_Hdistance,
            self.client,
            input_file_tag="decoded_fov",
            file_tag="microscope_stitched",
        )
        # ----------------------------------------------------------------

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
        assert self.metadata, self.logger.error(
            f"cannot process smFISH fovs because missing metadata attr"
        )
        assert self.analysis_parameters, self.logger.error(
            f"cannot process smFISH fovs because missing analysis_parameters attr"
        )
        assert self.running_functions, self.logger.error(
            f"cannot process smFISH fovs because missing running_functions attr"
        )
        assert self.grpd_fovs, self.logger.error(
            f"cannot process smFISH fovs because missing grpd_fovs attr"
        )
        assert self.client, self.logger.error(
            f"cannot process smFISH fovs because missing client attr"
        )
        assert self.tiles_org, self.logger.error(
            f"cannot process eel fovs because missing tiles organization attr"
        )

        fov_processing.processing_serial_fish_fov_graph(
            self.experiment_fpath,
            self.analysis_parameters,
            self.running_functions,
            self.tiles_org,
            self.metadata,
            self.grpd_fovs,
            self.save_intermediate_steps,
            self.preprocessed_image_tag,
            self.client,
            self.chunk_size,
            self.start_from_preprocessed_imgs,
        )

        # Removed the dots on the microscope stitched
        self.stitching_selected = "microscope_stitched"

        stitching.remove_duplicated_dots_serial_graph(
            self.experiment_fpath,
            self.data.dataset,
            self.tiles_org,
            self.same_dot_radius_duplicate_dots,
            self.stitching_selected,
            self.client,
        )

        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        stitching_selected = "microscope_stitched"
        io.simple_output_plotting_serial(
            self.experiment_fpath,
            stitching_selected,
            self.client,
            input_file_tag="microscope_stitched_cleaned",
            file_tag="cleaned_microscope_stitched",
        )
        # ----------------------------------------------------------------

    def microscope_stitched_remove_dots_eel_graph_step(self):
        """
        Function to remove the duplicated
        barcodes present in the overlapping regions of the tiles after stitching
        using the microscope coords

        Args:
        ----
        hamming_distance (int): Value to select the barcodes that are passing the
            screening (< hamming_distance). Default = 3
        same_dot_radius_duplicate_dots (int): Searching distance that define two dots as identical
            Default = 10
        stitching_selected (str): barcodes coords set where the duplicated dots will be
            removed

        The following attributes created by another step must be accessible:
        - dataset
        - tiles_org
        - client

        """
        assert self.client, self.logger.error(
            f"cannot remove duplicated dots because missing client attr"
        )
        assert isinstance(self.data.dataset, pd.DataFrame), self.logger.error(
            f"cannot remove duplicated dots because missing dataset attr"
        )

        assert isinstance(
            self.tiles_org, pysmFISH.stitching.organize_square_tiles
        ) | isinstance(
            self.tiles_org, pysmFISH.stitching.organize_square_tiles_old_room
        ), self.logger.error(
            f"cannot remove duplicated dots because tiles_org is missing attr"
        )

        # Removed the dots on the microscope stitched
        self.stitching_selected = "microscope_stitched"

        stitching.remove_duplicated_dots_graph(
            self.experiment_fpath,
            self.data.dataset,
            self.tiles_org,
            self.hamming_distance,
            self.same_dot_radius_duplicate_dots,
            self.stitching_selected,
            self.client,
        )

        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        selected_Hdistance = 3 / self.metadata["barcode_length"]
        stitching_selected = "microscope_stitched"
        io.simple_output_plotting(
            self.experiment_fpath,
            stitching_selected,
            selected_Hdistance,
            self.client,
            input_file_tag="microscope_stitched_cleaned",
            file_tag="cleaned_microscope_stitched",
        )
        # ----------------------------------------------------------------

    # def stitch_and_remove_dots_eel_graph_step(self):

    #     """
    #     Function to stitch the different fovs and remove the duplicated
    #     barcodes present in the overlapping regions of the tiles

    #     Args:
    #     ----
    #     hamming_distance (int): Value to select the barcodes that are passing the
    #         screening (< hamming_distance). Default = 3
    #     same_dot_radius_duplicate_dots (int): Searching distance that define two dots as identical
    #         Default = 10
    #     stitching_selected (str): barcodes coords set where the duplicated dots will be
    #         removed

    #     The following attributes created by another step must be accessible:
    #     - dataset
    #     - tiles_org
    #     - client

    #     """
    #     assert self.client, self.logger.error(
    #         f"cannot remove duplicated dots because missing client attr"
    #     )
    #     assert isinstance(self.data.dataset, pd.DataFrame), self.logger.error(
    #         f"cannot remove duplicated dots because missing dataset attr"
    #     )
    #     assert isinstance(
    #         self.tiles_org, pysmFISH.stitching.organize_square_tiles
    #     ), self.logger.error(
    #         f"cannot remove duplicated dots because tiles_org is missing attr"
    #     )

    #     self.adjusted_coords = stitching.stitching_graph(
    #         self.experiment_fpath,
    #         self.metadata["stitching_channel"],
    #         self.tiles_org,
    #         self.metadata,
    #         self.analysis_parameters["RegistrationReferenceHybridization"],
    #         self.client,
    #     )

    #     # Recalculate the overlapping regions after stitching
    #     self.tiles_org.tile_corners_coords_pxl = self.adjusted_coords
    #     self.tiles_org.determine_overlapping_regions()

    #     # Removed the dots on the global stitched
    #     self.stitching_selected = "global_stitched"
    #     stitching.remove_duplicated_dots_graph(
    #         self.experiment_fpath,
    #         self.data.dataset,
    #         self.tiles_org,
    #         self.hamming_distance,
    #         self.same_dot_radius_duplicate_dots,
    #         self.stitching_selected,
    #         self.client,
    #     )

    #     # ----------------------------------------------------------------
    #     # GENERATE OUTPUT FOR PLOTTING
    #     selected_Hdistance = 3 / self.metadata["barcode_length"]
    #     stitching_selected = "global_stitched"
    #     io.simple_output_plotting(
    #         self.experiment_fpath,
    #         stitching_selected,
    #         selected_Hdistance,
    #         self.client,
    #         input_file_tag="global_stitched_cleaned",
    #         file_tag="cleaned_global_stitched",
    #     )
    #     # ----------------------------------------------------------------

    #     # ----------------------------------------------------------------
    #     # GENERATE OUTPUT FOR PLOTTING
    #     selected_Hdistance = 3 / self.metadata["barcode_length"]
    #     stitching_selected = "global_stitched"
    #     io.simple_output_plotting(
    #         self.experiment_fpath,
    #         stitching_selected,
    #         selected_Hdistance,
    #         self.client,
    #         input_file_tag="global_stitched_removed",
    #         file_tag="removed_global_stitched",
    #     )
    #     # ----------------------------------------------------------------

    def microscope_stitched_remove_dots_serial_graph_step(self):
        """
        Function to remove the duplicated
        barcodes present in the overlapping regions of the tiles after stitching
        using the microscope coords

        Args:
        ----
        hamming_distance (int): Value to select the barcodes that are passing the
            screening (< hamming_distance). Default = 3
        same_dot_radius_duplicate_dots (int): Searching distance that define two dots as identical
            Default = 10
        stitching_selected (str): barcodes coords set where the duplicated dots will be
            removed

        The following attributes created by another step must be accessible:
        - dataset
        - tiles_org
        - client

        """
        assert self.client, self.logger.error(
            f"cannot remove duplicated dots because missing client attr"
        )
        assert isinstance(self.data.dataset, pd.DataFrame), self.logger.error(
            f"cannot remove duplicated dots because missing dataset attr"
        )

        assert isinstance(
            self.tiles_org, pysmFISH.stitching.organize_square_tiles
        ) | isinstance(
            self.tiles_org, pysmFISH.stitching.organize_square_tiles_old_room
        ), self.logger.error(
            f"cannot remove duplicated dots because tiles_org is missing attr"
        )

        # Removed the dots on the microscope stitched
        self.stitching_selected = "microscope_stitched"

        stitching.remove_duplicated_dots_serial_graph(
            self.experiment_fpath,
            self.data.dataset,
            self.tiles_org,
            self.same_dot_radius_duplicate_dots,
            self.stitching_selected,
            self.client,
        )

        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        stitching_selected = "microscope_stitched"
        io.simple_output_plotting_serial(
            self.experiment_fpath,
            stitching_selected,
            input_file_tag="microscope_stitched_cleaned",
            file_tag="cleaned_microscope_stitched",
        )
        # ----------------------------------------------------------------

    def stitch_and_remove_dots_eel_graph_step_global(self):

        """
        Function to stitch the different fovs and remove the duplicated
        barcodes present in the overlapping regions of the tiles
        Args:
        ----
        hamming_distance (int): Value to select the barcodes that are passing the
            screening (< hamming_distance). Default = 3
        same_dot_radius_duplicate_dots (int): Searching distance that define two dots as identical
            Default = 10
        stitching_selected (str): barcodes coords set where the duplicated dots will be
            removed
        The following attributes created by another step must be accessible:
        - dataset
        - tiles_org
        - client
        """
        assert self.client, self.logger.error(
            f"cannot remove duplicated dots because missing client attr"
        )
        assert isinstance(self.data.dataset, pd.DataFrame), self.logger.error(
            f"cannot remove duplicated dots because missing dataset attr"
        )
        assert isinstance(
            self.tiles_org, pysmFISH.stitching.organize_square_tiles
        ), self.logger.error(
            f"cannot remove duplicated dots because tiles_org is missing attr"
        )

        self.adjusted_coords = stitching.stitching_graph(
            self.experiment_fpath,
            self.metadata["stitching_channel"],
            self.tiles_org,
            self.metadata,
            self.analysis_parameters["RegistrationReferenceHybridization"],
            self.client,
        )

        # Recalculate the overlapping regions after stitching
        self.tiles_org.tile_corners_coords_pxl = self.adjusted_coords
        self.tiles_org.determine_overlapping_regions()

        # Removed the dots on the global stitched
        self.stitching_selected = "global_stitched"
        stitching.remove_duplicated_dots_graph(
            self.experiment_fpath,
            self.data.dataset,
            self.tiles_org,
            self.hamming_distance,
            self.same_dot_radius_duplicate_dots,
            self.stitching_selected,
            self.client,
        )

        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        selected_Hdistance = 3 / self.metadata["barcode_length"]
        stitching_selected = "global_stitched"
        io.simple_output_plotting(
            self.experiment_fpath,
            stitching_selected,
            selected_Hdistance,
            self.client,
            input_file_tag="global_stitched_cleaned",
            file_tag="cleaned_global_stitched",
        )
        # ----------------------------------------------------------------

        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        selected_Hdistance = 3 / self.metadata["barcode_length"]
        stitching_selected = "global_stitched"
        io.simple_output_plotting(
            self.experiment_fpath,
            stitching_selected,
            selected_Hdistance,
            self.client,
            input_file_tag="global_stitched_removed",
            file_tag="removed_global_stitched",
        )

    def stitch_and_remove_dots_serial_graph_step(self):

        """
        Function to stitch the different fovs and remove the duplicated
        barcodes present in the overlapping regions of the tiles

        Args:
        ----
        hamming_distance (int): Value to select the barcodes that are passing the
            screening (< hamming_distance). Default = 3
        same_dot_radius_duplicate_dots (int): Searching distance that define two dots as identical
            Default = 10
        stitching_selected (str): barcodes coords set where the duplicated dots will be
            removed

        The following attributes created by another step must be accessible:
        - dataset
        - tiles_org
        - client

        """
        assert self.client, self.logger.error(
            f"cannot remove duplicated dots because missing client attr"
        )
        assert isinstance(self.data.dataset, pd.DataFrame), self.logger.error(
            f"cannot remove duplicated dots because missing dataset attr"
        )
        assert isinstance(
            self.tiles_org, pysmFISH.stitching.organize_square_tiles
        ), self.logger.error(
            f"cannot remove duplicated dots because tiles_org is missing attr"
        )

        self.adjusted_coords = stitching.stitching_graph(
            self.experiment_fpath,
            self.metadata["stitching_channel"],
            self.tiles_org,
            self.metadata,
            self.analysis_parameters["RegistrationReferenceHybridization"],
            self.client,
        )

        # Recalculate the overlapping regions after stitching
        self.tiles_org.tile_corners_coords_pxl = self.adjusted_coords
        self.tiles_org.determine_overlapping_regions()

        # Removed the dots on the global stitched
        self.stitching_selected = "global_stitched"
        stitching.remove_duplicated_dots_serial_graph(
            self.experiment_fpath,
            self.data.dataset,
            self.tiles_org,
            self.same_dot_radius_duplicate_dots,
            self.stitching_selected,
            self.client,
        )

        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        stitching_selected = "global_stitched"
        io.simple_output_plotting_serial(
            self.experiment_fpath,
            stitching_selected,
            self.client,
            input_file_tag="global_stitched_cleaned",
            file_tag="cleaned_global_stitched",
        )
        # ----------------------------------------------------------------

        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        stitching_selected = "global_stitched"
        io.simple_output_plotting_serial(
            self.experiment_fpath,
            stitching_selected,
            self.client,
            input_file_tag="global_stitched_removed",
            file_tag="removed_global_stitched",
        )
        # ----------------------------------------------------------------

    def stitch_and_remove_dots_eel_graph_old_room_step(self):

        """
        Function to stitch the different fovs and remove the duplicated
        barcodes present in the overlapping regions of the tiles

        Args:
        ----
        hamming_distance (int): Value to select the barcodes that are passing the
            screening (< hamming_distance). Default = 3
        same_dot_radius_duplicate_dots (int): Searching distance that define two dots as identical
            Default = 10
        stitching_selected (str): barcodes coords set where the duplicated dots will be
            removed

        The following attributes created by another step must be accessible:
        - dataset
        - tiles_org
        - client

        """
        assert self.client, self.logger.error(
            f"cannot remove duplicated dots because missing client attr"
        )
        assert isinstance(self.data.dataset, pd.DataFrame), self.logger.error(
            f"cannot remove duplicated dots because missing dataset attr"
        )
        assert isinstance(
            self.tiles_org, pysmFISH.stitching.organize_square_tiles
        ), self.logger.error(
            f"cannot remove duplicated dots because tiles_org is missing attr"
        )

        self.adjusted_coords = stitching.stitching_graph(
            self.experiment_fpath,
            self.metadata["stitching_channel"],
            self.tiles_org,
            self.metadata,
            self.analysis_parameters["RegistrationReferenceHybridization"],
            self.client,
        )

        # Recalculate the overlapping regions after stitching
        self.tiles_org.tile_corners_coords_pxl = self.adjusted_coords
        self.tiles_org.determine_overlapping_regions()

        # Removed the dots on the global stitched
        self.stitching_selected = "global_stitched"
        stitching.remove_duplicated_dots_graph(
            self.experiment_fpath,
            self.data.dataset,
            self.tiles_org,
            self.hamming_distance,
            self.same_dot_radius_duplicate_dots,
            self.stitching_selected,
            self.client,
        )

        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        selected_Hdistance = 3 / self.metadata["barcode_length"]
        stitching_selected = "global_stitched"
        io.simple_output_plotting(
            self.experiment_fpath,
            stitching_selected,
            selected_Hdistance,
            self.client,
            input_file_tag="global_stitched_cleaned",
            file_tag="cleaned_global_stitched",
        )
        # ----------------------------------------------------------------

        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        selected_Hdistance = 3 / self.metadata["barcode_length"]
        stitching_selected = "global_stitched"
        io.simple_output_plotting(
            self.experiment_fpath,
            stitching_selected,
            selected_Hdistance,
            self.client,
            input_file_tag="global_stitched_removed",
            file_tag="removed_global_stitched",
        )
        # ----------------------------------------------------------------

    # def stitch_and_remove_dots_eel_graph_step(self):

    #     """
    #     Function to stitch the different fovs and remove the duplicated
    #     barcodes present in the overlapping regions of the tiles

    #     Args:
    #     ----
    #     hamming_distance (int): Value to select the barcodes that are passing the
    #         screening (< hamming_distance). Default = 3
    #     same_dot_radius_duplicate_dots (int): Searching distance that define two dots as identical
    #         Default = 10
    #     stitching_selected (str): barcodes coords set where the duplicated dots will be
    #         removed

    #     The following attributes created by another step must be accessible:
    #     - dataset
    #     - tiles_org
    #     - client

    #     """
    #     assert self.client, self.logger.error(
    #         f"cannot remove duplicated dots because missing client attr"
    #     )
    #     assert isinstance(self.data.dataset, pd.DataFrame), self.logger.error(
    #         f"cannot remove duplicated dots because missing dataset attr"
    #     )
    #     assert isinstance(
    #         self.tiles_org, pysmFISH.stitching.organize_square_tiles
    #     ), self.logger.error(
    #         f"cannot remove duplicated dots because tiles_org is missing attr"
    #     )

    #     self.adjusted_coords = stitching.stitching_graph_serial_nuclei(
    #         self.experiment_fpath,
    #         self.metadata["stitching_channel"],
    #         self.tiles_org,
    #         self.metadata,
    #         self.analysis_parameters["RegistrationReferenceHybridization"],
    #         self.client,
    #     )

    #     # Recalculate the overlapping regions after stitching
    #     self.tiles_org.tile_corners_coords_pxl = self.adjusted_coords
    #     self.tiles_org.determine_overlapping_regions()

    #     # Removed the dots on the global stitched
    #     self.stitching_selected = "global_stitched"
    #     stitching.remove_duplicated_dots_graph(
    #         self.experiment_fpath,
    #         self.data.dataset,
    #         self.tiles_org,
    #         self.hamming_distance,
    #         self.same_dot_radius_duplicate_dots,
    #         self.stitching_selected,
    #         self.client,
    #     )

    #     # ----------------------------------------------------------------
    #     # GENERATE OUTPUT FOR PLOTTING
    #     selected_Hdistance = 3 / self.metadata["barcode_length"]
    #     stitching_selected = "global_stitched"
    #     io.simple_output_plotting(
    #         self.experiment_fpath,
    #         stitching_selected,
    #         selected_Hdistance,
    #         self.client,
    #         input_file_tag="global_stitched_cleaned",
    #         file_tag="cleaned_global_stitched",
    #     )
    #     # ----------------------------------------------------------------

    #     # ----------------------------------------------------------------
    #     # GENERATE OUTPUT FOR PLOTTING
    #     selected_Hdistance = 3 / self.metadata["barcode_length"]
    #     stitching_selected = "global_stitched"
    #     io.simple_output_plotting(
    #         self.experiment_fpath,
    #         stitching_selected,
    #         selected_Hdistance,
    #         self.client,
    #         input_file_tag="global_stitched_removed",
    #         file_tag="removed_global_stitched",
    #     )
    #     # ----------------------------------------------------------------

    def stitch_and_remove_dots_eel_graph_step(self):

        """
        Function to stitch the different fovs and remove the duplicated
        barcodes present in the overlapping regions of the tiles


        Args:
        ----
        hamming_distance (int): Value to select the barcodes that are passing the
            screening (< hamming_distance). Default = 3
        same_dot_radius_duplicate_dots (int): Searching distance that define two dots as identical
            Default = 10
        stitching_selected (str): barcodes coords set where the duplicated dots will be
            removed

        The following attributes created by another step must be accessible:
        - dataset
        - tiles_org
        - client

        """

        assert isinstance(self.data.dataset, pd.DataFrame), self.logger.error(
            f"cannot remove duplicated dots because missing dataset attr"
        )
        assert isinstance(
            self.tiles_org, pysmFISH.stitching.organize_square_tiles
        ), self.logger.error(
            f"cannot remove duplicated dots because tiles_org is missing attr"
        )

        # Input parameters
        # Folder with parquet files from the Results folder, like: LBEXP20210718_EEL_Mouse_448_2_decoded_fov_670.parquet
        folder = os.path.join(self.experiment_fpath.as_posix(), "results")

        columns_to_load = [
            "dot_id",
            "r_px_microscope_stitched",
            "c_px_microscope_stitched",
            "r_px_original",
            "c_px_original",
            "channel",
            "hamming_distance",
            "decoded_genes",
            "round_num",
            "mapped_beads_type"
        ]

        # FOV numbers
        fovs = list(self.tiles_org.overlapping_regions.keys())

        # List of file names
        filenames = [
            os.path.join(
                folder, f"{self.metadata['experiment_name']}_decoded_fov_{i}.parquet"
            )
            for i in fovs
        ]

        # Dictionary with fov ID as keys a filename as value
        fov_filenames = dict(zip(fovs, filenames))

        # Load all dataframes in dictionary
        fov_df = {
            k: FOV_alignment.load_parquet(fov_filenames[k], columns_to_load)
            for k in fovs
        }

        selected_Hdistance = self.hamming_distance / self.metadata["barcode_length"]

        # Run cleaning
        fname_rna_merge, fname_large_beads = FOV_alignment.clean_microscope_stitched2(
            fovs,
            fov_df,
            self.tiles_org.overlapping_regions,
            self.tiles_org.tile_corners_coords_pxl,
            mode=self.fov_alignment_mode,
            bead_channel=self.metadata["stitching_channel"],
            max_hamming_dist=selected_Hdistance,
            matching_dot_radius=self.same_dot_radius_duplicate_dots,
            out_folder=folder,
            exp_name=self.metadata["experiment_name"],
            remove_distinct_genes=self.remove_distinct_genes,
            clip_size=self.clip_size,
            verbose=False,
        )  # Set to False in pipeline

        self.metadata["fname_rna_merge"] = fname_rna_merge
        self.metadata["fname_large_beads"] = fname_large_beads

    def stitch_and_remove_dots_eel_graph_step_bayesian(self):

        """
        Function to stitch the different fovs and remove the duplicated
        barcodes present in the overlapping regions of the tiles


        Args:
        ----
        hamming_distance (int): Value to select the barcodes that are passing the
            screening (< hamming_distance). Default = 3
        same_dot_radius_duplicate_dots (int): Searching distance that define two dots as identical
            Default = 10
        stitching_selected (str): barcodes coords set where the duplicated dots will be
            removed

        The following attributes created by another step must be accessible:
        - dataset
        - tiles_org
        - client

        """

        assert isinstance(self.data.dataset, pd.DataFrame), self.logger.error(
            f"cannot remove duplicated dots because missing dataset attr"
        )
        assert isinstance(
            self.tiles_org, pysmFISH.stitching.organize_square_tiles
        ), self.logger.error(
            f"cannot remove duplicated dots because tiles_org is missing attr"
        )

        # Input parameters
        # Folder with parquet files from the Results folder, like: LBEXP20210718_EEL_Mouse_448_2_decoded_fov_670.parquet
        folder = os.path.join(self.experiment_fpath.as_posix(), "results")

        columns_to_load = [
            "dot_id",
            "r_px_microscope_stitched",
            "c_px_microscope_stitched",
            "r_px_original",
            "c_px_original",
            "channel",
            "hamming_distance",
            "decoded_genes",
            "round_num",
            "raw_barcodes",
            "mapped_beads_type",
            "bit_1_intensity",	
            "bit_2_intensity",	
            "bit_3_intensity",
            "bit_4_intensity",
            "bit_5_intensity",
            "bit_6_intensity",
            "bit_7_intensity",
            "bit_8_intensity",
            "bit_9_intensity",
            "bit_10_intensity",	
            "bit_11_intensity",	
            "bit_12_intensity",	
            "bit_13_intensity",
            "bit_14_intensity",
            "bit_15_intensity",
            "bit_16_intensity",
        ]


        # FOV numbers
        fovs = list(self.tiles_org.overlapping_regions.keys())

        # List of file names
        filenames = [
            os.path.join(
                folder, f"{self.metadata['experiment_name']}_decoded_fov_{i}.parquet"
            )
            for i in fovs
        ]

        # Dictionary with fov ID as keys a filename as value
        fov_filenames = dict(zip(fovs, filenames))

        # Load all dataframes in dictionary
        fov_df = {
        }
       
        for k in fovs:
            try:
                fov_df[k] = pd.read_parquet(fov_filenames[k],columns=columns_to_load)
            except:
                print('Not existing FOV {}, filling empty data.'.format(k))
                fov_df[k] = pd.DataFrame(columns=columns_to_load)
        #print(fovs)
        
        selected_Hdistance = self.hamming_distance / self.metadata["barcode_length"]

        # Run cleaning
        fname_rna_merge, fname_large_beads = FOV_alignment.clean_microscope_stitched(
            #[x for x in range(2,4)],#
            fovs,
            fov_df,
            self.tiles_org.overlapping_regions,
            self.tiles_org.tile_corners_coords_pxl,
            mode=self.fov_alignment_mode,
            bead_channel=self.metadata["stitching_channel"],
            max_hamming_dist=selected_Hdistance,
            matching_dot_radius=self.same_dot_radius_duplicate_dots,
            out_folder=folder,
            exp_name=self.metadata["experiment_name"],
            remove_distinct_genes=self.remove_distinct_genes,
            clip_size=self.clip_size,
            verbose=False,
        )  # Set to False in pipeline

        self.metadata["fname_rna_merge"] = fname_rna_merge
        self.metadata["fname_large_beads"] = fname_large_beads

        df = pd.read_parquet(fname_rna_merge)
        df = self.process_bayesian(df)
        df.to_parquet(fname_rna_merge.split('.')[0]+'_bayesian.parquet' )

        df = df[df['probability'] >= 0.5]
        df = df[~df.decoded_genes.isna()]
        df.loc[
            :,
            [
            "r_px_microscope_stitched",
            "c_px_microscope_stitched",
            "hamming_distance",
            "decoded_genes",
            "channel",
            "round_num",
            "probability",
        ]].to_parquet(fname_rna_merge)

    def process_bayesian(self, df):
        from sklearn.preprocessing import normalize
        from pysmFISH.bayeseel import BayesEEL
        channels = df.channel.unique()

        dfs = []
        for ch in channels:
            print('Running bayesian calling on: ',ch)
            codebook_filename = self.metadata['list_all_codebooks'][np.where(self.metadata['list_all_channels']==ch)[0]][0]
            barcodes = pd.read_parquet(self.experiment_fpath/'codebook'/codebook_filename)
            df_i= df[df.channel == ch]
            in_decoded = df_i.loc[:,["bit_1_intensity","bit_2_intensity","bit_3_intensity","bit_4_intensity","bit_5_intensity","bit_6_intensity","bit_7_intensity","bit_8_intensity",
                "bit_9_intensity","bit_10_intensity","bit_11_intensity","bit_12_intensity",	"bit_13_intensity","bit_14_intensity","bit_15_intensity","bit_16_intensity",]].fillna(0)


            ids = (df_i.hamming_distance <= 2/16)
            sel = df_i[ids]
            known_barcodes = []
            for x in sel.raw_barcodes:
                known_barcodes.append(np.array([True if y == 1 else False for y in x]))
            known_barcodes = np.array(known_barcodes)
            in_decoded = in_decoded[ids].values

            print('Training BayesEEL: ', ch)

            BE = BayesEEL()
            BE_ = BE.fit(known_x=normalize(in_decoded),known_barcodes=known_barcodes)
            print('Transforming BayesEEL: ', ch)
            known_barcodes = []
            for x in df_i.raw_barcodes:
                known_barcodes.append(np.array([True if y == 1 else False for y in x]))
            known_barcodes = np.array(known_barcodes)
            
            in_decoded = df_i.loc[:,["bit_1_intensity","bit_2_intensity","bit_3_intensity","bit_4_intensity","bit_5_intensity","bit_6_intensity","bit_7_intensity","bit_8_intensity",
                "bit_9_intensity","bit_10_intensity","bit_11_intensity","bit_12_intensity",	"bit_13_intensity","bit_14_intensity","bit_15_intensity","bit_16_intensity",]].fillna(0)

            input_transform = in_decoded.values
            bs = []
            for b in barcodes.Code:
                bs.append(np.array([True if y == 1 else False for y in b]))
            bs = np.array(bs)

            idx, probabilities = BE_.transform(normalize(input_transform),bs)
            decoded_genes = np.array([barcodes.Gene.values[x] for x in idx])
            df_i['decoded_genes'] = decoded_genes
            print('All dots: ',df_i.shape[0])

            df_i['probability'] = probabilities
            dfs.append(df_i)
    
        return pd.concat(dfs)


    def processing_fresh_tissue_step(
        self,
        parsing=True,
        reprocessing=True,
        segment=True,
        tag_ref_beads="_ChannelEuropium_Cy3_",
        tag_nuclei="_ChannelCy3_",
    ):
        """
        This function create and run a processing graph that parse and filter the nuclei staining in fresh tissue
        and parse, filter and counts the Europium beads used for the registration with the smFISH images at high
        power.

        Args:
        ----
        tag_ref_beads (str): The tag reference of the .nd2 file containing the raw beads images. Default: '_ChannelEuropium_Cy3_'
        tag_ref_nuclei (str): The tag reference of the .nd2 file containing the raw images of nuclei. Default: '_ChannelCy3_'

        """
        assert self.client, self.logger.error(
            f"cannot process fresh tissue because missing client attr"
        )
        assert self.analysis_parameters, self.logger.error(
            f"cannot process fresh tissue because missing analysis_parameters attr"
        )
        assert self.running_functions, self.logger.error(
            f"cannot process fresh tissue because missing running_functions attr"
        )
        if not os.path.exists(Path(self.experiment_fpath)/ "fresh_tissue"/ "segmentation"):
            os.makedirs(Path(self.experiment_fpath)/ "fresh_tissue"/ "segmentation")
        if reprocessing:
            (
                self.ds_beads,
                self.ds_nuclei,
                self.nuclei_metadata,
            ) = fov_processing.process_fresh_sample_graph(
                self.experiment_fpath,
                self.running_functions,
                self.analysis_parameters,
                self.client,
                self.chunk_size,
                tag_ref_beads=tag_ref_beads,
                tag_nuclei=tag_nuclei,
                eel_metadata=self.metadata,
                fresh_tissue_segmentation_engine=self.fresh_tissue_segmentation_engine,
                diameter_size=self.diameter_size,
                parsing=parsing,
                save_steps_output=self.save_intermediate_steps,
            )

            pickle.dump(
                [
                    self.ds_beads,
                    self.ds_nuclei,
                    self.metadata,
                ],
                open(
                    Path(self.experiment_fpath)
                    / "fresh_tissue"
                    / "segmentation"
                    / "ds_tmp_data.pkl",
                    "wb",
                ),
            )

        (self.ds_beads, self.ds_nuclei, self.nuclei_metadata) = pickle.load(
             open(
                 Path(self.experiment_fpath)
                 / "fresh_tissue"
                 / "segmentation"
                 / "ds_tmp_data.pkl",
                 "rb",
             ),
        )
        # Segmentation

        if segment:
            fov_processing.segmentation_graph(
                self.ds_nuclei,
                self.chunk_size,
                self.experiment_fpath,
                self.fresh_tissue_segmentation_engine,
                self.diameter_size,
            )

        (
            self.nuclei_org_tiles,
            self.nuclei_adjusted_coords,
        ) = stitching.stitched_beads_on_nuclei_fresh_tissue(
            self.experiment_fpath,
            self.client,
            self.ds_nuclei,
            self.ds_beads,
            round_num=1,
        )

        pickle.dump(
            [
                self.ds_beads,
                self.ds_nuclei,
                self.nuclei_metadata,
                self.nuclei_org_tiles,
                self.nuclei_adjusted_coords,
            ],
            open(
                Path(self.experiment_fpath)
                / "fresh_tissue"
                / "segmentation"
                / "tmp_data.pkl",
                "wb",
            ),
        )

        (
            self.ds_beads,
            self.ds_nuclei,
            self.nuclei_metadata,
            self.nuclei_org_tiles,
            self.nuclei_adjusted_coords,
        ) = pickle.load(
            open(
                Path(self.experiment_fpath)
                / "fresh_tissue"
                / "segmentation"
                / "tmp_data.pkl",
                "rb",
            ),
        )
        
        segmentation_output_path = (
            Path(self.experiment_fpath) / "fresh_tissue" / "segmentation"
        )
        segmented_object_dict_recalculated = segmentation.create_label_image(
            self.experiment_fpath,
            segmentation_output_path,
            self.ds_nuclei,
            self.nuclei_org_tiles,
            self.nuclei_adjusted_coords,
            self.nuclei_metadata,
            self.client,
            self.min_overlapping_pixels_segmentation,
        )
        gc.collect()
    
    def processing_assign_dots(self):

        segmentation_output_path = (
            Path(self.experiment_fpath) / "fresh_tissue" / "segmentation"
        )
        
        (
            self.ds_beads,
            self.ds_nuclei,
            self.nuclei_metadata,
            self.nuclei_org_tiles,
            self.nuclei_adjusted_coords,
        ) = pickle.load(
            open(
                Path(self.experiment_fpath)
                / "fresh_tissue"
                / "segmentation"
                / "tmp_data.pkl",
                "rb",
            ),
        )

        segmented_object_dict_recalculated = pickle.load(
        open(
            segmentation_output_path / ("segmented_objects_dict_recalculated_ids.pkl"),
            "rb",
        ),
        )
        segmentation.register_assign(
            self.experiment_fpath,
            segmented_object_dict_recalculated,
            self.data.dataset,
            self.ds_nuclei.dataset,
            self.metadata,
            self.nuclei_metadata,
            self.pipeline_run_name,
            segmentation_output_path,
            self.max_expansion_radius,
            self.hamming_distance,
            self.bead_alignment_centering_mode,
            self.bead_alignment_radius
        )

    # --------------------------------
    # QC STEPS (some other included in the graph function)
    # --------------------------------

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
        assert self.client, self.logger.error(
            f"cannot run QC on registration because missing client attr"
        )
        assert self.analysis_parameters, self.logger.error(
            f"cannot run QC on registration because missing analysis_parameters attr"
        )
        assert self.metadata, self.logger.error(
            f"cannot process smFISH fovs because missing metadata attr"
        )
        assert isinstance(self.tile_corners_coords_pxl, np.ndarray), self.logger.error(
            f"cannot run QC on registration because missing tile_corners_coords_pxl attr"
        )

        qc_reg = qc_utils.QC_registration_error(
            self.client,
            self.experiment_fpath,
            self.analysis_parameters,
            self.metadata,
            self.tile_corners_coords_pxl,
        )

        qc_reg.run_qc()

    # --------------------------------
    # DATA REORGANIZATION
    # --------------------------------

    def transfer_data_after_processing(self):
        """
        Function use to clear space in the processing folder.
        - Tranfer parsed images zarr / filtered images zarr / dataset
        - Remove the log and tmp folders
        - Transfer the remaining data to cold storage
        Before transfering the data it is necessary to close the cluster and the client
        otherwise you won't be able to remove the logs folder.
        """

        data_organization.reorganize_processing_dir(
            self.experiment_fpath,
            self.storage_fpath,
            self.store_dataset,
            self.dataset_folder_storage_path,
            self.results_folder_storage_path,
        )

    # --------------------------------
    # RUNNING OPTIONS
    # --------------------------------

    def run_setup(self):
        start = datetime.now()

        self.create_folders_step()

        self.logger = logger_utils.json_logger(
            (self.experiment_fpath / "logs"), "pipeline_run"
        )
        self.logger.info(f"Start parsing")

        self.logger.info(f"Saved current git commit version")

        if self.run_type == "original":
            self.QC_check_experiment_yaml_file_step()
            self.logger.info(f"Checked config file")

        self.logger.info(
            f"{self.experiment_fpath.stem} timing: \
                Setup completed in {utils.nice_deltastring(datetime.now() - start)}."
        )

    def run_cluster_activation(self):
        start = datetime.now()
        self.processing_cluster_init_step()
        self.logger.info(f"Started dask processing cluster")
        self.logger.info(
            f"client dashboard {self.client.scheduler_info()['services']['dashboard']}"
        )
        self.logger.info(
            f"{self.experiment_fpath.stem} timing: \
                Cluester activation completed in {utils.nice_deltastring(datetime.now() - start)}."
        )

    def run_parsing(self):
        start = datetime.now()
        # Run parsing only if required
        self.logger.info(f"Parsing started")
        if self.parsing_type != "no_parsing":
            self.nikon_nd2_parsing_graph_step()
        self.logger.info(
            f"{self.experiment_fpath.stem} timing: \
                Parsing completed in {utils.nice_deltastring(datetime.now() - start)}."
        )

    def run_parsing_only(self):
        """
        Pipeline running the data organization and the parsing
        of the .nd2 files of the entire experiment
        """

        start = datetime.now()

        self.create_folders_step()

        self.logger = logger_utils.json_logger(
            (self.experiment_fpath / "logs"), "pipeline_run"
        )
        self.logger.info(f"Start parsing")

        self.logger.info(f"Saved current git commit version")

        self.QC_check_experiment_yaml_file_step()
        self.logger.info(f"Checked config file")

        self.processing_cluster_init_step()
        self.logger.info(f"Started dask processing cluster")
        self.logger.info(
            f"client dashboard {self.client.scheduler_info()['services']['dashboard']}"
        )

        # Run parsing only if required
        self.logger.info(f"Parsing started")
        if self.parsing_type != "no_parsing":
            self.nikon_nd2_parsing_graph_step()
        self.logger.info(
            f"{self.experiment_fpath.stem} timing: \
                Parsing completed in {utils.nice_deltastring(datetime.now() - start)}."
        )

        step_start = datetime.now()
        self.logger.info(f"Started creation of the dataset")
        self.prepare_processing_dataset_step()
        self.logger.info(
            f"{self.experiment_fpath.stem} timing:\
                Dataset creation completed in {utils.nice_deltastring(datetime.now() - step_start)}."
        )

    def run_required_steps(self):
        """
        Short pipeline used to make sure that the basic required
        step are run and will be included in more complex pipelines

        """
        start = datetime.now()
        step_start = datetime.now()
        self.logger.info(f"Started creation of the dataset")
        self.prepare_processing_dataset_step()
        self.logger.info(
            f"{self.experiment_fpath.stem} timing:\
                Dataset creation completed in {utils.nice_deltastring(datetime.now() - step_start)}."
        )
        self.create_analysis_config_file_from_dataset_step()
        self.logger.info(f"Created analysis_config.yaml file")
        self.determine_tiles_organization()
        self.logger.info(f"Determined the tiles organization")
        self.create_running_functions_step()
        self.logger.info(f"Created the running function dictionary")
        self.logger.info(
            f"{self.experiment_fpath.stem} timing:\
                Required steps completed in {utils.nice_deltastring(datetime.now() - start)}."
        )
        self.logger.info(f"")

    def run_full(self):
        """
        Full run from raw images from nikon or parsed images
        """
        start = datetime.now()
        self.run_setup()
        self.run_cluster_activation()
        self.run_parsing()
        self.run_required_steps()

        if self.resume:
            already_processed = (Path(self.experiment_fpath) / "results").glob(
                "*barcodes_max_array*.parquet"
            )
            already_done_fovs = []
            for fname in already_processed:
                fov_num = int(fname.stem.split("_")[-1])
                already_done_fovs.append(fov_num)
            not_processed_fovs = set(self.grpd_fovs.groups.keys()).difference(
                set(already_done_fovs)
            )
            self.data.dataset = self.data.dataset.loc[
                self.data.dataset.fov_num.isin(not_processed_fovs), :
            ]
            self.grpd_fovs = self.data.dataset.groupby("fov_num")

        if self.metadata["experiment_type"] == "eel-barcoded":
            step_start = datetime.now()
            self.processing_barcoded_eel_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    eel fov processing completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            step_start = datetime.now()
            self.QC_registration_error_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    QC registration completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            step_start = datetime.now()
            self.microscope_stitched_remove_dots_eel_graph_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    removal overlapping dots in microscope stitched {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            step_start = datetime.now()
            try:
                self.stitch_and_remove_dots_eel_graph_step()
            except:
                self.logger.info(f"Stitching using dots didn't work")
                pass
            else:
                self.logger.info(
                    f"{self.experiment_fpath.stem} timing: \
                    Stitching and removal of duplicated dots completed in {utils.nice_deltastring(datetime.now() - step_start)}."
                )

            step_start = datetime.now()
            self.processing_fresh_tissue_step()
            self.processing_assign_dots()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    Processing fresh tissue completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )

        elif self.metadata["experiment_type"] == "smfish-serial":
            step_start = datetime.now()
            self.processing_serial_fish_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    serial smfish fov processing completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )
        else:
            self.logger.error(
                f"the experiment type {self.metadata['experiment_type']} is unknown"
            )
            sys.exit(
                f"the experiment type {self.metadata['experiment_type']} is unknown"
            )

        step_start = datetime.now()
        # self.transfer_data_after_processing()
        # self.logger.info(f"{self.experiment_fpath.stem} timing: \
        #             data transfer after processing completed in {utils.nice_deltastring(datetime.now() - step_start)}.")

        self.logger.info(
            f"{self.experiment_fpath.stem} timing: \
                    Pipeline run completed in {utils.nice_deltastring(datetime.now() - start)}."
        )

        self.client.close()
        self.cluster.close()
        if self.processing_engine == "unmanaged_cluster":
            processing_cluster_setup.kill_process()

    def test_run_after_editing(self):
        """
        Full run from raw images from nikon or parsed images
        """

        start = datetime.now()
        if self.resume:
            already_processed = (Path(self.experiment_fpath) / "results").glob(
                "*decoded*.parquet"
            )
            already_done_fovs = []
            for fname in already_processed:
                fov_num = int(fname.stem.split("_")[-1])
                already_done_fovs.append(fov_num)
            not_processed_fovs = set(self.grpd_fovs.groups.keys()).difference(
                set(already_done_fovs)
            )
            self.data.dataset = self.data.dataset.loc[
                self.data.dataset.fov_num.isin(not_processed_fovs), :
            ]
            self.grpd_fovs = self.data.dataset.groupby("fov_num")

        if self.metadata["experiment_type"] == "eel-barcoded":
            step_start = datetime.now()
            self.processing_barcoded_eel_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    eel fov processing completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            # step_start = datetime.now()
            # self.QC_registration_error_step()
            # self.logger.info(f"{self.experiment_fpath.stem} timing: \
            #         QC registration completed in {utils.nice_deltastring(datetime.now() - step_start)}.")

        self.logger.info(
            f"{self.experiment_fpath.stem} timing: \
                    Pipeline run completed in {utils.nice_deltastring(datetime.now() - start)}."
        )

        # self.client.close()
        # self.cluster.close()

    def test_run_short(self):
        start = datetime.now()
        self.run_parsing_only()
        self.run_required_steps()
        if self.resume:
            already_processed = (Path(self.experiment_fpath) / "results").glob(
                "*decoded*.parquet"
            )
            already_done_fovs = []
            for fname in already_processed:
                fov_num = int(fname.stem.split("_")[-1])
                already_done_fovs.append(fov_num)
            not_processed_fovs = set(self.grpd_fovs.groups.keys()).difference(
                set(already_done_fovs)
            )
            self.data.dataset = self.data.dataset.loc[
                self.data.dataset.fov_num.isin(not_processed_fovs), :
            ]
            self.grpd_fovs = self.data.dataset.groupby("fov_num")

        if self.metadata["experiment_type"] == "eel-barcoded":
            step_start = datetime.now()
            self.processing_barcoded_eel_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    eel fov processing completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            step_start = datetime.now()
            self.QC_registration_error_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    QC registration completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )

        self.logger.info(
            f"{self.experiment_fpath.stem} timing: \
                    Pipeline run completed in {utils.nice_deltastring(datetime.now() - start)}."
        )

        self.client.close()
        self.cluster.close()

    def test_run_decoding(self):
        """
        Run analysis starting from the raw data files.
        Requires raw files
        """

        raw_files_path = list(
            (self.experiment_fpath / "results").glob("*_decoded_fov_*")
        )

        start = datetime.now()
        self.run_setup()
        self.run_cluster_activation()
        self.run_parsing()
        self.run_required_steps()
        if raw_files_path:

            step_start = datetime.now()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    eel fov processing from dots completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            step_start = datetime.now()
            self.rerun_decoding_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    eel fov processing from dots completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    Pipeline run completed in {utils.nice_deltastring(datetime.now() - start)}."
            )
        else:
            self.logger.error(
                f"{self.experiment_fpath.stem} missing files with raw counts"
            )
            raise FileNotFoundError

        self.client.close()
        self.cluster.close()

    def test_run_from_registration(self):
        """
        Run analysis starting from the raw data files.
        Requires raw files
        """

        raw_files_path = list(
            (self.experiment_fpath / "results").glob("*_raw_counts_*")
        )
        start = datetime.now()
        start = datetime.now()
        self.run_setup()
        self.run_cluster_activation()
        self.run_parsing()
        self.run_required_steps()

        if raw_files_path:

            if self.resume:
                already_processed = (Path(self.experiment_fpath) / "results").glob(
                    "*decoded*.parquet"
                )
                already_done_fovs = []
                for fname in already_processed:
                    fov_num = int(fname.stem.split("_")[-1])
                    already_done_fovs.append(fov_num)
                not_processed_fovs = set(self.grpd_fovs.groups.keys()).difference(
                    set(already_done_fovs)
                )
                self.data.dataset = self.data.dataset.loc[
                    self.data.dataset.fov_num.isin(not_processed_fovs), :
                ]
                self.grpd_fovs = self.data.dataset.groupby("fov_num")

            step_start = datetime.now()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    eel fov processing from dots completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            step_start = datetime.now()
            self.rerun_from_registration_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    eel fov processing from dots completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            step_start = datetime.now()
            self.QC_registration_error_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    QC registration completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            step_start = datetime.now()
            self.microscope_stitched_remove_dots_eel_graph_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    removal overlapping dots in microscope stitched {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            step_start = datetime.now()
            try:
                self.stitch_and_remove_dots_eel_graph_step()
            except:
                self.logger.info(f"Stitching using dots didn't work")
                pass
            else:
                self.logger.info(
                    f"{self.experiment_fpath.stem} timing: \
                    Stitching and removal of duplicated dots completed in {utils.nice_deltastring(datetime.now() - step_start)}."
                )

            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    Pipeline run completed in {utils.nice_deltastring(datetime.now() - start)}."
            )
        else:
            self.logger.error(
                f"{self.experiment_fpath.stem} missing files with raw counts"
            )
            raise FileNotFoundError

        self.client.close()
        self.cluster.close()

    def test_run_Lars_mouse_atlas(self):
        """
        The data were acquired before the rearrangement of the room and
        the fresh tissue data combined the DAPI/Europium in one single
        image.
        TODO: The fresh tissue processing need to be specifically implemented
            for this experiment

        """

        start = datetime.now()
        self.run_setup()
        self.run_cluster_activation()
        self.run_parsing()
        self.run_required_steps()
        self.determine_tiles_organization_before_room_reorganisation()

        if self.resume:
            already_processed = (Path(self.experiment_fpath) / "results").glob(
                "*barcodes_max_array*.parquet"
            )
            already_done_fovs = []
            for fname in already_processed:
                fov_num = int(fname.stem.split("_")[-1])
                already_done_fovs.append(fov_num)
            not_processed_fovs = set(self.grpd_fovs.groups.keys()).difference(
                set(already_done_fovs)
            )
            self.data.dataset = self.data.dataset.loc[
                self.data.dataset.fov_num.isin(not_processed_fovs), :
            ]
            self.grpd_fovs = self.data.dataset.groupby("fov_num")

        # Adjust the dataset because changes in the yaml files
        self.data.dataset.loc[
            self.data.dataset.channel == "Cy5", "processing_type"
        ] = "fish"
        self.data.dataset.loc[
            self.data.dataset.channel == "Europium", "processing_type"
        ] = "both-beads"

        if self.metadata["experiment_type"] == "eel-barcoded":
            step_start = datetime.now()
            self.processing_barcoded_eel_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    eel fov processing completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            step_start = datetime.now()
            self.QC_registration_error_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    QC registration completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            step_start = datetime.now()
            self.microscope_stitched_remove_dots_eel_graph_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                    removal overlapping dots in microscope stitched {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            step_start = datetime.now()
            try:
                self.stitch_and_remove_dots_eel_graph_step()
            except:
                self.logger.info(f"Stitching using dots didn't work")
                pass
            else:
                self.logger.info(
                    f"{self.experiment_fpath.stem} timing: \
                    Stitching and removal of duplicated dots completed in {utils.nice_deltastring(datetime.now() - step_start)}."
                )

        else:
            self.logger.error(
                f"the experiment type {self.metadata['experiment_type']} is unknown"
            )
            sys.exit(
                f"the experiment type {self.metadata['experiment_type']} is unknown"
            )

        step_start = datetime.now()
        # self.transfer_data_after_processing()
        # self.logger.info(f"{self.experiment_fpath.stem} timing: \
        #             data transfer after processing completed in {utils.nice_deltastring(datetime.now() - step_start)}.")

        self.logger.info(
            f"{self.experiment_fpath.stem} timing: \
                    Pipeline run completed in {utils.nice_deltastring(datetime.now() - start)}."
        )

        self.client.close()
        self.cluster.close()
        if self.processing_engine == "unmanaged_cluster":
            processing_cluster_setup.kill_process()

    def test_run_Alejandro_spinal_cord_atlas(self):

        """
        TODO: The fresh tissue processing need to be specifically implemented
            for this experiment

        """

        start = datetime.now()
        self.run_setup()
        self.run_cluster_activation()
        self.run_parsing()
        self.run_required_steps()
        self.determine_tiles_organization_before_room_reorganisation()

        if self.resume:
            already_processed = (Path(self.experiment_fpath) / "results").glob(
                "*barcodes_max_array*.parquet"
            )
            already_done_fovs = []
            for fname in already_processed:
                fov_num = int(fname.stem.split("_")[-1])
                already_done_fovs.append(fov_num)
            not_processed_fovs = set(self.grpd_fovs.groups.keys()).difference(
                set(already_done_fovs)
            )
            self.data.dataset = self.data.dataset.loc[
                self.data.dataset.fov_num.isin(not_processed_fovs), :
            ]
            self.grpd_fovs = self.data.dataset.groupby("fov_num")

        # Adjust the dataset because changes in the yaml files
        self.data.dataset.loc[:, "pipeline"] = "smfish-serial"
        self.data.dataset.loc[:, "experiment_type"] = "smfish-serial"
        self.data.dataset.loc[:, "stitching_type"] = "nuclei"
        self.data.dataset.loc[:, "stitching_channel"] = "DAPI"
        self.data.dataset.loc[
            self.data.dataset.channel == "DAPI", "processing_type"
        ] = "nuclei"
        self.data.dataset = self.data.dataset.loc[
            ~(self.data.dataset.channel == "Europium"), :
        ]
        self.data.dataset.loc[:, "machine"] = "ROBOFISH1"

        if self.metadata["experiment_type"] == "smfish-serial":
            step_start = datetime.now()
            self.processing_serial_fish_step()
            self.logger.info(
                f"{self.experiment_fpath.stem} timing: \
                        eel fov processing completed in {utils.nice_deltastring(datetime.now() - step_start)}."
            )

            step_start = datetime.now()
            try:
                self.stitch_and_remove_dots_serial_graph_step()
            except:
                self.logger.info(f"Stitching using dots didn't work")
                pass
            else:
                self.logger.info(
                    f"{self.experiment_fpath.stem} timing: \
                        Stitching and removal of duplicated dots completed in {utils.nice_deltastring(datetime.now() - step_start)}."
                )

        else:
            self.logger.error(
                f"the experiment type {self.metadata['experiment_type']} is unknown"
            )
            sys.exit(
                f"the experiment type {self.metadata['experiment_type']} is unknown"
            )

        self.client.close()
        self.cluster.close()
        if self.processing_engine == "unmanaged_cluster":
            processing_cluster_setup.kill_process()
