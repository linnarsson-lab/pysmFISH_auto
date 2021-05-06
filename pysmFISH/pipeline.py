from typing import *

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
    """
    Generale pipeline for running barcoded eel or serial smFISH experiments



    """

    def __init__(self, pipeline_run_name:str, experiment_fpath:str,
                run_type:str = 'new', parsing_type:str = 'original', **kwarg):

        self.logger = logger_utils.selected_logger()   
        self.pipeline_run_name = pipeline_run_name
        self.experiment_fpath = Path(experiment_fpath)
        self.run_type = run_type
        self.parsing_type = parsing_type

        # Collect some of the parameters. If missing a predefined value is assigned
        self.raw_data_folder_storage_path = kwarg.pop('raw_data_folder_storage_path', '/fish/rawdata')
        self.parsed_image_tag = kwarg.pop('parsed_image_tag','img_data')
        self.dataset_folder_storage_path = kwarg.pop('dataset_folder_storage_path','/fish/fish_datasets')
        self.save_intermediate_steps = kwarg.pop('save_intermediate_steps',False)
        self.store_dataset = kwarg.pop('store_dataset',True)

        # Parameters for processing in htcondor
        self.processing_env_config = {}
        self.processing_env_config['engine'] = kwarg.pop('engine','htcondor')
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

    def load_config_files(self):
        self.analysis_parameters = configuration_files.load_analysis_config_file(self.experiment_fpath)
        self.parsed_raw_data_fpath = self.experiment_fpath / (self.experiment_fpath.stem + '_' + self.parsed_image_tag + '.zarr') 
    
    def processing_cluster_init(self):
        # Start processing environment
        self.logger.info(f'Starting processing environment')
        self.cluster = processing_cluster_setup.start_processing_env(self.processing_env_config)
        self.client = Client(self.cluster)
        self.logger.debug(f'Dask dashboard info {self.client.scheduler_info()['services']}')


    def nikon_nd2_parsing_step(self):
        """
        explain all the steps that are run by this function
        """
        experiment_info = configuration_files.load_experiment_config_file(self.experiment_fpath)
        configuration_files.create_specific_analysis_config_file(self.experiment_fpath, experiment_info)
        
        # Create empty zarr file for the parse data
        parsed_raw_data_fpath = io.create_empty_zarr_file(experiment_fpath=self.experiment_fpath,
                                            tag=self.parsed_image_tag)
        if parsing_type == 'original':
            utils.collect_processing_files(self.experiment_fpath, experiment_info)
            utils.sort_data_into_folders(self.experiment_fpath, experiment_info)
            all_raw_nd2 = microscopy_file_parsers.nd2_raw_files_selector(self.experiment_fpath)

            parsing_futures = client.map(microscopy_file_parsers.nikon_nd2_autoparser_zarr,
                                    all_raw_nd2,
                                    parsed_raw_data_fpath=self.parsed_raw_data_fpath,
                                    experiment_info=experiment_info)

            # wait(parsing_futures)
            _ = self.client.gather(parsing_futures)
        
        else:
            # add error if not correct parsing type
            if parsing_type == 'reparsing_from_processing_folder':
                raw_files_fpath = self.experiment_fpath + '/raw_data'
                logger.info(f'raw_files_fpath {raw_files_fpath}')
            elif parsing_type == 'reparsing_from_storage':
                data_organization.microscopy_file_parsers.microscopy_file_parsers.nd2_raw_files_selector(self.storage_experiment_fpath, self.experiment_fpath)
                raw_files_fpath = self.storage_experiment_fpath + '/raw_data'
            
            all_raw_nd2 = microscopy_file_parsers.nd2_raw_files_selector_general(folder_fpath=self.raw_files_fpath)
            parsing_futures = client.map(microscopy_file_parsers.nikon_nd2_reparser_zarr,
                                    all_raw_nd2,
                                    parsed_raw_data_fpath=self.parsed_raw_data_fpath,
                                    experiment_info=experiment_info)

        _ = self.client.gather(parsing_futures)
        consolidated_grp = io.consolidate_zarr_metadata(self.parsed_raw_data_fpath)


    def prepare_processing_dataset(self):
        self.ds = data_models.Dataset()
        if self.dataset_path:
            self.ds.load_dataset(self.dataset_path)
        else:
            self.ds.create_full_dataset_from_zmetadata(self.parsed_raw_data_fpath)
        
        self.metadata = self.ds.collect_metadata(ds.dataset)
        self.grpd_fovs = self.ds.dataset.groupby('fov_num')


    def determine_tiles_organization(self):
        self.reference_round = analysis_parameters['RegistrationReferenceHybridization']
        self.tiles_org = stitching.organize_square_tiles(self.experiment_fpath,self.metadata,
                                    self.reference_round)
        self.tiles_org.run_tiles_organization()
        self.tile_corners_coords_pxl = self.tiles_org.tile_corners_coords_pxl


    def processing_barcoded_eel_fov(self):
        """ 
        This method create a processing graph XXXXXXXXX

        """
        dark_img = preprocessing.load_dark_image(self.experiment_fpath)
        
        # did this conversion to avoid to pass self to dask
        analysis_parameters = self.analysis_parameters
        running_functions = self.running_functions
        tile_corners_coords_pxl = self.tile_corners_coords_pxl
        
        
        dark_img = delayed(dark_img)
        analysis_parameters = delayed(analysis_parameters)
        running_functions = delayed(running_functions)
        tile_corners_coords_pxl = delayed(tile_corners_coords_pxl)

        codebook = configuration_files.load_codebook(experiment_fpath,metadata)
        codebook_df = delayed(codebook)
        
        all_processing = []
        
        for fov_num, group in self.grpd_fovs:
            all_counts_fov = []
            for index_value, fov_subdataset in group.iterrows():
                round_num = fov_subdataset.round_num
                channel = fov_subdataset.channel
                fov = fov_subdataset.fov_num
                experiment_name = fov_subdataset.experiment_name
                dask_delayed_name = 'filt_count_' +experiment_name + '_' + channel + \
                                '_round_' + str(round_num) + '_fov_' +str(fov) + '-' + tokenize()
                counts = delayed(fov_processing.single_fov_round_processing_eel, name=dask_delayed_name)(fov_subdataset,
                                            analysis_parameters,
                                            running_functions,
                                            dark_img,
                                            experiment_fpath,
                                            save_steps_output=save_intermediate_steps)
                all_counts_fov.append(counts)

            name = 'concat_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            all_counts_fov = delayed(pd.concat,name=name)(all_counts_fov,axis=0,ignore_index=True)
            
            name = 'register_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            registered_counts = delayed(fovs_registration.beads_based_registration,name=name)(all_counts_fov,
                                                analysis_parameters)

            name = 'decode_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            decoded = delayed(barcodes_analysis.extract_barcodes_NN_fast,name=name)(registered_counts, 
                                                                        analysis_parameters,codebook_df)                                                        
            
            name = 'stitch_to_mic_coords_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()  
            stitched_coords = delayed(stitching.stitch_using_microscope_fov_coords_new,name=name)(decoded[1],tile_corners_coords_pxl)
            
            name = 'save_df_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize() 
            saved_file = delayed(stitched_coords.to_parquet,name=name)(Path(experiment_fpath) / 'results'/ (experiment_name + \
                            '_decoded_fov_' + str(fov) + '.parquet'),index=False)
        
            all_processing.append(saved_file) 
        _ = dask.compute(*all_processing)
 
        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        selected_Hdistance = 3 / self.metadata['barcode_length']
        stitching_selected = 'microscope_stitched'
        io.simple_output_plotting(self.experiment_fpath, stitching_selected, selected_Hdistance, self.client)
        # ----------------------------------------------------------------  


    def processing_serial_fish_fov(self):
        """ 
        This method create a processing graph XXXXXXXXX

        """

        dark_img = preprocessing.load_dark_image(self.experiment_fpath)
        
        # did this conversion to avoid to pass self to dask
        analysis_parameters = self.analysis_parameters
        running_functions = self.running_functions
        tile_corners_coords_pxl = self.tile_corners_coords_pxl
        
        dark_img = delayed(dark_img)
        analysis_parameters = delayed(analysis_parameters)
        running_functions = delayed(running_functions)
        tile_corners_coords_pxl = delayed(tile_corners_coords_pxl)

        all_processing = []

        for fov_num, group in self.grpd_fovs:
            all_counts_fov = []
            all_nuclei_fov = []
            for index_value, fov_subdataset in group.iterrows():
                round_num = fov_subdataset.round_num
                channel = fov_subdataset.channel
                fov = fov_subdataset.fov_num
                experiment_name = fov_subdataset.experiment_name
                processing_type = fov_subdataset.processing_type

                if processing_type == 'nuclei':
                    dask_delayed_name = 'filt_' +experiment_name + '_' + channel + \
                                    '_round_' + str(round_num) + '_fov_' +str(fov) + '-' + tokenize()

                    out_nuclei = delayed(fov_processing.single_fov_round_processing_serial_nuclei,name=dask_delayed_name)(fov_subdataset,
                                            analysis_parameters,
                                            running_functions,
                                            dark_img,
                                            experiment_fpath,
                                            save_steps_output=save_intermediate_steps)
                    all_nuclei_fov.append(out_nuclei)

                else:
                    dask_delayed_name = 'filt_count_' +experiment_name + '_' + channel + \
                                    '_round_' + str(round_num) + '_fov_' +str(fov) + '-' + tokenize()
                    counts = delayed(fov_processing.single_fov_round_processing_eel,name=dask_delayed_name)(fov_subdataset,
                                                analysis_parameters,
                                                running_functions,
                                                dark_img,
                                                experiment_fpath,
                                                save_steps_output=save_intermediate_steps)
                    all_counts_fov.append(counts)
            
            name = 'concat_fish_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            all_counts_fov = delayed(pd.concat,name=name)(all_counts_fov,axis=0,ignore_index=True)
            
            name = 'create_nuclei_stack' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            filtered_nuclei_stack = delayed(utils.combine_filtered_images,name=name)(all_nuclei_fov,experiment_fpath,metadata)

            name = 'register_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            registered_counts = delayed(fovs_registration.nuclei_based_registration,name=name)(all_counts_fov,
                                                filtered_nuclei_stack,
                                                analysis_parameters)
                                                                                                
            name = 'stitch_to_mic_coords_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()  
            stitched_coords = delayed(stitching.stitch_using_microscope_fov_coords_new,name=name)(registered_counts,tile_corners_coords_pxl)
            
            name = 'save_df_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize() 
            saved_file = delayed(stitched_coords.to_parquet,name=name)(Path(experiment_fpath) / 'results'/ (experiment_name + \
                            '_decoded_fov_' + str(fov) + '.parquet'),index=False)
        
            all_processing.append(saved_file) 
        _ = dask.compute(*all_processing)

        # # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        stitching_selected = 'microscope_stitched'
        io.simple_output_plotting_serial(self.experiment_fpath, stitching_selected, self.client)
        # ----------------------------------------------------------------  
  


    # --------------------------------
    # RUNNING OPTIONS
    # --------------------------------
    
    def run_all(self):
        self.load_config_files()
        self.processing_cluster_init.start_processing_env()

        utils.create_folder_structure(self.experiment_fpath)
        self.logger.info(f'Folder structure completed')

        # Run parsing only if required
        if self.parsing_type != 'no_parsing':
            nikon_nd2_parsing_step()

        # Create or load processing dataset
        self.prepare_processing_dataset()

        # CREATE DARK IMAGES
        utils.create_dark_img(experiment_fpath,metadata)

        # CREATE RUNNING FUNCTIONS
        # used to select the sets of functions to run preprocessing and dots
        # calling according to the type of processing and sample
        running_functions = configuration_files.create_function_runner(experiment_fpath,metadata)

        # CREATE DARK IMAGES
        self.determine_tiles_organization()

        if self.metadata['experiment_type'] == 'eel-barcoded':
            self.processing_barcoded_eel_fov
        elif self.metadata['experiment_type'] == 'smfish-serial':
            self.processing_serial_fish_fov
        else:
            self.logger.error(f'the experiment type {self.metadata['experiment_type']} is unknown')
            sys.exit(f'the experiment type {self.metadata['experiment_type']} is unknown')