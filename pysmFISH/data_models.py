"""
Module containing classes and function used to organize data and metadata
"""
from typing import *
import time

import pandas as pd
import sys
import pickle
from pathlib import Path
import pandas as pd

from pysmFISH.io import open_consolidated_metadata
from pysmFISH.logger_utils import selected_logger


class Dataset():
    """Dataset class used to collect all the info related to the
    acquired images. It contains all the metadata needed during the
    processing.
    """
    
    def __init__(self):
        
        self.logger = selected_logger()
        
    def load_dataset(self, dataset_fpath: str):
        """Load a pre-existing dataset

        Args:
            dataset_fpath (str): Path to the existing dataset
        """
        self.dataset = pd.read_parquet(dataset_fpath)
        
    
    # TODO Add support for other file types
    def create_full_dataset_from_files(self, experiment_fpath: str, 
                                       experiment_info: dict,
                                       parsed_raw_data_fpath: str ,ftype: str = 'pkl'):
        
        """ Utility function that can be used to create the dataset from a 
        storage type different from the zarr structure used in pysmFISH. 
        processing. It requires a file of type ftype containing all the metadata

        

        Args:
            experiment_fpath (str): Path to the experiment to process
            experiment_info (dict): Dictionary with the configuration data
            parsed_raw_data_fpath (str): Path to the folder/zarr file with the parsed data
            ftype (str, optional): Path to the files containing the metadata. Defaults to 'pkl'.
        """    
        
        self.experiment_fpath = Path(experiment_fpath)
        self.experiment_info = experiment_info
        self.parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    
        date_tag = time.strftime("%y%m%d_%H_%M_%S")
        experiment_name = self.experiment_fpath.stem
        self.dataset_fpath = self.experiment_fpath / (date_tag + '_' + experiment_name + '_dataset.parquet')

        self.dataset = pd.DataFrame()
        all_pickle_list = list(self.parsed_raw_data_fpath.glob('*.' + ftype))
        if len(all_pickle_list):
            if ftype == 'pkl':
                for fdata in all_pickle_list:
                    single = pickle.load(open(fdata,'rb'))
                    fdata_loc = fdata.parent / fdata.stem
                    single['raw_data_location'] = fdata_loc.as_posix()
                    single_df = pd.DataFrame(single,index=[0])
                    self.dataset = pd.concat([self.dataset,single_df],axis=0,ignore_index=True)
                self.dataset.to_parquet(self.dataset_fpath, index=False)
            
            # TODO Add support for other file types
        else:
            self.logger.error(f'there are no files with the metadata dictionary')
            sys.exit(f'there are no files with the metadata dictionary')
    
    def create_full_dataset_from_zmetadata(self,
                                       parsed_raw_data_fpath: str):
        """Function used to create the full dataset starting from the consolidated
        metadata of the parsed zarr file

        Args:
            parsed_raw_data_fpath (str): Path to the zarr file with the parsed data
        """
        
        try:
            consolidated_metadata = open_consolidated_metadata(parsed_raw_data_fpath)
        except:
            self.logger.error(f'consolidated zarr metadata missing or broken')
            sys.exit(f'consolidated zarr metadata missing or broken')
        else:
            self.parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
            
            date_tag = time.strftime("%y%m%d_%H_%M_%S")
            experiment_name = self.parsed_raw_data_fpath.stem
            self.dataset_fpath = self.parsed_raw_data_fpath.parent / (date_tag + '_' + experiment_name + '_dataset.parquet')
            
            self.dataset = pd.DataFrame()
            for name, grp in consolidated_metadata.items():
                attrs_dict = dict(grp.attrs)
                fdata_loc = self.parsed_raw_data_fpath / attrs_dict['grp_name']
                attrs_dict['raw_data_location'] = fdata_loc.as_posix()
                attrs_df = pd.DataFrame(attrs_dict,index=[0])
                self.dataset = pd.concat([self.dataset,attrs_df],axis=0,ignore_index=True)
            self.dataset.to_parquet(self.dataset_fpath, index=False)
        
            
    def collect_metadata(self,df:pd.DataFrame)-> dict:
        """Function used to create a simplified metadata dictionary
        used in the processing. This create a smaller size variable
        compared to the dataset making it easier to send to all
        the workers during the processing

        Args:
            df (pd.DataFrame): Dataset DataFrame

        Returns:
            dict: Dictionary with the compiled metadata
        """
        self.metadata = {}
        self.metadata['list_all_fovs'] = df.fov_num.unique()
        self.metadata['list_all_channels'] = df.channel.unique()
        self.metadata['total_rounds'] = df.round_num.max()
        self.metadata['stitching_channel'] = df.iloc[0]['stitching_channel']
        self.metadata['img_width'] = df.iloc[0]['img_width']
        self.metadata['img_height'] = df.iloc[0]['img_height']
        self.metadata['img_zstack'] = df.iloc[0]['zstack']
        self.metadata['pixel_microns'] = df.iloc[0]['pixel_microns']
        self.metadata['experiment_name'] = df.iloc[0]['experiment_name']
        self.metadata['overlapping_percentage'] = df.iloc[0]['overlapping_percentage']
        self.metadata['machine'] = df.iloc[0]['machine']
        self.metadata['barcode_length'] = df.iloc[0]['barcode_length']
        self.metadata['processing_type'] = df.iloc[0]['processing_type']
        self.metadata['experiment_type'] = df.iloc[0]['experiment_type']
        self.metadata['pipeline'] =  df.iloc[0]['pipeline']
        self.metadata['stitching_type'] = df.iloc[0]['stitching_type']
        self.metadata['list_all_codebooks'] = df.codebook.unique()
        return self.metadata
        

    def grp_by_channel(self,df: pd.DataFrame)->pd.DataFrame:
        """Groupby the dataset by channel

        Args:
            df (pd.DataFrame): Dataset to groupby

        Returns:
            pd.DataFrame: Dataset grouped by channel
        """
        subset_df = self.dataset.groupby('channel')
        return subset_df
    
    def select_fovs_subset(self, df: pd.DataFrame, min_fov: int, max_fov: int)->pd.DataFrame:
        """Select a range of fovs between min and max values (extremites included)

        Args:
            df (pd.DataFrame): Dataset
            min_fov (int): lowest fov to select
            max_fov (int): highest fov to select

        Returns:
            pd.DataFrame: Sub-Dataset containing only the selected fovs 
        """
        subset_df = df.loc[(df.fov_num >= min_fov) & (df.fov_num <= max_fov),:]
        return subset_df
    
    
    def select_channel_subset(self,df: pd.DataFrame, channel: str)->pd.DataFrame:
        """Select all the fovs of a specific channel

        Args:
            df (pd.DataFrame): Dataset
            channel (str): Imaging channel of interest

        Returns:
            pd.DataFrame: Sub-Dataset containing the fovs of the
                    selected channel
        """
        subset_df = df.loc[(df.channel == channel),:]
        return subset_df
    
    
    def select_round_subset(self,df: pd.DataFrame, round_num: int)->pd.DataFrame:
        """Select all the fovs of a specific imaging round

        Args:
            df (pd.DataFrame): Dataset
            round_num (int): Imaging round of interest

        Returns:
            pd.DataFrame: Sub-Dataset containing the fovs of the
                    selected round
        """
        subset_df = df.loc[(df.round_num == round_num),:]
        return subset_df
    
    
    def select_specific_fov(self,df: pd.DataFrame, channel: str,round_num: int, fov_num: int)->pd.Series:
        """Select a specific fov according to fov number, channel and
        round number

        Args:
            df (pd.DataFrame): Dataset
            channel (str): Imaging channel of interest
            round_num (int): Imaging round of interest
            fov_num (int): Number of the fov of interest

        Returns:
            pd.Series: Sub-Dataset containing the fovs of the
                    selected round
        """
        subset_df = df.loc[(df.channel == channel) & 
                           (df.round_num == round_num) & 
                           (df.fov_num == fov_num),:]
        return subset_df.squeeze()

    def select_all_imgs_fov(self,df: pd.DataFrame, fovs_list: List)->pd.DataFrame:
        """Select a list of fovs

        Args:
            df (pd.DataFrame): Dataset
            fovs_list (List): List of fov number to select

        Returns:
            pd.DataFrame: Sub-Dataset containing the selected fovs
        """
        
        subset_df = df.loc[(df.fov_num.isin(fovs_list)),:]
        return subset_df
    
    # TODO additional utliti functions
    def load_to_numpy(self):
        pass
    
    def load_to_xarray(self):
        pass
    
    def merge_datasets(self):
        pass

    def save_dataset(self, df, fpath):
        df.to_parquet(fpath)


class Output_models():
    """
    Utility class containing the definition of the output dataframes
    """

    def __init__(self):

        self.output_variables_dot_calling = ['r_px_original','c_px_original',
                                        'dot_id','dot_intensity','selected_thr']

        self.output_variables_specific_registration = ['r_px_registered', 'c_px_registered',
                            'r_shift_px','c_shift_px','min_number_matching_dots_registration',
                            'reference_hyb', 'experiment_type','experiment_name','pxl_um',
                            'stitching_type', 'img_width_px','img_height_px','fov_acquisition_coords_x', 
                            'fov_acquisition_coords_y', 'fov_acquisition_coords_z']
        
        self.output_variables_specific_barcodes = ['barcodes_extraction_resolution', 'barcode_reference_dot_id',
                                'raw_barcodes','all_Hdistance_genes','zeroHdistance_genes','below2Hdistance_genes',
                                'below3Hdistance_genes','number_positive_bits','hamming_distance']
    
    
        self.dots_counts_dict = dict.fromkeys(self.output_variables_dot_calling)
        self.dots_counts_df = pd.DataFrame(columns=self.output_variables_dot_calling)

        # self.output_variables_registration = self.output_variables_dot_calling + self.output_variables_specific_registration
        # self.output_registration_df = pd.DataFrame(columns=self.output_variables_registration)

        # self.output_variables_barcodes = self.output_variables_registration + self.output_variables_specific_barcodes
        # self.barcode_analysis_df = pd.DataFrame(columns=self.output_variables_barcodes)

        self.barcode_analysis_df = pd.DataFrame(columns=self.output_variables_specific_barcodes)
