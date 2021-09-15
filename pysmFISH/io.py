"""
Module containing functions used to load and write data
"""

from typing import *
import zarr
import sys
import time
import numpy as np
import pandas as pd
from pathlib import Path
from dask import dataframe as dd
from dask.distributed import Client

from pysmFISH.utils import convert_from_uint16_to_float64
from pysmFISH.logger_utils import selected_logger


def create_empty_zarr_file(experiment_fpath:str,tag:str)-> str:
    """Function that create and empty zarr file 


    Args:
        experiment_fpath (str): location of the folder to be processed
        tag (str): string to add to the file name
    Returns:
        str: path of the created file
    """
    

    experiment_fpath = Path(experiment_fpath)
    experiment_name = experiment_fpath.stem
    zarr_fpath = experiment_fpath / (experiment_name + '_' + tag + '.zarr')
    
    store = zarr.DirectoryStore(zarr_fpath,'w')
    grp = zarr.group(store=store)
    return zarr_fpath


def consolidate_zarr_metadata(parsed_raw_data_fpath: str):
    """Function to consolidate all the zarr metadata in one unique
    json file for eady indexing and searching

    Args:
        parsed_raw_data_fpath (str): path to the zarr file for which
                    the metadata needs to be consolidated

    Returns:
        zarr groups instance with the consolidated metadata
    """

    logger = selected_logger()
    try:
        store = zarr.DirectoryStore(parsed_raw_data_fpath)
        consolidated_grp = zarr.consolidate_metadata(store)
    except:
        logger.error(f'cannot consolidate metadata of the parsed zarr file')
        sys.exit(f'cannot consolidate metadata of the parsed zarr file')
    else:
        return consolidated_grp


def open_consolidated_metadata(parsed_raw_data_fpath:str):
    """Load the consolidated json metadata file
    Args:
        parsed_raw_data_fpath (str): path to the zarr file with the
                consolidated metadata

    Returns:
        zarr groups instance with the consolidated metadata
    """
    logger = selected_logger()
    
    try:
        store = zarr.DirectoryStore(parsed_raw_data_fpath)
    except:
        logger.error(f'the metadata are not consolidated')
    else:
        consolidated_grp = zarr.open_consolidated(store)
        return consolidated_grp



def load_raw_images(zarr_grp_name:str,parsed_raw_data_fpath:str)->Tuple[np.ndarray,Dict]:
    """Function used to load a raw image and metadata from the 
    parsed raw file and the attrs for the filtering

    Args:
        zarr_grp_name (str): name to the group to process. The group contain the raw images and the 
            corresponding metadata.
            grp = experiment_name_channel_fov_X
            dataset = raw_data_fov_X
        parsed_raw_data_fpath (str): fpath to zarr store containing the parsed raw images

    Returns:
        Tuple[np.ndarray,Dict]: return the selected image and the corresponding metadata
    """
    logger = selected_logger()
    st = zarr.DirectoryStore(parsed_raw_data_fpath)
    root = zarr.group(store=st,overwrite=False)

    metadata = root[zarr_grp_name].attrs
    img = root[zarr_grp_name][metadata['fov_name']][...]

    return img, metadata



def load_general_zarr(fov_subdataset: pd.Series ,parsed_raw_data_fpath:str, tag:str)->Tuple[np.ndarray,Dict]:
    """Function used to load images stored in a zarr file (ex. preprocessed zarr)

    Args:
        fov_subdataset (pd.Series): Dataset metadata corresponding to as specific fov
        parsed_raw_data_fpath (str): path to the zarr file
        tag (str): string used to specify the zarr file

    Returns:
         Tuple[np.ndarray,Dict]: return the selected image and the corresponding metadata
    """
    logger = selected_logger()
    st = zarr.DirectoryStore(parsed_raw_data_fpath)
    root = zarr.group(store=st,overwrite=False)
    fov_name = tag + '_fov_' + str(fov_subdataset.fov_num)
    grp_name = fov_subdataset.experiment_name +'_' + fov_subdataset.channel + '_round_' + str(fov_subdataset.round_num) + '_fov_' + str(fov_subdataset.fov_num)
    img = root[grp_name][fov_name][...]
    img = convert_from_uint16_to_float64(img)
    metadata = root[grp_name].attrs
    return img, metadata



def simple_output_plotting(experiment_fpath: str, stitching_selected: str, 
                            selected_Hdistance: float, client, file_tag: str):
    """Utility function used to create a pandas dataframe with a simplified
    version of the eel analysis output that can be used for quick visualization

    Args:
        experiment_fpath (str): Path to the experiment to process
        stitching_selected (str): Define with stitched data will be selected
            for creating the simplified dataframe
        selected_Hdistance (float): Used to select the dots with hamming
            distance below this value
        client (Client): Dask client taking care of the processing 
        file_tag (str): tag to label the output file
    """

    experiment_fpath = Path(experiment_fpath)
    counts_dd = dd.read_parquet(experiment_fpath / 'results' / ('*' +file_tag+'*.parquet'),engine='pyarrow')

    date_tag = time.strftime("%y%m%d_%H_%M_%S")

    r_tag = 'r_px_' + stitching_selected
    c_tag = 'c_px_' + stitching_selected

    counts_dd_below  = counts_dd.loc[counts_dd.hamming_distance < selected_Hdistance, :]

    counts_df = counts_dd_below.loc[:,['fov_num',r_tag,c_tag, 'decoded_genes']].compute()

    counts_df=counts_df.dropna(subset=['decoded_genes'])
    fpath = experiment_fpath / 'results' / (date_tag + '_' + experiment_fpath.stem + '_data_summary_simple_plotting_'+file_tag+'.parquet')
    counts_df.to_parquet(fpath,index=False)


def simple_output_plotting_serial(experiment_fpath: str, stitching_selected: str, 
                                client,file_tag: str):
    """Utility function used to create a pandas dataframe with a simplified
    version of the serial analysis output that can be used for quick visualization

    Args:
        experiment_fpath (str): Path to the experiment to process
        stitching_selected (str): Define with stitched data will be selected
            for creating the simplified dataframe
        client (Client): Dask client taking care of the processing 
        file_tag (str): tag to label the output file

    """


    experiment_fpath = Path(experiment_fpath)
    counts_dd = dd.read_parquet(experiment_fpath / 'results' / ('*'+file_tag+'*.parquet'),engine='pyarrow')

    date_tag = time.strftime("%y%m%d_%H_%M_%S")

    r_tag = 'r_px_' + stitching_selected
    c_tag = 'c_px_' + stitching_selected


    counts_df = counts_dd.loc[:,['fov_num',r_tag,c_tag, 'target_name']].compute()

    counts_df=counts_df.dropna(subset=['target_name'])
    fpath = experiment_fpath / 'results' / (date_tag + '_' + experiment_fpath.stem + '_data_summary_simple_plotting_'+file_tag+'.parquet')
    counts_df.to_parquet(fpath,index=False)

