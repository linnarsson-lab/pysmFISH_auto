"""
Set of functions used to load and write data to the shoji database
that will contain all the analysis outputs
"""

from typing import *
import os
import zarr
import sys
import time
# import shoji
import numpy as np
from pathlib import Path
from dask import dataframe as dd
import pandas as pd

from pysmFISH.logger_utils import selected_logger
from pysmFISH.utils import convert_to_uint16



def create_empty_zarr_file(experiment_fpath:str,tag:str)-> str:

    """
    Function that create and empty zarr file 

    Args:
        experiment_fpath: str
            location of the folder to be processed
        tag: str
            tag to add to the file name
    Return:
        empty_fpath: str
            path of the created file

    """
    
    experiment_fpath = Path(experiment_fpath)
    experiment_name = experiment_fpath.stem
    zarr_fpath = experiment_fpath / (experiment_name + '_' + tag + '.zarr')
    
    store = zarr.DirectoryStore(zarr_fpath,'w')
    grp = zarr.group(store=store)
    return zarr_fpath


def consolidate_zarr_metadata(parsed_raw_data_fpath):
    """
    Function to consolidate all the zarr metadata in one unique
    json file for eady indexing and searching

    Args:
        parsed_raw_data_fpath: str
            path to the file with all the parsed images
    
    Returns:
        consolidated_grp: zarr group
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
    logger = selected_logger()
    
    try:
        store = zarr.DirectoryStore(parsed_raw_data_fpath)
    except:
        logger.error(f'the metadata are not consolidated')
    else:
        consolidated_grp = zarr.open_consolidated(store)
        return consolidated_grp


def load_raw_images(zarr_grp_name:str,parsed_raw_data_fpath:str)->np.ndarray:
    """
    Function used to load a raw image and metadata from the 
    parsed raw file and the attrs for the filtering
        parsed_raw_data_fpath: str
            fpath to zarr store containing the parsed raw images
        zarr_grp_name: str
            fpath to the group to process. The group contain the raw images and the 
            corresponding metadata

            grp = experiment_name_channel_fov_X
                dataset = raw_data_fov_X

    """
    logger = selected_logger()
    st = zarr.DirectoryStore(parsed_raw_data_fpath)
    root = zarr.group(store=st,overwrite=False)

    metadata = root[zarr_grp_name].attrs
    img = root[zarr_grp_name][metadata['fov_name']][...]
    return (img, dict(metadata))



def load_zarr_fov(zarr_fpath:str, fov:int):
    """
    Function used to load in memory the raw image stacks,
    it will also convert the data into float64

    Parameters:
    -----------
    zarr_fpah: str 
        path to the file containing the raw images
    fov: int
        number of the fov to load
    
    Returns:
    --------
    img_stack: float64
        raw image
    """

    logger = selected_logger()

    try:
        raw_store = zarr.DirectoryStore(zarr_fpath)
        raw_root = zarr.group(store=raw_store,overwrite=False)
    except:
        logger.error(f'cannot load the zarr file {zarr_fpath}')
        sys.exit(f'cannot load the zarr file {zarr_fpath}')
    
    try:
        img_stack = raw_root[fov][...]
    except:
        logger.error(f'cannot load image fov {fov} from {zarr_fpath}')
        sys.exit(f'cannot load image fov {fov} from {zarr_fpath}')

    return img_stack



def simple_output_plotting(experiment_fpath, stitching_selected, select_genes, client):

    experiment_fpath = Path(experiment_fpath)
    counts_dd = dd.read_parquet(experiment_fpath / 'results' / '*_cleaned_df*.parquet')

    date_tag = time.strftime("%y%m%d_%H_%M_%S")

    r_tag = 'r_px_' + stitching_selected
    c_tag = 'c_px_' + stitching_selected

    counts_df = counts_dd.loc[(counts_dd.dot_id == counts_dd.barcode_reference_dot_id),
                                ['fov_num',r_tag,c_tag, select_genes]]

    counts_df=counts_df.dropna(subset=[select_genes]).compute()
    fpath = experiment_fpath / 'results' / (date_tag + '_' + experiment_fpath.stem + '_data_summary_simple_plotting.parquet')
    counts_df.to_parquet(fpath,index=False)
# def connect_to_shoji_smfish_experiment(experiment_name: str):
    
#     logger = prefect_logging_setup(f'connect to shoji')
#     try:
#             db = shoji.connect()
#     except:
#         logger.error(f'Cannot connect to shoji DB')
#         err = signals.FAIL(f'Cannot connect to shoji DB')
#         raise err
#     else:
#         try:
#             ws = db.FISH[experiment_name]
#         except:
#             logger.error(f'experiment workspace missing')
#             err = signals.FAIL(f'experiment workspace missing')
#             raise err
#         else:
#             try:
#                 dots_ws = ws.dots_data
#             except:
#                 logger.error(f'dots_data workspace missing')
#                 err = signals.FAIL(f'dots_data workspace missing')
#                 raise err

#             else:
#                 try:
#                     images_properties_ws = ws.images_properties
#                 except:
#                     logger.error(f'images_properties_ws workspace missing')
#                     err = signals.FAIL(f'images_properties_ws workspace missing')
#                     raise err
#                 else:
#                     try:
#                         experiment_properties_ws = ws.experiment_properties
#                     except:
#                         logger.error(f'experiment_properties_ws workspace missing')
#                         err = signals.FAIL(f'experiment_properties_ws workspace missing')
#                         raise err
#                     else:
#                         try:
#                             analysis_parameters_ws = ws.analysis_parameters
#                         except:
#                             logger.error(f'analysis_parameters_ws workspace missing')
#                             err = signals.FAIL(f'analysis_parameters_ws workspace missing')
#                             raise err
                        
#                         else:
#                             return dots_ws, images_properties_ws, experiment_properties_ws, analysis_parameters_ws


# class load_analysis_parameters(Task):

#     """
#     Task to load all the possible parameters that can be
#     selected for the preprocessing 
#     NB: If you add additional parameters you are required to
#         modify the code and add them to the dictionary
    
#     Args:
#     -----
#         experiment_name: str
#                 name of the experiment. Used to access the parameters in the shoji db
#     """
    
#     def run(self, experiment_name):
#         """
#         Task to load all the possible parameters that can be
#         selected for the preprocessing 
#         NB: If you add additional parameters you are required to
#             modify the code and add them to the dictionary
        
#         Args:
#         -----
#             experiment_name: str
#                     name of the experiment. Used to access the parameters in the shoji db
        
#         Returns:
#         --------
#             analysis_parameters: dict
#                 dictionary with all the parameters that will be used for the analysis
#         """
    
#         dots_ws, images_properties_ws, experiment_properties_ws, analysis_parameters_ws = connect_to_shoji_smfish_experiment(experiment_name)
#         # Collect all parameters
#         analysis_parameters = {}
#         analysis_parameters['fish'] = {}
#         analysis_parameters['small-beads'] = {}
#         analysis_parameters['large-beads'] = {}
#         analysis_parameters['staining'] = {}
#         analysis_parameters['fresh-nuclei'] = {}
#         analysis_parameters['BarcodesExtractionResolution'] = analysis_parameters_ws[:].BarcodesExtractionResolution            
#         analysis_parameters['RegistrationReferenceHybridization'] = analysis_parameters_ws[:].RegistrationReferenceHybridization

#         analysis_parameters['fish']['PreprocessingFishFlatFieldKernel'] = analysis_parameters_ws[:].PreprocessingFishFlatFieldKernel
#         analysis_parameters['fish']['PreprocessingFishFilteringSmallKernel'] = analysis_parameters_ws[:].PreprocessingFishFilteringSmallKernel
#         analysis_parameters['fish']['PreprocessingFishFilteringLaplacianKernel'] = analysis_parameters_ws[:].PreprocessingFishFilteringLaplacianKernel
#         analysis_parameters['fish']['CountingFishMinObjDistance'] = analysis_parameters_ws[:].CountingFishMinObjDistance
#         analysis_parameters['fish']['CountingFishMaxObjSize'] = analysis_parameters_ws[:].CountingFishMinObjSize
#         analysis_parameters['fish']['CountingFishMinObjSize'] = analysis_parameters_ws[:].CountingFishMinObjSize
#         analysis_parameters['fish']['CountingFishNumPeaksPerLabel'] = analysis_parameters_ws[:].CountingFishNumPeaksPerLabel
        
#         analysis_parameters['small-beads']['PreprocessingBeadsRegistrationFlatFieldKernel'] = analysis_parameters_ws[:].PreprocessingSmallBeadsRegistrationFlatFieldKernel
#         analysis_parameters['small-beads']['PreprocessingBeadsRegistrationFilteringSmallKernel'] = analysis_parameters_ws[:].PreprocessingSmallBeadsRegistrationFilteringSmallKernel
#         analysis_parameters['small-beads']['PreprocessingBeadsRegistrationFilteringLaplacianKernel'] = analysis_parameters_ws[:].PreprocessingSmallBeadsRegistrationFilteringLaplacianKernel 
#         analysis_parameters['small-beads']['CountingBeadsRegistrationMinObjDistance'] = analysis_parameters_ws[:].CountingSmallBeadsRegistrationMinObjDistance
#         analysis_parameters['small-beads']['CountingBeadsRegistrationMinObjSize'] = analysis_parameters_ws[:].CountingSmallBeadsRegistrationMinObjSize
#         analysis_parameters['small-beads']['CountingBeadsRegistrationMaxObjSize'] = analysis_parameters_ws[:].CountingSmallBeadsRegistrationMaxObjSize
#         analysis_parameters['small-beads']['CountingBeadsRegistrationNumPeaksPerLabel'] = analysis_parameters_ws[:].CountingSmallBeadsRegistrationNumPeaksPerLabel 
        
#         analysis_parameters['large-beads']['PreprocessingBeadsRegistrationFlatFieldKernel'] = analysis_parameters_ws[:].PreprocessingLargeBeadsRegistrationFlatFieldKernel
#         analysis_parameters['large-beads']['PreprocessingBeadsRegistrationFilteringSmallKernel'] = analysis_parameters_ws[:].PreprocessingLargeBeadsRegistrationFilteringSmallKernel
#         analysis_parameters['large-beads']['PreprocessingBeadsRegistrationFilteringLaplacianKernel'] = analysis_parameters_ws[:].PreprocessingLargeBeadsRegistrationFilteringLaplacianKernel
#         analysis_parameters['large-beads']['CountingBeadsRegistrationMinObjDistance'] = analysis_parameters_ws[:].CountingLargeBeadsRegistrationMinObjDistance
#         analysis_parameters['large-beads']['CountingBeadsRegistrationMinObjSize'] = analysis_parameters_ws[:].CountingLargeBeadsRegistrationMinObjSize
#         analysis_parameters['large-beads']['CountingBeadsRegistrationMaxObjSize'] = analysis_parameters_ws[:].CountingLargeBeadsRegistrationMaxObjSize
#         analysis_parameters['large-beads']['CountingBeadsRegistrationNumPeaksPerLabel'] = analysis_parameters_ws[:].CountingLargeBeadsRegistrationNumPeaksPerLabel
        
#         analysis_parameters['staining']['PreprocessingStainingFlatFieldKernel'] = analysis_parameters_ws[:].PreprocessingStainingFlatFieldKernel
        
#         analysis_parameters['fresh-nuclei']['PreprocessingFreshNucleiLargeKernelSize'] =  analysis_parameters_ws[:].PreprocessingFreshNucleiLargeKernelSize

#         return analysis_parameters





# @task(name = 'load-preprocessing-parameters')
# def load_analysis_parameters(experiment_name:str):
#     """
#     Function to load all the possible parameters that can be
#     selected for the preprocessing 
#     NB: If you add additional parameters you are required to
#         modify the code and add them to the dictionary
#     experiment_name: str
#             name of the experiment. Used to access the parameters in the shoji db
#     """
#     logger = prefect_logging_setup(f'load-analysis-parameters')
    
#     dots_ws, images_properties_ws, experiment_properties_ws, analysis_parameters_ws = connect_to_shoji_smfish_experiment(experiment_name)
#     # Collect all parameters
#     analysis_parameters = {}
#     analysis_parameters['fish'] = {}
#     analysis_parameters['small-beads'] = {}
#     analysis_parameters['large-beads'] = {}
#     analysis_parameters['staining'] = {}
#     analysis_parameters['fresh-nuclei'] = {}
#     analysis_parameters['BarcodesExtractionResolution'] = analysis_parameters_ws[:].BarcodesExtractionResolution            
#     analysis_parameters['RegistrationReferenceHybridization'] = analysis_parameters_ws[:].RegistrationReferenceHybridization

#     analysis_parameters['fish']['PreprocessingFishFlatFieldKernel'] = analysis_parameters_ws[:].PreprocessingFishFlatFieldKernel
#     analysis_parameters['fish']['PreprocessingFishFilteringSmallKernel'] = analysis_parameters_ws[:].PreprocessingFishFilteringSmallKernel
#     analysis_parameters['fish']['PreprocessingFishFilteringLaplacianKernel'] = analysis_parameters_ws[:].PreprocessingFishFilteringLaplacianKernel
#     analysis_parameters['fish']['CountingFishMinObjDistance'] = analysis_parameters_ws[:].CountingFishMinObjDistance
#     analysis_parameters['fish']['CountingFishMaxObjSize'] = analysis_parameters_ws[:].CountingFishMinObjSize
#     analysis_parameters['fish']['CountingFishMinObjSize'] = analysis_parameters_ws[:].CountingFishMinObjSize
#     analysis_parameters['fish']['CountingFishNumPeaksPerLabel'] = analysis_parameters_ws[:].CountingFishNumPeaksPerLabel
    
#     analysis_parameters['small-beads']['PreprocessingBeadsRegistrationFlatFieldKernel'] = analysis_parameters_ws[:].PreprocessingSmallBeadsRegistrationFlatFieldKernel
#     analysis_parameters['small-beads']['PreprocessingBeadsRegistrationFilteringSmallKernel'] = analysis_parameters_ws[:].PreprocessingSmallBeadsRegistrationFilteringSmallKernel
#     analysis_parameters['small-beads']['PreprocessingBeadsRegistrationFilteringLaplacianKernel'] = analysis_parameters_ws[:].PreprocessingSmallBeadsRegistrationFilteringLaplacianKernel 
#     analysis_parameters['small-beads']['CountingBeadsRegistrationMinObjDistance'] = analysis_parameters_ws[:].CountingSmallBeadsRegistrationMinObjDistance
#     analysis_parameters['small-beads']['CountingBeadsRegistratiohMinObjSize'] = analysis_parameters_ws[:].CountingSmallBeadsRegistrationMinObjSize
#     analysis_parameters['small-beads']['CountingBeadsRegistrationMaxObjSize'] = analysis_parameters_ws[:].CountingSmallBeadsRegistrationMaxObjSize
#     analysis_parameters['small-beads']['CountingBeadsRegistrationNumPeaksPerLabel'] = analysis_parameters_ws[:].CountingSmallBeadsRegistrationNumPeaksPerLabel 
    
#     analysis_parameters['large-beads']['PreprocessingBeadsRegistrationFlatFieldKernel'] = analysis_parameters_ws[:].PreprocessingLargeBeadsRegistrationFlatFieldKernel
#     analysis_parameters['large-beads']['PreprocessingBeadsRegistrationFilteringSmallKernel'] = analysis_parameters_ws[:].PreprocessingLargeBeadsRegistrationFilteringSmallKernel
#     analysis_parameters['large-beads']['PreprocessingBeadsRegistrationFilteringLaplacianKernel'] = analysis_parameters_ws[:].PreprocessingLargeBeadsRegistrationFilteringLaplacianKernel
#     analysis_parameters['large-beads']['CountingBeadsRegistrationMinObjDistance'] = analysis_parameters_ws[:].CountingLargeBeadsRegistrationMinObjDistance
#     analysis_parameters['large-beads']['CountingBeadsRegistrationMinObjSize'] = analysis_parameters_ws[:].CountingLargeBeadsRegistrationMinObjSize
#     analysis_parameters['large-beads']['CountingBeadsRegistrationMaxObjSize'] = analysis_parameters_ws[:].CountingLargeBeadsRegistrationMaxObjSize
#     analysis_parameters['large-beads']['CountingBeadsRegistrationNumPeaksPerLabel'] = analysis_parameters_ws[:].CountingLargeBeadsRegistrationNumPeaksPerLabel
    
#     analysis_parameters['staining']['PreprocessingStainingFlatFieldKernel'] = analysis_parameters_ws[:].PreprocessingStainingFlatFieldKernel
    
#     analysis_parameters['fresh-nuclei']['PreprocessingFreshNucleiLargeKernelSize'] =  analysis_parameters_ws[:].PreprocessingFreshNucleiLargeKernelSize

#     return analysis_parameters


# @task(name = 'save-img-and-metadata')
# def save_images_metadata(img_metadata:Tuple):
#     """
#     Function used to store the images metadata in the shoji database

#     Args:
#         img_metadata : Tuple 
#         contanining (image,metadata) 
#             Filtered image uint16 and image_metadata dictionary with the metadata collected duriring the 
#             parsing or acquisition of the images
#     """

#     logger = prefect.utilities.logging.get_logger('save-image-metadata')
#     metadata = img_metadata[1]
#     experiment_name = metadata['experiment_name']
#     img = convert_to_uint16(img_metadata[0])
#     fov_acquisition_coords = np.array((metadata['fov_acquisition_coords_x'],
#                                        metadata['fov_acquisition_coords_y'],
#                                        metadata['fov_acquisition_coords_z']),dtype=np.float64)
#     img_shape = np.array((metadata['img_height'],
#                           metadata['img_width']),np.uint16)
    
    
#     try:
#         dots_ws, images_properties_ws, experiment_properties_ws, analysis_parameters_ws = connect_to_shoji_smfish_experiment(experiment_name)
#     except:
#         logger.error(f'cannot connect to shoji db')
#         signals.FAIL(f'cannot connect to shoji db')

#     else:
#         images_properties_ws.fov.append({
#             'GroupName': np.array([metadata['grp_name']],dtype=np.object).reshape(1,1,1),
#             'FovName' : np.array([metadata['fov_name']],dtype=np.object),
#             'FovNumber' : np.array([metadata['fov_num']],dtype=np.uint16),
#             'AcquistionChannel' : np.array([metadata['channel']],dtype=np.object),
#             'TargetName' : np.array([metadata['target_name']],dtype=np.object).reshape(1,1,1),
#             'ImageShape' : np.array(img_shape,dtype=np.uint16).reshape(1,1,1,2),
#             'PixelMicrons' : np.array([metadata['pixel_microns']],dtype=np.float64).reshape(1,1,1),
#             'HybridizationNumber' : np.array([metadata['hybridization_num']],dtype=np.uint8),
#             'PreprocessedImage': img[None][None][None],
#             'FovCoords': np.array(fov_acquisition_coords,dtype=np.float64).reshape(1,1,1,3),
#             'RegistrationShift' : np.array([0,0],dtype=np.float64).reshape(1,1,1,2),
#             'RegistrationError' : np.array([0],dtype=np.float64).reshape(1,1,1),
#             'StitchingShift' : np.array([0,0],dtype=np.float64).reshape(1,1,1,2),
#             'StitchingError' : np.array([0],dtype=np.float64).reshape(1,1,1),
#             'FieldsOfView': np.array([metadata['fields_of_view']],dtype=np.uint16)
#         })


# @task(name = 'save-dots-identification')
# def save_dots_data(filtered_img_meta:Tuple):
#     """
#     Function used to store the dots relative data in the shoji database

#     Args:
#         filtered_img_meta: Tuple
#         countains (counts_dict, img_metadata)
#         counts_dict : contains data relative to the counted dots
#             DotsCoordsFOV
#             DotID
#             FovNumber
#             DotIntensity
#             SelectedThreshold
#             DotChannel
#         img_metadata
#     """

#     logger = prefect_logging_setup(f'save-dots-identification')
#     counts_dict, metadata = filtered_img_meta
#     experiment_name = metadata['experiment_name']
    
#     dots_ws, images_properties_ws, experiment_properties_ws, analysis_parameters_ws = connect_to_shoji_smfish_experiment(experiment_name)
    
#     DotCoordsFOV = counts_dict['DotsCoordsFOV'].astype(np.float64)
#     DotCoordsFOV = DotCoordsFOV.reshape(DotCoordsFOV.shape[0],1,1,2)
    
#     DotIntensity = counts_dict['DotIntensity'].astype(np.float64)
#     DotIntensity = DotIntensity[:,np.newaxis,np.newaxis]
    
#     SelectedThreshold = counts_dict['SelectedThreshold'].astype(np.float64)
#     SelectedThreshold = SelectedThreshold[:,np.newaxis,np.newaxis]
    
#     DotChannel = counts_dict['DotChannel'].astype(np.object)
    
#     ProcessingType = np.repeat(metadata['processing_type'],DotIntensity.shape[0])
#     ProcessingType = ProcessingType.astype(np.object)
#     ProcessingType = ProcessingType[:,np.newaxis,np.newaxis]
    
#     HybridizationNumber = np.repeat(metadata['hybridization_num'],DotIntensity.shape[0])
#     HybridizationNumber = HybridizationNumber.astype(np.uint8)
    
#     DotsCoordsRegisteredFOV = np.zeros([DotIntensity.shape[0],1,1,2],dtype=np.float64)
    
#     DotsCoordsStitched = np.zeros([DotIntensity.shape[0],1,1,2],dtype=np.float64)

#     BarcodeReferenceDotID = np.repeat(str(),DotIntensity.shape[0])
#     BarcodeReferenceDotID = BarcodeReferenceDotID.astype(np.object)

#     RawBarcode = np.zeros([DotIntensity.shape[0],16],dtype=np.bool)
    
#     GeneID = np.repeat(str(),DotIntensity.shape[0])
#     GeneID = GeneID.astype(np.object)
#     # GeneID = GeneID.reshape(DotIntensity.shape[0],1)
    

#     HammingDistanceRawBarcode = np.zeros([DotIntensity.shape[0]],dtype=np.float64)

#     dots_ws.dots.append({
#         'DotCoordsFOV': DotCoordsFOV,
#         'DotID' : counts_dict['DotID'].astype(np.object),          
#         'FovNumber' :  counts_dict['FovNumber'].astype(np.uint16),                
#         'HybridizationNumber' :  HybridizationNumber,
#         'DotIntensity' : DotIntensity,            
#         'SelectedThreshold' : SelectedThreshold,            
#         'DotChannel' : DotChannel,               
#         'ProcessingType' : ProcessingType,               
#         'DotsCoordsRegisteredFOV' : DotsCoordsRegisteredFOV,
#         'DotsCoordsStitched' : DotsCoordsStitched,
#         'BarcodeReferenceDotID' : BarcodeReferenceDotID,
#         'RawBarcode' : RawBarcode,
#         'GeneID': GeneID,
#         'HammingDistanceRawBarcode' : HammingDistanceRawBarcode

#     })  

#     del counts_dict
#     del metadata



