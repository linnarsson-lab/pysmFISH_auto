""" 
Collection of functions to run filtering 
"""
from typing import *
import pickle
import sys
import time

import dask
import numpy as np
import scipy.ndimage as nd
from skimage import filters, morphology, measure
from pathlib import Path
from dask.distributed import Client


# pysmFISH imports
from pysmFISH.io import load_raw_images
from pysmFISH.dots_calling import osmFISH_peak_based_detection
from pysmFISH.dots_calling import osmFISH_dots_thr_selection, osmFISH_dots_mapping
from pysmFISH.dots_calling import osmFISH_barcoded_peak_based_detection_masked_thr
from pysmFISH.utils import convert_from_uint16_to_float64
from pysmFISH.data_models import Output_models

from pysmFISH.logger_utils import selected_logger


def load_dark_image(experiment_fpath:str)->np.ndarray:
    """
    Function used to load the dark image  previously created and 
    saved in the extra folder present in the experiment

    Parameters:
    -----------
    experiment_fpath: str
        path to the experiment

    """

    logger = selected_logger()

    search = '*dark_img.npy'
    dark_img_folder = Path(experiment_fpath) / 'extra_processing_data'
    dark_img_fpath = list(dark_img_folder.glob(search))

    try:
        dark_img = np.load(dark_img_fpath[0])
    except:
        logger.error(f'cannot load the dark image file')
        sys.exit(f'cannot load the dark image file')
    else:
        logger.info(f'loaded {dark_img_fpath[0].stem} dark image')
        dark_img = convert_from_uint16_to_float64(dark_img)
        return dark_img



def standard_not_norm_preprocessing(
        zarr_grp_name,
        parsed_raw_data_fpath,
        processing_parameters,
        dark_img):

    """
   Standard function for filtering the fov image
   The image is not normalized
    
    Args:
    -----
        zarr_grp_name: str
            group representing the image to process
        parsed_raw_data_fpath: str
            path to the zarr file containing the parsed images
        processing_parameters: dict
            dictionary with the parameters used to process the images
    """

    logger = selected_logger()
    
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    experiment_fpath = parsed_raw_data_fpath.parent
    FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']
    FilteringSmallKernel=processing_parameters['PreprocessingFishFilteringSmallKernel']
    LaplacianKernel=processing_parameters['PreprocessingFishFilteringLaplacianKernel']
    

    try:
        raw_fish_images_meta = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)
    except:
        logger.error(f'cannot load {zarr_grp_name} raw fish image')
        sys.exit(f'cannot load {zarr_grp_name} raw fish image')
    else:
        logger.info(f'loaded {zarr_grp_name} raw fish image')
        img = raw_fish_images_meta[0]
        img_metadata = raw_fish_images_meta[1]
        img = convert_from_uint16_to_float64(img)

        img -= dark_img
        img[img<0] = 0

        background = filters.gaussian(img,FlatFieldKernel,preserve_range=False)
        img /= background
        img = nd.gaussian_laplace(img,LaplacianKernel)
        img = -img # the peaks are negative so invert the signal
        img[img<=0] = 0 # All negative values set to zero also = to avoid -0.0 issues
        img = np.abs(img) # to avoid -0.0 issues

        img = img.max(axis=0)

        return (img, img_metadata)



def filter_remove_large_objs(
        zarr_grp_name,
        parsed_raw_data_fpath,
        processing_parameters,
        dark_img):

    """
    Function to:
    - preprocess the fish images
    - count the dots and save the data
    
    Args:
    -----
        zarr_grp_name: str
            group representing the image to process
        parsed_raw_data_fpath: str
            path to the zarr file containing the parsed images
        processing_parameters: dict
            dictionary with the parameters used to process the images
    """

    logger = selected_logger()
    
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    experiment_fpath = parsed_raw_data_fpath.parent
    FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']
    FilteringSmallKernel=processing_parameters['PreprocessingFishFilteringSmallKernel']
    LaplacianKernel=processing_parameters['PreprocessingFishFilteringLaplacianKernel']

    LargeObjRemovalPercentile = processing_parameters['LargeObjRemovalPercentile']
    LargeObjRemovalMinObjSize = processing_parameters['LargeObjRemovalMinObjSize']
    LargeObjRemovalSelem = processing_parameters['LargeObjRemovalSelem']


    try:
        raw_fish_images_meta = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)
    except:
        logger.error(f'cannot load {zarr_grp_name} raw fish image')
        sys.exit(f'cannot load {zarr_grp_name} raw fish image')
    else:
        logger.info(f'loaded {zarr_grp_name} raw fish image')

        img = raw_fish_images_meta[0]
        img_metadata = raw_fish_images_meta[1]
        img = convert_from_uint16_to_float64(img)

        img -= dark_img
        img[img<0] = 0

        background = filters.gaussian(img,FlatFieldKernel,preserve_range=False)
        img /= background
        img = nd.gaussian_laplace(img,LaplacianKernel)
        img = -img # the peaks are negative so invert the signal
        img[img<=0] = 0 # All negative values set to zero also = to avoid -0.0 issues
        img = np.abs(img) # to avoid -0.0 issues
        
        img = img.max(axis=0)


        mask = np.zeros_like(img)
        idx=  img > np.percentile(img,LargeObjRemovalPercentile)
        mask[idx] = 1
    
        labels = nd.label(mask)

        properties = measure.regionprops(labels[0])    
        for ob in properties:
            if ob.area < LargeObjRemovalMinObjSize:
                mask[ob.coords[:,0],ob.coords[:,1]]=0

        mask = np.logical_not(mask)

        masked_img = img*mask

        
        return img, masked_img, img_metadata


    def large_beads_preprocessing(zarr_grp_name,
        parsed_raw_data_fpath,
        processing_parameters,
        dark_img):

        logger = selected_logger()
        
        parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
        experiment_fpath = parsed_raw_data_fpath.parent
        FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']

        raw_fish_images_meta = load_raw_images(zarr_grp_name,
                                        parsed_raw_data_fpath)

        img = raw_fish_images_meta[0]
        img_metadata = raw_fish_images_meta[1]
        img = convert_from_uint16_to_float64(img)
        img -= dark_img
        img[img<0] = 0
        img = np.abs(img) # to avoid -0.0 issues

        img = img.max(axis=0)

        img /= filters.gaussian(img,FlatFieldKernel,preserve_range=False)


        return img, img_metadata

