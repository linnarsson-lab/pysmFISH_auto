""" 
Collection of functions to run filtering 
"""
from typing import *
import pickle
import sys
import time
import zarr

import dask
import numpy as np
import scipy.ndimage as nd
from skimage import filters, morphology, measure
from skimage import img_as_float64, img_as_uint
from pathlib import Path
from dask.distributed import Client


import cupy as cp
from cupyx.scipy import ndimage as ndx


# pysmFISH imports
from pysmFISH.io import load_raw_images
from pysmFISH.io import load_zarr_fov
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
        # dark_img = convert_from_uint16_to_float64(dark_img)
        dark_img = img_as_float64(dark_img)
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

        return img, img_metadata



def standard_norm_preprocessing(
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

        img_mean_z = img.mean(axis=(1,2))
        img_mean_z = img_mean_z[:,np.newaxis,np.newaxis]
        img_std_z = img.std(axis=(1,2))
        img_std_z = img_std_z[:,np.newaxis,np.newaxis]
        img_nn= (img - img_mean_z)/ img_std_z
        img_nn = img_nn.max(axis=0)
        img_nn[img_nn<=0] = 0 # All negative values set to zero

        return img_nn, img_metadata



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


def filter_remove_large_objs_gpu(
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
        img = img_as_float64(img)

        img -= dark_img
        img[img<0] = 0

        img = cp.array(img)

        background = ndx.gaussian_filter(img,FlatFieldKernel)
        img /= background
        img = ndx.gaussian_laplace(img,LaplacianKernel)
        img = -img # the peaks are negative so invert the signal
        img[img<=0] = 0 # All negative values set to zero also = to avoid -0.0 issues
        img = cp.abs(img) # to avoid -0.0 issues
        
        img = img.max(axis=0)

        mask = cp.zeros_like(img)
        idx=  img > cp.percentile(img,LargeObjRemovalPercentile)
        mask[idx] = 1
    
        labels = ndx.label(mask)

        labels = labels.get()
        mask = mask.get()
        img = img.get()

        properties = measure.regionprops(labels[0])    
        for ob in properties:
            if ob.area < LargeObjRemovalMinObjSize:
                mask[ob.coords[:,0],ob.coords[:,1]]=0

        mask = np.logical_not(mask)

        masked_img = img*mask

        
        return img, masked_img, img_metadata


def large_beads_preprocessing_gpu(zarr_grp_name,
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
    img = img_as_float64(img)
    img -= dark_img
    img[img<0] = 0
    img = np.abs(img) # to avoid -0.0 issues

    img = img.max(axis=0)

    img = cp.array(img)

    img /= ndx.gaussian_filter(img,FlatFieldKernel,preserve_range=False)

    img = cp.asnumpy(img)
    return img, img_metadata



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


def both_beads_preprocessing(zarr_grp_name,
        parsed_raw_data_fpath,
        processing_parameters,
        dark_img):

    """
    Function used to process only large beads in both-beads condition
    Used for testing experiment
    """

    logger = selected_logger()
    
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    experiment_fpath = parsed_raw_data_fpath.parent
    FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']

    raw_fish_images_meta = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)

    logger.info(f'loaded {zarr_grp_name} raw fish image')
    img = raw_fish_images_meta[0]
    img_metadata = raw_fish_images_meta[1]
    img = convert_from_uint16_to_float64(img)
    img -= dark_img
    img[img<0] = 0
    img = np.abs(img) # to avoid -0.0 issues

    img = img.max(axis=0)

    img /= filters.gaussian(img,FlatFieldKernel,preserve_range=False)


    return img, img_metadata


def nuclei_registration_filtering(zarr_grp_name,
        parsed_raw_data_fpath,
        processing_parameters,
        dark_img):
        
        """
        This function remove the background from large structures like nuclei
        For the sigma I seleced a value quite bigger than
        the nuclei size in order to remove them from the 
        image. I used what showed on the gaussian filter code page and on this
        link on stackoverflow: 
        http://stackoverflow.com/questions/25216382/gaussian-filter-in-scipy

        Parameters:
        -----------

        img_stack: np.array float64
            raw image to filter
        
        large_kernel_size_sigma: list
            list with the kernel size used to remove large objects in the image
            to identify the background

        Returns
        -----------

        img_stack: np.array float64 
            img stack after filtering 

        """

        parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
        experiment_fpath = parsed_raw_data_fpath.parent
        FlatFieldKernel=processing_parameters['PreprocessingNucleiFlatFieldKernel']
    
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
            # img -= filters.gaussian(img,FilteringSmallKernel,preserve_range=False)
            img[img<0] = 0

            background = filters.gaussian(img,FlatFieldKernel,preserve_range=False)
            img -= background
            img[img<=0] = 0 # All negative values set to zero also = to avoid -0.0 issues
            img = np.abs(img) # to avoid -0.0 issues
            img /= background
            
            img = img.max(axis=0)
        
            return img, img_metadata


def fresh_nuclei_filtering(
        parsed_raw_data_fpath,
        filtered_raw_data_fpath,
        fov,
        processing_parameters):
        
        """
        This function remove the background from large structures like nuclei
        For the sigma I seleced a value quite bigger than
        the nuclei size in order to remove them from the 
        image. I used what showed on the gaussian filter code page and on this
        link on stackoverflow: 
        http://stackoverflow.com/questions/25216382/gaussian-filter-in-scipy

        Arguments
        -----------

        img_stack: np.array float64
            3D numpy array with the image
        
        Returns
        -----------

        filtered_image: np.array float64 
            2D flattened image 

        """
        PreprocessingFreshNucleiLargeKernelSize = processing_parameters['fresh-nuclei']['PreprocessingFreshNucleiLargeKernelSize']

        filtered_store = zarr.DirectoryStore(filtered_raw_data_fpath)
        filtered_root = zarr.group(store=filtered_store,overwrite=False)

        img_stack = load_zarr_fov(parsed_raw_data_fpath,fov)
        img_stack = img_as_float64(img_stack)

        # Clean the image from the background
        img_stack = img_stack-filters.gaussian(img_stack,sigma=PreprocessingFreshNucleiLargeKernelSize)
        # Remove the negative values        
        img_stack[img_stack<0] = 0
        # Flatten the image
        flattened_img = np.amax(img_stack,axis=0)

        flattened_img = img_as_uint(flattened_img)
        filtered_root.create_dataset(fov, data=flattened_img, shape=flattened_img.shape, chunks=(None,None),overwrite=True)
        
