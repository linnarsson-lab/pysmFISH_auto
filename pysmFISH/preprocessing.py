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
    min_distance=processing_parameters['CountingFishMinObjDistance']
    min_obj_size=processing_parameters['CountingFishMinObjSize']
    max_obj_size=processing_parameters['CountingFishMaxObjSize']
    num_peaks_per_label=processing_parameters['CountingFishNumPeaksPerLabel']


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
        # img = np.abs(img) # to avoid -0.0 issues

        img = img.max(axis=0)

        return (img, img_metadata)