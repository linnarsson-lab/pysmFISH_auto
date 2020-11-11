'''
Functions used to preprocess the raw images
'''

from typing import *
import numpy as np
import scipy.ndimage as nd
from skimage import filters
from pathlib import Path

from pysmFISH.utils import convert_from_uint16_to_float64

from prefect import task
from prefect.engine import signals
from prefect.utilities.logging import get_logger


from pysmFISH.logger_utils import prefect_logging_setup


# @task(name='load-dark-image')
def load_dark_image(experiment_fpath:str)->np.ndarray:
    """
    Function used to load the dark image  previously created and 
    saved in the extra folder present in the experiment

    Parameters:
    -----------
    experiment_fpath: str
        path to the experiment

    """

    logger = get_logger()

    search = '*dark_img.npy'
    dark_img_folder = Path(experiment_fpath) / 'extra_processing_data'
    dark_img_fpath = list(dark_img_folder.glob(search))

    try:
        dark_img = np.load(dark_img_fpath[0])
    except:
        logger.error(f'cannot load the dark image file')
        err = signals.FAIL(f'cannot load the dark image file')
        raise err
    else:
        logger.info(f'loaded {dark_img_fpath[0].stem} dark image')
        dark_img = convert_from_uint16_to_float64(dark_img)
        return dark_img


# @task(name='preprocessing-raw-fish-images')
def preprocessing_dot_raw_image(img_meta:tuple,dark_img:np.ndarray,
                            FlatFieldKernel:np.ndarray,FilteringSmallKernel:np.ndarray, 
                            LaplacianKernel:np.ndarray )->np.ndarray:
    """
    Function that preprocess the images in order to enhance the dots

    Args:
        img_meta: tuple
            tuple containing (image np.ndarray and metadata dict)
        dark_img: np.ndarray
            image for correction of dark current
        FlatFieldKernel: np.ndarray
            Kernel to remove the large object and calculating the background image
            for correction of vignetting
        FilteringSmallKernel: np.ndarray
            Kernel to remove dots signal to calculate the background of the image
        LaplacianKernel: np.ndarray
            Kernel used to enhance the dots signal

    """
    logger = get_logger()
    img = img_meta[0]
    img_metadata = img_meta[1]
    img = convert_from_uint16_to_float64(img)
    img -= dark_img
    img = np.amax(img, axis=0)
    background = filters.gaussian(img,FlatFieldKernel,preserve_range=True)
    img -= filters.gaussian(img,FilteringSmallKernel,preserve_range=True)
    img[img<0] = 0
    img /= background
    img = nd.gaussian_laplace(img,LaplacianKernel)
    img = -img
    img[img<0] = 0
    img = (img - np.mean(img)) / np.std(img)
    img[img<0] = 0 
    logger.info(f'max {img.max()}')
    
    return (img, img_metadata)


@task(name='preprocessing-images')
def test_preprocessing_large_scale(img_meta:tuple,
                    experiment_fpath:str,
                    FlatFieldKernel:np.ndarray,FilteringSmallKernel:np.ndarray, 
                    LaplacianKernel:np.ndarray):

    dark_img = load_dark_image(experiment_fpath)
    img_meta = preprocessing_dot_raw_image(img_meta,dark_img,
                            FlatFieldKernel,FilteringSmallKernel, 
                            LaplacianKernel)
    return img_meta