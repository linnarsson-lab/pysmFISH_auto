'''
Functions used to preprocess the raw images
'''

from typing import *
import logging
import sys
import numpy as np
from skimage import filters, measure, morphology
import scipy.ndimage as nd
from pathlib import Path

from pysmFISH.utils import convert_from_uint16_to_float64


def load_dark_image(experiment_fpath:str)->np.ndarray:
    """
    Function used to load the dark image  previously created and 
    saved in the extra folder present in the experiment

    Parameters:
    -----------
    experiment_fpath: str
        path to the experiment

    """

    logger = logging.getLogger(__name__)

    dark_img_fpath = Path(experiment_fpath) / 'extra_processing_data' / 'dark_img.npy'

    try:
        dark_img = np.load(dark_img_fpath)
    except:
        logger.error(f'cannot load the dark image file')
        sys.exit(f'cannot load the dark image file')
    
    dark_img = convert_from_uint16_to_float64(dark_img)
    return dark_img

def dark_image_removal(img_stack:np.ndarray,dark_image:np.ndarray)-> np.ndarray:
    """
    Function to remove the dark image with the background
    noise of the detector

    Parameters:
    ----------
    img_stack: np.ndarray float64
        raw image stack to process

    dark_image: np.ndarray float64
        image containing the background noise of the detector

    Returns:
    --------

    img_stack: np.ndarray float64
        image stack after removal of the detector noise

    """

    img_stack = img_stack - dark_image
    
    return img_stack



def flat_field_correction(img_stack:np.ndarray,filtering_element:np.ndarray)-> np.ndarray:
    """
    Function to correct the vignetting caused by missmatching between 
    light guide shape and camera chip

    Parameters:
    ----------
    img_stack: np.ndarray float64
        raw image stack to process

    filtering_element: np.ndarray float64
        vector with the sigma of the gaussian filter used to smooth
        the image

    Returns:
    --------

    img_stack: np.ndarray float64
        image stack after removal of the vignetting

    """

    # test if it is better to preserve_rangebool = True
    img_stack_filtered = filters.gaussian(img_stack,sigma=filtering_element,preserve_range=True)
    img_stack = img_stack / img_stack_filtered

    return img_stack


def flat_field_correction_filt_staining(img_stack:np.ndarray,filtering_element:np.ndarray)-> np.ndarray:
    """
    Function to correct the vignetting caused by missmatching between 
    light guide shape and camera chip

    Parameters:
    ----------
    img_stack: np.ndarray float64
        raw image stack to process

    filtering_element: np.ndarray float64
        vector with the sigma of the gaussian filter used to smooth
        the image

    Returns:
    --------

    img_stack: np.ndarray float64
        image stack after removal of the vignetting

    """

    # test if it is better to preserve_rangebool = True
    img_stack_filtered = filters.gaussian(img_stack,sigma=filtering_element,preserve_range=True)
    img_stack = img_stack / img_stack_filtered

    # Clean the image from the background
    img_stack = img_stack-img_stack_filtered
    # Remove the negative values        
    img_stack[img_stack<0] = 0

    return img_stack



def normalization_nn(img_stack):
    """
    Normalization of the cleaned images to reduce the 
    effect of edges and improve selection of the dots
    """
    img_mean_z = img_stack.mean(axis=(1,2))
    img_mean_z = img_mean_z[:,np.newaxis,np.newaxis]
    img_std_z = img_stack.std(axis=(1,2))
    img_std_z = img_std_z[:,np.newaxis,np.newaxis]
    img_nn= (img_stack - img_mean_z)/ img_std_z
    img_nn = img_nn.max(axis=0)
    return img_nn




def smFISH_filtering(img_stack:np.ndarray,small_kernel_size_sigma:list,kernel_size_sigma_laplacian:list)->np.ndarray:
    
    """
    functio used to remove the background from the smFISH and enhance the dots.

    Parameters:
    -----------
    img_stack: np.array float64
        raw image to filter
    small_kernel_size_sigma: list 
        list with the kernel size used to remove the dots to identify the background
    kernel_size_sigma_laplacian: list
        list with the kernel used to enhanche the dots signal

    Returns
    -----------
    img_stack: np.array float64 
        img stack after filtering
    """

    # Use a gaussian with kernel bigger than the dots to estimate the background
    # and remove it from the image
    img_stack = img_stack-filters.gaussian(img_stack,small_kernel_size_sigma,preserve_range=True)
    img_stack[img_stack<0] = 0

    # Enhance the dots
    img_stack = nd.gaussian_laplace(img_stack,kernel_size_sigma_laplacian)
    img_stack = -img_stack # the peaks are negative so invert the signal
    img_stack[img_stack<0] = 0 # All negative values set to zero 

    # Flatten the image
    # flattened_img = np.amax(img_stack,axis=0)
    # flattened_img = convert_to_uint16(flattened_img)

    return img_stack



def large_structures_filtering(img_stack:np.ndarray,large_kernel_size_sigma:list)->np.ndarray:
        
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

        # Clean the image from the background
        img_stack = img_stack-filters.gaussian(img_stack,sigma=large_kernel_size_sigma,preserve_range=True)
        # Remove the negative values        
        img_stack[img_stack<0] = 0
        # Flatten the image
        # flattened_img = np.amax(self.img_stack,axis=0)
        # flattened_img = convert_to_uint16(self.flattened_
        return img_stack

        

def removal_large_objects(img_stack:np.ndarray, percentile_level:int,min_obj_size:int,selem_size:int)->np.ndarray:

    """
    function used to remove large bright contaminants in the image.
    Can be used to remove lipofuscin from the image

    Parameters:
    -----------
    img_stack: np.ndarray float64
        image stack to process
    percentile_level: int
        value used to identify the bright objects to remove
    min_obj_size: int
        size in px above which the identified objects are removed
    selem_size: int
        size of the structuring element use to enlarge the mask region

    Returns:
    --------
    mask: np.ndarray float64
        mask used to filter out the objects before normalization

    """

    # img_max = np.amax(img_stack,axis=0)
    mask = np.zeros_like(img_stack)
    idx=  img_stack > np.percentile(img_stack,percentile_level)
    mask[idx] = 1
    # mask[~idx] = 0
    labels = nd.label(mask)

    properties = measure.regionprops(labels[0])    
    for ob in properties:
        if ob.area < min_obj_size:
            mask[ob.coords[:,0],ob.coords[:,1]]=0


    mask = morphology.binary_dilation(np.logical_not(mask), selem=morphology.disk(selem_size))
    mask = np.logical_not(mask)

    return mask


def removal_large_objects_max_proj(img:np.ndarray, percentile_level:int,min_obj_size:int)->np.ndarray:

    """
    function used to remove large bright contaminants in the image.
    Can be used to remove lipofuscin from the image

    Parameters:
    -----------
    img: 2D np.ndarray float64
        image stack to process
    percentile_level: int
        value used to identify the bright objects to remove
    min_obj_size: int
        size in px above which the identified objects are removed

    Returns:
    --------
    mask: np.ndarray float64
        mask used to filter out the objects before normalization

    """

    mask = np.zeros_like(img)
    idx=  img > np.percentile(img,percentile_level)
    mask[idx] = 1
   
    labels = nd.label(mask)

    properties = measure.regionprops(labels[0])    
    for ob in properties:
        if ob.area < min_obj_size:
            mask[ob.coords[:,0],ob.coords[:,1]]=0

    mask = np.logical_not(mask)

    return mask



def mask_image(img_stack:np.ndarray, mask:np.ndarray)->np.ndarray:
    """
    function used to mask an image stack and remove or select objects
    ex. for the removal of lipofuscin in the filtered image stack

    Parameters:
    -----------
    img_stack: np.ndarray float64
        image stack to process

    mask: np.ndarray bool
        mask used to remove or select objects from the image
        of interest
    
    
    Returns:
    --------
    img_stack: np.ndarray float64
        masked image stack
    """

    img_stack = img_stack*mask

    return img_stack


def combine_rounds(prepreocessed_list:list):
    """
    This function is used to collect all the flattened images
    corresponding to the different rounds of a fov.

    The output is an image stack in whick the number of the z planes
    correspond to the nuber of the round.

    Parameters:
    ----------
    prepreocessed_list: list
        list of tuples (flattened_img, fov_name, round_name)

    Returns:
    --------

    img_stack: np.ndarray float64
        image stack with the filtered images
    
    """
    img_z = len(prepreocessed_list)
    img_r,img_c = prepreocessed_list[0][0].shape

    img_stack = np.zeros([img_z,img_r,img_c])

    for img, fov_name, round_name in prepreocessed_list:
        round_number = int(round_name.split('_')[-1])-1
        img_stack[round_number,:,:] = img

    return img_stack