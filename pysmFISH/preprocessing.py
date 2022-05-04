""" 
Collection of functions to run filtering 

All the filtering functions will return a tuple of different length with 
in position 0 the image that will be used for downstream dots calling and 
last position the filtered image that will be saved.
Ex. 
standard_not_norm_preprocessing -> (img,): in this case img will be fed to the
counting algorithm and also saved

filter_remove_large_objs -> (masked_img, img): in this case masked_img will be
fed to the counting algorithm and img saved

Hypotetical_filt -> (processed, X, X, X , img): in this case processed will be
fed to the counting algorithm and img saved


"""
from typing import *
import sys

import numpy as np
import scipy.ndimage as nd
from skimage import filters, morphology, measure
from skimage import img_as_float64
from pathlib import Path


# pysmFISH imports
from pysmFISH.io import load_raw_images
from pysmFISH.utils import convert_from_uint16_to_float64
from pysmFISH.logger_utils import selected_logger
try:
    from csbdeep.utils import normalize
except:
    pass

def load_dark_image(experiment_fpath:str)->np.ndarray:
    """Function used to load the dark image  previously created and 
    saved in the extra folder present in the experiment

    Args:
        experiment_fpath (str): Path to the experiment to process

    Returns:
        np.ndarray: Dark image
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
        zarr_grp_name: str,
        parsed_raw_data_fpath: str,
        processing_parameters: dict,
        dark_img: np.ndarray)-> Tuple[Tuple[np.ndarray,],dict]:
    
    """Standard function for filtering the fov image
    The image is not normalized

    Args:
        zarr_grp_name (str): group name of the image to process
        parsed_raw_data_fpath (str): path to the zarr file containing the parsed images
        processing_parameters (dict): dictionary with the parameters used to process the images
        dark_img (np.ndarray): Dark image used to remove camera dark noise
    Returns:
        Tuple[Tuple[np.ndarray,],dict]: ((filtered_image,),metadata)
    """
 

    logger = selected_logger()
    
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']
    LaplacianKernel=processing_parameters['PreprocessingFishFilteringLaplacianKernel']
    
    try:
        img, metadata = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)
    except:
        logger.error(f'cannot load {zarr_grp_name} raw fish image')
        sys.exit(f'cannot load {zarr_grp_name} raw fish image')
    else:
        logger.info(f'loaded {zarr_grp_name} raw fish image')
        img = convert_from_uint16_to_float64(img)

        img -= dark_img
        img[img<0] = 0
        # (1,100,100)
        background = filters.gaussian(img,FlatFieldKernel,preserve_range=False)
        img /= background
        img = nd.gaussian_laplace(img,LaplacianKernel)
        img = -img # the peaks are negative so invert the signal
        img[img<=0] = 0 # All negative values set to zero also = to avoid -0.0 issues
        img = np.abs(img) # to avoid -0.0 issues

        img = img.max(axis=0)

        return ((img,),metadata)



def without_flat_field_correction(
        zarr_grp_name: str,
        parsed_raw_data_fpath: str,
        processing_parameters: dict,
        dark_img: np.ndarray)-> Tuple[Tuple[np.ndarray,],dict]:
    
    """Standard function for filtering the fov image where the
    flat field correction step is removed.
    The image is not normalized


    Args:
        zarr_grp_name (str): group name of the image to process
        parsed_raw_data_fpath (str): path to the zarr file containing the parsed images
        processing_parameters (dict): dictionary with the parameters used to process the images
        dark_img (np.ndarray): Dark image used to remove camera dark noise

    Returns:
        Tuple[Tuple[np.ndarray,],dict]: ((filtered_image,),metadata)
    """
  

    logger = selected_logger()
    
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']
    LaplacianKernel=processing_parameters['PreprocessingFishFilteringLaplacianKernel']
    

    try:
        img, metadata = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)
    except:
        logger.error(f'cannot load {zarr_grp_name} raw fish image')
        sys.exit(f'cannot load {zarr_grp_name} raw fish image')
    else:
        logger.info(f'loaded {zarr_grp_name} raw fish image')
        img = convert_from_uint16_to_float64(img)

        img -= dark_img
        img[img<0] = 0

        background = filters.gaussian(img,LaplacianKernel,preserve_range=False)
        img -= background
        img[img<=0] = 0
        img = nd.gaussian_laplace(img,LaplacianKernel)
        img = -img # the peaks are negative so invert the signal
        img[img<=0] = 0 # All negative values set to zero also = to avoid -0.0 issues
        img = np.abs(img) # to avoid -0.0 issues

        img = img.max(axis=0)

        return ((img,),metadata)


def standard_norm_preprocessing(
        zarr_grp_name: str,
        parsed_raw_data_fpath: str,
        processing_parameters: dict,
        dark_img: np.ndarray)-> Tuple[Tuple[np.ndarray,],dict]:
    
    """Function to preprocess and normalize the intensity of the
        fov. Normalization: (img_stack - img_stack_z_mean)/ img_stack_z_std

    Args:
        zarr_grp_name (str): group name of the image to process
        parsed_raw_data_fpath (str): path to the zarr file containing the parsed images
        processing_parameters (dict): dictionary with the parameters used to process the images
        dark_img (np.ndarray): Dark image used to remove camera dark noise

    Returns:
        Tuple[Tuple[np.ndarray,],dict]: ((filtered_normalized_image,),metadata)
    """
   

    logger = selected_logger()
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']
    LaplacianKernel=processing_parameters['PreprocessingFishFilteringLaplacianKernel']


    try:
        img, metadata = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)
    except:
        logger.error(f'cannot load {zarr_grp_name} raw fish image')
        sys.exit(f'cannot load {zarr_grp_name} raw fish image')
    else:
        logger.info(f'loaded {zarr_grp_name} raw fish image')
            
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

        return ((img_nn,),metadata)



def filter_remove_large_objs(
        zarr_grp_name: str,
        parsed_raw_data_fpath: str,
        processing_parameters: dict,
        dark_img: np.ndarray)-> Tuple[Tuple[np.ndarray,np.ndarray],dict]:

    """Function used to mask large objects (ex. lipofuscin) present in the image

    Args:
        zarr_grp_name (str): group name of the image to process
        parsed_raw_data_fpath (str): path to the zarr file containing the parsed images
        processing_parameters (dict): dictionary with the parameters used to process the images
        dark_img (np.ndarray): Dark image used to remove camera dark noise

    Returns:
        Tuple[Tuple[np.ndarray,],dict]: ((masked_image, filtered_image),metadata)
    """
  

    logger = selected_logger()
    
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']
    LaplacianKernel=processing_parameters['PreprocessingFishFilteringLaplacianKernel']
    LargeObjRemovalPercentile = processing_parameters['LargeObjRemovalPercentile']
    LargeObjRemovalMinObjSize = processing_parameters['LargeObjRemovalMinObjSize']
    LargeObjRemovalSelem = processing_parameters['LargeObjRemovalSelem']

    try:
        img, metadata = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)
    except:
        logger.error(f'cannot load {zarr_grp_name} raw fish image')
        sys.exit(f'cannot load {zarr_grp_name} raw fish image')
    else:
        logger.info(f'loaded {zarr_grp_name} raw fish image')

        img = img_as_float64(img)
        img -= dark_img
        img[img<0] = 0

        background = filters.gaussian(img,(1,5,5),preserve_range=False)
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

        mask = morphology.binary_dilation(mask, selem=morphology.disk(LargeObjRemovalSelem))
        mask = np.logical_not(mask)

        masked_img = img*mask

        
        return ((masked_img,img),metadata)


def filter_remove_large_objs_no_flat(
        zarr_grp_name: str,
        parsed_raw_data_fpath: str,
        processing_parameters: dict,
        dark_img: np.ndarray)-> Tuple[Tuple[np.ndarray,np.ndarray],dict]:

    """Function used to mask large objects (ex. lipofuscin) present in the image
        without flat field correction

    Args:
        zarr_grp_name (str): group name of the image to process
        parsed_raw_data_fpath (str): path to the zarr file containing the parsed images
        processing_parameters (dict): dictionary with the parameters used to process the images
        dark_img (np.ndarray): Dark image used to remove camera dark noise

    Returns:
        Tuple[Tuple[np.ndarray,],dict]: ((masked_image, filtered_image),metadata)
    """
  

    logger = selected_logger()
    
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']
    LaplacianKernel=processing_parameters['PreprocessingFishFilteringLaplacianKernel']
    LargeObjRemovalPercentile = processing_parameters['LargeObjRemovalPercentile']
    LargeObjRemovalMinObjSize = processing_parameters['LargeObjRemovalMinObjSize']
    LargeObjRemovalSelem = processing_parameters['LargeObjRemovalSelem']

    try:
        img, metadata = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)
    except:
        logger.error(f'cannot load {zarr_grp_name} raw fish image')
        sys.exit(f'cannot load {zarr_grp_name} raw fish image')
    else:
        logger.info(f'loaded {zarr_grp_name} raw fish image')

        img = img_as_float64(img)
        img -= dark_img
        img[img<0] = 0

        background = filters.gaussian(img,(1,5,5),preserve_range=False)
        img /= background
        img = nd.gaussian_laplace(img,LaplacianKernel)
        img = -img # the peaks are negative so invert the signal
        img[img<=0] = 0 # All negative values set to zero also = to avoid -0.0 issues
        img = np.abs(img) # to avoid -0.0 issues
        img = img.max(axis=0)
        #img = normalize(img,clip=True, dtype=np.float64)

        mask = np.zeros_like(img)
        idx=  img > np.percentile(img,LargeObjRemovalPercentile)
        mask[idx] = 1
    
        labels = nd.label(mask)

        properties = measure.regionprops(labels[0])    
        for ob in properties:
            if ob.area < LargeObjRemovalMinObjSize:
                mask[ob.coords[:,0],ob.coords[:,1]]=0

        mask = morphology.binary_dilation(mask, selem=morphology.disk(LargeObjRemovalSelem))
        mask = np.logical_not(mask)

        masked_img = img*mask

        
        return ((masked_img,img),metadata)



def large_beads_preprocessing(zarr_grp_name: str,
        parsed_raw_data_fpath: str,
        processing_parameters: dict,
        dark_img: np.ndarray)-> Tuple[Tuple[np.ndarray,],dict]:

    """Function used filter the reference image containing large beads 

    Args:
        zarr_grp_name (str): group name of the image to process
        parsed_raw_data_fpath (str): path to the zarr file containing the parsed images
        processing_parameters (dict): dictionary with the parameters used to process the images
        dark_img (np.ndarray): Dark image used to remove camera dark noise

    Returns:
        Tuple[Tuple[np.ndarray,],dict]: ((filtered_image,),metadata)
    """

    logger = selected_logger()
    
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']

    img, metadata = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)

    img = convert_from_uint16_to_float64(img)
    img -= dark_img
    img[img<0] = 0
    img = np.abs(img) # to avoid -0.0 issues

    img = img.max(axis=0)

    img /= filters.gaussian(img,FlatFieldKernel,preserve_range=False)


    return ((img,), metadata)


def both_beads_preprocessing(zarr_grp_name: str,
        parsed_raw_data_fpath: str,
        processing_parameters: dict,
        dark_img: np.ndarray)-> Tuple[Tuple[np.ndarray,],dict]:

    """Function used filter the reference image containing large beads 

    Args:
        zarr_grp_name (str): group name of the image to process
        parsed_raw_data_fpath (str): path to the zarr file containing the parsed images
        processing_parameters (dict): dictionary with the parameters used to process the images
        dark_img (np.ndarray): Dark image used to remove camera dark noise

    Returns:
        Tuple[Tuple[np.ndarray,],dict]: ((filtered_image,),metadata)
    """

    # logger = selected_logger()
    
    # parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    # experiment_fpath = parsed_raw_data_fpath.parent
    # FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']

    # img, metadata = load_raw_images(zarr_grp_name,
    #                         parsed_raw_data_fpath)

    # logger.info(f'loaded {zarr_grp_name} raw fish image')
    # img = convert_from_uint16_to_float64(img)
    # img -= dark_img
    # img[img<0] = 0
    # img = np.abs(img) # to avoid -0.0 issues

    # img = img.max(axis=0)

    # img /= filters.gaussian(img,FlatFieldKernel,preserve_range=False)


    # return ((img,),metadata)

    logger = selected_logger()
    
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']
    LaplacianKernel=processing_parameters['PreprocessingFishFilteringLaplacianKernel']
    

    try:
        img, metadata = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)
    except:
        logger.error(f'cannot load {zarr_grp_name} raw fish image')
        sys.exit(f'cannot load {zarr_grp_name} raw fish image')
    else:
        logger.info(f'loaded {zarr_grp_name} raw fish image')
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

        return ((img,),metadata)


def nuclei_registration_filtering(zarr_grp_name: str,
        parsed_raw_data_fpath: str,
        processing_parameters: dict)-> Tuple[Tuple[np.ndarray,],dict]:
    
    """This function remove the background from large structures like nuclei
        For the sigma I seleced a value quite bigger than
        the nuclei size in order to remove them from the 
        image. I used what showed on the gaussian filter code page and on this
        link on stackoverflow: 
        http://stackoverflow.com/questions/25216382/gaussian-filter-in-scipy

    Args:
        zarr_grp_name (str): group name of the image to process
        parsed_raw_data_fpath (str): path to the zarr file containing the parsed images
        processing_parameters (dict): dictionary with the parameters used to process the images

    Returns:
        Tuple[Tuple[np.ndarray,],dict]: ((filtered_image,),metadata)
    """
        

    logger = selected_logger()
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    experiment_fpath = parsed_raw_data_fpath.parent
    FlatFieldKernel=processing_parameters['PreprocessingNucleiFlatFieldKernel']

    try:
        img, metadata = load_raw_images(zarr_grp_name,
                                parsed_raw_data_fpath)
    except:
        logger.error(f'cannot load {zarr_grp_name} raw fish image')
        sys.exit(f'cannot load {zarr_grp_name} raw fish image')
    else:
        logger.info(f'loaded {zarr_grp_name} raw fish image')

        
        img = convert_from_uint16_to_float64(img)


        # Clean the image from the background
        img = img-filters.gaussian(img,sigma=FlatFieldKernel)
        # Remove the negative values        
        img[img<0] = 0
        # Flatten the image
        flattened_img = np.amax(img,axis=0)


        # img -= dark_img
        # # img -= filters.gaussian(img,FilteringSmallKernel,preserve_range=False)
        # img[img<0] = 0

        # background = filters.gaussian(img,FlatFieldKernel,preserve_range=False)
        # img -= background
        # img[img<=0] = 0 # All negative values set to zero also = to avoid -0.0 issues
        # img = np.abs(img) # to avoid -0.0 issues
        # img /= background
        
        # img = img.max(axis=0)
    
        return ((flattened_img,),metadata)



def fresh_nuclei_filtering(zarr_grp_name: str,
        parsed_raw_data_fpath: str,
        processing_parameters: dict)-> Tuple[Tuple[np.ndarray,],dict]:
        
        """This function remove the background and filter the nuclei in the
        low power magnification images used for identifying the cells position
        after running eel.
        For the sigma I seleced a value quite bigger than
        the nuclei size in order to remove them from the 
        image. I used what showed on the gaussian filter code page and on this
        link on stackoverflow: 
        http://stackoverflow.com/questions/25216382/gaussian-filter-in-scipy

    Args:
        zarr_grp_name (str): group name of the image to process
        parsed_raw_data_fpath (str): path to the zarr file containing the parsed images
        processing_parameters (dict): dictionary with the parameters used to process the images

    Returns:
        Tuple[Tuple[np.ndarray,],dict]: ((filtered_image,),metadata)
    """

        logger = selected_logger()
        parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
        experiment_fpath = parsed_raw_data_fpath.parent
        FlatFieldKernel = processing_parameters['PreprocessingFreshNucleiLargeKernelSize']
    
        try:
            img, metadata = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)
        except:
            logger.error(f'cannot load {zarr_grp_name} raw fish image')
            sys.exit(f'cannot load {zarr_grp_name} raw fish image')
        else:
            logger.info(f'loaded {zarr_grp_name} raw fish image')

            
            img = convert_from_uint16_to_float64(img)


            # Clean the image from the background
            img = img-filters.gaussian(img,sigma=FlatFieldKernel)
            # Remove the negative values        
            img[img<0] = 0
            # Flatten the image
            flattened_img = np.amax(img,axis=0)


            # img -= dark_img
            # # img -= filters.gaussian(img,FilteringSmallKernel,preserve_range=False)
            # img[img<0] = 0

            # background = filters.gaussian(img,FlatFieldKernel,preserve_range=False)
            # img -= background
            # img[img<=0] = 0 # All negative values set to zero also = to avoid -0.0 issues
            # img = np.abs(img) # to avoid -0.0 issues
            # img /= background
            
            # img = img.max(axis=0)
        
            return ((flattened_img,),metadata)


def fresh_tissue_beads_preprocessing(zarr_grp_name: str,
        parsed_raw_data_fpath: str,
        processing_parameters: dict,
        dark_img: np.ndarray)-> Tuple[Tuple[np.ndarray,],dict]:

    """Function used filter the reference image containing large beads 

    Args:
        zarr_grp_name (str): group name of the image to process
        parsed_raw_data_fpath (str): path to the zarr file containing the parsed images
        processing_parameters (dict): dictionary with the parameters used to process the images
        dark_img (np.ndarray): Dark image used to remove camera dark noise

    Returns:
        Tuple[Tuple[np.ndarray,],dict]: ((filtered_image,),metadata)
    """

    logger = selected_logger()
    
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']

    img, metadata = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)

    img = convert_from_uint16_to_float64(img)
    img = img.max(axis=0)

    img -= filters.gaussian(img,FlatFieldKernel,preserve_range=False)
    img[img<0] = 0

    return ((img,), metadata)

    