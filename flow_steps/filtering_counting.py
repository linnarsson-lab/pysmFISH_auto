""" 
Collection of functions to run filtering of the images and counting
"""
from typing import *
import pickle
import sys

import dask
import numpy as np
import scipy.ndimage as nd
from skimage import filters
from pathlib import Path
from dask.distributed import Client

# pysmFISH imports
from pysmFISH.io import load_raw_images
from pysmFISH.dots_calling import osmFISH_peak_based_detection
from pysmFISH.utils import convert_from_uint16_to_float64

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
        err = sys.exit(f'cannot load the dark image file')
    else:
        logger.info(f'loaded {dark_img_fpath[0].stem} dark image')
        dark_img = convert_from_uint16_to_float64(dark_img)
        return dark_img



def single_fish_filter_count_standard(
        zarr_grp_name,
        parsed_raw_data_fpath,
        FlatFieldKernel,
        FilteringSmallKernel, 
        LaplacianKernel,
        min_distance,
        min_obj_size,
        max_obj_size,
        num_peaks_per_label):


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
        FlatFieldKernel: np.ndarray
            size of the kernel use to remove the backgroun (ex. [10,10])
        FilteringSmallKernel: np.ndarray
            size of the kernel use to smooth the dots (ex. [8,8])
        LaplacianKernel: np.ndarray
            size of the kernel use enhance dots signal (ex. [1,1])
        min_distance: int
            minimum distance between two dots
        min_obj_size: int
            minimum size of a dot
        max_obj_size: int
            max size of a dot or a cluster of dots (depending on num_peaks_per_label)
        num_peaks_per_label
            max number of peaks called in a masked object
    
    """

    logger = selected_logger()
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    experiment_fpath = parsed_raw_data_fpath.parent
    
    try:
        raw_fish_images_meta = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)
    except:
        logger.error(f'cannot load {zarr_grp_name} raw fish image')
        sys.exit(f'cannot load {zarr_grp_name} raw fish image')
    else:
        logger.info(f'loaded {zarr_grp_name} raw fish image')
        try:
            # This may chnaged if the image will be store in shoji
            dark_img = load_dark_image(experiment_fpath)
        except:
            logger.error(f'cannot load dark reference fish image')
            sys.exit(f'cannot load dark reference fish image')
        else:
            logger.info('loaded dark reference image')

            img = raw_fish_images_meta[0]
            img_metadata = raw_fish_images_meta[1]
            img = convert_from_uint16_to_float64(img)
            img -= dark_img
            img = np.amax(img, axis=0)
            background = filters.gaussian(img,FlatFieldKernel,preserve_range=False)
            img -= filters.gaussian(img,FilteringSmallKernel,preserve_range=False)
            # img[img<0] = 0
            img /= background
            img = nd.gaussian_laplace(img,LaplacianKernel)
            # if np.all(img < 0):
            #     # This line was included to flip the peaks in the laplacian processing
            #     # I need to find the logic why it is not necessary anymore
            #     # The selection for np.all is not correct but I just want to keep it and see if
            #     # it works before remove it...find conditions when it will break
            #     logger.debug(f'image values are negative after laplacian. Values flipped')
            #     img = -img
            # img[img<0] = 0
            img = -img
            # img[img<0] = 0

            img = (img - np.mean(img)) / np.std(img)
            img[img<0] = 0 
        
            

            fish_counts = osmFISH_peak_based_detection((img, img_metadata),
                                                    min_distance,
                                                    min_obj_size,
                                                    max_obj_size,
                                                    num_peaks_per_label)
            
          
            fname = experiment_fpath / 'tmp' / (zarr_grp_name + '_filtered.pkl')
            pickle.dump((img, img_metadata),open(fname,'wb'))
            
            # save_dots_data(fish_counts)
            fname = experiment_fpath / 'tmp' / (zarr_grp_name + '_dots.pkl')
            pickle.dump(fish_counts,open(fname,'wb'))



def filtering_counting_runner(cluster,
                            running_function,
                            parsed_raw_data_fpath,
                            grp_name,
                            sorted_images_list,
                            processing_parameters):
    """
    Function used to run the filtering and counting of the images.
    the passed preprocessing function regulate the type of preprocessing 
    and counting defined
    
    Args:
        cluster: dask-obj 
            cluster used to run the processing
        running_function: python-obj
            filtering and couning function used to process the data
        parsed_raw_data_fpath: str
            path to the parsed file
        grp_name: str
            name of the grp to process (ex. fish, beads, staining)
        sorted_images_list: list
            list of the zarr grps of a specific type (ex. fish, beads, staining)
        processing_parameters: dict
            parameters used to define the processing conditions
    
    """

    assert len(sorted_images_list) > 0, sys.exit(f'missing grp of files to process')

    if grp_name == 'fish':
        FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']
        FilteringSmallKernel=processing_parameters['PreprocessingFishFilteringSmallKernel']
        LaplacianKernel=processing_parameters['PreprocessingFishFilteringLaplacianKernel']
        min_distance=processing_parameters['CountingFishMinObjDistance']
        min_obj_size=processing_parameters['CountingFishMinObjSize']
        max_obj_size=processing_parameters['CountingFishMaxObjSize']
        num_peaks_per_label=processing_parameters['CountingFishNumPeaksPerLabel']
    elif grp_name == 'beads':
        FlatFieldKernel=processing_parameters['PreprocessingBeadsRegistrationFlatFieldKernel']
        FilteringSmallKernel=processing_parameters['PreprocessingBeadsRegistrationFilteringSmallKernel']
        LaplacianKernel=processing_parameters['PreprocessingBeadsRegistrationFilteringLaplacianKernel']
        min_distance=processing_parameters['CountingBeadsRegistrationMinObjDistance']
        min_obj_size=processing_parameters['CountingBeadsRegistrationMinObjSize']
        max_obj_size=processing_parameters['CountingBeadsRegistrationMaxObjSize']
        num_peaks_per_label=processing_parameters['CountingBeadsRegistrationNumPeaksPerLabel']
    elif grp_name == 'staining':
        pass

    client = Client(cluster)
    L = client.map(running_function,sorted_images_list,
                                                parsed_raw_data_fpath=parsed_raw_data_fpath,
                                                FlatFieldKernel=FlatFieldKernel,
                                                FilteringSmallKernel=FilteringSmallKernel, 
                                                LaplacianKernel=LaplacianKernel,
                                                min_distance=min_distance,
                                                min_obj_size=min_obj_size,
                                                max_obj_size=max_obj_size,
                                                num_peaks_per_label=num_peaks_per_label)
    _ = client.gather(L)

    # all_futures = []
    # for processing_file in sorted_images_list:
    #     future = dask.delayed(running_function)(processing_file,
    #                                             parsed_raw_data_fpath=parsed_raw_data_fpath,
    #                                             FlatFieldKernel=FlatFieldKernel,
    #                                             FilteringSmallKernel=FilteringSmallKernel, 
    #                                             LaplacianKernel=LaplacianKernel,
    #                                             min_distance=min_distance,
    #                                             min_obj_size=min_obj_size,
    #                                             max_obj_size=max_obj_size,
    #                                             num_peaks_per_label=num_peaks_per_label)
    #     all_futures.append(future)
    
    # _ = dask.compute(*all_futures)