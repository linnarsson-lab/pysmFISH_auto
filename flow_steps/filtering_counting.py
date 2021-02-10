""" 
Collection of functions to run filtering of the images and counting
"""
from typing import *
import pickle
import sys

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
        err = sys.exit(f'cannot load the dark image file')
    else:
        logger.info(f'loaded {dark_img_fpath[0].stem} dark image')
        dark_img = convert_from_uint16_to_float64(dark_img)
        return dark_img



def single_fish_filter_count_standard(
        zarr_grp_name,
        parsed_raw_data_fpath,
        processing_parameters):

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
            img -= filters.gaussian(img,FilteringSmallKernel,preserve_range=False)
            img[img<0] = 0

            background = filters.gaussian(img,FlatFieldKernel,preserve_range=False)
            img /= background
            img = nd.gaussian_laplace(img,LaplacianKernel)
            img = -img # the peaks are negative so invert the signal
            img[img<0] = 0 # All negative values set to zero 

            img_mean_z = img.mean(axis=(1,2))
            img_mean_z = img_mean_z[:,np.newaxis,np.newaxis]
            img_std_z = img.std(axis=(1,2))
            img_std_z = img_std_z[:,np.newaxis,np.newaxis]
            img_nn= (img - img_mean_z)/ img_std_z
            img_nn = img_nn.max(axis=0)
            img_nn[img_nn<=0] = 0 # All negative values set to zero also = to avoid -0.0 issues


            fish_counts = osmFISH_peak_based_detection((img_nn, img_metadata),
                                                    min_distance,
                                                    min_obj_size,
                                                    max_obj_size,
                                                    num_peaks_per_label)
            
          
            fname = experiment_fpath / 'tmp' / 'filtered_images' / (zarr_grp_name + '_filtered.pkl')
            pickle.dump((img_nn, img_metadata),open(fname,'wb'))
            
            # save_dots_data(fish_counts)
            fname = experiment_fpath / 'tmp' / 'raw_counts' / (zarr_grp_name + '_dots.pkl')
            pickle.dump(fish_counts,open(fname,'wb'))


def single_fish_filter_count_standard_not_norm(
        zarr_grp_name,
        parsed_raw_data_fpath,
        processing_parameters):

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
            img -= filters.gaussian(img,FilteringSmallKernel,preserve_range=False)
            img[img<0] = 0

            background = filters.gaussian(img,FlatFieldKernel,preserve_range=False)
            img /= background
            img = nd.gaussian_laplace(img,LaplacianKernel)
            img = -img # the peaks are negative so invert the signal
            img[img<=0] = 0 # All negative values set to zero also = to avoid -0.0 issues

            img = img.max(axis=0)

            fish_counts = osmFISH_peak_based_detection((img, img_metadata),
                                                    min_distance,
                                                    min_obj_size,
                                                    max_obj_size,
                                                    num_peaks_per_label)
            
          
            fname = experiment_fpath / 'tmp' / 'filtered_images' / (zarr_grp_name + '_filtered.pkl')
            pickle.dump((img, img_metadata),open(fname,'wb'))
            

            # save_dots_data(fish_counts)
            fname = experiment_fpath / 'tmp' / 'raw_counts' / (zarr_grp_name + '_dots.pkl')
            pickle.dump(fish_counts,open(fname,'wb'))


def single_fish_filter_count_avoid_large_obj(
        zarr_grp_name,
        parsed_raw_data_fpath,
        processing_parameters):

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
    min_distance=processing_parameters['CountingFishMinObjDistance']
    min_obj_size=processing_parameters['CountingFishMinObjSize']
    max_obj_size=processing_parameters['CountingFishMaxObjSize']
    num_peaks_per_label=processing_parameters['CountingFishNumPeaksPerLabel']

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
            img -= filters.gaussian(img,FilteringSmallKernel,preserve_range=False)
            img[img<0] = 0

            background = filters.gaussian(img,FlatFieldKernel,preserve_range=False)
            img /= background
            img = nd.gaussian_laplace(img,LaplacianKernel)
            img = -img # the peaks are negative so invert the signal
            img[img<=0] = 0 # All negative values set to zero also = to avoid -0.0 issues

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


            fish_counts = osmFISH_barcoded_peak_based_detection_masked_thr((img, img_metadata),
                                                    masked_img,
                                                    min_distance,
                                                    min_obj_size,
                                                    max_obj_size,
                                                    num_peaks_per_label)
            
          
            fname = experiment_fpath / 'tmp' / 'filtered_images' / (zarr_grp_name + '_filtered.pkl')
            pickle.dump((img, img_metadata),open(fname,'wb'))
            

            # save_dots_data(fish_counts)
            fname = experiment_fpath / 'tmp' / 'raw_counts' / (zarr_grp_name + '_dots.pkl')
            pickle.dump(fish_counts,open(fname,'wb'))




def filtering_counting_large_beads(zarr_grp_name,
        parsed_raw_data_fpath,
        processing_parameters):

    logger = selected_logger()
    
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    experiment_fpath = parsed_raw_data_fpath.parent
    min_distance=processing_parameters['CountingFishMinObjDistance']
    min_obj_size=processing_parameters['CountingFishMinObjSize']
    max_obj_size=processing_parameters['CountingFishMaxObjSize']
    num_peaks_per_label=processing_parameters['CountingFishNumPeaksPerLabel']


    img = raw_fish_images_meta[0]
    img_metadata = raw_fish_images_meta[1]
    img = convert_from_uint16_to_float64(img)

    img -= dark_img
    img[img<0] = 0

    img = img.max(axis=0)

    fish_counts = osmFISH_barcoded_peak_based_detection_masked_thr((img, img_metadata),
                                                    masked_img,
                                                    min_distance,
                                                    min_obj_size,
                                                    max_obj_size,
                                                    num_peaks_per_label)
            
          
    fname = experiment_fpath / 'tmp' / 'filtered_images' / (zarr_grp_name + '_filtered.pkl')
    pickle.dump((img, img_metadata),open(fname,'wb'))
    

    # save_dots_data(fish_counts)
    fname = experiment_fpath / 'tmp' / 'raw_counts' / (zarr_grp_name + '_dots.pkl')
    pickle.dump(fish_counts,open(fname,'wb'))


def both_beads_filt_count_mask(
        zarr_grp_name,
        parsed_raw_data_fpath,
        processing_parameters):

    """
    Function to:
    - preprocess the fish images
    - count the dots and save the data
    The processing try to catch small beads in presence of large beads.
    The large beads are masked before calculating the thr for the counting
    but the counting is run in the non masked image
    
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
    
    data_models = Output_models()
    counts_dict = data_models.dots_counts_dict
    fill_value = np.nan

    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    experiment_fpath = parsed_raw_data_fpath.parent
    FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']
    FilteringSmallKernel=processing_parameters['PreprocessingFishFilteringSmallKernel']
    LaplacianKernel=processing_parameters['PreprocessingFishFilteringLaplacianKernel']
    min_distance=processing_parameters['CountingFishMinObjDistance']
    min_obj_size=processing_parameters['CountingFishMinObjSize']
    max_obj_size=processing_parameters['CountingFishMaxObjSize']
    num_peaks_per_label=processing_parameters['CountingFishNumPeaksPerLabel']

    LargeObjRemovalPercentile = processing_parameters['LargeObjRemovalPercentile']
    LargeObjRemovalMinObjSize = processing_parameters['LargeObjRemovalMinObjSize']
    LargeObjRemovalSelem = processing_parameters['LargeObjRemovalSelem']

    parameters_dict = {
    'min_distance': min_distance,
    'min_obj_size': min_obj_size,
    'max_obj_size': max_obj_size,
    'num_peaks_per_label':num_peaks_per_label}

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
            hybridization_num = img_metadata['hybridization_num']

            fov = img_metadata['fov_num']
            # Initialise an empty version of the counts dict
            counts_dict['r_px_original'] = np.array([fill_value])
            counts_dict['c_px_original'] = np.array([fill_value])
            counts_dict['dot_id'] = np.array([fill_value])
            counts_dict['fov_num'] = np.array(fov)
            counts_dict['round_num'] = np.array([img_metadata['hybridization_num']])
            counts_dict['dot_intensity'] = np.array([fill_value])
            counts_dict['selected_thr'] = np.array([fill_value])
            counts_dict['dot_channel'] = np.array([img_metadata['channel']])
            counts_dict['target_name'] = np.array([img_metadata['target_name']])
        
            img = convert_from_uint16_to_float64(img)
            img -= dark_img
            img = np.amax(img, axis=0)
            background = filters.gaussian(img,FlatFieldKernel,preserve_range=False)
            img /= background
            img -= background
            img = nd.gaussian_laplace(img,LaplacianKernel)
            img = -img
            
            mask = np.zeros_like(img)
            idx=  img > np.percentile(img,LargeObjRemovalPercentile)
            mask[idx] = 1
            # mask[~idx] = 0
            labels = nd.label(mask)

            properties = measure.regionprops(labels[0])    
            for ob in properties:
                if ob.area > LargeObjRemovalMinObjSize:
                    mask[ob.coords[:,0],ob.coords[:,1]]=0


            mask = morphology.binary_dilation(np.logical_not(mask), selem=morphology.disk(LargeObjRemovalSelem))
            masked_image = mask*img
            
            thr_calculation = osmFISH_dots_thr_selection(img,parameters_dict)
            thr_calculation.counting_graph()
            thr_calculation.thr_identification()

            if not np.isnan(thr_calculation.selected_thr):
                dots = osmFISH_dots_mapping(masked_image,thr_calculation.selected_thr,parameters_dict)
                if isinstance(dots.selected_peaks,np.ndarray):
                    # Peaks have been identified
                    total_dots = dots.selected_peaks.shape[0]
                    dot_id_array = np.array([str(fov)+'_'+str(hybridization_num)+'_'+ img_metadata['channel'] +'_'+str(nid) for nid in range(total_dots)])
                    fov_array = np.repeat(fov,total_dots)
                    thr_array = np.repeat(thr_calculation.selected_thr,total_dots)
                    channel_array = np.repeat(img_metadata['channel'],total_dots)
                    hybridization_num_array = np.repeat(img_metadata['hybridization_num'],total_dots)
                    target_name_array = np.repeat(img_metadata['target_name'],total_dots)

                    counts_dict['r_px_original'] = dots.selected_peaks[:,0]
                    counts_dict['c_px_original'] = dots.selected_peaks[:,1]
                    counts_dict['dot_id'] = dot_id_array
                    counts_dict['fov_num'] = fov_array
                    counts_dict['round_num'] = hybridization_num_array
                    counts_dict['dot_intensity'] = dots.intensity_array
                    counts_dict['selected_thr'] = thr_array
                    counts_dict['dot_channel'] = channel_array
                    counts_dict['target_name'] = target_name_array
                else:
                    logger.info(f' fov {fov} does not have counts (mapping)')
                        
            else:
                logger.info(f' fov {fov} does not have counts (thr)')
        
        fname = experiment_fpath / 'tmp' / 'filtered_images' / (zarr_grp_name + '_filtered.pkl')
        pickle.dump((img, img_metadata),open(fname,'wb'))
        
        # save_dots_data(fish_counts)
        fname = experiment_fpath / 'tmp' / 'raw_counts' / (zarr_grp_name + '_dots.pkl')
        pickle.dump((counts_dict,img_metadata),open(fname,'wb'))




def filtering_counting_runner(client,
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
        client: dask-client
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