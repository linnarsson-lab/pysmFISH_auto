import numpy as np
import scipy.ndimage as nd
from skimage import filters, morphology, measure
from pathlib import Path

# pysmFISH imports
from pysmFISH.io import load_raw_images
from pysmFISH.utils import convert_from_uint16_to_float64

from flow_steps.filtering_counting import load_dark_image
from pysmFISH.logger_utils import selected_logger

def test_fun(zarr_grp_name,
        parsed_raw_data_fpath,
        processing_parameters):

    min_distance=processing_parameters['CountingFishMinObjDistance']
    min_obj_size=processing_parameters['CountingFishMinObjSize']
    max_obj_size=processing_parameters['CountingFishMaxObjSize']
    num_peaks_per_label=processing_parameters['CountingFishNumPeaksPerLabel']
    FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']


    logger = selected_logger()
    parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    experiment_fpath = parsed_raw_data_fpath.parent

    FlatFieldKernel=processing_parameters['PreprocessingFishFlatFieldKernel']
    raw_fish_images_meta = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)

    img = raw_fish_images_meta[0]
    img_metadata = raw_fish_images_meta[1]
    img = convert_from_uint16_to_float64(img)
    dark_img = load_dark_image(experiment_fpath)
    img -= dark_img
    img[img<0] = 0
    img = np.abs(img) # to avoid -0.0 issues

    img = img.max(axis=0)

    img /= filters.gaussian(img,FlatFieldKernel,preserve_range=False)

    return img