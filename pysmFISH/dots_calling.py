from typing import *
import numpy as np
import argparse
from toolz.itertoolz import get
import zarr
import re
import sys
import logging
import pickle
import pandas as pd
from sympy import Point, Line
from skimage import feature, measure, morphology, img_as_float
from skimage.filters import rank_order
from scipy import ndimage as nd
from pathlib import Path

from pysmFISH.utils import convert_from_uint16_to_float64
from pysmFISH.data_models import Output_models

from pysmFISH.logger_utils import selected_logger


class osmFISH_dots_thr_selection():
    """
    Class used to automatically define the threshold used to call the
    signal peaks. This function calculate the threshold without masking large object 
    and contamination. 
    
    This is the original class used in the osmFISH paper.
    """
    
    
    def __init__(self, img:np.ndarray, parameters_dict:Dict, min_int:float=False, max_int:float=False,min_peaks:int=False):
        """Initialize the class

        Args:
            img (np.ndarray): Image to process
            parameters_dict (Dict): Parameters used to define the peaks.
            min_int (float, optional): Minimum intensity value to use for the binning of
                    the signal intensities. Defaults to False.
            max_int (float, optional): Maximum intensity value to use for the binning of
                    the signal intensities. Defaults to False.
            min_peaks (int, optional): Minimum number of peaks required for the
                    calculation of the counting threshold. Defaults to False.
        """
        self.img = img
        self.parameters_dict = parameters_dict
        self.min_int = min_int
        self.max_int = max_int
        self.min_peaks = min_peaks

        if self.min_peaks == False:
            self.min_peaks = 3
        
        self.min_distance = self.parameters_dict['min_distance']
        
        self.fill_value = np.nan

        # List with the total peaks calculated for each threshold
        self.total_peaks = []
        self.thr_used = []

    def counting_graph(self):
        """Function used for the construction of the number of peaks(Y) / thresholds(X)
        graph used to define the threshold.
        """

        binning = 100
        # Define the range of thr to be tested
        if self.min_int and self.max_int:
            self.thr_array = np.linspace(self.min_int,self.max_int,num=binning)
        elif self.min_int:
            self.thr_array = np.linspace(self.min_int,self.img.max(),num=binning)
        elif self.max_int:
            self.thr_array = np.linspace(np.min(self.img[np.nonzero(self.img)]),self.max_int,num=binning)
        else:
            self.thr_array = np.linspace(np.min(self.img[np.nonzero(self.img)]),self.img.max(),num=binning)
    
        # Calculate the number of peaks for each threshold. In this calculation
        # the size of the objects is not considered
        self.peak_counter_min = 0
        self.peak_counter_max = 0
        for vl, thr in enumerate(self.thr_array):
            # The border is excluded from the counting
            self.peaks = feature.peak_local_max(self.img,min_distance=self.min_distance,\
                threshold_abs=thr,exclude_border=False, indices=True,\
                num_peaks=np.inf, footprint=None,labels=None)
            
            self.number_peaks = len(self.peaks)

            # Stop the counting when the number of peaks detected falls below 3
            if self.number_peaks<=self.min_peaks:
                self.stop_thr = thr # Move in the upper loop so you will stop at the previous thr
                break
            else:
                self.total_peaks.append(len(self.peaks))
                self.thr_used.append(thr)

    def thr_identification(self):  
        """Function that use the number of peaks / thresholds graph to define the threshold
        to used for the counting.
        - calculate the gradient of the number of peaks / threshold function
        - remove the initial minimum point
        - calculate the segment that join the extremities of the gradient. This version
          of the code uses sympy.
        - Calculate the thr corresponding to the point of max distance from the segment
        """
        # Consider the case of no detectected peaks or if there is only one Thr
        # that create peaks (list total_peaks have only one element and )
        # if np.array(total_peaks).sum()>0 or len(total_peaks)>1:
        if len(self.total_peaks)>1:

            # Trim the threshold array in order to match the stopping point
            # used the [0][0] to get the first number and then take it out from list
            # thr_array = thr_array[:np.where(thr_array==stop_thr)[0][0]]
            self.thr_array = np.array(self.thr_used)

            # Calculate the gradient of the number of peaks distribution
            grad = np.gradient(self.total_peaks)
            
            # Restructure the data in order to avoid to consider the min_peak in the
            # calculations

            # Coord of the gradient min_peak
            grad_min_peak_coord = np.argmin(grad)
            
            # Trim the data to remove the peak.
            self.trimmed_thr_array = self.thr_array[grad_min_peak_coord:]
            self.trimmed_grad = grad[grad_min_peak_coord:]

            if self.trimmed_thr_array.shape>(1,):

                # Trim the coords array in order to maintain the same length of the 
                # tr and pk
                self.trimmed_total_peaks = self.total_peaks[grad_min_peak_coord:]

                # To determine the threshold we will determine the Thr with the biggest
                # distance to the segment that join the end points of the calculated
                # gradient

                # Distances list
                distances = []

                # Calculate the coords of the end points of the gradient
                p1 = Point(self.trimmed_thr_array[0],self.trimmed_grad[0])
                p2 = Point(self.trimmed_thr_array[-1],self.trimmed_grad[-1])
                
                # Create a line that join the points
                s = Line(p1,p2)
                allpoints = np.arange(0,len(self.trimmed_thr_array))
                
                # Calculate the distance between all points and the line
                for p in allpoints:
                    dst = s.distance(Point(self.trimmed_thr_array[p],self.trimmed_grad[p]))
                    distances.append(dst.evalf())

                # Remove the end points from the lists
                self.trimmed_thr_array = self.trimmed_thr_array[1:-1]
                self.trimmed_grad = self.trimmed_grad[1:-1]
                self.trimmed_total_peaks = self.trimmed_total_peaks[1:-1]
                self.trimmed_distances = distances[1:-1]
            
                # Determine the coords of the selected Thr
                # Converted trimmed_distances to array because it crashed
                # on Sanger.
                if self.trimmed_distances: # Most efficient way will be to consider the length of Thr list
                    thr_idx = np.argmax(np.array(self.trimmed_distances))
                    self.selected_thr = self.trimmed_thr_array[thr_idx]
                    # The selected threshold usually causes oversampling of the number of dots
                    # I added a stringency parameter (int n) to use to select the Thr+n 
                    # for the counting. It selects a stringency only if the trimmed_thr_array
                    # is long enough. Also consider the case in which the stringency in negative
                else:
                    self.selected_thr = self.fill_value
                    self.trimmed_thr_array = self.fill_value
            else:
                self.selected_thr = self.fill_value
                self.trimmed_thr_array = self.fill_value
        else:
            self.selected_thr = self.fill_value
            self.trimmed_thr_array = self.fill_value



class osmFISH_dots_mapping():
    """Function used to count the peaks after identification of the threshold
    and masking of large objects.

    This is the original class used in the osmFISH paper. 
    """

    def __init__(self,img: np.ndarray,thr: float,parameters_dict: dict):
        """Class initialization

        Args:
            img (np.ndarray): Image to process
            thr (float): Precalculate threshold for masking the image 
            parameters_dict (dict): Parameters used to define the peaks.
        """
       # Calculate the selected peaks after removal of the big and small objects
        
        self.img = img
        self.thr = thr
        # make an error if selected Thr <0
        self.parameters_dict = parameters_dict
        
        self.min_distance = self.parameters_dict['min_distance']
        self.min_obj_size = self.parameters_dict['min_obj_size']
        self.max_obj_size = self.parameters_dict['max_obj_size']
        self.num_peaks_per_label = self.parameters_dict['num_peaks_per_label']

        self.fill_value = np.nan

        # Threshold the image using the selected threshold
        img_mask = self.img>self.thr
    
        labels = nd.label(img_mask)[0]
        
        properties = measure.regionprops(labels)
            
        for ob in properties:
            if ob.area<self.min_obj_size or ob.area>self.max_obj_size:
                img_mask[ob.coords[:,0],ob.coords[:,1]]=0
        
        labels = nd.label(img_mask)[0]

        # Collect the properties of the labels after size selection
        properties = measure.regionprops(labels,intensity_image=self.img)

        self.selected_peaks = feature.peak_local_max(self.img, min_distance=self.min_distance, 
                                threshold_abs=self.thr, exclude_border=True, 
                                footprint=None, labels=labels,num_peaks_per_label=self.num_peaks_per_label)                            
        
        
        # # calling peak_local_max without Labels argument
        # selected_peaks_mask = feature.peak_local_max(self.img, min_distance=self.min_distance, 
        #                 threshold_abs=self.thr, exclude_border=True,
        #                 footprint=None,num_peaks=np.inf,indices=False).astype(int)                         
                              
        # # instead, make sure the selected peaks does not meet zeros at labels (background)
        # labels_mask = (labels > 0).astype(int)
        # selected_peaks_mask = selected_peaks_mask * labels_mask
        # self.selected_peaks = np.vstack(np.where(selected_peaks_mask)).T
        
        
        if self.selected_peaks.size:
            self.intensity_array = self.img[self.selected_peaks[:,0],self.selected_peaks[:,1]]
        else:
            self.intensity_array = np.nan


def peak_thrs_local_max_fast(image: np.ndarray, min_distance: int=1, 
                        threshold_abs: float=None,threshold_rel:float=None, 
                        exclude_border: int=True, indices: bool=True,
                        num_peaks: int=np.inf, footprint: np.ndarray=None, 
                        labels: np.ndarray=None)->np.ndarray:
    """Function after modification:
    returns the coordinates for a range of thresholds

    Peaks are the local maxima in a region of `2 * min_distance + 1`
    (i.e. peaks are separated by at least `min_distance`).

    If peaks are flat (i.e. multiple adjacent pixels have identical
    intensities), the coordinates of all such pixels are returned.

    If both `threshold_abs` and `threshold_rel` are provided, the maximum
    of the two is chosen as the minimum intensity threshold of peaks.

    Notes
    -----
    The peak local maximum function returns the coordinates of local peaks
    (maxima) in an image. A maximum filter is used for finding local maxima.
    This operation dilates the original image. After comparison of the dilated
    and original image, this function returns the coordinates or a mask of the
    peaks where the dilated image equals the original image.

    Examples
    --------
    >>> img1 = np.zeros((7, 7))
    >>> img1[3, 4] = 1
    >>> img1[3, 2] = 1.5
    >>> img1
    array([[ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
           [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
           [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
           [ 0. ,  0. ,  1.5,  0. ,  1. ,  0. ,  0. ],
           [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
           [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
           [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ]])

    >>> peak_local_max(img1, min_distance=1)
    array([[3, 2],
           [3, 4]])

    >>> peak_local_max(img1, min_distance=2)
    array([[3, 2]])

    >>> img2 = np.zeros((20, 20, 20))
    >>> img2[10, 10, 10] = 1
    >>> peak_local_max(img2, exclude_border=0)
    array([[10, 10, 10]])

    Args:
        image (np.ndarray): Input image.
        min_distance (int, optional): Minimum number of pixels separating peaks in a region of `2 *
                    min_distance + 1` (i.e. peaks are separated by at least
                    `min_distance`). Defaults to 1.
        threshold_abs (float, optional): Minimum intensity of peaks. By default, the absolute threshold is
                    the minimum intensity of the image. Defaults to None.
        threshold_rel (float, optional): Minimum intensity of peaks, calculated as `max(image) * threshold_rel`. 
                    Defaults to None.
        exclude_border (int, optional): If nonzero, `exclude_border` excludes peaks from
                    within `exclude_border`-pixels of the border of the image.. Defaults to True.
        indices (bool, optional): If True, the output will be an array representing peak
                    coordinates.  If False, the output will be a boolean array shaped as
                    `image.shape` with peaks present at True elements.. Defaults to True.
        num_peaks (int, optional): Maximum number of peaks. When the number of peaks exceeds `num_peaks`,
                    return `num_peaks` peaks based on highest peak intensity.. Defaults to np.inf.
        footprint (np.ndarray, optional): If provided, `footprint == 1` represents the local region within which
                    to search for peaks at every point in `image`.  Overrides
                    `min_distance` (also for `exclude_border`).. Defaults to None.
        labels (np.ndarray, optional): If provided, each unique region `labels == value` represents a unique
                    region to search for peaks. Zero is reserved for background.. Defaults to None.

    Returns:
        np.ndarray: If `indices = True`  : (row, column, ...) coordinates of peaks.
                    If `indices = False` : Boolean array shaped like `image`, with peaks
                    represented by True values.
    """

    if type(exclude_border) == bool:
        exclude_border = min_distance if exclude_border else 0

    out = np.zeros_like(image, dtype=np.bool)

    # In the case of labels, recursively build and return an output
    # operating on each label separately
    if labels is not None:
        label_values = np.unique(labels)
        # Reorder label values to have consecutive integers (no gaps)
        if np.any(np.diff(label_values) != 1):
            mask = labels >= 1
            labels[mask] = 1 + rank_order(labels[mask])[0].astype(labels.dtype)
        labels = labels.astype(np.int32)

        # New values for new ordering
        label_values = np.unique(labels)
        for label in label_values[label_values != 0]:
            maskim = (labels == label)
            out += feature.peak_local_max(image * maskim, min_distance=min_distance,
                                  threshold_abs=threshold_abs,
                                  threshold_rel=threshold_rel,
                                  exclude_border=exclude_border,
                                  indices=False, num_peaks=np.inf,
                                  footprint=footprint, labels=None)

        if indices is True:
            return np.transpose(out.nonzero())
        else:
            return out.astype(np.bool)

    if np.all(image == image.flat[0]):
        if indices is True:
            return np.empty((0, 2), np.int)
        else:
            return out

    # Non maximum filter
    if footprint is not None:
        image_max = nd.maximum_filter(image, footprint=footprint,
                                       mode='constant')
    else:
        size = 2 * min_distance + 1
        image_max = nd.maximum_filter(image, size=size, mode='constant')
    mask = image == image_max

    if exclude_border:
        # zero out the image borders
        for i in range(mask.ndim):
            mask = mask.swapaxes(0, i)
            remove = (footprint.shape[i] if footprint is not None
                      else 2 * exclude_border)
            mask[:remove // 2] = mask[-remove // 2:] = False
            mask = mask.swapaxes(0, i)

    # find top peak candidates above a threshold
    thresholds = []
    if threshold_abs is None:
        threshold_abs = image.min()
    thresholds.append(threshold_abs)
    if threshold_rel is not None:
        thresholds.append(threshold_rel * image.max())
    if thresholds:
        mask_original = mask  # save the local maxima's of the image
        thrs_coords = {}  # dictionary holds the coordinates correspond for each threshold
        for threshold in thresholds[0]:
            mask = mask_original
            mask &= image > threshold

            # get coordinates of peaks
            coordinates = np.transpose(mask.nonzero())

            if coordinates.shape[0] > num_peaks:
                intensities = image.flat[np.ravel_multi_index(coordinates.transpose(),
                                                              image.shape)]
                idx_maxsort = np.argsort(intensities)[::-1]
                coordinates = coordinates[idx_maxsort][:num_peaks]

            if indices is True:
                thrs_coords[threshold] = coordinates
            else:
                nd_indices = tuple(coordinates.T)
                out[nd_indices] = True
                return out
    if thresholds and thrs_coords:
        return thrs_coords


def osmFISH_peak_based_detection_fast(ImgStack: np.ndarray, 
                            fov_subdataset: pd.Series,
                            parameters_dict: dict,
                            dimensions: int=2,
                            stringency:int =0,
                            min_int:float=False,     
                            max_int:float=False,
                            min_peaks:int=False)->pd.DataFrame:
    
    """This function is used to calculate the threshold to use for the dots
    counting in a 3D image. It is based on the function used in the osmFISH
    paper but doesn’t require simpy. It integrate the peak_thrs_local_max_fast
    and the calculation of the peaks on the masked image.

    Args:
        ImgStack (np.ndarray): preprocessed image used to count the dots
        fov_subdataset (pd.Series): Series with the metadata info relative to the image to process
        parameters_dict (dict): Parameters used to define the peaks.
        dimensions (int, optional): Image dimension (2 for 2D or 3 for 3D). Defaults to 2.
        stringency (int, optional): Select a thr before or after the one calculated
                    automatically. Defaults to 0.
        min_int (float, optional): Minimum intensity value to use for the binning of
                    the signal intensities. Defaults to False.
        max_int (float, optional): Maximum intensity value to use for the binning of
                    the signal intensities. Defaults to False.
        min_peaks (int, optional): Minimum number of peaks required for the
                    calculation of the counting threshold. Defaults to False.

    Returns:
        pd.DataFrame: counts data
    """



    
    logger = selected_logger()
    
    
    if min_peaks == False:
        min_peaks = 3

    fill_value = np.nan

    # List with the total peaks calculated for each threshold
    thr_used = []

    binning = 100
    # Define the range of thr to be tested
    if min_int and max_int:
        ThrArray = np.linspace(min_int,max_int,num=binning)
    elif min_int:
        ThrArray = np.linspace(min_int,ImgStack.max(),num=binning)
    elif max_int:
        ThrArray = np.linspace(np.min(ImgStack[np.nonzero(ImgStack)]),max_int,num=binning)
    else:
        ThrArray = np.linspace(np.min(ImgStack[np.nonzero(ImgStack)]),ImgStack.max(),num=binning)
    

    fov = fov_subdataset.fov_num
    round_num = fov_subdataset.round_num
    channel = fov_subdataset.channel
    target_name = fov_subdataset.target_name
    
    fill_value = np.nan
    data_models = Output_models()
    counts_dict = data_models.dots_counts_dict
    
    # Initialise an empty version of the counts dict
    counts_dict['r_px_original'] = np.array([fill_value])
    counts_dict['c_px_original'] = np.array([fill_value])
    counts_dict['dot_id'] = np.array([fill_value])
    counts_dict['dot_intensity'] = np.array([fill_value])
    counts_dict['selected_thr'] = np.array([fill_value])
    
    min_distance = parameters_dict['CountingFishMinObjDistance']
    min_obj_size = parameters_dict['CountingFishMinObjSize']
    max_obj_size = parameters_dict['CountingFishMaxObjSize']
    #num_peaks_per_label = self.parameters_dict['num_peaks_per_label']
    
    fov_subdataset_df = pd.DataFrame(fov_subdataset).T
    
    # List of ndarrays with the coords of the peaks calculated for each threshold
    PeaksCoords = []
    Selected_Peaks2_mask = None
    # Determine if working in 2D or 3D 
    if dimensions == 2:
        if len(ImgStack.shape) > 2:
            ImgStack = np.amax(ImgStack, axis=0)

    # Define the Thr array
#     ThrArray = np.linspace(ImgStack.min(), ImgStack.max(), num=binning)


    # Calculate the number of peaks for each threshold
    # Exclude border beacause of artefact of image processing
    thrs_peaks = peak_thrs_local_max_fast(ImgStack, min_distance=min_distance,
                                     threshold_abs=ThrArray, exclude_border=True, indices=True,
                                     num_peaks=np.inf, footprint=None, labels=None)
    lists = sorted(thrs_peaks.items())  # sorted by key, return a list of tuples. tuple[0]: threshold, tuple[1]: coords
    x, PeaksCoords = zip(*lists)  # unpack a list of pairs into two tuples
    TotalPeaks = []
    for j in range(len(PeaksCoords)):
        TotalPeaks += (len(PeaksCoords[j]),)  # get number of peaks
    # print("Thresholds distribution %.3f seconds" % (timings['thrs_dist']))

    # Consider the case of no detectected peaks or if there is only one Thr
    # that create peaks (list TotalPeaks have only one element and )
    # if np.array(TotalPeaks).sum()>0 or len(TotalPeaks)>1:
    if len(TotalPeaks) > 3:

        # Trim the threshold array in order to match the stopping point
        # used the [0][0] to get the first number and then take it out from list
        # ThrArray = ThrArray[:np.where(ThrArray == StopThr)[0][0]]

        # Trim and convert to types as Simone's
        TotalPeaks = np.array(TotalPeaks)
        TotalPeaks = list(TotalPeaks[TotalPeaks > min_peaks])

        ThrArray = ThrArray[:len(TotalPeaks)]

        PeaksCoords = np.array(PeaksCoords)
        PeaksCoords = PeaksCoords[:len(TotalPeaks)]
        PeaksCoords = list(PeaksCoords)
        
        if len(TotalPeaks) > 3:
            # Calculate the gradient of the number of peaks distribution
            # grad = np.gradient(TotalPeaks)
            grad = np.gradient(TotalPeaks,edge_order=1)

            # Restructure the data in order to avoid to consider the min_peak in the
            # calculations

            # Coord of the gradient min_peak
            grad_min_peak_coord = np.argmin(grad)

            # Trim the data to remove the peak.
            trimmed_thr_array = ThrArray[grad_min_peak_coord:]
            trimmed_grad = grad[grad_min_peak_coord:]

            if trimmed_thr_array.shape > (1,):

                # Trim the coords array in order to maintain the same length of the 
                # tr and pk
                Trimmed_PeaksCoords = PeaksCoords[grad_min_peak_coord:]
                trimmed_total_peaks = TotalPeaks[grad_min_peak_coord:]

                # To determine the threshold we will determine the Thr with the biggest
                # distance to the segment that join the end points of the calculated
                # # gradient

                # Calculate the coords of the end points of the gradient
                p1 = np.array([trimmed_thr_array[0],trimmed_grad[0]])
                p2 = np.array([trimmed_thr_array[-1],trimmed_grad[-1]])

                # Create a line that join the points
                allpoints = np.arange(0,len(trimmed_thr_array))
                allpoints_coords = np.array([trimmed_thr_array[allpoints],trimmed_grad[allpoints]]).T

                distances = []
                for point in allpoints_coords:
                    distances.append(np.linalg.norm(np.cross(p2-p1, p1-point))/np.linalg.norm(p2-p1))


                # Remove the end points from the lists
                trimmed_thr_array = trimmed_thr_array[1:-1]
                trimmed_grad = trimmed_grad[1:-1]
                trimmed_total_peaks = trimmed_total_peaks[1:-1]
                trimmed_distances = distances[1:-1]

                # Determine the coords of the selected Thr
                # Converted Trimmed_distances to array because it crashed
                # on Sanger.
                if trimmed_distances:  # Most efficient way will be to consider the length of Thr list
                    Thr_idx = np.argmax(np.array(trimmed_distances))
                    Calculated_Thr = trimmed_thr_array[Thr_idx]
                    # The selected threshold usually causes oversampling of the number of dots
                    # I added a stringency parameter (int n) to use to select the Thr+n 
                    # for the counting. It selects a stringency only if the Trimmed_ThrArray
                    # is long enough
                    if Thr_idx + stringency < len(trimmed_thr_array):
                        Selected_Thr = trimmed_thr_array[Thr_idx + stringency]
                        Selected_Peaks = Trimmed_PeaksCoords[Thr_idx + stringency]
                    else:
                        Selected_Thr = trimmed_thr_array[Thr_idx]
                        Selected_Peaks = Trimmed_PeaksCoords[Thr_idx]

                    # Calculate the selected peaks after removal of the big and small objects
                    # Threshold the image using the selected threshold
                    # if Selected_Thr > 0:
                    #     ImgMask = ImgStack > Selected_Thr

                    ImgMask = ImgStack > Selected_Thr


                    Labels = nd.label(ImgMask)[0]

                    Properties = measure.regionprops(Labels)

                    for ob in Properties:
                        if ob.area < min_obj_size or ob.area > max_obj_size:
                            ImgMask[ob.coords[:, 0], ob.coords[:, 1]] = 0

                    Labels = nd.label(ImgMask)[0]

                    # # calling peak_local_max without Labels argument
                    # Selected_Peaks2_mask = feature.peak_local_max(ImgStack, min_distance=min_distance, 
                    #                 threshold_abs=Selected_Thr, exclude_border=True, indices=False,
                    #                 footprint=None,num_peaks=np.inf).astype(int)                         
                    
                
                    # # instead, make sure the selected peaks does not meet zeros at labels (background)
                    # Labels_mask = (Labels > 0).astype(int)
                    # Selected_Peaks2_mask = Selected_Peaks2_mask * Labels_mask
                    # Selected_Peaks2 = np.vstack(np.where(Selected_Peaks2_mask)).T
                    
                    Selected_Peaks2 = feature.peak_local_max(ImgStack, min_distance=min_distance, 
                                    threshold_abs=Selected_Thr, exclude_border=True, indices=True,
                                    footprint=None,labels=Labels,num_peaks=np.inf).astype(int)        
    
                    
                    if Selected_Peaks2.size:
                        # Intensity counting of the max peaks
    #                     Selected_peaks_coords = np.where(Selected_Peaks2)
    #                     Selected_peaks_int = ImgStack[Selected_peaks_coords[0], Selected_peaks_coords[1]]
                        Selected_peaks_int = ImgStack[Selected_Peaks2[:, 0], Selected_Peaks2[:, 1]]
                                            
                        # Peaks have been identified
                        total_dots = Selected_Peaks2.shape[0]
                        dot_id_array = np.array([str(fov)+'_'+str(round_num)+'_'+ channel +'_'+str(nid) for nid in range(total_dots)])
                        thr_array = np.repeat(Selected_Thr,total_dots)
                        channel_array = np.repeat(channel,total_dots)

                        counts_dict['r_px_original']  = Selected_Peaks2[:,0]
                        counts_dict['c_px_original'] = Selected_Peaks2[:,1]
                        counts_dict['dot_id'] = dot_id_array
                        counts_dict['dot_intensity'] = Selected_peaks_int
                        counts_dict['selected_thr'] = thr_array
                    else:
                        logger.info(f' fov {fov} does not have counts (mapping)')

                else:
                    logger.info(f' fov {fov} Trimmed distance equal to zero')
            else:
                logger.info(f' fov {fov} calculated Thr array to small for selection of Thr')

        else:
            logger.info(f' fov {fov} does not have counts for calculating Thr')

    else:
        logger.info(f' fov {fov} does not have counts for calculating Thr')

        
        
    counts_df = pd.DataFrame(counts_dict)
    fov_subdataset_df = pd.DataFrame(fov_subdataset).T
    fov_subdataset_df = pd.concat([fov_subdataset_df]*counts_df.shape[0],axis=0).sort_index().reset_index(drop=True)
    counts_df = pd.concat([counts_df,fov_subdataset_df],axis=1)
    
    
    return counts_df




def beads_peak_based_detection(img: np.ndarray, 
                            fov_subdataset: pd.Series,
                            processing_parameters: Dict)->pd.DataFrame:
    """Counts the peaks in the reference images with small and large beads.
    It first identify the large beads and then mask them and identify the
    small ones.

    Args:
        img (np.ndarray): Reference image with large and small beads
        fov_subdataset (pd.Series): Series with the metadata info relative to the image to process
        processing_parameters (Dict): Parameters used to define the peaks.

    Returns:
        pd.DataFrame: Beads counts
    """

    stitching_type = fov_subdataset.stitching_type      

    LargeObjRemovalPercentile = processing_parameters['LargeObjRemovalPercentile']
    LargeObjRemovalMinObjSize = processing_parameters['LargeObjRemovalMinObjSize']
    LargeObjRemovalSelem = processing_parameters['LargeObjRemovalSelem']


    if stitching_type == 'both-beads':

        large_beads_counts_df = osmFISH_peak_based_detection_fast(img,fov_subdataset,processing_parameters)
        large_beads_counts_df['mapped_beads_type'] = 'large'
        

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

        processing_parameters_small = {
                'CountingFishMinObjDistance': 5,
                'CountingFishMaxObjSize': 20,
                'CountingFishMinObjSize': 2,
                'CountingFishNumPeaksPerLabel': 1}

        small_beads_counts_df = osmFISH_peak_based_detection_fast(masked_img,fov_subdataset,processing_parameters_small)
        small_beads_counts_df['mapped_beads_type'] = 'small'

        counts_df = pd.concat([large_beads_counts_df,small_beads_counts_df], axis=0, copy=False)  

    elif stitching_type == 'small-beads':
        counts_df = osmFISH_peak_based_detection_fast(img,fov_subdataset,processing_parameters)
        counts_df['mapped_beads_type'] = 'small'

    elif stitching_type == 'large-beads':
        counts_df = osmFISH_peak_based_detection_fast(img,fov_subdataset,processing_parameters)
        counts_df['mapped_beads_type'] = 'large'
        

    return counts_df


def osmFISH_peak_based_detection(ImgStack: np.ndarray,
                                        fov_subdataset: pd.Series,
                                        parameters_dict: dict):
    
    """
    This funtion apply the same peak based detection strategy used for 
    dots calling in the osmFISH paper
    
    Args:
    -----------
    img_meta: tuple
        tuple containing (image np.ndarray and metadata dict)
    min_distance: np.float64
        minimum distance between two peaks
    min_obj_size: np.uint16
        minimum object size of the objects that will be processed for peak detection
        objects below this value are discharged
    max_obj_size: np.uint16
        maximum object size of the objects that will be processed for peak detection
        objects above this value are discharged
    num_peaks_per_label: np.uint16
        Max number of peaks detected in each segmented object. Use None for max detection

    """

    
    
    logger = selected_logger()


    fov = fov_subdataset.fov_num
    round_num = fov_subdataset.round_num
    channel = fov_subdataset.channel
    target_name = fov_subdataset.target_name
    
    fill_value = np.nan
    data_models = Output_models()
    counts_dict = data_models.dots_counts_dict
    
    # Initialise an empty version of the counts dict
    counts_dict['r_px_original'] = np.array([fill_value])
    counts_dict['c_px_original'] = np.array([fill_value])
    counts_dict['dot_id'] = np.array([fill_value])
    counts_dict['dot_intensity'] = np.array([fill_value])
    counts_dict['selected_thr'] = np.array([fill_value])
    
    
    fov_subdataset_df = pd.DataFrame(fov_subdataset).T
    

    counting_parameters_dict = {
                            'min_distance': parameters_dict['CountingFishMinObjDistance'],
                            'min_obj_size': parameters_dict['CountingFishMinObjSize'],
                            'max_obj_size':  parameters_dict['CountingFishMaxObjSize'],
                            'num_peaks_per_label':  parameters_dict['CountingFishNumPeaksPerLabel'],
                                }
    fill_value = np.nan
    counts = osmFISH_dots_thr_selection(ImgStack,counting_parameters_dict)
    counts.counting_graph()
    counts.thr_identification()
    data_models = Output_models()
    counts_dict = data_models.dots_counts_dict

                    
    if not np.isnan(counts.selected_thr):
            dots = osmFISH_dots_mapping(ImgStack,counts.selected_thr,counting_parameters_dict)
            if isinstance(dots.selected_peaks,np.ndarray):
                # Peaks have been identified
                total_dots = dots.selected_peaks.shape[0]
                dot_id_array = np.array([str(fov)+'_'+str(round_num)+'_'+ channel +'_'+str(nid) for nid in range(total_dots)])
                thr_array = np.repeat(counts.selected_thr,total_dots)

                counts_dict['r_px_original']  = dots.selected_peaks[:,0]
                counts_dict['c_px_original'] = dots.selected_peaks[:,1]
                counts_dict['dot_id'] = dot_id_array
                counts_dict['dot_intensity'] = dots.intensity_array
                counts_dict['selected_thr'] = thr_array
                
            else:
                logger.info(f' fov {fov} does not have counts (mapping)')
                
    else:
        logger.info(f' fov {fov} does not have counts (thr)')
    
    
    counts_df = pd.DataFrame(counts_dict)
    fov_subdataset_df = pd.DataFrame(fov_subdataset).T
    fov_subdataset_df = pd.concat([fov_subdataset_df]*counts_df.shape[0],axis=0).sort_index().reset_index(drop=True)
    counts_df = pd.concat([counts_df,fov_subdataset_df],axis=1)
    
    
    return counts_df


# TODO remove unused functions below

# def osmFISH_peak_based_detection(img_meta:Tuple[np.ndarray, Dict],
#                                         min_distance: np.float64,
#                                         min_obj_size: np.uint16,
#                                         max_obj_size: np.uint16,
#                                         num_peaks_per_label: np.uint16):
    
#     """
#     This funtion apply the same peak based detection strategy used for 
#     dots calling in the osmFISH paper
    
#     Args:
#     -----------
#     img_meta: tuple
#         tuple containing (image np.ndarray and metadata dict)
#     min_distance: np.float64
#         minimum distance between two peaks
#     min_obj_size: np.uint16
#         minimum object size of the objects that will be processed for peak detection
#         objects below this value are discharged
#     max_obj_size: np.uint16
#         maximum object size of the objects that will be processed for peak detection
#         objects above this value are discharged
#     num_peaks_per_label: np.uint16
#         Max number of peaks detected in each segmented object. Use None for max detection

#     """

    
    
#     logger = selected_logger()

#     img = img_meta[0]
#     img_metadata = img_meta[1]
#     fov = img_metadata['fov_num']
#     hybridization = img_metadata['hybridization_num']
#     target_name = img_metadata['target_name']
    
#     logger.info(f'logging osmFISH_peak_based_detection fov {fov}')
#     hybridization_num = img_metadata['hybridization_num']

#     counting_parameters_dict = {
#                             'min_distance': min_distance,
#                             'min_obj_size': min_obj_size,
#                             'max_obj_size': max_obj_size,
#                             'num_peaks_per_label': num_peaks_per_label,
#                                 }
#     fill_value = np.nan
#     counts = osmFISH_dots_thr_selection(img,counting_parameters_dict)
#     counts.counting_graph()
#     counts.thr_identification()
#     data_models = Output_models()
#     counts_dict = data_models.dots_counts_dict

#     # Initialise an empty version of the counts dict
#     counts_dict['r_px_original'] = np.array([fill_value])
#     counts_dict['c_px_original'] = np.array([fill_value])
#     counts_dict['dot_id'] = np.array([fill_value])
#     counts_dict['fov_num'] = np.array(fov)
#     counts_dict['round_num'] = np.array([img_metadata['hybridization_num']])
#     counts_dict['dot_intensity'] = np.array([fill_value])
#     counts_dict['selected_thr'] = np.array([fill_value])
#     counts_dict['dot_channel'] = np.array([img_metadata['channel']])
#     counts_dict['target_name'] = np.array([img_metadata['target_name']])
                    
#     if not np.isnan(counts.selected_thr):
#             dots = osmFISH_dots_mapping(img,counts.selected_thr,counting_parameters_dict)
#             if isinstance(dots.selected_peaks,np.ndarray):
#                 # Peaks have been identified
#                 total_dots = dots.selected_peaks.shape[0]
#                 dot_id_array = np.array([str(fov)+'_'+str(hybridization_num)+'_'+ img_metadata['channel'] +'_'+str(nid) for nid in range(total_dots)])
#                 fov_array = np.repeat(fov,total_dots)
#                 thr_array = np.repeat(counts.selected_thr,total_dots)
#                 channel_array = np.repeat(img_metadata['channel'],total_dots)
#                 hybridization_num_array = np.repeat(img_metadata['hybridization_num'],total_dots)
#                 target_name_array = np.repeat(img_metadata['target_name'],total_dots)

#                 counts_dict['r_px_original']  = dots.selected_peaks[:,0]
#                 counts_dict['c_px_original'] = dots.selected_peaks[:,1]
#                 counts_dict['dot_id'] = dot_id_array
#                 counts_dict['fov_num'] = fov_array
#                 counts_dict['round_num'] = hybridization_num_array
#                 counts_dict['dot_intensity'] = dots.intensity_array
#                 counts_dict['selected_thr'] = thr_array
#                 counts_dict['dot_channel'] = channel_array
#                 counts_dict['target_name'] = target_name_array
#             else:
#                 logger.info(f' fov {fov} does not have counts (mapping)')
                
#     else:
#         logger.info(f' fov {fov} does not have counts (thr)')
    
#     return (counts_dict, img_metadata)




# def osmFISH_peak_based_detection_test(img:np.ndarray,
#                                     fov_subdataset,
#                                     processing_parameters:Dict):
    
#     """
#     This funtion apply the same peak based detection strategy used for 
#     dots calling in the osmFISH paper
    
#     Args:
#     -----------
#     img_meta: tuple
#         tuple containing (image np.ndarray and metadata dict)
#     min_distance: np.float64
#         minimum distance between two peaks
#     min_obj_size: np.uint16
#         minimum object size of the objects that will be processed for peak detection
#         objects below this value are discharged
#     max_obj_size: np.uint16
#         maximum object size of the objects that will be processed for peak detection
#         objects above this value are discharged
#     num_peaks_per_label: np.uint16
#         Max number of peaks detected in each segmented object. Use None for max detection

#     """

    
    
#     logger = selected_logger()

#     fov = fov_subdataset.fov_num
#     round_num = fov_subdataset.round_num
#     channel = fov_subdataset.channel
#     target_name = fov_subdataset.target_name

    
#     logger.info(f'logging osmFISH_peak_based_detection fov {fov}')

#     counting_parameters_dict = {
#                             'min_distance': processing_parameters['CountingFishMinObjDistance'],
#                             'min_obj_size': processing_parameters['CountingFishMinObjSize'],
#                             'max_obj_size': processing_parameters['CountingFishMaxObjSize'],
#                             'num_peaks_per_label': processing_parameters['CountingFishNumPeaksPerLabel'],
#                                 }


#     fill_value = np.nan
#     counts = osmFISH_dots_thr_selection(img,counting_parameters_dict)
#     counts.counting_graph()
#     counts.thr_identification()
#     data_models = Output_models()
#     counts_dict = data_models.dots_counts_dict

#     # Initialise an empty version of the counts dict
#     counts_dict['r_px_original'] = np.array([fill_value])
#     counts_dict['c_px_original'] = np.array([fill_value])
#     counts_dict['dot_id'] = np.array([fill_value])
#     counts_dict['dot_intensity'] = np.array([fill_value])
#     counts_dict['selected_thr'] = np.array([fill_value])
                    
#     if not np.isnan(counts.selected_thr):
#             dots = osmFISH_dots_mapping(img,counts.selected_thr,counting_parameters_dict)
#             if isinstance(dots.selected_peaks,np.ndarray):
#                 # Peaks have been identified
#                 total_dots = dots.selected_peaks.shape[0]
#                 dot_id_array = np.array([str(fov)+'_'+str(round_num)+'_'+ channel +'_'+str(nid) for nid in range(total_dots)])
#                 fov_array = np.repeat(fov,total_dots)
#                 thr_array = np.repeat(counts.selected_thr,total_dots)
#                 channel_array = np.repeat(channel,total_dots)
#                 hybridization_num_array = np.repeat(round_num,total_dots)
#                 target_name_array = np.repeat(target_name,total_dots)

#                 counts_dict['r_px_original']  = dots.selected_peaks[:,0]
#                 counts_dict['c_px_original'] = dots.selected_peaks[:,1]
#                 counts_dict['dot_id'] = dot_id_array
#                 counts_dict['dot_intensity'] = dots.intensity_array
#                 counts_dict['selected_thr'] = thr_array
#             else:
#                 logger.info(f' fov {fov} does not have counts (mapping)')
                
#     else:
#         logger.info(f' fov {fov} does not have counts (thr)')
    
#     counts_df = pd.DataFrame(counts_dict)
#     fov_subdataset_df = pd.DataFrame(fov_subdataset).T
#     fov_subdataset_df = pd.concat([fov_subdataset_df]*counts_df.shape[0],axis=0).sort_index().reset_index(drop=True)
#     counts_df = pd.concat([counts_df,fov_subdataset_df],axis=1)
#     return counts_df


# def osmFISH_barcoded_peak_based_detection_masked_thr_test(img_meta:Tuple[np.ndarray, Dict],
#                 masked_img:np.ndarray,
#                 processing_parameters:Dict):

#     # Need to add all the cases with no counts
#     """
#     This class apply the same strategy used for dots calling in osmFISH
#     on files that are structured for barcoded analysis using a masked image
#     for the selection of the thr

#     Attributes:
#     -----------
#     fov_name: str
#         name of the fov that is processed ex. fov_1
#     img_stack: np.ndarray
#         image stack containing in which each layer 
#         correspond to a round
#     parameters_dic:
#         dict with the parameters used for the counting

#     """

#     logger = selected_logger()

#     img = img_meta[0]
#     img_metadata = img_meta[1]
#     fov = img_metadata['fov_num']
#     hybridization = img_metadata['hybridization_num']
#     target_name = img_metadata['target_name']
    
#     logger.info(f'logging osmFISH_peak_based_detection fov {fov}')
#     hybridization_num = img_metadata['hybridization_num']

#     counting_parameters_dict = {
#                             'min_distance': processing_parameters['CountingFishMinObjDistance'],
#                             'min_obj_size': processing_parameters['CountingFishMinObjSize'],
#                             'max_obj_size': processing_parameters['CountingFishMaxObjSize'],
#                             'num_peaks_per_label': processing_parameters['CountingFishNumPeaksPerLabel'],
#                                 }

#     fill_value = np.nan
        

#     counts = osmFISH_dots_thr_selection(masked_img,counting_parameters_dict)
#     counts.counting_graph()
#     counts.thr_identification()

#     data_models = Output_models()
#     counts_dict = data_models.dots_counts_dict

#     # Initialise an empty version of the counts dict
#     counts_dict['r_px_original'] = np.array([fill_value])
#     counts_dict['c_px_original'] = np.array([fill_value])
#     counts_dict['dot_id'] = np.array([fill_value])
#     counts_dict['fov_num'] = np.array(fov)
#     counts_dict['round_num'] = np.array([img_metadata['hybridization_num']])
#     counts_dict['dot_intensity'] = np.array([fill_value])
#     counts_dict['selected_thr'] = np.array([fill_value])
#     counts_dict['dot_channel'] = np.array([img_metadata['channel']])
#     counts_dict['target_name'] = np.array([img_metadata['target_name']])
                    
#     if not np.isnan(counts.selected_thr):
#             dots = osmFISH_dots_mapping(img,counts.selected_thr,counting_parameters_dict)
#             if isinstance(dots.selected_peaks,np.ndarray):
#                 # Peaks have been identified
#                 total_dots = dots.selected_peaks.shape[0]
#                 dot_id_array = np.array([str(fov)+'_'+str(hybridization_num)+'_'+ img_metadata['channel'] +'_'+str(nid) for nid in range(total_dots)])
#                 fov_array = np.repeat(fov,total_dots)
#                 thr_array = np.repeat(counts.selected_thr,total_dots)
#                 channel_array = np.repeat(img_metadata['channel'],total_dots)
#                 hybridization_num_array = np.repeat(img_metadata['hybridization_num'],total_dots)
#                 target_name_array = np.repeat(img_metadata['target_name'],total_dots)

#                 counts_dict['r_px_original']  = dots.selected_peaks[:,0]
#                 counts_dict['c_px_original'] = dots.selected_peaks[:,1]
#                 counts_dict['dot_id'] = dot_id_array
#                 counts_dict['fov_num'] = fov_array
#                 counts_dict['round_num'] = hybridization_num_array
#                 counts_dict['dot_intensity'] = dots.intensity_array
#                 counts_dict['selected_thr'] = thr_array
#                 counts_dict['dot_channel'] = channel_array
#                 counts_dict['target_name'] = target_name_array
#             else:
#                 logger.info(f' fov {fov} does not have counts (mapping)')
                
#     else:
#         logger.info(f' fov {fov} does not have counts (thr)')
    
#     return (counts_dict, img_metadata)





# def osmFISH_barcoded_peak_based_detection_masked_thr(img_meta:Tuple[np.ndarray, Dict],
#                 masked_img:np.ndarray,
#                 min_distance: np.float64,
#                 min_obj_size: np.uint16,
#                 max_obj_size: np.uint16,
#                 num_peaks_per_label: np.uint16):

#     # Need to add all the cases with no counts
#     """
#     This class apply the same strategy used for dots calling in osmFISH
#     on files that are structured for barcoded analysis using a masked image
#     for the selection of the thr

#     Attributes:
#     -----------
#     fov_name: str
#         name of the fov that is processed ex. fov_1
#     img_stack: np.ndarray
#         image stack containing in which each layer 
#         correspond to a round
#     parameters_dic:
#         dict with the parameters used for the counting

#     """

#     logger = selected_logger()

#     img = img_meta[0]
#     img_metadata = img_meta[1]
#     fov = img_metadata['fov_num']
#     hybridization = img_metadata['hybridization_num']
#     target_name = img_metadata['target_name']
    
#     logger.info(f'logging osmFISH_peak_based_detection fov {fov}')
#     hybridization_num = img_metadata['hybridization_num']

#     counting_parameters_dict = {
#                             'min_distance': min_distance,
#                             'min_obj_size': min_obj_size,
#                             'max_obj_size': max_obj_size,
#                             'num_peaks_per_label': num_peaks_per_label,
#                                 }

#     fill_value = np.nan
        

#     counts = osmFISH_dots_thr_selection(masked_img,counting_parameters_dict)
#     counts.counting_graph()
#     counts.thr_identification()

#     data_models = Output_models()
#     counts_dict = data_models.dots_counts_dict

#     # Initialise an empty version of the counts dict
#     counts_dict['r_px_original'] = np.array([fill_value])
#     counts_dict['c_px_original'] = np.array([fill_value])
#     counts_dict['dot_id'] = np.array([fill_value])
#     counts_dict['fov_num'] = np.array(fov)
#     counts_dict['round_num'] = np.array([img_metadata['hybridization_num']])
#     counts_dict['dot_intensity'] = np.array([fill_value])
#     counts_dict['selected_thr'] = np.array([fill_value])
#     counts_dict['dot_channel'] = np.array([img_metadata['channel']])
#     counts_dict['target_name'] = np.array([img_metadata['target_name']])
                    
#     if not np.isnan(counts.selected_thr):
#             dots = osmFISH_dots_mapping(img,counts.selected_thr,counting_parameters_dict)
#             if isinstance(dots.selected_peaks,np.ndarray):
#                 # Peaks have been identified
#                 total_dots = dots.selected_peaks.shape[0]
#                 dot_id_array = np.array([str(fov)+'_'+str(hybridization_num)+'_'+ img_metadata['channel'] +'_'+str(nid) for nid in range(total_dots)])
#                 fov_array = np.repeat(fov,total_dots)
#                 thr_array = np.repeat(counts.selected_thr,total_dots)
#                 channel_array = np.repeat(img_metadata['channel'],total_dots)
#                 hybridization_num_array = np.repeat(img_metadata['hybridization_num'],total_dots)
#                 target_name_array = np.repeat(img_metadata['target_name'],total_dots)

#                 counts_dict['r_px_original']  = dots.selected_peaks[:,0]
#                 counts_dict['c_px_original'] = dots.selected_peaks[:,1]
#                 counts_dict['dot_id'] = dot_id_array
#                 counts_dict['fov_num'] = fov_array
#                 counts_dict['round_num'] = hybridization_num_array
#                 counts_dict['dot_intensity'] = dots.intensity_array
#                 counts_dict['selected_thr'] = thr_array
#                 counts_dict['dot_channel'] = channel_array
#                 counts_dict['target_name'] = target_name_array
#             else:
#                 logger.info(f' fov {fov} does not have counts (mapping)')
                
#     else:
#         logger.info(f' fov {fov} does not have counts (thr)')
    
#     return (counts_dict, img_metadata)


# def osmFISH_dots_thr_selection_np(img:np.ndarray, min_distance:int, 
#                                min_int:float=False, 
#                                max_int:float=False,
#                                min_peaks:int=False):

#     if min_peaks == False:
#         min_peaks = 3


#     fill_value = 9999

#     # List with the total peaks calculated for each threshold
#     total_peaks = []
#     thr_used = []

#     binning = 100
#     # Define the range of thr to be tested
#     if min_int and max_int:
#         thr_array = np.linspace(min_int,max_int,num=binning)
#     elif min_int:
#         thr_array = np.linspace(min_int,img.max(),num=binning)
#     elif max_int:
#         thr_array = np.linspace(np.min(img[np.nonzero(img)]),max_int,num=binning)
#     else:
#         thr_array = np.linspace(np.min(img[np.nonzero(img)]),img.max(),num=binning)

#     # Calculate the number of peaks for each threshold. In this calculation
#     # the size of the objects is not considered
#     peak_counter_min = 0
#     peak_counter_max = 0
#     for vl, thr in enumerate(thr_array):
#         # The border is excluded from the counting
#         peaks = feature.peak_local_max(img,min_distance=min_distance,\
#             threshold_abs=thr,exclude_border=False, indices=True,\
#             num_peaks=np.inf, footprint=None,labels=None)

#         number_peaks = len(peaks)

#         # Stop the counting when the number of peaks detected falls below 3
#         if number_peaks<=min_peaks:
#             stop_thr = thr # Move in the upper loop so you will stop at the previous thr
#             break
#         else:
#             total_peaks.append(len(peaks))
#             thr_used.append(thr)

#     # Consider the case of no detectected peaks or if there is only one Thr
#     # that create peaks (list total_peaks have only one element and )
#     # if np.array(total_peaks).sum()>0 or len(total_peaks)>1:
#     if len(total_peaks)>1:

#         # Trim the threshold array in order to match the stopping point
#         # used the [0][0] to get the first number and then take it out from list
#         # thr_array = thr_array[:np.where(thr_array==stop_thr)[0][0]]
#         thr_array = np.array(thr_used)

#         # Calculate the gradient of the number of peaks distribution
#         grad = np.gradient(total_peaks)

#         # Restructure the data in order to avoid to consider the min_peak in the
#         # calculations

#         # Coord of the gradient min_peak
#         grad_min_peak_coord = np.argmin(grad)

#         # Trim the data to remove the peak.
#         trimmed_thr_array = thr_array[grad_min_peak_coord:]
#         trimmed_grad = grad[grad_min_peak_coord:]

#         if trimmed_thr_array.shape>(1,):

#             # Trim the coords array in order to maintain the same length of the 
#             # tr and pk
#             trimmed_total_peaks = total_peaks[grad_min_peak_coord:]

#             # To determine the threshold we will determine the Thr with the biggest
#             # distance to the segment that join the end points of the calculated
#             # gradient

            
#             # Calculate the coords of the end points of the gradient
#             p1 = np.array([trimmed_thr_array[0],trimmed_grad[0]])
#             p2 = np.array([trimmed_thr_array[-1],trimmed_grad[-1]])

#             # Create a line that join the points
#             allpoints = np.arange(0,len(trimmed_thr_array))
#             allpoints_coords = np.array([trimmed_thr_array[allpoints],trimmed_grad[allpoints]]).T

#             distances = []
#             for point in allpoints_coords:
#                 distances.append(np.linalg.norm(np.cross(p2-p1, p1-point))/np.linalg.norm(p2-p1))

# #             # Calculate the distance between all points and the line
# #             for p in allpoints:
# #                 dst = s.distance(Point(trimmed_thr_array[p],trimmed_grad[p]))
# #                 distances.append(dst.evalf())

#             # Remove the end points from the lists
#             trimmed_thr_array = trimmed_thr_array[1:-1]
#             trimmed_grad = trimmed_grad[1:-1]
#             trimmed_total_peaks = trimmed_total_peaks[1:-1]
#             trimmed_distances = distances[1:-1]

#             # Determine the coords of the selected Thr
#             # Converted trimmed_distances to array because it crashed
#             # on Sanger.
#             if trimmed_distances: # Most efficient way will be to consider the length of Thr list
#                 thr_idx = np.argmax(np.array(trimmed_distances))
#                 selected_thr = trimmed_thr_array[thr_idx]
#                 # The selected threshold usually causes oversampling of the number of dots
#                 # I added a stringency parameter (int n) to use to select the Thr+n 
#                 # for the counting. It selects a stringency only if the trimmed_thr_array
#                 # is long enough. Also consider the case in which the stringency in negative
#             else:
#                 selected_thr = fill_value
#                 trimmed_thr_array = fill_value
#         else:
#             selected_thr = fill_value
#             trimmed_thr_array = fill_value
#     else:
#         selected_thr = fill_value
#         trimmed_thr_array = fill_value
    
#     return selected_thr
