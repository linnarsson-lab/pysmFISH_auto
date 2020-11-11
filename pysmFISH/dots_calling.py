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
from skimage import feature, measure, img_as_float
from scipy import ndimage as nd
from pathlib import Path

from pysmFISH.utils import load_pipeline_config_file, convert_from_uint16_to_float64


from prefect import task
from prefect.engine import signals
from prefect.utilities.logging import get_logger


from pysmFISH.logger_utils import prefect_logging_setup, simple_writing_logger



class osmFISH_dots_thr_selection():

    '''
    Class used to calculate and optimise the thr for the identification of the dots.
    

    Attributes:
    -----------
    img

    Returns:
    ---------
    '''
    
    def __init__(self, img:np.ndarray, parameters_dict:Dict, min_int:float=False, max_int:float=False,min_peaks:int=False):
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

        # Define the range of thr to be tested
        if self.min_int and self.max_int:
            self.thr_array = np.linspace(self.min_int,self.max_int,num=100)
        elif self.min_int:
            self.thr_array = np.linspace(self.min_int,self.img.max(),num=100)
        elif self.max_int:
            self.thr_array = np.linspace(self.img.min(),self.max_int,num=100)
        else:
            self.thr_array = np.linspace(self.img.min(),self.img.max(),num=100)
    
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

    def __init__(self,img,thr,parameters_dict):
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
                                threshold_abs=self.thr, exclude_border=False, indices=True, num_peaks=np.inf, 
                                footprint=None, labels=labels,num_peaks_per_label=self.num_peaks_per_label)                            
        
        if self.selected_peaks.size:
            self.intensity_array = self.img[self.selected_peaks[:,0],self.selected_peaks[:,1]]
        else:
            self.intensity_array = np.nan


# @task(name='peak-based-detection')
def osmFISH_peak_based_detection(img_meta:tuple,
                                        min_distance: np.float64,
                                        min_obj_size: np.uint16,
                                        max_obj_size: np.uint16,
                                        num_peaks_per_label: np.uint16):
    
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

    
    
    logger = get_logger('logging osmFISH_peak_based_detection')

    img = img_meta[0]
    img_metadata = img_meta[1]
    fov = img_metadata['fov_num']
    
    logger.info(f'logging osmFISH_peak_based_detection fov {fov}')
    hybridization_num = img_metadata['hybridization_num']

    counting_parameters_dict = {
                            'min_distance': min_distance,
                            'min_obj_size': min_obj_size,
                            'max_obj_size': max_obj_size,
                            'num_peaks_per_label': num_peaks_per_label,
                                }
    fill_value = np.nan
    counts = osmFISH_dots_thr_selection(img,counting_parameters_dict)
    counts.counting_graph()
    counts.thr_identification()

    if not np.isnan(counts.selected_thr):
            dots = osmFISH_dots_mapping(img,counts.selected_thr,counting_parameters_dict)
            if isinstance(dots.selected_peaks,np.ndarray):
                total_dots = dots.selected_peaks.shape[0]
                dot_id_array = np.array([str(fov)+'_'+str(hybridization_num)+'_'+ img_metadata['channel'] +'_'+str(nid) for nid in range(total_dots)])
                fov_array = np.repeat(fov,total_dots)
                thr_array = np.repeat(counts.selected_thr,total_dots)
                channel_array = np.repeat(img_metadata['channel'],total_dots)
                counts_dict = {
                    'DotsCoordsFOV': dots.selected_peaks,
                    'DotID': dot_id_array,
                    'FovNumber': fov_array,
                    'DotIntensity': dots.intensity_array,
                    'SelectedThreshold':thr_array,
                    'DotChannel':channel_array}
            else:
                logger.info(f' fov {fov} does not have counts (mapping)')
                counts_dict = {
                    'DotsCoordsFOV': np.array([fill_value, fill_value]),
                    'DotID': np.array([fill_value]),
                    'FovNumber': np.array(fov),
                    'DotIntensity': np.array([fill_value]),
                    'SelectedThreshold':np.array([fill_value]),
                    'DotChannel':np.array([fill_value])}
    else:
        logger.info(f' fov {fov} does not have counts (thr)')
        counts_dict = {
                    'DotsCoordsFOV': np.array([fill_value, fill_value]),
                    'DotID': np.array([fill_value]),
                    'FovNumber': np.array(fov),
                    'DotIntensity': np.array([fill_value]),
                    'SelectedThreshold':np.array([fill_value]),
                    'DotChannel':np.array([fill_value])}
    
    return (counts_dict, img_metadata)