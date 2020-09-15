from typing import *
import numpy as np
import argparse
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


class osmFISH_barcoded_peak_based_detection():

    # Need to add all the cases with no counts
    """
    This class apply the same strategy used for dots calling in osmFISH
    on files that are structured for barcoded analysis

    Attributes:
    -----------
    fov_name: str
        name of the fov that is processed ex. fov_1
    img_stack: np.ndarray
        image stack containing in which each layer 
        correspond to a round
    parameters_dic:
        dict with the parameters used for the counting

    """

    def __init__(self, fov_name:str, img_stack:np.ndarray, counting_parameters_dict:Dict):
        self.fov_name = fov_name
        self.img_stack = img_stack
        self.counting_parameters_dict = counting_parameters_dict
        self.fov = int(self.fov_name.split('_')[-1])

        self.fill_value = np.nan
        self.logger = logging.getLogger(__name__)
    
    def call_dots(self):
        try:
            number_layers = self.img_stack.shape[0]
        except:
            self.logger.error(f'image is not a stack')

        column_names = ['r_px_original','c_px_original','dot_id','fov_num','round_num','dot_intensity_norm','dot_intensity_not','selected_thr' ]
        self.output = pd.DataFrame(columns = column_names)

        if number_layers > 1:
            self.max_thr = []
            # First round for calcolating thr range
            for layer_idx in range(number_layers):
                counts = osmFISH_dots_thr_selection(self.img_stack[layer_idx,:,:],self.counting_parameters_dict)
                counts.counting_graph()
                counts.thr_identification()
                if not np.isnan(counts.selected_thr):
                        dots = osmFISH_dots_mapping(self.img_stack[layer_idx,:,:],counts.selected_thr,self.counting_parameters_dict)
                        if isinstance(dots.selected_peaks,np.ndarray):
                            total_dots = dots.selected_peaks.shape[0]
                            dot_id_array = np.array([str(self.fov)+'_'+str(layer_idx+1)+'_'+str(nid) for nid in range(total_dots)])
                            fov_array = np.repeat(self.fov,total_dots)
                            round_array = np.repeat(layer_idx+1,total_dots)
                            thr_array = np.repeat(counts.selected_thr,total_dots)
                            intensity_array_non_normalized = self.img_stack[layer_idx,:,:][dots.selected_peaks[:,0],dots.selected_peaks[:,1]]
                            self.output_round_dict = {'r_px_original': dots.selected_peaks[:,0],
                                'c_px_original': dots.selected_peaks[:,1],
                                'dot_id': dot_id_array,
                                'fov_num': fov_array,
                                'round_num': round_array,
                                'dot_intensity_norm': dots.intensity_array,
                                'dot_intensity_not': intensity_array_non_normalized,
                                'selected_thr':thr_array}
                        else:
                            self.output_round_dict = {'r_px_original': np.array([self.fill_value]),
                                'c_px_original': np.array([self.fill_value]),
                                'dot_id': np.array([self.fill_value]),
                                'fov_num': np.array(self.fov),
                                'round_num': np.array([layer_idx+1]),
                                'dot_intensity_norm': np.array([self.fill_value]),
                                'dot_intensity_not': np.array([self.fill_value]),
                                'selected_thr':np.array([self.fill_value])}
                else:
                    self.output_round_dict = {'r_px_original': np.array([self.fill_value]),
                            'c_px_original': np.array([self.fill_value]),
                            'dot_id': np.array([self.fill_value]),
                            'fov_num': np.array(self.fov),
                            'round_num': np.array([layer_idx+1]),
                            'dot_intensity_norm': np.array([self.fill_value]),
                            'dot_intensity_not': np.array([self.fill_value]),
                            'selected_thr':np.array([self.fill_value])}
                self.output = pd.concat([self.output,pd.DataFrame(self.output_round_dict)],axis=0,ignore_index=True)
        else:

            self.logger.error(f'this image in not suitable for barcoded experiments')



class osmFISH_serial_peak_detection():

    # Need to add all the cases with no counts
    """
    This class calculate the number of dots after normalization of the
    image

    Attributes:
    -----------
    fov_name: str
        name of the fov that is processed ex. fov_1
    img_stack: np.ndarray
        image stack containing in which each layer 
        correspond to a round
    parameters_dic:
        dict with the parameters used for the counting
   
         
    """

    def __init__(self, fov_name:str, img_stack:np.ndarray, parameters_dict:Dict):
        self.fov_name = fov_name
        self.img_stack = img_stack
        self.parameters_dict = parameters_dict
        self.fov = int(self.fov_name.split('_')[-1])

        self.fill_value = np.nan
        self.logger = logging.getLogger(__name__)
    
    def call_dots(self):
        
        try:
            number_layers = self.img_stack.shape[0]
        except:
            self.logger.error(f'image is not a stack')


        column_names = ['r_px_original','c_px_original','dot_id','fov_num','round_num','dot_intensity','selected_thr' ]
        self.output = pd.DataFrame(columns = column_names)

        counts = osmFISH_dots_thr_selection( self.img_stack,self.parameters_dict)
        counts.counting_graph()
        counts.thr_identification()
        if not np.isnan(counts.selected_thr):
            dots = osmFISH_dots_mapping(self.img_stack,counts.selected_thr,self.parameters_dict)
            if isinstance(dots.selected_peaks,np.ndarray):
                total_dots = dots.selected_peaks.shape[0]
                dot_id_array = np.array([str(self.fov)+'_'+str(nid) for nid in range(total_dots)])
                fov_array = np.repeat(self.fov,total_dots)
                thr_array = np.repeat(counts.selected_thr,total_dots)
                intensity_array= self.img_stack[dots.selected_peaks[:,0],dots.selected_peaks[:,1]]
                self.output_round_dict = {'r_px_original': dots.selected_peaks[:,0],
                    'c_px_original': dots.selected_peaks[:,1],
                    'dot_id': dot_id_array,
                    'fov_num': fov_array,
                    'dot_intensity': dots.intensity_array,
                    'selected_thr':thr_array}
            else:
                self.output_round_dict = {'r_px_original': np.array([self.fill_value]),
                    'c_px_original': np.array([self.fill_value]),
                    'dot_id': np.array([self.fill_value]),
                    'fov_num': np.array(self.fov),
                    'dot_intensity': np.array([self.fill_value]),
                    'selected_thr':np.array([self.fill_value])}
        else:
            self.output_round_dict = {'r_px_original': np.array([self.fill_value]),
                    'c_px_original': np.array([self.fill_value]),
                    'dot_id': np.array([self.fill_value]),
                    'fov_num': np.array(self.fov),
                    'dot_intensity': np.array([self.fill_value]),
                    'selected_thr':np.array([self.fill_value])}
        self.output = pd.concat([self.output,pd.DataFrame(self.output_round_dict)],axis=0,ignore_index=True)



class osmFISH_serial_peak_detection_defined_thr():

    # Need to add all the cases with no counts
    """
    This class calculate the number of dots after normalization of the
    image

    Attributes:
    -----------
    fov_name: str
        name of the fov that is processed ex. fov_1
    img_stack: np.ndarray
        image stack containing in which each layer 
        correspond to a round
    parameters_dic:
        dict with the parameters used for the counting
   
         
    """

    def __init__(self, fov_name:str, img_stack:np.ndarray, parameters_dict:Dict, thr:float):
        self.fov_name = fov_name
        self.img_stack = img_stack
        self.parameters_dict = parameters_dict
        self.thr = thr
        self.fov = int(self.fov_name.split('_')[-1])

        self.fill_value = np.nan
        self.logger = logging.getLogger(__name__)
    
    def call_dots(self):
        
        try:
            number_layers = self.img_stack.shape[0]
        except:
            self.logger.error(f'image is not a stack')


        column_names = ['r_px_original','c_px_original','dot_id','fov_num','round_num','dot_intensity','selected_thr' ]
        self.output = pd.DataFrame(columns = column_names)

        dots = osmFISH_dots_mapping(self.img_stack,self.thr,self.parameters_dict)
        if isinstance(dots.selected_peaks,np.ndarray):
            total_dots = dots.selected_peaks.shape[0]
            dot_id_array = np.array([str(self.fov)+'_'+str(nid) for nid in range(total_dots)])
            fov_array = np.repeat(self.fov,total_dots)
            thr_array = np.repeat(self.thr,total_dots)
            intensity_array= self.img_stack[dots.selected_peaks[:,0],dots.selected_peaks[:,1]]
            self.output_round_dict = {'r_px_original': dots.selected_peaks[:,0],
                'c_px_original': dots.selected_peaks[:,1],
                'dot_id': dot_id_array,
                'fov_num': fov_array,
                'dot_intensity': dots.intensity_array,
                'selected_thr':thr_array}
        else:
            self.output_round_dict = {'r_px_original': np.array([self.fill_value]),
                'c_px_original': np.array([self.fill_value]),
                'dot_id': np.array([self.fill_value]),
                'fov_num': np.array(self.fov),
                'dot_intensity': np.array([self.fill_value]),
                'selected_thr':np.array([self.fill_value])}
        self.output = pd.concat([self.output,pd.DataFrame(self.output_round_dict)],axis=0,ignore_index=True)