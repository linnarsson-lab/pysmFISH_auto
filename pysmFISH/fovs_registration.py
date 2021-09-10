"""
group of class or functions use to register the fov between
different rounds.
"""
from typing import *
import logging
import zarr
import pickle
import pandas as pd
import numpy as np
import sys
from pathlib import Path
from sklearn.neighbors import NearestNeighbors
from skimage import transform
from skimage.measure import ransac
from scipy.spatial import distance
from scipy.ndimage import fourier_shift

from skimage.feature import register_translation
# from skimage.registration import phase_cross_correlation UPDATE SKIMAGEN
from skimage import filters
from skimage.registration import phase_cross_correlation


from sklearn.neighbors import NearestNeighbors

import itertools
import math
import operator
from scipy.optimize import minimize

import gc


from pysmFISH.logger_utils import selected_logger
from pysmFISH.errors import Registration_errors

from pysmFISH.data_models import Output_models
from pysmFISH import utils

def create_fake_image(img_shape: np.ndarray,coords: np.ndarray)->np.ndarray:
    """Function used to create an image from dots counts. The image
    will be used for FFT based registration. The dots are mapped in
    the image and the signal is enhanced using a gaussian. This
    increase the features that can be used in the FFT based
    registration.

    Args:
        img_shape (np.ndarray): Shape of the image to create. It matches
            the shape of the images to process.
        coords (np.ndarray): Coordinates of the peaks detected.

    Returns:
        np.ndarray: synthetic image used for registration
    """
    gaussian_sigma = 5
    img = np.zeros(img_shape,dtype=np.float64)
    img[coords[:,0].astype(int),coords[:,1].astype(int)] = 1000
    img = filters.gaussian(img,sigma=gaussian_sigma)
    return img


def register_images(img: np.ndarray, shift: np.ndarray)->np.ndarray:
    """Function to create a new image shifted according
    a predefined shift.

    Args:
        img (np.ndarray): Image to shift
        shift (np.ndarray): Shift

    Returns:
        np.ndarray: Shifted image
    """
    offset_image = fourier_shift(np.fft.fftn(img), shift)
    offset_image = np.fft.ifftn(offset_image).real
    return offset_image


def combine_register_filtered_images(output_dict: dict, metadata: dict, 
                                    all_rounds_shifts:dict)->np.ndarray:
    """Function used to register the image throughout all rounds.

    Args:
        output_dict (dict): dict containg output of the filtering ((img,),metadata) organised
                            by channe and round.
        metadata (dict): Metadata that characterize the acquired images
        registered_counts (pd.DataFrame): Counts after registration.

    Returns:
        np.ndarray: Image stack with all the rounds registered.
    """
    registered_img_stack = {}
    if output_dict:
        for channel, all_rounds_data in output_dict.items():
            img_stack = np.zeros([int(metadata['total_rounds']),int(metadata['img_width']),int(metadata['img_height'])])
            for round_num, filt_out in all_rounds_data.items():
                img = filt_out[0][0]
                shift = all_rounds_shifts[int(round_num)]
                if shift.size:
                    shifted_img = register_images(img,shift)
                    img_stack[round_num-1,:,:] = shifted_img
            registered_img_stack[channel] = img_stack
    return registered_img_stack


def beads_based_registration(all_counts_fov: pd.DataFrame,
                            analysis_parameters: Dict)->pd.DataFrame:
    """[summary]

    Args:
        all_counts_fov (pd.DataFrame): [description]
        analysis_parameters (Dict): [description]

    Returns:
        pd.DataFrame: [description]
    """
    
    # Used index to avoid to remake the output dataframe
    
    stitching_channel = all_counts_fov['stitching_channel'].iloc[0]
    
    stitching_channel_df = all_counts_fov.loc[(all_counts_fov.channel == stitching_channel) &
                                              (all_counts_fov.mapped_beads_type == 'large')  , :]
    
    if stitching_channel_df.shape[0] < 2:
        stitching_channel_df = all_counts_fov.loc[(all_counts_fov.channel == stitching_channel) , :]

    
    fish_df = all_counts_fov.loc[all_counts_fov.channel != stitching_channel,:]
    
    reference_round_num = analysis_parameters['RegistrationReferenceHybridization']
    registration_tollerance_pxl = analysis_parameters['RegistrationTollerancePxl']
    
    registration_errors = Registration_errors()
    
    
    # Determine if there are any round with missing counts in the registration
    if stitching_channel_df[stitching_channel_df['dot_id'].isnull()].empty :
    
        all_counts_fov['r_px_registered'] = np.nan
        all_counts_fov['c_px_registered'] = np.nan
        all_counts_fov['r_shift_px'] = np.nan
        all_counts_fov['c_shift_px'] = np.nan
        all_counts_fov['min_number_matching_dots_registration'] = np.nan

        # Register stitching channel
        ref_counts_df = stitching_channel_df.loc[stitching_channel_df.round_num == reference_round_num,:]

        all_counts_fov.loc[ref_counts_df.index,'r_px_registered'] =  \
                ref_counts_df.loc[ref_counts_df.round_num == reference_round_num,'r_px_original']
        all_counts_fov.loc[ref_counts_df.index,'c_px_registered'] =  \
                ref_counts_df.loc[ref_counts_df.round_num == reference_round_num,'c_px_original']
        all_counts_fov.loc[ref_counts_df.index,'r_shift_px'] =  0
        all_counts_fov.loc[ref_counts_df.index,'c_shift_px'] =  0
        all_counts_fov.loc[ref_counts_df.index,'min_number_matching_dots_registration'] =  1000

        # Register fish
        fish_ref_round = fish_df.loc[fish_df.round_num == reference_round_num,:]
        
        all_counts_fov.loc[fish_ref_round.index,'r_px_registered'] =  \
                fish_ref_round.loc[fish_ref_round.index,'r_px_original']
        all_counts_fov.loc[fish_ref_round.index,'c_px_registered'] =  \
                fish_ref_round.loc[fish_ref_round.index,'c_px_original']
        all_counts_fov.loc[fish_ref_round.index,'r_shift_px'] =  0
        all_counts_fov.loc[fish_ref_round.index,'c_shift_px'] =  0
        all_counts_fov.loc[fish_ref_round.index,'min_number_matching_dots_registration'] =  1000
        
        # Create reference fake image for registration
        img_width = ref_counts_df.iloc[0]['img_width']
        img_height = ref_counts_df.iloc[0]['img_height']
        ref_coords = ref_counts_df.loc[:,['r_px_original', 'c_px_original']].to_numpy()
        img_ref = create_fake_image((img_width, img_height),ref_coords)

        all_rounds = all_counts_fov.round_num.unique()
        all_rounds = all_rounds[all_rounds != reference_round_num]
        for tran_round_num in all_rounds:

           
            tran_counts_df = stitching_channel_df.loc[stitching_channel_df.round_num == tran_round_num,:]
            img_width = tran_counts_df.iloc[0]['img_width']
            img_height = tran_counts_df.iloc[0]['img_height']
            tran_coords = tran_counts_df.loc[:,['r_px_original', 'c_px_original']].to_numpy()
            img_tran = create_fake_image((img_width, img_height),tran_coords)


            shift, error, diffphase = register_translation(img_ref, img_tran)
            registered_tran_coords = tran_coords + shift
            min_num_matching_dots = identify_matching_register_dots_NN(ref_coords,
                                                        registered_tran_coords,
                                                        registration_tollerance_pxl)


            # Register stitching channel
            all_counts_fov.loc[tran_counts_df.index,'r_px_registered'] =  registered_tran_coords[:,0]
            all_counts_fov.loc[tran_counts_df.index,'c_px_registered'] =  registered_tran_coords[:,1]
            all_counts_fov.loc[tran_counts_df.index,'r_shift_px'] =  shift[0]
            all_counts_fov.loc[tran_counts_df.index,'c_shift_px'] =  shift[1]
            all_counts_fov.loc[tran_counts_df.index,
                        'min_number_matching_dots_registration'] =  min_num_matching_dots
            
            
            # Register fish
            fish_ref_round = fish_df.loc[fish_df.round_num == tran_round_num,:]
            all_counts_fov.loc[fish_ref_round.index,'r_px_registered'] =  fish_ref_round['r_px_original'] + shift[0] 
            all_counts_fov.loc[fish_ref_round.index,'c_px_registered'] =  fish_ref_round['c_px_original'] + shift[1]
            all_counts_fov.loc[fish_ref_round.index,'r_shift_px'] =  shift[0] 
            all_counts_fov.loc[fish_ref_round.index,'c_shift_px'] =  shift[1]
            all_counts_fov.loc[fish_ref_round.index,
                        'min_number_matching_dots_registration'] =  min_num_matching_dots
            
            
    else:
        
        all_counts_fov['r_px_registered'] = np.nan
        all_counts_fov['c_px_registered'] = np.nan
        all_counts_fov['r_shift_px'] = np.nan
        all_counts_fov['c_shift_px'] = np.nan
        all_counts_fov['min_number_matching_dots_registration'] = registration_errors.missing_counts_reg_channel
     
    return all_counts_fov


def beads_based_registration_stitching_channel(stitching_channel_df: pd.DataFrame,
                analysis_parameters: Dict)-> Tuple[pd.DataFrame, np.ndarray, pd.DataFrame]:
    """Registration of the reference channel that contained reference beads.

    Args:
        stitching_channel_df (pd.DataFrame): Coordinates of the beads
        analysis_parameters (Dict): Processing parameters

    Returns:
        Tuple[pd.DataFrame, np.ndarray, pd.DataFrame]: (stitching_channel_df, 
                        all_rounds_shifts, all_rounds_matching_dots)
    
    """
                
    # Used index to avoid to remake the output dataframe
    

    reference_round_num = analysis_parameters['RegistrationReferenceHybridization']
    registration_tollerance_pxl = analysis_parameters['RegistrationTollerancePxl']
    
    registration_errors = Registration_errors()
    
    all_rounds_shifts = {}
    all_rounds_matching_dots = {}
    # Determine if there are any round with missing counts in the registration
    if stitching_channel_df[stitching_channel_df['dot_id'].isnull()].empty :
    
        stitching_channel_df['r_px_registered'] = np.nan
        stitching_channel_df['c_px_registered'] = np.nan
        stitching_channel_df['r_shift_px'] = np.nan
        stitching_channel_df['c_shift_px'] = np.nan
        stitching_channel_df['min_number_matching_dots_registration'] = np.nan

        # Register stitching channel
        ref_counts_df = stitching_channel_df.loc[stitching_channel_df.round_num == reference_round_num,:]

        stitching_channel_df.loc[ref_counts_df.index,'r_px_registered'] =  \
                ref_counts_df.loc[ref_counts_df.round_num == reference_round_num,'r_px_original']
        stitching_channel_df.loc[ref_counts_df.index,'c_px_registered'] =  \
                ref_counts_df.loc[ref_counts_df.round_num == reference_round_num,'c_px_original']
        stitching_channel_df.loc[ref_counts_df.index,'r_shift_px'] =  0
        stitching_channel_df.loc[ref_counts_df.index,'c_shift_px'] =  0
        stitching_channel_df.loc[ref_counts_df.index,'min_number_matching_dots_registration'] =  1000

        all_rounds_shifts[reference_round_num] = np.array([0,0])
        all_rounds_matching_dots[reference_round_num] = 1000

        # Create reference fake image for registration
        img_width = ref_counts_df.iloc[0]['img_width']
        img_height = ref_counts_df.iloc[0]['img_height']
        ref_coords = ref_counts_df.loc[:,['r_px_original', 'c_px_original']].to_numpy()
        img_ref = create_fake_image((img_width, img_height),ref_coords)

        all_rounds = stitching_channel_df.round_num.unique()
        all_rounds = all_rounds[all_rounds != reference_round_num]
        for tran_round_num in all_rounds:

           
            tran_counts_df = stitching_channel_df.loc[stitching_channel_df.round_num == tran_round_num,:]
            img_width = tran_counts_df.iloc[0]['img_width']
            img_height = tran_counts_df.iloc[0]['img_height']
            tran_coords = tran_counts_df.loc[:,['r_px_original', 'c_px_original']].to_numpy()
            img_tran = create_fake_image((img_width, img_height),tran_coords)


            shift, error, diffphase = register_translation(img_ref, img_tran)
            registered_tran_coords = tran_coords + shift
            min_num_matching_dots = identify_matching_register_dots_NN(ref_coords,
                                                        registered_tran_coords,
                                                        registration_tollerance_pxl)

            all_rounds_shifts[tran_round_num] = shift
            all_rounds_matching_dots[tran_round_num] = min_num_matching_dots
            

            # Register stitching channel
            stitching_channel_df.loc[tran_counts_df.index,'r_px_registered'] =  registered_tran_coords[:,0]
            stitching_channel_df.loc[tran_counts_df.index,'c_px_registered'] =  registered_tran_coords[:,1]
            stitching_channel_df.loc[tran_counts_df.index,'r_shift_px'] =  shift[0]
            stitching_channel_df.loc[tran_counts_df.index,'c_shift_px'] =  shift[1]
            stitching_channel_df.loc[tran_counts_df.index,
                        'min_number_matching_dots_registration'] =  min_num_matching_dots
            
    else:
        
        stitching_channel_df['r_px_registered'] = np.nan
        stitching_channel_df['c_px_registered'] = np.nan
        stitching_channel_df['r_shift_px'] = np.nan
        stitching_channel_df['c_shift_px'] = np.nan
        stitching_channel_df['min_number_matching_dots_registration'] = registration_errors.missing_counts_reg_channel
        all_rounds_shifts = np.nan
     
    return stitching_channel_df, all_rounds_shifts, all_rounds_matching_dots




def beads_based_registration_fish(all_counts_fov: pd.DataFrame,
                                        all_rounds_shifts: np.ndarray,
                                        all_rounds_matching_dots: pd.DataFrame,
                                        analysis_parameters: Dict)->pd.DataFrame:
    """Adjust the coords of the fish dots using the shift calculated in the 
    reference channel.

    Args:
        all_counts_fov (pd.DataFrame): Coords of all peaks
        all_rounds_shifts (np.ndarray): Calculated shift from the registration 
                                of the reference channel
        all_rounds_matching_dots (pd.DataFrame): Marching dots between rounds
        analysis_parameters (Dict): Processing parameters

    Returns:
        pd.DataFrame: All registered counts
    """
    
    registration_errors = Registration_errors()
    
    # Determine if there are any round with missing counts in the registration
    if all_counts_fov[all_counts_fov['dot_id'].isnull()].empty :
    
        all_counts_fov['r_px_registered'] = np.nan
        all_counts_fov['c_px_registered'] = np.nan
        all_counts_fov['r_shift_px'] = np.nan
        all_counts_fov['c_shift_px'] = np.nan
        all_counts_fov['min_number_matching_dots_registration'] = np.nan

        if isinstance(all_rounds_shifts,dict):

            for round_num, shift in all_rounds_shifts.items():
                    
                all_counts_fov.loc[all_counts_fov.round_num == round_num,'r_px_registered'] =  all_counts_fov['r_px_original'] + shift[0] 
                all_counts_fov.loc[all_counts_fov.round_num == round_num,'c_px_registered'] =  all_counts_fov['c_px_original'] + shift[1]
                all_counts_fov.loc[all_counts_fov.round_num == round_num,'r_shift_px'] =  shift[0] 
                all_counts_fov.loc[all_counts_fov.round_num == round_num,'c_shift_px'] =  shift[1]
                all_counts_fov.loc[all_counts_fov.round_num == round_num,
                            'min_number_matching_dots_registration'] =  all_rounds_matching_dots[round_num]
        else:

            all_counts_fov['r_px_registered'] = np.nan
            all_counts_fov['c_px_registered'] = np.nan
            all_counts_fov['r_shift_px'] = np.nan
            all_counts_fov['c_shift_px'] = np.nan
            all_counts_fov['min_number_matching_dots_registration'] = registration_errors.missing_counts_reg_channel
                 
    else:
        
        all_counts_fov['r_px_registered'] = np.nan
        all_counts_fov['c_px_registered'] = np.nan
        all_counts_fov['r_shift_px'] = np.nan
        all_counts_fov['c_shift_px'] = np.nan
        all_counts_fov['min_number_matching_dots_registration'] = registration_errors.missing_counts_reg_channel
     
    return all_counts_fov


# TODO MUST ADD THE CONSIDERATION OF POTENTIAL ERRORS
def nuclei_based_registration(all_counts_fov: pd.DataFrame,
                            img_stack: np.ndarray,
                            analysis_parameters: Dict)->pd.DataFrame:
    """Function used to run the registration of a single fov
    through all hybridization. The registration is done using the dots
    coords determined in each image

    Args:
        all_counts_fov (pd.DataFrame): Coords of all peaks
        img_stack (np.ndarray): Image stack of the nuclei from different rounds
        analysis_parameters (Dict): Processing parameters

    Returns:
        pd.DataFrame: All registered counts
    """
   
    
    logger = selected_logger()
    
    # Used index to avoid to remake the output dataframe
    
    stitching_channel = all_counts_fov['stitching_channel'].iloc[0]
    stitching_channel_df = all_counts_fov.loc[all_counts_fov.channel == stitching_channel, :]
    fish_df = all_counts_fov.loc[all_counts_fov.channel != stitching_channel,:]
    
    reference_round_num = analysis_parameters['RegistrationReferenceHybridization']
    ref_round_num_img = reference_round_num -1
    registration_tollerance_pxl = analysis_parameters['RegistrationTollerancePxl']
    
    registration_errors = Registration_errors()
    
    
    # Determine if there are any round with missing counts in the registration
    if stitching_channel_df[stitching_channel_df['dot_id'].isnull()].empty :
    
        all_counts_fov['r_px_registered'] = np.nan
        all_counts_fov['c_px_registered'] = np.nan
        all_counts_fov['r_shift_px'] = np.nan
        all_counts_fov['c_shift_px'] = np.nan
        all_counts_fov['min_number_matching_dots_registration'] = np.nan
        all_counts_fov['RMS'] = np.nan


        # Register fish and enter ref round info
        fish_ref_round = fish_df.loc[fish_df.round_num == reference_round_num,:]

        all_counts_fov.loc[fish_ref_round.index,'r_px_registered'] =  \
                fish_ref_round.loc[fish_ref_round.index,'r_px_original']
        all_counts_fov.loc[fish_ref_round.index,'c_px_registered'] =  \
                fish_ref_round.loc[fish_ref_round.index,'c_px_original']
        all_counts_fov.loc[fish_ref_round.index,'r_shift_px'] =  0
        all_counts_fov.loc[fish_ref_round.index,'c_shift_px'] =  0
        all_counts_fov.loc[fish_ref_round.index,'min_number_matching_dots_registration'] =  1000
        all_counts_fov.loc[fish_ref_round.index,'RMS'] =  0

        all_rounds = np.arange(img_stack.shape[0])
        all_rounds = all_rounds[all_rounds != ref_round_num_img]

        ref_img = img_stack[ref_round_num_img,:,:]
    
        for r_num in all_rounds:
            # Register fish
            fish_ref_round = fish_df.loc[fish_df.round_num == (r_num+1),:]

            shift, error, diffphase = phase_cross_correlation(ref_img, img_stack[r_num,:,:],return_error=True)    
            
            all_counts_fov.loc[fish_ref_round.index,'r_px_registered'] =  fish_ref_round['r_px_original'] + shift[0] 
            all_counts_fov.loc[fish_ref_round.index,'c_px_registered'] =  fish_ref_round['c_px_original'] + shift[1]
            all_counts_fov.loc[fish_ref_round.index,'r_shift_px'] =  shift[0] 
            all_counts_fov.loc[fish_ref_round.index,'c_shift_px'] =  shift[1]
            all_counts_fov.loc[fish_ref_round.index,
                        'min_number_matching_dots_registration'] =  error
            all_counts_fov.loc[fish_ref_round.index,'RMS'] =  0

    else:
        
        all_counts_fov['r_px_registered'] = np.nan
        all_counts_fov['c_px_registered'] = np.nan
        all_counts_fov['r_shift_px'] = np.nan
        all_counts_fov['c_shift_px'] = np.nan
        all_counts_fov['min_number_matching_dots_registration'] = registration_errors.missing_counts_reg_channel
        all_counts_fov['RMS'] = np.nan
     
    return all_counts_fov


# TODO Adjust and port this function in case of triangulation
# class triangles_based_registration():
#     """
#     Class used to register the different rounds by searaching and
#     matching all possible triangles formed by the dots in the reference
#     and translated image. This function run only a registration to the reference
#     round
    
#     The calculation of the triangle is based on list processing and may 
#     be improved in ported to numpy.
#     https://stackoverflow.com/questions/43126580/match-set-of-x-y-points-to-another-set-that-is-scaled-rotated-translated-and

#     """

#     def __init__(self, counts, experiment_fpath, channel_name, roi_number, fov_name):
#         self.counts = counts
#         self.experiment_fpath = Path(experiment_fpath)
#         self.channel_name = channel_name
#         self.roi_number = roi_number
#         self.fov_name = fov_name
        
#         self.logger = logging.getLogger(__name__)
        
#         self.pipeline_config_fpath = self.experiment_fpath / 'pipeline_config'
#         self.experiment_config_fpath =  self.pipeline_config_fpath / 'experiment.yaml'
#         self.experiment_config = load_pipeline_config_file(self.experiment_config_fpath)
        
#         searching_key = '*roi_' +str(self.roi_number) + '_' + self.channel_name + '_images_config.yaml'        
        
#         try:
#             self.image_config_fpath = list(self.pipeline_config_fpath.glob(searching_key))[0]
#         except:
#             self.logger.error(f'the reference beads image_config file is missing {searching_key}')
#             sys.exit(f'the reference beads image_config file')

        
#         # Load registration parameters
#         self.image_config = load_pipeline_config_file(self.image_config_fpath)
#         self.registration_parameters = self.image_config[fov_name]['fov_analysis_parameters']['rounds_registration']
#         self.chunk_size = self.registration_parameters['chunk_size']
#         self.min_dots_chunk = self.registration_parameters['min_dots_chunk']
#         self.min_error_triangles = self.registration_parameters['min_error_triangles']
#         self.percent_padding = self.registration_parameters['percent_padding']
#         self.reference_round = self.registration_parameters['reference_round']
#         self.collect_all_chunks = self.registration_parameters['collect_all_chunks']
#         self.reference_round_name = 'round_' + str(self.reference_round)
       
#         # The top lef coords are 0,0 because we are using relative coords
#         self.tl_coords = (0,0)

#         self.img_dimensions = (self.image_config[fov_name]['rounds'][self.reference_round_name]['shape']['height'],
#                                self.image_config[fov_name]['rounds'][self.reference_round_name]['shape']['width'])    
                    

#     @staticmethod
#     def combine_coords(counts, round_num):
#         data_reference = counts.loc[counts['round_num'] == round_num]
#         r_px = data_reference.r_px_original.to_list()
#         c_px = data_reference.c_px_original.to_list()
#         coords = np.array(list(zip(r_px,c_px)))
#         position_idx = data_reference.index
#         return coords, position_idx


#     @staticmethod
#     def obj_fun(pars,x,src):
#         tx, ty = pars
#         H = np.array([[1, 0, tx],\
#             [0, 1, ty]])
#         src1 = np.c_[src,np.ones(src.shape[0])]
#         return np.sum( (x - src1.dot(H.T)[:,:2])**2 )

#     @staticmethod
#     def apply_transform(pars, src):
#         tx, ty = pars
#         H = np.array([[1, 0, tx],\
#             [0, 1, ty]])
#         src1 = np.c_[src,np.ones(src.shape[0])]
#         return src1.dot(H.T)[:,:2]

#     @staticmethod
#     def distance(x1,y1,x2,y2):
#         return math.sqrt((x2 - x1)**2 + (y2 - y1)**2 )

#     @staticmethod
#     def list_subtract(list1,list2):
#         return np.absolute(np.array(list1)-np.array(list2))

#     def tri_sides(self,set_x, set_x_tri):

#         triangles = []
#         for i in range(len(set_x_tri)):

#             point1 = set_x_tri[i][0]
#             point2 = set_x_tri[i][1]
#             point3 = set_x_tri[i][2]

#             point1x, point1y = set_x[point1][0], set_x[point1][1]
#             point2x, point2y = set_x[point2][0], set_x[point2][1]
#             point3x, point3y = set_x[point3][0], set_x[point3][1] 

#             len1 = self.distance(point1x,point1y,point2x,point2y)
#             len2 = self.distance(point1x,point1y,point3x,point3y)
#             len3 = self.distance(point2x,point2y,point3x,point3y)

#             # you need to normalize in case the ref and the tran
#             # are warped
#             #min_side = min(len1,len2,len3)
#             #len1/=min_side
#             #len2/=min_side
#             #len3/=min_side
#             t=[len1,len2,len3]
#             t.sort()
#             triangles.append(t)

#         return triangles


#     def identify_matching_coords(self,set_A, set_B, threshold):
#         match_A_pts = []
#         match_B_pts = []
#         set_A_tri = list(itertools.combinations(range(len(set_A)), 3))
#         set_B_tri = list(itertools.combinations(range(len(set_B)), 3))
#         A_triangles = self.tri_sides(set_A, set_A_tri)
#         B_triangles = self.tri_sides(set_B, set_B_tri)
#         sums = []
#         for i in range(len(A_triangles)):
#             for j in range(len(B_triangles)):
#                 k = sum(self.list_subtract(A_triangles[i], B_triangles[j]))
#                 if k < threshold:
#                     sums.append([i,j,k])
#         # sort by smallest sum
#         sums = sorted(sums, key=operator.itemgetter(2))
#         if len(sums):
#             match_A = set_A_tri[sums[0][0]]
#             match_B = set_B_tri[sums[0][1]]
#             for i in range(3):
#                 match_A_pts.append(set_A[match_A[i]])
#                 match_B_pts.append(set_B[match_B[i]])
#         return (match_A_pts,match_B_pts)



#     def calculate_chunks(self):
#         self.chunks = chunking(self.img_dimensions,self.chunk_size,
#                           self.registration_parameters['percent_padding'],self.tl_coords)   
#         self.chunks.block_chunking()
#         self.Coords_Padded_Chunks_list = self.chunks.Coords_Padded_Chunks_list
        
#     def calculate_dots_chunks(self,coords,chunk_coords):  
#         r_tl = chunk_coords[0]
#         r_br = chunk_coords[1]
#         c_tl = chunk_coords[2]
#         c_br = chunk_coords[3]

#         # Select only the coords in the trimmed region
#         coords_in_chunk = coords[((r_tl < coords[:,0]) & (coords[:,0]<r_br)\
#                     & (c_tl <coords[:,1]) &(coords[:,1]<c_br)),: ]
#         return coords_in_chunk


#     def optimize_chunking(self,ref_coords, tran_coords):       
#         self.enough_dots = False
#         chunk_size = self.chunk_size
#         while chunk_size < min(self.img_dimensions):
#             chunks = chunking(self.img_dimensions, chunk_size, self.percent_padding, self.tl_coords)
#             chunks.block_chunking()
#             Coords_Padded_Chunks_list = chunks.Coords_Padded_Chunks_list
#             ref_max_number_dots = []
#             tran_max_number_dots = []
#             ref_total = []
#             tran_total = []
#             for chunk_coords in Coords_Padded_Chunks_list:
#                 ref_coords_in_chunk = self.calculate_dots_chunks(ref_coords,chunk_coords)
#                 tran_coords_in_chunk = self.calculate_dots_chunks(tran_coords,chunk_coords)
#                 if ref_coords_in_chunk.shape[0] > self.min_dots_chunk and tran_coords_in_chunk.shape[0] > self.min_dots_chunk:
#                         self.enough_dots = True
#                         break
#             if self.enough_dots:
#                 break
#             else:
#                 self.enough_dots = False
#                 chunk_size += 200

#         if self.enough_dots:
#             # Collect the ref and tran coords from the chunks with enough dots
#             self.ref_tran_screening_list = []
#             for chunk_coords in Coords_Padded_Chunks_list:
#                 ref_coords_in_chunk = self.calculate_dots_chunks(ref_coords,chunk_coords)
#                 tran_coords_in_chunk = self.calculate_dots_chunks(tran_coords,chunk_coords)
#                 if ref_coords_in_chunk.shape[0] > self.min_dots_chunk and tran_coords_in_chunk.shape[0] > self.min_dots_chunk:
#                     self.ref_tran_screening_list.append((ref_coords_in_chunk,tran_coords_in_chunk,chunk_coords))
#             self.chunk_size = chunk_size


#     def rounds_to_register(self):
#         # remember that round_num is not pythonic
#         self.rounds_list = list(self.counts['round_num'].unique())
#         if self.reference_round == False:
#             logger.error(f'missing reference round')
#             sys.exit(f'missing reference round')
#         elif self.reference_round <0 or self.reference_round > max(self.rounds_list):
#             logger.error(f'selected reference round {self.reference_round} is out of range')
#             sys.exit(f'selected reference round {self.reference_round} is out of range')
#         else:
#             self.rounds_list = list(self.counts['round_num'].unique()) 
#             self.processing_combinations = []
#             self.rounds_list.remove(self.reference_round)
#             for el in self.rounds_list:
#                 self.processing_combinations.append((self.reference_round,el))

#     def cpl_registration(self,ref_coords,tran_coords,cpl):
#         self.optimize_chunking(ref_coords, tran_coords)
#         self.completed_registration = False
#         if self.enough_dots:
#             match_ref_pts_all = []
#             match_tran_pts_all = []
#             if self.collect_all_chunks:
#             # Collect all matching dots in all chunked regions with number of dots above threshold
#                 for ref_coords_in_chunk,tran_coords_in_chunk, chunk_coords in self.ref_tran_screening_list:
#                         match_ref_pts, match_tran_pts = self.identify_matching_coords(ref_coords_in_chunk,tran_coords_in_chunk,self.min_error_triangles)
#                         if len(match_ref_pts) and len(match_tran_pts):
#                             match_ref_pts_all.append(match_ref_pts)
#                             match_tran_pts_all.append(match_tran_pts)
#                 match_ref_pts_all = [pts for grp in match_ref_pts_all for pts in grp]
#                 match_tran_pts_all = [pts for grp in match_tran_pts_all for pts in grp]
#             else:
#                 ref_coords_in_chunk,tran_coords_in_chunk, chunk_coords = self.ref_tran_screening_list[0]
#                 match_ref_pts_all, match_tran_pts_all = self.identify_matching_coords(ref_coords_in_chunk,tran_coords_in_chunk,self.min_error_triangles)

#             if len(match_ref_pts_all):
#                 match_ref_pts_all = np.vstack(match_ref_pts_all)
#                 match_tran_pts_all = np.vstack(match_tran_pts_all)
#                 minimization_output = minimize(self.obj_fun,[0,0],args=(match_ref_pts_all,match_tran_pts_all), method='Nelder-Mead')
#                 if minimization_output.success:
#                     self.tran_registered_coords = self.apply_transform(minimization_output.x, tran_coords)
#                     self.transformation_matrix = minimization_output.x
#                     self.completed_registration = True
#                 else:
#                     self.logger.info(f'chunk {chunk_coords} of {cpl} failed minimization of distances')
#             else:
#                 self.logger.info(f'chunk {chunk_coords} of {cpl} did not find matching triangles')

#         else:
#             self.logger.info(f'cannot register rounds {cpl} not enough dots')
#             self.tran_registered_coords = tran_coords
#             self.transformation_matrix = np.empty([1,2])
#             self.transformation_matrix[:] = np.nan

#         if not self.completed_registration:
#             self.logger.info(f'was not possible to register {cpl} ')
#             self.tran_registered_coords = tran_coords
#             self.transformation_matrix = np.empty([1,2])
#             self.transformation_matrix[:] = np.nan


#     def register(self):
#         all_registration_matrix = {}
#         r_px_col_name = 'r_px_registered'
#         c_px_col_name = 'c_px_registered'
#         self.counts[r_px_col_name] = np.nan
#         self.counts[c_px_col_name] = np.nan
#         self.counts['r_transformation_registration'] = np.nan
#         self.counts['c_transformation_registration'] = np.nan
#         self.rounds_to_register()
#         for cpl in self.processing_combinations:
#             ref_round, tran_round = cpl
#             ref_coords, ref_position_idx = self.combine_coords(self.counts,ref_round)
#             tran_coords, tran_position_idx = self.combine_coords(self.counts,tran_round)
#             self.cpl_registration(ref_coords,tran_coords,cpl)
#             all_registration_matrix[cpl] = self.transformation_matrix
#             self.ref_registered_coords = ref_coords
#             self.counts.loc[ref_position_idx, r_px_col_name] = self.ref_registered_coords[:,0]
#             self.counts.loc[ref_position_idx, c_px_col_name] = self.ref_registered_coords[:,1]
#             self.counts.loc[ref_position_idx, 'r_transformation_registration'] = np.ones(self.ref_registered_coords.shape[0])
#             self.counts.loc[ref_position_idx, 'c_transformation_registration'] = np.ones(self.ref_registered_coords.shape[0])
#             self.counts.loc[tran_position_idx, r_px_col_name] = self.tran_registered_coords[:,0]
#             self.counts.loc[tran_position_idx, c_px_col_name] = self.tran_registered_coords[:,1]
#             all_transf = np.tile(self.transformation_matrix,(self.tran_registered_coords.shape[0],1))
#             self.counts.loc[tran_position_idx, 'r_transformation_registration'] = all_transf[:,0]
#             self.counts.loc[tran_position_idx, 'c_transformation_registration'] = all_transf[:,1]
       
#         return self.counts, all_registration_matrix




# TODO Remove not used functions
########## ----------------------------------------





# @task(name='fft_registration_beads')# ADD LOGGER
def fft_registration_beads(reference_coords:np.ndarray, translated_coords:np.ndarray,
                       img_shape: np.ndarray, fov_num:int,
                       hybridization_num_translated:int):

    # CHECK THE VALUES AND CATCH THE ERROR TO BLOCK FOV FOR REGISTRATION     
    if np.any(np.isnan(reference_coords)):
        logger.error(f'missing reference round for registration for fov {fov_num}') 
        signals.SKIP(f'missing reference round for registration for fov {fov_num}')
        shift = np.array([np.nan,np.nan])
        error = 0
        tran_registered_coords = translated_coords
        
        # RETURN VALUES TO WRITE ON DB
    elif np.any(np.isnan(translated_coords)):
        logger.error(f'missing registration round {hybridization_num_translated} for registration for fov {fov_num}')   
        shift = np.array([np.nan,np.nan])
        error = 0
        tran_registered_coords = translated_coords
        
    else:
        img_ref = create_fake_image(img_shape,reference_coords)    
        img_tran = create_fake_image(img_shape,translated_coords)
        shift, error, diffphase = register_translation(img_ref, img_tran)
        tran_registered_coords = translated_coords + shift
        tran_registered_coords = tran_registered_coords.astype(int)

    return tran_registered_coords, shift, error, fov_num, hybridization_num_translated 


def identify_matching_register_dots_NN(ref_dots_coords_fov,tran_registered_coords,registration_tollerance_pxl):
    # put in the wrapping function
    #if (ref_dots_coords_fov.shape[0] >0) and (tran_registered_coords.shape[0] >0):
            
    # initialize network
    nn = NearestNeighbors(1, metric="euclidean")
    nn.fit(ref_dots_coords_fov)

    # Get the nn
    dists, indices = nn.kneighbors(tran_registered_coords, return_distance=True)

    # select only the nn that are below pxl distance
    idx_selected_coords_compare = np.where(dists <= registration_tollerance_pxl)[0]

    number_matching_dots = idx_selected_coords_compare.shape[0]
    
    return number_matching_dots


def calculate_shift_hybridization_fov(processing_files:List,analysis_parameters:dict, save=True):
    """
    Function used to run the registration of a single fov
    through all hybridization. The registration is done using the dots
    coords determined in each image

    Args:
        processing_files: List
            list of the files with the counts to register
        analysis_parameters: Dict
            dict that contains the settings for the analysis 
    """
    
    logger = selected_logger()
    all_rounds_shifts = {}
    
    registration_errors = Registration_errors()
    data_models = Output_models()
    output_registration_df = data_models.output_registration_df
    status = 'SUCCESS'

    reference_hybridization = analysis_parameters['RegistrationReferenceHybridization']
    round_num = reference_hybridization
    reference_hybridization_str = 'Hybridization' + str(reference_hybridization).zfill(2)
    registration_tollerance_pxl = analysis_parameters['RegistrationTollerancePxl']

    # collect info used to generate fname when files are corrupted   
    experiment_fpath = processing_files[0].parent.parent.parent
    channel = (processing_files[0].stem).split('_')[-4]
    experiment_name = experiment_fpath.stem
    fov = (processing_files[0].stem).split('_')[-2]
    file_tags = {'experiment_fpath':experiment_fpath,
                'experiment_name':experiment_name,
                'channel':channel,
                'fov':fov}

    fname = experiment_fpath / 'tmp' / 'registered_counts' / (experiment_name + '_' + channel + '_registered_fov_' + fov + '.parquet')
    shift_fname = experiment_fpath / 'tmp' / 'registered_counts' / (experiment_name + '_' + channel + '_all_rounds_shifts_fov_' + fov + '.pkl')
    # Load reference hybridization data
    try:
        ref_fpath = [fpath for fpath in processing_files if reference_hybridization_str in fpath.as_posix()][0]
    except:
        logger.error(f'missing reference hyb file for fov {fov}')
        output_registration_df = output_registration_df.append({'min_number_matching_dots_registration':registration_errors.missing_file_reg_channel, 
                                            'fov_num':int(fov),'dot_channel':channel,'round_num':int(round_num)},ignore_index=True)
        status = 'FAILED'
        all_rounds_shifts[round_num] = np.array([np.nan,np.nan])
    else:
        try:
            ref_counts,ref_img_metadata = pickle.load(open(ref_fpath, 'rb'))
        except:
            logger.error(f'cannot open the reference hyb file for fov {fov}')
            output_registration_df = output_registration_df.append({'min_number_matching_dots_registration':registration_errors.cannot_load_file_reg_channel,
                                            'fov_num':int(fov),'dot_channel':channel,'round_num':int(round_num) },ignore_index=True)
            status = 'FAILED'
            all_rounds_shifts[round_num] = np.array([np.nan,np.nan])
        else:
            # Check if there are dots detected in the reference round
            if np.any(np.isnan(ref_counts['r_px_original'])) or (ref_counts['r_px_original'].shape[0]<registration_tollerance_pxl):
                logger.error(f'There are no dots in there reference hyb for fov {fov} ')
                output_registration_df = output_registration_df.append({'min_number_matching_dots_registration':registration_errors.missing_counts_reg_channel,
                                            'fov_num':int(fov),'dot_channel':channel,'round_num':int(round_num) },ignore_index=True)
                status = 'FAILED'
                all_rounds_shifts[round_num] = np.array([np.nan,np.nan])
            else:
                ref_counts_df = pd.DataFrame(ref_counts)
                ref_counts_df['r_px_registered'] = ref_counts_df['r_px_original']
                ref_counts_df['c_px_registered'] = ref_counts_df['c_px_original']
                ref_counts_df['r_shift_px'] = 0
                ref_counts_df['c_shift_px'] = 0
                ref_counts_df['min_number_matching_dots_registration'] = 1000
                all_rounds_shifts[round_num] = np.array([0,0])

                output_registration_df = pd.concat([output_registration_df,ref_counts_df],axis=0,ignore_index=True)
                
                tran_processing_files =processing_files.copy()
                tran_processing_files.remove(ref_fpath)
                img_width = ref_img_metadata['img_width']
                img_height = ref_img_metadata['img_height']

                img_shape = (img_width, img_height)
                ref_coords = ref_counts_df.loc[:,['r_px_original', 'c_px_original']].to_numpy()
                img_ref = create_fake_image(img_shape,ref_coords)
        

                for fpath in tran_processing_files:
                    try:
                        tran_counts, tran_img_metadata = pickle.load(open(fpath, 'rb'))
                    except:
                        logger.error(f'cannot open {fpath.stem} file for fov {fov}')
                        # If there is an error in the opening reset the df
                        round_num = int((fpath.stem).split('_')[-5].split('Hybridization')[-1])
                        output_registration_df = data_models.output_registration_df
                        output_registration_df = output_registration_df.append({'min_number_matching_dots_registration':registration_errors.cannot_load_file_reg_channel,
                                            'fov_num':int(fov) ,'dot_channel':channel,'round_num':int(round_num)},ignore_index=True)
                        status = 'FAILED'
                        all_rounds_shifts[round_num] = np.array([np.nan,np.nan])
                        break
                    else:
                        tran_counts_df = pd.DataFrame(tran_counts)
                        tran_coords = tran_counts_df.loc[:,['r_px_original', 'c_px_original']].to_numpy()
                        if np.any(np.isnan(tran_coords)) or tran_coords.shape[0]<registration_tollerance_pxl:
                            round_num = int((fpath.stem).split('_')[-5].split('Hybridization')[-1])
                            # If dots are missing, reset the df
                            output_registration_df = data_models.output_registration_df
                            output_registration_df = output_registration_df.append({'min_number_matching_dots_registration':registration_errors.missing_counts_reg_channel,
                                            'fov_num':int(fov),'dot_channel':channel,'round_num':int(round_num) },ignore_index=True)
                            status = 'FAILED'
                            all_rounds_shifts[round_num] = np.array([np.nan,np.nan])
                            logger.error(f' {fpath.stem} file for fov {fov} has no counts')
                            break
                        else:
                            
                            round_num = tran_counts['round_num'][0]
                            ref_img_metadata['reference_hyb'] = str(reference_hybridization)
                            img_tran = create_fake_image(img_shape,tran_coords)
                            shift, error, diffphase = register_translation(img_ref, img_tran)
                            all_rounds_shifts[round_num] = shift
                            registered_tran_coords = tran_coords + shift
                            min_num_matching_dots = identify_matching_register_dots_NN(ref_coords,
                                                                                    registered_tran_coords,
                                                                                    registration_tollerance_pxl)
                            
                            tran_counts_df['r_px_registered'] = registered_tran_coords[:,0]
                            tran_counts_df['c_px_registered'] = registered_tran_coords[:,1]
                            tran_counts_df['r_shift_px'] = shift[0]
                            tran_counts_df['c_shift_px'] = shift[1]
                            tran_counts_df['min_number_matching_dots_registration'] = min_num_matching_dots

                            output_registration_df = pd.concat([output_registration_df,tran_counts_df],axis=0,ignore_index=True)
                            output_registration_df.attrs[fov] = ref_img_metadata

                # Save the dataframe
                
                output_registration_df['reference_hyb'] = reference_hybridization
                output_registration_df['experiment_type'] = ref_img_metadata['experiment_type']
                output_registration_df['experiment_name'] = ref_img_metadata['experiment_name']
                output_registration_df['pxl_um'] = ref_img_metadata['pixel_microns']
                output_registration_df['stitching_type'] = ref_img_metadata['stitching_type']
                output_registration_df['img_width_px'] = ref_img_metadata['img_width']
                output_registration_df['img_height_px'] = ref_img_metadata['img_height']
                output_registration_df['fov_acquisition_coords_x'] = ref_img_metadata['fov_acquisition_coords_x']
                output_registration_df['fov_acquisition_coords_y'] = ref_img_metadata['fov_acquisition_coords_y']
                output_registration_df['fov_acquisition_coords_z'] = ref_img_metadata['fov_acquisition_coords_z']
        
        # Save extra metadata in the
        if save:
            output_registration_df.to_parquet(fname,index=False)
            pickle.dump(all_rounds_shifts,open(shift_fname,'wb'))
        return output_registration_df, all_rounds_shifts, file_tags, status



def calculate_shift_hybridization_fov_test(fov:int,
                                            counts_output:Dict,
                                            analysis_parameters:Dict, 
                                            experiment_info:Dict):
    """
    Function used to run the registration of a single fov
    through all hybridization. The registration is done using the dots
    coords determined in each image

    Args:
        processing_files: List
            list of the files with the counts to register
        analysis_parameters: Dict
            dict that contains the settings for the analysis 
    """
    
    logger = selected_logger()
    all_rounds_shifts = {}
    
    registration_errors = Registration_errors()
    data_models = Output_models()
    output_registration_df = data_models.output_registration_df
    status = 'SUCCESS'

    reference_hybridization = analysis_parameters['RegistrationReferenceHybridization']
    round_num = reference_hybridization
    reference_hybridization_str = 'Hybridization' + str(reference_hybridization).zfill(2)
    channel = experiment_info['StitchingChannel']
    registration_tollerance_pxl = analysis_parameters['RegistrationTollerancePxl']
    counts_zarr_names = list(counts_output['registration'][channel].keys())

    try:
        ref_zarr_name = [el for el in counts_zarr_names if reference_hybridization_str in el][0]
    except:
        logger.error(f'missing registration counts fov {fov}')
        output_registration_df = output_registration_df.append({'min_number_matching_dots_registration':registration_errors.missing_file_reg_channel, 
                                            'fov_num':int(fov),'dot_channel':channel,'round_num':int(round_num)},ignore_index=True)
        status = 'FAILED'
        all_rounds_shifts[round_num] = np.array([np.nan,np.nan])

    else:
        ref_counts, ref_img_metadata = counts_output['registration'][channel][ref_zarr_name]
        # Check if there are dots detected in the reference round
        if np.any(np.isnan(ref_counts['r_px_original'])) or (ref_counts['r_px_original'].shape[0]<registration_tollerance_pxl):
            logger.error(f'There are no dots in there reference hyb for fov {fov} ')
            output_registration_df = output_registration_df.append({'min_number_matching_dots_registration':registration_errors.missing_counts_reg_channel,
                                        'fov_num':int(fov),'dot_channel':channel,'round_num':int(round_num) },ignore_index=True)
            status = 'FAILED'
            all_rounds_shifts[round_num] = np.array([np.nan,np.nan])
        else:
            ref_counts_df = pd.DataFrame(ref_counts)
            ref_counts_df['r_px_registered'] = ref_counts_df['r_px_original']
            ref_counts_df['c_px_registered'] = ref_counts_df['c_px_original']
            ref_counts_df['r_shift_px'] = 0
            ref_counts_df['c_shift_px'] = 0
            ref_counts_df['min_number_matching_dots_registration'] = 1000
            all_rounds_shifts[round_num] = np.array([0,0])

            output_registration_df = pd.concat([output_registration_df,ref_counts_df],axis=0,ignore_index=True)
            
            tran_processing_files = counts_zarr_names.copy()
            tran_processing_files.remove(ref_zarr_name)
            img_width = ref_img_metadata['img_width']
            img_height = ref_img_metadata['img_height']

            img_shape = (img_width, img_height)
            ref_coords = ref_counts_df.loc[:,['r_px_original', 'c_px_original']].to_numpy()
            img_ref = create_fake_image(img_shape,ref_coords)
    

            for zarr_name in tran_processing_files:
                round_num = int(zarr_name.split('_')[-4].split('Hybridization')[-1])
                try:
                    tran_counts, tran_img_metadata = counts_output['registration'][channel][zarr_name]
                except:
                    logger.error(f'cannot open {zarr_name} file for fov {fov}')
                    # If there is an error in the opening reset the df
                    output_registration_df = data_models.output_registration_df
                    output_registration_df = output_registration_df.append({'min_number_matching_dots_registration':registration_errors.cannot_load_file_reg_channel,
                                        'fov_num':int(fov) ,'dot_channel':channel,'round_num':int(round_num)},ignore_index=True)
                    status = 'FAILED'
                    all_rounds_shifts[round_num] = np.array([np.nan,np.nan])
                    break
                else:
                    tran_counts_df = pd.DataFrame(tran_counts)
                    tran_coords = tran_counts_df.loc[:,['r_px_original', 'c_px_original']].to_numpy()
                    if np.any(np.isnan(tran_coords)) or tran_coords.shape[0]<registration_tollerance_pxl:
                        # If dots are missing, reset the df
                        output_registration_df = data_models.output_registration_df
                        output_registration_df = output_registration_df.append({'min_number_matching_dots_registration':registration_errors.missing_counts_reg_channel,
                                        'fov_num':int(fov),'dot_channel':channel,'round_num':int(round_num) },ignore_index=True)
                        status = 'FAILED'
                        all_rounds_shifts[round_num] = np.array([np.nan,np.nan])
                        logger.error(f' {zarr_name} file for fov {fov} has no counts')
                        break
                    else:
                        
                        round_num = tran_counts['round_num'][0]
                        ref_img_metadata['reference_hyb'] = str(reference_hybridization)
                        img_tran = create_fake_image(img_shape,tran_coords)
                        shift, error, diffphase = register_translation(img_ref, img_tran)
                        all_rounds_shifts[round_num] = shift
                        registered_tran_coords = tran_coords + shift
                        min_num_matching_dots = identify_matching_register_dots_NN(ref_coords,
                                                                                registered_tran_coords,
                                                                                registration_tollerance_pxl)
                        
                        tran_counts_df['r_px_registered'] = registered_tran_coords[:,0]
                        tran_counts_df['c_px_registered'] = registered_tran_coords[:,1]
                        tran_counts_df['r_shift_px'] = shift[0]
                        tran_counts_df['c_shift_px'] = shift[1]
                        tran_counts_df['min_number_matching_dots_registration'] = min_num_matching_dots

                        output_registration_df = pd.concat([output_registration_df,tran_counts_df],axis=0,ignore_index=True)

            # Save the dataframe
            
            output_registration_df['reference_hyb'] = reference_hybridization
            output_registration_df['experiment_type'] = ref_img_metadata['experiment_type']
            output_registration_df['experiment_name'] = ref_img_metadata['experiment_name']
            output_registration_df['pxl_um'] = ref_img_metadata['pixel_microns']
            output_registration_df['stitching_type'] = ref_img_metadata['stitching_type']
            output_registration_df['img_width_px'] = ref_img_metadata['img_width']
            output_registration_df['img_height_px'] = ref_img_metadata['img_height']
            output_registration_df['fov_acquisition_coords_x'] = ref_img_metadata['fov_acquisition_coords_x']
            output_registration_df['fov_acquisition_coords_y'] = ref_img_metadata['fov_acquisition_coords_y']
            output_registration_df['fov_acquisition_coords_z'] = ref_img_metadata['fov_acquisition_coords_z']
    
    
    return output_registration_df, all_rounds_shifts, status


def register_fish(processing_files:List,analysis_parameters:Dict,
                        registered_reference_channel_df,all_rounds_shifts:Dict,file_tags:Dict,status:str,
                        save=True):

    logger = selected_logger()
    registration_errors = Registration_errors()
    data_models = Output_models()
    fov = file_tags['fov']
    channel = (processing_files[0].stem).split('_')[-4]
    file_tags['channel'] = channel
    registered_fish_df = data_models.output_registration_df

    if status == 'FAILED':
        error = registered_reference_channel_df['min_number_matching_dots_registration'].values[0]
        round_num = registered_reference_channel_df['round_num'].values[0]
        registered_fish_df = registered_fish_df.append({'min_number_matching_dots_registration':error,
                                                           'fov_num':int(fov),'dot_channel':channel,'round_num':int(round_num) },ignore_index=True)
    elif status == 'SUCCESS':
        reference_hybridization = registered_reference_channel_df.attrs[fov]['reference_hyb']
        reference_hybridization_str = 'Hybridization' + str(reference_hybridization).zfill(2)
        for fpath in processing_files:
            round_num = int((fpath.stem).split('_')[-5].split('Hybridization')[-1])
            try:
                fish_counts, fish_img_metadata = pickle.load(open(fpath,'rb'))
            except:
                logger.error(f'cannot open the processing files')
                
                registered_fish_df = registered_fish_df.append({'min_number_matching_dots_registration':registration_errors.cannot_load_file_fish_channel,
                                                'fov_num':int(fov),'dot_channel':channel,'round_num':round_num },ignore_index=True)
                status = 'FAILED'
                break
            else:                
                if np.any(np.isnan(fish_counts['r_px_original'])):
                    logger.error(f'There are no dots in there reference hyb for fov {fov} ')
                    registered_fish_df = registered_fish_df.append({'min_number_matching_dots_registration':registration_errors.missing_counts_fish_channel,
                                                'fov_num':int(fov),'dot_channel':channel,'round_num':round_num },ignore_index=True)
                    status = 'FAILED'
                    break

                else:
                    if reference_hybridization_str in fpath.as_posix():
                        fish_img_metadata['reference_hyb'] = reference_hybridization
                        registered_fish_df.attrs[fov] = fish_img_metadata
                    fish_counts_df = pd.DataFrame(fish_counts)
                    
                    subset_df = registered_reference_channel_df.loc[registered_reference_channel_df['round_num'] == round_num, :]
                    subset_df = subset_df.reset_index()
                        
                    shift = all_rounds_shifts[round_num]
                    fish_counts_df['r_px_registered'] = fish_counts['r_px_original'] + shift[0]
                    fish_counts_df['c_px_registered'] = fish_counts['c_px_original'] + shift[1]
                    fish_counts_df['r_shift_px'] = shift[0]
                    fish_counts_df['c_shift_px'] = shift[1]
                    fish_counts_df['min_number_matching_dots_registration'] = subset_df.loc[0,'min_number_matching_dots_registration'] 
                    # fish_counts_df['reference_hyb'] = reference_hybridization
                    # fish_counts_df['experiment_type'] = subset_df.loc[0,'experiment_type']
                    # fish_counts_df['experiment_name'] = subset_df.loc[0,'experiment_name']
                    # fish_counts_df['pxl_um'] = subset_df.loc[0,'pxl_um']
                    # fish_counts_df['stitching_type'] = subset_df.loc[0,'stitching_type']
                    # fish_counts_df['fov_acquisition_coords_x'] = subset_df.loc[0,'fov_acquisition_coords_x']
                    # fish_counts_df['fov_acquisition_coords_y'] = subset_df.loc[0,'fov_acquisition_coords_y']
                    # fish_counts_df['fov_acquisition_coords_z'] = subset_df.loc[0,'fov_acquisition_coords_z']

                    registered_fish_df = pd.concat([registered_fish_df,fish_counts_df],axis=0,ignore_index=True)
                    status = 'SUCCESS'

    if save:
        fname = file_tags['experiment_fpath'] / 'tmp' / 'registered_counts' / (file_tags['experiment_name'] + '_' + file_tags['channel'] + '_registered_fov_' + file_tags['fov'] + '.parquet')
        registered_fish_df.to_parquet(fname,index=False)
    return registered_fish_df, file_tags, status


def register_fish_test(fov:int,
                        channel:str,
                        counts_output:Dict,
                        registered_reference_channel_df,
                        all_rounds_shifts:Dict,
                        analysis_parameters:Dict,
                        status:str):

    logger = selected_logger()
    registration_errors = Registration_errors()
    data_models = Output_models()

    registered_fish_df = data_models.output_registration_df

    counts_zarr_names = list(counts_output['fish'][channel].keys())

    if status == 'FAILED':
        error = registered_reference_channel_df['min_number_matching_dots_registration'].values[0]
        round_num = registered_reference_channel_df['round_num'].values[0]
        registered_fish_df = registered_fish_df.append({'min_number_matching_dots_registration':error,
                                                           'fov_num':int(fov),'dot_channel':channel,'round_num':int(round_num) },ignore_index=True)
    elif status == 'SUCCESS':
        reference_hybridization = analysis_parameters['RegistrationReferenceHybridization']
        reference_hybridization_str = 'Hybridization' + str(reference_hybridization).zfill(2)
        
        for zarr_name in counts_zarr_names:
            round_num = int(zarr_name.split('_')[-4].split('Hybridization')[-1])
            try:
                fish_counts, fish_img_metadata = counts_output['fish'][channel][zarr_name]
            except:
                logger.error(f'missing the counts in {zarr_name}')
                
                registered_fish_df = registered_fish_df.append({'min_number_matching_dots_registration':registration_errors.cannot_load_file_fish_channel,
                                                'fov_num':int(fov),'dot_channel':channel,'round_num':round_num },ignore_index=True)
                status = 'FAILED'
                break
            else:                
                if np.any(np.isnan(fish_counts['r_px_original'])):
                    logger.error(f'There are no dots in there reference hyb for fov {fov} ')
                    registered_fish_df = registered_fish_df.append({'min_number_matching_dots_registration':registration_errors.missing_counts_fish_channel,
                                                'fov_num':int(fov),'dot_channel':channel,'round_num':round_num },ignore_index=True)
                    status = 'FAILED'
                    break

                else:
                    if reference_hybridization_str in zarr_name:
                        fish_img_metadata['reference_hyb'] = reference_hybridization
                        registered_fish_df.attrs[fov] = fish_img_metadata
                    fish_counts_df = pd.DataFrame(fish_counts)
                    
                    subset_df = registered_reference_channel_df.loc[registered_reference_channel_df['round_num'] == round_num, :]
                    subset_df = subset_df.reset_index()
                        
                    shift = all_rounds_shifts[round_num]
                    fish_counts_df['r_px_registered'] = fish_counts['r_px_original'] + shift[0]
                    fish_counts_df['c_px_registered'] = fish_counts['c_px_original'] + shift[1]
                    fish_counts_df['r_shift_px'] = shift[0]
                    fish_counts_df['c_shift_px'] = shift[1]
                    fish_counts_df['min_number_matching_dots_registration'] = subset_df.loc[0,'min_number_matching_dots_registration'] 
                    fish_counts_df['reference_hyb'] = reference_hybridization
                    fish_counts_df['experiment_type'] = subset_df.loc[0,'experiment_type']
                    fish_counts_df['experiment_name'] = subset_df.loc[0,'experiment_name']
                    fish_counts_df['pxl_um'] = subset_df.loc[0,'pxl_um']
                    fish_counts_df['stitching_type'] = subset_df.loc[0,'stitching_type']
                    fish_counts_df['fov_acquisition_coords_x'] = subset_df.loc[0,'fov_acquisition_coords_x']
                    fish_counts_df['fov_acquisition_coords_y'] = subset_df.loc[0,'fov_acquisition_coords_y']
                    fish_counts_df['fov_acquisition_coords_z'] = subset_df.loc[0,'fov_acquisition_coords_z']
                    fish_counts_df['img_width_px'] = fish_img_metadata['img_width']
                    fish_counts_df['img_height_px'] = fish_img_metadata['img_height']




                    registered_fish_df = pd.concat([registered_fish_df,fish_counts_df],axis=0,ignore_index=True)
                    status = 'SUCCESS'

    return registered_fish_df, status






def calculate_shift_hybridization_fov_nuclei(processing_files:List,analysis_parameters:dict, save=True):
    """
    Function used to run the registration of a single fov
    through all hybridization. The registration is done using the dots
    coords determined in each image

    Args:
        processing_files: List
            list of the files with the filtered images of the nuclei
        analysis_parameters: Dict
            dict that contains the settings for the analysis 
    """
    
    logger = selected_logger()
    all_rounds_shifts = {}
    all_rounds_shifts_RMS = {}
    registration_errors = Registration_errors()
    data_models = Output_models()
    output_registration_df = data_models.output_registration_df
    status = 'SUCCESS'

    reference_hybridization = analysis_parameters['RegistrationReferenceHybridization']
    round_num = reference_hybridization
    reference_hybridization_str = 'Hybridization' + str(reference_hybridization).zfill(2)
    registration_tollerance_pxl = analysis_parameters['RegistrationTollerancePxl']

    # collect info used to generate fname when files are corrupted   
    experiment_fpath = processing_files[0].parent.parent.parent
    channel = (processing_files[0].stem).split('_')[-4]
    experiment_name = experiment_fpath.stem
    fov = (processing_files[0].stem).split('_')[-2]
    file_tags = {'experiment_fpath':experiment_fpath,
                'experiment_name':experiment_name,
                'channel':channel,
                'fov':fov}

    shift_fname = experiment_fpath / 'tmp' / 'registered_counts' / (experiment_name + '_' + channel + '_nuclei_all_rounds_shifts_fov_' + fov + '.pkl')
    shift_error_fname = experiment_fpath / 'tmp' / 'registered_counts' / (experiment_name + '_' + channel + '_nuclei_all_rounds_shifts_errors_fov_' + fov + '.pkl')
    
    # Load reference hybridization data
    try:
        ref_fpath = [fpath for fpath in processing_files if reference_hybridization_str in fpath.as_posix()][0]
    except:
        logger.error(f'missing reference hyb file for fov {fov}')
        status = 'FAILED'
        all_rounds_shifts[round_num] = np.array([np.nan,np.nan])
        all_rounds_shifts_RMS[round_num] = 1
    else:
        try:
            ref_img,ref_img_metadata = pickle.load(open(ref_fpath, 'rb'))
        except:
            logger.error(f'cannot open the reference hyb file for fov {fov}')
            status = 'FAILED'
            all_rounds_shifts[round_num] = np.array([np.nan,np.nan])
            all_rounds_shifts_RMS[round_num] = 1
        else:   
                tran_processing_files =processing_files.copy()
                tran_processing_files.remove(ref_fpath)
                
                for fpath in tran_processing_files:
                    try:
                        tran_img, tran_img_metadata = pickle.load(open(fpath, 'rb'))
                    except:
                        logger.error(f'cannot open {fpath.stem} file for fov {fov}')
                        # If there is an error in the opening reset the df
                        round_num = int((fpath.stem).split('_')[-5].split('Hybridization')[-1])
                        status = 'FAILED'
                        all_rounds_shifts[round_num] = np.array([np.nan,np.nan])
                        all_rounds_shifts_RMS[round_num] = 1
                        break
                    else:                    
                            round_num = int((fpath.stem).split('_')[-5].split('Hybridization')[-1])
                            ref_img_metadata['reference_hyb'] = str(reference_hybridization)
                            
                            shift, error, diffphase = phase_cross_correlation(ref_img, tran_img,return_error=True,)
                            all_rounds_shifts[round_num] = shift
                            all_rounds_shifts_RMS[round_num] = error
        
        # Save extra metadata in the
        if save:
            pickle.dump(all_rounds_shifts,open(shift_fname,'wb'))
            pickle.dump(all_rounds_shifts_RMS,open(shift_error_fname,'wb'))

        return all_rounds_shifts, all_rounds_shifts_RMS, file_tags, status






def calculate_shift_hybridization_fov_nuclei_test(
                                            img_stack:np.ndarray,
                                            analysis_parameters:Dict):
    """
    Function used to run the registration of a single fov
    through all hybridization. The registration is done using the dots
    coords determined in each image

    Args:
        processing_files: List
            list of the files with the counts to register
        analysis_parameters: Dict
            dict that contains the settings for the analysis 


    MUST ADD THE CONSIDERATION OF POTENTIAL ERRORS
    """
    
    logger = selected_logger()
    
    registration_errors = Registration_errors()
  
    ref_round_num = analysis_parameters['RegistrationReferenceHybridization'] - 1
    registration_tollerance_pxl = analysis_parameters['RegistrationTollerancePxl']

    all_rounds = np.arange(img_stack.shape[0])
    all_rounds = all_rounds[all_rounds != ref_round_num]

    ref_img = img_stack[ref_round_num,:,:]

    all_rounds_reg = []
    ref_df = pd.DataFrame({'round_num':ref_round_num + 1,
                                            'r_shift_px': 0,
                                            'c_shift_px':0,
                                            'min_number_matching_dots_registration':1000,
                                            'RMS':0})

    all_rounds_reg.append(ref_df)
    for r_num in all_rounds:
        shift, error, diffphase = phase_cross_correlation(ref_img, img_stack[r_num,:,:],return_error=True)    
        tran_df = pd.DataFrame({'round_num':[r_num + 1],
                                        'r_shift_px': [shift[0]],
                                        'c_shift_px':[shift[1]],
                                        'min_number_matching_dots_registration':error,
                                        'RMS':error})

        all_rounds_reg.append(tran_df)

    output_registration_df = pd.concat([all_rounds_reg],axis=0,ignore_index=True)


    return output_registration_df





def register_fish_on_nuclei_test(fov:int,
                        counts_output:Dict,
                        registered_reference_channel_df,
                        all_rounds_shifts:Dict,
                        analysis_parameters:Dict,
                        status:str):

    logger = selected_logger()
    registration_errors = Registration_errors()
    data_models = Output_models()

    registered_fish_df = data_models.output_registration_df

    if status == 'FAILED':
        error = registered_reference_channel_df['min_number_matching_dots_registration'].values[0]
        round_num = registered_reference_channel_df['round_num'].values[0]
        registered_fish_df = registered_fish_df.append({'min_number_matching_dots_registration':error,
                                                           'fov_num':int(fov),'dot_channel':channel,'round_num':int(round_num) },ignore_index=True)
    elif status == 'SUCCESS':
        reference_hybridization = analysis_parameters['RegistrationReferenceHybridization']
        reference_hybridization_str = 'Hybridization' + str(reference_hybridization).zfill(2)
        
        for channel, all_counts_dict in counts_output['fish'].items():
            for zarr_name, fish_counts_tpl in all_counts_dict.items():
                round_num = int(zarr_name.split('_')[-4].split('Hybridization')[-1])
                
                fish_counts, fish_img_metadata = fish_counts_tpl

                fish_counts_df = pd.DataFrame(fish_counts)
                shift = all_rounds_shifts[round_num]
                fish_counts_df['r_px_registered'] = fish_counts['r_px_original'] + shift[0]
                fish_counts_df['c_px_registered'] = fish_counts['c_px_original'] + shift[1]
                fish_counts_df['r_shift_px'] = shift[0]
                fish_counts_df['c_shift_px'] = shift[1]
                # fish_counts_df['min_number_matching_dots_registration'] = fish_counts['min_number_matching_dots_registration']
                fish_counts_df['reference_hyb'] = reference_hybridization
                # fish_counts_df['experiment_type'] = fish_counts['experiment_type']
                # fish_counts_df['experiment_name'] = fish_counts['experiment_name']
                # fish_counts_df['pxl_um'] = fish_counts['pxl_um']
                # fish_counts_df['stitching_type'] = fish_counts['stitching_type']
                # fish_counts_df['fov_acquisition_coords_x'] = fish_counts['fov_acquisition_coords_x']
                # fish_counts_df['fov_acquisition_coords_y'] = fish_counts['fov_acquisition_coords_y']
                # fish_counts_df['fov_acquisition_coords_z'] = fish_counts['fov_acquisition_coords_z']
                # fish_counts_df['img_width_px'] = fish_img_metadata['img_width']
                # fish_counts_df['img_height_px'] = fish_img_metadata['img_height']
                
                registered_fish_df = pd.concat([registered_fish_df,fish_counts_df],axis=0,ignore_index=True)
        status = 'SUCCESS'

    return registered_fish_df, status


def register_fish_on_nuclei(processing_files:List,analysis_parameters:Dict,
                        registered_reference_channel_df,all_rounds_shifts:Dict,file_tags:Dict,status:str,
                        save=True):
    """
    The only difference to the other function is the channel name definition
    """
    logger = selected_logger()
    registration_errors = Registration_errors()
    data_models = Output_models()
    fov = file_tags['fov']
    channel = 'all_channels'
    file_tags['channel'] = channel
    registered_fish_df = data_models.output_registration_df

    if status == 'FAILED':
        error = registered_reference_channel_df['min_number_matching_dots_registration'].values[0]
        round_num = registered_reference_channel_df['round_num'].values[0]
        registered_fish_df = registered_fish_df.append({'min_number_matching_dots_registration':error,
                                                           'fov_num':int(fov),'dot_channel':channel,'round_num':int(round_num) },ignore_index=True)
    elif status == 'SUCCESS':
        reference_hybridization = registered_reference_channel_df.attrs[fov]['reference_hyb']
        reference_hybridization_str = 'Hybridization' + str(reference_hybridization).zfill(2)
        for fpath in processing_files:
            round_num = int((fpath.stem).split('_')[-5].split('Hybridization')[-1])
            try:
                fish_counts, fish_img_metadata = pickle.load(open(fpath,'rb'))
            except:
                logger.error(f'cannot open the processing files')
                
                registered_fish_df = registered_fish_df.append({'min_number_matching_dots_registration':registration_errors.cannot_load_file_fish_channel,
                                                'fov_num':int(fov),'dot_channel':channel,'round_num':round_num },ignore_index=True)
                status = 'FAILED'
                break
            else:                
                if np.any(np.isnan(fish_counts['r_px_original'])):
                    logger.error(f'There are no dots in there reference hyb for fov {fov} ')
                    registered_fish_df = registered_fish_df.append({'min_number_matching_dots_registration':registration_errors.missing_counts_fish_channel,
                                                'fov_num':int(fov),'dot_channel':channel,'round_num':round_num },ignore_index=True)
                    status = 'FAILED'
                    break

                else:
                    if reference_hybridization_str in fpath.as_posix():
                        fish_img_metadata['reference_hyb'] = reference_hybridization
                        registered_fish_df.attrs[fov] = fish_img_metadata
                    fish_counts_df = pd.DataFrame(fish_counts)
                    
                    subset_df = registered_reference_channel_df.loc[registered_reference_channel_df['round_num'] == round_num, :]
                    subset_df = subset_df.reset_index()
                        
                    shift = all_rounds_shifts[round_num]
                    fish_counts_df['r_px_registered'] = fish_counts['r_px_original'] + shift[0]
                    fish_counts_df['c_px_registered'] = fish_counts['c_px_original'] + shift[1]
                    fish_counts_df['r_shift_px'] = shift[0]
                    fish_counts_df['c_shift_px'] = shift[1]
                    fish_counts_df['min_number_matching_dots_registration'] = subset_df.loc[0,'min_number_matching_dots_registration'] 
                    registered_fish_df = pd.concat([registered_fish_df,fish_counts_df],axis=0,ignore_index=True)
                    status = 'SUCCESS'

    if save:
        fname = file_tags['experiment_fpath'] / 'tmp' / 'registered_counts' / (file_tags['experiment_name'] + '_' + file_tags['channel'] +'_registered_fov_' + file_tags['fov'] + '.parquet')
        registered_fish_df.to_parquet(fname,index=False)
    return registered_fish_df, file_tags, status

def create_registration_grps(experiment_fpath:str,registration_channel:str, fovs:List, save=True):
    """
    Function to create groups of files that need to be registered together
    Args:
        experiment_fpath: str
            Path to the exp directory 
        registration_channel : str
            Name of the channel used for registration of the fovs
        fovs: List
            List of the keys used for grouping
    Returns:
        all_grps: List
            List of tuples containing the file names of reference
            channel and fish grouped by fov
    """
    logger = selected_logger()
    tmp_fpath = Path(experiment_fpath) / 'tmp/raw_counts'
    all_files = set(tmp_fpath.glob('*'))
    registration_files = list(tmp_fpath.glob('*' + registration_channel + '*'))
    fish_files = all_files.difference(registration_files)

    all_grps = []
    for fov in fovs:
        search_key_reg = registration_channel + '_fov_' + str(fov) + '_dots.pkl'
        search_key_fish =  '_fov_' + str(fov) + '_dots.pkl'
        grp_reg = [reg_f for reg_f in registration_files if search_key_reg in reg_f.as_posix()]
        grp_fish = [fish_f for fish_f in fish_files if search_key_fish in fish_f.as_posix()]
        logger.debug(f'fov {fov} for registration on {registration_channel} has {len(grp_reg)} counts files')
        all_grps.append((grp_reg, grp_fish))
    if save:
        fname = Path(experiment_fpath) / 'tmp' / 'registration_groups.pkl'
        pickle.dump(all_grps, open(fname,'wb'))
    return all_grps



































#         counts_df = data[0]




# def registration_fish_hybridization(reference_coords:np.ndarray,shift:np.ndarray,
#                                     fov_num:int, hybridization_num_translated:int):
#     """
#     Function for the registration of the fish counts using
#     the shift calculated by aligning the registration channel

#     Parameters:
#     -----------

#     all_rounds_shifts: dict
#         dictionary containing the round number as key and 
#         shift as item
#     counts: pandas dataframe
#         pandas dataframe containing the fish counting 
#         output

#     """
#     logger = prefect.utilities.logging.get_logger("registration_fish")

#     if len(reference_coords):
#         if np.any(np.isnan(shift)):
#             registered_coords = reference_coords
#             logger.info(f'registration of {fov_num} of hybridization {hybridization_num_translated} failed')
#         else:
#             registered_coords = reference_coords + shift
#     else:
#         logger.info(f'no counts in fov {fov_num} of hybridization {hybridization_num_translated}')
#         registered_coords=np.array([np.nan,np.nan])
    
#     return registered_coords, shift



# @task(name='create-combinations-registration-fov')
# def hybridizations_registration_grps(experiment_info):

#     logger = prefect_logging_setup("create-combination-registration-fov")

#     experiment_name = experiment_info['EXP_number']
#     reference_hybridization_number = 1
#     reference_fov = 355
   
#     dots_ws, images_properties_ws, experiment_properties_ws, analysis_parameters_ws = connect_to_shoji_smfish_experiment(experiment_name)

#     stitching_channel = experiment_properties_ws[:].StitchingChannel
#     fields_of_view = images_properties_ws.FieldsOfView[(images_properties_ws.FovNumber == reference_fov) &
#                                                 (images_properties_ws.AcquistionChannel == stitching_channel) &
#                                                 (images_properties_ws.HybridizationNumber == reference_hybridization_number)]


#     fields_of_view = fields_of_view.reshape(fields_of_view.shape[1],)
    
#     fields_of_view = [355,73,122]
    
#     number_hybridizations = np.unique( images_properties_ws.HybridizationNumber[(images_properties_ws.AcquistionChannel == stitching_channel) & 
#                                                                             (images_properties_ws.FovNumber == reference_fov)])
    
#     # I will keep registration also for the reference hyb with itself to avoid extra function to rewrite the coords of
#     # unprocessed hyb
#     # number_hybridizations = number_hybridizations[number_hybridizations > reference_hybridization_number]
    
#     processing_combinations = list(itertools.product([stitching_channel],fields_of_view,[reference_hybridization_number],number_hybridizations))
    
#     return processing_combinations


# @task(name='calculate-shift-fov')
# def calculate_shift_hybridization_fov(processing_combination,experiment_info):
#     logger = prefect_logging_setup("calculate-shift-fov")
#     # Add try exect for determine if the ws are there
#     experiment_name = experiment_info['EXP_number']
    
#     dots_ws, images_properties_ws, experiment_properties_ws, analysis_parameters_ws = connect_to_shoji_smfish_experiment(experiment_name)
    
#     img_shape = images_properties_ws.ImageShape[(images_properties_ws.AcquistionChannel == processing_combination[0]) & 
#                                             (images_properties_ws.FovNumber == processing_combination[1]) & 
#                                             (images_properties_ws.HybridizationNumber == processing_combination[2])]
#     img_shape = img_shape.reshape(2,)
    

#     ref_dots_coords_fov = dots_ws.DotCoordsFOV[(dots_ws.FovNumber == processing_combination[1]) &
#                                 (dots_ws.DotChannel == processing_combination[0]) &
#                                 (dots_ws.HybridizationNumber == processing_combination[2])]
#     ref_dots_coords_fov = ref_dots_coords_fov.reshape(ref_dots_coords_fov.shape[0],2)
    

#     tran_dots_coords_fov = dots_ws.DotCoordsFOV[(dots_ws.FovNumber == processing_combination[1]) &
#                                     (dots_ws.DotChannel == processing_combination[0]) &
#                                     (dots_ws.HybridizationNumber == processing_combination[3])]
#     tran_dots_coords_fov = tran_dots_coords_fov.reshape(tran_dots_coords_fov.shape[0],2)
    
#     tran_registered_coords, shift, error, fov_num, hybridization_num_translated = \
#                             fft_registration_beads(ref_dots_coords_fov, tran_dots_coords_fov,
#                         img_shape=img_shape, fov_num=processing_combination[1],
#                         hybridization_num_translated=processing_combination[3])
    
    
#     # Write out the registration shift on shoji
#     shift = shift.reshape(1,1,1,2)
#     stitching_channel = experiment_properties_ws[:].StitchingChannel
#     images_properties_ws.RegistrationShift[(images_properties_ws.FovNumber == processing_combination[1]) & 
#                 (images_properties_ws.HybridizationNumber == processing_combination[3]) &
#                 (images_properties_ws.AcquistionChannel == stitching_channel)] = shift
    
#     # Translate the coords of all the dots for all channels and save data in shoji
#     dots_ws.DotsCoordsRegisteredFOV[(dots_ws.FovNumber == processing_combination[1]) &
#             (dots_ws.HybridizationNumber == processing_combination[3])] =  dots_ws.DotCoordsFOV[(dots_ws.FovNumber == processing_combination[1]) &
#             (dots_ws.HybridizationNumber == processing_combination[3])] + shift


#     # Calculate error shift and write it out on shoji
#     if (ref_dots_coords_fov.shape[0] >0) and (tran_registered_coords.shape[0] >0):
#         pxl = analysis_parameters_ws.BarcodesExtractionResolution[:]
#         number_matching_dots = identify_matching_register_dots_NN(ref_dots_coords_fov,tran_registered_coords,pxl)
#         number_matching_dots = np.array(number_matching_dots, dtype=np.float64).reshape(1,1,1)
#         all_errors = images_properties_ws.RegistrationError[(images_properties_ws.FovNumber == processing_combination[1]) & 
#                 (images_properties_ws.HybridizationNumber == processing_combination[3])] 
#         number_matching_dots =np.repeat(number_matching_dots,all_errors.shape[0]).reshape(all_errors.shape[0],1,1)
#         images_properties_ws.RegistrationError[(images_properties_ws.FovNumber == processing_combination[1]) & 
#                 (images_properties_ws.HybridizationNumber == processing_combination[3])] = number_matching_dots
#     else:
#         all_errors = images_properties_ws.RegistrationError[(images_properties_ws.FovNumber == processing_combination[1]) & 
#                 (images_properties_ws.HybridizationNumber == processing_combination[3])] 
#         number_matching_dots =np.repeat(0,all_errors.shape[0]).reshape(all_errors.shape[0],1,1)
#         images_properties_ws.RegistrationError[(images_properties_ws.FovNumber == processing_combination[1]) & 
#                 (images_properties_ws.HybridizationNumber == processing_combination[3])] = number_matching_dots




# @task(name='registration-fish-round')
# def registration_fish_hybridization(reference_coords:np.ndarray,shift:np.ndarray,
#                                     fov_num:int, hybridization_num_translated:int):
#     """
#     Function for the registration of the fish counts using
#     the shift calculated by aligning the registration channel

#     Parameters:
#     -----------

#     all_rounds_shifts: dict
#         dictionary containing the round number as key and 
#         shift as item
#     counts: pandas dataframe
#         pandas dataframe containing the fish counting 
#         output

#     """
#     logger = prefect.utilities.logging.get_logger("registration_fish")

#     if len(reference_coords):
#         if np.any(np.isnan(shift)):
#             registered_coords = reference_coords
#             logger.info(f'registration of {fov_num} of hybridization {hybridization_num_translated} failed')
#         else:
#             registered_coords = reference_coords + shift
#     else:
#         logger.info(f'no counts in fov {fov_num} of hybridization {hybridization_num_translated}')
#         registered_coords=np.array([np.nan,np.nan])
    
#     return registered_coords, shift























############
# utility function to use in stitching
def determine_overlap_region(self):
        """Determine the overlap between two neighbouring tiles

        Parameters:
        -----------

        ind1: int
            Index (flattened) of tile 1
        ind2: int
            Index (flattened) of tile 2

        micData: object
            MicroscopeData object containing coordinates

        Returns:
        --------

        overlap1: np.array
            Overlapping part of tile_1
        overlap2: np.array
            Overlapping part of tile_2
        plot_order: np.array
            Numpy array of ones. The shape of this array is
            used for plotting the overlaps in well fitting
            subplots.
        """
        
        
        if np.ma.is_masked(self.micData.tile_set.flat[:][self.ind1]):
            tile_1 = False
        else:
            tile_1 = True
            fnum_tile_1 = self.micData.tile_set.flat[:][self.ind1] + self.micData.tile_nr.min()
            
        
        if np.ma.is_masked(self.micData.tile_set.flat[:][self.ind2]):
            tile_2 = False
        else:
            tile_2 = True
            fnum_tile_2 = self.micData.tile_set.flat[:][self.ind2] + self.micData.tile_nr.min()
            

        if (tile_1 and tile_2):
            self.tiles_num = (fnum_tile_1, fnum_tile_2)
            self.tile1_x_coords = self.micData.x_coords[self.micData.tile_set.flat[:][self.ind1]]
            self.tile2_x_coords = self.micData.x_coords[self.micData.tile_set.flat[:][self.ind2]]
            self.tile1_y_coords = self.micData.y_coords[self.micData.tile_set.flat[:][self.ind1]]
            self.tile2_y_coords = self.micData.y_coords[self.micData.tile_set.flat[:][self.ind2]]

            tile_1_fpath = [fpath for fpath in self.counting_files_list if 'pos_'+str(fnum_tile_1)+'.' in str(fpath)][0]
            self.tile_1_store = zarr.DirectoryStore(tile_1_fpath)
            self.tile_1_root = zarr.group(store=self.tile_1_store, overwrite=False)
            self.tile_1_counts = self.tile_1_root['stringency_raw_counts']['coords_original'][...]
            tile_1_ref_coords = np.array([self.tile1_y_coords,self.tile1_x_coords])
            tile_1_ref_coords = tile_1_ref_coords[:,np.newaxis]
            self.tile_1_adj_coords = self.tile_1_counts + tile_1_ref_coords


            tile_2_fpath = [fpath for fpath in self.counting_files_list if 'pos_'+str(fnum_tile_2)+'.' in str(fpath)][0]
            self.tile_2_store = zarr.DirectoryStore(tile_2_fpath)
            self.tile_2_root = zarr.group(store=self.tile_2_store, overwrite=False)
            self.tile_2_counts = self.tile_2_root['stringency_raw_counts']['coords_original'][...]
            tile_2_ref_coords = np.array([self.tile2_y_coords,self.tile2_x_coords])
            tile_2_ref_coords = tile_2_ref_coords[:,np.newaxis]
            self.tile_2_adj_coords = self.tile_2_counts + tile_2_ref_coords


            
            if self.tile1_y_coords > self.tile2_y_coords:
                r_tl = self.tile1_y_coords
                r_br = self.tile2_y_coords + self.img_size

                r_bl = self.tile2_y_coords + self.img_size
                r_tr = self.tile1_y_coords
                
            else:
                r_tl = self.tile2_y_coords
                r_br = self.tile1_y_coords + self.img_size
                
                r_bl = self.tile1_y_coords + self.img_size
                r_tr = self.tile2_y_coords

            if self.tile1_x_coords > self.tile2_x_coords:
                c_tl = self.tile1_x_coords
                c_br = self.tile2_x_coords + self.img_size
                
                c_tr = self.tile2_x_coords + self.img_size
                c_bl = self.tile1_x_coords
                
            else:
                c_tl = self.tile2_x_coords
                c_br = self.tile1_x_coords + self.img_size
                
                c_bl = self.tile2_x_coords
                c_tr = self.tile1_x_coords + self.img_size


            self.tl_coords = np.array([r_tl,c_tl])
            self.br_coords = np.array([r_br,c_br])
            self.tr_coords = np.array([r_tr,c_tr])
            self.bl_coords = np.array([r_bl,c_bl])
            num_r = np.abs(self.tl_coords[1] - self.tr_coords[1])
            num_c = np.abs(self.tr_coords[0] - self.bl_coords[0])
            self.region_dimensions = (num_c,num_r)


        else:
            self.tl_coords = self.br_coords = self.tr_coords = self.bl_coords = None
            
###################################




class chunking():
    """
    utility class to create processing chunks
    """

    def __init__(self, region_dimensions, chunk_size, percent_padding, tl_coords):
        self.region_dimensions = region_dimensions
        self.chunk_size = chunk_size
        self.percent_padding = percent_padding
        self.tl_coords = tl_coords
    

    @staticmethod
    def block_chunks_calculator(dimension,chunk_size):
        """
        Helper function to calculate the size of the chunks created according
        the length of the vector and the chunk size.

        Parameters:
        -----------

        dimension: int
            Length of the vector to Chunk
        chunkSize: int 
            Dimension of the Chunks

        Returns:
        -----------

        chunks_sizes: np.array 
            Array of the sizes of the created chunks. It deals with conditions 
            when the expected chunks size do not fit an even number of times in the 
            dimension
        """
        number_even_chunks=int(dimension//chunk_size)
        total_size_even_chunks=number_even_chunks*chunk_size
        odd_tile_size=dimension-total_size_even_chunks
        chunk_sizes=[]
        chunks_sizes=list(np.repeat(chunk_size,number_even_chunks-1))
        if odd_tile_size < chunk_size:
            chunks_sizes.append(chunk_size+odd_tile_size)
        else:
            chunks_sizes.append(odd_tile_size)
        return tuple(chunks_sizes)
    
    def block_chunking(self):
        """
        Function used to generate the coords of the images according to the
        chunking 

        Parameters:
        -----------

        PercentPadding: float 
            Percent of overlapping between the different images (Ex. 0.2).
        ChunkSize: int 
            Dimension of the Chunks.

        Returns:
        -----------

        Coords_Chunks_list: list 
            List of np.array with the coords of the images without padding
        Coords_Padded_Chunks_list: list 
            List of np.array with the coords of the images with padding

        Notes:
        ------

        For both lists each np.array contains the coords in the following order:
        [row_tl,row_br,col_tl,col_br]

        """
        num_r,num_c = self.region_dimensions
        pixel_padding = int(self.chunk_size*self.percent_padding)
        self.starting_position = self.tl_coords

        # Calculate the size of the chunks
        r_chunks_size = self.block_chunks_calculator(num_r,self.chunk_size)
        
        c_chunks_size = self.block_chunks_calculator(num_c,self.chunk_size)
        
        # Calculate the total numbers of chunks
        nr_chunks = len(r_chunks_size)
        
        nc_chunks = len(c_chunks_size)
       


        # Coords top left corner (tl)
        if nr_chunks == 1:
            r_coords_tl = self.starting_position[0]
        else:  
            r_coords_tl = np.arange(self.starting_position[0],(self.starting_position[0]+self.chunk_size*(nr_chunks)),self.chunk_size)
        
        
        if nc_chunks == 1:
            c_coords_tl = self.starting_position[1]
        else:
            c_coords_tl = np.arange(self.starting_position[1],(self.starting_position[1]+self.chunk_size*(nc_chunks)),self.chunk_size)

        
        # Coords of all the tl in the image
        r_coords_tl_all,c_coords_tl_all = np.meshgrid(r_coords_tl,c_coords_tl,indexing='ij')
        self.coords_all_to_test = [r_coords_tl_all,c_coords_tl_all]
        # Calculate all the br coords
        r_coords_br_all = r_coords_tl_all.copy()
        c_coords_br_all = c_coords_tl_all.copy()

        for c in np.arange(0,r_coords_tl_all.shape[1]):
            r_coords_br_all[:,c] = r_coords_br_all[:,c]+r_chunks_size

        for r in np.arange(0,r_coords_tl_all.shape[0]):
             c_coords_br_all[r,:] = c_coords_br_all[r,:]+c_chunks_size

        # Calculate the padded coords
        r_coords_tl_all_padded = r_coords_tl_all-pixel_padding
        c_coords_tl_all_padded = c_coords_tl_all-pixel_padding
        r_coords_br_all_padded = r_coords_br_all+pixel_padding
        c_coords_br_all_padded = c_coords_br_all+pixel_padding

        # Correct for coords out of the image (where tl<0,br>Img.shape)
        r_coords_tl_all_padded[r_coords_tl_all_padded<0] = r_coords_tl_all[r_coords_tl_all_padded<0]
        c_coords_tl_all_padded[c_coords_tl_all_padded<0] = c_coords_tl_all[c_coords_tl_all_padded<0]
        r_coords_br_all_padded[r_coords_br_all_padded>num_r] = r_coords_br_all[r_coords_br_all_padded>num_r]
        c_coords_br_all_padded[c_coords_br_all_padded>num_c] = c_coords_br_all[c_coords_br_all_padded>num_c]

        # The coords list are generated as:
        # row_tl,row_br,col_tl,col_br


        # Create a list for the padded coords
        self.Coords_Padded_Chunks_list = list()
        for r in np.arange(0,r_coords_tl_all_padded.shape[0]):
            for c in np.arange(0,r_coords_tl_all_padded.shape[1]):
                self.Coords_Padded_Chunks_list.append(np.array([r_coords_tl_all_padded[r][c],\
                                                           r_coords_br_all_padded[r][c],\
                                                           c_coords_tl_all_padded[r][c],\
                                                           c_coords_br_all_padded[r][c]])) 
    


class reference_beads_registration():

    def __init__(self,ref_hyb_coords, comp_hyb_coords, Coords_Padded_Chunks_list, n_neighbors,
                min_acceptable_distance, min_samples, residual_threshold, max_trials, matching_radius):

        self.ref_hyb_coords = ref_hyb_coords
        self.comp_hyb_coords = comp_hyb_coords
        self.Coords_Padded_Chunks_list = Coords_Padded_Chunks_list

        self.n_neighbors = n_neighbors
        self.min_acceptable_distance = min_acceptable_distance
        self.min_samples = min_samples
        self.residual_threshold = residual_threshold
        self.max_trials = max_trials
        self.matching_radius = matching_radius

        self.logger = logging.getLogger(__name__)

    @staticmethod
    def calculate_min_distances(selected_ref):
        ref_dist = distance.cdist(selected_ref.T, selected_ref.T)
        ref_dist = np.triu(ref_dist)
        ref_dist = ref_dist[:-1,:]
        mref = np.ma.masked_where(ref_dist==0,ref_dist)
        mref_min = mref.min(axis=1)
        return np.sort(mref_min.data)
    
    @staticmethod
    # CALCULATE ERROR OVER THE DIAGONAL
    def errors(transformed_coords, ref_good):
        diagonals = np.sqrt((ref_good[:,0]- transformed_coords[:,0])**2 + (ref_good[:,1]- transformed_coords[:,1])**2)
        diagonal_mean = np.mean(diagonals,axis=0)
        diagonal_median = np.median(diagonals,axis=0)
        err1 = np.sum((diagonals-diagonal_mean),axis=0)/len(diagonals)
        err2 = np.sum((diagonals-diagonal_mean)**2,axis=0)/len(diagonals)
        diag_std = np.std(diagonals)
        diag_sem = diag_std/len(diagonals)
        
        # delta = np.abs(ref_good - transformed_coords)
        # delta_mean = np.mean(delta,axis=0)
        # err1 = np.sum((delta-delta_mean),axis=0)/len(transformed_coords)
        # err2 = np.sum((delta-delta_mean)**2,axis=0)/len(transformed_coords)
        return err1,err2, diagonal_mean, diagonal_median, diag_std, diag_sem

    def calculate_NN_roi(self,chunk_coords):    
        r_tl = chunk_coords[0]
        r_br = chunk_coords[1]
        c_tl = chunk_coords[2]
        c_br = chunk_coords[3]

        # Select only the coords in the trimmed region
        ref_trimmed = self.ref_hyb_coords[:,((r_tl < self.ref_hyb_coords[0,:]) & (self.ref_hyb_coords[0,:]<r_br)\
                                      & (c_tl <self.ref_hyb_coords[1,:]) &(self.ref_hyb_coords[1,:]<c_br)) ]
        tran_trimmed = self.comp_hyb_coords[:,((r_tl < self.comp_hyb_coords[0,:]) & (self.comp_hyb_coords[0,:]<r_br)\
                                          & (c_tl <self.comp_hyb_coords[1,:]) &(self.comp_hyb_coords[1,:]<c_br)) ]
        
        # Add check if there are dots in the timmed region
        if ref_trimmed.size >= 2*self.n_neighbors and tran_trimmed.size >= 2*self.n_neighbors:
            nbrs_ref = NearestNeighbors(n_neighbors=self.n_neighbors, algorithm='ball_tree',radius=self.matching_radius).fit(ref_trimmed.T)
            nbrs_tr = NearestNeighbors(n_neighbors=self.n_neighbors, algorithm='ball_tree',radius=self.matching_radius).fit(tran_trimmed.T)
        else:
            nbrs_ref = nbrs_tr = np.nan
            
        return ref_trimmed, tran_trimmed, nbrs_ref, nbrs_tr
   

    def calculate_NN_overlapping_region(self, coords_overlapping_region):
        
        # Add check if there are dots in the timmed region
        if coords_overlapping_region.size >= 2*self.n_neighbors and coords_overlapping_region.size >= 2*self.n_neighbors:
            nbrs = NearestNeighbors(n_neighbors=self.n_neighbors, algorithm='ball_tree',radius=self.matching_radius).fit(coords_overlapping_region.T)
        else:
            nbrs = np.nan
            
        return nbrs

    
    def calculate_registration(self):
        
        self.c = []
        passing_num = 0
        matching_points = False
        if self.ref_hyb_coords.size and self.comp_hyb_coords.size:
            for chunk_coords in self.Coords_Padded_Chunks_list:
                ref_trimmed, tran_trimmed, nbrs_ref, nbrs_tr= self.calculate_NN_roi(chunk_coords)
                if ref_trimmed.size >= 2*self.n_neighbors and tran_trimmed.size >= 2*self.n_neighbors:
                    if nbrs_ref != np.nan and nbrs_tr != np.nan:
                        # create dots id list
                        trans_dot_id_list = np.arange(tran_trimmed.shape[1])
                        for tran_dot_id in trans_dot_id_list:
                            searching=tran_trimmed[:,tran_dot_id]
                            searching = searching[np.newaxis,:]
                            dist, idx  = nbrs_tr.kneighbors(searching,n_neighbors=self.n_neighbors)
                            selected_tran = tran_trimmed[:,idx[0]]
                            tran_dist = self.calculate_min_distances(selected_tran)
                            for id_r in np.arange(ref_trimmed.shape[1]):
                                searching_r=ref_trimmed[:,id_r]
                                searching_r = searching_r[np.newaxis,:]
                                dist_r, idx_r  = nbrs_ref.kneighbors(searching_r,n_neighbors=self.n_neighbors)
                                if idx_r[0].shape[0] == idx[0].shape[0]: 
                                    selected_ref = ref_trimmed[:,idx_r[0]]
                                    ref_dist = self.calculate_min_distances(selected_ref)
                                    if np.all(np.abs(ref_dist - tran_dist)<self.min_acceptable_distance):
                                        ref = ref_trimmed[:,idx_r[0]]
                                        ref_srt = ref[:,ref[0,:].argsort()]
                                        tran = tran_trimmed[:,idx[0]]
                                        tran_srt = tran[:,tran[0,:].argsort()]
                                        cpls = np.concatenate((ref_srt.T,tran_srt.T),axis=1)
                                        self.c.append(cpls)
                                        if passing_num == 0:
                                            matching_points = cpls
                                            passing_num += 1
                                        else:
                                            matching_points = np.concatenate((matching_points,cpls),axis=0)
            if isinstance(matching_points,np.ndarray):
                matching_points = np.unique(matching_points,axis=0)
                self.ref = matching_points[:,0:2]
                self.tran = matching_points[:,2:]
                if (self.min_samples < self.ref.shape[0]) and (self.min_samples < self.tran.shape[0]):
                    self.model, self.inliers = ransac((self.tran, self.ref), transform.SimilarityTransform, min_samples=self.min_samples,
                                    residual_threshold=self.residual_threshold, max_trials=self.max_trials)
    #                 self.model = transform.estimate_transform('Affine',self.tran, self.ref)
                    self.missing_pos = False
                else:
                    self.missing_pos = True 
            else:
                self.missing_pos = True 
        else: 
            self.missing_pos = True
        
        if not self.missing_pos:
            if self.model:
                self.tr_good = self.tran[self.inliers]
                self.ref_good = self.ref[self.inliers]
                self.transformed_coords = transform.matrix_transform(self.tr_good, self.model.params)
                self.transformed_all_coords = transform.matrix_transform(self.comp_hyb_coords.T, self.model.params)
                self.delta = self.ref_good - self.tr_good
                self.translation_mean = np.mean(self.delta,axis=0)
                self.translation_median = np.median(self.delta,axis=0)
                self.translation_diagonal = np.sqrt((self.ref_good[:,0] - self.tr_good[:,0])**2 + (self.ref_good[:,1] - self.tr_good[:,1])**2)
                self.translational_diagonal_mean = np.mean(self.translation_diagonal)
                self.translational_diagonal_median = np.median(self.translation_diagonal)
                self.translation_diagonal_std = np.std(self.translation_diagonal)
                self.translation_diagonal_sem = np.std(self.translation_diagonal)/ len(self.translation_diagonal)
                self.err1,self.err2, self.diagonal_mean, self.diagonal_median, self.diag_std, self.diag_sem = self.errors(self.transformed_coords, self.ref_good)
                self.used_points = [self.ref_good,self.tr_good]
            else:
                self.missing_pos = True


    def deploy(self):
        self.registration_data = {}
        if self.ref_hyb_coords.size >= 2*self.n_neighbors:
            self.ref_nbrs = self.calculate_NN_overlapping_region(self.ref_hyb_coords)
            if self.comp_hyb_coords.size >= 2*self.n_neighbors:
                self.comp_nbrs = self.calculate_NN_overlapping_region(self.comp_hyb_coords)
                if (self.ref_nbrs != np.nan) and (self.comp_nbrs != np.nan):
                    self.calculate_registration()
                else:
                    self.logger.info(f'no neighbors identified')
                    self.missing_pos = True
            else:
                self.logger.info(f'comp region \
                        does not contain enough dots for registration')
                self.missing_pos = True
                
                # ADJUST WHEN SAVING THE DATA
                self.registration_data = {'missing_pos':self.missing_pos}
        else:
            self.logger.info(f'reference region\
                            does not contain enough dots for registration')
            self.missing_pos = True
            
        if self.missing_pos:
            self.registration_data = {'missing_pos':self.missing_pos}
        else:
            self.registration_data = {'translation_diagonal_mean':self.translational_diagonal_mean, 
                                    'translation_diagonal_median': self.translational_diagonal_median, 
                                    'translation_diagonal_std': self.translation_diagonal_std,
                                    'translation_diagonal_sem':self.translation_diagonal_sem,
                                    'used_points': self.used_points,
                                    'model_params':self.model.params, 
                                    'err1':self.err1, 
                                    'err2':self.err2, 
                                    'diagonal_mean':self.diagonal_mean,
                                    'diagonal_median':self.translational_diagonal_median,
                                    'diag_std':self.diag_std, 
                                    'diag_sem':self.diag_sem,
                                    'transformed_coords':self.transformed_coords,
                                    'transformed_all_coords':self.transformed_all_coords,
                                    'missing_pos':self.missing_pos}



class reference_beads_registration_couple():

    """
    class used to register a couple of rounds. It is used to monitor the outcome
    because it saves a lot of information useful for troubleshooting. Once the
    parameters are well define the corresponding non test function is used in the
    data processing
    """

    def __init__(self,ref_hyb_coords, comp_hyb_coords, Coords_Padded_Chunks_list, n_neighbors,
                min_acceptable_distance, min_samples, residual_threshold, max_trials, matching_radius):

        self.ref_hyb_coords = ref_hyb_coords
        self.comp_hyb_coords = comp_hyb_coords
        self.Coords_Padded_Chunks_list = Coords_Padded_Chunks_list

        self.n_neighbors = n_neighbors
        self.min_acceptable_distance = min_acceptable_distance
        self.min_samples = min_samples
        self.residual_threshold = residual_threshold
        self.max_trials = max_trials
        self.matching_radius = matching_radius

        self.logger = logging.getLogger(__name__)

    @staticmethod
    def calculate_min_distances(selected_ref):
        ref_dist = distance.cdist(selected_ref, selected_ref)
        ref_dist = np.triu(ref_dist)
        ref_dist = ref_dist[:-1,:]
        mref = np.ma.masked_where(ref_dist==0,ref_dist)
        mref_min = mref.min(axis=1)
        return np.sort(mref_min.data)
    
    @staticmethod
    # CALCULATE ERROR OVER THE DIAGONAL
    def errors(transformed_coords, ref_good):
        diagonals = np.sqrt((ref_good[:,0]- transformed_coords[:,0])**2 + (ref_good[:,1]- transformed_coords[:,1])**2)
        diagonal_mean = np.mean(diagonals,axis=0)
        diagonal_median = np.median(diagonals,axis=0)
        err1 = np.sum((diagonals-diagonal_mean),axis=0)/len(diagonals)
        err2 = np.sum((diagonals-diagonal_mean)**2,axis=0)/len(diagonals)
        diag_std = np.std(diagonals)
        diag_sem = diag_std/len(diagonals)
        
        # delta = np.abs(ref_good - transformed_coords)
        # delta_mean = np.mean(delta,axis=0)
        # err1 = np.sum((delta-delta_mean),axis=0)/len(transformed_coords)
        # err2 = np.sum((delta-delta_mean)**2,axis=0)/len(transformed_coords)
        return err1,err2, diagonal_mean, diagonal_median, diag_std, diag_sem

    def calculate_NN_roi(self,chunk_coords):    
        r_tl = chunk_coords[0]
        r_br = chunk_coords[1]
        c_tl = chunk_coords[2]
        c_br = chunk_coords[3]

        # Select only the coords in the trimmed region
        ref_trimmed = self.ref_hyb_coords[((r_tl < self.ref_hyb_coords[:,0]) & (self.ref_hyb_coords[:,0]<r_br)\
                                      & (c_tl <self.ref_hyb_coords[:,1]) &(self.ref_hyb_coords[:,1]<c_br)),:]
        tran_trimmed = self.comp_hyb_coords[((r_tl < self.comp_hyb_coords[:,0]) & (self.comp_hyb_coords[:,0]<r_br)\
                                          & (c_tl <self.comp_hyb_coords[:,1]) &(self.comp_hyb_coords[:,1]<c_br)),:]
        
        # Add check if there are dots in the timmed region
        if ref_trimmed.size >= 2*self.n_neighbors and tran_trimmed.size >= 2*self.n_neighbors:
            nbrs_ref = NearestNeighbors(n_neighbors=self.n_neighbors, algorithm='ball_tree',radius=self.matching_radius).fit(ref_trimmed)
            nbrs_tr = NearestNeighbors(n_neighbors=self.n_neighbors, algorithm='ball_tree',radius=self.matching_radius).fit(tran_trimmed)
        else:
            nbrs_ref = nbrs_tr = np.nan
            
        return ref_trimmed, tran_trimmed, nbrs_ref, nbrs_tr
   

    def calculate_NN_overlapping_region(self, coords_overlapping_region):
        
        # Add check if there are dots in the timmed region
        if coords_overlapping_region.size >= 2*self.n_neighbors and coords_overlapping_region.size >= 2*self.n_neighbors:
            nbrs = NearestNeighbors(n_neighbors=self.n_neighbors, algorithm='ball_tree',radius=self.matching_radius).fit(coords_overlapping_region)
        else:
            nbrs = np.nan
            
        return nbrs

    
    def calculate_registration(self):
        self.monitor = []
        self.c = []
        passing_num = 0
        self.matching_points = False
        if self.ref_hyb_coords.size and self.comp_hyb_coords.size:
            for chunk_coords in self.Coords_Padded_Chunks_list:
                ref_trimmed, tran_trimmed, nbrs_ref, nbrs_tr= self.calculate_NN_roi(chunk_coords)
                #if ref_trimmed.shape[0] >= 2*self.n_neighbors and tran_trimmed.shape[0] >= 2*self.n_neighbors:
                if ref_trimmed.shape[0] >= 20 and tran_trimmed.shape[0] >= 20: 
                    self.monitor.append((ref_trimmed,tran_trimmed))
                    if nbrs_ref != np.nan and nbrs_tr != np.nan:
                        # create dots id list
                        trans_dot_id_list = np.arange(tran_trimmed.shape[0])
                        for tran_dot_id in trans_dot_id_list:
                            searching=tran_trimmed[tran_dot_id,:]
                            searching = searching[np.newaxis,:]
                            dist, idx  = nbrs_tr.kneighbors(searching,n_neighbors=self.n_neighbors)
                            selected_tran = tran_trimmed[idx[0],:]
                            tran_dist = self.calculate_min_distances(selected_tran)
                            for id_r in np.arange(ref_trimmed.shape[0]):
                                self.searching_r=ref_trimmed[id_r,:]
                                self.searching_r = self.searching_r[np.newaxis,:]
                                dist_r, self.idx_r  = nbrs_ref.kneighbors(self.searching_r,n_neighbors=self.n_neighbors)
                                if self.idx_r[0].shape[0] == idx[0].shape[0]: 
                                    selected_ref = ref_trimmed[self.idx_r[0],:]
                                    ref_dist = self.calculate_min_distances(selected_ref)
                                    if np.all(np.abs(ref_dist - tran_dist)<self.min_acceptable_distance):
                                        ref = ref_trimmed[self.idx_r[0],:]
                                        ref_srt = ref[ref[0,:].argsort(),:]
                                        tran = tran_trimmed[idx[0],:]
                                        tran_srt = tran[tran[0,:].argsort(),:]
                                        cpls = np.concatenate((ref_srt,tran_srt),axis=1)
                                        self.c.append(cpls)
                                        if passing_num == 0:
                                            self.matching_points = cpls
                                            passing_num += 1
                                        else:
                                            self.matching_points = np.concatenate((self.matching_points,cpls),axis=0)
            if isinstance(self.matching_points,np.ndarray):
                self.matching_points = np.unique(self.matching_points,axis=0)
                self.ref = self.matching_points[:,0:2]
                self.tran = self.matching_points[:,2:]
                if (self.min_samples < self.ref.shape[0]) and (self.min_samples < self.tran.shape[0]):
                    self.model, self.inliers = ransac((self.tran, self.ref), transform.SimilarityTransform, min_samples=self.min_samples,
                                    residual_threshold=self.residual_threshold, max_trials=self.max_trials)
    #                 self.model = transform.estimate_transform('Affine',self.tran, self.ref)
                    self.missing_pos = False
                else:
                    self.missing_pos = True 
            else:
                self.missing_pos = True 
        else: 
            self.missing_pos = True
        
        if not self.missing_pos:
            if self.model:
                self.tr_good = self.tran[self.inliers]
                self.ref_good = self.ref[self.inliers]
                self.transformed_coords = transform.matrix_transform(self.tr_good, self.model.params)
                self.transformed_all_coords = transform.matrix_transform(self.comp_hyb_coords, self.model.params)
                self.delta = self.ref_good - self.tr_good
                self.translation_mean = np.mean(self.delta,axis=0)
                self.translation_median = np.median(self.delta,axis=0)
                self.translation_diagonal = np.sqrt((self.ref_good[:,0] - self.tr_good[:,0])**2 + (self.ref_good[:,1] - self.tr_good[:,1])**2)
                self.translational_diagonal_mean = np.mean(self.translation_diagonal)
                self.translational_diagonal_median = np.median(self.translation_diagonal)
                self.translation_diagonal_std = np.std(self.translation_diagonal)
                self.translation_diagonal_sem = np.std(self.translation_diagonal)/ len(self.translation_diagonal)
                self.err1,self.err2, self.diagonal_mean, self.diagonal_median, self.diag_std, self.diag_sem = self.errors(self.transformed_coords, self.ref_good)
                self.used_points = [self.ref_good,self.tr_good]
            else:
                self.missing_pos = True


    def deploy(self):
        self.registration_data = {}
        if self.ref_hyb_coords.size >= 2*self.n_neighbors:
            self.ref_nbrs = self.calculate_NN_overlapping_region(self.ref_hyb_coords)
            if self.comp_hyb_coords.size >= 2*self.n_neighbors:
                self.comp_nbrs = self.calculate_NN_overlapping_region(self.comp_hyb_coords)
                if (self.ref_nbrs != np.nan) and (self.comp_nbrs != np.nan):
                    self.calculate_registration()
                else:
                    self.logger.info(f'no neighbors identified')
                    self.missing_pos = True
            else:
                self.logger.info(f'comp region \
                        does not contain enough dots for registration')
                self.missing_pos = True
                
                # ADJUST WHEN SAVING THE DATA
                self.registration_data = {'missing_pos':self.missing_pos}
        else:
            self.logger.info(f'reference region\
                            does not contain enough dots for registration')
            self.missing_pos = True
            
        if self.missing_pos:
            self.registration_data = {'missing_pos':self.missing_pos}
        else:
            self.registration_data = {'translation_diagonal_mean':self.translational_diagonal_mean, 
                                    'translation_diagonal_median': self.translational_diagonal_median, 
                                    'translation_diagonal_std': self.translation_diagonal_std,
                                    'translation_diagonal_sem':self.translation_diagonal_sem,
                                    'used_points': self.used_points,
                                    'model_params':self.model.params, 
                                    'err1':self.err1, 
                                    'err2':self.err2, 
                                    'diagonal_mean':self.diagonal_mean,
                                    'diagonal_median':self.translational_diagonal_median,
                                    'diag_std':self.diag_std, 
                                    'diag_sem':self.diag_sem,
                                    'transformed_coords':self.transformed_coords,
                                    'transformed_all_coords':self.transformed_all_coords,
                                    'missing_pos':self.missing_pos}





