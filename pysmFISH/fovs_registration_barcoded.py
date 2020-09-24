
"""
group of class or functions use to register the fov
"""
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

from pysmFISH.utils import load_pipeline_config_file
from skimage.feature import register_translation
# from skimage.registration import phase_cross_correlation UPDATE SKIMAGEN
from skimage import filters

import itertools
import math
import operator
from scipy.optimize import minimize

import prefect
from prefect import task
from prefect.engine import signals

from pysmFISH.logger_utils import prefect_logging_setup



def create_fake_image(img_shape,coords):
    gaussian_sigma = 5
    img = np.zeros(img_shape,dtype=np.float64)
    img[coords[:,0].astype(int),coords[:,1].astype(int)] = 1000
    img = filters.gaussian(img,sigma=gaussian_sigma)
    return img


@task(name='fft_registration_beads')# ADD LOGGER
def fft_registration_beads(reference_coords:np.ndarray, translated_coords:np.ndarray,
                       img_width: int, img_height:int, fov_num:int,
                       hybridization_num_translated:int):

    logger = prefect.utilities.logging.get_logger("fft_registration_beads")
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
        img_shape = [img_width,img_height]
        img_ref = create_fake_image(img_shape,reference_coords)    
        img_tran = create_fake_image(img_shape,translated_coords)
        shift, error, diffphase = register_translation(img_ref, img_tran)
        tran_registered_coords = translated_coords + shift
        tran_registered_coords = tran_registered_coords.astype(int)

    return tran_registered_coords, shift, error, fov_num, hybridization_num_translated 


@task(name='registration-fish-round')
def registration_fish_hybridization(reference_coords:np.ndarray,shift:np.ndarray,
                                    fov_num:int, hybridization_num_translated:int):
    """
    Function for the registration of the fish counts using
    the shift calculated by aligning the registration channel

    Parameters:
    -----------

    all_rounds_shifts: dict
        dictionary containing the round number as key and 
        shift as item
    counts: pandas dataframe
        pandas dataframe containing the fish counting 
        output

    """
    logger = prefect.utilities.logging.get_logger("registration_fish")

    if len(reference_coords):
        if np.any(np.isnan(shift)):
            registered_coords = reference_coords
            logger.info(f'registration of {fov_num} of hybridization {hybridization_num_translated} failed')
        else:
            registered_coords = reference_coords + shift
    else:
        logger.info(f'no counts in fov {fov_num} of hybridization {hybridization_num_translated}')
        registered_coords=np.array([np.nan,np.nan])
    
    return registered_coords, shift

























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





class triangles_based_registration():
    """
    Class used to register the different rounds by searaching and
    matching all possible triangles formed by the dots in the reference
    and translated image. This function run only a registration to the reference
    round
    
    The calculation of the triangle is based on list processing and may 
    be improved in ported to numpy.
    https://stackoverflow.com/questions/43126580/match-set-of-x-y-points-to-another-set-that-is-scaled-rotated-translated-and

    """

    def __init__(self, counts, experiment_fpath, channel_name, roi_number, fov_name):
        self.counts = counts
        self.experiment_fpath = Path(experiment_fpath)
        self.channel_name = channel_name
        self.roi_number = roi_number
        self.fov_name = fov_name
        
        self.logger = logging.getLogger(__name__)
        
        self.pipeline_config_fpath = self.experiment_fpath / 'pipeline_config'
        self.experiment_config_fpath =  self.pipeline_config_fpath / 'experiment.yaml'
        self.experiment_config = load_pipeline_config_file(self.experiment_config_fpath)
        
        searching_key = '*roi_' +str(self.roi_number) + '_' + self.channel_name + '_images_config.yaml'        
        
        try:
            self.image_config_fpath = list(self.pipeline_config_fpath.glob(searching_key))[0]
        except:
            self.logger.error(f'the reference beads image_config file is missing {searching_key}')
            sys.exit(f'the reference beads image_config file')

        
        # Load registration parameters
        self.image_config = load_pipeline_config_file(self.image_config_fpath)
        self.registration_parameters = self.image_config[fov_name]['fov_analysis_parameters']['rounds_registration']
        self.chunk_size = self.registration_parameters['chunk_size']
        self.min_dots_chunk = self.registration_parameters['min_dots_chunk']
        self.min_error_triangles = self.registration_parameters['min_error_triangles']
        self.percent_padding = self.registration_parameters['percent_padding']
        self.reference_round = self.registration_parameters['reference_round']
        self.collect_all_chunks = self.registration_parameters['collect_all_chunks']
        self.reference_round_name = 'round_' + str(self.reference_round)
       
        # The top lef coords are 0,0 because we are using relative coords
        self.tl_coords = (0,0)

        self.img_dimensions = (self.image_config[fov_name]['rounds'][self.reference_round_name]['shape']['height'],
                               self.image_config[fov_name]['rounds'][self.reference_round_name]['shape']['width'])    
                    

    @staticmethod
    def combine_coords(counts, round_num):
        data_reference = counts.loc[counts['round_num'] == round_num]
        r_px = data_reference.r_px_original.to_list()
        c_px = data_reference.c_px_original.to_list()
        coords = np.array(list(zip(r_px,c_px)))
        position_idx = data_reference.index
        return coords, position_idx


    @staticmethod
    def obj_fun(pars,x,src):
        tx, ty = pars
        H = np.array([[1, 0, tx],\
            [0, 1, ty]])
        src1 = np.c_[src,np.ones(src.shape[0])]
        return np.sum( (x - src1.dot(H.T)[:,:2])**2 )

    @staticmethod
    def apply_transform(pars, src):
        tx, ty = pars
        H = np.array([[1, 0, tx],\
            [0, 1, ty]])
        src1 = np.c_[src,np.ones(src.shape[0])]
        return src1.dot(H.T)[:,:2]

    @staticmethod
    def distance(x1,y1,x2,y2):
        return math.sqrt((x2 - x1)**2 + (y2 - y1)**2 )

    @staticmethod
    def list_subtract(list1,list2):
        return np.absolute(np.array(list1)-np.array(list2))

    def tri_sides(self,set_x, set_x_tri):

        triangles = []
        for i in range(len(set_x_tri)):

            point1 = set_x_tri[i][0]
            point2 = set_x_tri[i][1]
            point3 = set_x_tri[i][2]

            point1x, point1y = set_x[point1][0], set_x[point1][1]
            point2x, point2y = set_x[point2][0], set_x[point2][1]
            point3x, point3y = set_x[point3][0], set_x[point3][1] 

            len1 = self.distance(point1x,point1y,point2x,point2y)
            len2 = self.distance(point1x,point1y,point3x,point3y)
            len3 = self.distance(point2x,point2y,point3x,point3y)

            # you need to normalize in case the ref and the tran
            # are warped
            #min_side = min(len1,len2,len3)
            #len1/=min_side
            #len2/=min_side
            #len3/=min_side
            t=[len1,len2,len3]
            t.sort()
            triangles.append(t)

        return triangles


    def identify_matching_coords(self,set_A, set_B, threshold):
        match_A_pts = []
        match_B_pts = []
        set_A_tri = list(itertools.combinations(range(len(set_A)), 3))
        set_B_tri = list(itertools.combinations(range(len(set_B)), 3))
        A_triangles = self.tri_sides(set_A, set_A_tri)
        B_triangles = self.tri_sides(set_B, set_B_tri)
        sums = []
        for i in range(len(A_triangles)):
            for j in range(len(B_triangles)):
                k = sum(self.list_subtract(A_triangles[i], B_triangles[j]))
                if k < threshold:
                    sums.append([i,j,k])
        # sort by smallest sum
        sums = sorted(sums, key=operator.itemgetter(2))
        if len(sums):
            match_A = set_A_tri[sums[0][0]]
            match_B = set_B_tri[sums[0][1]]
            for i in range(3):
                match_A_pts.append(set_A[match_A[i]])
                match_B_pts.append(set_B[match_B[i]])
        return (match_A_pts,match_B_pts)



    def calculate_chunks(self):
        self.chunks = chunking(self.img_dimensions,self.chunk_size,
                          self.registration_parameters['percent_padding'],self.tl_coords)   
        self.chunks.block_chunking()
        self.Coords_Padded_Chunks_list = self.chunks.Coords_Padded_Chunks_list
        
    def calculate_dots_chunks(self,coords,chunk_coords):  
        r_tl = chunk_coords[0]
        r_br = chunk_coords[1]
        c_tl = chunk_coords[2]
        c_br = chunk_coords[3]

        # Select only the coords in the trimmed region
        coords_in_chunk = coords[((r_tl < coords[:,0]) & (coords[:,0]<r_br)\
                    & (c_tl <coords[:,1]) &(coords[:,1]<c_br)),: ]
        return coords_in_chunk


    def optimize_chunking(self,ref_coords, tran_coords):       
        self.enough_dots = False
        chunk_size = self.chunk_size
        while chunk_size < min(self.img_dimensions):
            chunks = chunking(self.img_dimensions, chunk_size, self.percent_padding, self.tl_coords)
            chunks.block_chunking()
            Coords_Padded_Chunks_list = chunks.Coords_Padded_Chunks_list
            ref_max_number_dots = []
            tran_max_number_dots = []
            ref_total = []
            tran_total = []
            for chunk_coords in Coords_Padded_Chunks_list:
                ref_coords_in_chunk = self.calculate_dots_chunks(ref_coords,chunk_coords)
                tran_coords_in_chunk = self.calculate_dots_chunks(tran_coords,chunk_coords)
                if ref_coords_in_chunk.shape[0] > self.min_dots_chunk and tran_coords_in_chunk.shape[0] > self.min_dots_chunk:
                        self.enough_dots = True
                        break
            if self.enough_dots:
                break
            else:
                self.enough_dots = False
                chunk_size += 200

        if self.enough_dots:
            # Collect the ref and tran coords from the chunks with enough dots
            self.ref_tran_screening_list = []
            for chunk_coords in Coords_Padded_Chunks_list:
                ref_coords_in_chunk = self.calculate_dots_chunks(ref_coords,chunk_coords)
                tran_coords_in_chunk = self.calculate_dots_chunks(tran_coords,chunk_coords)
                if ref_coords_in_chunk.shape[0] > self.min_dots_chunk and tran_coords_in_chunk.shape[0] > self.min_dots_chunk:
                    self.ref_tran_screening_list.append((ref_coords_in_chunk,tran_coords_in_chunk,chunk_coords))
            self.chunk_size = chunk_size


    def rounds_to_register(self):
        # remember that round_num is not pythonic
        self.rounds_list = list(self.counts['round_num'].unique())
        if self.reference_round == False:
            logger.error(f'missing reference round')
            sys.exit(f'missing reference round')
        elif self.reference_round <0 or self.reference_round > max(self.rounds_list):
            logger.error(f'selected reference round {self.reference_round} is out of range')
            sys.exit(f'selected reference round {self.reference_round} is out of range')
        else:
            self.rounds_list = list(self.counts['round_num'].unique()) 
            self.processing_combinations = []
            self.rounds_list.remove(self.reference_round)
            for el in self.rounds_list:
                self.processing_combinations.append((self.reference_round,el))

    def cpl_registration(self,ref_coords,tran_coords,cpl):
        self.optimize_chunking(ref_coords, tran_coords)
        self.completed_registration = False
        if self.enough_dots:
            match_ref_pts_all = []
            match_tran_pts_all = []
            if self.collect_all_chunks:
            # Collect all matching dots in all chunked regions with number of dots above threshold
                for ref_coords_in_chunk,tran_coords_in_chunk, chunk_coords in self.ref_tran_screening_list:
                        match_ref_pts, match_tran_pts = self.identify_matching_coords(ref_coords_in_chunk,tran_coords_in_chunk,self.min_error_triangles)
                        if len(match_ref_pts) and len(match_tran_pts):
                            match_ref_pts_all.append(match_ref_pts)
                            match_tran_pts_all.append(match_tran_pts)
                match_ref_pts_all = [pts for grp in match_ref_pts_all for pts in grp]
                match_tran_pts_all = [pts for grp in match_tran_pts_all for pts in grp]
            else:
                ref_coords_in_chunk,tran_coords_in_chunk, chunk_coords = self.ref_tran_screening_list[0]
                match_ref_pts_all, match_tran_pts_all = self.identify_matching_coords(ref_coords_in_chunk,tran_coords_in_chunk,self.min_error_triangles)

            if len(match_ref_pts_all):
                match_ref_pts_all = np.vstack(match_ref_pts_all)
                match_tran_pts_all = np.vstack(match_tran_pts_all)
                minimization_output = minimize(self.obj_fun,[0,0],args=(match_ref_pts_all,match_tran_pts_all), method='Nelder-Mead')
                if minimization_output.success:
                    self.tran_registered_coords = self.apply_transform(minimization_output.x, tran_coords)
                    self.transformation_matrix = minimization_output.x
                    self.completed_registration = True
                else:
                    self.logger.info(f'chunk {chunk_coords} of {cpl} failed minimization of distances')
            else:
                self.logger.info(f'chunk {chunk_coords} of {cpl} did not find matching triangles')

        else:
            self.logger.info(f'cannot register rounds {cpl} not enough dots')
            self.tran_registered_coords = tran_coords
            self.transformation_matrix = np.empty([1,2])
            self.transformation_matrix[:] = np.nan

        if not self.completed_registration:
            self.logger.info(f'was not possible to register {cpl} ')
            self.tran_registered_coords = tran_coords
            self.transformation_matrix = np.empty([1,2])
            self.transformation_matrix[:] = np.nan


    def register(self):
        all_registration_matrix = {}
        r_px_col_name = 'r_px_registered'
        c_px_col_name = 'c_px_registered'
        self.counts[r_px_col_name] = np.nan
        self.counts[c_px_col_name] = np.nan
        self.counts['r_transformation_registration'] = np.nan
        self.counts['c_transformation_registration'] = np.nan
        self.rounds_to_register()
        for cpl in self.processing_combinations:
            ref_round, tran_round = cpl
            ref_coords, ref_position_idx = self.combine_coords(self.counts,ref_round)
            tran_coords, tran_position_idx = self.combine_coords(self.counts,tran_round)
            self.cpl_registration(ref_coords,tran_coords,cpl)
            all_registration_matrix[cpl] = self.transformation_matrix
            self.ref_registered_coords = ref_coords
            self.counts.loc[ref_position_idx, r_px_col_name] = self.ref_registered_coords[:,0]
            self.counts.loc[ref_position_idx, c_px_col_name] = self.ref_registered_coords[:,1]
            self.counts.loc[ref_position_idx, 'r_transformation_registration'] = np.ones(self.ref_registered_coords.shape[0])
            self.counts.loc[ref_position_idx, 'c_transformation_registration'] = np.ones(self.ref_registered_coords.shape[0])
            self.counts.loc[tran_position_idx, r_px_col_name] = self.tran_registered_coords[:,0]
            self.counts.loc[tran_position_idx, c_px_col_name] = self.tran_registered_coords[:,1]
            all_transf = np.tile(self.transformation_matrix,(self.tran_registered_coords.shape[0],1))
            self.counts.loc[tran_position_idx, 'r_transformation_registration'] = all_transf[:,0]
            self.counts.loc[tran_position_idx, 'c_transformation_registration'] = all_transf[:,1]
       
        return self.counts, all_registration_matrix



