from typing import *
import logging
import shutil
import itertools
import math
import sys
import operator
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from pathlib import Path
from sklearn.neighbors import NearestNeighbors
from scipy.optimize import minimize

from pynndescent import NNDescent

from pysmFISH.utils import load_pipeline_config_file
from pysmFISH.logger_utils import selected_logger


class organize_square_tiles():

    """
    Class designed to determine the tile organization and identify the coords of the
    overlapping regions between the tiles

    Parameters:
    -----------

    experiment_fpath: str
        location of the experiment to process
    stitching_channel: str
        channel to use for the stitching step
    roi_num: int
        region of interest to process

    """

    def __init__(self, experiment_fpath:str, metadata:Dict,round_num:int):
        """
        round_num = int
            reference channel
        """
        
        self.logger = selected_logger()
        self.experiment_fpath = Path(experiment_fpath)
        self.metadata = metadata
        self.round_num = round_num
        
        self.experiment_name = self.metadata['experiment_name']
        self.stitching_channel = self.metadata['stitching_channel']
        self.overlapping_percentage = int(self.metadata['overlapping_percentage']) / 100
         
        self.pixel_size = self.metadata['pixel_microns']
        self.img_width = self.metadata['img_width']
        self.img_height = self.metadata['img_height']
        
        
        if  self.img_width ==  self.img_height:
            self.img_size = self.img_width
        else:
            self.logger.error(f'the images to stitch are not square')
            sys.exit(f'the images to stitch are not square')
            
    
    def extract_microscope_coords(self): 
        
        fname = self.experiment_fpath / 'microscope_tiles_coords'/(self.experiment_name + '_Hybridization'+  \
                            str(self.round_num).zfill(2) + '_' + self.stitching_channel + '_fovs_coords.npy')
        coords = np.load(fname) 
        self.x_coords = coords[:,1]
        self.y_coords = coords[:,2]
    
    def normalize_coords(self):

        if self.metadata['machine'] == 'ROBOFISH2':
            # RobofishII has stage with reference point
            # in the center (0,0)
            # consider that we get the top-right corner of the image as well

            # consider that we get the top-right corner of the image as well
            y_min = np.amin(self.y_coords)
            x_min = np.amin(self.x_coords)
            x_max = np.amax(self.x_coords)
            y_max = np.amax(self.y_coords)

            # Put the coords to zero
    
            if x_max >=0 :
                self.x_coords = self.x_coords - x_min
            else:
                self.x_coords + np.abs(x_min)
            
            if y_max>0:
                self.y_coords = self.y_coords - y_min
            else:
                self.y_coords = self.y_coords + np.abs(y_min)

            # flip y_axis
            self.y_coords = self.y_coords - self.y_coords.max()
            self.y_coords = - self.y_coords


            # change the coords from x,y to r,c
            adjusted_coords = np.zeros([self.x_coords.shape[0],2])
            adjusted_coords[:,0] = self.y_coords
            adjusted_coords[:,1] = self.x_coords

        elif self.metadata['machine'] == 'ROBOFISH1':
            # The current system has stage ref coords BOTTOM-RIGH
           
            # Normalize to (0,0) still BOTTOM-RIGHT
            y_min = np.amin(self.y_coords)
            x_min = np.amin(self.x_coords)

            self.x_coords = self.x_coords - x_min
            self.y_coords = self.y_coords - y_min

            # flip axis to move (0,0) on TOP-LEF
            self.x_coords = self.x_coords - self.x_coords.max()
            self.x_coords = - self.x_coords

            self.y_coords = self.y_coords - self.y_coords.max()
            self.y_coords = - self.y_coords

            # change the coords from x,y to r,c
            adjusted_coords = np.zeros([self.x_coords.shape[0],2])
            adjusted_coords[:,0] = self.y_coords
            adjusted_coords[:,1] = self.x_coords
        
        elif self.metadata['machine'] == 'NOT_DEFINED':
            self.logger.error(f'Need to define the specs for stitching NOT_DEFINED machine')
            sys.exit(f'Need to define the specs for stitching NOT_DEFINED machine')
            
        else:
            self.logger.error(f'define the right machine used to collected the data')
            sys.exit(f'define the right machine used to collected the data')
        
        self.tile_corners_coords_pxl = adjusted_coords / self.pixel_size
    
    
    def save_graph_original_coords(self):
        # Turn interactive plotting off
        saving_fpath = self.experiment_fpath / 'output_figures' / 'microscope_space_tiles_organization.png'
        plt.ioff()
        # Create image type axes
        labels = [str(nr) for nr in np.arange(self.x_coords.shape[0])]
        fig = plt.figure(figsize=(20,10))
        plt.plot(self.x_coords,self.y_coords,'or')

        for label, x, y in zip(labels, self.x_coords,self.y_coords):
            plt.annotate(
                label,
                xy=(x,y), xytext=(-2, 2),
                textcoords='offset points', ha='center', va='bottom',fontsize=12)
        plt.tight_layout()
        plt.savefig(saving_fpath)
    
    
    def save_graph_image_space_coords(self):
        # Turn interactive plotting off
        saving_fpath = self.experiment_fpath / 'output_figures' / 'image_space_tiles_organization.png'
        plt.ioff()
        # Create image type axes
        labels = [str(nr) for nr in np.arange(self.tile_corners_coords_pxl.shape[0])]
        fig = plt.figure(figsize=(20,10))
        plt.gca().invert_yaxis()
        plt.plot(self.tile_corners_coords_pxl[:,1],self.tile_corners_coords_pxl[:,0],'or')

        for label, x, y in zip(labels, self.tile_corners_coords_pxl[:,1],self.tile_corners_coords_pxl[:,0]):
            plt.annotate(
                label,
                xy=(x,y), xytext=(-2, 2),
                textcoords='offset points', ha='center', va='bottom',fontsize=12)
        plt.tight_layout()
        plt.savefig(saving_fpath)
        
    
    def identify_adjacent_tiles(self):
        shift_percent_tolerance = 0.05
        searching_radius = self.img_size - (self.img_size*self.overlapping_percentage) + (self.img_size*shift_percent_tolerance)
        nn = NearestNeighbors(n_neighbors=5,radius=searching_radius, metric='euclidean')
        nn.fit(self.tile_corners_coords_pxl)
        self.dists, self.indices = nn.kneighbors(self.tile_corners_coords_pxl, return_distance=True)


    def determine_overlapping_regions(self):
        # remember that overlapping region can be an empty dictionary
        self.overlapping_regions = {}
        self.overlapping_order ={}
        for idx in np.arange(self.indices.shape[0]):
            self.overlapping_regions[idx] = {}
            self.overlapping_order[idx] = {}
        for idx in np.arange(self.indices.shape[0]):
            # Determine the indices that identify the correct adjacent
            processing_indices = self.indices[idx,:]
            processing_dists = self.dists[idx,:]
            ref_tile = processing_indices[0]
            self.overlapping_regions[ref_tile] = {}
            self.overlapping_order[ref_tile] = {}
            trimmed_indices = processing_indices[1:]
            trimmed_dists = processing_dists[1:]

            idx_adj = np.where(trimmed_dists < self.img_size)
            adj_tiles_id = trimmed_indices[idx_adj]
            adj_cpls = [(ref_tile, adj_tile) for adj_tile in adj_tiles_id]
            
    
            # remove pairs that are already selected
            only_new_cpls = [cpl for cpl in adj_cpls if (cpl[1],cpl[0]) not in self.overlapping_regions[cpl[1]].keys()]
            

            if self.metadata['machine'] == 'ROBOFISH2':
                # If tile coords are top left
                for cpl in only_new_cpls:

                    tile1_r_coords = self.tile_corners_coords_pxl[cpl[0]][0]
                    tile2_r_coords = self.tile_corners_coords_pxl[cpl[1]][0]
                    tile1_c_coords = self.tile_corners_coords_pxl[cpl[0]][1]
                    tile2_c_coords = self.tile_corners_coords_pxl[cpl[1]][1]

                    if tile1_r_coords > tile2_r_coords:
                        r_tl = tile1_r_coords
                        r_br = tile2_r_coords + self.img_size

                        r_bl = tile2_c_coords + self.img_size
                        r_tr = tile1_c_coords

                        row_order = ('bottom','top')

                    else:
                        r_tl = tile2_r_coords
                        r_br = tile1_r_coords + self.img_size

                        r_bl = tile1_r_coords + self.img_size
                        r_tr = tile2_r_coords

                        row_order = ('top','bottom')

                    if tile1_c_coords > tile2_c_coords:
                        c_tl = tile1_c_coords
                        c_br = tile2_c_coords + self.img_size

                        c_tr = tile2_c_coords + self.img_size
                        c_bl = tile1_c_coords

                        col_order = ('right','left')

                    else:
                        c_tl = tile2_c_coords
                        c_br = tile1_c_coords + self.img_size

                        c_bl = tile2_c_coords
                        c_tr = tile1_c_coords + self.img_size

                        col_order = ('left','right')

                

            elif self.metadata['machine'] == 'ROBOFISH1':
                # If tile coords are bottom right
                for cpl in only_new_cpls:

                    tile1_r_coords = self.tile_corners_coords_pxl[cpl[0]][0]
                    tile2_r_coords = self.tile_corners_coords_pxl[cpl[1]][0]
                    tile1_c_coords = self.tile_corners_coords_pxl[cpl[0]][1]
                    tile2_c_coords = self.tile_corners_coords_pxl[cpl[1]][1]

                    if tile1_r_coords > tile2_r_coords:
                        r_tl = tile1_r_coords - self.img_size
                        r_br = tile2_r_coords

                        r_bl = tile2_c_coords
                        r_tr = tile1_c_coords - self.img_size

                        row_order = ('bottom','top')

                    else:
                        r_tl = tile2_r_coords - self.img_size
                        r_br = tile1_r_coords 

                        r_bl = tile1_r_coords 
                        r_tr = tile2_r_coords - self.img_size

                        row_order = ('top','bottom')

                    if tile1_c_coords > tile2_c_coords:
                        c_tl = tile1_c_coords - self.img_size
                        c_br = tile2_c_coords 

                        c_tr = tile2_c_coords
                        c_bl = tile1_c_coords - self.img_size

                        col_order = ('right','left')

                    else:
                        c_tl = tile2_c_coords - self.img_size
                        c_br = tile1_c_coords 

                        c_bl = tile2_c_coords - self.img_size
                        c_tr = tile1_c_coords 

                        col_order = ('left','right')
            else:
                pass

            self.overlapping_regions[ref_tile][cpl] = [r_tl, r_br, c_tl, c_br]
            self.overlapping_order[ref_tile][cpl] = {'row_order':row_order,'column_order':col_order}

    def run_tiles_organization(self):
        self.extract_microscope_coords()
        self.save_graph_original_coords()
        self.normalize_coords()
        self.save_graph_image_space_coords()
        self.identify_adjacent_tiles()
        self.determine_overlapping_regions()
        fname = self.experiment_fpath / 'results' / 'microscope_tile_corners_coords_pxl.npy'
        np.save(fname,self.tile_corners_coords_pxl)

def stitch_using_microscope_fov_coords(decoded_df_fpath,tile_corners_coords_pxl):
    decoded_df_fpath = Path(decoded_df_fpath)
    decoded_df = pd.read_parquet(decoded_df_fpath)
    fov = int((decoded_df_fpath.stem).split('_')[-1])
    r_microscope_coords = tile_corners_coords_pxl[fov,0]
    c_microscope_coords = tile_corners_coords_pxl[fov,1]
    if decoded_df['r_px_registered'].isnull().values.all():
        decoded_df['r_px_microscope_stitched'] = np.nan
        decoded_df['c_px_microscope_stitched'] = np.nan
    else:
        decoded_df['r_px_microscope_stitched'] =  r_microscope_coords - decoded_df['r_px_registered']
        decoded_df['c_px_microscope_stitched'] =  c_microscope_coords - decoded_df['c_px_registered']

        # decoded_df['r_px_microscope_stitched'] =  r_microscope_coords + decoded_df['r_px_registered']
        # decoded_df['c_px_microscope_stitched'] =  c_microscope_coords + decoded_df['c_px_registered']
    
    decoded_df.to_parquet(decoded_df_fpath)
    return decoded_df


def stitch_using_microscope_fov_coords_test(decoded_df,fov,tile_corners_coords_pxl):
    r_microscope_coords = tile_corners_coords_pxl[fov,0]
    c_microscope_coords = tile_corners_coords_pxl[fov,1]
    if decoded_df['r_px_registered'].isnull().values.all():
        decoded_df['r_px_microscope_stitched'] = np.nan
        decoded_df['c_px_microscope_stitched'] = np.nan
    else:
        decoded_df['r_px_microscope_stitched'] =  r_microscope_coords - decoded_df['r_px_registered']
        decoded_df['c_px_microscope_stitched'] =  c_microscope_coords - decoded_df['c_px_registered']

        # decoded_df['r_px_microscope_stitched'] =  r_microscope_coords + decoded_df['r_px_registered']
        # decoded_df['c_px_microscope_stitched'] =  c_microscope_coords + decoded_df['c_px_registered']
    return decoded_df


def stitch_using_microscope_fov_coords_new(decoded_df,tile_corners_coords_pxl):
    if decoded_df['r_px_registered'].isnull().values.all():
        decoded_df['r_px_microscope_stitched'] = np.nan
        decoded_df['c_px_microscope_stitched'] = np.nan
    else:
        fov = decoded_df.iloc[0]['fov_num']
        r_microscope_coords = tile_corners_coords_pxl[fov,0]
        c_microscope_coords = tile_corners_coords_pxl[fov,1]
        decoded_df['r_px_microscope_stitched'] =  r_microscope_coords - decoded_df['r_px_registered']
        decoded_df['c_px_microscope_stitched'] =  c_microscope_coords - decoded_df['c_px_registered']

        # decoded_df['r_px_microscope_stitched'] =  r_microscope_coords + decoded_df['r_px_registered']
        # decoded_df['c_px_microscope_stitched'] =  c_microscope_coords + decoded_df['c_px_registered']
    return decoded_df

# REMOVE OVERLAPPING DOTS ACCORDING TO GENE
# preprocessing and removal part to put in the flow file

# all_files = (Path(experiment_fpath) / 'tmp' / 'registered_counts').glob('*decoded*.parquet')
# counts_dd_list = [dd.read_parquet(counts_file) for counts_file in all_files]
# counts_dd = dd.concat(counts_dd_list, axis=0)
# counts_dd = counts_dd.loc[counts_dd.dot_id == counts_dd.barcode_reference_dot_id,['barcode_reference_dot_id',
#                                                                                     r_tag, c_tag, select_genes, 
#                                                                                     'fov_num']]
# counts_df = counts_dd.dropna(subset=[select_genes]).compute()
# grpd = counts_df.groupby(select_genes)

# all_futures = []

# for gene, count_df in grpd:
#     future = client.submit(remove_overlapping_dots_from_gene,
#                             experiment_fpath = experiment_fpath,
#                             counts_df=counts_df,
#                             unfolded_overlapping_regions_dict=corrected_overlapping_regions_dict,
#                             stitching_selected=stitching_selected,
#                             gene = gene,
#                             same_dot_radius = same_dot_radius)
    
#     all_futures.append(future)

def get_dots_in_overlapping_regions(counts_df, unfolded_overlapping_regions_dict, 
                       stitching_selected, gene):    
    r_tag = 'r_px_' + stitching_selected
    c_tag = 'c_px_' + stitching_selected
    
    ref_tiles_df = pd.DataFrame(columns=counts_df.columns)
    comp_tiles_df = pd.DataFrame(columns=counts_df.columns)
    
    
    grpd_df = counts_df.groupby('fov_num')
    list_fov = list(grpd_df.groups.keys())
    
    for cpl, chunk_coords in unfolded_overlapping_regions_dict.items():
        
        if (cpl[0] in list_fov) and (cpl[1] in list_fov):
            r_tl = chunk_coords[0]
            r_br = chunk_coords[1]
            c_tl = chunk_coords[2]
            c_br = chunk_coords[3]

            barcoded_ref_df = grpd_df.get_group(cpl[0])
            barcoded_comp_df = grpd_df.get_group(cpl[1])

            overlapping_ref_df = barcoded_ref_df.loc[(barcoded_ref_df[r_tag] > r_tl) & (barcoded_ref_df[r_tag] < r_br) 
                                               & (barcoded_ref_df[c_tag] > c_tl) & (barcoded_ref_df[c_tag] < c_br),:]


            overlapping_comp_df = barcoded_comp_df.loc[(barcoded_comp_df[r_tag] > r_tl) & (barcoded_comp_df[r_tag] < r_br) 
                                               & (barcoded_comp_df[c_tag] > c_tl) & (barcoded_comp_df[c_tag] < c_br),:]


            ref_tiles_df = ref_tiles_df.append(overlapping_ref_df)
            comp_tiles_df = comp_tiles_df.append(overlapping_comp_df)
        
    return ref_tiles_df, comp_tiles_df


def identify_duplicated_dots(ref_tiles_df,comp_tiles_df,stitching_selected,same_dot_radius):
    
    r_tag = 'r_px_' + stitching_selected
    c_tag = 'c_px_' + stitching_selected

    overlapping_ref_coords = ref_tiles_df.loc[:, [r_tag,c_tag]].to_numpy()
    overlapping_comp_coords = comp_tiles_df.loc[:, [r_tag,c_tag]].to_numpy()
    dots_ids = comp_tiles_df.loc[:, ['dot_id']].to_numpy()
    index = NNDescent(overlapping_ref_coords,metric='euclidean',n_neighbors=1)
    indices, dists = index.query(overlapping_comp_coords,k=1)
    idx_dists = np.where(dists < same_dot_radius)[0]
    dots_id_to_remove = dots_ids[idx_dists]
    dots_id_to_remove = list(dots_id_to_remove.reshape(dots_id_to_remove.shape[0],))
    return dots_id_to_remove


def remove_overlapping_dots_from_gene(experiment_fpath,counts_df,unfolded_overlapping_regions_dict,
                                    stitching_selected,gene,same_dot_radius):

    experiment_fpath = Path(experiment_fpath)
    ref_tiles_df, comp_tiles_df = get_dots_in_overlapping_regions(counts_df,unfolded_overlapping_regions_dict, 
                       stitching_selected, gene)
    dots_id_to_remove = identify_duplicated_dots(ref_tiles_df,comp_tiles_df,stitching_selected,same_dot_radius)
    cleaned_df = counts_df.loc[~counts_df.barcode_reference_dot_id.isin(dots_id_to_remove), :]
    fpath = experiment_fpath / 'results' / (experiment_fpath.stem + '_' + gene +'_counts.parquet')
    cleaned_df.to_parquet(fpath,index=False)

    pass


# REMOVED OVERLAPPING DOTS ACCORDING TO FOV (MUCH FASTER THAN FOR GENE)
# EXPECIALLY FOR LARGE AREAS WITH A LOT OF COUNTS

def get_all_dots_in_overlapping_regions(counts_df, chunk_coords, 
                       stitching_selected):    
    
    r_tag = 'r_px_' + stitching_selected
    c_tag = 'c_px_' + stitching_selected
    
    subset_df = counts_df.loc[counts_df.dot_id == counts_df.barcode_reference_dot_id, :]
    
    
    r_tl = chunk_coords[0]
    r_br = chunk_coords[1]
    c_tl = chunk_coords[2]
    c_br = chunk_coords[3]
    
    overlapping_ref_df = subset_df.loc[(subset_df[r_tag] > r_tl) & (subset_df[r_tag] < r_br) 
                                               & (subset_df[c_tag] > c_tl) & (subset_df[c_tag] < c_br),:]
        
    return overlapping_ref_df

def identify_duplicated_dots_NNDescend(ref_tiles_df,comp_tiles_df,stitching_selected,same_dot_radius):
    
    r_tag = 'r_px_' + stitching_selected
    c_tag = 'c_px_' + stitching_selected

    overlapping_ref_coords = ref_tiles_df.loc[:, [r_tag,c_tag]].to_numpy()
    overlapping_comp_coords = comp_tiles_df.loc[:, [r_tag,c_tag]].to_numpy()
    dots_ids = comp_tiles_df.loc[:, ['dot_id']].to_numpy()
    index = NNDescent(overlapping_ref_coords,metric='euclidean',n_neighbors=1)
    indices, dists = index.query(overlapping_comp_coords,k=1)
    idx_dists = np.where(dists < same_dot_radius)[0]
    dots_id_to_remove = dots_ids[idx_dists]
    dots_id_to_remove = list(dots_id_to_remove.reshape(dots_id_to_remove.shape[0],))
    return dots_id_to_remove

def identify_duplicated_dots_sklearn(ref_tiles_df,comp_tiles_df, stitching_selected,same_dot_radius):
    
    nn = NearestNeighbors(n_neighbors=1,radius=same_dot_radius, metric='euclidean')
    
    r_tag = 'r_px_' + stitching_selected
    c_tag = 'c_px_' + stitching_selected
    
    overlapping_ref_coords = ref_tiles_df.loc[:, [r_tag,c_tag]].to_numpy()
    overlapping_comp_coords = comp_tiles_df.loc[:, [r_tag,c_tag]].to_numpy()
    dots_ids = comp_tiles_df.loc[:, ['dot_id']].to_numpy()
    nn.fit(overlapping_ref_coords)
    dists, indices = nn.kneighbors(overlapping_comp_coords, return_distance=True)
    idx_dists = np.where(dists < same_dot_radius)[0]
    dots_id_to_remove = dots_ids[idx_dists]
    dots_id_to_remove = list(dots_id_to_remove.reshape(dots_id_to_remove.shape[0],))
    
    return dots_id_to_remove

def remove_overlapping_dots_fov(cpl, chunk_coords, experiment_fpath,
                                    stitching_selected,
                                    select_genes,same_dot_radius):

    logger = selected_logger()
    all_dots_id_to_remove = []
    experiment_fpath = Path(experiment_fpath)
    
    try:
        counts1_fpath = list((experiment_fpath / 'tmp' / 'registered_counts').glob('*decoded*_fov_' + str(cpl[0]) + '.parquet'))[0]
    except:
        logger.error(f'count file missing for fov {cpl[0]}')
    
    else:
        try:
            counts2_fpath = list((experiment_fpath / 'tmp' / 'registered_counts').glob('*decoded*_fov_' + str(cpl[1]) + '.parquet'))[0]
        except:
            logger.error(f'count file missing for fov {cpl[1]}')
        else:
    
            counts1_df = pd.read_parquet(counts1_fpath)
            counts2_df = pd.read_parquet(counts2_fpath)
            
            overlap_count1 = get_all_dots_in_overlapping_regions(counts1_df, chunk_coords, 
                            stitching_selected)
            overlap_count2 = get_all_dots_in_overlapping_regions(counts2_df, chunk_coords, 
                            stitching_selected)
            
            count1_grp = overlap_count1.groupby(select_genes)
            count2_grp = overlap_count2.groupby(select_genes)
            
            for gene, over_c1_df in count1_grp:
                try:
                    over_c2_df = count2_grp.get_group(gene)
                except:
                    pass
                else:
                    dots_id_to_remove = identify_duplicated_dots_sklearn(over_c1_df,over_c2_df,
                                                                stitching_selected,same_dot_radius)
                    if len(dots_id_to_remove):
                        all_dots_id_to_remove.append(dots_id_to_remove)
            all_dots_id_to_remove = [el for tg in all_dots_id_to_remove for el in tg]
            return all_dots_id_to_remove

def clean_from_duplicated_dots(fov, dots_id_to_remove,experiment_fpath):
    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)
    try:
        fname = list((experiment_fpath / 'tmp' / 'registered_counts').glob('*_decoded_fov_' + str(fov) + '.parquet'))[0]
    except:
        logger.error(f'missing decoded file for fov {fov}')
    else:
        save_name = fname.stem.split('_decoded_fov_')[0] + '_cleaned_df_fov_' + str(fov) + '.parquet'
        save_name = experiment_fpath / 'results' / save_name
        if len(dots_id_to_remove):
            try:
                counts_df = pd.read_parquet(fname)
                logger.error(f'loaded {fname}')
                
            except:
                logger.error(f'missing {fname}')
            else:
                cleaned_df = counts_df.loc[~counts_df.barcode_reference_dot_id.isin(dots_id_to_remove), :]
                cleaned_df.to_parquet(save_name,index=False)
                logger.error(f'saved {fname}')
        else:
            try:
                _ = shutil.copy2(fname.as_posix(),save_name.as_posix())
                logger.error(f'copied {fname}')
            except:
                logger.error(f'cannot copy {fname} to {save_name}')



class r_c_chunking():
    """
    Utility class used to chunk and arbitrary region and obtain the coords if the chunks.
    The chunking can be different between row and 
    columns

    Parameters:
    -----------

    region_dimensions: np.ndarray
        number of rows and columns of the region to chunk
    r_chunk_size: float
        size of the chunks along the rows
    c_chunk_size: float
        size of the chunks along the columns
    tl_coords: np.ndarray
        coordinate of the top left corner of the region to chunk
        to use to calculate the coords of the chunks

    """

    def __init__(self, region_dimensions, r_chunk_size, c_chunk_size, tl_coords):
        self.region_dimensions = region_dimensions
        self.r_chunk_size = r_chunk_size
        self.c_chunk_size = c_chunk_size
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

        Notes:
        ------

        For both lists each np.array contains the coords in the following order:
        [row_tl,row_br,col_tl,col_br]

        """
        num_r,num_c = self.region_dimensions
        self.starting_position = self.tl_coords
        self.end_position = self.tl_coords + self.region_dimensions

        # Calculate the size of the chunks
        r_chunks_size = self.block_chunks_calculator(num_r,self.r_chunk_size)
        
        c_chunks_size = self.block_chunks_calculator(num_c,self.c_chunk_size)
        
        # Calculate the total numbers of chunks
        nr_chunks = len(r_chunks_size)
        
        nc_chunks = len(c_chunks_size)
       


        # Coords top left corner (tl)
        if nr_chunks == 1:
            r_coords_tl = self.starting_position[0]
        else: 
            r_coords_tl = [self.starting_position[0]]
            for i in np.arange(1,nr_chunks):
                r_coords_tl.append(r_coords_tl[i-1] + self.r_chunk_size )
            r_coords_tl = np.array(r_coords_tl)
            # r_coords_tl = np.arange(self.starting_position[0],(self.starting_position[0]+self.r_chunk_size*(nr_chunks)),self.r_chunk_size)
            
        if nc_chunks == 1:
            c_coords_tl = self.starting_position[1]
        else:
            c_coords_tl = [self.starting_position[1]]
            for i in np.arange(1,nc_chunks):
                c_coords_tl.append(c_coords_tl[i-1] + self.c_chunk_size )
            c_coords_tl = np.array(c_coords_tl)

            # c_coords_tl = np.arange(self.starting_position[1],(self.starting_position[1]+self.c_chunk_size*(nc_chunks)),self.c_chunk_size)

        
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

        
        # The coords list are generated as:
        # row_tl,row_br,col_tl,col_br


        # Create a list for the padded coords
        self.coords_chunks_list = list()
        for r in np.arange(0,r_coords_tl_all.shape[0]):
            for c in np.arange(0,r_coords_tl_all.shape[1]):
                self.coords_chunks_list.append(np.array([r_coords_tl_all[r][c],\
                                                           r_coords_br_all[r][c],\
                                                           c_coords_tl_all[r][c],\
                                                           c_coords_br_all[r][c]])) 
    


class triangles_based_dots_stitching():
    """
    Class used to register the different rounds by searaching and
    matching all possible triangles formed by the dots in the reference
    and translated image. This function run only a registration to the reference
    round
    
    The calculation of the triangle is based on list processing and may 
    be improved in ported to numpy.
    https://stackoverflow.com/questions/43126580/match-set-of-x-y-points-to-another-set-that-is-scaled-rotated-translated-and

    """

    
    def __init__(self, ref_overlapping_counts, comp_overlapping_counts, chunk_coords):
        self.ref_overlapping_counts = ref_overlapping_counts
        self.comp_overlapping_counts = comp_overlapping_counts
        self.chunk_coords = chunk_coords
        
        self.r_tl = self.chunk_coords[0]
        self.r_br = self.chunk_coords[1]
        self.c_tl = self.chunk_coords[2]
        self.c_br = self.chunk_coords[3]
        
        
        num_r = np.abs(np.abs(self.r_tl) - np.abs(self.r_br))
        num_c = np.abs(np.abs(self.c_tl) - np.abs(self.c_br))
        self.overlapping_region_dimensions = np.array([num_r,num_c])
        
        if num_r > num_c:
            self.chunk_search_ax = 'r'
            self.r_chunk_size = num_c
            self.c_chunk_size = num_c
            self.max_chunk_size = num_r
        else:
            self.chunk_search_ax = 'c'
            self.r_chunk_size = num_r
            self.c_chunk_size = num_r
            self.max_chunk_size = num_c
        
        self.min_dots_chunk = 6
        self.min_error_triangles = 1
        self.max_dots = 12
        self.logger = logging.getLogger(__name__)
        self.tl_coords = np.array([self.r_tl, self.c_tl])
          

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
        self.chunks = r_c_chunking(self.overlapping_region_dimensions,self.r_chunk_size,
                        self.c_chunk_size,self.tl_coords)   
        self.chunks.block_chunking()
        self.coords_chunks_list = self.chunks.coords_chunks_list
        
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
        if self.chunk_search_ax == 'c':
            chunk_size = self.c_chunk_size
        else:
            chunk_size = self.r_chunk_size
        
        while chunk_size < self.max_chunk_size:
            chunks = r_c_chunking(self.overlapping_region_dimensions,self.r_chunk_size,
                        self.c_chunk_size,self.tl_coords)
            chunks.block_chunking()
            coords_chunks_list = chunks.coords_chunks_list
            ref_max_number_dots = []
            tran_max_number_dots = []
            ref_total = []
            tran_total = []
            for chunk_coords in coords_chunks_list:
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
                if self.chunk_search_ax == 'c':
                    self.c_chunk_size += 200
                else:
                    self.r_chunk_size += 200
                                  
        if self.enough_dots:
            # Collect the ref and tran coords from the chunks with enough dots
            self.ref_tran_screening_list = []
            for chunk_coords in coords_chunks_list:
                ref_coords_in_chunk = self.calculate_dots_chunks(ref_coords,chunk_coords)
                tran_coords_in_chunk = self.calculate_dots_chunks(tran_coords,chunk_coords)
                if ref_coords_in_chunk.shape[0] > self.min_dots_chunk and tran_coords_in_chunk.shape[0] > self.min_dots_chunk:
                    self.ref_tran_screening_list.append((ref_coords_in_chunk,tran_coords_in_chunk,chunk_coords))

    def register(self,ref_coords,tran_coords):
        self.optimize_chunking(ref_coords, tran_coords)
        self.completed_registration = False
        if self.enough_dots:
            match_ref_pts_all = []
            match_tran_pts_all = []
            # Collect all matching dots in all chunked regions with number of dots above threshold
            for ref_coords_in_chunk,tran_coords_in_chunk, chunk_coords in self.ref_tran_screening_list:
                    match_ref_pts, match_tran_pts = self.identify_matching_coords(ref_coords_in_chunk,tran_coords_in_chunk,self.min_error_triangles)
                    if len(match_ref_pts) and len(match_tran_pts):
                        match_ref_pts_all.append(match_ref_pts)
                        match_tran_pts_all.append(match_tran_pts)
                    if len(match_ref_pts_all) > self.max_dots:
                        break
            match_ref_pts_all = [pts for grp in match_ref_pts_all for pts in grp]
            match_tran_pts_all = [pts for grp in match_tran_pts_all for pts in grp]

            if len(match_ref_pts_all):
                match_ref_pts_all = np.vstack(match_ref_pts_all)
                match_tran_pts_all = np.vstack(match_tran_pts_all)
                minimization_output = minimize(self.obj_fun,[0,0],args=(match_ref_pts_all,match_tran_pts_all), method='Nelder-Mead')
                if minimization_output.success:
                    self.tran_registered_coords = self.apply_transform(minimization_output.x, tran_coords)
                    self.transformation_matrix = minimization_output.x
                    self.completed_registration = True
                else:
                    self.logger.info(f'chunk {chunk_coords} failed minimization of distances')
            else:
                self.logger.info(f'chunk {chunk_coords} did not find matching triangles')

        else:
            self.logger.info(f'cannot register rounds not enough dots')
            self.tran_registered_coords = tran_coords
            self.transformation_matrix = np.empty([1,2])
            self.transformation_matrix[:] = np.nan

        if not self.completed_registration:
            self.logger.info(f'was not possible to register ')
            self.tran_registered_coords = tran_coords
            self.transformation_matrix = np.empty([1,2])
            self.transformation_matrix[:] = np.nan



def register_adj_tiles(experiment_fpath, roi_num, stitching_channel, idx_reference_tile,overlapping_regions,tile_corners_coords_pxl):
    stitching_shift = {}
    # tmp = []
    experiment_fpath = Path(experiment_fpath)
    counts_fpath = experiment_fpath / 'counts'/ ('roi_' + str(roi_num)) / stitching_channel
    search_key = '*_fov_' + str(idx_reference_tile) + '.parquet'
    # Adde crosscheck for error
    ref_counts_fpath = list(counts_fpath.glob(search_key))[0]
    ref_counts_df = pd.read_parquet(ref_counts_fpath)

    ref_tile_coords = tile_corners_coords_pxl[idx_reference_tile]
    ref_counts_selected = ref_counts_df.loc[ref_counts_df.round_num == 1, ['r_px_registered', 'c_px_registered']]
    ref_counts_selected.loc[:, ['r_px_registered', 'c_px_registered']] += ref_tile_coords

    for cpl, chunk_coords in overlapping_regions[idx_reference_tile].items():
        # Add crosscheck for error
        search_key = '*_fov_' + str(cpl[1]) + '.parquet'
        comp_counts_fpath = list(counts_fpath.glob(search_key))[0]
        comp_counts_df = pd.read_parquet(comp_counts_fpath)
        comp_counts_selected = comp_counts_df.loc[comp_counts_df.round_num == 1, ['r_px_registered', 'c_px_registered']]
        comp_tile_coords = tile_corners_coords_pxl[cpl[1]]
        comp_counts_selected.loc[:, ['r_px_registered', 'c_px_registered']] += comp_tile_coords

        r_tl = chunk_coords[0]
        r_br = chunk_coords[1]
        c_tl = chunk_coords[2]
        c_br = chunk_coords[3]

        overlapping_ref_coords_df = ref_counts_selected.loc[(ref_counts_selected.r_px_registered > r_tl) & (ref_counts_selected.r_px_registered < r_br) 
                                                   & (ref_counts_selected.c_px_registered > c_tl) & (ref_counts_selected.c_px_registered < c_br),['r_px_registered','c_px_registered']]

        overlapping_comp_coords_df = comp_counts_selected.loc[(comp_counts_selected.r_px_registered > r_tl) & (comp_counts_selected.r_px_registered < r_br) 
                                                   & (comp_counts_selected.c_px_registered > c_tl) & (comp_counts_selected.c_px_registered < c_br),['r_px_registered','c_px_registered']]

        overlapping_ref_coords = overlapping_ref_coords_df.to_numpy()
        overlapping_comp_coords = overlapping_comp_coords_df.to_numpy()

        tr_st = triangles_based_dots_stitching(overlapping_ref_coords, overlapping_comp_coords, chunk_coords)
        tr_st.optimize_chunking(overlapping_ref_coords,overlapping_comp_coords)
        tr_st.register(overlapping_ref_coords,overlapping_comp_coords)
        stitching_shift[cpl] = tr_st.transformation_matrix
        # tmp.append([overlapping_ref_coords, overlapping_comp_coords, tr_st.tran_registered_coords])
        
    return stitching_shift



if __name__ ==  '__main__':

    exp = '/Users/simone/Documents/local_data_storage/test_micdata/LBEXP20200325_oPool11'
    stitching_channel = 'Europium'
    roi_num = 0
    to = organize_square_tiles(exp,stitching_channel,roi_num)
    # to.save_graph()
    to.normalize_coords()
    to.identify_adjacent_tiles()
    to.determine_overlapping_regions()
    idx_reference_tile = 0
    experiment_fpath = '/Users/simone/Documents/local_data_storage/dots_analysis/LBEXP20200325_oPool11/'
    stitching_shift, tmp = register_adj_tiles(experiment_fpath, roi_num, stitching_channel, idx_reference_tile,to.overlapping_regions,to.tile_corners_coords_pxl)
    to.save_graph()