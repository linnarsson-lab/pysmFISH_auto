"""
Function and classes used to identify barcodes
"""
from typing import *
import pandas as pd
import numpy as np
import pickle
import logging
from sklearn.neighbors import NearestNeighbors
# from pynndescent import NNDescent
from pathlib import Path
from itertools import groupby

from pysmFISH.logger_utils import selected_logger
from pysmFISH.data_models import Output_models
from pysmFISH.errors import Registration_errors


class simplify_barcodes_reference():
    """Utility Class use to convert excels files with codebook info
    in smaller size pandas dataframe/parquet files to pass to dask
    workers during the processing. This utility function must be
    run before running the experiment analysis. The pipeline
    require the output of this function.
    """
    
    def __init__(self, barcode_fpath: str):
        """Class initialization

        Args:
            barcode_fpath (str): Path to the xlsx file with the codebook
        """
        
        self.barcode_fpath = Path(barcode_fpath)
        self.barcode_fname = self.barcode_fpath.stem

    @staticmethod
    def format_codeword(codeword: str):
        """[summary]

        Args:
            codeword (str): codeword representing a gene

        Returns:
            byte: codeword converted in byte representation
        """
        str_num = codeword.split('[')[-1].split(']')[0]
        converted_codeword = np.array([int(el) for el in list(str_num)]).astype(np.int8)
        converted_codeword = converted_codeword.tobytes()
        return converted_codeword

    def convert_codebook(self):
        used_gene_codebook_df = pd.read_excel(self.barcode_fpath)
        # used_gene_codebook_df = pd.read_parquet(self.barcode_fpath)
        self.codebook_df = used_gene_codebook_df.loc[:,['Barcode','Gene']]
        self.codebook_df.rename(columns = {'Barcode':'Code'}, inplace = True)
        self.codebook_df.Code = self.codebook_df.Code.apply(lambda x: self.format_codeword(x))
        self.codebook_df.to_parquet(self.barcode_fpath.parent / (self.barcode_fname + '.parquet'))



def dots_hoods(coords: np.ndarray,pxl: int)->np.ndarray:
    """Function that calculate the coords of the peaks searching
    neighborhood for identifying the barcodes.

    Args:
        coords (np.ndarray): coords of the identified peaks
        pxl (int): size of the neighborhood in pixel

    Returns:
        np.ndarray: coords that define the neighborhood (r_tl,r_br,c_tl,c_tr)
    """
    r_tl = coords[:,0]-pxl
    r_br = coords[:,0]+pxl
    c_tl = coords[:,1]-pxl
    c_tr = coords[:,1]+pxl
    r_tl = r_tl[:,np.newaxis]
    r_br = r_br[:,np.newaxis]
    c_tl = c_tl[:,np.newaxis]
    c_tr = c_tr[:,np.newaxis]
    chunks_coords = np.hstack((r_tl,r_br,c_tl,c_tr))
    chunks_coords = chunks_coords.astype(int)
    return chunks_coords


def extract_dots_images(barcoded_df: pd.DataFrame,registered_img_stack: np.ndarray,
                experiment_fpath: str, metadata: dict):
    """Function used to extract the images corresponding to a barcode
    after running the decoding identification. It can save the images
    but to avoid increasing too much the space occupied by a processed
    experiment an array with the maximum intensity value of the pxl in
    each round is calculated and saved 

    Args:
        barcoded_df (pd.DataFrame): Dataframe with decoded barcodes 
                    for a specific field of view.
        registered_img_stack (np.ndarray): Preprocessed image of a single field of view
                    the imaging round correspond to the z-stack position
        experiment_fpath (str): Path to the folder of the experiment to process
        metadata (dict): Overall experiment info
    """
    if isinstance(registered_img_stack, np.ndarray) and (barcoded_df.shape[0] >1):
        experiment_fpath = Path(experiment_fpath)
          
        round_intensity_labels = ['bit_' + str(el) +'_intensity' for el in np.arange(1,metadata['total_rounds']+1)]

        barcodes_names = barcoded_df['barcode_reference_dot_id'].values
        coords = barcoded_df.loc[:, ['r_px_registered', 'c_px_registered']].to_numpy()
        barcodes_extraction_resolution = barcoded_df['barcodes_extraction_resolution'].values[0]

        chunks_coords = dots_hoods(coords,barcodes_extraction_resolution)
        chunks_coords[chunks_coords<0]=0
        chunks_coords[chunks_coords>registered_img_stack.shape[1]]= registered_img_stack.shape[1]
    
        for idx in np.arange(chunks_coords.shape[0]):
            selected_region = registered_img_stack[:,chunks_coords[idx,0]:chunks_coords[idx,1]+1,chunks_coords[idx,2]:chunks_coords[idx,3]+1]
            if selected_region.size >0:
                max_array = selected_region.max(axis=(1,2))
                barcoded_df.loc[barcoded_df.dot_id == barcodes_names[idx],round_intensity_labels] = max_array


        # for channel in channels:
            # all_regions[channel] = {}
            # all_max[channel] = {}
            # img_stack = registered_img_stack[channel]
            # trimmed_df_channel = trimmed_df.loc[trimmed_df.channel == channel]
            # if trimmed_df_channel.shape[0] >0:

            #     barcodes_names = trimmed_df_channel['barcode_reference_dot_id'].values
            #     coords = trimmed_df_channel.loc[:, ['r_px_registered', 'c_px_registered']].to_numpy()
            #     barcodes_extraction_resolution = trimmed_df_channel['barcodes_extraction_resolution'].values[0]

            #     chunks_coords = dots_hoods(coords,barcodes_extraction_resolution)
            #     chunks_coords[chunks_coords<0]=0
            #     chunks_coords[chunks_coords>img_stack.shape[1]]= img_stack.shape[1]
            
            
            #     for idx in np.arange(chunks_coords.shape[0]):
            #         selected_region = img_stack[:,chunks_coords[idx,0]:chunks_coords[idx,1]+1,chunks_coords[idx,2]:chunks_coords[idx,3]+1]
            #         if selected_region.size >0:
            #             max_array = selected_region.max(axis=(1,2))
            #             # all_regions[channel][barcodes_names[idx]]= selected_region
            #             all_max[channel][barcodes_names[idx]]= max_array
            #             barcoded_df.loc[barcoded_df.dot_id == barcodes_names[idx],round_intensity_labels] = max_array

        # fpath = experiment_fpath / 'tmp' / 'combined_rounds_images' / (experiment_name + '_' + channel + '_img_dict_fov_' + str(fov) + '.pkl')
        # pickle.dump(all_regions,open(fpath,'wb'))
        # fpath = experiment_fpath / 'results' / (experiment_name + '_barcodes_max_array_dict_fov_' + str(fov) + '.pkl')
        # pickle.dump(all_max,open(fpath,'wb'))
    return barcoded_df









def identify_flipped_bits(codebook: pd.DataFrame, gene: str, 
                    raw_barcode: ByteString)-> Tuple[ByteString, ByteString]:
    """Utility function used to identify the position of the bits that are
    flipped after the nearest neighbors and the definition of the
    acceptable hamming distance for a single dot.

    Args:
        codebook (pd.DataFrame): Codebook used for the decoding
        gene (str): Name of the gene identified
        raw_barcode (ByteString): identifide barcode from the images

    Returns:
        Tuple[ByteString, ByteString]: (flipped_position, flipping_direction)
    """
    gene_barcode_str =codebook.loc[codebook.Gene == gene, 'Code'].values[0]
    gene_barcode = np.frombuffer(gene_barcode_str, np.int8)
    raw_barcode = np.frombuffer(raw_barcode, np.int8)
    flipped_positions = np.where(raw_barcode != gene_barcode)[0].astype(np.int8)
    flipping_directions = (gene_barcode[flipped_positions] - raw_barcode[flipped_positions]).astype(np.int8)
    # flipped_positions = flipped_positions.tobytes()
    # flipping_directions = flipping_directions.tobytes()
    return flipped_positions,flipping_directions


def define_flip_direction(codebook_dict: dict,experiment_fpath: str, 
            output_df: pd.DataFrame):
    """Function used to determinethe the position of the bits that are
    flipped after the nearest neighbors and the definition of the
    acceptable hamming distance for fov.

    Args:
        codebook (dict): Codebooks used for the decoding
        experiment_fpath (str): Path to the folder of the experiment to process
        output_df (pd.DataFrame): Dataframe with the decoded results for 
                    the specific fov.
    """
    if output_df.shape[0] > 1:
        correct_hamming_distance = 0
        selected_hamming_distance = 3 / output_df.iloc[0].barcode_length
        experiment_fpath = Path(experiment_fpath)
        experiment_name = experiment_fpath.stem
        channels = codebook_dict.keys()
        all_evaluated = []
        for channel in channels: 
            codebook = codebook_dict[channel]
            fov = output_df.fov_num.values[0]
            trimmed_df = output_df.loc[(output_df.dot_id == output_df.barcode_reference_dot_id) &
                                (output_df.channel == channel) &
                                (output_df['hamming_distance'] > correct_hamming_distance) &
                                (output_df['hamming_distance'] < selected_hamming_distance),
                                    ['barcode_reference_dot_id', 'decoded_genes', 'raw_barcodes','hamming_distance']]
            trimmed_df = trimmed_df.dropna(subset=['decoded_genes'])
            trimmed_df.loc[:,('flip_and_direction')] = trimmed_df.apply(lambda x: identify_flipped_bits(codebook,x.decoded_genes,x.raw_barcodes),axis=1)
            trimmed_df['flip_position'] = trimmed_df['flip_and_direction'].apply(lambda x: x[0])
            trimmed_df['flip_direction'] = trimmed_df['flip_and_direction'].apply(lambda x: x[1])
            trimmed_df.drop(columns=['flip_and_direction'],inplace=True)
            all_evaluated.append(trimmed_df)
        
        all_evaluated = pd.concat(all_evaluated,axis=0,ignore_index=True,inplace=True)

        
        fpath = experiment_fpath / 'results' / (experiment_name + '_' + channel + '_df_flip_direction_fov' + str(fov) + '.parquet')
        all_evaluated.to_parquet(fpath)
        # return trimmed_df


def chunk_dfs(dataframes_list: list, chunk_size: int):
    """ 
    Functions modified from
    https://stackoverflow.com/questions/45217120/how-to-efficiently-join-merge-concatenate-large-data-frame-in-pandas
    yields n dataframes at a time where n == chunksize 
    """
    dfs = []
    for f in dataframes_list:
        dfs.append(f)
        if len(dfs) == chunk_size:
            yield dfs
            dfs  = []
    if dfs:
        yield dfs


def merge_with_concat(dfs: list)->pd.DataFrame:
    """Utility function used to merge dataframes

    Args:
        dsf (list): List with the dataframe to merge

    Returns:
        pd.DataFrame: Merged dataframe
    """                                           

#     dfs = (df.set_index(col, drop=True) for df in dfs)
    merged = pd.concat(dfs, axis=0, join='outer', copy=False)
    return merged


"""
    Class used to extract the barcodes from the registered
    counts using nearest neighbour

    Parameters:
    -----------
    counts: pandas.DataFrame
        pandas file with the fov counts after
        registration
    analysis_parameters: dict
        parameters for data processing 
    codebook_df: pandas.DataFrame
        pandas file with the codebook used to
        deconvolve the barcode

    NB: if there is a problem with the registration the barcode assigned 
        will be 0*barcode_length
    
    """


def extract_barcodes_NN_fast_multicolor(registered_counts_df: pd.DataFrame, analysis_parameters: Dict,
                codebook_df: pd.DataFrame, metadata:dict)-> Tuple[pd.DataFrame,pd.DataFrame]:
    """Function used to extract the barcodes from the registered
    counts using nearest neighbour. if there is a problem with the registration the barcode assigned 
    will be 0*barcode_length

    Args:
        registered_counts_df (pd.Dataframe): Fov counts after registration
        analysis_parameters (Dict): Parameters for data processing 
        codebook_df (pd.DataFrame): codebook used to deconvolve the barcode
    Returns:
        Tuple[pd.DataFrame,pd.DataFrame]: (barcoded_round, all_decoded_dots_df)
    """

    logger = selected_logger()

    barcodes_extraction_resolution = analysis_parameters['BarcodesExtractionResolution']
    RegistrationMinMatchingBeads = analysis_parameters['RegistrationMinMatchingBeads']
    barcode_length = metadata['barcode_length']
    registration_errors = Registration_errors()

    stitching_channel = metadata['stitching_channel']
    
 
    registered_counts_df.dropna(subset=['dot_id'],inplace=True)
    # Starting level for selection of dots
    dropping_counts = registered_counts_df.copy(deep=True)

    all_decoded_dots_list = []
    barcoded_round = []
    if registered_counts_df['r_px_registered'].isnull().values.any():
        
        all_decoded_dots_df = pd.DataFrame(columns = registered_counts_df.columns)
        all_decoded_dots_df['decoded_genes'] = np.nan
        all_decoded_dots_df['hamming_distance'] = np.nan
        all_decoded_dots_df['number_positive_bits'] = np.nan
        all_decoded_dots_df['barcode_reference_dot_id'] = np.nan
        all_decoded_dots_df['raw_barcodes'] = np.nan
        all_decoded_dots_df['barcodes_extraction_resolution'] = barcodes_extraction_resolution
        # Save barcoded_round and all_decoded_dots_df
        return registered_counts_df, all_decoded_dots_df
        
    else:
        for ref_round_number in np.arange(1,barcode_length+1):

            #ref_round_number = 1
            reference_round_df = dropping_counts.loc[dropping_counts.round_num == ref_round_number,:]
            # Step one (all dots not in round 1)
            compare_df = dropping_counts.loc[dropping_counts.round_num!=ref_round_number,:]

            if (not reference_round_df.empty):
                if not compare_df.empty:
                    nn = NearestNeighbors(1, metric="euclidean")
                    nn.fit(reference_round_df[['r_px_registered','c_px_registered']])
                    dists, indices = nn.kneighbors(compare_df[['r_px_registered','c_px_registered']], return_distance=True)

                    # select only the nn that are below barcodes_extraction_resolution distance
                    idx_distances_below_resolution = np.where(dists <= barcodes_extraction_resolution)[0]

                    comp_idx = idx_distances_below_resolution
                    ref_idx = indices[comp_idx].flatten()

                    # Subset the dataframe according to the selected points
                    # The reference selected will have repeated points
                    comp_selected_df = compare_df.iloc[comp_idx]
                    ref_selected_df = reference_round_df.iloc[ref_idx]

                    # The size of ref_selected_df w/o duplicates may be smaller of reference_round_df if 
                    # some of the dots in reference_round_df have no neighbours

                    # Test approach where we get rid of the single dots
                    comp_selected_df.loc[:,'barcode_reference_dot_id'] = ref_selected_df['dot_id'].values
                    ref_selected_df_no_duplicates = ref_selected_df.drop_duplicates()
                    ref_selected_df_no_duplicates.loc[:,'barcode_reference_dot_id'] = ref_selected_df_no_duplicates['dot_id'].values

                    # Collect singletons
                    # Remeber that this method works only because there are no duplicates inside the dataframes
                    # https://stackoverflow.com/questions/48647534/python-pandas-find-difference-between-two-data-frames
                    if reference_round_df.shape[0] > ref_selected_df_no_duplicates.shape[0]:
                        singletons_df = pd.concat([reference_round_df,ref_selected_df_no_duplicates]).drop_duplicates(keep=False)
                        singletons_df.loc[:,'barcode_reference_dot_id'] = singletons_df['dot_id'].values
                        barcoded_round = pd.concat([comp_selected_df, ref_selected_df_no_duplicates,singletons_df], axis=0,ignore_index=False)
                    else:
                        barcoded_round = pd.concat([comp_selected_df, ref_selected_df_no_duplicates], axis=0,ignore_index=False)

                    barcoded_round = pd.concat([comp_selected_df, ref_selected_df_no_duplicates,singletons_df], axis=0,ignore_index=False)
                    barcoded_round_grouped = barcoded_round.groupby('barcode_reference_dot_id')

                    compare_df = compare_df.drop(comp_selected_df.index)
                    dropping_counts = compare_df

                else:
                    # Collecting singleton of last bit
                    reference_round_df.loc[:,'barcode_reference_dot_id'] = reference_round_df['dot_id'].values
                    barcoded_round_grouped = reference_round_df.groupby('barcode_reference_dot_id')
                    ref_selected_df_no_duplicates = reference_round_df

                for brdi, grp in barcoded_round_grouped:
                    barcode = np.zeros([barcode_length],dtype=np.int8)
                    barcode[grp.round_num.values.astype(np.int8)-1] = 1
                    #hamming_dist, index_gene = nn_sklearn.kneighbors(barcode.reshape(1, -1), return_distance=True)
                    #gene= codebook_df.loc[index_gene.reshape(index_gene.shape[0]),'Gene'].tolist()

                    barcode = barcode.tostring()
                    ref_selected_df_no_duplicates.loc[ref_selected_df_no_duplicates.barcode_reference_dot_id == brdi,'raw_barcodes'] = barcode
                    #ref_selected_df_no_duplicates.loc[ref_selected_df_no_duplicates.barcode_reference_dot_id == brdi,'decoded_gene_name'] = gene
                    #ref_selected_df_no_duplicates.loc[ref_selected_df_no_duplicates.barcode_reference_dot_id == brdi,'hamming_distance'] = hamming_dist.flatten()[0]

                    #fish_counts.loc[grp.index,'barcode_reference_dot_id'] = brdi
                    #fish_counts.loc[grp.index,'raw_barcodes'] = barcode
                    #dists, index = nn_sklearn.kneighbors(all_barcodes, return_distance=True)

                all_decoded_dots_list.append(ref_selected_df_no_duplicates)

        if all_decoded_dots_list:
            all_decoded_dots_df = pd.concat(all_decoded_dots_list,ignore_index=False)
                    
                
            codebook_df = convert_str_codebook(codebook_df,'Code')
            codebook_array = make_codebook_array(codebook_df,'Code')
            nn_sklearn = NearestNeighbors(n_neighbors=1, metric="hamming")
            nn_sklearn.fit(codebook_array)

            all_barcodes = np.vstack(all_decoded_dots_df.raw_barcodes.map(lambda x: np.frombuffer(x, np.int8)).values)
            dists_arr, index_arr = nn_sklearn.kneighbors(all_barcodes, return_distance=True)
            genes=codebook_df.loc[index_arr.reshape(index_arr.shape[0]),'Gene'].tolist()

            all_decoded_dots_df.loc[:,'decoded_genes'] = genes
            all_decoded_dots_df.loc[:,'hamming_distance'] = dists_arr
            all_decoded_dots_df.loc[:,'number_positive_bits'] = all_barcodes.sum(axis=1)

            all_decoded_dots_df['barcodes_extraction_resolution'] = barcodes_extraction_resolution
        else:
            all_decoded_dots_df = pd.DataFrame(columns = registered_counts_df.columns)
            all_decoded_dots_df['decoded_genes'] = np.nan
            all_decoded_dots_df['hamming_distance'] = np.nan
            all_decoded_dots_df['number_positive_bits'] = np.nan
            all_decoded_dots_df['barcode_reference_dot_id'] = np.nan
            all_decoded_dots_df['raw_barcodes'] = np.nan
            all_decoded_dots_df['barcodes_extraction_resolution'] = barcodes_extraction_resolution

        
        # Save barcoded_round and all_decoded_dots_df
        return barcoded_round, all_decoded_dots_df



# TODO Remove all the functions below
######## -------------------------------------------------------------------

class extract_barcodes_NN():
    """
    Class used to extract the barcodes from the registered
    counts using nearest neighbour

    Parameters:
    -----------
    counts: pandas.DataFrame
        pandas file with the fov counts after
        registration
    analysis_parameters: dict
        parameters for data processing
    experiment_config: Dict
        dictionary with the experimental data 
    codebook_df: pandas.DataFrame
        pandas file with the codebook used to
        deconvolve the barcode

    NB: if there is a problem with the registration the barcode assigned 
        will be 0*barcode_length
    
    """

    def __init__(self, counts, analysis_parameters:Dict,experiment_config:Dict,codebook_df,file_tags,status:str):
        
        self.barcodes_extraction_resolution = analysis_parameters['BarcodesExtractionResolution']
        self.RegistrationMinMatchingBeads = analysis_parameters['RegistrationMinMatchingBeads']
        self.barcode_length = experiment_config['Barcode_length']
        self.counts = counts
        self.logger = selected_logger()
        self.codebook_df = codebook_df
        self.file_tags = file_tags
        self.status = status
        self.registration_errors = Registration_errors()



        
    @staticmethod
    def barcode_nn(counts_df, ref_round_number, barcodes_extraction_resolution):
        column_names = list(counts_df.columns.values)
        column_names = column_names.append('barcode_reference_dot_id')
        barcoded_df = pd.DataFrame(columns=column_names)

        reference_array = counts_df.loc[counts_df.round_num == ref_round_number, ['r_px_registered','c_px_registered']].to_numpy()
        reference_round_df = counts_df.loc[counts_df.round_num == ref_round_number,:].reset_index(drop=True)
        # Step one (all dots not in round 1)
        coords_compare = counts_df.loc[counts_df.round_num != ref_round_number, ['r_px_registered','c_px_registered']].to_numpy()
        compare_df = counts_df.loc[counts_df.round_num != ref_round_number,:].reset_index(drop=True)

        if (reference_array.shape[0] >0) and (coords_compare.shape[0] >0):
            # initialize network
            nn = NearestNeighbors(1, metric="euclidean")
            nn.fit(reference_array)

            # Get the nn
            dists, indices = nn.kneighbors(coords_compare, return_distance=True)

            # select only the nn that are below barcodes_extraction_resolution distance
            idx_selected_coords_compare = np.where(dists <= barcodes_extraction_resolution)[0]

            compare_selected_df = compare_df.loc[idx_selected_coords_compare,:]
            compare_selected_df['barcode_reference_dot_id'] = np.nan

            # ref_idx = indices[idx_selected_coords_compare]
            # compare_selected_df.loc[compare_selected_df.index.isin(idx_selected_coords_compare),'barcode_reference_dot_id'] = reference_round_df.loc[ref_idx,'dot_id'].values[0]
            
            for idx in idx_selected_coords_compare:
                ref_idx = indices[idx]
                compare_selected_df.loc[idx,'barcode_reference_dot_id'] = reference_round_df.loc[ref_idx,'dot_id'].values[0]

            reference_round_df['barcode_reference_dot_id'] = reference_round_df.dot_id

            barcoded_df = barcoded_df.append([compare_selected_df, reference_round_df], ignore_index=True)

            compare_df = compare_df.drop(compare_selected_df.index)
            compare_df = compare_df.reset_index(drop=True)
    
        return compare_df, barcoded_df

    @staticmethod
    def convert_str_codebook(codebook_df,column_name):
        codebook_df[column_name] = codebook_df[column_name].map(lambda x: np.frombuffer(x, np.int8))
        return codebook_df

    @staticmethod
    def make_codebook_array(codebook_df,column_name):
        codebook_array = np.zeros((len(codebook_df[column_name]),codebook_df[column_name][0].shape[0]))
        for idx, el in enumerate(codebook_df[column_name]):
            row = codebook_df[column_name][idx]
            row = row[np.newaxis,:]
            codebook_array[idx,:] = row
        return codebook_array


    def run_extraction(self):

        data_models = Output_models()
        registration_errors = Registration_errors()
        fov = self.file_tags['fov']
        channel = self.file_tags['channel']
        self.barcoded_fov_df = data_models.barcode_analysis_df
        self.barcoded_fov_df.attrs = self.counts.attrs

        if self.status == 'FAILED':
            error = self.counts['min_number_matching_dots_registration'].values[0]
            round_num = self.counts['round_num'].values[0]
            self.barcoded_fov_df = self.barcoded_fov_df.append({'min_number_matching_dots_registration':error,
                                                           'fov_num':int(fov),'dot_channel':channel,'round_num': round_num },ignore_index=True)
        elif self.status == 'SUCCESS':

            if (min(self.counts.loc[:,'min_number_matching_dots_registration']) < self.RegistrationMinMatchingBeads):
                round_num = self.counts['round_num'].values[0]
                self.barcoded_fov_df = self.barcoded_fov_df.append({'min_number_matching_dots_registration':registration_errors.registration_below_extraction_resolution, 
                                            'fov_num':int(fov),'dot_channel':channel,'round_num': round_num},ignore_index=True)
                self.status = 'FAILED'
            else:
                hd_2 = 2 / self.barcode_length
                hd_3 = 3 / self.barcode_length
                # barcode_length = len(self.counts['round_num'].unique())
                rounds = np.arange(1,self.barcode_length+1)
                self.codebook_df = self.convert_str_codebook(self.codebook_df,'Code')
                codebook_array = self.make_codebook_array(self.codebook_df,'Code')
                nn_sklearn = NearestNeighbors(n_neighbors=1, metric="hamming")
                nn_sklearn.fit(codebook_array)

                # remove points with np.NAN
                # self.counts = self.counts.dropna()
                for round_num in rounds:
                    compare_df, barcoded_df = self.barcode_nn(self.counts, round_num, self.barcodes_extraction_resolution)
                    self.barcoded_fov_df = self.barcoded_fov_df.append(barcoded_df, ignore_index=True)
                    self.counts = compare_df

                self.counts['barcode_reference_dot_id'] = self.counts.dot_id
                self.barcoded_fov_df = self.barcoded_fov_df.append(self.counts, ignore_index=True)
                self.barcoded_fov_df['barcodes_extraction_resolution'] = self.barcodes_extraction_resolution
                self.grpd = self.barcoded_fov_df.groupby('barcode_reference_dot_id')
                # self.all_barcodes = {}
                # for name, group in self.grpd:
                #     rounds_num = group.round_num.values
                #     dot_ids = group.dot_id.values
                #     rounds_num = rounds_num.astype(int)
                #     barcode = np.zeros([self.barcode_length],dtype=np.int8)
                #     barcode[(rounds_num-1)] += 1

                #     dists_arr, index_arr = nn_sklearn.kneighbors(barcode.reshape(1, -1), return_distance=True)
                #     gene=self.codebook_df.loc[index_arr.reshape(index_arr.shape[0]),'Gene'].tolist()[0]
                #     self.barcoded_fov_df.loc[self.barcoded_fov_df.barcode_reference_dot_id == name,'raw_barcodes'] = barcode.tostring()
                #     self.barcoded_fov_df.loc[self.barcoded_fov_df.barcode_reference_dot_id == name,'all_Hdistance_genes'] = gene
                #     if dists_arr[0][0] == 0:
                #         self.barcoded_fov_df.loc[self.barcoded_fov_df.barcode_reference_dot_id == name,'0Hdistance_genes'] = gene
                #     elif dists_arr[0][0] < hd_2:
                #         self.barcoded_fov_df.loc[self.barcoded_fov_df.barcode_reference_dot_id == name,'below2Hdistance_genes'] = gene
                #     elif dists_arr[0][0] < hd_3:
                #         self.barcoded_fov_df.loc[self.barcoded_fov_df.barcode_reference_dot_id == name,'below3Hdistance_genes'] = gene


                barcode_reference_dot_id_list = []
                num_unique_dots = np.unique(self.barcoded_fov_df.loc[:,'barcode_reference_dot_id']).shape[0]
                
                # There are no dots is the df
                if num_unique_dots > 0:
                
                    all_barcodes = np.zeros([num_unique_dots,self.barcode_length],dtype=np.int8)
                    for idx, (name, group) in enumerate(self.grpd):
                        barcode_reference_dot_id_list.append(name)
                        barcode = np.zeros([self.barcode_length],dtype=np.int8)
                        rounds_num = group.round_num.values
                        rounds_num = rounds_num.astype(int)
                        barcode[(rounds_num-1)] += 1
                        all_barcodes[idx,:] = barcode

                    dists_arr, index_arr = nn_sklearn.kneighbors(all_barcodes, return_distance=True)
                    genes=self.codebook_df.loc[index_arr.reshape(index_arr.shape[0]),'Gene'].tolist()
                    for idx,name in enumerate(barcode_reference_dot_id_list):
                        barcode = all_barcodes[idx,:]
                        gene = genes[idx]
                        hd = dists_arr[idx][0]

                        cols = ['raw_barcodes','all_Hdistance_genes','number_positive_bits','hamming_distance'] # will add last column depending on hd
                        writing_data = [barcode.tostring(),gene,barcode.sum(),hd]
                        if hd == 0:
                            cols = cols + ['zeroHdistance_genes']
                            writing_data = writing_data + [gene]

                        if hd < hd_2:
                            cols = cols + ['below2Hdistance_genes']
                            writing_data = writing_data + [gene]

                        if hd < hd_3:
                            cols = cols + ['below3Hdistance_genes']
                            writing_data = writing_data + [gene]

                        self.barcoded_fov_df.loc[self.barcoded_fov_df.barcode_reference_dot_id == name,cols] = writing_data


                        # self.barcoded_fov_df.loc[self.barcoded_fov_df.barcode_reference_dot_id == name,'raw_barcodes'] = barcode.tostring()
                        # self.barcoded_fov_df.loc[self.barcoded_fov_df.barcode_reference_dot_id == name,'all_Hdistance_genes'] = gene
                        # self.barcoded_fov_df.loc[self.barcoded_fov_df.barcode_reference_dot_id == name,'number_positive_bits'] = barcode.sum()
                        # self.barcoded_fov_df.loc[self.barcoded_fov_df.barcode_reference_dot_id == name,'hamming_distance'] = hd

                        # if hd == 0:
                        #     self.barcoded_fov_df.loc[self.barcoded_fov_df.barcode_reference_dot_id == name,'0Hdistance_genes'] = gene
                        # elif hd < hd_2:
                        #     self.barcoded_fov_df.loc[self.barcoded_fov_df.barcode_reference_dot_id == name,'below2Hdistance_genes'] = gene
                        # elif hd < hd_3:
                        #     self.barcoded_fov_df.loc[self.barcoded_fov_df.barcode_reference_dot_id == name,'below3Hdistance_genes'] = gene
   
        fname = self.file_tags['experiment_fpath'] / 'tmp' / 'registered_counts' / (self.file_tags['experiment_name'] + '_' + self.file_tags['channel'] + '_decoded_fov_' + self.file_tags['fov'] + '.parquet')
        self.barcoded_fov_df.to_parquet(fname,index=False)
   


class extract_barcodes_NN_test():
    """
    Class used to extract the barcodes from the registered
    counts using nearest neighbour

    Parameters:
    -----------
    counts: pandas.DataFrame
        pandas file with the fov counts after
        registration
    analysis_parameters: dict
        parameters for data processing
    experiment_config: Dict
        dictionary with the experimental data 
    codebook_df: pandas.DataFrame
        pandas file with the codebook used to
        deconvolve the barcode

    NB: if there is a problem with the registration the barcode assigned 
        will be 0*barcode_length
    
    """

    def __init__(self, fov, channel, counts, analysis_parameters:Dict,experiment_config:Dict,codebook_df,status:str):
        
        self.barcodes_extraction_resolution = analysis_parameters['BarcodesExtractionResolution']
        self.RegistrationMinMatchingBeads = analysis_parameters['RegistrationMinMatchingBeads']
        self.barcode_length = experiment_config['Barcode_length']
        self.fov = fov
        self.channel = channel
        self.counts = counts
        self.logger = selected_logger()
        self.codebook_df = codebook_df
        self.status = status
        self.registration_errors = Registration_errors()



        
    @staticmethod
    def barcode_nn(counts_df, ref_round_number, barcodes_extraction_resolution):
        column_names = list(counts_df.columns.values)
        column_names = column_names.append('barcode_reference_dot_id')
        barcoded_df = pd.DataFrame(columns=column_names)

        reference_array = counts_df.loc[counts_df.round_num == ref_round_number, ['r_px_registered','c_px_registered']].to_numpy()
        reference_round_df = counts_df.loc[counts_df.round_num == ref_round_number,:].reset_index(drop=True)
        # Step one (all dots not in round 1)
        coords_compare = counts_df.loc[counts_df.round_num != ref_round_number, ['r_px_registered','c_px_registered']].to_numpy()
        compare_df = counts_df.loc[counts_df.round_num != ref_round_number,:].reset_index(drop=True)

        if (reference_array.shape[0] >0) and (coords_compare.shape[0] >0):
            # initialize network
            nn = NearestNeighbors(1, metric="euclidean")
            nn.fit(reference_array)

            # Get the nn
            dists, indices = nn.kneighbors(coords_compare, return_distance=True)

            # select only the nn that are below barcodes_extraction_resolution distance
            idx_selected_coords_compare = np.where(dists <= barcodes_extraction_resolution)[0]

            compare_selected_df = compare_df.loc[idx_selected_coords_compare,:]
            compare_selected_df['barcode_reference_dot_id'] = np.nan

            # ref_idx = indices[idx_selected_coords_compare]
            # compare_selected_df.loc[compare_selected_df.index.isin(idx_selected_coords_compare),'barcode_reference_dot_id'] = reference_round_df.loc[ref_idx,'dot_id'].values[0]
            
            for idx in idx_selected_coords_compare:
                ref_idx = indices[idx]
                compare_selected_df.loc[idx,'barcode_reference_dot_id'] = reference_round_df.loc[ref_idx,'dot_id'].values[0]

            reference_round_df['barcode_reference_dot_id'] = reference_round_df.dot_id

            barcoded_df = barcoded_df.append([compare_selected_df, reference_round_df], ignore_index=True)

            compare_df = compare_df.drop(compare_selected_df.index)
            compare_df = compare_df.reset_index(drop=True)
    
        return compare_df, barcoded_df

    @staticmethod
    def convert_str_codebook(codebook_df,column_name):
        codebook_df[column_name] = codebook_df[column_name].map(lambda x: np.frombuffer(x, np.int8))
        return codebook_df

    @staticmethod
    def make_codebook_array(codebook_df,column_name):
        codebook_array = np.zeros((len(codebook_df[column_name]),codebook_df[column_name][0].shape[0]))
        for idx, el in enumerate(codebook_df[column_name]):
            row = codebook_df[column_name][idx]
            row = row[np.newaxis,:]
            codebook_array[idx,:] = row
        return codebook_array



    def run_extraction(self):

        data_models = Output_models()
        registration_errors = Registration_errors()
        self.barcoded_fov_df = data_models.barcode_analysis_df
        self.barcoded_fov_df.attrs = self.counts.attrs

        if self.status == 'FAILED':
            error = self.counts['min_number_matching_dots_registration'].values[0]
            round_num = self.counts['round_num'].values[0]
            self.barcoded_fov_df = self.barcoded_fov_df.append({'min_number_matching_dots_registration':error,
                                                           'fov_num':int(self.fov),'dot_channel':self.channel,'round_num': round_num },ignore_index=True)
        elif self.status == 'SUCCESS':

            if (min(self.counts.loc[:,'min_number_matching_dots_registration']) < self.RegistrationMinMatchingBeads):
                round_num = self.counts['round_num'].values[0]
                self.barcoded_fov_df = self.barcoded_fov_df.append({'min_number_matching_dots_registration':registration_errors.registration_below_extraction_resolution, 
                                            'fov_num':int(self.fov),'dot_channel':self.channel,'round_num': round_num},ignore_index=True)
                self.status = 'FAILED'
            else:
                hd_2 = 2 / self.barcode_length
                hd_3 = 3 / self.barcode_length
                # barcode_length = len(self.counts['round_num'].unique())
                rounds = np.arange(1,self.barcode_length+1)
                self.codebook_df = self.convert_str_codebook(self.codebook_df,'Code')
                codebook_array = self.make_codebook_array(self.codebook_df,'Code')
                nn_sklearn = NearestNeighbors(n_neighbors=1, metric="hamming")
                nn_sklearn.fit(codebook_array)

                # remove points with np.NAN
                # self.counts = self.counts.dropna()
                for round_num in rounds:
                    compare_df, barcoded_df = self.barcode_nn(self.counts, round_num, self.barcodes_extraction_resolution)
                    self.barcoded_fov_df = self.barcoded_fov_df.append(barcoded_df, ignore_index=True)
                    self.counts = compare_df

                self.counts['barcode_reference_dot_id'] = self.counts.dot_id
                self.barcoded_fov_df = self.barcoded_fov_df.append(self.counts, ignore_index=True)
                self.barcoded_fov_df['barcodes_extraction_resolution'] = self.barcodes_extraction_resolution
                self.grpd = self.barcoded_fov_df.groupby('barcode_reference_dot_id')

                barcode_reference_dot_id_list = []
                num_unique_dots = np.unique(self.barcoded_fov_df.loc[:,'barcode_reference_dot_id']).shape[0]
                
                # There are no dots is the df
                if num_unique_dots > 0:
                
                    all_barcodes = np.zeros([num_unique_dots,self.barcode_length],dtype=np.int8)
                    for idx, (name, group) in enumerate(self.grpd):
                        barcode_reference_dot_id_list.append(name)
                        barcode = np.zeros([self.barcode_length],dtype=np.int8)
                        rounds_num = group.round_num.values
                        rounds_num = rounds_num.astype(int)
                        barcode[(rounds_num-1)] += 1
                        all_barcodes[idx,:] = barcode

                    dists_arr, index_arr = nn_sklearn.kneighbors(all_barcodes, return_distance=True)
                    genes=self.codebook_df.loc[index_arr.reshape(index_arr.shape[0]),'Gene'].tolist()
                    for idx,name in enumerate(barcode_reference_dot_id_list):
                        barcode = all_barcodes[idx,:]
                        gene = genes[idx]
                        hd = dists_arr[idx][0]

                        cols = ['raw_barcodes','all_Hdistance_genes','number_positive_bits','hamming_distance'] # will add last column depending on hd
                        writing_data = [barcode.tostring(),gene,barcode.sum(),hd]
                        if hd == 0:
                            cols = cols + ['zeroHdistance_genes']
                            writing_data = writing_data + [gene]

                        if hd < hd_2:
                            cols = cols + ['below2Hdistance_genes']
                            writing_data = writing_data + [gene]

                        if hd < hd_3:
                            cols = cols + ['below3Hdistance_genes']
                            writing_data = writing_data + [gene]

                        self.barcoded_fov_df.loc[self.barcoded_fov_df.barcode_reference_dot_id == name,cols] = writing_data






class extract_barcodes_NN_new():
    """
    Class used to extract the barcodes from the registered
    counts using nearest neighbour

    Parameters:
    -----------
    counts: pandas.DataFrame
        pandas file with the fov counts after
        registration
    analysis_parameters: dict
        parameters for data processing
    experiment_config: Dict
        dictionary with the experimental data 
    codebook_df: pandas.DataFrame
        pandas file with the codebook used to
        deconvolve the barcode

    NB: if there is a problem with the registration the barcode assigned 
        will be 0*barcode_length
    
    """

    def __init__(self, registered_counts, analysis_parameters:Dict,experiment_config:Dict,codebook_df):
        
        self.counts_df = registered_counts
        self.analysis_parameters = analysis_parameters
        self.experiment_config = experiment_config
        self.codebook_df = codebook_df
        
        self.logger = selected_logger()
        
        self.barcodes_extraction_resolution = analysis_parameters['BarcodesExtractionResolution']
        self.RegistrationMinMatchingBeads = analysis_parameters['RegistrationMinMatchingBeads']
        self.barcode_length = self.counts_df.loc[0]['barcode_length']
        self.registration_errors = Registration_errors()

        self.stitching_channel = self.counts_df['stitching_channel'].iloc[0]
        
    @staticmethod
    def convert_str_codebook(codebook_df,column_name):
        codebook_df[column_name] = codebook_df[column_name].map(lambda x: np.frombuffer(x, np.int8))
        return codebook_df

    @staticmethod
    def make_codebook_array(codebook_df,column_name):
        codebook_array = np.zeros((len(codebook_df[column_name]),codebook_df[column_name][0].shape[0]))
        for idx, el in enumerate(codebook_df[column_name]):
            row = codebook_df[column_name][idx]
            row = row[np.newaxis,:]
            codebook_array[idx,:] = row
        return codebook_array

    
    @staticmethod
    def barcode_nn(counts_df, ref_round_number, barcodes_extraction_resolution):
        column_names = list(counts_df.columns.values)
        column_names = column_names.append('barcode_reference_dot_id')
        barcoded_df = pd.DataFrame(columns=column_names)

        
        reference_array = counts_df.loc[counts_df.round_num == ref_round_number, ['r_px_registered','c_px_registered']].to_numpy()
        reference_round_df = counts_df.loc[counts_df.round_num == ref_round_number,:].reset_index(drop=True)
    
        # Step one (all dots not in round 1)
        coords_compare = counts_df.loc[counts_df.round_num != ref_round_number, ['r_px_registered','c_px_registered']].to_numpy()
        compare_df = counts_df.loc[counts_df.round_num != ref_round_number,:].reset_index(drop=True)
        
        if (reference_array.shape[0] >0) and (coords_compare.shape[0] >0):
            # initialize network
            nn = NearestNeighbors(1, metric="euclidean")
            nn.fit(reference_array)

            # Get the nn
            dists, indices = nn.kneighbors(coords_compare, return_distance=True)

            # select only the nn that are below barcodes_extraction_resolution distance
            idx_selected_coords_compare = np.where(dists <= barcodes_extraction_resolution)[0]
            compare_selected_df = compare_df.loc[idx_selected_coords_compare,:]
            compare_selected_df['barcode_reference_dot_id'] = np.nan
            
            for k,v in groupby(idx_selected_coords_compare):
                if len(list(v)) > 3:
                    print("key: '{}'--> group: {}".format(k, len(list(v))))
            # ref_idx = indices[idx_selected_coords_compare].squeeze()
            # compare_selected_df.loc[compare_selected_df.index.isin(idx_selected_coords_compare),'barcode_reference_dot_id'] = reference_round_df.loc[ref_idx,'dot_id'].values[0]

            for idx in idx_selected_coords_compare:
                ref_idx = indices[idx]
                compare_selected_df.loc[idx,'barcode_reference_dot_id'] = reference_round_df.loc[ref_idx,'dot_id'].values[0]

            reference_round_df['barcode_reference_dot_id'] = reference_round_df.dot_id

            barcoded_df = barcoded_df.append([compare_selected_df, reference_round_df], ignore_index=True)

            compare_df = compare_df.drop(compare_selected_df.index)
            compare_df = compare_df.reset_index(drop=True)
            
        return compare_df, barcoded_df

        
    def run_extraction(self):

        data_models = Output_models()
        registration_errors = Registration_errors()
        self.barcoded_spec = data_models.barcode_analysis_df

        if not self.counts_df[self.counts_df['dot_id'].isnull()].empty:
            print('shitty FOV')
            self.all_combine_df = pd.concat([self.counts_df,self.barcoded_spec],axis=1)

        elif (min(self.counts_df['min_number_matching_dots_registration']) < self.RegistrationMinMatchingBeads):

            self.counts_df['min_number_matching_dots_registration'] = registration_errors.registration_below_extraction_resolution
            self.all_combine_df = pd.concat([self.counts_df,self.barcoded_spec],axis=1)

        else:
            self.counts_df = pd.concat([self.counts_df,self.barcoded_spec],axis=1)
            self.fish_counts = self.counts_df.loc[self.counts_df.channel != self.stitching_channel,:]
            
            hd_2 = 2 / self.barcode_length
            hd_3 = 3 / self.barcode_length
            rounds = np.arange(1,self.barcode_length+1)
            self.codebook_df = self.convert_str_codebook(self.codebook_df,'Code')
            codebook_array = self.make_codebook_array(self.codebook_df,'Code')
            nn_sklearn = NearestNeighbors(n_neighbors=1, metric="hamming")
            nn_sklearn.fit(codebook_array)
            
    
            self.index_df = pd.DataFrame(index=self.fish_counts.index,
                                        columns = ['barcode_reference_dot_id','raw_barcodes'])
            
            barcoded_df_list = []
            self.barcoded_fov_df = pd.DataFrame()
            for round_num in rounds:
                compare_df, barcoded_df = self.barcode_nn(self.fish_counts, round_num, self.barcodes_extraction_resolution)
                barcoded_df_list.append(barcoded_df)
                self.fish_counts = compare_df
            self.barcoded_fov_df = pd.concat(barcoded_df_list, ignore_index=True)
                
            self.fish_counts['barcode_reference_dot_id'] = self.fish_counts.dot_id
            self.barcoded_fov_df = self.barcoded_fov_df.append(self.fish_counts, ignore_index=True)
            self.barcoded_fov_df['barcodes_extraction_resolution'] = self.barcodes_extraction_resolution        
            self.grpd = self.barcoded_fov_df.groupby('barcode_reference_dot_id')

            barcode_reference_dot_id_list = []
            num_unique_dots = np.unique(self.barcoded_fov_df.loc[:,'barcode_reference_dot_id']).shape[0]

            # There are no dots is the df
            if num_unique_dots > 0:

                all_barcodes = np.zeros([num_unique_dots,self.barcode_length],dtype=np.int8)
                for idx, (name, group) in enumerate(self.grpd):
                    barcode_reference_dot_id_list.append(name)
                    barcode = np.zeros([self.barcode_length],dtype=np.int8)
                    rounds_num = group.round_num.values
                    rounds_num = rounds_num.astype(int)
                    barcode[(rounds_num-1)] += 1
                    all_barcodes[idx,:] = barcode

                dists_arr, index_arr = nn_sklearn.kneighbors(all_barcodes, return_distance=True)
                genes=self.codebook_df.loc[index_arr.reshape(index_arr.shape[0]),'Gene'].tolist()
                
                self.all_combine_df_list = []
                self.all_combine_df = pd.DataFrame()
                
                for idx, (name, group) in enumerate(self.grpd):
                    barcode = all_barcodes[idx,:]
                    gene = genes[idx]
                    hd = dists_arr[idx][0]

                    group_df = group
                    
                    group_df['raw_barcodes'] = barcode.tostring()
                    group_df['all_Hdistance_genes'] = gene
                    group_df['number_positive_bits'] = barcode.sum()
                    group_df['hamming_distance'] = hd
                    group_df['zeroHdistance_genes'] = np.nan
                    group_df['below2Hdistance_genes'] = np.nan
                    group_df['below3Hdistance_genes'] = np.nan
                    
                    if hd == 0:
                        group_df['zeroHdistance_genes'] = gene
                    if hd < hd_2:
                        group_df['below2Hdistance_genes'] = gene
                    if hd < hd_3:
                        group_df['below3Hdistance_genes'] = gene
                    
                    self.all_combine_df_list.append(group_df)
                
                # self.all_combine_df = pd.concat(self.all_combine_df_list, axis=0,copy=False) 
                chunk_size = 500
                self.all_combine_df = merge_with_concat((merge_with_concat(dfs) for dfs in chunk_dfs(self.all_combine_df_list, chunk_size)))
  
            else:
                self.all_combine_df = self.fish_counts   # add the missing column and control error for stitching


def decoder_fun(registered_counts, analysis_parameters,experiment_info,codebook_df):
    dc = extract_barcodes_NN_new(registered_counts, analysis_parameters,experiment_info,codebook_df)
    dc.run_extraction()
    return dc.all_combine_df


def convert_str_codebook(codebook_df,column_name):
    codebook_df[column_name] = codebook_df[column_name].map(lambda x: np.frombuffer(x, np.int8))
    return codebook_df

def make_codebook_array(codebook_df,column_name):
    codebook_array = np.zeros((len(codebook_df[column_name]),codebook_df[column_name][0].shape[0]))
    for idx, el in enumerate(codebook_df[column_name]):
        row = codebook_df[column_name][idx]
        row = row[np.newaxis,:]
        codebook_array[idx,:] = row
    return codebook_array






#                 chunk_size = 100
#                 self.all_combine_df = merge_with_concat((merge_with_concat(dfs) for dfs in chunk_dfs(self.all_combine_df_list, chunk_size)))
                

def extract_barcodes_NN_fast(registered_counts_df, analysis_parameters:Dict,codebook_df):
    """
    Class used to extract the barcodes from the registered
    counts using nearest neighbour

    Parameters:
    -----------
    counts: pandas.DataFrame
        pandas file with the fov counts after
        registration
    analysis_parameters: dict
        parameters for data processing 
    codebook_df: pandas.DataFrame
        pandas file with the codebook used to
        deconvolve the barcode

    NB: if there is a problem with the registration the barcode assigned 
        will be 0*barcode_length
    
    """

    logger = selected_logger()

    barcodes_extraction_resolution = analysis_parameters['BarcodesExtractionResolution']
    RegistrationMinMatchingBeads = analysis_parameters['RegistrationMinMatchingBeads']
    barcode_length = registered_counts_df.loc[0]['barcode_length']
    registration_errors = Registration_errors()

    stitching_channel = registered_counts_df['stitching_channel'].iloc[0]
    
    fish_counts = registered_counts_df.loc[registered_counts_df.channel != stitching_channel,:]
    stitching_channel_counts = registered_counts_df.loc[registered_counts_df.channel == stitching_channel,:]
    
    fish_counts.dropna(subset=['dot_id'],inplace=True)
    # Starting level for selection of dots
    dropping_counts = fish_counts.copy(deep=True)

    all_decoded_dots_list = []

    if fish_counts['r_px_registered'].isnull().values.any():
        
        all_decoded_dots_df = pd.DataFrame(columns = fish_counts.columns)
        all_decoded_dots_df['decoded_genes'] = np.nan
        all_decoded_dots_df['hamming_distance'] = np.nan
        all_decoded_dots_df['number_positive_bits'] = np.nan
        all_decoded_dots_df['barcode_reference_dot_id'] = np.nan
        all_decoded_dots_df['raw_barcodes'] = np.nan
        all_decoded_dots_df['barcodes_extraction_resolution'] = barcodes_extraction_resolution
        # Save barcoded_round and all_decoded_dots_df
        return fish_counts, all_decoded_dots_df
        
    else:
    
        for ref_round_number in np.arange(1,barcode_length+1):

            #ref_round_number = 1
            reference_round_df = dropping_counts.loc[dropping_counts.round_num == ref_round_number,:]
            # Step one (all dots not in round 1)
            compare_df = dropping_counts.loc[dropping_counts.round_num != ref_round_number,:]

            if (not reference_round_df.empty) and (not compare_df.empty):
                nn = NearestNeighbors(1, metric="euclidean")
                nn.fit(reference_round_df[['r_px_registered','c_px_registered']])
                dists, indices = nn.kneighbors(compare_df[['r_px_registered','c_px_registered']], return_distance=True)

                # select only the nn that are below barcodes_extraction_resolution distance
                idx_distances_below_resolution = np.where(dists <= barcodes_extraction_resolution)[0]

                comp_idx = idx_distances_below_resolution
                ref_idx = indices[comp_idx].flatten()

                # Subset the dataframe according to the selected points
                # The reference selected will have repeated points
                comp_selected_df = compare_df.iloc[comp_idx]
                ref_selected_df = reference_round_df.iloc[ref_idx]

               # The size of ref_selected_df w/o duplicates may be smaller of reference_round_df if 
                # some of the dots in reference_round_df have no neighbours

                # Test approach where we get rid of the single dots
                comp_selected_df.loc[:,'barcode_reference_dot_id'] = ref_selected_df['dot_id'].values
                ref_selected_df_no_duplicates = ref_selected_df.drop_duplicates()
                ref_selected_df_no_duplicates.loc[:,'barcode_reference_dot_id'] = ref_selected_df_no_duplicates['dot_id'].values

                barcoded_round = pd.concat([comp_selected_df, ref_selected_df_no_duplicates], axis=0,ignore_index=False)
                barcoded_round_grouped = barcoded_round.groupby('barcode_reference_dot_id')
                for brdi, grp in barcoded_round_grouped:
                    barcode = np.zeros([barcode_length],dtype=np.int8)
                    barcode[grp.round_num.values.astype(np.int8)-1] = 1
                    #hamming_dist, index_gene = nn_sklearn.kneighbors(barcode.reshape(1, -1), return_distance=True)
                    #gene= codebook_df.loc[index_gene.reshape(index_gene.shape[0]),'Gene'].tolist()

                    barcode = barcode.tostring()
                    ref_selected_df_no_duplicates.loc[ref_selected_df_no_duplicates.barcode_reference_dot_id == brdi,'raw_barcodes'] = barcode
                    #ref_selected_df_no_duplicates.loc[ref_selected_df_no_duplicates.barcode_reference_dot_id == brdi,'decoded_gene_name'] = gene
                    #ref_selected_df_no_duplicates.loc[ref_selected_df_no_duplicates.barcode_reference_dot_id == brdi,'hamming_distance'] = hamming_dist.flatten()[0]

                    #fish_counts.loc[grp.index,'barcode_reference_dot_id'] = brdi
                    #fish_counts.loc[grp.index,'raw_barcodes'] = barcode
                    #dists, index = nn_sklearn.kneighbors(all_barcodes, return_distance=True)

                all_decoded_dots_list.append(ref_selected_df_no_duplicates)

                all_decoded_dots_df = pd.concat(all_decoded_dots_list,ignore_index=False)
                
                compare_df = compare_df.drop(comp_selected_df.index)
                dropping_counts = compare_df

        codebook_df = convert_str_codebook(codebook_df,'Code')
        codebook_array = make_codebook_array(codebook_df,'Code')
        nn_sklearn = NearestNeighbors(n_neighbors=1, metric="hamming")
        nn_sklearn.fit(codebook_array)

        all_barcodes = np.vstack(all_decoded_dots_df.raw_barcodes.map(lambda x: np.frombuffer(x, np.int8)).values)
        dists_arr, index_arr = nn_sklearn.kneighbors(all_barcodes, return_distance=True)
        genes=codebook_df.loc[index_arr.reshape(index_arr.shape[0]),'Gene'].tolist()

        all_decoded_dots_df.loc[:,'decoded_genes'] = genes
        all_decoded_dots_df.loc[:,'hamming_distance'] = dists_arr
        all_decoded_dots_df.loc[:,'number_positive_bits'] = all_barcodes.sum(axis=1)

        all_decoded_dots_df = pd.concat([all_decoded_dots_df,stitching_channel_counts])
        all_decoded_dots_df['barcodes_extraction_resolution'] = barcodes_extraction_resolution

        # Save barcoded_round and all_decoded_dots_df
        return barcoded_round, all_decoded_dots_df



#                 chunk_size = 100
#                 self.all_combine_df = merge_with_concat((merge_with_concat(dfs) for dfs in chunk_dfs(self.all_combine_df_list, chunk_size)))
                




   # ---------------------------------------------------




class extract_barcodes():
    """
    Class used to extract the barcodes from the registered
    counts

    Parameters:
    -----------
    counts: pandas dataframe
        contains all the counts for a fov
    pxl: int
        size of the hood to search for positive
        barcodes
    """

    def __init__(self, counts, pxl:int,save=True):
        self.counts = counts
        self.pxl = pxl
        self.save = save

        self.logger = logging.getLogger(__name__)

    @staticmethod
    def combine_coords(counts, round_num, counts_type='registered'):
        data_reference = counts.loc[counts['round_num'] == round_num]
        if counts_type == 'registered':
            r_px = data_reference.r_px_registered.to_list()
            c_px = data_reference.c_px_registered.to_list()
        elif counts_type == 'original':
            r_px = data_reference.r_px_original.to_list()
            c_px = data_reference.c_px_original.to_list()
        coords = np.array(list(zip(r_px,c_px)))
        position_idx = data_reference.index
        return coords, position_idx
    
    @staticmethod
    def dots_hoods(coords,pxl):
        r_tl = coords[:,0]-pxl
        r_br = coords[:,0]+pxl
        c_tl = coords[:,1]-pxl
        c_tr = coords[:,1]+pxl
        r_tl = r_tl[:,np.newaxis]
        r_br = r_br[:,np.newaxis]
        c_tl = c_tl[:,np.newaxis]
        c_tr = c_tr[:,np.newaxis]
        chunks_coords = np.hstack((r_tl,r_br,c_tl,c_tr))
        chunks_coords = chunks_coords.astype(int)
        return chunks_coords

    @staticmethod
    def barcode_detection(all_coords_image,chunk_coords,pxl):
        selected_region = all_coords_image[:,chunk_coords[0]:chunk_coords[1]+1,chunk_coords[2]:chunk_coords[3]+1]
        barcode = np.sum(selected_region, axis=(1,2))
        return barcode, selected_region

    def run_extraction(self):
        # add check to make sure that the counts are there
        column_names = ['r_px_original','c_px_original','dot_id','fov_num','round_num','dot_intensity_norm',
                        'dot_intensity_not','selected_thr' ,'channel','r_shift','c_shift','r_px_registered',
                        'c_px_registered','pixel_microns']
        column_names_output = column_names.copy()
        column_names_output.append('pxl_hood_size')
        column_names_output.append('raw_barcodes')
        self.barcodes_binary_images = {}
        if len(self.counts['r_px_registered']) > 1:
            r_num = int(self.counts['r_px_registered'].max())
            c_num = int(self.counts['c_px_registered'].max())
            rounds_num = len(self.counts['round_num'].unique())
            all_coords_image = np.zeros([rounds_num,r_num+1,c_num+1],dtype=np.int8)
            for round_num in np.arange(1,rounds_num+1):
                coords, ref_position_idx = self.combine_coords(self.counts, round_num, counts_type='registered')
                if np.any(np.isnan(coords)):
                    stop = True
                    break
                else:
                    stop = False
                    coords = coords.astype(np.int)
                    all_coords_image[round_num-1,coords[:,0],coords[:,1]] = 1

            if stop:
                self.logger.error(f'missing round {coords} no barcodes')
                self.output = pd.DataFrame(np.nan,index=[0],columns = column_names_output)
            else:
                self.output = pd.DataFrame(columns = column_names_output)
                self.output['pxl_hood_size'] = np.nan
                self.output['raw_barcodes'] = np.nan
                # self.raw_barcodes = []
                for round_num in np.arange(1,rounds_num+1):
                    coords, ref_position_idx = self.combine_coords(self.counts, round_num, counts_type='registered')
                    dot_ids = self.counts.loc[ref_position_idx,'dot_id']
                    dot_ids = dot_ids.reset_index()
                    chunks_coords = self.dots_hoods(coords,self.pxl)
                    for row in np.arange(coords.shape[0]):
                        chunk = chunks_coords[row,:]
                        dot_id = dot_ids.loc[row,'dot_id']
                        selected_dot_coords = coords[row,:].astype(int)
                        barcode,selected_region = self.barcode_detection(all_coords_image,chunk,self.pxl)
                        # Add only the barcodes that are not all negatives (can happen because I am blanking the img)
                        if barcode.sum():
                            self.barcodes_binary_images[dot_id] = {}
                            self.barcodes_binary_images[dot_id]['img'] = selected_region
                            self.barcodes_binary_images[dot_id]['raw_barcode'] = barcode
                            # self.raw_barcodes.append(barcode)
                            # Blank the image
                            all_coords_image[:,selected_dot_coords[0],selected_dot_coords[1]] = 0
                            dot_data = self.counts[(self.counts['round_num'] == round_num) & (self.counts['r_px_registered'] == selected_dot_coords[0]) & (self.counts['c_px_registered'] == selected_dot_coords[1])] 
                            barcode_st = ''
                            barcode_st = barcode_st.join([",".join(item) for item in barcode.astype(str)])
                            dot_data.insert(0,'raw_barcodes',barcode_st)
                            dot_data.insert(0,'pxl_hood_size',self.pxl)
                            self.output = pd.concat([self.output,dot_data],axis=0,ignore_index=True)
        
                # self.output['pxl_hood_size'] = self.pxl
                # self.output['raw_barcodes'] = self.raw_barcodes
        else:
            if self.save:
                self.output = pd.DataFrame(np.nan,index=[0],columns = column_names_output)





# class extract_barcodes_NN_descend():
#     """
#     Class used to extract the barcodes from the registered
#     counts using nearest neighbour

#     Parameters:
#     -----------
#     counts: pandas dataframe
#         contains all the counts for a fov
#     pxl: int
#         size of the hood to search for positive
#         barcodes

#     NB: If a round is missing everything get processed
#         but the missing rounds get 0 in the barcode
#     """

#     def __init__(self, counts, pxl:int,total_rounds:int):
#         self.counts = counts
#         self.pxl = pxl
#         self.total_rounds = total_rounds

#         #self.logger = logging.getLogger(__name__)

        
#     @staticmethod
#     def barcode_nn(counts_df, ref_round_number, pxl):
#         column_names = list(counts_df.columns.values)
#         column_names = column_names.append('barcode_reference_dot_id')
#         barcoded_df = pd.DataFrame(columns=column_names)

#         reference_array = counts_df.loc[counts_df.round_num == ref_round_number, ['r_px_registered','c_px_registered']].to_numpy()
#         reference_round_df = counts_df.loc[counts_df.round_num == ref_round_number,:].reset_index(drop=True)
#         # Step one (all dots not in round 1)
#         coords_compare = counts_df.loc[counts_df.round_num != ref_round_number, ['r_px_registered','c_px_registered']].to_numpy()
#         compare_df = counts_df.loc[counts_df.round_num != ref_round_number,:].reset_index(drop=True)

#         if (reference_array.shape[0] >0) and (coords_compare.shape[0] >0):
#             # initialize network
    
#             index = NNDescent(reference_array,metric='euclidean',n_neighbors=1)        
#             # Get the nn
#             indices, dists = index.query(coords_compare,k=1)
#             # select only the nn that are below pxl distance
#             idx_selected_coords_compare = np.where(dists <= pxl)[0]

#             compare_selected_df = compare_df.loc[idx_selected_coords_compare,:]
#             compare_selected_df['barcode_reference_dot_id'] = np.nan

#             for idx in idx_selected_coords_compare:
#                 ref_idx = indices[idx]
#                 compare_selected_df.loc[idx,'barcode_reference_dot_id'] = reference_round_df.loc[ref_idx,'dot_id'].values[0]

#             reference_round_df['barcode_reference_dot_id'] = reference_round_df.dot_id

#             barcoded_df = barcoded_df.append([compare_selected_df, reference_round_df], ignore_index=True)

#             compare_df = compare_df.drop(compare_selected_df.index)
#             compare_df = compare_df.reset_index(drop=True)
    
#         return compare_df, barcoded_df


#     def run_extraction(self):

#         # count all barcodes
#         columns_names = ['r_px_original', 'c_px_original', 'dot_id', 'fov_num',
#         'round_num', 'dot_intensity_norm', 'dot_intensity_not',
#         'selected_thr', 'channel', 'r_shift', 'c_shift', 'r_px_registered',
#         'c_px_registered', 'pixel_microns', 'pxl_hood_size', 'barcode_reference_dot_id']

#         self.barcoded_fov_df = pd.DataFrame(columns=columns_names)
        
#         rounds = np.arange(1,self.total_rounds+1)
        
#          # Check that counts not NaN
#         if not self.counts['r_px_original'].isnull().all():

#             # remove points with np.NAN
#             self.counts = self.counts.dropna()
#             for round_num in np.arange(1,self.total_rounds):
#                 compare_df, barcoded_df = self.barcode_nn(self.counts, round_num, self.pxl)
#                 self.barcoded_fov_df = self.barcoded_fov_df.append(barcoded_df, ignore_index=True)
#                 self.counts = compare_df

#             self.counts['barcode_reference_dot_id'] = self.counts.dot_id
#             self.barcoded_fov_df = self.barcoded_fov_df.append(self.counts, ignore_index=True)
#             self.barcoded_fov_df['pxl_hood_size'] = self.pxl
#             self.grpd = self.barcoded_fov_df.groupby('barcode_reference_dot_id')
#             self.all_barcodes = {}
#             for name, group in self.grpd:
#                 rounds_num = group.round_num.values
#                 dot_ids = group.dot_id.values
#                 rounds_num = rounds_num.astype(int)
#                 barcode = np.zeros([1,self.total_rounds],dtype=int)
#                 barcode[:,(rounds_num-1)] += 1
#                 barcode_st = ''
#                 barcode_st = barcode_st.join([",".join(item) for item in barcode.astype(str)])
#                 self.all_barcodes[name] = {}
#                 self.all_barcodes[name]['dot_ids'] = dot_ids
#                 self.all_barcodes[name]['barcode'] = barcode
#                 self.all_barcodes[name]['barcode_str'] = barcode_st
#         else:
#             self.all_barcodes[name] = {}


# def convert_str_codebook(codebook_df,column_name):
#     codebook_df[column_name] = codebook_df[column_name].map(lambda x: np.frombuffer(x, np.int8))
#     return codebook_df


# def make_codebook_array(codebook_df,column_name):
#     codebook_array = np.zeros((len(codebook_df[column_name]),codebook_df[column_name][0].shape[0]))
#     for idx, el in enumerate(codebook_df[column_name]):
#         row = codebook_df[column_name][idx]
#         row = row[np.newaxis,:]
#         codebook_array[idx,:] = row
#     return codebook_array

# def nn_decoding_barcodes(counts_df, codebook_df):

#     if counts_df.empty:
#         counts_df['gene'] = np.nan
#         counts_df['hamming_distance_barcode'] = np.nan
#     else:
#         codebook_barcode_array_df = convert_str_codebook(codebook_df,'Code')
#         codebook_array = make_codebook_array(codebook_df,'Code')
#         nn_sklearn = NearestNeighbors(n_neighbors=1, metric="hamming")
#         nn_sklearn.fit(codebook_array)

#         grpd_counts = counts_df.groupby('barcode_reference_dot_id')

#         barcodes_dots_ids = list(grpd_counts.groups.keys())
#         subgroup_counts = counts_df.loc[counts_df.dot_id.isin(barcodes_dots_ids),:].copy()
#         subgroup_counts = subgroup_counts.reset_index()

#         subgroup_counts_array_df = convert_str_codebook(subgroup_counts,'raw_barcodes')
#         subgroup_counts_array = make_codebook_array(subgroup_counts_array_df,'raw_barcodes')

#         dists_arr, index_arr = nn_sklearn.kneighbors(subgroup_counts_array, return_distance=True)

#         all_genes=codebook_df.loc[index_arr.reshape(index_arr.shape[0]),'Gene'].tolist()

#         subgroup_counts_array_df['hamming_distance_barcode'] = dists_arr
#         subgroup_counts_array_df['gene'] = all_genes

#         counts_df['gene'] = np.nan
#         counts_df['hamming_distance_barcode'] = np.nan
#         for ref in subgroup_counts_array_df.dot_id.values:  
#             counts_df.loc[counts_df.barcode_reference_dot_id == ref, ['gene']] = subgroup_counts_array_df.loc[subgroup_counts_array_df.dot_id == ref, ['gene']].values[0][0]
#             counts_df.loc[counts_df.barcode_reference_dot_id == ref, ['hamming_distance_barcode']] = subgroup_counts_array_df.loc[subgroup_counts_array_df.dot_id == ref, ['hamming_distance_barcode']].values[0][0]        

#     return counts_df




# if __name__ == '__main__':
#     import pandas as pd
#     import pickle
#     counts = pd.read_csv('/Users/simone/Documents/local_data_storage/LBEXP20200325_oPool11/counts/roi_0/Cy5/LBEXP20200325_oPool11_roi_0_Cy5_fov_90.csv')
#     pxl = 0.5
#     total_rounds = 16
#     bc = extract_barcodes_NN(counts,pxl,total_rounds)
#     bc.run_extraction()
#     pickle.dump(bc.barcodes_binary_images,open('/Users/simone/Downloads/binary_images.pkl','wb'))
#     print(f'cane')