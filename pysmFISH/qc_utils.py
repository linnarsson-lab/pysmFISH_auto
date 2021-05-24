"""
set of useful function to determine if different process steps
that terminated without error generated the correct output.
Ex. missing positions after parsing the data
"""

from typing import *
import re
import sys
import pandas as pd
import numpy as np
from pathlib import Path

from dask import dataframe as dd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from pysmFISH.logger_utils import selected_logger
from pysmFISH.configuration_files import load_experiment_config_file


class QC_registration_error():

    def __init__(self, client, experiment_fpath, analysis_parameters, tiles_coords):
    # def __init__(self, client, experiment_fpath, analysis_parameters, tiles_coords, img_width, img_height):
        self.client = client
        self.experiment_fpath = Path(experiment_fpath)
        self.analysis_parameters = analysis_parameters
        self.tiles_coords = tiles_coords
        # self.img_width = img_width
        # self.img_height = img_height
        matplotlib.use("Agg")

    def create_error_df(self):
        all_counts_folder = self.experiment_fpath / 'results'
        search_key = '*_decoded_*'
        self.error_output_df= pd.DataFrame()
        all_counts_dd = dd.read_parquet(all_counts_folder / search_key)
        registration_error_df = all_counts_dd.groupby('fov_num').agg({'min_number_matching_dots_registration': ['min']}).compute()
        for idx, row in registration_error_df.itertuples():
            search_key = '*decoded_fov_' + str(idx) +'.parquet'
            fname = list(all_counts_folder.glob(search_key))[0]
            fov_data_df = pd.read_parquet(fname)
            matching_rounds = fov_data_df.loc[(fov_data_df.fov_num == idx) &
                        (fov_data_df.min_number_matching_dots_registration == row),
                        ['fov_num','min_number_matching_dots_registration','round_num']]
            self.error_output_df = pd.concat([self.error_output_df,matching_rounds.iloc[0]],axis=1)
        self.error_output_df = self.error_output_df.T
        # #     RegistrationMinMatchingBeads = analysis_parameters['RegistrationMinMatchingBeads']
        #     registration_error_df.columns = ["_".join(x) for x in registration_error_df.columns.ravel()]
        self.error_output_df.to_parquet(self.experiment_fpath / 'tmp' / 'registration_error.parquet')

    def plot_error(self):
        plt.ioff()
        scale_value = 5
        self.tiles_coords = self.tiles_coords / scale_value
        RegistrationMinMatchingBeads = self.analysis_parameters['RegistrationMinMatchingBeads']
        fovs = self.error_output_df['fov_num'].values.astype(int)
        rounds_num = self.error_output_df['round_num'].values.astype(int)
        min_errors = self.error_output_df['min_number_matching_dots_registration'].values.astype(int)

        fig = plt.figure(figsize=(30,20))
        ax = fig.add_subplot(111)

        r_coords = self.tiles_coords[:,0]
        c_coords = self.tiles_coords[:,1]
        r_coords_min = r_coords.min()
        c_coords_min = c_coords.min()
        to_zero_r_coords = r_coords - r_coords_min
        to_zero_c_coords = c_coords - c_coords_min


        errors_normalized = (min_errors -min(min_errors)) / (max(min_errors -min(min_errors)))
        threshold = (RegistrationMinMatchingBeads-min(min_errors)) / (max(min_errors -min(min_errors)))
        nodes = [0,threshold, threshold, 1.0]

        colors = ["black", "black", "blue", "magenta"]
        cmap = LinearSegmentedColormap.from_list("", list(zip(nodes, colors)))
        cmap.set_under("black")

    
        sc = ax.scatter(to_zero_c_coords,to_zero_r_coords,c=errors_normalized,cmap=cmap, s = 1000)
        plt.gca().invert_yaxis()
        for fov, round_num, match, x, y in zip(fovs,rounds_num, min_errors, to_zero_c_coords,to_zero_r_coords):
            
            ax.annotate(
                fov,
                xy=(x,y), xytext=(-0, 15),
                textcoords='offset points', ha='center', va='bottom',fontsize=12)
            
            ax.annotate(round_num, (x, y), color='white', weight='bold', 
                            fontsize=10, ha='center', va='center')
            
            ax.annotate('m' + str(match), xy=(x, y), xytext=(0, -10), color='white', weight='bold', 
                            textcoords='offset points', fontsize=5, ha='center', va='center')
    

        ax.set_aspect('equal')
        ax.axis('off')
        plt.tight_layout()

        plt.savefig(self.experiment_fpath / 'output_figures' / 'registration_error.png',dpi=200,pad_inches=0)


        # Old plotting with squares
        # RegistrationMinMatchingBeads = self.analysis_parameters['RegistrationMinMatchingBeads']
        # # create colormap for error in the registration
        # errors_normalized = (min_errors -min(min_errors)) / (max(min_errors -min(min_errors)))
        # threshold = (RegistrationMinMatchingBeads-min(min_errors)) / (max(min_errors -min(min_errors)))
        # nodes = [0,threshold, threshold, 1.0]

        # colors = ["black", "black", "blue", "magenta"]
        # cmap = LinearSegmentedColormap.from_list("", list(zip(nodes, colors)))
        # cmap.set_under("black")
        

        # sc = ax.scatter(to_zero_c_coords,to_zero_r_coords,c=errors_normalized,cmap=cmap, s = 1000)
        # for fov, round_num, match, x, y in zip(fovs,rounds_num, min_errors, to_zero_c_coords,to_zero_r_coords):
            
        #     ax.annotate(
        #         fov,
        #         xy=(x,y), xytext=(-0, 15),
        #         textcoords='offset points', ha='center', va='bottom',fontsize=12)
            
        #     ax.annotate(round_num, (x, y), color='white', weight='bold', 
        #                     fontsize=10, ha='center', va='center')
            
        #     ax.annotate('m' + str(match), xy=(x, y), xytext=(0, -10), color='white', weight='bold', 
        #                     textcoords='offset points', fontsize=5, ha='center', va='center')
            

        # ax.set_aspect('equal')
        # ax.axis('off')
        # plt.gca().invert_yaxis()
        # plt.tight_layout()

        # plt.savefig(self.experiment_fpath / 'output_figures' / 'registration_error.png',dpi=200,pad_inches=0)

    def run_qc(self):
        self.create_error_df()
        self.plot_error()


def QC_check_experiment_yaml_file(experiment_fpath:str):
    """
    This function is used to check that the parameter
    needed for processing are included in the experiment_name.yaml
    file and have the expected values
    
    Args:
        experiment_fpath: str
            str with the path to the folder of the experiment to process
    """
    experiment_fpath = Path(experiment_fpath)
    logger = selected_logger()
    
    try:
        experiment_info = load_experiment_config_file(experiment_fpath)
    except:
        logger.error(f'missing experiment info yaml file')
        sys.exit(f'missing experiment info yaml file')
    else:

        codebooks_list = (experiment_fpath.parent / 'codebooks').glob('*')
        probe_sets_list = (experiment_fpath.parent / 'probes_sets').glob('*')


        experiment_info_keys = list(experiment_info.keys())

        required_keys = ['Stitching_type',
                            'Experiment_type',
                            'Barcode_length',
                            'Barcode',
                            'Codebook',
                            'Machine',
                            'Operator',
                            'Overlapping_percentage',
                            'Probe_FASTA_name',
                            'Species',
                            'Start_date',
                            'Strain',
                            'Tissue',
                            'Pipeline']

        for key in required_keys:
            if key not in experiment_info_keys:
                logger.error(f'{key} field is missing in the experiment file')
                sys.exit(f'{key} field is missing in the experiment file')
        

        if 'serial' not in experiment_info['Experiment_type']:
            if not experiment_info['Codebook']:
                logger.error(f"The experiment type {experiment_info['Experiment_type']} needs a codebook")
                sys.exit(f"The experiment type {experiment_info['Experiment_type']} needs a codebook")
            else:
                present = [x.name for x in codebooks_list if experiment_info['Codebook'] == x.name]
                if not present:
                    logger.error(f'Specified codebook is missing from the database')
                    sys.exit(f'Specified codebook is missing from the database')

            if experiment_info['Barcode'] not in ['True', 'False']:
                logger.error(f'Value corresponding to Barcode keyword must be True or False ')
                sys.exit(f'Value corresponding to Barcode keyword must be True or False ')
        
            if experiment_info['Barcode'] == 'True':
                if 'Barcode_length' not in experiment_info_keys:
                    logger.error(f'Barcode_length keyword in the experiment file')
                    sys.exit(f'Barcode_length keyword in the experiment file')
                elif experiment_info['Barcode_length'] != 16:
                    logger.error(f'Wrong barcode length')
                    sys.exit(f'Barcode_length keyword in the experiment file')
            
        if experiment_info['Machine'] not in ['ROBOFISH1', 'ROBOFISH2', 'NOT-DEFINED']:
            logger.error(f'Wrong machine name')
            sys.exit(f'Wrong machine name')
        
        if experiment_info['Stitching_type'] not in ['small-beads', 'large-beads','both-beads', 'nuclei']:
            logger.error(f'Wrong Stitching_type selected in the experiment file')
            sys.exit(f'Wrong Stitching_type selected in the experiment file')
        
        if experiment_info['Experiment_type'] not in ['smfish-serial', 'smfish-barcoded', 'eel-barcoded']:
            logger.error(f'Wrong Experiment_type selected in the experiment file')
            sys.exit(f'Wrong Experiment_type selected in the experiment file')

        if not experiment_info['Probe_FASTA_name']:
            logger.error(f'Experiment require the probes name')
            sys.exit(f'Probes keyword in the experiment file')
        else:
            present = [x.name for x in probe_sets_list if experiment_info['Probe_FASTA_name'] == x.name]
            if not present:
                logger.error(f'Specified probes set is missing from the database')
                sys.exit(f'Specified probes set is missing from the database')


def QC_matching_nd2_metadata_robofish(all_raw_files:list):
    """
    This function is used to check that each of the nd2 files
    generated by the microscope has a matching pkl metadata
    file generated by robofish

    Args:
        all_raw_files: list
            list with all the paths of the nd2 files to process

    """

    logger = selected_logger()
    experiment_fpath = all_raw_files[0].parent
    all_info_files = list(experiment_fpath.glob('*.pkl'))
  
    if len(all_info_files) == 0:
        logger.error(f"no .pkl files in the folder")
        sys.exit(f"no .pkl files in the folder")

    # Determine if there are multiple metadata files with same number
    all_codes = []
    for meta_file_path in all_info_files:
        # collect the count code
        count_code = re.search(r'(Count)\d{5}', meta_file_path.stem)
        assert count_code, logger.error(f'{meta_file_path.stem} does not contain the CountXXXXX code')
        count_code = count_code.group()
        all_codes. append(count_code)
    
    if all_codes:
        all_codes_counts_dict = {i:all_codes.count(i) for i in all_codes}
        all_codes_counts_array = np.array(list(all_codes_counts_dict.values()))
        if np.any(all_codes_counts_array>1):
            for count_code,value in all_codes_counts_dict.items():
                if value >1:
                    logger.error(f' multiple pkl files with {count_code}')
                    sys.exit(f' multiple pkl files with {count_code}')
            
            logger.error(f'fix naming of the files with the repeated codes')
            sys.exit(f'fix naming of the files with the repeated codes')

    missing_pkl = []
    for nd2_file_path in all_raw_files:
        # collect the count code
        try:
            count_code = re.search(r'(Count)\d{5}', nd2_file_path.stem).group()
        except:
            count_code = None
            logger.error(f'{nd2_file_path.stem} does not contain the CountXXXXX code')
            sys.exit(f'{nd2_file_path.stem} does not contain the CountXXXXX code')
        else:
            try:
                info_file = [info_file for info_file in all_info_files if count_code  in info_file.stem][0]
            except IndexError:
                logger.error(f'{nd2_file_path.stem} does not have the corresponding pkl file')
                missing_pkl.append(nd2_file_path.stem)

    if missing_pkl:
        logger.error(f'collect the missing pkl for {missing_pkl} before parsing')
        sys.exit(f'collect the missing pkl for {missing_pkl} before parsing')