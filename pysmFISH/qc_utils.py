"""
set of useful function to determine if different process steps
that terminated without error generated the correct output.
Ex. missing positions after parsing the data
"""

from typing import *
from prefect import task
from prefect.engine import signals
import re
import sys
import numpy as np
from pathlib import Path

from dask import dataframe as dd

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatch

from pysmFISH.logger_utils import json_logger
from pysmFISH.configuration_files import load_experiment_config_file


class QC_registration_error():

    def __init__(self, client, experiment_fpath, analysis_parameters, tiles_coords, img_width, img_height):
        self.client = client
        self.experiment_fpath = Path(experiment_fpath)
        self.analysis_parameters = analysis_parameters
        self.tiles_coords = tiles_coords
        self.img_width = img_width
        self.img_height = img_height

    def create_error_df(self):
        all_counts_folder = self.experiment_fpath / 'tmp' / 'registered_counts'
        search_key = '*decoded*'
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
            self.error_output_df = pd.concat([error_output_df,matching_rounds.iloc[0]],axis=1)
        self.error_output_df = error_output_df.T
        # #     RegistrationMinMatchingBeads = analysis_parameters['RegistrationMinMatchingBeads']
        #     registration_error_df.columns = ["_".join(x) for x in registration_error_df.columns.ravel()]

    def plot_error(self):
    
        scale_value = 5
        self.tiles_coords = self.tiles_coords / scale_value
        
        fovs = self.error_output_df['fov_num'].values
        rounds_num = self.error_output_df['round_num'].values
        min_errors = self.error_output_df['min_number_matching_dots_registration'].values
        
        
        plt.ioff()
        fig = plt.figure(figsize=(30,20))
        ax = fig.add_subplot(111)

        r_coords = self.tiles_coords[:,0]
        c_coords = self.tiles_coords[:,1]
        r_coords_min = r_coords.min()
        c_coords_min = c_coords.min()
        to_zero_r_coords = r_coords - r_coords_min
        to_zero_c_coords = c_coords - c_coords_min

        ax.plot(to_zero_c_coords,to_zero_r_coords,'xk')


        width = self.img_width / scale_value
        height = self.img_height / scale_value


        # create colormap for error in the registration
        errors_normalized = (min_errors -min(min_errors)) / (max(min_errors -min(min_errors)))
        threshold = (RegistrationMinMatchingBeads-min(min_errors)) / (max(min_errors -min(min_errors)))
        nodes = [0,threshold, threshold, 1.0]

        colors = ["black", "black", "blue", "magenta"]
        cmap = LinearSegmentedColormap.from_list("", list(zip(nodes, colors)))
        cmap.set_under("black")
        
        
        rectangles = {}
        # to change after I corrected the error with the '10'
        for idx,fov in enumerate(fovs):
        
            round_num = rounds_num[idx]
            min_error = min_errors[idx]
            
            tag = 'fov_' + str(fov) + '\n' + 'match: '+ str(min_error) + '\n' + 'round: ' + str(round_num)
            rectangles.update({tag : mpatch.Rectangle( xy=(to_zero_c_coords[idx],to_zero_r_coords[idx] ),  # point of origin.
                                            width=width,
                                            height=height,
                                            linewidth=0.5,
                                            edgecolor='k',
                                            fill=True,
                                            alpha = 0.8,
                                            facecolor=cmap(errors_normalized[idx]))})

        for r in rectangles:
            ax.add_artist(rectangles[r])
            rx, ry = rectangles[r].get_xy()
            cx = rx + rectangles[r].get_width()/2.0
            cy = ry + rectangles[r].get_height()/2.0

            ax.annotate(r, (cx, cy), color='white', weight='bold', 
                        fontsize=10, ha='center', va='center')


        ax.set_xlim((-width, to_zero_c_coords.max()+ width + (0.05*width)))
        ax.set_ylim((-height, to_zero_r_coords.max()+ height + (0.05*height)))
        ax.set_aspect('equal')
        ax.axis('off')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(experiment_fpath / 'output_figures' / 'registration_error.png',dpi=200,pad_inches=0)
        plt.show()

    def run_qc(self):
        self.create_error_df()
        self.plot_error()


def check_experiment_yaml_file(experiment_fpath:str):
    """
    This function is used to check that the parameter
    needed for processing are included in the experiment_name.yaml
    file and have the expected values
    
    Args:
        experiment_fpath: str
            str with the path to the folder of the experiment to process
    """
    experiment_fpath = Path(experiment_fpath)
    logger = json_logger((experiment_fpath / 'logs'),'qc_experiment_yaml')
    
    try:
        experiment_info = load_experiment_config_file(experiment_fpath)
    except:
        logger.error(f'missing experiment info yaml file')
        sys.exit(f'missing experiment info yaml file')
    else:

        codebooks_list = (experiment_fpath.parent / 'codebooks').glob('*')
        probe_sets_list = (experiment_fpath.parent / 'probes_sets').glob('*')


        experiment_info_keys = list(experiment_info.keys())

        if 'Codebook' not in experiment_info_keys:
            logger.error(f'Codebook keyword in the experiment file')

        if 'Barcode' not in experiment_info_keys:
            logger.error(f'Barcode keyword in the experiment file')
        elif type(experiment_info['Barcode']) != bool:
            logger.error(f'Value corresponding to Barcode keyword must be bool ')
            sys.exit(f'Value corresponding to Barcode keyword must be bool ')
        
        if type(experiment_info['Barcode']):
            if 'Barcode_length' not in experiment_info_keys:
                logger.error(f'Barcode_length keyword in the experiment file')
                sys.exit(f'Barcode_length keyword in the experiment file')
            elif experiment_info['Barcode_length'] != 16:
                logger.error(f'Wrong barcode length')
                sys.exit(f'Barcode_length keyword in the experiment file')
            
            if 'Codebook' not in experiment_info_keys:
                logger.error(f'Codebook keyword in the experiment file')
                sys.exit(f'Codebook keyword in the experiment file')
            else:
                present = [x.name for x in codebooks_list if experiment_info['Codebook'] == x.name]
                if not present:
                    logger.error(f'Specified codebook is missing from the database')
                    sys.exit(f'Specified codebook is missing from the database')

        if 'Machine' not in experiment_info_keys:
            logger.error(f'Machine keyword in the experiment file')
        elif experiment_info['Machine'] not in ['ROBOFISH1', 'ROBOFISH2', 'NOT-DEFINED']:
            logger.error(f'Wrong machine name')
            sys.exit(f'Wrong machine name')
        

        if 'Overlapping_percentage' not in experiment_info_keys:
            logger.error(f'Overlapping_percentage keyword in the experiment file')
            sys.exit(f'Overlapping_percentage keyword in the experiment file')

        if 'Species' not in experiment_info_keys:
            logger.error(f'Species keyword in the experiment file')
        elif experiment_info['Species'] not in ['Mus Musculus', 'Homo Sapiens']:
            logger.error(f'Unknown Species selected')

        if 'roi' not in experiment_info_keys:
            logger.error(f'roi keyword in the experiment file')

        if 'StitchingChannel' not in experiment_info_keys:
            logger.error(f'StitchingChannel keyword in the experiment file')
            sys.exit(f'StitchingChannel keyword in the experiment file')

        if 'Stitching_type' not in experiment_info_keys:
            logger.error(f'Stitching_type keyword in the experiment file')
            sys.exit(f'Stitching_type keyword in the experiment file')
        elif experiment_info['Stitching_type'] not in ['small-beads', 'large-beads','both-beads', 'nuclei',]:
            logger.error(f'Wrong Stitching_type selected in the experiment file')
            sys.exit(f'Wrong Stitching_type selected in the experiment file')

        if 'EXP_name' not in experiment_info_keys:
            logger.error(f'EXP_name keyword in the experiment file')
            sys.exit(f'EXP_name keyword in the experiment file')
        
        if 'Experiment_type' not in experiment_info_keys:
            logger.error(f'Experiment_type keyword in the experiment file')
            sys.exit(f'Experiment_type keyword in the experiment file')
        elif experiment_info['Experiment_type'] not in ['smfish-serial', 'smfish-barcoded', 'eel-barcoded']:
            logger.error(f'Wrong Experiment_type selected in the experiment file')
            sys.exit(f'Wrong Experiment_type selected in the experiment file')

        if 'Probe_FASTA_name' not in experiment_info_keys:
            logger.error(f'Probes keyword in the experiment file')
            sys.error(f'Probes keyword in the experiment file')
        else:
            present = [x.name for x in probe_sets_list if experiment_info['Probe_FASTA_name'] == x.name]
            if not present:
                logger.error(f'Specified probes set is missing from the database')
                sys.exit(f'Specified probes set is missing from the database')


def qc_matching_nd2_metadata_robofish(all_raw_files:list):
    """
    This function is used to check that each of the nd2 files
    generated by the microscope has a matching pkl metadata
    file generated by robofish

    Args:
        all_raw_files: list
            list with all the paths of the nd2 files to process

    """

    logger = json_logger((experiment_fpath / 'logs'),'qc_matching_pkl')
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