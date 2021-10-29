"""
Set of utility functions that do not belong to any
specific module.
"""


from typing import *
import shutil
import sys
import git
import os
import click
import ctypes
import numpy as np
from skimage import img_as_float64
from pathlib import Path
from datetime import datetime


from pysmFISH.logger_utils import selected_logger
from pysmFISH.configuration_files import create_general_analysis_config_file


def create_dir(dir_path: str):
    """Create a directory

    The directory is created only if it is not present in the folder.
    If the directory is already present it is not overwritten

    Args:
        dir_path (str): Create a directory and change it access to 777.
    """
    try:
        os.stat(dir_path)
    except:
        os.mkdir(dir_path)
        os.chmod(dir_path,0o777)


def get_git_hash() -> str:
    """Function used to get the hash of the current commit used for the 
    processing of the data

    Returns:
        version (str): hash of the commit used for processing the data
    """
    logger = selected_logger()
    repo = git.Repo(search_parent_directories=True)
    version = repo.head.object.hexsha
    logger.debug(f"current code version: {version}")
    return version


def create_processing_env(processing_folder_path:str):
    logger = selected_logger()
    
    processing_folder_path = Path(processing_folder_path)
    config_path = processing_folder_path / 'config_db'
    codebooks_path = processing_folder_path / 'codebooks'
    probes_path = processing_folder_path / 'probes_sets'
    
    # Create the directories that will contain the files required for the processing
    create_dir(config_path)
    create_dir(codebooks_path)
    create_dir(probes_path)

    # Create the analysis master files
    create_general_analysis_config_file(config_path)


def nice_deltastring(delta) -> str:
    """Function that provide a nice visualization of execution time

    From Sten Linnarsson lab
    https://github.com/linnarsson-lab/cytograph-shoji/blob/6389e8864c755f056ab7c9b51892650e5ed4f040/cytograph/pipeline/workflow.py#L12
    
    Args:
        delta (datetime.timedelta): execution time. Ex. End time - start time

    Returns:
        str: time difference in ms
    """

    result = []
    s = delta.total_seconds()
    h = s // 60 // 60
    if h >= 1:
        result.append(f"{int(h)}h")
        s -= h * 60 * 60
    m = s // 60
    if m >= 1:
        result.append(f"{int(m)}m")
        s -= m * 60
    if s >= 1:
        result.append(f"{int(s)}s")
    if len(result) > 0:
        return " ".join(result)
    else:
        return f"{delta.microseconds // 1000} ms"


def free_space(hd_path:str, min_free_space:int) -> bool:
    """Function used to determine if there is enough space in the
        HD where the experiment will be processed
    Args:
        hd_path (str): pathway of the target HD where the processing will be run
        min_free_space (int): minimum space required for processing in Gb (ex: 1000)

    Returns:
        bool: True if there is enough free space to process the experiment
    """

    logger = selected_logger()
    total, used, free = shutil.disk_usage(hd_path)
    free_space_giga = free // (2**30)
    if free_space_giga <= min_free_space:
        logger.info(f'Free space in the HD: {free_space_giga} Gb data cannot be transferred,\
                        not enough space on the HD')
        return False
    else:
        logger.info(f'Free space in the HD: {free_space_giga} Gb data can be transferred')
        return True


def end_processing_file(path_destination:str, completion_pattern:str='processing_completed.txt'):
    """Function used to create the file that signal that the
    analysis is completed

    Args:
        path_destination (str): Path where to save the file
        completion_pattern (str, optional): Tag indicating the completion. Defaults to 'processing_completed.txt'.
    """
    
    logger = selected_logger()
    
    fname = Path(path_destination) / completion_pattern
    try:
        open(fname, 'a').close()
    except OSError:
        logger.error('Failed to create the processing termination file')
    else:
        logger.info('processing termination files created')






def create_folder_structure(experiment_fpath:str,run_type:str):
    """Create the folder structure where to sort the files generated by ROBOFISH

    FOLDER STRUCTURE
    - original_robofish_logs: contains all the original robofish logs.
    - extra_files: contains the extra files acquired during imaging.
    - extra_processing_data: contains extra files used in the analysis 
            like the dark images for flat field correction.
    - pipeline_config: contains all the configuration files.
    - raw_data: contains the renamed .nd2 files and the corresponding 
            pickle configuration files.
    - output_figures: contains the reports and visualizations
    - probes: contains the fasta file with the probes used in the experiment

    Args:
        experiment_fpath (str): path to the experiment to process
        run_type (str): type of processing run. It can be:
                - original: if processing a new experiment
                - re-run: if it is a reprocessing
    """

    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)

    if run_type == 'new':

        folders_list = ['raw_data',
                        'original_robofish_logs',
                        'extra_processing_data',
                        'extra_files',
                        'pipeline_config',
                        'output_figures',
                        'probes',
                        'logs',
                        'results',
                        'microscope_tiles_coords',
                        'notebooks']
    else:
        folders_list = ['extra_processing_data',
                    'pipeline_config',
                    'output_figures',
                    'probes',
                    'logs',
                    'results',
                    'microscope_tiles_coords',
                    'notebooks']
    
    for folder_name in folders_list:
        try:
            os.stat(experiment_fpath / folder_name )
            logger.info(f'{folder_name} already exist')
        except FileNotFoundError:
            os.mkdir(experiment_fpath / folder_name)
            os.chmod(experiment_fpath / folder_name,0o777)


def collect_processing_files(experiment_fpath:str, experiment_info:Dict):
    """Gather codebooks and probes from the storage folders

    Args:
        experiment_fpath (str): path to the experiment to process
        experiment_info (Dict): content of the configuration file (experiment_name_config.yaml)
            present in the experiment folder
    Todo:
        Modify the selection and the loading of multiple codebooks / probes
    """
    logger = selected_logger()

    experiment_fpath = Path(experiment_fpath)

    try:
        machine = experiment_info['Machine']
    except NameError:
        machine = 'NOT_DEFINED'


    for idx, probes in experiment_info['Probes_FASTA'].items():
        if probes != 'None':
            probes_fpath = experiment_fpath.parent / 'probes_sets' / probes
            try:
                shutil.copy(probes_fpath, (experiment_fpath / 'probes'/ probes))
            except FileNotFoundError:
                logger.error('missing probes set file')
                sys.exit('missing probes set file')
        else:
            
            if 'barcoded' in experiment_info['Experiment_type']:
                codebooks_folder = experiment_fpath.parent / 'codebooks'
                
                for idx, codebook in experiment_info['Codebooks'].items():
                    if codebook != 'None':
                        codebook_fpath = codebooks_folder / codebook
                        
                        # Create codebook folder in the experiment folder
                        try:
                            os.stat(experiment_fpath / 'codebook' )
                            logger.info(f'codebook folder already exist')
                        except FileNotFoundError:
                            os.mkdir(experiment_fpath / 'codebook')
                            os.chmod(experiment_fpath / 'codebook',0o777)
                        try:
                            shutil.copy(codebook_fpath, (experiment_fpath / 'codebook'))
                        except FileNotFoundError:
                            logger.error('codebook is missing')
                            sys.exit('codebook is missing')


def sort_data_into_folders(experiment_fpath:str,experiment_info:Dict):
    """Sort the files created by ROBOFISH in subfolders

    Args:
        experiment_fpath (str): path to the experiment to process
        experiment_info (Dict): content of the configuration file (experiment_name_config.yaml)
            present in the experiment folder
    """
    
    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)

    # Move the robofish generated logs
    robofish_logs = list(experiment_fpath.glob('*.log'))
    if len(robofish_logs):
        for log in robofish_logs:
            shutil.move(log.as_posix(), (experiment_fpath / 'original_robofish_logs').as_posix())
            logger.debug(f'moved {log.stem} into original_robofish_logs')
    else:
        logger.debug(f'there are no original robofish logs in the experiment folder')

    # Move the crosses images
    crosses_files = list(experiment_fpath.glob('*_cross.nd2'))
    if len(crosses_files):
        for cross in crosses_files:
            shutil.move(cross.as_posix(), (experiment_fpath / 'extra_files').as_posix())
            logger.debug(f'moved {cross.stem} into extra_files')
    else:
        logger.debug(f'there are no crosses images in the experiment folders')
    
    # Move the coords files
    coords_files = list(experiment_fpath.glob('*.xml'))
    if len(coords_files):
        for coord in coords_files:
            shutil.move(coord.as_posix(), (experiment_fpath / 'extra_files').as_posix())
            logger.debug(f'moved {coord.stem} into extra_files')
    else:
        logger.debug(f'there are no coords files in the experiment folders')
    
    # Move the fresh nuclei file for eel
    if 'eel' in experiment_info['Experiment_type']:
        beads_files = list(experiment_fpath.glob('*ChannelEuropium_Cy3*'))
        nuclei_files = list(experiment_fpath.glob('*ChannelCy3_Nuclei_*'))
        
        if len(nuclei_files) or (len(beads_files)):
            try:
                os.stat(experiment_fpath / 'fresh_tissue' )
                logger.debug(f'fresh_tissue already exist')
            except FileNotFoundError:
                os.mkdir(experiment_fpath / 'fresh_tissue')
                os.chmod(experiment_fpath / 'fresh_tissue',0o777)

                os.mkdir(experiment_fpath / 'fresh_tissue' / 'results')
                os.chmod(experiment_fpath / 'fresh_tissue' / 'results',0o777)

                os.mkdir(experiment_fpath / 'fresh_tissue' / 'output_figures')
                os.chmod(experiment_fpath / 'fresh_tissue' / 'output_figures',0o777)
            
            for nuclei in nuclei_files:
                shutil.move(nuclei.as_posix(), (experiment_fpath / 'fresh_tissue').as_posix())
                logger.debug(f'moved {nuclei.stem} into fresh tissue')

            for beads in beads_files:
                shutil.move(beads.as_posix(), (experiment_fpath / 'fresh_tissue').as_posix())
                logger.debug(f'moved {beads.stem} into fresh tissue')
        
        large_views = list(experiment_fpath.glob('*ChannelDAPI*'))
        if len(large_views):
            for large_view in large_views:
                shutil.move(large_view.as_posix(), (experiment_fpath / 'extra_files').as_posix())
                logger.debug(f'moved {large_view.stem} into extra_files')
        else:
            logger.debug(f'The experiment does not have large view images of fresh nuclei')

    # move remaining .nd2 files
    all_nd2 = list(experiment_fpath.glob('*.nd2'))
    remaining = [el for el in all_nd2 if 'Count' not in el.stem]
    if len(remaining):
        for fremaining in remaining:
            shutil.move(fremaining.as_posix(), (experiment_fpath / 'extra_files').as_posix())
            logger.debug(f'moved {fremaining.stem} into extra_files')
    else:
        logger.debug(f'There are no extra files to be moved')




def copytree(src:str, dst:str, symlinks:bool=False, ignore:str=None):
    """ Function used to copy an entire directory tree into another directory

    Code copied from:
    https://stackoverflow.com/questions/1868714/how-do-i-copy-an-entire-directory-of-files-into-an-existing-directory-using-pyth

    Args:
        src (str): path of the directory to copy
        dst (str): path of the destination
        symlinks (bool, optional): Set to true if you want to copy symlinks as well. Defaults to False.
        ignore (str, optional): Files to ignore. Defaults to None.
    """
    if not os.path.exists(dst):
        os.makedirs(dst)
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            copytree(s, d, symlinks, ignore)
        else:
            if not os.path.exists(d) or os.stat(s).st_mtime - os.stat(d).st_mtime > 1:
                shutil.copy2(s, d)



def convert_to_uint16_full_float64_range(image: np.ndarray)->np.ndarray:
    """Utility function to convert images to full range float64
    Args:
        image (np.ndarray): Image as numpy array

    Returns:
        image (np.ndarray): Image as float64 numpy array
    """
    # Clip the values above 1
    image = image / np.finfo(np.float64).max
    # image[image > 1] = 1
    # in the current processing step the images are float
    # but not anymore between (0,1)
    # Scale to the max of the uint16
    image *= np.iinfo(np.uint16).max
    # Round and convert to integer
    image = np.uint16(np.rint(image))
    return image


def convert_from_uint16_to_float64(image: np.ndarray)->np.ndarray:
    """Utility function to convert images to float64
    Now a similar function is included in scikit-image but was not
    there when I started the developing

    Args:
        image (np.ndarray): Image as numpy array

    Returns:
        image (np.ndarray): Image as float64 numpy array
    """
    image = image / np.iinfo(np.uint16).max
    image = img_as_float64(image)
    return image


def convert_to_uint16(image:np.ndarray)-> np.ndarray:
    """ Utility function to convert images to uint16
    
    Consider that the original images were uint16
    that is why i normalize for np.iinfo(np.uint16).max*1.0
    otherwise the values will be to small

    Args:
        image (np.ndarray): Image as numpy array, probably float64 after
        processing

    Returns:
        image: (np.ndarray): Image as numpy array converted to uint16
    """
    
 
    image = image / (np.iinfo(np.uint16).max*1.0)
    image *= np.iinfo(np.uint16).max
    # Round and convert to integer
    image = image.astype(np.uint16)
    return image



class OptionEatAll(click.Option):
    """ Option class to ingest undefined number of options.

    Code copied from:
    https://stackoverflow.com/questions/48391777/nargs-equivalent-for-options-in-click

    """

    def __init__(self, *args, **kwargs):
        self.save_other_options = kwargs.pop('save_other_options', True)
        nargs = kwargs.pop('nargs', -1)
        assert nargs == -1, 'nargs, if set, must be -1 not {}'.format(nargs)
        super(OptionEatAll, self).__init__(*args, **kwargs)
        self._previous_parser_process = None
        self._eat_all_parser = None

    def add_to_parser(self, parser, ctx):

        def parser_process(value, state):
            # method to hook to the parser.process
            done = False
            value = [value]
            if self.save_other_options:
                # grab everything up to the next option
                while state.rargs and not done:
                    for prefix in self._eat_all_parser.prefixes:
                        if state.rargs[0].startswith(prefix):
                            done = True
                    if not done:
                        value.append(state.rargs.pop(0))
            else:
                # grab everything remaining
                value += state.rargs
                state.rargs[:] = []
            value = tuple(value)

            # call the actual process
            self._previous_parser_process(value, state)

        retval = super(OptionEatAll, self).add_to_parser(parser, ctx)
        for name in self.opts:
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break
        return retval


def trim_memory() -> int:
     libc = ctypes.CDLL("libc.so.6")
     return libc.malloc_trim(0)