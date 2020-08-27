from typing import *
import logging
import yaml
import sys
import os
import click
import datetime
import numpy as np
from collections import OrderedDict
from skimage import img_as_float64
from pathlib import Path













# to avoid reference for nested structures
# https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml (comment)
yaml.SafeDumper.ignore_aliases = lambda *args : True

def convert_to_uint16_full_float64_range(image):
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

def convert_to_uint16(image):
    """
    Consider that the original images were uint16
    that is why i normalize for np.iinfo(np.uint16).max*1.0
    otherwise the values will be to small
    """
 
    image = image / (np.iinfo(np.uint16).max*1.0)
    image *= np.iinfo(np.uint16).max
    # Round and convert to integer
    image = image.astype(np.uint16)
    return image

def convert_from_uint16_to_float64(image):
    image = image / np.iinfo(np.uint16).max
    image = img_as_float64(image)
    return image


def load_pipeline_config_file(pipeline_config_fpath):

    try:
        pipeline_config_dict = OrderedDict(yaml.safe_load(open(pipeline_config_fpath, 'rb')))
        return pipeline_config_dict
    except FileExistsError:
        logging.exception(f'{pipeline_config_fpath} missing')
        sys.exit(f'{pipeline_config_fpath} missing')
    except NameError:
        logging.exception(f'{pipeline_config_fpath} wrong pathway name')
        sys.exit(f'{pipeline_config_fpath} wrong pathway name')



def load_running_analysis_config_file(experiment_fpath):
    experiment_fpath = Path(experiment_fpath)
    config_fpath = experiment_fpath /  'pipeline_config'
    running_analysis_config_fpath = list(config_fpath.glob('*_analysis_run.yaml'))
    if len(running_analysis_config_fpath) == 0:
        logging.exception(f'missing running analysis config file')
        sys.exit(f'missing running analysis config file')
    elif len(running_analysis_config_fpath) > 1:
        logging.exception(f'too many running analysis config file')
        sys.exit(f'too many running analysis config file')
    else:
        running_analysis_config_fpath = running_analysis_config_fpath[0]
        running_analysis_config = load_pipeline_config_file(running_analysis_config_fpath)
        return running_analysis_config

def create_dir(dir_path):
    try:
        os.stat(dir_path)
    except:
        os.mkdir(dir_path)
        os.chmod(dir_path,0o777)



def modify_images_config(experiment_fpath:str,keyword:str):
    """
    Function used to batch modify the values of the
    processing parameters in the images_config.yaml
    files

    Parameters:
    -----------
    experiment_fpath: str
        path of the experiment to process
    keyword: str
        selection criteria used to process a subset of
        configuration files

    """
    experiment_fpath = Path(experiment_fpath)
    pipeline_config_dpath = experiment_fpath / 'pipeline_config'

    for config_fpath in pipeline_config_dpath.glob('*_images_config.yaml'):
        if keyword in config_fpath.name:
            img_config_dict = load_pipeline_config_file(config_fpath)
            for fov_name, fov_data in img_config_dict.items():
                for round_name, round_data in fov_data.items():
                    # round_data['analysis_parameters'] = {}
                    # round_data['analysis_parameters']['preprocessing'] = {}
                    # round_data['analysis_parameters']['preprocessing']['flat_field_kernel'] = (2,100,100)
                    # round_data['analysis_parameters']['preprocessing']['filtering_small_kernel'] = (1,8,8)
                    # round_data['analysis_parameters']['preprocessing']['filtering_laplacian_kernel'] = (0.2,0.5,0.5)
                    # round_data['analysis_parameters']['preprocessing']['large_obj_removal_percentile'] = 99
                    # round_data['analysis_parameters']['preprocessing']['large_obj_removal_min_obj_size']= 50
                    # round_data['analysis_parameters']['preprocessing']['large_obj_removal_selem'] = 3
                    round_data['analysis_parameters']['counting'] = {}
                    round_data['analysis_parameters']['counting']['min_distance'] = 4
                    round_data['analysis_parameters']['counting']['min_obj_size'] = 4
                    round_data['analysis_parameters']['counting']['max_obj_size'] = 200
                    round_data['analysis_parameters']['beads_registration'] = {}

            with open(config_fpath, 'w') as new_config:
                    yaml.safe_dump(dict(img_config_dict), new_config,default_flow_style=False,sort_keys=False)



def create_stringency_selection_fov_yaml(pipeline_config_fpath):
    
    pipeline_config_dict = load_pipeline_config_file(pipeline_config_fpath)
    fovs = (pipeline_config_dict['channel_coords_dict'].keys())
    stringency = pipeline_config_dict['pipeline_required_parameters']['counting_settings']['stringency']
    pipeline_config_dict['pipeline_required_parameters']['counting_settings']['stringency'] = {}
    pipeline_config_dict['pipeline_required_parameters']['counting_settings']['stringency'] = {k:int(stringency) for k in fovs}
    with open(pipeline_config_fpath, 'w') as new_config:
            yaml.safe_dump(dict(pipeline_config_dict), new_config,default_flow_style=False,sort_keys=False)

def create_specific_fov_counting_thr_yaml(pipeline_config_fpath):
    pipeline_config_dict = load_pipeline_config_file(pipeline_config_fpath)
    fovs = (pipeline_config_dict['channel_coords_dict'].keys())
    stringency = pipeline_config_dict['pipeline_required_parameters']['counting_settings']['thr_fov'] = {}
    pipeline_config_dict['pipeline_required_parameters']['counting_settings']['thr_fov'] = {k:None for k in fovs}
    with open(pipeline_config_fpath, 'w') as new_config:
            yaml.safe_dump(dict(pipeline_config_dict), new_config,default_flow_style=False,sort_keys=False)

def select_dir(processing_directory,keyword):
    list_hyb_fpaths = []
    processed_dirs = next(os.walk(processing_directory))[1]
    
    if keyword:
        for dir_name in processed_dirs:
            if keyword in dir_name:
                dir_path = processing_directory / dir_name
                list_hyb_fpaths.append(dir_path)
    else:
        for dir_name in processed_dirs:
            dir_path = processing_directory / dir_name
            list_hyb_fpaths.append(dir_path)

    return list_hyb_fpaths



class OptionEatAll(click.Option):
    '''
    https://stackoverflow.com/questions/48391777/nargs-equivalent-for-options-in-click

    Option class to ingest undefined number of options.
    '''

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