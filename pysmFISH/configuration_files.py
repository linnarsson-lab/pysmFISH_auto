"""
create configuration files needed to run the pipeline

"""
from typing import *
import yaml
import sys
from pathlib import Path
from collections import OrderedDict

from prefect import task
from prefect import Task
from prefect.engine import signals

from pysmFISH.logger_utils import selected_logger


# to avoid reference for nested structures
# https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml (comment)
yaml.SafeDumper.ignore_aliases = lambda *args : True



def create_processing_env_config_file(experiment_fpath:str):
    """
    This function creates a configuration files with the parameters requested for the 
    setup of the processing cluster. 
    The current file is defaulted to htcondor

    Args:
        experiment_fpath: str
            path to the experiment that will be processed
    """

    processing_env_config_fpath = experiment_fpath / 'pipeline_config' / 'processing_env_config.yaml'
    processing_env_config = OrderedDict()

    processing_env_config['htcondor'] = {}

    processing_env_config['htcondor']['smfish-serial'] = {}
    processing_env_config['htcondor']['smfish-serial']['cores'] = 1
    processing_env_config['htcondor']['smfish-serial']['memory'] = '10GB'
    processing_env_config['htcondor']['smfish-serial']['disk'] = '0.1GB'
    processing_env_config['htcondor']['smfish-serial']['local_directory'] = '/tmp'

    processing_env_config['htcondor']['smfish-barcoded'] = {}
    processing_env_config['htcondor']['smfish-barcoded']['cores'] = 1
    processing_env_config['htcondor']['smfish-barcoded']['memory'] = '10GB'
    processing_env_config['htcondor']['smfish-barcoded']['disk'] = '0.1GB'
    processing_env_config['htcondor']['smfish-barcoded']['local_directory'] = '/tmp'

    processing_env_config['htcondor']['eel-barcoded'] = {}
    processing_env_config['htcondor']['eel-barcoded']['cores'] = 1
    processing_env_config['htcondor']['eel-barcoded']['memory'] = '10GB'
    processing_env_config['htcondor']['eel-barcoded']['disk'] = '0.1GB'
    processing_env_config['htcondor']['eel-barcoded']['local_directory'] = '/tmp'

    with open(processing_env_config_fpath, 'w') as new_config:
            yaml.safe_dump(dict(processing_env_config), new_config,default_flow_style=False,sort_keys=False)




def create_general_analysis_config_file(experiment_fpath:str):
    """
    This function creates a basic standard configuration files with the parameters used for running
    all available analysis. It will be stored in the config_db folder. The data required for the analysis
    will be extracted from the files using the experiment_info file. 
   
    Args:
        experiment_fpath: str
            path to the experiment that will be processed

    """
    
    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)
    analysis_config_fpath = experiment_fpath.parent / 'config_db' / 'analysis_config.yaml'
    analysis_parameters = OrderedDict()

    analysis_parameters['eel-barcoded'] = {}
    
    analysis_parameters['eel-barcoded']['ROBOFISH1'] = {}
    analysis_parameters['eel-barcoded']['ROBOFISH1']['fish'] = {}
    analysis_parameters['eel-barcoded']['ROBOFISH1']['small-beads'] = {}
    analysis_parameters['eel-barcoded']['ROBOFISH1']['large-beads'] = {}
    analysis_parameters['eel-barcoded']['ROBOFISH1']['staining'] = {}
    analysis_parameters['eel-barcoded']['ROBOFISH1']['fresh-nuclei'] = {}
    analysis_parameters['eel-barcoded']['ROBOFISH1']['BarcodesExtractionResolution'] = 3         
    analysis_parameters['eel-barcoded']['ROBOFISH1']['RegistrationReferenceHybridization'] = 1

    analysis_parameters['eel-barcoded']['ROBOFISH1']['fish']['PreprocessingFishFlatFieldKernel'] = (100,100)
    analysis_parameters['eel-barcoded']['ROBOFISH1']['fish']['PreprocessingFishFilteringSmallKernel'] = (8,8)
    analysis_parameters['eel-barcoded']['ROBOFISH1']['fish']['PreprocessingFishFilteringLaplacianKernel'] = (0.5,0.5)
    analysis_parameters['eel-barcoded']['ROBOFISH1']['fish']['CountingFishMinObjDistance'] = 2
    analysis_parameters['eel-barcoded']['ROBOFISH1']['fish']['CountingFishMaxObjSize'] = 200
    analysis_parameters['eel-barcoded']['ROBOFISH1']['fish']['CountingFishMinObjSize'] = 2
    analysis_parameters['eel-barcoded']['ROBOFISH1']['fish']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['ROBOFISH1']['fish']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['ROBOFISH1']['fish']['LargeObjRemovalPercentile'] = 99
    analysis_parameters['eel-barcoded']['ROBOFISH1']['fish']['LargeObjRemovalMinObjSize'] = 50
    analysis_parameters['eel-barcoded']['ROBOFISH1']['fish']['LargeObjRemovalSelem'] = 3

    analysis_parameters['eel-barcoded']['ROBOFISH1']['small-beads']['PreprocessingFishFlatFieldKernel'] = (100,100)
    analysis_parameters['eel-barcoded']['ROBOFISH1']['small-beads']['PreprocessingFishFilteringSmallKernel'] = (8,8)
    analysis_parameters['eel-barcoded']['ROBOFISH1']['small-beads']['PreprocessingFishFilteringLaplacianKernel'] = (0.5,0.5)
    analysis_parameters['eel-barcoded']['ROBOFISH1']['small-beads']['CountingFishMinObjDistance'] = 2
    analysis_parameters['eel-barcoded']['ROBOFISH1']['small-beads']['CountingFishMaxObjSize'] = 200
    analysis_parameters['eel-barcoded']['ROBOFISH1']['small-beads']['CountingFishMinObjSize'] = 2
    analysis_parameters['eel-barcoded']['ROBOFISH1']['small-beads']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['ROBOFISH1']['small-beads']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['ROBOFISH1']['small-beads']['LargeObjRemovalPercentile'] = 99
    analysis_parameters['eel-barcoded']['ROBOFISH1']['small-beads']['LargeObjRemovalMinObjSize'] = 50
    analysis_parameters['eel-barcoded']['ROBOFISH1']['small-beads']['LargeObjRemovalSelem'] = 3

    analysis_parameters['eel-barcoded']['ROBOFISH1']['large-beads']['PreprocessingFishFlatFieldKernel'] = (100,100)
    analysis_parameters['eel-barcoded']['ROBOFISH1']['large-beads']['PreprocessingFishFilteringSmallKernel'] = (8,8)
    analysis_parameters['eel-barcoded']['ROBOFISH1']['large-beads']['PreprocessingFishFilteringLaplacianKernel'] = (0.5,0.5)
    analysis_parameters['eel-barcoded']['ROBOFISH1']['large-beads']['CountingFishMinObjDistance'] = 2
    analysis_parameters['eel-barcoded']['ROBOFISH1']['large-beads']['CountingFishMaxObjSize'] = 200
    analysis_parameters['eel-barcoded']['ROBOFISH1']['large-beads']['CountingFishMinObjSize'] = 2
    analysis_parameters['eel-barcoded']['ROBOFISH1']['large-beads']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['ROBOFISH1']['large-beads']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['ROBOFISH1']['large-beads']['LargeObjRemovalPercentile'] = 99
    analysis_parameters['eel-barcoded']['ROBOFISH1']['large-beads']['LargeObjRemovalMinObjSize'] = 50
    analysis_parameters['eel-barcoded']['ROBOFISH1']['large-beads']['LargeObjRemovalSelem'] = 3

    analysis_parameters['eel-barcoded']['ROBOFISH1']['staining']['PreprocessingStainingFlatFieldKernel'] = (2,100, 100)
    
    analysis_parameters['eel-barcoded']['ROBOFISH1']['fresh-nuclei']['PreprocessingFreshNucleiLargeKernelSize'] =  (5,50,50)

    analysis_parameters['eel-barcoded']['ROBOFISH2'] = {}
    analysis_parameters['eel-barcoded']['ROBOFISH2']['fish'] = {}
    analysis_parameters['eel-barcoded']['ROBOFISH2']['small-beads'] = {}
    analysis_parameters['eel-barcoded']['ROBOFISH2']['large-beads'] = {}
    analysis_parameters['eel-barcoded']['ROBOFISH2']['staining'] = {}
    analysis_parameters['eel-barcoded']['ROBOFISH2']['fresh-nuclei'] = {}
    analysis_parameters['eel-barcoded']['ROBOFISH2']['BarcodesExtractionResolution'] = 3         
    analysis_parameters['eel-barcoded']['ROBOFISH2']['RegistrationReferenceHybridization'] = 1

    analysis_parameters['eel-barcoded']['ROBOFISH2']['fish']['PreprocessingFishFlatFieldKernel'] = (100,100)
    analysis_parameters['eel-barcoded']['ROBOFISH2']['fish']['PreprocessingFishFilteringSmallKernel'] = (8,8)
    analysis_parameters['eel-barcoded']['ROBOFISH2']['fish']['PreprocessingFishFilteringLaplacianKernel'] = (0.5,0.5)
    analysis_parameters['eel-barcoded']['ROBOFISH2']['fish']['CountingFishMinObjDistance'] = 2
    analysis_parameters['eel-barcoded']['ROBOFISH2']['fish']['CountingFishMaxObjSize'] = 200
    analysis_parameters['eel-barcoded']['ROBOFISH2']['fish']['CountingFishMinObjSize'] = 2
    analysis_parameters['eel-barcoded']['ROBOFISH2']['fish']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['ROBOFISH2']['fish']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['ROBOFISH2']['fish']['LargeObjRemovalPercentile'] = 99
    analysis_parameters['eel-barcoded']['ROBOFISH2']['fish']['LargeObjRemovalMinObjSize'] = 50
    analysis_parameters['eel-barcoded']['ROBOFISH2']['fish']['LargeObjRemovalSelem'] = 3

    analysis_parameters['eel-barcoded']['ROBOFISH2']['small-beads']['PreprocessingFishFlatFieldKernel'] = (100,100)
    analysis_parameters['eel-barcoded']['ROBOFISH2']['small-beads']['PreprocessingFishFilteringSmallKernel'] = (8,8)
    analysis_parameters['eel-barcoded']['ROBOFISH2']['small-beads']['PreprocessingFishFilteringLaplacianKernel'] = (0.5,0.5)
    analysis_parameters['eel-barcoded']['ROBOFISH2']['small-beads']['CountingFishMinObjDistance'] = 2
    analysis_parameters['eel-barcoded']['ROBOFISH2']['small-beads']['CountingFishMaxObjSize'] = 200
    analysis_parameters['eel-barcoded']['ROBOFISH2']['small-beads']['CountingFishMinObjSize'] = 2
    analysis_parameters['eel-barcoded']['ROBOFISH2']['small-beads']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['ROBOFISH2']['small-beads']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['ROBOFISH2']['small-beads']['LargeObjRemovalPercentile'] = 99
    analysis_parameters['eel-barcoded']['ROBOFISH2']['small-beads']['LargeObjRemovalMinObjSize'] = 50
    analysis_parameters['eel-barcoded']['ROBOFISH2']['small-beads']['LargeObjRemovalSelem'] = 3

    analysis_parameters['eel-barcoded']['ROBOFISH2']['large-beads']['PreprocessingFishFlatFieldKernel'] = (100,100)
    analysis_parameters['eel-barcoded']['ROBOFISH2']['large-beads']['PreprocessingFishFilteringSmallKernel'] = (8,8)
    analysis_parameters['eel-barcoded']['ROBOFISH2']['large-beads']['PreprocessingFishFilteringLaplacianKernel'] = (0.5,0.5)
    analysis_parameters['eel-barcoded']['ROBOFISH2']['large-beads']['CountingFishMinObjDistance'] = 2
    analysis_parameters['eel-barcoded']['ROBOFISH2']['large-beads']['CountingFishMaxObjSize'] = 200
    analysis_parameters['eel-barcoded']['ROBOFISH2']['large-beads']['CountingFishMinObjSize'] = 2
    analysis_parameters['eel-barcoded']['ROBOFISH2']['large-beads']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['ROBOFISH2']['large-beads']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['ROBOFISH2']['large-beads']['LargeObjRemovalPercentile'] = 99
    analysis_parameters['eel-barcoded']['ROBOFISH2']['large-beads']['LargeObjRemovalMinObjSize'] = 50
    analysis_parameters['eel-barcoded']['ROBOFISH2']['large-beads']['LargeObjRemovalSelem'] = 3

    analysis_parameters['eel-barcoded']['ROBOFISH2']['staining']['PreprocessingStainingFlatFieldKernel'] = (2,100,100)
    
    analysis_parameters['eel-barcoded']['ROBOFISH2']['fresh-nuclei']['PreprocessingFreshNucleiLargeKernelSize'] =  (5,50,50)


    analysis_parameters['eel-barcoded']['NOT_DEFINED'] = {}
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['fish'] = {}
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['small-beads'] = {}
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['large-beads'] = {}
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['staining'] = {}
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['fresh-nuclei'] = {}
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['BarcodesExtractionResolution'] = 3         
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['RegistrationReferenceHybridization'] = 1

    analysis_parameters['eel-barcoded']['NOT_DEFINED']['fish']['PreprocessingFishFlatFieldKernel'] = (100,100)
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['fish']['PreprocessingFishFilteringSmallKernel'] = (8,8)
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['fish']['PreprocessingFishFilteringLaplacianKernel'] = (0.5,0.5)
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['fish']['CountingFishMinObjDistance'] = 2
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['fish']['CountingFishMaxObjSize'] = 200
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['fish']['CountingFishMinObjSize'] = 2
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['fish']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['fish']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['fish']['LargeObjRemovalPercentile'] = 99
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['fish']['LargeObjRemovalMinObjSize'] = 50
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['fish']['LargeObjRemovalSelem'] = 3

    analysis_parameters['eel-barcoded']['NOT_DEFINED']['small-beads']['PreprocessingFishFlatFieldKernel'] = (100,100)
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['small-beads']['PreprocessingFishFilteringSmallKernel'] = (8,8)
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['small-beads']['PreprocessingFishFilteringLaplacianKernel'] = (0.5,0.5)
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['small-beads']['CountingFishMinObjDistance'] = 2
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['small-beads']['CountingFishMaxObjSize'] = 200
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['small-beads']['CountingFishMinObjSize'] = 2
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['small-beads']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['small-beads']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['small-beads']['LargeObjRemovalPercentile'] = 99
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['small-beads']['LargeObjRemovalMinObjSize'] = 50
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['small-beads']['LargeObjRemovalSelem'] = 3

    analysis_parameters['eel-barcoded']['NOT_DEFINED']['large-beads']['PreprocessingFishFlatFieldKernel'] = (100,100)
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['large-beads']['PreprocessingFishFilteringSmallKernel'] = (8,8)
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['large-beads']['PreprocessingFishFilteringLaplacianKernel'] = (0.5,0.5)
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['large-beads']['CountingFishMinObjDistance'] = 2
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['large-beads']['CountingFishMaxObjSize'] = 200
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['large-beads']['CountingFishMinObjSize'] = 2
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['large-beads']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['large-beads']['CountingFishNumPeaksPerLabel'] = 1
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['large-beads']['LargeObjRemovalPercentile'] = 99
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['large-beads']['LargeObjRemovalMinObjSize'] = 50
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['large-beads']['LargeObjRemovalSelem'] = 3

    analysis_parameters['eel-barcoded']['NOT_DEFINED']['staining']['PreprocessingStainingFlatFieldKernel'] = (2,100,100)
    
    analysis_parameters['eel-barcoded']['NOT_DEFINED']['fresh-nuclei']['PreprocessingFreshNucleiLargeKernelSize'] =  (5,50,50)

    try:
        with open(analysis_config_fpath, 'w') as new_config:
                yaml.safe_dump(dict(analysis_parameters), new_config,default_flow_style=False,sort_keys=False)
    except:
        logger.error(f'cannot save the analysis_config_file')


def create_specific_analysis_config_file(experiment_fpath:str, experiment_info:Dict):
    """
    Select the analysis parameters according to the processing machine. If the
    machine is not defined in the experiment_info dictionary a generic processing set up is defined

    Args:
        experiment_fpath: str
            path to the experiment that will be processed

        experiment_info: dict
            dictionary containing all the info generated by robofish
    """

    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)
    general_analysis_config_fpath = experiment_fpath.parent / 'config_db' / 'analysis_config.yaml'
    analysis_config_fpath = experiment_fpath / 'pipeline_config' / 'analysis_config.yaml'
    analysis_config = OrderedDict()

    general_analysis_config = yaml.safe_load(open(general_analysis_config_fpath, 'rb'))

    try:
        machine = experiment_info['Machine']
    except NameError:
        machine = 'NOT_DEFINED'

    try:
        experiment_type = experiment_info['Experiment_type']
        analysis_config = general_analysis_config[experiment_type][machine]
    except:
        logger.error(f'Unidentified experiment type in the config.yaml file')
    else:
        try:
            with open(analysis_config_fpath, 'w') as new_config:
                yaml.safe_dump(dict(analysis_config), new_config,default_flow_style=False,sort_keys=False)
        except:
            logger.error(f'cannot save the analysis_config_file')


def load_experiment_config_file(experiment_fpath:str):
    """
    Function that load the experiment general information generated by the
    machines.

    Args:
        experiment_fpath: str
            location of the folder to be processed
    Return:
        experiment_info: ordered dict
            ordered dict with the parsed info generated by the instrument

    """
    logger = selected_logger()  
    experiment_fpath = Path(experiment_fpath)
    experiment_name = experiment_fpath.stem
    experiment_name = experiment_name.split('_auto')[0]
    search_key = experiment_name + '_config.yaml'
    
    try:
        experiment_info_fpath = list(experiment_fpath.glob(search_key))[0]
    except:
        logger.error(f'No experiment info file in {experiment_fpath}')
        sys.exit(f'No experiment info file in {experiment_fpath}')
    else:
        try:
            experiment_info = yaml.safe_load(open(experiment_info_fpath, 'rb'))
            return experiment_info
        except:
            logger.error(f'Experiment info file has the wrong name in {experiment_fpath}')
            sys.exit(f'Experiment info file has the wrong name in {experiment_fpath}')


def load_processing_env_config_file(config_db_fpath:str):
    """
    Function used to load the parameters used for setting up 
    the processing cluster

    Args:
        config_db_fpath; str
            path to the folder containing the data_transfer_config.yaml
    """
    # logger = prefect_logging_setup('load-processing-env-config')
    config_db_fpath = Path(config_db_fpath)
    processing_env_config_fpath = config_db_fpath / 'processing_env_config.yaml'
    try:
        processing_env_config = OrderedDict(yaml.safe_load(open(processing_env_config_fpath, 'rb')))
    except (FileExistsError,NameError,FileNotFoundError) as e:
        # logger.debug(f'{processing_env_config_fpath} missing')
        try:
            create_processing_env_config_file(experiment_fpath)
        except:
            # logger.error('cannot create the processing config file')
            # signals.FAIL('cannot create the processing config file')
            print('cane')
        else:
            processing_env_config = OrderedDict(yaml.safe_load(open(processing_env_config_fpath, 'rb')))
            return processing_env_config
    else:
        return processing_env_config


def load_analysis_config_file(experiment_fpath:str):
    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)
    analysis_config_fpath = experiment_fpath / 'pipeline_config' / 'analysis_config.yaml'
    try:
        analysis_config = OrderedDict(yaml.safe_load(open(analysis_config_fpath, 'rb')))
    except (FileExistsError,NameError,FileNotFoundError) as e:
        logger.debug(f'{analysis_config_fpath} missing')
        signals.FAIL(f'{analysis_config_fpath} missing')
    else:
        return analysis_config


# def load_transferring_config(config_db_fpath:str)->Dict:
#     """
#     Function used to load the parameters used for transfering
#     the data from the storage to the processing HD

#     Args:
#         config_db_fpath; str
#             path to the folder containing the data_transfer_config.yaml
#     """
    
#     logger = prefect_logging_setup(f'load-transfer-config')
#     config_db_fpath = Path(config_db_fpath)

#     transfer_config_fpath = config_db_fpath / 'data_transfer_config.yaml'
#     try:
#          transfer_config = OrderedDict(yaml.safe_load(open(transfer_config_fpath, 'rb')))
#     except (FileExistsError,NameError,FileNotFoundError) as e:
#         logger.debug(f'{transfer_config_fpath} missing')
#         signals.FAIL(f'{transfer_config_fpath} missing')
#     else:
#         return  transfer_config

class load_transferring_config(Task):
    """
    Class used to load the parameters required for transfering
    the data from the storage to the processing HD

    Args:
        config_db_fpath; str
            path to the folder containing the data_transfer_config.yaml
    """
    def run(self, config_db_fpath:str)->Dict:
        """
        Function used to load the parameters used for transfering
        the data from the storage to the processing HD

        Args:
        -----
            config_db_fpath; str
                path to the folder containing the data_transfer_config.yaml

        Returns:
        --------
            transfer_config: Dict
                dictionary with all the info necessary for transferring the data

        """
        
        config_db_fpath = Path(config_db_fpath)

        transfer_config_fpath = config_db_fpath / 'data_transfer_config.yaml'
        try:
            transfer_config = OrderedDict(yaml.safe_load(open(transfer_config_fpath, 'rb')))
        except (FileExistsError,NameError,FileNotFoundError) as e:
            self.logger.error(f'{transfer_config_fpath} missing')
            signals.FAIL(f'{transfer_config_fpath} missing')
        else:
            return  transfer_config




# @task(name='load_experiments_info_file')
# def load_experiment_config_file(experiment_fpath:str):
#     """
#     Function that load the experiment general information generated by the
#     machines.

#     Args:
#         experiment_fpath: str
#             location of the folder to be processed
#     Return:
#         experiment_info: ordered dict
#             ordered dict with the parsed info generated by the instrument

#     """
    
#     logger = prefect_logging_setup('load_experiments_info_file')

#     experiment_fpath = Path(experiment_fpath)
#     experiment_name = experiment_fpath.stem
#     experiment_name = experiment_name.split('_auto')[0]
#     search_key = experiment_name + '_config.yaml'
    
#     try:
#         experiment_info_fpath = list(experiment_fpath.glob(search_key))[0]
#     except:
#         logger.error(f'No experiment info file in {experiment_fpath}')
#         fail_signal = signals.FAIL("No experiment info file in the folder")
#         raise fail_signal
    
#     try:
#         experiment_info = OrderedDict(yaml.safe_load(open(experiment_info_fpath, 'rb')))
#         return experiment_info
#     except:
#         logger.error(f'Experiment info file has the wrong name in {experiment_fpath}')
#         fail_signal = signals.FAIL("Experiment info file has the wrong name in the folder")
#         raise fail_signal