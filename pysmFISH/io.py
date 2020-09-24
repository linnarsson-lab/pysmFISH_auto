"""
Set of functions used to load and write data to the shoji database
that will contain all the analysis outputs
"""

from typing import *
import os
import shoji
import numpy as np
from pathlib import Path


from prefect import task
from prefect.engine import signals

from pysmFISH.logger_utils import prefect_logging_setup


@task(name = 'load-preprocessing-parameters')
def load_analysis_parameters(experiment_name:str):
    """
    Function to load all the possible parameters that can be
    selected for the preprocessing 
    NB: If you add additional parameters you are required to
        modify the code and add them to the dictionary
    experiment_name: str
            name of the experiment. Used to access the parameters in the shoji db
    """
    logger = prefect_logging_setup(f'load-analysis-parameters')
    
    try:
        db = shoji.connect()
    except:
        logger.error(f'cannot connect to shoji db')
        err = signals.FAIL(f'cannot connect to shoji db')
        raise err
    else:
        try:
            ws = db.FISH[experiment_name]
        except:
            logger.error(f'cannot connect to the experiment working space')
            err = signals.FAIL(f'cannot connect to the experiment working space')
            raise err
        else:
            try:
                analysis_parameters_ws = db.FISH[experiment_name]['analysis_parameters']
            except:
                logger.error(f'the analysis paramenters workspace is missing')
                err = signals.FAIL(f'the analysis paramenters workspace is missing')
                raise err
            else:
                
                # Collect all parameters
                analysis_parameters = {}
                analysis_parameters['fish'] = {}
                analysis_parameters['small-beads'] = {}
                analysis_parameters['large-beads'] = {}
                analysis_parameters['staining'] = {}
                analysis_parameters['fresh-nuclei'] = {}
                analysis_parameters['BarcodesExtractionResolution'] = analysis_parameters_ws[:].BarcodesExtractionResolution            
                analysis_parameters['RegistrationReferenceHybridization'] = analysis_parameters_ws[:].RegistrationReferenceHybridization

                analysis_parameters['fish']['PreprocessingFishFlatFieldKernel'] = analysis_parameters_ws[:].PreprocessingFishFlatFieldKernel
                analysis_parameters['fish']['PreprocessingFishFilteringSmallKernel'] = analysis_parameters_ws[:].PreprocessingFishFilteringSmallKernel
                analysis_parameters['fish']['PreprocessingFishFilteringLaplacianKernel'] = analysis_parameters_ws[:].PreprocessingFishFilteringLaplacianKernel
                analysis_parameters['fish']['CountingFishMinObjDistance'] = analysis_parameters_ws[:].CountingFishMinObjDistance
                analysis_parameters['fish']['CountingFishMaxObjSize'] = analysis_parameters_ws[:].CountingFishMinObjSize
                analysis_parameters['fish']['CountingFishMinObjSize'] = analysis_parameters_ws[:].CountingFishMinObjSize
                analysis_parameters['fish']['CountingFishNumPeaksPerLabel'] = analysis_parameters_ws[:].CountingFishNumPeaksPerLabel
                
                analysis_parameters['small-beads']['PreprocessingBeadsRegistrationFlatFieldKernel'] = analysis_parameters_ws[:].PreprocessingSmallBeadsRegistrationFlatFieldKernel
                analysis_parameters['small-beads']['PreprocessingBeadsRegistrationFilteringSmallKernel'] = analysis_parameters_ws[:].PreprocessingSmallBeadsRegistrationFilteringSmallKernel
                analysis_parameters['small-beads']['PreprocessingBeadsRegistrationFilteringLaplacianKernel'] = analysis_parameters_ws[:].PreprocessingSmallBeadsRegistrationFilteringLaplacianKernel 
                analysis_parameters['small-beads']['CountingBeadsRegistrationMinObjDistance'] = analysis_parameters_ws[:].CountingSmallBeadsRegistrationMinObjDistance
                analysis_parameters['small-beads']['CountingBeadsRegistratiohMinObjSize'] = analysis_parameters_ws[:].CountingSmallBeadsRegistrationMinObjSize
                analysis_parameters['small-beads']['CountingBeadsRegistrationMaxObjSize'] = analysis_parameters_ws[:].CountingSmallBeadsRegistrationMaxObjSize
                analysis_parameters['small-beads']['CountingBeadsRegistrationNumPeaksPerLabel'] = analysis_parameters_ws[:].CountingSmallBeadsRegistrationNumPeaksPerLabel 
                
                analysis_parameters['large-beads']['PreprocessingBeadsRegistrationFlatFieldKernel'] = analysis_parameters_ws[:].PreprocessingLargeBeadsRegistrationFlatFieldKernel
                analysis_parameters['large-beads']['PreprocessingBeadsRegistrationFilteringSmallKernel'] = analysis_parameters_ws[:].PreprocessingLargeBeadsRegistrationFilteringSmallKernel
                analysis_parameters['large-beads']['PreprocessingBeadsRegistrationFilteringLaplacianKernel'] = analysis_parameters_ws[:].PreprocessingLargeBeadsRegistrationFilteringLaplacianKernel
                analysis_parameters['large-beads']['CountingBeadsRegistrationMinObjDistance'] = analysis_parameters_ws[:].CountingLargeBeadsRegistrationMinObjDistance
                analysis_parameters['large-beads']['CountingBeadsRegistrationMinObjSize'] = analysis_parameters_ws[:].CountingLargeBeadsRegistrationMinObjSize
                analysis_parameters['large-beads']['CountingBeadsRegistrationMaxObjSize'] = analysis_parameters_ws[:].CountingLargeBeadsRegistrationMaxObjSize
                analysis_parameters['large-beads']['CountingBeadsRegistrationNumPeaksPerLabel'] = analysis_parameters_ws[:].CountingLargeBeadsRegistrationNumPeaksPerLabel
                
                analysis_parameters['staining']['PreprocessingStainingFlatFieldKernel'] = analysis_parameters_ws[:].PreprocessingStainingFlatFieldKernel
                
                analysis_parameters['fresh-nuclei']['PreprocessingFreshNucleiLargeKernelSize'] =  analysis_parameters_ws[:].PreprocessingFreshNucleiLargeKernelSize

                return analysis_parameters