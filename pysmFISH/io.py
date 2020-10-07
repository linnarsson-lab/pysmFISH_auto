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
from pysmFISH.utils import convert_to_uint16

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


@task(name = 'save-img-and-metadata')
def save_images_metadata(img_metadata:Tuple):
    """
    Function used to store the images metadata in the shoji database

    Args:
        img_metadata : Tuple 
        contanining (image,metadata) 
            Filtered image uint16 and image_metadata dictionary with the metadata collected duriring the 
            parsing or acquisition of the images
    """

    logger = prefect_logging_setup(f'save-images-metadata')
    metadata = img_metadata[1]
    experiment_name = metadata['experiment_name']
    img = convert_to_uint16(img_metadata[0])
    fov_acquisition_coords = np.array((metadata['fov_acquisition_coords_x'],
                                       metadata['fov_acquisition_coords_y'],
                                       metadata['fov_acquisition_coords_z']),dtype=np.float64)
    img_shape = np.array((metadata['img_height'],
                          metadata['img_width']),np.uint16)
    try:
        db = shoji.connect()
    except:
        logger.error(f'Cannot connect to shoji DB')
        err = signals.FAIL(f'Cannot connect to shoji DB')
        raise err
    else:
        try:
             ws = db.FISH[experiment_name]
        except:
            logger.error(f'experiment workspace missing')
            err = signals.FAIL(f'experiment workspace missing')
            raise err

        else:
            try:
                images_properties_ws = db.FISH[experiment_name]['images_properties']
            except:
                logger.error(f'image properties workspace missing')
                err = signals.FAIL(f'image properties workspace missing')
                raise err
        
            else:
                images_properties_ws.fov.append({
                    'GroupName': np.array([metadata['grp_name']],dtype=np.object).reshape(1,1,1),
                    'FovName' : np.array([metadata['fov_name']],dtype=np.object),
                    'FovNumber' : np.array([metadata['fov_num']],dtype=np.uint16),
                    'AcquistionChannel' : np.array([metadata['channel']],dtype=np.object),
                    'TargetName' : np.array([metadata['target_name']],dtype=np.object).reshape(1,1,1),
                    'ImageShape' : np.array(img_shape,dtype=np.uint16).reshape(1,1,1,2),
                    'PixelMicrons' : np.array([metadata['pixel_microns']],dtype=np.float64).reshape(1,1,1),
                    'HybridizationNumber' : np.array([metadata['hybridization_num']],dtype=np.uint8),
                    'PreprocessedImage': img[None][None][None],
                    'FovCoords': np.array(fov_acquisition_coords,dtype=np.float64).reshape(1,1,1,3),
                    'RegistrationShift' : np.array([0,0],dtype=np.float64).reshape(1,1,1,2),
                    'RegistrationError' : np.array([0],dtype=np.float64).reshape(1,1,1),
                    'StitchingShift' : np.array([0,0],dtype=np.float64).reshape(1,1,1,2),
                    'StitchingError' : np.array([0],dtype=np.float64).reshape(1,1,1),
                    'FieldsOfView': np.array([metadata['fields_of_view']],dtype=np.uint16)

                })


@task(name = 'save-dots-identification')
def save_dots_data(filtered_img_meta:Tuple):
    """
    Function used to store the dots relative data in the shoji database

    Args:
        filtered_img_meta: Tuple
        countains (counts_dict, img_metadata)
        counts_dict : contains data relative to the counted dots
            DotsCoordsFOV
            DotID
            FovNumber
            DotIntensity
            SelectedThreshold
            DotChannel
        img_metadata
    """

    logger = prefect_logging_setup(f'save-dots-identification')
    counts_dict, metadata = filtered_img_meta
    experiment_name = metadata['experiment_name']
    try:
        db = shoji.connect()
    except:
        logger.error(f'Cannot connect to shoji DB')
        err = signals.FAIL(f'Cannot connect to shoji DB')
        raise err
    else:
        try:
             ws = db.FISH[experiment_name]
        except:
            logger.error(f'experiment workspace missing')
            err = signals.FAIL(f'experiment workspace missing')
            raise err

        else:
            try:
                dots_data_ws = db.FISH[experiment_name]['dots_data']
            except:
                logger.error(f'image properties workspace missing')
                err = signals.FAIL(f'image properties workspace missing')
                raise err
        
            else:
                DotCoordsFOV = counts_dict['DotsCoordsFOV'].astype(np.float64)
                DotCoordsFOV = DotCoordsFOV.reshape(DotCoordsFOV.shape[0],1,1,2)
                
                DotIntensity = counts_dict['DotIntensity'].astype(np.float64)
                DotIntensity = DotIntensity[:,np.newaxis,np.newaxis]
                
                SelectedThreshold = counts_dict['SelectedThreshold'].astype(np.float64)
                SelectedThreshold = SelectedThreshold[:,np.newaxis,np.newaxis]
                
                DotChannel = counts_dict['DotChannel'].astype(np.object)
                
                ProcessingType = np.repeat(metadata['processing_type'],DotIntensity.shape[0])
                ProcessingType = ProcessingType.astype(np.object)
                ProcessingType = ProcessingType[:,np.newaxis,np.newaxis]
                
                HybridizationNumber = np.repeat(metadata['hybridization_num'],DotIntensity.shape[0])
                HybridizationNumber = HybridizationNumber.astype(np.uint8)
                
                DotsCoordsRegisteredFOV = np.zeros([DotIntensity.shape[0],1,1,2],dtype=np.float64)
                
                DotsCoordsStitched = np.zeros([DotIntensity.shape[0],1,1,2],dtype=np.float64)

                BarcodeReferenceDotID = np.repeat(str(),DotIntensity.shape[0])
                BarcodeReferenceDotID = BarcodeReferenceDotID.astype(np.object)

                RawBarcode = np.zeros([DotIntensity.shape[0],16],dtype=np.bool)
                
                GeneID = np.repeat(str(),DotIntensity.shape[0])
                GeneID = GeneID.astype(np.object)
                # GeneID = GeneID.reshape(DotIntensity.shape[0],1)
                

                HammingDistanceRawBarcode = np.zeros([DotIntensity.shape[0]],dtype=np.float64)

                dots_data_ws.dots.append({
                    'DotCoordsFOV': DotCoordsFOV,
                    'DotID' : counts_dict['DotID'].astype(np.object),          
                    'FovNumber' :  counts_dict['FovNumber'].astype(np.uint16),                
                    'HybridizationNumber' :  HybridizationNumber,
                    'DotIntensity' : DotIntensity,            
                    'SelectedThreshold' : SelectedThreshold,            
                    'DotChannel' : DotChannel,               
                    'ProcessingType' : ProcessingType,               
                    'DotsCoordsRegisteredFOV' : DotsCoordsRegisteredFOV,
                    'DotsCoordsStitched' : DotsCoordsStitched,
                    'BarcodeReferenceDotID' : BarcodeReferenceDotID,
                    'RawBarcode' : RawBarcode,
                    'GeneID': GeneID,
                    'HammingDistanceRawBarcode' : HammingDistanceRawBarcode

                })  

                del counts_dict
                del metadata