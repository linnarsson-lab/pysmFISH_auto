from typing import *
import yaml
import time

import numpy as np
import pandas as pd
import sys
import pickle
from pathlib import Path
import pandas as pd

from pysmFISH.io import open_consolidated_metadata
from pysmFISH.logger_utils import selected_logger

class Dataset():
    
    def __init__(self):
        
        self.logger = selected_logger()
        
    def load_dataset(self, dataset_fpath):
        self.dataset = pd.read_parquet(dataset_fpath)
        
    def create_full_dataset_from_files(self, experiment_fpath, 
                                       experiment_info,
                                       parsed_raw_data_fpath,ftype='pkl'): 
        """
        Starting point if you want to work with another type of storage
        but want to use the same dataset structure for downstream
        processing
        """
        
        self.experiment_fpath = Path(experiment_fpath)
        self.experiment_info = experiment_info
        self.parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
    
        date_tag = time.strftime("%y%m%d_%H_%M_%S")
        experiment_name = self.experiment_fpath.stem
        self.dataset_fpath = self.experiment_fpath / (date_tag + '_' + experiment_name + '_dataset.parquet')

        self.dataset = pd.DataFrame()
        all_pickle_list = list(self.parsed_raw_data_fpath.glob('*.' + ftype))
        if len(all_pickle_list):
            for fdata in all_pickle_list:
                single = pickle.load(open(fdata,'rb'))
                fdata_loc = fdata.parent / fdata.stem
                single['raw_data_location'] = fdata_loc.as_posix()
                single_df = pd.DataFrame(single,index=[0])
                self.dataset = pd.concat([self.dataset,single_df],axis=0,ignore_index=True)
            self.dataset.to_parquet(self.dataset_fpath, index=False)
        
        else:
            self.logger.error(f'there are no files with the metadata dictionary')
            sys.exit(f'there are no files with the metadata dictionary')
    
    def create_full_dataset_from_zmetadata(self,experiment_fpath, 
                                       experiment_info,
                                       parsed_raw_data_fpath):  
        
        try:
            consolidated_metadata = open_consolidated_metadata(parsed_raw_data_fpath)
        except:
            self.logger.error(f'consolidated zarr metadata missing or broken')
            sys.exit(f'consolidated zarr metadata missing or broken')
        else:
            self.experiment_fpath = Path(experiment_fpath)
            self.experiment_info = experiment_info
            self.parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
            
            date_tag = time.strftime("%y%m%d_%H_%M_%S")
            experiment_name = self.experiment_fpath.stem
            self.dataset_fpath = self.experiment_fpath / (date_tag + '_' + experiment_name + '_dataset.parquet')
            
            self.dataset = pd.DataFrame()
            for name, grp in consolidated_metadata.items():
                attrs_dict = dict(grp.attrs)
                fdata_loc = self.parsed_raw_data_fpath / attrs_dict['grp_name']
                attrs_dict['raw_data_location'] = fdata_loc.as_posix()
                attrs_df = pd.DataFrame(attrs_dict,index=[0])
                self.dataset = pd.concat([self.dataset,attrs_df],axis=0,ignore_index=True)
            self.dataset.to_parquet(self.dataset_fpath, index=False)
        
            
    def collect_info(self):
        self.list_all_fovs = self.dataset.fov_num.unique()
        self.list_all_channels = self.dataset.channel.unique()
        self.total_rounds = self.dataset.round_num.max()
        self.stitching_channel = self.dataset.iloc[0]['stitching_channel']
        # may not be needed
        self.img_width = self.dataset.iloc[0]['img_width']
        self.img_height = self.dataset.iloc[0]['img_height']
        self.img_zstack = self.dataset.iloc[0]['zstack']
        
    
    def grp_by_channel(self,df):
        subset_df = self.dataset.groupby('channel')
    
    
    def select_fovs_subset(self, df, min_fov, max_fov):
        subset_df = df.loc[(df.fov_num >= min_fov) & (df.fov_num <= max_fov),:]
        return subset_df
    
    
    def select_channel_subset(self,df, channel):
        subset_df = df.loc[(df.channel == channel),:]
        return subset_df
    
    
    def select_round_subset(self,df, round_num):
        subset_df = df.loc[(df.round_num == round_num),:]
        return subset_df
    
    
    def select_specific_fov(self,df, channel,round_num, fov_num):
        subset_df = df.loc[(df.channel == channel) & 
                           (df.round_num == round_num) & 
                           (df.fov_num == fov_num),:]
        return subset_df.squeeze()

    def select_all_imgs_fov(self,df, fovs_list):
        """
        fov is a list
        """
        subset_df = df.loc[(df.fov_num.isin(fovs_list)),:]
        return subset_df
    
    def load_to_numpy(self):
        pass
    
    def load_to_xarray(self):
        pass
    
    


class Output_models():
    """
    Class containing the definition of the output data
    """

    def __init__(self):

        self.output_variables_dot_calling = ['r_px_original','c_px_original',
                                        'dot_id','fov_num','round_num',
                                        'dot_intensity','selected_thr','dot_channel',
                                        'target_name']

        self.output_variables_specific_registration = ['r_px_registered', 'c_px_registered',
                            'r_shift_px','c_shift_px','min_number_matching_dots_registration',
                            'reference_hyb', 'experiment_type','experiment_name','pxl_um',
                            'stitching_type', 'img_width_px','img_height_px','fov_acquisition_coords_x', 
                            'fov_acquisition_coords_y', 'fov_acquisition_coords_z']
        
        self.output_variables_specific_barcodes = ['barcodes_extraction_resolution', 'barcode_reference_dot_id',
                                'raw_barcodes','all_Hdistance_genes','zeroHdistance_genes','below2Hdistance_genes',
                                'below3Hdistance_genes','number_positive_bits','hamming_distance']
    
    
        self.dots_counts_dict = dict.fromkeys(self.output_variables_dot_calling)
        self.dots_counts_df = pd.DataFrame(columns=self.output_variables_dot_calling)

        self.output_variables_registration = self.output_variables_dot_calling + self.output_variables_specific_registration
        self.output_registration_df = pd.DataFrame(columns=self.output_variables_registration)

        self.output_variables_barcodes = self.output_variables_registration + self.output_variables_specific_barcodes
        self.barcode_analysis_df = pd.DataFrame(columns=self.output_variables_barcodes)



















# class shoji_db_fish(Task):
#     """
#     Task used to create the shoji database that will be used to store
#     the output of the analysis. It will also contain the standard settings
#     for automated analysis and the parameter obtained from the experiment.yaml
#     file generated by the machine

#     NB: The settings for the analysis are  adjusted according to the machine 
#         used for acquisition. In our case the machines have different components.

#     Args:
#     -----
#         experiment_info: dict
#             dictionary containing all the info generated by robofish
#     """
    
#     def run(self, experiment_info):
#         """
#         Function used to create the shoji database that will be used to store
#         the output of the analysis. It will also contain the standard settings
#         for automated analysis and the parameter obtained from the experiment.yaml
#         file generated by the machine

#         NB: The settings for the analysis are  adjusted according to the machine 
#             used for acquisition. In our case the machines have different components.

#         Args:
#         -----
#             experiment_info: dict
#                 dictionary containing all the info generated by robofish
#         """
#         experiment_name = experiment_info['EXP_name']
#         try:
#             machine = experiment_info['Machine']
#         except NameError:
#             machine = 'NOT_DEFINED'

#         try:
#             db = shoji.connect()
#         except:
#             self.logger.error(f'Cannot connect to shoji DB')
#             err = signals.FAIL(f'Cannot connect to shoji DB')
#             raise err
#         else:
#             if 'FISH' not in db:
#                 db.FISH = shoji.Workspace()
            
#             if experiment_name not in db.FISH:
#                 db.FISH[experiment_name] = shoji.Workspace()
#                 ws = db.FISH[experiment_name]
#             else:
#                 ws = db.FISH[experiment_name]
#                 self.logger.info('Experiment already present in the database')
            
#             if 'experiment_properties' not in db.FISH[experiment_name]:
#                 db.FISH[experiment_name]['experiment_properties'] = shoji.Workspace()
#                 experiment_properties_ws = db.FISH[experiment_name]['experiment_properties']
#             else:
#                 experiment_properties_ws = db.FISH[experiment_name]['experiment_properties']
            
#             if 'analysis_parameters' not in db.FISH[experiment_name]:
#                 db.FISH[experiment_name]['analysis_parameters'] = shoji.Workspace()
#                 analysis_parameters_ws = db.FISH[experiment_name]['analysis_parameters']
#             else:
#                 analysis_parameters_ws = db.FISH[experiment_name]['analysis_parameters']
            
#             if 'images_properties' not in db.FISH[experiment_name]:
#                 db.FISH[experiment_name]['images_properties'] = shoji.Workspace()
#                 images_properties_ws = db.FISH[experiment_name]['images_properties']
#             else:
#                 images_properties_ws = db.FISH[experiment_name]['images_properties']
            
#             if 'dots_data' not in db.FISH[experiment_name]:
#                 db.FISH[experiment_name]['dots_data'] = shoji.Workspace()
#                 dots_data_ws = db.FISH[experiment_name]['dots_data']
#             else:
#                 dots_data_ws = db.FISH[experiment_name]['dots_data']
            
            
#             # Create the experiment properties workspace
#             try:
#                 experiment_properties_ws.Age = shoji.Tensor("uint8", dims=(), inits=np.array(int(experiment_info['Age']),dtype=np.uint8))
#             except ValueError:
#                 experiment_properties_ws.Age = shoji.Tensor("uint8", dims=(), inits=None)

#             experiment_properties_ws.Barcode =                       shoji.Tensor("bool", dims=(), inits=np.array((experiment_info['Barcode'] == 'True'),dtype=np.bool))
#             experiment_properties_ws.BarcodeLength =                 shoji.Tensor("uint8",dims=(), inits=np.array(experiment_info['Barcode_length'],dtype= np.uint8))
#             experiment_properties_ws.ChamberExp =                    shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Chamber_EXP'],dtype=object)) 
#             experiment_properties_ws.Chemistry =                     shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Chemistry'],dtype=object)) 
#             experiment_properties_ws.CodebookName =                  shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Codebook'],dtype=object)) 
#             experiment_properties_ws.ProbeSetName =                  shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Probe_FASTA_name'],dtype=object)) 
#             experiment_properties_ws.ExperimentType =                shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Experiment_type'],dtype=object))
#             experiment_properties_ws.ExperimentName =                shoji.Tensor("string",dims=(),inits=np.array(experiment_info['EXP_name'],dtype=object))
#             experiment_properties_ws.Machine =                       shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Machine'],dtype=object))
#             experiment_properties_ws.Program =                       shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Program'],dtype=object))
#             experiment_properties_ws.Operator =                      shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Operator'],dtype=object))
#             experiment_properties_ws.Sample =                        shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Sample'],dtype=object))
#             experiment_properties_ws.SectionID =                     shoji.Tensor("string",dims=(),inits=np.array(experiment_info['SectionID'],dtype=object))
#             experiment_properties_ws.Species =                       shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Species'],dtype=object))
#             experiment_properties_ws.Strain =                        shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Strain'],dtype=object))
#             experiment_properties_ws.Tissue =                        shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Tissue'],dtype=object))
#             experiment_properties_ws.Roi =                           shoji.Tensor("string",dims=(),inits=np.array(experiment_info['roi'],dtype=object))
#             experiment_properties_ws.RegionImaged =                  shoji.Tensor("string",dims=(),inits=np.array(experiment_info['RegionImaged'],dtype=object))
#             experiment_properties_ws.SampleOrientation =             shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Orrientation'],dtype=object))
#             experiment_properties_ws.TilesOverlappingPercentage =    shoji.Tensor("float32",dims=(),inits=np.array(np.float32(experiment_info['Overlapping_percentage']),dtype=np.float32))
#             experiment_properties_ws.StitchingChannel =              shoji.Tensor("string",dims=(),inits=np.array(experiment_info['StitchingChannel'],dtype=object))
#             experiment_properties_ws.StitchingType =                 shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Stitching_type'],dtype=object))
#             experiment_properties_ws.DataGenerationDate =            shoji.Tensor("string",dims=(),inits=np.array(experiment_info['Start_date'],dtype=object))



#             analysis_parameters_ws.RegistrationReferenceHybridization = shoji.Tensor("uint8", dims=(), inits=np.array(1,dtype=np.uint8))

#             if machine == 'ROBOFISH1':
#                 analysis_parameters_ws.PreprocessingFishFlatFieldKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([100,100], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingFishFilteringSmallKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([8,8], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingFishFilteringLaplacianKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([1,1], dtype=np.uint16))

#                 analysis_parameters_ws.PreprocessingSmallBeadsRegistrationFlatFieldKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([100,100], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingSmallBeadsRegistrationFilteringSmallKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([8,8], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingSmallBeadsRegistrationFilteringLaplacianKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([1,1], dtype=np.uint16))

#                 analysis_parameters_ws.PreprocessingLargeBeadsRegistrationFlatFieldKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([100,100], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingLargeBeadsRegistrationFilteringSmallKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([8,8], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingLargeBeadsRegistrationFilteringLaplacianKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([1,1], dtype=np.uint16)) 

#                 analysis_parameters_ws.CountingFishMinObjDistance =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingFishMinObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingFishMaxObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(200,dtype=np.uint16))
#                 analysis_parameters_ws.CountingFishNumPeaksPerLabel =  shoji.Tensor("uint16", dims=(), inits=np.array(1,dtype=np.uint16))

#                 analysis_parameters_ws.CountingSmallBeadsRegistrationMinObjDistance =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingSmallBeadsRegistrationMinObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingSmallBeadsRegistrationMaxObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(200,dtype=np.uint16))
#                 analysis_parameters_ws.CountingSmallBeadsRegistrationNumPeaksPerLabel =  shoji.Tensor("uint16", dims=(), inits=np.array(1,dtype=np.uint16))

#                 analysis_parameters_ws.CountingLargeBeadsRegistrationMinObjDistance =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingLargeBeadsRegistrationMinObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingLargeBeadsRegistrationMaxObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(200,dtype=np.uint16))
#                 analysis_parameters_ws.CountingLargeBeadsRegistrationNumPeaksPerLabel =  shoji.Tensor("uint16", dims=(), inits=np.array(1,dtype=np.uint16))

#                 analysis_parameters_ws.BarcodesExtractionResolution =   shoji.Tensor("uint8", dims=(), inits=np.array(3,dtype=np.uint8))                                                                   

#                 analysis_parameters_ws.PreprocessingStainingFlatFieldKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([100,100], dtype=np.uint16))

#                 analysis_parameters_ws.PreprocessingFreshNucleiLargeKernelSize = shoji.Tensor("uint16", dims=(2,), inits=np.array([50,50], dtype=np.uint16))

#             elif machine == 'ROBOFISH2':

#                 analysis_parameters_ws.PreprocessingFishFlatFieldKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([100,100], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingFishFilteringSmallKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([8,8], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingFishFilteringLaplacianKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([1,1], dtype=np.uint16)) 

#                 analysis_parameters_ws.PreprocessingSmallBeadsRegistrationFlatFieldKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([100,100], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingSmallBeadsRegistrationFilteringSmallKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([8,8], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingSmallBeadsRegistrationFilteringLaplacianKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([1,1], dtype=np.uint16)) 

#                 analysis_parameters_ws.PreprocessingLargeBeadsRegistrationFlatFieldKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([100,100], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingLargeBeadsRegistrationFilteringSmallKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([8,8], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingLargeBeadsRegistrationFilteringLaplacianKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([1,1], dtype=np.uint16))

#                 analysis_parameters_ws.CountingFishMinObjDistance =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingFishMinObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingFishMaxObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(200,dtype=np.uint16))
#                 analysis_parameters_ws.CountingFishNumPeaksPerLabel =  shoji.Tensor("uint16", dims=(), inits=np.array(1,dtype=np.uint16))

#                 analysis_parameters_ws.CountingSmallBeadsRegistrationMinObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingSmallBeadsRegistrationMinObjDistance =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingSmallBeadsRegistrationMaxObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(200,dtype=np.uint16))
#                 analysis_parameters_ws.CountingSmallBeadsRegistrationNumPeaksPerLabel =  shoji.Tensor("uint16", dims=(), inits=np.array(1,dtype=np.uint16))

#                 analysis_parameters_ws.CountingLargeBeadsRegistrationMinObjDistance =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingLargeBeadsRegistrationMinObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingLargeBeadsRegistrationMaxObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(200,dtype=np.uint16))
#                 analysis_parameters_ws.CountingLargeBeadsRegistrationNumPeaksPerLabel =  shoji.Tensor("uint16", dims=(), inits=np.array(1,dtype=np.uint16))

#                 analysis_parameters_ws.BarcodesExtractionResolution =   shoji.Tensor("uint8", dims=(), inits=np.array(3,dtype=np.uint8))                                                              

#                 analysis_parameters_ws.PreprocessingStainingFlatFieldKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([100,100], dtype=np.uint16))

#                 analysis_parameters_ws.PreprocessingFreshNucleiLargeKernelSize = shoji.Tensor("uint16", dims=(2,), inits=np.array([50,50], dtype=np.uint16))

#             elif machine == 'NOT_DEFINED':

#                 analysis_parameters_ws.PreprocessingFishFlatFieldKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([100,100], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingFishFilteringSmallKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([8,8], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingFishFilteringLaplacianKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([1,1], dtype=np.uint16)) 

#                 analysis_parameters_ws.PreprocessingSmallBeadsRegistrationFlatFieldKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([100,100], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingSmallBeadsRegistrationFilteringSmallKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([8,8], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingSmallBeadsRegistrationFilteringLaplacianKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([1,1], dtype=np.uint16)) 

#                 analysis_parameters_ws.PreprocessingLargeBeadsRegistrationFilteringSmallKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([8,8], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingLargeBeadsRegistrationFlatFieldKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([100,100], dtype=np.uint16))
#                 analysis_parameters_ws.PreprocessingLargeBeadsRegistrationFilteringLaplacianKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([1,1], dtype=np.uint16))

#                 analysis_parameters_ws.CountingFishMinObjDistance =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingFishMinObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingFishMaxObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(200,dtype=np.uint16))
#                 analysis_parameters_ws.CountingFishNumPeaksPerLabel =  shoji.Tensor("uint16", dims=(), inits=np.array(1,dtype=np.uint16))

#                 analysis_parameters_ws.CountingSmallBeadsRegistrationMinObjDistance =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingSmallBeadsRegistrationMinObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingSmallBeadsRegistrationMaxObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(200,dtype=np.uint16))
#                 analysis_parameters_ws.CountingSmallBeadsRegistrationNumPeaksPerLabel =  shoji.Tensor("uint16", dims=(), inits=np.array(1,dtype=np.uint16))

#                 analysis_parameters_ws.CountingLargeBeadsRegistrationMinObjDistance =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingLargeBeadsRegistrationMinObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(2,dtype=np.uint16))
#                 analysis_parameters_ws.CountingLargeBeadsRegistrationMaxObjSize =  shoji.Tensor("uint16", dims=(), inits=np.array(200,dtype=np.uint16))
#                 analysis_parameters_ws.CountingLargeBeadsRegistrationNumPeaksPerLabel =  shoji.Tensor("uint16", dims=(), inits=np.array(1,dtype=np.uint16))

#                 analysis_parameters_ws.BarcodesExtractionResolution =   shoji.Tensor("uint8", dims=(), inits=np.array(3,dtype=np.uint8))                                                                    

#                 analysis_parameters_ws.PreprocessingStainingFlatFieldKernel = shoji.Tensor("uint16", dims=(2,), inits=np.array([100,100], dtype=np.uint16))

#                 analysis_parameters_ws.PreprocessingFreshNucleiLargeKernelSize = shoji.Tensor("uint16", dims=(2,), inits=np.array([50,50], dtype=np.uint16))


#             # Create dimension and tensors for storage of image related properties   
#             images_properties_ws.fov =                       shoji.Dimension(shape=None)

#             images_properties_ws.hybridization =             shoji.Dimension(shape=1)      
#             images_properties_ws.acquisitioncoords=          shoji.Dimension(shape=3)      
#             images_properties_ws.registrationshiftcoords=    shoji.Dimension(shape=2)
#             images_properties_ws.channel =                   shoji.Dimension(shape=1)
#             images_properties_ws.imageshaperc =              shoji.Dimension(shape=2)

#             # Rank-1 tensors
#             images_properties_ws.FovName = shoji.Tensor("string",dims=('fov',))
#             images_properties_ws.AcquistionChannel = shoji.Tensor("string",dims=('fov',))
#             images_properties_ws.FovNumber = shoji.Tensor("uint16",dims=('fov',))
#             images_properties_ws.HybridizationNumber = shoji.Tensor("uint8",dims=('fov',))

#             # Higher ranking tensors
#             images_properties_ws.FieldsOfView = shoji.Tensor("uint16",dims=('fov',None))
#             images_properties_ws.GroupName = shoji.Tensor("string",dims=('fov','hybridization','channel'))
#             images_properties_ws.TargetName = shoji.Tensor("string",dims=('fov','hybridization','channel'))
#             images_properties_ws.ImageShape = shoji.Tensor("uint16",dims=('fov','hybridization','channel','imageshaperc'))
#             images_properties_ws.PixelMicrons = shoji.Tensor("float64",dims=('fov','hybridization','channel'))
#             images_properties_ws.PreprocessedImage = shoji.Tensor("uint16",dims=('fov','hybridization','channel',None,None))
#             images_properties_ws.FovCoords = shoji.Tensor("float64",dims=('fov','hybridization','channel','acquisitioncoords'))
#             images_properties_ws.RegistrationShift = shoji.Tensor("float64",dims=('fov','hybridization','channel','registrationshiftcoords'))
#             images_properties_ws.RegistrationError = shoji.Tensor("float64",dims=('fov','hybridization','channel'))
#             images_properties_ws.StitchingShift = shoji.Tensor("float64",dims=('fov','hybridization','channel','registrationshiftcoords'))
#             images_properties_ws.StitchingError = shoji.Tensor("float64",dims=('fov','hybridization','channel'))


#             # Create dimension and tensors for dots acquisition
#             # Dimenstion of the tensors
#             dots_data_ws.hybridization =    shoji.Dimension(shape=1)   
#             dots_data_ws.fov =              shoji.Dimension(shape=1)   
#             dots_data_ws.dots =             shoji.Dimension(shape=None)   
#             dots_data_ws.rc =               shoji.Dimension(shape=2)      
#             dots_data_ws.bits =             shoji.Dimension(shape=16)  # Depend on the barcodes used
#             dots_data_ws.gene =             shoji.Dimension(shape=1) 
        
#             # Rank-1 tensors
#             dots_data_ws.DotID =                           shoji.Tensor("string", dims=('dots',))
#             dots_data_ws.FovNumber =                       shoji.Tensor("uint16", dims=('dots',))
#             dots_data_ws.HybridizationNumber =             shoji.Tensor("uint8", dims=('dots',))
#             dots_data_ws.BarcodeReferenceDotID =           shoji.Tensor("string", dims=('dots',))
#             dots_data_ws.DotChannel =                      shoji.Tensor("string", dims=('dots',))
#             dots_data_ws.GeneID =                          shoji.Tensor("string", dims=('dots',))
#             dots_data_ws.HammingDistanceRawBarcode =       shoji.Tensor("float64", dims=('dots',))

#             # add rank-1 tensors for the fitering of the barcodes and genes

#             # Higher ranking tensors
#             dots_data_ws.DotCoordsFOV =                    shoji.Tensor("float64", dims=('dots','fov','hybridization','rc'))
#             dots_data_ws.DotIntensity =                    shoji.Tensor("float64", dims=('dots','fov','hybridization'))
#             dots_data_ws.SelectedThreshold =               shoji.Tensor("float64", dims=('dots','fov','hybridization'))
#             dots_data_ws.ProcessingType =                  shoji.Tensor("string", dims=('dots','fov','hybridization'))
#             dots_data_ws.DotsCoordsRegisteredFOV =         shoji.Tensor("float64", dims=('dots','fov','hybridization','rc'))
#             dots_data_ws.DotsCoordsStitched =              shoji.Tensor("float64", dims=('dots','fov','hybridization','rc'))
#             dots_data_ws.RawBarcode =                      shoji.Tensor("bool", dims=('dots','bits'))