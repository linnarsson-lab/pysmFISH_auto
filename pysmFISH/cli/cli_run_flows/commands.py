# import click
# import time
# from pathlib import Path

# from pysmFISH.logger_utils import selected_logger



# from flows.general_flow_single_fov import general_flow


# @click.group('flows_runner')
# def flows_runner():
#     """
#     group of commands used to setup a processing env
#     """
#     pass

# @flows_runner.command('general-flow-runer')

# @click.option('--experiment_fpath', type=str, help='Path to the folder \
#                 where the experiments will be processed')
# @click.option('--run_type', type=str,default='new', help='Type of run \
#                 select between new/re-run')
# @click.option('--parsing_type', type=str, default='original', help='key to select the type of data parsing to run \
#                 - original: parse the data out from robofish system \
#                 - reparsing_from_processing_folder: parse the raw data stored in the \
#                                 experiment folder in the processing HD \
#                 - reparsing_from_storage: parse the raw data stored in the \
#                                 experiment folder in the storage HD \
#                 - no_parsing: skip parsing step')


# @click.option('--fresh_nuclei_processing', type=bool,default=False, 
#                 help='True if you want to process the fresh nuclei staining \
#                     used for eel segmentation')

# @click.option('--raw_data_folder_storage_path', type=str,default='/fish/rawdata', 
#             help='where to store the processed experiment')

# @click.option('--dataset_folder_storage_path', type=str,default='/fish/fish_datasets', 
#             help='where to store the parsed dataset')

# @click.option('--save_intermediate_steps', type=bool,default=False, 
#             help='True if you want to save intermediate processing results')

# @click.option('--store_dataset', type=bool,default=True, 
#             help='True if you want to save the parsed dataset')

            
# def run_general_flow(experiment_fpath:str, 
#                 run_type:str='new', 
#                 parsing_type:str='original',
#                 fresh_nuclei_processing:bool = False,
#                 raw_data_folder_storage_path:str = '/fish/rawdata',
#                 dataset_folder_storage_path:str = '/fish/fish_datasets',
#                 save_intermediate_steps:bool = False,
#                 store_dataset:bool = True):

#     logger = selected_logger()
#     start_time = time.time()
#     click.echo('    ')
#     click.echo('--------------------------------------------------------------')
#     click.echo('Run general processing flow')
#     click.echo('--------------------------------------------------------------')
    
#     general_flow(experiment_fpath, 
#                 run_type = run_type, 
#                 parsing_type=parsing_type,
#                 fresh_nuclei_processing= fresh_nuclei_processing,
#                 raw_data_folder_storage_path= raw_data_folder_storage_path,
#                 dataset_folder_storage_path= dataset_folder_storage_path,
#                 save_intermediate_steps = save_intermediate_steps,
#                 store_dataset=store_dataset)
