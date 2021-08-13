"""
Module with uility function used to browse throught the different steps 
of the processing
"""
import napari
import sys
import pandas as pd
from pathlib import Path

from pysmFISH import data_models
from pysmFISH import io
from pysmFISH.logger_utils import selected_logger

def visualize_raw_counts(experiment_fpath:str,dataset_name:str,
                        fov_num:int,channel:str):
    """Utility function to visualize the raw count for a FOV.
    The preprocessed images for each round are combined with the
    corresponding counts and plot on different napari layers

    Args:
        experiment_fpath (str): Path to the experiment to process
        dataset_name (str): Name of the dataset
        fov_num (int): Number of the field of view to visualize
        channel (str): Channel to visualize
    """
    
    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)
    preprocessed_img_fpath = experiment_fpath / (experiment_fpath.stem + '_preprocessed_img_data.zarr')
    dataset_fpath = experiment_fpath / dataset_name

    # Load the dataset
    data = data_models.Dataset()
    try:
        data.load_dataset(dataset_fpath)
    except:
        logger.error(f"missing dataset")
        sys.exit(f"missing dataset")
    else:

        fov_dataset = data.dataset.loc[(data.dataset.fov_num == fov_num) &
                                    (data.dataset.channel == channel),:]

        counts_fname = fov_dataset.experiment_name.iloc[0] + '_raw_counts_channel_' + channel + '_fov_' +str(fov_num) + '.parquet'
        counts_fpath = experiment_fpath / 'results' / counts_fname
    
        try:
            all_counts_df = pd.read_parquet(counts_fpath)
        except FileNotFoundError:
            logger.error(f"missing file with raw counts")
            sys.exit(f"missing file with raw counts")
        else:
            vw = napari.Viewer(title='raw counts fov ' + str(fov_num) + 'channel ' + channel)
            for round_num in fov_dataset.round_num.unique():
                fov_subdataset = fov_dataset.loc[fov_dataset.round_num == round_num,:].iloc[0]
                preprocessed_img = io.load_general_zarr(fov_subdataset,
                                                        preprocessed_img_fpath,
                                                        tag='preprocessed_data')
                counts_coords = all_counts_df.loc[all_counts_df.round_num == round_num,['r_px_original','c_px_original']].to_numpy()
                
                vw.add_image(preprocessed_img,name='img round ' + str(round_num), colormap='inferno',
                            visible=False)
                vw.add_points(counts_coords, name= 'counts round '+ str(round_num), face_color='cyan',
                            symbol='+', size=5, opacity=0.5,visible=False)

                    # contrast_limits=(preprocessed_img.min(), preprocessed_img.max())