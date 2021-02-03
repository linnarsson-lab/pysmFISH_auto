import click
import time
from pathlib import Path

from pysmFISH.logger_utils import selected_logger
from pysmFISH.visualization import viz_napari_counts



@click.group('viz')
def viz_utilities():
    """
    group of commands used to visualise data
    """
    pass

@viz_utilities.command('viz-napari-counts')

@click.option('--data_folder_path', type=str, help='Path to processed data. \
                plot the dataframe with decoded in the name')


@click.option('--microscope_tiles_coords_fpath', type=str, help='Path to tiles coords.')


@click.option('--select_genes', type=str, help='level of selection according to hamming distance. \
                the default value is below3Hdistance_genes \
                Possible alternatives: \
                below3Hdistance_genes \
                below2Hdistance_genes \
                zeroHdistance_genes \
                all_Hdistance_genes')

@click.option('--stitching_selected', type=str, help='Select the type of stitching to plot \
                the default value is microscope_stitched \
                Possible alternatives: \
                microscope_stitched \
                XXXXXXXXXXXXXXXXXXX \
                XXXXXXXXXXXXXXXXXXX')

@click.option('--all_genes_visible', type=bool, default=False,help='Select if make all genes visible \
                the default value is False')


def vizcounts(data_folder_path:str,microscope_tiles_coords_fpath:str, 
                select_genes:str,stitching_selected:str,all_genes_visible:bool):
    logger = selected_logger()
    click.echo('    ')
    click.echo('--------------------------------------------------------------')
    click.echo('visualize the stitched counts')
    click.echo('--------------------------------------------------------------')
    viz_napari_counts(data_folder_path,microscope_tiles_coords_fpath,
                        select_genes,stitching_selected,all_genes_visible)
