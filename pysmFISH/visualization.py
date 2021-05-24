# import napari
# from dask.distributed import Client
# from dask import dataframe as dd
# from pathlib import Path
# import pandas as pd
# import numpy as np
# import random


# def viz_napari_counts(data_folder_path,
#                       microscope_tiles_coords_fpath,
#                       select_genes='below3Hdistance_genes',
#                       stitching_selected='microscope_stitched', all_genes_visible=False):
    
#     data_folder_path = Path(data_folder_path)
#     tag = '*_decoded_*'
#     r = lambda: random.randint(0,255)
    
#     client = Client()
    
#     microscope_tiles_coords = np.load(microscope_tiles_coords_fpath)
#     fovs = np.arange(microscope_tiles_coords.shape[0])
    
#     r_coords_name = 'r_px_' + stitching_selected
#     c_coords_name = 'c_px_' + stitching_selected
#     decoded_files = dd.read_parquet(data_folder_path / tag)
#     decoded_files = decoded_files.loc[decoded_files.dot_id == decoded_files.barcode_reference_dot_id,:]
# #     decoded_files = decoded_files.dropna(subset=[select_genes])
#     selected = decoded_files.loc[:,[r_coords_name,c_coords_name,select_genes]].compute()
#     gene_grp = selected.groupby(select_genes)
#     client.close()
#     with napari.gui_qt():
#         # create a Viewer and add an image here
#         vw = napari.Viewer(title='stitched counts ')
        
#         _ = vw.add_points(microscope_tiles_coords,
#                             properties={'label': fovs},
#                             text = 'label',
#                             edge_color = 'magenta',
#                             face_color = 'magenta',
#                             symbol = '+',
#                             size = 30,
#                             visible= False
#                             )
        
#         for idx, (gene, coords) in enumerate(gene_grp):
#             coords = coords.loc[:,[r_coords_name,c_coords_name]].to_numpy()
#             col = '#'+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
#             # col = "#%02X%02X%02X" % (r(),r(),r())
#             _ = vw.add_points(coords,name=gene,size=4,symbol='o',visible=all_genes_visible,edge_color=col, face_color=col)