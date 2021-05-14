import pickle
import gc
import dask
import sys
import pandas as pd
import numpy as np
import zarr
from pathlib import Path
from dask import delayed
from dask import dataframe as dd
from dask.base import tokenize


import pysmFISH
from pysmFISH import io
from pysmFISH import fovs_registration
from pysmFISH import barcodes_analysis
from pysmFISH import stitching
from pysmFISH import preprocessing
from pysmFISH import configuration_files
from pysmFISH import utils
from pysmFISH import microscopy_file_parsers
from pysmFISH import data_models
from pysmFISH.logger_utils import selected_logger

def combine_steps(*args):
    pass


def single_fov_round_processing_eel(fov_subdataset,
                                   analysis_parameters,
                                   running_functions,
                                   dark_img,
                                   experiment_fpath,
                                   preprocessed_zarr_fpath,
                                   save_steps_output=False):

    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)

    # Path of directory where to save the intermediate results
    filtered_img_path = experiment_fpath / 'tmp' / 'filtered_images'
    raw_counts_path = experiment_fpath / 'tmp' / 'raw_counts'
    
    experiment_name = fov_subdataset.experiment_name
    pipeline = fov_subdataset.pipeline
    processing_type = fov_subdataset.processing_type
    zarr_grp_name = fov_subdataset.grp_name
    
    raw_data_location = Path(fov_subdataset.raw_data_location)
    parsed_raw_data_fpath = raw_data_location.parent

    if processing_type == 'fish':
        processing_parameters = analysis_parameters['fish']
        filtering_fun = running_functions['fish_channels_preprocessing']
        counting_fun = running_functions['fish_channels_dots_calling']
        

    elif processing_type != 'staining':
        processing_parameters = analysis_parameters[processing_type]
        filtering_fun = running_functions['reference_channels_preprocessing']
        counting_fun = running_functions['reference_channels_dots_calling']


    filt_out = getattr(pysmFISH.preprocessing,filtering_fun)(
                                                    zarr_grp_name,
                                                    parsed_raw_data_fpath,
                                                    processing_parameters,
                                                    dark_img)

    counts = getattr(pysmFISH.dots_calling,counting_fun)(
                                                        filt_out[0][0],
                                                        fov_subdataset,
                                                        processing_parameters)                                              

    if save_steps_output:

        # Save the file as zarr
        store = zarr.DirectoryStore(preprocessed_zarr_fpath)
        root = zarr.group(store=store,overwrite=False)
        tag_name = experiment_name + '_' + fov_subdataset.channel + '_round_' + str(fov_subdataset.round_num) + '_fov_' + str(fov_subdataset.fov_num)
        dgrp = root.create_group(tag_name,overwrite=True)
        for k, v in filt_out[1].items():
            dgrp.attrs[k] = v
        fov_name = 'preprocessed_data_fov_' + str(fov_subdataset.fov_num)
        dgrp.attrs['fov_name'] = fov_name
        img = utils.convert_to_uint16(filt_out[0][-1])
        dset = dgrp.create_dataset(fov_name, data=img, shape=img.shape, chunks=None,overwrite=True)

        # counts.to_parquet(raw_counts_path / (fname + '.parquet'),index=False)

    # return counts, (fov_subdataset.channel,fov_subdataset.round_num,img)
    return counts, filt_out


def single_fov_round_processing_serial_nuclei(fov_subdataset,
                                   analysis_parameters,
                                   running_functions,
                                   dark_img,
                                   experiment_fpath,
                                   preprocessed_image_tag,
                                   save_steps_output=False):

    """
    Process the nuclei for registration
    """
    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)

    # Path of directory where to save the intermediate results
    filtered_img_path = experiment_fpath / 'tmp' / 'filtered_images'
    raw_counts_path = experiment_fpath / 'tmp' / 'raw_counts'
    
    experiment_name = fov_subdataset.experiment_name
    pipeline = fov_subdataset.pipeline
    processing_type = fov_subdataset.processing_type
    zarr_grp_name = fov_subdataset.grp_name
    
    raw_data_location = Path(fov_subdataset.raw_data_location)
    parsed_raw_data_fpath = raw_data_location.parent

    processing_parameters = analysis_parameters[processing_type]

    img = getattr(pysmFISH.preprocessing,running_functions['reference_channels_preprocessing'])(
                                                                        zarr_grp_name,
                                                                        parsed_raw_data_fpath,
                                                                        processing_parameters)

    if save_steps_output:

        # Save the file as zarr
        store = zarr.DirectoryStore(preprocessed_zarr_fpath)
        root = zarr.group(store=store,overwrite=False)
        tag_name = experiment_name + '_' + fov_subdataset.channel + '_round_' + str(fov_subdataset.round_num) + '_fov_' + str(fov_subdataset.fov_num)
        dgrp = root.create_group(tag_name,overwrite=True)
        for k, v in filt_out[1].items():
            dgrp.attrs[k] = v
        fov_name = 'preprocessed_data_fov_' + str(fov_subdataset.fov_num)
        dgrp.attrs['fov_name'] = fov_name
        img = utils.convert_to_uint16(img)
        dset = dgrp.create_dataset(fov_name, data=img, shape=img.shape, chunks=None,overwrite=True)

    return (img,fov_subdataset)


def processing_barcoded_eel_fov_graph(experiment_fpath,analysis_parameters,
                                    running_functions, tile_corners_coords_pxl,metadata,
                                    grpd_fovs,save_intermediate_steps, 
                                    preprocessed_image_tag, client ):
        """ 
        This method create a processing graph XXXXXXXXX

        """
        experiment_fpath = Path(experiment_fpath)
        io.create_empty_zarr_file(experiment_fpath.as_posix(), preprocessed_image_tag)
        preprocessed_zarr_fpath = experiment_fpath / (experiment_fpath.stem + '_' + preprocessed_image_tag + '.zarr')

        utils.create_dark_img(experiment_fpath,metadata)

        dark_img = preprocessing.load_dark_image(experiment_fpath)
        
        # did this conversion to avoid to pass self to dask
        # analysis_parameters = analysis_parameters
        # running_functions = running_functions
        # tile_corners_coords_pxl = tile_corners_coords_pxl
        
        
        dark_img = delayed(dark_img)
        analysis_parameters = delayed(analysis_parameters)
        running_functions = delayed(running_functions)
        tile_corners_coords_pxl = delayed(tile_corners_coords_pxl)

        codebook = configuration_files.load_codebook(experiment_fpath,metadata)
        codebook_df = delayed(codebook)
        
        all_processing = []
        all_filtered_images = []
        all_fovs = list(grpd_fovs.groups.keys())
        chunks = [all_fovs[x:x+30] for x in range(0, len(all_fovs), 30)]
        for chunk in chunks:
            for fov_num in chunk:
                group = grpd_fovs.get_group(fov_num)
        # for fov_num, group in grpd_fovs:
                all_counts_fov = []
                for index_value, fov_subdataset in group.iterrows():
                    round_num = fov_subdataset.round_num
                    channel = fov_subdataset.channel
                    fov = fov_subdataset.fov_num
                    experiment_name = fov_subdataset.experiment_name
                    dask_delayed_name = 'filt_count_' +experiment_name + '_' + channel + \
                                    '_round_' + str(round_num) + '_fov_' +str(fov) + '-' + tokenize()
                    fov_out = delayed(single_fov_round_processing_eel, name=dask_delayed_name,nout=2)(fov_subdataset,
                                                analysis_parameters,
                                                running_functions,
                                                dark_img,
                                                experiment_fpath,
                                                preprocessed_zarr_fpath,
                                                save_steps_output=save_intermediate_steps,
                                                dask_key_name=dask_delayed_name)
                    counts, filt_out = fov_out[0], fov_out[1]
                    all_counts_fov.append(counts)
                    if channel != fov_subdataset.stitching_channel:
                        all_filtered_images.append(filt_out)
                
                # Modify for channels name

                name = 'concat_' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize()
                all_counts_fov = delayed(pd.concat,name=name)(all_counts_fov,axis=0,ignore_index=True)
                
                name = 'register_' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize()
                registered_counts = delayed(fovs_registration.beads_based_registration,name=name)(all_counts_fov,
                                                    analysis_parameters)

                name = 'decode_' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize()
                decoded = delayed(barcodes_analysis.extract_barcodes_NN_fast,name=name)(registered_counts, 
                                                                            analysis_parameters,codebook_df)                                                        
                
                name = 'stitch_to_mic_coords_' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize()  
                stitched_coords = delayed(stitching.stitch_using_microscope_fov_coords_new,name=name)(decoded[1],tile_corners_coords_pxl)
                
                name = 'save_df_' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize() 
                saved_file = delayed(stitched_coords.to_parquet,name=name)(Path(experiment_fpath) / 'results'/ (experiment_name + \
                                '_decoded_fov_' + str(fov) + '.parquet'),index=False)
            
                name = 'register_and_combine_filt_imgs' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize() 
                
                combined_images = delayed(fovs_registration.combine_register_filtered_images,name=name)(all_filtered_images,stitched_coords,
                                                                                                fov_subdataset.stitching_channel)


                name = 'extract_barcodes_int' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize()
                barcodes_max = delayed(barcodes_analysis.extract_dots_images,name=name)(stitched_coords,combined_images,
                                                                                        experiment_fpath)


                name = 'extract_bit_flip_direction' +experiment_name + '_' + channel + '_' \
                                    + '_fov_' +str(fov) + '-' + tokenize()
                flip_direction = delayed(barcodes_analysis.define_flip_direction,name=name)(codebook,
                                                                                        experiment_fpath,
                                                                                        stitched_coords)

                end = delayed(combine_steps)(saved_file,barcodes_max,flip_direction)

                all_processing.append(end)
            
            _ = dask.compute(*all_processing)

        io.consolidate_zarr_metadata(preprocessed_zarr_fpath)
 
        # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        selected_Hdistance = 3 / metadata['barcode_length']
        stitching_selected = 'microscope_stitched'
        io.simple_output_plotting(experiment_fpath, stitching_selected, selected_Hdistance, client)
        # ----------------------------------------------------------------  


def processing_serial_fish_fov_graph(experiment_fpath,analysis_parameters,
                                    running_functions, tile_corners_coords_pxl,metadata,
                                    grpd_fovs,save_intermediate_steps, 
                                    preprocessed_zarr_fpath, client):
        """ 
        This method create a processing graph XXXXXXXXX

        """

        experiment_fpath = Path(experiment_fpath)
        io.create_empty_zarr_file(experiment_fpath, preprocessed_image_tag)
        preprocessed_zarr_fpath = experiment_fpath / (experiment_fpath.stem + '_' + preprocessed_image_tag + '.zarr')

        utils.create_dark_img(experiment_fpath,metadata)


        dark_img = preprocessing.load_dark_image(experiment_fpath)
        
        # did this conversion to avoid to pass self to dask
        analysis_parameters = analysis_parameters
        running_functions = running_functions
        tile_corners_coords_pxl = tile_corners_coords_pxl
        
        dark_img = delayed(dark_img)
        analysis_parameters = delayed(analysis_parameters)
        running_functions = delayed(running_functions)
        tile_corners_coords_pxl = delayed(tile_corners_coords_pxl)

        all_processing = []

        for fov_num, group in grpd_fovs:
            all_counts_fov = []
            all_nuclei_fov = []
            for index_value, fov_subdataset in group.iterrows():
                round_num = fov_subdataset.round_num
                channel = fov_subdataset.channel
                fov = fov_subdataset.fov_num
                experiment_name = fov_subdataset.experiment_name
                processing_type = fov_subdataset.processing_type

                if processing_type == 'nuclei':
                    dask_delayed_name = 'filt_' +experiment_name + '_' + channel + \
                                    '_round_' + str(round_num) + '_fov_' +str(fov) + '-' + tokenize()

                    out_nuclei = delayed(fov_processing.single_fov_round_processing_serial_nuclei,name=dask_delayed_name)(fov_subdataset,
                                            analysis_parameters,
                                            running_functions,
                                            dark_img,
                                            experiment_fpath,
                                            preprocessed_zarr_fpath,
                                            save_steps_output=save_intermediate_steps,
                                            dask_key_name=dask_delayed_name)
                    all_nuclei_fov.append(out_nuclei)

                else:
                    dask_delayed_name = 'filt_count_' +experiment_name + '_' + channel + \
                                    '_round_' + str(round_num) + '_fov_' +str(fov) + '-' + tokenize()
                    counts = delayed(fov_processing.single_fov_round_processing_eel,name=dask_delayed_name)(fov_subdataset,
                                                analysis_parameters,
                                                running_functions,
                                                dark_img,
                                                experiment_fpath,
                                                preprocessed_zarr_fpath,
                                                save_steps_output=save_intermediate_steps,
                                                dask_key_name=dask_delayed_name)
                    all_counts_fov.append(counts)
            
            name = 'concat_fish_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            all_counts_fov = delayed(pd.concat,name=name)(all_counts_fov,axis=0,ignore_index=True)
            
            name = 'create_nuclei_stack' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            filtered_nuclei_stack = delayed(utils.combine_filtered_images,name=name)(all_nuclei_fov,experiment_fpath,metadata)

            name = 'register_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()
            registered_counts = delayed(fovs_registration.nuclei_based_registration,name=name)(all_counts_fov,
                                                filtered_nuclei_stack,
                                                analysis_parameters)
                                                                                                
            name = 'stitch_to_mic_coords_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize()  
            stitched_coords = delayed(stitching.stitch_using_microscope_fov_coords_new,name=name)(registered_counts,tile_corners_coords_pxl)
            
            name = 'save_df_' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize() 
            saved_file = delayed(stitched_coords.to_parquet,name=name)(Path(experiment_fpath) / 'results'/ (experiment_name + \
                            '_decoded_fov_' + str(fov) + '.parquet'),index=False)
        
            all_processing.append(saved_file) 
        _ = dask.compute(*all_processing)

        io.consolidate_zarr_metadata(preprocessed_zarr_fpath)

        # # ----------------------------------------------------------------
        # GENERATE OUTPUT FOR PLOTTING
        stitching_selected = 'microscope_stitched'
        io.simple_output_plotting_serial(experiment_fpath, stitching_selected, client)
        # ----------------------------------------------------------------  
  







def single_fov_fresh_tissue_beads(processing_tag,
                                   fov_subdataset,
                                   analysis_parameters,
                                   running_functions,
                                   dark_img,
                                   experiment_fpath,
                                   preprocessed_zarr_fpath,
                                   save_steps_output=True):

    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)

    experiment_name = fov_subdataset.experiment_name
    zarr_grp_name = fov_subdataset.grp_name

    parsed_raw_data_fpath = experiment_fpath / 'fresh_tissue'/ (experiment_name +'_img_data.zarr')
    
    processing_parameters = analysis_parameters['fresh-tissue'][processing_tag]

    if processing_tag == 'beads':
        filtering_fun = running_functions['fresh_sample_reference_preprocessing']
        counting_fun = running_functions['fresh_sample_reference_dots_calling']


        filt_out = getattr(pysmFISH.preprocessing,filtering_fun)(
                                                        zarr_grp_name,
                                                        parsed_raw_data_fpath,
                                                        processing_parameters,
                                                        dark_img)

        counts = getattr(pysmFISH.dots_calling,counting_fun)(
                                                            filt_out[0][0],
                                                            fov_subdataset,
                                                            processing_parameters)                                              

    elif processing_tag == 'nuclei':
        filtering_fun = running_functions['fresh_sample_nuclei_preprocessing']
        filt_out = getattr(pysmFISH.preprocessing,filtering_fun)(
                                                        zarr_grp_name,
                                                        parsed_raw_data_fpath,
                                                        processing_parameters)


    if save_steps_output:

        # Save the file as zarr
        store = zarr.DirectoryStore(preprocessed_zarr_fpath)
        root = zarr.group(store=store,overwrite=False)
        tag_name = experiment_name + '_fresh_tissue_' + processing_tag + '_fov_' + str(fov_subdataset.fov_num)
        dgrp = root.create_group(tag_name,overwrite=True)
        for k, v in filt_out[1].items():
            dgrp.attrs[k] = v
        fov_name = 'preprocessed_data_fov_' + str(fov_subdataset.fov_num)
        dgrp.attrs['fov_name'] = fov_name
        img = utils.convert_to_uint16(filt_out[0][-1])
        dset = dgrp.create_dataset(fov_name, data=img, shape=img.shape, chunks=None,overwrite=True)


    if processing_tag == 'beads':
        return counts, filt_out

    elif processing_tag == 'nuclei':
        return filt_out


def process_fresh_sample_graph(experiment_fpath, running_functions, 
                                analysis_parameters, client, tag_ref_beads, 
                                tag_nuclei,parsing=True):
    logger = selected_logger()
    all_parsing = []
    
    if parsing:
        presence_nuclei = 0
        presence_beads = 0
        all_fresh_tissue_fpath = list((Path(experiment_fpath) / 'fresh_tissue').glob('*.nd2'))
        if all_fresh_tissue_fpath:
            for fpath in all_fresh_tissue_fpath:
                if tag_ref_beads in fpath.stem:
                    parsed_beads_fpath = experiment_fpath / 'fresh_tissue'/ (fpath.stem +'_img_data.zarr')
                    parsing_future = client.submit(microscopy_file_parsers.nikon_nd2_parser_simple_mfov,fpath)
                    all_parsing.append(parsing_future)
                    presence_beads = 1
                elif tag_nuclei in fpath.stem:
                    parsed_nuclei_fpath = experiment_fpath / 'fresh_tissue'/ (fpath.stem +'_img_data.zarr')
                    parsing_future = client.submit(microscopy_file_parsers.nikon_nd2_parser_simple_mfov,fpath)
                    all_parsing.append(parsing_future)
                    presence_nuclei = 1

            if presence_beads and presence_nuclei:
                _ = client.gather(all_parsing)
                io.consolidate_zarr_metadata(parsed_beads_fpath)
                io.consolidate_zarr_metadata(parsed_nuclei_fpath)
            else:
                if not presence_beads:
                    logger.error(f'missing fresh-tissue beads file')
                    sys.exit(f'missing fresh-tissue beads file')
                elif not presence_nuclei:
                    logger.error(f'missing fresh-tissue nuclei file')
                    sys.exit(f'missing fresh-tissue nuclei file')
                else:
                    logger.error(f'missing fresh-tissue beads and nuclei files')
                    sys.exit(f'missing fresh-tissue beads and nuclei files')

    else:
        all_fresh_tissue_fpath = list((Path(experiment_fpath) / 'fresh_tissue').glob('*.zarr'))
        if all_fresh_tissue_fpath:
            for fpath in all_fresh_tissue_fpath:
                if tag_ref_beads in fpath.stem:
                    parsed_beads_fpath = fpath
                    presence_beads = 1
                elif tag_nuclei in fpath.stem:
                    parsed_nuclei_fpath = fpath
                    presence_nuclei = 1
              
        if not presence_beads:
            logger.error(f'missing fresh-tissue beads parsed file')
            sys.exit(f'missing fresh-tissue beads parsed file')
        elif not presence_nuclei:
            logger.error(f'missing fresh-tissue nuclei parsed file')
            sys.exit(f'missing fresh-tissue nuclei parsed file')
        elif (not presence_nuclei) & (not presence_beads    ):
            logger.error(f'missing fresh-tissue beads and nuclei parsed files')
            sys.exit(f'missing fresh-tissue beads and nuclei parsed files')  

  
    # Create dataset
    
    ds_beads = data_models.Dataset()
    ds_nuclei = data_models.Dataset()
    ds_beads.create_full_dataset_from_zmetadata(parsed_beads_fpath)
    ds_nuclei.create_full_dataset_from_zmetadata(parsed_nuclei_fpath)

    beads_grpd_fovs = ds_beads.dataset.groupby('fov_num')
    nuclei_grpd_fovs = ds_nuclei.dataset.groupby('fov_num')
    
    # In this case I fake a dark image. It must be collected from the 
    # robofish system
    img_width = ds_nuclei.dataset.iloc[0].img_width
    img_height = ds_nuclei.dataset.iloc[0].img_height
    dark_img = np.zeros([img_width,img_height])

    base_path = experiment_fpath / 'fresh_tissue'
    # nuclei_base = parsed_nuclei_fpath.stem.split('_img_data.zarr')[0]
    # beads_base = parsed_beads_fpath.stem.split('_img_data.zarr')[0]
    nuclei_filtered_fpath = base_path /  (base_path.stem + '_nuclei_preprocessed_img_data.zarr')
    io.create_empty_zarr_file(base_path.as_posix(), tag='nuclei_preprocessed_img_data')
    beads_filtered_fpath = base_path /  (base_path.stem + '_beads_preprocessed_img_data.zarr')
    io.create_empty_zarr_file(base_path.as_posix(), tag='beads_preprocessed_img_data')
    

    all_counts_beads = []
    all_processing_nuclei = []
    processing_tag = 'beads'
    for fov_num, group in beads_grpd_fovs:
        for index_value, fov_subdataset in group.iterrows():
            round_num = fov_subdataset.round_num
            channel = fov_subdataset.channel
            fov = fov_subdataset.fov_num
            experiment_name = fov_subdataset.experiment_name
            dask_delayed_name = 'filt_count_beads_fov'+ str(fov) + '_' + tokenize()
            fov_out = delayed(single_fov_fresh_tissue_beads, name=dask_delayed_name)(
                                            processing_tag,
                                            fov_subdataset,
                                            analysis_parameters,
                                            running_functions,
                                            dark_img,
                                            experiment_fpath,
                                            preprocessed_zarr_fpath=beads_filtered_fpath,
                                            save_steps_output=True,
                                            dask_key_name=dask_delayed_name)
            counts, filt_out = fov_out[0], fov_out[1]
            all_counts_beads.append(counts)
    

    name = 'concat_all_counts_beads_fresh_tissue'+ '-' + tokenize()
    all_counts_fov = delayed(pd.concat,name=name)(all_counts_beads,axis=0,ignore_index=True)

    # Add registration and recalculation of all the coords

    name = 'save_df_beads_fresh_tissue' +experiment_name + '_' + channel + '_' \
                                + '_fov_' +str(fov) + '-' + tokenize() 
    saved_file = delayed(all_counts_fov.to_parquet,name=name)(Path(experiment_fpath) / 'results'/ (experiment_name + \
                    '_counts_beads_fresh_tissue.parquet'),index=False)

    processing_tag='nuclei'
    for fov_num, group in nuclei_grpd_fovs:
        for index_value, fov_subdataset in group.iterrows():
            round_num = fov_subdataset.round_num
            channel = fov_subdataset.channel
            fov = fov_subdataset.fov_num
            experiment_name = fov_subdataset.experiment_name
            dask_delayed_name = 'filt_nuclei_fov'+ str(fov) + '_' + tokenize()
            fov_out = delayed(single_fov_fresh_tissue_beads, name=dask_delayed_name)(
                                            processing_tag,
                                            fov_subdataset,
                                            analysis_parameters,
                                            running_functions,
                                            dark_img,
                                            experiment_fpath,
                                            preprocessed_zarr_fpath=nuclei_filtered_fpath,
                                            save_steps_output=True,
                                            dask_key_name=dask_delayed_name)
            
            all_processing_nuclei.append(fov_out)

 
    end = delayed(combine_steps)(saved_file,all_processing_nuclei)
      
    _ = dask.compute(end)
    

